#!/usr/bin/env python3
"""Driver: drug -> human protein TARGETS edges. Unblocks M2 through M6."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.aliases import AliasIndex
from kgav.emit import Emit
from kgav.provenance import chembl_release, write_source
from kgav.host_targets import (
    DEFAULT_MIN_PCHEMBL,
    DEFAULT_MIN_PHASE,
    DEFAULT_PHASE_GATE,
    PHASE_GATE_POLICIES,
    gate_census,
    host_protein_ids,
    ingest_host_targets,
    known_inchikeys,
    load_mechanisms,
    query_host_activities,
)
from kgav.pipeline import layer_dir, layers_before
from kgav.schema import load_schema, report_violations


def build_parser(root: Path | None = None) -> argparse.ArgumentParser:
    """Extracted from main() so a test can read the defaults.

    --phase-gate carries a decision (fix-list B1 option D), and a default
    that only a 7-minute rebuild can observe is a default nothing checks.
    """
    root = Path(__file__).resolve().parents[1] if root is None else root
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path,
                    default=root / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--spine", type=Path, default=layer_dir(root, "spine"))
    ap.add_argument("--host", type=Path, default=layer_dir(root, "host"))
    ap.add_argument("--out", type=Path, default=layer_dir(root, "hosttargets"))
    ap.add_argument("--releases", type=Path, default=root / "data/releases")
    ap.add_argument("--min-pchembl", type=float, default=DEFAULT_MIN_PCHEMBL)
    ap.add_argument("--min-phase", type=int, default=DEFAULT_MIN_PHASE)
    ap.add_argument("--phase-gate", choices=PHASE_GATE_POLICIES,
                    default=DEFAULT_PHASE_GATE,
                    help="floor (default): the phase floor for everyone, "
                         "which is how v0.1 is built. relaxed: the floor OR "
                         "the graph already holds the compound -- measured "
                         "and declined (fix-list B1), kept runnable, NOT "
                         "evaluated; its admission rule correlates with the "
                         "labels.")
    ap.add_argument("--gate-census", action="store_true",
                    help="count what each SQL gate discards. A full scan of "
                         "~24.5M activity rows, about 76s; off by default so "
                         "a rebuild does not pay for it.")
    return ap


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    args = build_parser(root).parse_args()

    print("loading host layer")
    index = AliasIndex.from_releases(args.spine, args.host)
    proteins = host_protein_ids(index)
    print(f"  {len(proteins):,} host proteins in the graph")

    print("\nloading drug_mechanism action types")
    mechanisms = load_mechanisms(args.db)
    print(f"  {len(mechanisms):,} (compound, target) pairs with a stated action")

    # WHAT THE FOUR SQL GATES DISCARD (fix-list E7). Behind a flag because
    # it is a full scan of 24.5M activity rows -- 76s against an 8s ingest
    # -- and the answer only moves when ChEMBL does. The figures are
    # recorded in host_targets.gate_census' docstring.
    if args.gate_census:
        print("\n  gate census (full scan, ~76s):")
        c = gate_census(args.db, args.min_pchembl, args.min_phase)
        total = 0
        for k in sorted(c):
            if k.startswith("dropped_") or k == "kept":
                print(f"    {k:<42} {c[k]:>12,}")
                total += c[k]
        if total != c["rows_scanned"]:
            print(f"    ! ACCOUNTING DOES NOT CLOSE: {c['rows_scanned']:,} "
                  f"scanned, {total:,} accounted for")
            return 1
        print(f"    {'(sums to the scan)':<42} {c['rows_scanned']:>12,}")

    # WHICH COMPOUNDS THE GRAPH ALREADY HAS, read from the layers built
    # before this one -- which pipeline.layers_before names, so this cannot
    # drift when a layer moves. AND IT WILL NOT FIND THEM ALL from here:
    # 660 of the 1,327 compounds that gain a route come from ncats, which
    # is built after this layer, so --phase-gate relaxed sees the 667
    # chembl ones. hosttargets was moved after ncats to fix that, and moved
    # back when B1 was declined; the count printed below is the only honest
    # way to say which half ran.
    known = None
    if args.phase_gate == "relaxed":
        before = layers_before("hosttargets")
        known = known_inchikeys(args.releases, before)
        print("\nphase gate: RELAXED for compounds already in the graph "
              "-- not the shipped policy (fix-list B1)")
        print(f"  {len(known):,} compound(s) from {len(before)} earlier "
              f"layer(s): {', '.join(before)}")
    else:
        print(f"\nphase gate: floor only (max_phase >= {args.min_phase}) "
              f"-- how v0.1 was built")

    print(f"\nquerying activities (pChEMBL >= {args.min_pchembl}, "
          f"max_phase >= {args.min_phase}"
          f"{' OR already in the graph' if known else ''})")
    rows = query_host_activities(args.db, args.min_pchembl, args.min_phase,
                                 known=known)
    print(f"  {len(rows):,} activities")
    # WHAT THE RELAXATION ADMITTED, COUNTED. A policy whose effect is not
    # reported is indistinguishable from one that did nothing -- the rule
    # L10 settled for dose_series_policy.
    if known:
        relaxed_rows = [r for r in rows
                        if r.get("max_phase") is None
                        or r["max_phase"] < args.min_phase]
        relaxed_cmpds = {r["inchikey"] for r in relaxed_rows if r.get("inchikey")}
        print(f"  {len(relaxed_rows):,} of them admitted ONLY by the "
              f"relaxation, over {len(relaxed_cmpds):,} compound(s)")

    em = Emit()
    s = ingest_host_targets(em, rows, mechanisms, proteins, "infores:chembl")

    print(f"\n{s['edges']:,} TARGETS edges over {len(em.nodes):,} compounds")
    # WHICH KEY DECIDED EACH NODE ID (fix-list L27), same four counters and
    # the same reason as ingest_chembl.py: they were populated from the day
    # the identity repair landed and printed by neither driver.
    identity = {k: s[k] for k in sorted(s) if k.startswith("identity_")}
    if identity:
        print("  which key decided the node id:")
        for k, n in identity.items():
            print(f"    {k.removeprefix('identity_'):<12} {n:>8,}")
        print(f"    {'(total)':<12} {sum(identity.values()):>8,}")
    print("  direction (from drug_mechanism where stated):")
    for k in ("direction_inhibitor", "direction_agonist", "direction_modulator",
              "direction_unknown"):
        if s[k]:
            print(f"    {k.replace('direction_',''):<12} {s[k]:>7,}")
    # WHAT KIND OF ASSAY, PER EDGE. "binding" was a constant on all 18,954.
    print("  assay_type (from ChEMBL's own assay_type column):")
    for k in sorted(k for k in s if k.startswith("assay_type_")):
        print(f"    {k.removeprefix('assay_type_'):<12} {s[k]:>7,}")
    # AND THE LIABILITY PANELS, BY NAME. A and T coded rows are ADMET and
    # toxicity screens -- CYP3A4, CYP2D6, CYP2C9, hERG -- which become
    # TARGETS edges indistinguishable from therapeutic target engagement
    # and feed M2 through M6. Not dropped: inhibiting CYP3A4 at 100 nM is
    # a true statement, and whether a repurposing route may walk through
    # it is a judgement for whoever reads this number.
    liab = s["row_chembl_assay_type_A"] + s["row_chembl_assay_type_T"]
    if liab:
        print(f"\n  ! {liab:,} source row(s) are ADMET or toxicity panels "
              f"(ChEMBL A/T): mostly\n    CYP and hERG liability screens. "
              f"They are emitted as TARGETS edges and\n    reach M2-M6 like "
              f"any other. assay_type records which, nothing filters them.")

    # PUBLICATION COVERAGE, PRINTED. Every TARGETS edge in v0.1 carried no
    # publication, because the query selected docs.year and not
    # docs.pubmed_id. Nothing reported the absence, so it took an audit of
    # the assembled release to find it. A run now states its own coverage.
    cited, uncited = s["with_publication"], s["no_publication"]
    if cited or uncited:
        print(f"\n  publications: {cited:,} of {cited + uncited:,} edges "
              f"cite at least one PMID ({cited / (cited + uncited):.1%})")
        if s["duplicate_added_publication"]:
            print(f"    {s['duplicate_added_publication']:,} citation(s) came "
                  f"from a row that collapsed into an existing edge")

    if s["relation_imputed"]:
        print(f"\n  {s['relation_imputed']:,} row(s) had no usable "
              f"standard_relation and were recorded as exact")

    print("\n  dropped:")
    for k in ("no_inchikey", "target_not_in_graph", "unconvertible_units"):
        if s[k]:
            print(f"    {k:<22} {s[k]:>7,}")
    # NOT UNDER "dropped" ANY MORE. A duplicate row is collapsed into the
    # edge the first row created and contributes its PMID to it; calling
    # that a drop is what made losing its citation look like nothing.
    if s["duplicate"]:
        print(f"\n  collapsed into an existing edge: {s['duplicate']:,}")

    # How many TARGETS land on a coronavirus dependency factor -- this is the
    # number that decides whether M2 produces paths at all.
    dep: set[str] = set()
    orcs = layer_dir(root, "orcs") / "edges.jsonl"
    if orcs.exists():
        for line in orcs.read_text().splitlines():
            if not line.strip():
                continue
            e = json.loads(line)
            if (e["predicate"] == "HOST_FACTOR_FOR"
                    and e["qualifiers"].get("direction") == "dependency"):
                dep.add(e["subject"])
        hit = {e["object"] for e in em.edges if e["object"] in dep}
        m2 = sum(1 for e in em.edges if e["object"] in dep)
        print(f"\nM2 reachability: {len(hit):,} targeted proteins are coronavirus "
              f"dependency factors\n  {m2:,} TARGETS edges feed a single-hop M2 path")

    all_nodes = list(index.nodes.values()) + list(em.nodes.values())
    violations = load_schema().validate_batch(all_nodes, em.edges)
    if report_violations(violations):
        return 1

    # DATE PRECISION, REPORTED (fix-list C10). A year-only date is written
    # as a year; the placeholder is counted apart, because "undated" spread
    # across seven sources is seven gaps with seven remedies.
    if s["dated_year_only"] or s["undated_placeholder"]:
        tot = s["dated_year_only"] + s["undated_placeholder"]
        print(f"\n  dates: {s['dated_year_only']:,} at YEAR precision, "
              f"{s['undated_placeholder']:,} undated "
              f"({s['undated_placeholder'] / tot:.1%})")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in em.nodes.values()))
    (args.out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in em.edges))
    # WHICH RELEASE OF THE SOURCE (fix-list V2/V3/V5). Read from the
    # file on disk, never from a register kept by hand: a release named in
    # a path or a filename is evidence, one typed into a config is a claim.
    write_source(args.out, chembl=chembl_release(args.db))
    print(f"{len(em.nodes):,} compounds, {len(em.edges):,} edges -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
