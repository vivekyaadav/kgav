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
    gate_census,
    host_protein_ids,
    ingest_host_targets,
    load_mechanisms,
    query_host_activities,
)
from kgav.schema import load_schema, report_violations


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path,
                    default=root / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--spine", type=Path, default=root / "data/releases/v0.1-spine")
    ap.add_argument("--host", type=Path, default=root / "data/releases/v0.1-host")
    ap.add_argument("--out", type=Path, default=root / "data/releases/v0.1-hosttargets")
    ap.add_argument("--min-pchembl", type=float, default=6.0)
    ap.add_argument("--min-phase", type=int, default=2)
    ap.add_argument("--gate-census", action="store_true",
                    help="count what each SQL gate discards. A full scan of "
                         "~24.5M activity rows, about 76s; off by default so "
                         "a rebuild does not pay for it.")
    args = ap.parse_args()

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

    print(f"\nquerying activities (pChEMBL >= {args.min_pchembl}, "
          f"max_phase >= {args.min_phase})")
    rows = query_host_activities(args.db, args.min_pchembl, args.min_phase)
    print(f"  {len(rows):,} activities")

    em = Emit()
    s = ingest_host_targets(em, rows, mechanisms, proteins, "infores:chembl")

    print(f"\n{s['edges']:,} TARGETS edges over {len(em.nodes):,} compounds")
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
    orcs = root / "data/releases/v0.1-orcs/edges.jsonl"
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
