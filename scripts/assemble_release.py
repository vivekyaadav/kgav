#!/usr/bin/env python3
"""Day 8 driver: assemble layers, validate the whole graph, write a manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.assemble import (
    LAYER_PRECEDENCE,
    assemble,
    connectivity_report,
    orphan_detail,
)
from kgav.schema import load_schema, print_capped, report_violations


#: Viral nodes that legitimately have no edges, each with the reason it is
#: allowed to. ANY OTHER VIRAL ORPHAN FAILS THE BUILD.
#:
#: A host protein with no edges is expected -- the proteome is ingested
#: whole. A viral one means the spine built a node and no layer joined it,
#: which is a failed join, and until the host/viral split landed it was
#: invisible inside "Protein 634" in every release.
#:
#: The entry below is NOT deleted, and the reason is the spine's own, from
#: the comment on the pp1a/pp1ab alias stub: "the loser's id must remain
#: resolvable -- VirHostNet and IntAct cite pp1a chain ids freely, and
#: without this alias those interactions silently fail to join". The
#: uncleaved pp1a chain id is precisely the kind an external database
#: cites, so the node earns its place by being resolvable rather than by
#: being connected. Declared, so it is a decision on the record instead of
#: an anomaly nobody can read.
EXPECTED_VIRAL_ORPHANS = {
    "UniProtKB:PRO_0000449634":
        "SARS-CoV-2 uncleaved pp1a. No gene symbol, so the spine writes no "
        "ENCODED_BY, and no experiment attaches to the uncleaved form -- "
        "evidence lands on the nsps inside it. Kept so an external citation "
        "of the pp1a chain id resolves.",
}


def undeclared_viral_orphans(orphans: list[dict]) -> list[dict]:
    """Viral orphans nobody has accounted for. These fail the build."""
    return [n for n in orphans
            if (n.get("properties") or {}).get("is_viral")
            and n["id"] not in EXPECTED_VIRAL_ORPHANS]


def _rel(p: Path, root: Path) -> str:
    try:
        return str(p.relative_to(root))
    except ValueError:
        return str(p)


def _hash(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--releases", type=Path, default=root / "data" / "releases")
    ap.add_argument("--out", type=Path, default=root / "data" / "releases" / "v0.1")
    ap.add_argument("--version", default="v0.1")
    args = ap.parse_args()

    declared = {name: args.releases / f"v0.1-{name}" for name in LAYER_PRECEDENCE}
    layers = {k: v for k, v in declared.items() if v.exists()}
    if not layers:
        raise SystemExit(f"no v0.1-* release directories in {args.releases}")

    # A LISTED LAYER WITH NO DIRECTORY IS THE QUIETEST FAILURE HERE. The graph
    # assembles, the schema passes, every count looks reasonable, and the
    # metapaths that needed that layer reach nothing -- which is how the folds
    # layer went missing long enough for M8 to look like a modelling problem.
    # Said out loud now, because a skipped layer is almost always an ingest
    # that failed earlier in the same session.
    for name, path in declared.items():
        if name not in layers:
            print(f"  ! {name}: {path} not found — SKIPPED. Any metapath "
                  f"needing it will\n    reach nothing, and nothing further "
                  f"will say so.")

    schema_early = load_schema()
    identity = {ec.predicate: ec.identity_qualifiers
                for ec in schema_early.edge_classes if ec.identity_qualifiers}
    print(f"assembling {len(layers)} layers: {', '.join(layers)}")
    if identity:
        print(f"  identity-defining qualifiers: {identity}")
    # Printed because a guard nobody can see is a guard nobody can trust. If
    # this line is absent, two screens of one compound are silently merging
    # and LAYER_PRECEDENCE is deciding which measurement survives.
    src_ident = schema_early.source_identity_predicates()
    if src_ident:
        print(f"  identified partly by source: {sorted(src_ident)}")
    derived = {ec.predicate: schema_early.derived_absence_qualifiers(ec.predicate)
               for ec in schema_early.edge_classes}
    derived = {p: d for p, d in derived.items() if d}
    if derived:
        print(f"  derived at assembly (never ingested): {derived}")
    # SPELLING IS NOT DISAGREEMENT. Supplied here rather than inside
    # assemble because deciding it needs rdkit, the optional `chem` extra;
    # without it the comparator is simply absent and every difference is
    # reported as a conflict, which is the previous behaviour.
    equivalent = {}
    try:
        from kgav.normalize_chem import same_structure
        same_structure("C", "C")     # import rdkit NOW, not mid-assembly
    except ImportError:
        print("  rdkit absent: two spellings of one SMILES will be counted "
              "as a conflict")
    else:
        equivalent["smiles"] = same_structure
    a = assemble(layers, identity, schema_early.derive_qualifiers,
                 schema_early.source_identity_predicates(), equivalent)

    print(f"\n{len(a.nodes):,} nodes, {len(a.edges):,} edges")
    for c, n in a.by_class().most_common():
        print(f"  {c:<16} {n:>8,}")
    print()
    for p, n in a.by_predicate().most_common():
        print(f"  {p:<32} {n:>8,}")

    if a.stats:
        print("\nmerges and conflicts:")
        for k, v in sorted(a.stats.items()):
            print(f"  {k:<34} {v:>7,}")
    if a.conflicts:
        # "first 10 conflicts" with no total read the same at 10 and at 400.
        print(f"\n  {len(a.conflicts):,} property conflict(s):")
        print_capped(a.conflicts, 10, unit="more conflict(s)")

    dangling = a.dangling()
    if dangling:
        print(f"\nDANGLING — {len(dangling):,} edges reference a node that does not exist")
        print_capped(
            dangling, 8,
            lambda k: f"{k[0]} -{k[1]}-> {k[2]}   "
                      f"(missing: {k[0] if k[0] not in a.nodes else k[2]})",
            indent="  ", unit="more dangling edge(s)")
        return 1

    untouched_ids = set(a.nodes) - {x for k in a.edges for x in (k[0], k[2])}
    # The KEY is the id. Reading it from the stored dict instead would make
    # the gate depend on every layer having written an "id" field into its
    # own node body, which is true today and is not the authority.
    orphan_rows = [{**a.nodes[i], "id": i} for i in untouched_ids]
    detail = orphan_detail(a)
    if detail:
        print("\nnodes with no edges — a HOST protein with none is expected, "
              "the proteome is\ningested whole; a VIRAL one means the spine "
              "built it and no layer joined it:")
        # "<- a join failed" printed directly above "DECLARED ... kept
        # deliberately" is a report contradicting itself one line apart,
        # which is the thing this release keeps being audited for. The flag
        # only fires where the taxon has an UNDECLARED orphan.
        open_taxa = {(n.get("properties") or {}).get("taxon_id")
                     for n in undeclared_viral_orphans(orphan_rows)}
        viral_orphans = 0
        for c, by in sorted(detail.items(),
                            key=lambda kv: -sum(kv[1].values())):
            print(f"  {c:<16} {sum(by.values()):>8,}")
            for label, n in by.most_common():
                flag = ""
                if label.startswith("viral"):
                    viral_orphans += n
                    tax = label.removeprefix("viral ")
                    flag = ("   <- a join failed" if tax in open_taxa
                            else "   <- declared below")
                print(f"    {label:<30} {n:>8,}{flag}")
        if not viral_orphans:
            print("  No viral node is orphaned, so every count above is the "
                  "benign case.")
        else:
            declared = [(nid, why)
                        for nid, why in sorted(EXPECTED_VIRAL_ORPHANS.items())
                        if nid in untouched_ids]
            for nid, why in declared:
                print(f"\n  DECLARED  {nid}")
                for line in textwrap.wrap(why, 70):
                    print(f"            {line}")

    undeclared = undeclared_viral_orphans(orphan_rows)
    if undeclared:
        print(f"\n{len(undeclared)} VIRAL node(s) with no edge and no "
              f"declared reason. A viral orphan is a\nfailed join: the spine "
              f"built it and no layer connected it. Run\nscripts/"
              f"diagnose_orphans.py to see which of the three failures it "
              f"is, then\neither fix the join or add it to "
              f"EXPECTED_VIRAL_ORPHANS with the reason.")
        print_capped(undeclared, 20,
                     lambda n: f"{n['id']}   "
                               f"{(n.get('properties') or {}).get('taxon_id')}",
                     unit="more undeclared viral orphan(s)")
        return 1

    print("\nconnectivity by predicate:")
    for p, combos in sorted(connectivity_report(a).items()):
        for combo, n in combos.most_common(3):
            print(f"  {p:<32} {combo:<34} {n:>8,}")

    schema = schema_early
    print("\nvalidating assembled graph")
    violations = schema.validate_batch(a.nodes.values(), a.edges.values())
    if report_violations(violations):
        return 1

    # Hub exclusion is class-scoped. A global cut would exclude the virus
    # nodes (SARS-CoV-2 has degree 10,679) and the viral replicase -- but those
    # are where M1-M6 END. High degree there is structural, not uninformative.
    lim = schema.retrieval_limits
    cut = lim.get("hub_degree_cut", 460)
    classes = set(lim.get("hub_exclude_classes") or ["Protein"])
    include_viral = bool(lim.get("hub_exclude_viral", False))
    deg = a.degree()
    hubs = [nid for nid, k in deg.items()
            if k > cut and a.nodes[nid]["class"] in classes
            and (include_viral or not a.nodes[nid]["properties"].get("is_viral"))]
    print(f"\nhub exclusion: {sorted(classes)} with degree > {cut}"
          f"{'' if include_viral else ', viral excluded from the cut'}")
    print(f"  {len(hubs)} nodes excluded as path connectors")
    def _hub(nid):
        props = a.nodes[nid].get("properties", {})
        label = props.get("gene_symbol") or props.get("label") or nid
        return f"{str(label)[:24]:<26} {deg[nid]:>7,}"

    print_capped(sorted(hubs, key=lambda n: -deg[n]), 10, _hub,
                 unit="more excluded hub(s)")

    a.write(args.out)
    manifest = a.manifest(
        {name: {"path": _rel(p, root),
                "nodes_sha": _hash(p / "nodes.jsonl") if (p / "nodes.jsonl").exists() else None,
                "edges_sha": _hash(p / "edges.jsonl") if (p / "edges.jsonl").exists() else None}
         for name, p in layers.items()},
        schema.version)
    manifest["release"] = args.version
    manifest["hub_cut_degree"] = cut
    (args.out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nwrote {args.out} + MANIFEST.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
