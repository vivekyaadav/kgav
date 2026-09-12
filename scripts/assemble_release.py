#!/usr/bin/env python3
"""Day 8 driver: assemble layers, validate the whole graph, write a manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.assemble import (
    LAYER_PRECEDENCE,
    assemble,
    connectivity_report,
    orphan_nodes,
)
from kgav.schema import load_schema


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

    layers = {name: args.releases / f"v0.1-{name}" for name in LAYER_PRECEDENCE}
    layers = {k: v for k, v in layers.items() if v.exists()}
    if not layers:
        raise SystemExit(f"no v0.1-* release directories in {args.releases}")

    schema_early = load_schema()
    identity = {ec.predicate: ec.identity_qualifiers
                for ec in schema_early.edge_classes if ec.identity_qualifiers}
    print(f"assembling {len(layers)} layers: {', '.join(layers)}")
    if identity:
        print(f"  identity-defining qualifiers: {identity}")
    a = assemble(layers, identity)

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
        print(f"\n  first {min(len(a.conflicts), 10)} conflicts:")
        for c in a.conflicts[:10]:
            print(f"    {c}")

    dangling = a.dangling()
    if dangling:
        print(f"\nDANGLING — {len(dangling):,} edges reference a node that does not exist")
        for s, p, o in dangling[:8]:
            missing = s if s not in a.nodes else o
            print(f"  {s} -{p}-> {o}   (missing: {missing})")
        return 1

    orphans = orphan_nodes(a)
    if orphans:
        print("\nnodes with no edges (not an error; a large count in a class that "
              "should be connected means a join failed):")
        for c, n in orphans.most_common():
            print(f"  {c:<16} {n:>8,}")

    print("\nconnectivity by predicate:")
    for p, combos in sorted(connectivity_report(a).items()):
        for combo, n in combos.most_common(3):
            print(f"  {p:<32} {combo:<34} {n:>8,}")

    schema = schema_early
    print("\nvalidating assembled graph")
    violations = schema.validate_batch(a.nodes.values(), a.edges.values())
    if violations:
        print(f"SCHEMA FAIL — {len(violations):,} violations")
        from collections import Counter
        for code, n in Counter(v.code for v in violations).most_common():
            print(f"  {code:<28} {n:>7,}")
        for v in violations[:6]:
            print(f"    {v}")
        return 1
    print("schema: PASS")

    cut, hubs = a.hub_threshold(schema.retrieval_limits.get("hub_percentile", 99.9))
    print(f"\nhub cut at p{schema.retrieval_limits.get('hub_percentile')}: "
          f"degree > {cut} ({len(hubs)} nodes)")
    deg = a.degree()
    for nid in sorted(hubs, key=lambda n: -deg[n])[:10]:
        props = a.nodes[nid].get("properties", {})
        label = props.get("gene_symbol") or props.get("label") or nid
        print(f"  {str(label)[:24]:<26} {a.nodes[nid]['class']:<14} {deg[nid]:>7,}")

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
