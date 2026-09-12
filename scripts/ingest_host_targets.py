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
from kgav.host_targets import (
    host_protein_ids,
    ingest_host_targets,
    load_mechanisms,
    query_host_activities,
)
from kgav.schema import load_schema


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
    args = ap.parse_args()

    print("loading host layer")
    index = AliasIndex.from_releases(args.spine, args.host)
    proteins = host_protein_ids(index)
    print(f"  {len(proteins):,} host proteins in the graph")

    print("\nloading drug_mechanism action types")
    mechanisms = load_mechanisms(args.db)
    print(f"  {len(mechanisms):,} (compound, target) pairs with a stated action")

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
    print("\n  dropped:")
    for k in ("no_inchikey", "target_not_in_graph", "unconvertible_units", "duplicate"):
        if s[k]:
            print(f"    {k:<22} {s[k]:>7,}")

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
    if violations:
        print(f"\nSCHEMA FAIL — {len(violations)} violations")
        for v in violations[:8]:
            print(f"  {v}")
        return 1
    print("\nschema: PASS")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in em.nodes.values()))
    (args.out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in em.edges))
    print(f"{len(em.nodes):,} compounds, {len(em.edges):,} edges -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
