#!/usr/bin/env python3
"""Day 7 driver: chemistry from ChEMBL 37."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.aliases import AliasIndex
from kgav.chembl import (
    build_viral_target_index,
    ingest_activities,
    query_activities,
)
from kgav.emit import Emit
from kgav.schema import load_schema


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path,
                    default=root / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--config", type=Path, default=root / "config" / "viruses.yaml")
    ap.add_argument("--spine", type=Path, default=root / "data/releases/v0.1-spine")
    ap.add_argument("--host", type=Path, default=root / "data/releases/v0.1-host")
    ap.add_argument("--out", type=Path, default=root / "data/releases/v0.1-chembl")
    args = ap.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    taxon_map: dict[str, str] = {}
    for v in cfg["viruses"]:
        taxon_map[str(v["taxon"])] = str(v["taxon"])
        for iso in v.get("isolate_taxa") or []:
            taxon_map[str(iso)] = str(v["taxon"])

    print("loading spine + host")
    index = AliasIndex.from_releases(args.spine, args.host)
    viral_targets = build_viral_target_index(index, taxon_map)
    n_targets = sum(len(v) for v in viral_targets.values())
    print(f"  {len(index):,} nodes | {n_targets} viral chain targets across "
          f"{len(viral_targets)} viruses")

    print(f"\nquerying {args.db.name}")
    rows = query_activities(args.db, sorted(taxon_map))
    print(f"  {len(rows):,} activities")

    em = Emit()
    s = ingest_activities(em, rows, viral_targets, taxon_map, "infores:chembl")

    print(f"\n{s['protein_edges']:,} INHIBITS + {s['organism_edges']:,} "
          f"HAS_ANTIVIRAL_ACTIVITY_AGAINST")
    print(f"  {s['censored']:,} censored values kept with their relation "
          f"(hard negatives for Day 12)")
    print("\n  target resolution from assay description:")
    for k in sorted(k for k in s if k.startswith("resolved_")):
        print(f"    {k.replace('resolved_',''):<10} {s[k]:>7,}")
    print(f"    {'unresolved':<10} {s['unresolved']:>7,}")
    print("\n  dropped:")
    for k in ("no_inchikey", "unmapped_taxon", "unconvertible_units",
              "no_target_node", "duplicate"):
        if s[k]:
            print(f"    {k:<22} {s[k]:>7,}")

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

    per_virus: Counter = Counter()
    for e in em.edges:
        if e["predicate"] == "HAS_ANTIVIRAL_ACTIVITY_AGAINST":
            per_virus[e["object"]] += 1
    if per_virus:
        print("\ncell-based antiviral edges per virus:")
        for tid, n in per_virus.most_common():
            lbl = (index.nodes.get(tid, {}).get("properties", {}) or {}).get("label", tid)
            print(f"  {lbl:<14} {n:>7,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
