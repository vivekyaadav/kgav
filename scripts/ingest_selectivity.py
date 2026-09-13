#!/usr/bin/env python3
"""Build the selectivity layer from same-document EC50/CC50 pairs."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.emit import Emit
from kgav.schema import load_schema
from kgav.selectivity import SI_THRESHOLD, ingest_selectivity, query_pairs


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path,
                    default=root / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--config", type=Path, default=root / "config" / "viruses.yaml")
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--out", type=Path, default=root / "data/releases/v0.1-selectivity")
    ap.add_argument("--threshold", type=float, default=SI_THRESHOLD)
    args = ap.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    taxon_map: dict[str, str] = {}
    for v in cfg["viruses"]:
        taxon_map[str(v["taxon"])] = str(v["taxon"])
        for iso in v.get("isolate_taxa") or []:
            taxon_map[str(iso)] = str(v["taxon"])

    print("pairing antiviral EC50/IC50 with CC50/TC50 within the same document")
    rows = query_pairs(args.db, sorted(taxon_map))
    print(f"  {len(rows):,} matched (compound, document, virus) rows")

    em = Emit()
    s = ingest_selectivity(em, rows, taxon_map, "infores:chembl", args.threshold)

    total = s["selective"] + s["cytotoxic"]
    print(f"\n{s['edges']:,} compounds with a verified selectivity index")
    if total:
        print(f"  {s['selective']:,} selective (SI >= {args.threshold:g})  "
              f"{100*s['selective']/total:.0f}%")
        print(f"  {s['cytotoxic']:,} cytotoxic (SI < {args.threshold:g})  "
              f"{100*s['cytotoxic']/total:.0f}%")
    for k in ("no_inchikey", "unmapped_taxon", "si_uncomputable"):
        if s[k]:
            print(f"  dropped {k}: {s[k]:,}")

    per_virus: Counter = Counter()
    for e in em.edges:
        per_virus[e["object"]] += 1
    if per_virus:
        print("\nper virus:")
        for tid, n in per_virus.most_common():
            print(f"  {tid:<22} {n:>6,}")

    nodes = [json.loads(x) for x in (args.release / "nodes.jsonl").read_text().splitlines()
             if x.strip()]
    known = {n["id"] for n in nodes}
    dangling = [e for e in em.edges if e["subject"] not in known or e["object"] not in known]
    if dangling:
        # A compound paired here but absent from the release means the
        # chemistry layer filtered it out; emitting the edge would dangle.
        print(f"\n  {len(dangling):,} edges reference compounds absent from the "
              f"release and are dropped")
        em.edges = [e for e in em.edges if e not in dangling]

    violations = load_schema().validate_batch(nodes, em.edges)
    if violations:
        print(f"\nSCHEMA FAIL — {len(violations)} violations")
        for v in violations[:8]:
            print(f"  {v}")
        return 1
    print("\nschema: PASS")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("")
    (args.out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in em.edges))
    print(f"{len(em.edges):,} edges -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
