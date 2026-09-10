#!/usr/bin/env python3
"""Validate a KG release against the frozen schema.

    python scripts/validate_release.py data/releases/v0.1/

Exits non-zero if any violation is found. This is the CI gate.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.schema import load_schema


def _read(base: Path, stem: str) -> list[dict]:
    pq, jl = base / f"{stem}.parquet", base / f"{stem}.jsonl"
    if pq.exists():
        import pandas as pd
        df = pd.read_parquet(pq)
        for col in ("properties", "qualifiers", "publications"):
            if col in df.columns:
                df[col] = df[col].apply(lambda v: json.loads(v) if isinstance(v, str) else v)
        return df.to_dict("records")
    if jl.exists():
        with open(jl, encoding="utf-8") as fh:
            return [json.loads(line) for line in fh if line.strip()]
    raise SystemExit(f"no {stem}.parquet or {stem}.jsonl in {base}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("release_dir", type=Path)
    ap.add_argument("--schema", type=Path, default=None)
    ap.add_argument("--max-examples", type=int, default=5)
    args = ap.parse_args()

    schema = load_schema(args.schema)
    nodes = _read(args.release_dir, "nodes")
    edges = _read(args.release_dir, "edges")

    print(f"schema v{schema.version} | {len(nodes):,} nodes | {len(edges):,} edges")
    violations = schema.validate_batch(nodes, edges)

    if not violations:
        print("\nPASS — no violations")
        print(f"held-out predicates (strip from training graph): {sorted(schema.held_out_predicates())}")
        return 0

    counts = Counter(v.code for v in violations)
    print(f"\nFAIL — {len(violations):,} violations across {len(counts)} codes\n")
    for code, n in counts.most_common():
        print(f"  {code:<28} {n:>8,}")
        for v in [x for x in violations if x.code == code][: args.max_examples]:
            print(f"      {v.ref}: {v.message[:140]}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
