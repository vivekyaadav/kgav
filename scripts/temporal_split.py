#!/usr/bin/env python3
"""Build the temporal split, run the leakage audits, and evaluate on it."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.baselines import (
    Graph,
    combine,
    degree_ranking,
    dwpc_scores,
    hits_at_k,
    mrr,
    rank,
)
from kgav.schema import load_schema
from kgav.temporal import audit, build_split, write_split

MIN_TRUSTWORTHY_EXPECTATION = 1.0


def evaluate(scores: dict[str, float], pos: set[str], k: int) -> dict:
    r = rank(scores)
    p = len(pos & set(scores))
    exp = k * p / len(r) if r else 0.0
    h = hits_at_k(r, pos, k)
    return {"pool": len(r), "pos": p, "hits": h, "expected": exp,
            "lift": (h / exp) if exp else 0.0, "mrr": mrr(r, pos),
            "trustworthy": exp >= MIN_TRUSTWORTHY_EXPECTATION}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--out", type=Path, default=root / "data/releases/v0.1-train2021")
    ap.add_argument("--results", type=Path, default=root / "data/results")
    ap.add_argument("--cutoff", type=int, default=2021)
    ap.add_argument("--undated", default="include",
                    choices=["include", "exclude", "computed_only"])
    ap.add_argument("-k", type=int, default=100)
    args = ap.parse_args()

    s = load_schema()
    held = s.held_out_predicates()
    print(f"cutoff {args.cutoff} | undated policy: {args.undated}")

    split = build_split(args.release, args.cutoff, held, args.undated)
    print(f"\ntraining graph: {len(split.train_edges):,} edges")
    for key in ("edge_train", "edge_dropped_future", "undated_kept",
                "undated_kept_computed", "undated_dropped",
                "label_train", "label_test", "label_undated", "label_censored"):
        if split.stats[key]:
            print(f"  {key:<24} {split.stats[key]:>8,}")

    print("\nleakage audits:")
    problems = audit(split, args.release, held)
    fatal = [p for p in problems if p.startswith(("L3", "L5", "label leak"))]
    for p in problems:
        print(f"  {'FAIL' if p in fatal else 'note'}  {p}")
    if not problems:
        print("  all clean")
    if fatal:
        print("\nrefusing to evaluate on a leaking split")
        return 1

    write_split(split, args.release, args.out)
    print(f"\nwrote training graph -> {args.out}")

    g = Graph.load(args.out, skip_predicates=held)
    print(f"loaded {len(g.nodes):,} nodes | {len(g.drugs()):,} drugs")

    dwpc = dwpc_scores(g, s.metapaths, s.retrieval_limits)
    deg = degree_ranking(g)
    results: dict[str, dict] = {}

    for virus, test_pos in sorted(split.test_labels.items(),
                                  key=lambda kv: -len(kv[1])):
        label = (g.nodes.get(virus, {}).get("properties", {}) or {}).get("label", virus)
        print(f"\n=== {label}: {len(test_pos):,} post-{args.cutoff} compounds")
        reachable = test_pos & set(g.nodes)
        print(f"    {len(reachable):,} exist in the pre-{args.cutoff} graph "
              f"(ceiling on recall)")
        if not reachable:
            print("    unevaluable")
            continue
        per: dict[str, dict] = {}
        for name in sorted(dwpc):
            sc = {d: by[virus] for d, by in dwpc[name].items() if virus in by}
            if sc:
                per[name] = evaluate(sc, test_pos, args.k)
        per["COMBINED"] = evaluate(combine(dwpc, virus), test_pos, args.k)
        per["degree"] = evaluate(deg, test_pos, args.k)
        print(f"    {'scorer':<10} {'pool':>7} {'pos':>5} {'hits':>5} "
              f"{'exp':>7} {'lift':>7} {'mrr':>7}")
        for name, m in per.items():
            flag = "" if m["trustworthy"] else " (!)"
            print(f"    {name:<10} {m['pool']:>7,} {m['pos']:>5,} {m['hits']:>5} "
                  f"{m['expected']:>7.1f} {m['lift']:>7.2f} {m['mrr']:>7.4f}{flag}")
        results[label] = per

    args.results.mkdir(parents=True, exist_ok=True)
    (args.results / f"temporal_{args.cutoff}_{args.undated}.json").write_text(
        json.dumps({"cutoff": args.cutoff, "undated_policy": args.undated,
                    "stats": dict(split.stats), "audits": problems,
                    "results": results}, indent=2))
    print(f"\nwrote {args.results}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
