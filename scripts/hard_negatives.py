#!/usr/bin/env python3
"""Evaluate against MEASURED negatives, cross-sectionally and temporally.

The question every previous run could not answer: does the graph rank actives
above compounds that were assayed and found INACTIVE, rather than merely above
compounds nobody tested?
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.baselines import Graph, combine, degree_ranking, dwpc_scores
from kgav.labels import build_labels, evaluate_against_negatives
from kgav.schema import load_schema

MIN_CLASS_SIZE = 10


def run(g: Graph, schema, labels, tag: str, k_min: int = MIN_CLASS_SIZE) -> dict:
    dwpc = dwpc_scores(g, schema.metapaths, schema.retrieval_limits)
    deg = degree_ranking(g)
    out: dict[str, dict] = {}

    viruses = sorted(labels.positives, key=lambda v: -len(labels.positives[v]))
    for virus in viruses:
        pos, neg = labels.positives.get(virus, set()), labels.negatives.get(virus, set())
        label = (g.nodes.get(virus, {}).get("properties", {}) or {}).get("label", virus)
        print(f"\n=== {label} [{tag}]")
        print(f"    {len(pos):,} measured active | {len(neg):,} measured INACTIVE")
        if len(pos) < k_min or len(neg) < k_min:
            print(f"    skipped: need at least {k_min} of each class")
            continue

        scorers = {name: {d: by[virus] for d, by in dwpc[name].items() if virus in by}
                   for name in sorted(dwpc)}
        scorers["COMBINED"] = combine(dwpc, virus)
        scorers["degree"] = deg

        print(f"    {'scorer':<10} {'pos':>5} {'neg':>6} {'AUC':>7} "
              f"{'cov_pos':>8} {'cov_neg':>8}")
        per: dict[str, dict] = {}
        for name, sc in scorers.items():
            if not sc:
                continue
            m = evaluate_against_negatives(sc, pos, neg)
            per[name] = m
            if not m["evaluable"]:
                print(f"    {name:<10} {m['n_pos']:>5} {m['n_neg']:>6} "
                      f"{'-':>7}   one class empty in this pool")
                continue
            print(f"    {name:<10} {m['n_pos']:>5} {m['n_neg']:>6} {m['auc']:>7.3f} "
                  f"{m['coverage_pos']:>7.1%} {m['coverage_neg']:>7.1%}")
        out[label] = per
    return out


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--train", type=Path, default=root / "data/releases/v0.1-train2021")
    ap.add_argument("--results", type=Path, default=root / "data/results")
    ap.add_argument("--cutoff", type=int, default=2021)
    ap.add_argument("--inactive-above-nm", type=float, default=10_000.0)
    ap.add_argument("--selectivity", default="all",
                    choices=["all", "selective-only", "verified-only"],
                    help="all: every measured active (the original protocol). "
                         "selective-only: actives with a verified SI >= 10. "
                         "verified-only: actives with ANY verified SI, "
                         "selective or not -- isolates the effect of "
                         "restricting to compounds with matched data from the "
                         "effect of the selectivity filter itself.")
    ap.add_argument("--si-threshold", type=float, default=10.0)
    args = ap.parse_args()

    schema = load_schema()
    held = schema.held_out_predicates()

    all_labels = build_labels(args.release, inactive_above_nm=args.inactive_above_nm)

    if args.selectivity != "all":
        from kgav.selectivity import selectivity_index_by_compound
        si = selectivity_index_by_compound(args.release)
        print(f"selectivity filter '{args.selectivity}': "
              f"{len(si):,} (compound, virus) pairs have a verified index")
        for virus, pos in list(all_labels.positives.items()):
            if args.selectivity == "verified-only":
                keep = {d for d in pos if (d, virus) in si}
            else:
                keep = {d for d in pos if si.get((d, virus), -1) >= args.si_threshold}
            dropped = len(pos) - len(keep)
            all_labels.positives[virus] = keep
            if dropped:
                print(f"  {virus}: {len(pos):,} -> {len(keep):,} actives "
                      f"({dropped:,} excluded)")
    print(f"labels from the full release (inactive if censored above "
          f"{args.inactive_above_nm:,.0f} nM):")
    for key in ("active", "inactive", "undecidable", "ambiguous_dropped"):
        if all_labels.stats[key]:
            print(f"  {key:<20} {all_labels.stats[key]:>7,}")

    print("\n" + "=" * 66)
    print("CROSS-SECTIONAL: full graph, all measured compounds")
    print("=" * 66)
    g_full = Graph.load(args.release, skip_predicates=held)
    cross = run(g_full, schema, all_labels, "cross-sectional")

    temporal: dict = {}
    if args.train.exists():
        print("\n" + "=" * 66)
        print(f"TEMPORAL: pre-{args.cutoff} graph, post-{args.cutoff} measurements")
        print("=" * 66)
        test_labels = build_labels(args.release, year_from=args.cutoff + 1,
                                   inactive_above_nm=args.inactive_above_nm)
        # The selectivity filter must apply here too. Filtering only the
        # cross-sectional labels would report a temporal number computed on a
        # different positive set from the one it is being compared against.
        if args.selectivity != "all":
            from kgav.selectivity import selectivity_index_by_compound
            si_t = selectivity_index_by_compound(args.release)
            for virus, pos in list(test_labels.positives.items()):
                if args.selectivity == "verified-only":
                    keep = {d for d in pos if (d, virus) in si_t}
                else:
                    keep = {d for d in pos
                            if si_t.get((d, virus), -1) >= args.si_threshold}
                test_labels.positives[virus] = keep
        print(f"post-{args.cutoff} labels: {test_labels.stats['active']:,} active, "
              f"{test_labels.stats['inactive']:,} inactive")
        g_train = Graph.load(args.train, skip_predicates=held)
        temporal = run(g_train, schema, test_labels, f"post-{args.cutoff}")
    else:
        print(f"\n{args.train} not found -- run temporal_split.py first")

    args.results.mkdir(parents=True, exist_ok=True)
    (args.results / "hard_negatives.json").write_text(json.dumps(
        {"inactive_above_nm": args.inactive_above_nm,
         "selectivity_filter": args.selectivity,
         "label_stats": dict(all_labels.stats),
         "cross_sectional": cross, "temporal": temporal}, indent=2))
    print(f"\nwrote {args.results / 'hard_negatives.json'}")
    print("\nAUC 0.5 is chance. Unlike Hits@k this is prevalence-independent, so "
          "it is comparable\nacross viruses whose base rates differ by two orders "
          "of magnitude.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
