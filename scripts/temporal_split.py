#!/usr/bin/env python3
"""Build the temporal split, run the leakage audits, and evaluate on it."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.baselines import (
    Graph,
    combine,
    confound_warning,
    degree_ranking,
    degree_relative,
    dwpc_scores,
    lift_row,
)
from kgav.calibration import MIN_REACHED_POSITIVES
from kgav.provenance import write_refusal, write_results
from kgav.schema import load_schema
from kgav.temporal import audit, build_split, write_split

# THIS SCRIPT USED TO CARRY ITS OWN COPY OF lift_row, and the copy was the
# version from before the pool-size bug was fixed:
#
#     "trustworthy": exp >= MIN_TRUSTWORTHY_EXPECTATION      <- local copy
#     "trustworthy": exp >= MIN_TRUSTWORTHY_EXPECTATION and pool > k
#
# baselines.lift_row's own docstring records why the second condition exists:
# when the pool is no larger than k, ranked[:k] IS the whole pool, so lift is
# exactly pool/k and carries no information. The fix landed in baselines and
# never reached this copy, so the temporal output printed SARS-CoV-2 M7 0.91
# (pool 91), SARS-CoV M7 0.34 (pool 34) and MERS-CoV M7 0.14 (pool 14)
# unflagged, while compare_metapaths.py flagged the same channel on the same
# release. One concept, two implementations, and only one of them was fixed.
evaluate = lift_row

DEFAULT_OVERLAP_POLICY = "drop-from-test"


def _vs(m: dict) -> str:
    """Lift relative to the degree baseline, which is the comparison that
    carries information; p alone is against a random ranking."""
    v = m.get("lift_vs_degree")
    if v is None:
        return "      -"
    return f"{v:>7.2f}"


def _p(m: dict) -> str:
    """p shown only where it answers something: None means the question was
    empty (no positives, or k covering the whole pool)."""
    v = m.get("p_value")
    return "      -" if v is None else (f"{v:>7.4f}" if v >= 0.0001
                                        else "< .0001")


def _results_name(args) -> str:
    """The results filename, carrying any NON-DEFAULT policy.

    WHY. --overlap-policy keep exists to reproduce the refusal, and it wrote
    to the same path as the real run -- so running the control destroyed the
    result it was a control for, replacing a full temporal evaluation with a
    refusal record. Same shape as data/results vs data/results-corrected:
    two different runs, one filename.

    The default keeps its existing name, so nothing already on disk is
    renamed and no reader has to learn a new path.
    """
    name = f"temporal_{args.cutoff}_{args.undated}"
    if args.overlap_policy != DEFAULT_OVERLAP_POLICY:
        name += f"_overlap-{args.overlap_policy}"
    return name + ".json"


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
    ap.add_argument("--overlap-policy", default=DEFAULT_OVERLAP_POLICY,
                    choices=["drop-from-test", "keep"],
                    help="a compound measured on BOTH sides of the cutoff is "
                         "already known, so by default it leaves the test "
                         "set; 'keep' reproduces the refusal")
    args = ap.parse_args()

    s = load_schema()
    held = s.held_out_predicates()
    print(f"cutoff {args.cutoff} | undated policy: {args.undated}")

    # held = what to strip from the training graph; labels = what to partition
    # into positives and negatives. Not the same set since schema 0.11.0.
    split = build_split(args.release, args.cutoff, held, args.undated,
                        label_predicates=s.evaluation_label_predicates(),
                        overlap_policy=args.overlap_policy)
    for v in split.violations:
        print(f"  already known: {v}")
    print(f"\ntraining graph: {len(split.train_edges):,} edges")
    for key in ("edge_train", "edge_dropped_future", "undated_kept",
                "undated_kept_computed", "undated_dropped",
                "held_out_not_a_label",
                "label_train_active", "label_train_inactive",
                "label_test_active", "label_test_inactive",
                "label_undecidable", "label_undated",
                "train_ambiguous_dropped", "test_ambiguous_dropped",
                "test_positive_already_known_dropped",
                "test_negative_already_known_dropped"):
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
        # RECORD THE REFUSAL WHERE THE RESULTS WOULD HAVE GONE. Returning
        # without writing left the last successful run's numbers on disk, and
        # those numbers came from a graph that no longer exists. check_results
        # called them STALE and said "regenerate it" -- advice that cannot be
        # followed, because this branch is why they cannot be regenerated.
        args.results.mkdir(parents=True, exist_ok=True)
        out = args.results / _results_name(args)
        write_refusal(out, {"cutoff": args.cutoff,
                            "undated_policy": args.undated,
                            "overlap_policy": args.overlap_policy,
                            "stats": dict(split.stats), "audits": problems},
                      fatal, args.release, s.version)
        print(f"recorded the refusal -> {out}")
        print("  (the previous run's numbers are gone: they described an "
              "older graph)")
        return 1

    write_split(split, args.release, args.out)
    print(f"\nwrote training graph -> {args.out}")

    g = Graph.load(args.out, skip_predicates=held, symmetry=s.symmetry())
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
              f"({len(reachable) / len(test_pos):.1%} ceiling on recall)")
        # THE CEILING HAS TO GATE, NOT JUST ANNOTATE. It was a `note` in the
        # audit and a printed line here, and the only refusal was at exactly
        # zero -- so HCoV-NL63 was evaluated on 2 reachable compounds out of 5
        # and its lift reported alongside SARS-CoV-2's. Every metric on a test
        # set that is 85% unreachable by construction is dominated by the
        # unreachable part, which is the same failure the cross-sectional side
        # gates with MIN_REACHED_POSITIVES. One constant, imported, so the two
        # protocols cannot drift to different floors.
        if len(reachable) < MIN_REACHED_POSITIVES:
            print(f"    NOT EVALUABLE — {len(reachable)} reachable, below the "
                  f"{MIN_REACHED_POSITIVES} needed for a metric to describe "
                  f"the scorer rather than the ceiling")
            results[label] = {"not_evaluable": {
                "reachable": len(reachable), "test_positives": len(test_pos),
                "floor": MIN_REACHED_POSITIVES}}
            continue
        per: dict[str, dict] = {}
        for name in sorted(dwpc):
            sc = {d: by[virus] for d, by in dwpc[name].items() if virus in by}
            if sc:
                per[name] = evaluate(sc, test_pos, args.k)
        per["COMBINED"] = evaluate(combine(dwpc, virus), test_pos, args.k)
        per["degree"] = evaluate(deg, test_pos, args.k)
        degree_relative(per)
        print(f"    {'scorer':<10} {'pool':>7} {'pos':>5} {'hits':>5} "
              f"{'exp':>7} {'lift':>7} {'p':>7} {'vs deg':>7} {'mrr':>7}")
        for name, m in per.items():
            flag = "" if m["trustworthy"] else " (!)"
            print(f"    {name:<10} {m['pool']:>7,} {m['pos']:>5,} {m['hits']:>5} "
                  f"{m['expected']:>7.1f} {m['lift']:>7.2f} {_p(m)} "
                  f"{_vs(m)} {m['mrr']:>7.4f}{flag}")
        if (warn := confound_warning(per)):
            print(f"    ! {warn}")
        results[label] = per

    args.results.mkdir(parents=True, exist_ok=True)
    write_results(args.results / _results_name(args),
                  {"cutoff": args.cutoff, "undated_policy": args.undated,
                   "overlap_policy": args.overlap_policy,
                   "reachability_floor": MIN_REACHED_POSITIVES,
                   "stats": dict(split.stats), "audits": problems,
                   "already_known": split.violations,
                   "results": results}, args.release, s.version)
    print(f"\nwrote {args.results}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
