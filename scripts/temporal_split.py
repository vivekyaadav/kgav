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
    hits_cell,
    lift_cell,
    lift_row,
    p_cell,
    vs_degree_cell,
    undetermined_note,
)
from kgav.calibration import MIN_REACHED_POSITIVES
from kgav.labels import evaluate_against_negatives
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


def _auc(m: dict) -> str:
    """The tie-safe AUC with its interval, or why it is not reported.

    A channel reaching fewer than MIN_REACHED_POSITIVES measured actives has
    an AUC describing the imputed floor rather than the channel, which is
    the same rule docs/evaluation.md applies cross-sectionally. Printed as
    "not evaluable (n)" rather than blank, because an empty cell reads as a
    missing value instead of a refusal.
    """
    if m.get("auc") is None:
        return "—"
    if not m.get("auc_evaluable"):
        return f"not evaluable ({m.get('auc_reached_pos', 0)})"
    return (f"{m['auc']:.3f} [{m['auc_ci_lo']:.2f},{m['auc_ci_hi']:.2f}]")


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

    write_split(split, args.release, args.out, schema_version=s.version,
                policy={"undated": args.undated,
                        "overlap_policy": args.overlap_policy})
    print(f"\nwrote training graph -> {args.out}")
    print(f"  stamped against {args.release.name}; "
          f"policy undated={args.undated}, overlap={args.overlap_policy}")

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
        # KEPT, NOT REBUILT. The AUC below must describe the SAME scorer the
        # lift describes, and a second expression producing "the same" score
        # dict is the two-copies failure this project has hit four times.
        scorers: dict[str, dict[str, float]] = {}
        for name in sorted(dwpc):
            sc = {d: by[virus] for d, by in dwpc[name].items() if virus in by}
            if sc:
                scorers[name] = sc
                per[name] = evaluate(sc, test_pos, args.k)
        scorers["COMBINED"] = combine(dwpc, virus)
        scorers["degree"] = deg
        per["COMBINED"] = evaluate(scorers["COMBINED"], test_pos, args.k)
        per["degree"] = evaluate(scorers["degree"], test_pos, args.k)

        # THE TIE-SAFE METRIC, FINALLY CALLED. Hits@k reads a cut, so where a
        # tie block straddles k the figure is decided by dict insertion
        # order: M1 here returns 0-8 of its 8 positives. The AUC does not
        # read a cut -- labels.evaluate_against_negatives counts tied pairs
        # as half, "which matters because DWPC produces many exact ties at
        # percentile boundaries" -- so it is invariant to that order.
        #
        # This protocol already held both halves and used neither:
        # split.test_negatives is populated (3,630 inactives on v0.1) and
        # only test_pos was ever passed. A tie-sensitive lift was computed
        # against positives alone while the measured negatives sat unused.
        test_neg = split.test_negatives.get(virus, set())
        for name, row in per.items():
            a = evaluate_against_negatives(scorers[name], test_pos, test_neg)
            # The per-channel gate is the SAME constant the virus-level gate
            # above uses, applied to the positives this channel actually
            # reaches rather than to the virus's total. M1 reaches 8 and is
            # refused; COMBINED reaches enough and is not.
            reached = round(a.get("coverage_pos", 0.0) * (a.get("n_pos") or 0))
            row.update({
                "auc": a.get("auc"),
                "auc_ci_lo": a.get("auc_ci_lo"), "auc_ci_hi": a.get("auc_ci_hi"),
                "auc_n_pos": a.get("n_pos"), "auc_n_neg": a.get("n_neg"),
                "auc_reached_pos": reached,
                "auc_coverage_pos": a.get("coverage_pos"),
                "auc_coverage_neg": a.get("coverage_neg"),
                "auc_evaluable": bool(a.get("evaluable")
                                      and reached >= MIN_REACHED_POSITIVES),
            })
        degree_relative(per)
        print(f"    {'scorer':<10} {'pool':>7} {'pos':>5} {'hits':>5} "
              f"{'exp':>7} {'lift':>11} {'p':>7} {'vs deg':>7} {'mrr':>7} "
              f"{'AUC [95% CI]':>22}")
        for name, m in per.items():
            flag = "" if m["trustworthy"] else " (!)"
            # The hits column carries the RANGE where a tie block straddles
            # k, because the single figure there is one draw from it. On
            # v0.1 M1 that range is 2-12 at both k=50 and k=100.
            if not m.get("hits_determined", True):
                flag += " ~tied"
            print(f"    {name:<10} {m['pool']:>7,} {m['pos']:>5,} "
                  f"{hits_cell(m, 5)} "
                  f"{m['expected']:>7.1f} {lift_cell(m)} {p_cell(m)} "
                  f"{vs_degree_cell(m)} {m['mrr']:>7.4f} {_auc(m):>22}{flag}")
        if (warn := confound_warning(per)):
            print(f"    ! {warn}")
        if (note := undetermined_note(per)):
            print(f"    ! {note}")
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
