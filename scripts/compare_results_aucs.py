#!/usr/bin/env python3
"""Every AUC in data/results, this release against an earlier commit.

WHY THIS EXISTS. A rebuild moves a thousand figures, and the only question
that decides whether the rebuild was a repair or a regression is whether the
scores went down. That question was answered by hand, once, by reading commit
messages -- and a repair that quietly cost 0.02 of AUC on one cell would have
been invisible beside 106 failing prose assertions.

WHAT IT COMPARES. Every numeric leaf whose key is `auc` or begins with `auc_`,
in every results file, keyed by (file, path) so a cell is matched to itself
rather than to a cell of the same name elsewhere. `auc_n_pos` and `auc_n_neg`
are sample sizes rather than scores and are excluded; `auc_reached_pos` is a
count and is reported separately for the same reason.

WHAT IT REFUSES TO DO. It does not say whether a drop is acceptable. A drop
can be the correct consequence of a repair -- deduplicating a salt onto its
parent removes a compound that was being scored twice, and the honest AUC is
the lower one -- so the script reports the move and names the cell, and the
reading is a human's. What it will not let happen is the move going unseen.

A VANISHED CELL IS THE FINDING THIS WAS WRITTEN FOR. The L27 identity repair
left every AUC within 0.003 and took FOUR (run, virus) cells below the
reachability floor: MERS-CoV and SARS-CoV under temporal `computed_only` and
`exclude`, each at 9 reachable test positives against a floor of 10. No
figure moved to say so. All four read CONFOUNDED with no readable significant
channel, so nothing claimable was lost -- but that was established after the
fact, by reading the superseded files out of git, which is what this does by
command.

    python scripts/compare_results_aucs.py --against 5151de8
    make auc-delta ARGS="--against HEAD~2 --fail-below -0.01"
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
RESULTS = ROOT / "data" / "results"

from kgav.calibration import MIN_REACHED_POSITIVES  # noqa: E402

#: Leaves named like an AUC that are not one. Sample sizes and reach counts
#: live under the same prefix because they describe the AUC's population.
NOT_A_SCORE = ("auc_n_pos", "auc_n_neg", "auc_reached_pos",
               "auc_reached_neg")

#: A RATIO IS NOT A SCORE, AND THE FIRST VERSION OF THIS SCRIPT COMPARED THEM
#: TOGETHER. `auc_skill` is (auc - floor) / (ceiling - floor), so on a cell
#: whose span is a few thousandths it swings by most of its range for an AUC
#: move of 0.001: the largest "fall" in the first run was -0.9286 on
#: MERS-CoV M5's skill, where the AUC itself moved 0.0048. Gating on the
#: mixture would have refused every rebuild forever while telling a reader
#: nothing, which is the shape of a check that cannot pass. Skill is reported
#: in its own section and never gated.
RATIOS = ("auc_skill",)


def leaves(obj, path=()):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from leaves(v, path + (k,))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from leaves(v, path + (str(i),))
    else:
        yield path, obj


def _cells(obj, path=()):
    """Every dict that holds an `auc`, with its path. The scorer's own
    readability flags are siblings of the score, so the cell is the unit."""
    if isinstance(obj, dict):
        if "auc" in obj and isinstance(obj.get("auc"), (int, float)):
            yield path, obj
        for k, v in obj.items():
            yield from _cells(v, path + (k,))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _cells(v, path + (str(i),))


def aucs(directory: Path) -> tuple[dict[tuple[str, ...], float],
                                   dict[tuple[str, ...], float],
                                   set[tuple[str, ...]],
                                   set[tuple[str, ...]]]:
    """(scores, ratios, reported, claiming), keyed by (file, *json path).

    `scores` are AUCs and the quantities on the same 0..1 scale: the
    corrected interval bounds, the imputation floor and the ceiling. `ratios`
    are `auc_skill`, kept apart for the reason on RATIOS above.

    `reported` is the set of CELL paths this project is willing to interpret.
    THAT DISTINCTION IS THE WHOLE POINT OF THE GATE. Every material AUC fall
    the L27 repair produced is in a cell whose own flags refuse it --
    MERS-CoV M2 under temporal `include` fell 0.5300 -> 0.4288 on SIX
    negatives, where reaching one more negative moves the AUC by a tenth, and
    the file already says `trustworthy: false`. Gating on all cells would
    refuse a repair over numbers the project does not publish; gating on none
    would miss a real regression. So the gate reads the scorer's own verdict
    on whether the number means anything.
    """
    scores: dict[tuple[str, ...], float] = {}
    ratios: dict[tuple[str, ...], float] = {}
    reported: set[tuple[str, ...]] = set()
    claiming: set[tuple[str, ...]] = set()
    for f in sorted(directory.glob("*.json")):
        try:
            d = json.loads(f.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        for p, v in leaves(d):
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                continue
            key = p[-1]
            if key in NOT_A_SCORE:
                continue
            if key in RATIOS:
                ratios[(f.name,) + p] = float(v)
            elif key == "auc" or key.startswith("auc_"):
                scores[(f.name,) + p] = float(v)
        for cp, cell in _cells(d):
            if cell.get("trustworthy") is False:
                continue
            if cell.get("hits_determined") is False:
                continue
            if cell.get("readable") is False:
                continue
            if cell.get("measurable") is False:
                continue
            # AN ABSENT FLAG IS NOT A PASS, and the first version treated it
            # as one. `hard_negatives_row1.json` carries no `trustworthy`, so
            # every cell in it was admitted -- including MERS-CoV COMBINED at
            # 6.4% coverage of 47 positives, THREE reached compounds, whose
            # AUC moved 0.0438 and triggered the gate. An AUC over three
            # reached positives is what MIN_REACHED_POSITIVES exists to
            # refuse on both protocols, so the same constant decides here,
            # imported rather than re-typed.
            cov, npos = cell.get("coverage_pos"), cell.get("n_pos")
            if cov is None:
                cov, npos = cell.get("auc_coverage_pos"), cell.get("auc_n_pos")
            if cov is not None and npos:
                if round(cov * npos) < MIN_REACHED_POSITIVES:
                    continue
            reported.add((f.name,) + cp)
            if _claims(cell):
                claiming.add((f.name,) + cp)
    return scores, ratios, reported, claiming


def _claims(cell: dict) -> bool:
    """Whether this cell asserts a result: an interval that excludes chance.

    Two spellings in this project, both read here. The structure scorer
    records `clears_corrected_alpha` directly; the hard-negative and temporal
    protocols record the interval and the claim is that it excludes 0.5. A
    cell that asserts nothing in either release cannot be a result that got
    weaker, which is what makes this the unit the gate should judge.
    """
    if cell.get("clears_corrected_alpha") is True:
        return True
    for lo_k, hi_k in (("auc_ci_lo_corrected", "auc_ci_hi_corrected"),
                       ("auc_ci_lo", "auc_ci_hi")):
        lo, hi = cell.get(lo_k), cell.get(hi_k)
        if lo is not None and hi is not None:
            if lo > 0.5 or hi < 0.5:
                return True
    return False


def at_commit(commit: str, into: Path) -> Path:
    """data/results as of `commit`, extracted into a scratch directory."""
    tar = subprocess.run(["git", "archive", commit, "data/results"],
                         cwd=ROOT, capture_output=True)
    if tar.returncode != 0:
        raise SystemExit(f"git archive {commit} failed: "
                         f"{tar.stderr.decode()[:200]}")
    subprocess.run(["tar", "-x", "-C", str(into)], input=tar.stdout,
                   check=True)
    return into / "data" / "results"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--against", default="HEAD~1",
                    help="the commit whose data/results to compare with")
    # TWO DIRECTORIES, NO GIT. For comparing a candidate regeneration against
    # the shipped one before committing either -- and it is what makes the
    # gate testable, which matters because a gate that has only ever passed
    # is indistinguishable from one that cannot fail.
    ap.add_argument("--against-dir", type=Path, default=None,
                    help="compare with this results directory instead of a "
                         "commit")
    ap.add_argument("--new-dir", type=Path, default=None,
                    help="treat this directory as the new side instead of "
                         "data/results")
    ap.add_argument("--tol", type=float, default=5e-5,
                    help="moves smaller than this read as unchanged")
    ap.add_argument("--top", type=int, default=10,
                    help="how many of the largest moves to name per direction")
    ap.add_argument("--fail-below", type=float, default=None,
                    help="exit 1 if any cell fell by more than this (a "
                         "negative number, e.g. -0.01)")
    args = ap.parse_args(argv)

    if args.against_dir is not None:
        old, old_ratio, old_rep, old_claim = aucs(args.against_dir)
        source = str(args.against_dir)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            old, old_ratio, old_rep, old_claim = aucs(
                at_commit(args.against, Path(tmp)))
        source = args.against
    new_side = args.new_dir or RESULTS
    new, new_ratio, new_rep, new_claim = aucs(new_side)

    shared = sorted(set(old) & set(new))
    print(f"{source} -> {new_side}")
    print(f"  {len(old):,} score leaves before, {len(new):,} after, "
          f"{len(shared):,} comparable (AUC, interval bounds, floor, ceiling)")

    up = [k for k in shared if new[k] - old[k] > args.tol]
    dn = [k for k in shared if old[k] - new[k] > args.tol]
    print(f"  rose {len(up):,}   fell {len(dn):,}   "
          f"unchanged {len(shared) - len(up) - len(dn):,}")

    for keys, label in ((dn, "FELL"), (up, "ROSE")):
        if not keys:
            continue
        print(f"\n  {label}, largest {min(args.top, len(keys))}:")
        for k in sorted(keys, key=lambda k: -abs(new[k] - old[k]))[:args.top]:
            print(f"    {new[k] - old[k]:+.4f}  {old[k]:.4f} -> {new[k]:.4f}"
                  f"   {'/'.join(k)}")

    # A CELL THAT STOPPED EXISTING, which no delta can show. Grouped, because
    # one refused (run, virus) removes every channel under it at once and a
    # flat list of 48 paths reads as 48 losses.
    gone = sorted(set(old) - set(new))
    if gone:
        print(f"\n  {len(gone):,} leaf/leaves EXIST BEFORE AND NOT AFTER "
              f"-- a measurement the release no longer supports:")
        for group, n in Counter(k[:3] for k in gone).most_common():
            print(f"    {n:>4}  {' / '.join(group)}")
    fresh = sorted(set(new) - set(old))
    if fresh:
        print(f"\n  {len(fresh):,} leaf/leaves are new:")
        for k in fresh[:args.top]:
            print(f"    {new[k]:.4f}  {'/'.join(k)}")

    # SKILL, REPORTED AND NOT GATED. Beside each move, the AUC move of the
    # same cell -- which is the whole reason skill is not a score: a reader
    # seeing -0.93 needs the 0.0048 next to it.
    r_shared = sorted(set(old_ratio) & set(new_ratio))
    r_moved = [k for k in r_shared
               if abs(new_ratio[k] - old_ratio[k]) > args.tol]
    print(f"\n  auc_skill (a RATIO over each cell's own span, not a score): "
          f"{len(r_moved):,} of {len(r_shared):,} moved")
    for k in sorted(r_moved, key=lambda k: -abs(new_ratio[k] - old_ratio[k])
                    )[:args.top]:
        auc_key = k[:-1] + ("auc",)
        auc_note = (f"   [its AUC {new[auc_key] - old[auc_key]:+.4f}]"
                    if auc_key in old and auc_key in new else "")
        print(f"    {new_ratio[k] - old_ratio[k]:+.4f}  "
              f"{old_ratio[k]:.4f} -> {new_ratio[k]:.4f}   "
              f"{'/'.join(k)}{auc_note}")

    # THE AUC ITSELF, IN CELLS BOTH RELEASES REPORT. The interval bounds,
    # floor and ceiling move with coverage and are reported above; what a
    # verdict rests on is `auc`, and only where the scorer says the number
    # means something.
    judged = [k for k in shared
              if k[-1] == "auc" and k[:-1] in old_rep and k[:-1] in new_rep]
    j_up = sum(1 for k in judged if new[k] - old[k] > args.tol)
    j_dn = [k for k in judged if old[k] - new[k] > args.tol]
    print(f"\n  AUC in cells BOTH releases report: {len(judged):,} cells, "
          f"{j_up:,} rose, {len(j_dn):,} fell")
    if j_dn:
        print("    largest falls:")
        for k in sorted(j_dn, key=lambda k: new[k] - old[k])[:args.top]:
            print(f"    {new[k] - old[k]:+.4f}  {old[k]:.4f} -> {new[k]:.4f}"
                  f"   {'/'.join(k)}")

    # A COMPARATOR IS NOT A RESULT, and there are exactly two of them. The
    # project's claims are made on a controlled scorer and are compared
    # against (1) `degree`, the null hypothesis, and (2) the same scorer with
    # its control removed -- `structure_max_no_pub_exclusion`, published as
    # the leak magnitude. Both falling means the argument got STRONGER: a
    # weaker null and a smaller leak. Three of the seven cells the first
    # version of this gate refused were degree rows and the last was the
    # uncontrolled variant, so without this category the gate reports a
    # repair as a regression. They are counted on their own line and never
    # gated. The category is stated once rather than grown per failing cell,
    # which is the difference between a definition and a tuned floor.
    # `in` ON A TUPLE IS MEMBERSHIP, NOT SUBSTRING, and the first version of
    # this line read `"no_pub_exclusion" in k`. The path element is
    # `structure_max_no_pub_exclusion`, so that half never matched and the
    # gate went on refusing the one cell this category was written for --
    # while the `degree` half worked, because that element is exactly
    # "degree". A half-matching filter is the worst kind: it looks applied.
    base = [k for k in judged
            if any(part == "degree" or "no_pub_exclusion" in part
                   for part in k)]
    if base:
        d = [new[k] - old[k] for k in base]
        print(f"\n  COMPARATORS (degree, and the control-free variant -- "
              f"not results): {len(base)} cells, "
              f"{sum(1 for x in d if x < -args.tol)} fell, "
              f"{sum(1 for x in d if x > args.tol)} rose, "
              f"largest fall {min(d):+.4f}. A comparator falling is the "
              f"argument strengthening.")

    # WHAT A VERDICT RESTS ON: a cell whose interval excludes chance in
    # EITHER release. A cell asserting nothing in both cannot be a result
    # that weakened; a cell that STOPPED asserting is caught whatever its
    # magnitude, which is the failure a floor alone would miss.
    claim = [k for k in judged if k not in set(base)
             and (k[:-1] in old_claim or k[:-1] in new_claim)]
    lost = [k for k in claim if k[:-1] in old_claim and k[:-1] not in new_claim]
    gained = [k for k in claim if k[:-1] not in old_claim and k[:-1] in new_claim]
    print(f"\n  cells carrying a claim (interval excludes chance in either "
          f"release): {len(claim)}")
    print(f"    stopped carrying it: {len(lost)}    started: {len(gained)}")
    for k in lost[:args.top]:
        print(f"    LOST  {old[k]:.4f} -> {new[k]:.4f}   {'/'.join(k)}")
    for k in gained[:args.top]:
        print(f"    NEW   {old[k]:.4f} -> {new[k]:.4f}   {'/'.join(k)}")

    if args.fail_below is not None:
        worst = min((new[k] - old[k] for k in claim), default=0.0)
        bad = [k for k in claim if new[k] - old[k] < args.fail_below]
        print(f"\n  largest fall in a claim-carrying AUC {worst:+.4f} "
              f"against a floor of {args.fail_below:+.4f}")
        if bad or lost:
            if bad:
                print(f"  REFUSED: {len(bad)} claim-carrying AUC(s) fell "
                      f"further than the floor allows:")
                for k in sorted(bad, key=lambda k: new[k] - old[k])[:args.top]:
                    print(f"    {new[k] - old[k]:+.4f}  {'/'.join(k)}")
            if lost:
                print(f"  REFUSED: {len(lost)} cell(s) stopped carrying a "
                      f"claim, which no floor can excuse.")
            return 1
        print("  PASS: no claim-carrying AUC fell further than the floor "
              "allows, and none stopped carrying its claim.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
