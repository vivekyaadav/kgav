#!/usr/bin/env python
"""Does any channel carry signal beyond degree.

THE PROTOCOL IS docs/prereg_partial_degree.md, written and committed before
this file. Nothing below departs from it except where marked GAP, and each
GAP names what the registration failed to say and which reading was taken.

Every verdict in docs/evaluation.md compares a channel against a degree
baseline. diagnose_scorer.py measured that both scorings of a channel are
degree with a sign -- legacy correlates +0.3 to +0.7 with the compound's
degree, pdp correlates negatively on 25 of 36 cells, and four cells read
exactly -1.000. A comparison of two quantities each monotone in a third
cannot say which explains the labels.

So this puts both in one model, per (virus, channel, scoring):

    label ~ b_c * channel_percentile + b_d * degree_percentile + intercept

and tests b_c by likelihood ratio against the degree-only model. b_c is the
estimand: what the channel contributes once degree is accounted for.

Registered and not chosen here: lambda, alpha, the floor, the correction,
the feature transform, the test. They are constants below and
tests/test_prereg_partial_degree.py asserts them against the document.

USAGE

  python scripts/partial_degree.py
  python scripts/partial_degree.py --temporal          # the secondary split
  python scripts/partial_degree.py --negatives all     # the sensitivity

Read-only. Writes nothing, recommends nothing, ranks no drug.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import (  # noqa: E402
    SCORING_LEGACY,
    SCORING_PDP,
    Graph,
    channel_percentiles,
    degree_ranking,
    dwpc_scores,
)
from kgav.calibration import MIN_REACHED_POSITIVES  # noqa: E402
from kgav.labels import build_labels  # noqa: E402
from kgav.schema import load_schema  # noqa: E402

PREREG = "docs/prereg_partial_degree.md"

# Registered in that document, before this file existed.
LAMBDA = 1.0
ALPHA = 0.05
MAX_ITER = 100
TOL = 1e-10


# --------------------------------------------------------------- arithmetic
def solve(a: list[list[float]], b: list[float]) -> list[float] | None:
    """Gaussian elimination with partial pivoting. None if singular.

    Three parameters at most, so a hand-rolled solve is shorter than the
    argument for taking on a dependency -- and pyproject declares pyyaml,
    pandas and pyarrow, nothing numeric beyond what pandas pulls in.
    """
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-14:
            return None
        m[col], m[piv] = m[piv], m[col]
        for r in range(n):
            if r == col:
                continue
            f = m[r][col] / m[col][col]
            for c in range(col, n + 1):
                m[r][c] -= f * m[col][c]
    return [m[i][n] / m[i][i] for i in range(n)]


def loglik(x: list[list[float]], y: list[int], beta: list[float]) -> float:
    """The DATA log-likelihood, with no penalty term.

    The penalty is a fitting device; including it in the statistic would
    make the ratio depend on lambda twice.
    """
    total = 0.0
    for row, yi in zip(x, y):
        eta = sum(b * v for b, v in zip(beta, row))
        # log(1 + exp(eta)) computed without overflowing on a separated cell
        total += yi * eta - (eta + math.log1p(math.exp(-eta)) if eta > 0
                             else math.log1p(math.exp(eta)))
    return total


def fit(x: list[list[float]], y: list[int], lam: float = LAMBDA
        ) -> tuple[list[float], bool]:
    """Ridge logistic by IRLS. Column 0 must be the intercept, unpenalised.

    Returns (beta, converged). Separation is expected on small cells: the
    penalty bounds the coefficients, and a cell that still does not converge
    is reported NOT FITTED rather than given a number.
    """
    k = len(x[0])
    beta = [0.0] * k
    pen = [0.0] + [lam] * (k - 1)
    for _ in range(MAX_ITER):
        grad = [0.0] * k
        hess = [[0.0] * k for _ in range(k)]
        for row, yi in zip(x, y):
            eta = sum(b * v for b, v in zip(beta, row))
            p = 1.0 / (1.0 + math.exp(-eta)) if eta > -700 else 0.0
            w = p * (1.0 - p)
            for i in range(k):
                grad[i] += (yi - p) * row[i]
                for j in range(k):
                    hess[i][j] += w * row[i] * row[j]
        for i in range(k):
            grad[i] -= pen[i] * beta[i]
            hess[i][i] += pen[i]
        step = solve(hess, grad)
        if step is None:
            return beta, False
        beta = [b + s for b, s in zip(beta, step)]
        if max(abs(s) for s in step) < TOL:
            return beta, True
    return beta, False


def chi2_sf_1df(stat: float) -> float:
    """P(X > stat) for one degree of freedom, exactly.

    chi2 with 1 df is the square of a standard normal, so the upper tail is
    erfc(sqrt(stat/2)) -- in the standard library, no scipy and no series
    approximation to get subtly wrong.
    """
    if stat <= 0:
        return 1.0
    return math.erfc(math.sqrt(stat / 2.0))


def standardise(values: list[float]) -> list[float] | None:
    """Centre and scale. None when the column is constant, which is not a
    feature: a channel giving every drug the same percentile cannot be
    distinguished from the intercept, and fitting it would return a
    coefficient whose size is set by the penalty alone."""
    n = len(values)
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / n
    if var <= 0:
        return None
    sd = var ** 0.5
    return [(v - mean) / sd for v in values]


# ------------------------------------------------------------------- the cell
class Cell:
    __slots__ = ("virus", "chan", "scoring", "status", "n", "pos",
                 "b_c", "b_d", "stat", "p", "p_adj")

    def __init__(self, virus, chan, scoring, status, n=0, pos=0,
                 b_c=None, b_d=None, stat=None, p=None):
        self.virus, self.chan, self.scoring = virus, chan, scoring
        self.status, self.n, self.pos = status, n, pos
        self.b_c, self.b_d, self.stat, self.p = b_c, b_d, stat, p
        self.p_adj = None


def run_cell(virus, chan, scoring, pool, positives, negatives) -> Cell:
    """One (virus, channel, scoring) cell of the registered measurement."""
    labelled = sorted((set(positives) | set(negatives)) & set(pool))
    pos = [d for d in labelled if d in positives]
    neg = [d for d in labelled if d in negatives]
    if len(pos) < MIN_REACHED_POSITIVES or len(neg) < MIN_REACHED_POSITIVES:
        return Cell(virus, chan, scoring, "NOT RUN", len(labelled), len(pos))
    return Cell(virus, chan, scoring, "fitted", len(labelled), len(pos))


def evaluate(cell: Cell, chan_pct: list[float], deg_pct: list[float],
             y: list[int]) -> Cell:
    """Fit both models and take the likelihood ratio. One df: the full model
    adds exactly the channel."""
    c = standardise(chan_pct)
    d = standardise(deg_pct)
    if c is None or d is None:
        cell.status = "NO SPREAD"
        return cell
    full_x = [[1.0, ci, di] for ci, di in zip(c, d)]
    red_x = [[1.0, di] for di in d]
    bf, ok_f = fit(full_x, y)
    br, ok_r = fit(red_x, y)
    if not (ok_f and ok_r):
        cell.status = "NOT FITTED"
        return cell
    stat = 2.0 * (loglik(full_x, y, bf) - loglik(red_x, y, br))
    cell.b_c, cell.b_d = bf[1], bf[2]
    cell.stat = stat
    cell.p = chi2_sf_1df(max(stat, 0.0))
    return cell


# ----------------------------------------------------------------- reporting
def fmt(v, spec="+.3f"):
    return "    n/a" if v is None else format(v, spec)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path, default=ROOT / "data/releases/v0.1")
    ap.add_argument("--negatives", choices=("measured", "all"),
                    default="measured")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    schema = load_schema()
    g = Graph.load(args.release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    labels = build_labels(args.release)
    degree = degree_ranking(g)

    print(f"protocol: {PREREG}")
    print("Nothing below departs from it except where marked GAP.\n")
    print(f"release {args.release.name} | {len(g.nodes):,} nodes | "
          f"{len(g.drugs()):,} drugs")
    print(f"lambda = {LAMBDA}  alpha = {ALPHA}  "
          f"floor = {MIN_REACHED_POSITIVES} positives AND "
          f"{MIN_REACHED_POSITIVES} negatives")
    print()
    print("GAP: the registration says \"10 negatives\" and never defines a")
    print("negative. Taken here as MEASURED inactive (labels.negatives), the")
    print("set hard_negatives.py already evaluates on, because an unlabelled")
    print("compound is not an inactive one. --negatives all is the")
    print(f"sensitivity. Running: {args.negatives}.\n")

    both = {s: dwpc_scores(g, schema.metapaths, schema.retrieval_limits,
                           scoring=s)
            for s in (SCORING_LEGACY, SCORING_PDP)}

    cells: list[Cell] = []
    for virus in sorted(labels.positives):
        pos = labels.positives.get(virus, set())
        for chan in sorted(schema.metapaths):
            for scoring in (SCORING_LEGACY, SCORING_PDP):
                pool = {d: hits[virus]
                        for d, hits in both[scoring].get(chan, {}).items()
                        if virus in hits}
                if not pool:
                    continue
                neg = (labels.negatives.get(virus, set())
                       if args.negatives == "measured"
                       else set(pool) - set(pos))
                cell = run_cell(virus, chan, scoring, pool, pos, neg)
                if cell.status == "fitted":
                    pct = channel_percentiles(pool)
                    dpct = channel_percentiles(
                        {d: degree.get(d, 0.0) for d in pool})
                    rows = sorted((set(pos) | set(neg)) & set(pool))
                    cell = evaluate(cell, [pct[d] for d in rows],
                                    [dpct[d] for d in rows],
                                    [1 if d in pos else 0 for d in rows])
                cells.append(cell)

    tested = [c for c in cells if c.p is not None]
    for c in tested:                       # Bonferroni over cells ACTUALLY run
        c.p_adj = min(1.0, c.p * len(tested))

    print(f"{len(cells)} cells, {len(tested)} fitted and tested. "
          f"Bonferroni factor {len(tested)}.\n")
    print(f"  {'virus':<18} {'chan':<4} {'scoring':<20} {'n':>5} {'pos':>4} "
          f"{'b_c':>8} {'b_d':>8} {'LR':>8} {'p':>9} {'p_adj':>9}")
    for c in cells:
        if c.p is None:
            print(f"  {c.virus:<18} {c.chan:<4} {c.scoring:<20} "
                  f"{c.n:5d} {c.pos:4d}  {c.status}")
            continue
        print(f"  {c.virus:<18} {c.chan:<4} {c.scoring:<20} {c.n:5d} "
              f"{c.pos:4d} {fmt(c.b_c)} {fmt(c.b_d)} {fmt(c.stat, '8.3f')} "
              f"{fmt(c.p, '9.4f')} {fmt(c.p_adj, '9.4f')}")

    print()
    negative_statistics(tested)
    why_not_tested(cells)
    return verdict(tested)


def negative_statistics(tested: list[Cell]) -> None:
    """A likelihood ratio cannot be negative, and some of these are.

    NOT NOISE, AND NOT A CODING ERROR. The registration says the two models
    are "both penalised identically", meaning the same lambda -- but over
    DIFFERENT numbers of coefficients. The full model maximises
    loglik - lambda/2 * (b_c^2 + b_d^2) and the reduced one
    loglik - lambda/2 * b_d^2, so admitting b_c lets the fit trade b_d away
    against the penalty and land on a LOWER data log-likelihood than the
    reduced model reached. 2 * (ll_full - ll_red) then goes negative.

    Reproduced in 146 of 300 synthetic null cells shaped like the real ones
    -- n about 3,000, 15 positives, channel and degree correlated at 0.95 --
    with the most negative at -0.601.

    THE DIRECTION MATTERS. p is set to 1 for these, so the defect can only
    fail to reject. P1 is a null and its non-rejection survives; what the
    defect costs is power, exactly where channel and degree are collinear,
    which is every cell worth testing.

    The fix belongs in a future registration, not here: take the ratio of
    the PENALISED objectives, which cannot be negative because the full
    model evaluated at b_c = 0 IS the reduced model.
    """
    bad = [c for c in tested if c.stat is not None and c.stat < 0]
    if not bad:
        return
    worst = min(c.stat for c in bad)
    print(f"  {len(bad)} of {len(tested)} cells returned a NEGATIVE "
          f"likelihood ratio (most negative {worst:+.3f}).")
    print("  The two models carry the same lambda over different numbers of")
    print("  coefficients, so the full fit can trade b_d against the penalty")
    print("  and reach a lower DATA log-likelihood than the reduced one.")
    print("  p is set to 1 for those, which is conservative: the defect can")
    print("  only fail to reject, so P1's non-rejection survives it. What it")
    print("  costs is power, where channel and degree are collinear.")
    print()


def why_not_tested(cells: list[Cell]) -> None:
    """The reasons, counted. An untested cell is evidence about the graph --
    NO SPREAD in particular is the M1 degeneracy showing up again -- and a
    verdict that printed only "nothing to test" would lose it."""
    reasons: dict[str, int] = {}
    for c in cells:
        if c.p is None:
            reasons[c.status] = reasons.get(c.status, 0) + 1
    if reasons:
        print("  not tested, by reason:")
        for status, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
            note = {
                "NOT RUN": "below the registered floor",
                "NO SPREAD": "the channel gives every drug one percentile, "
                             "so it cannot be told from the intercept",
                "NOT FITTED": "separated, and the penalty did not bound it",
            }.get(status, "")
            print(f"    {status:<11} {n:4d}  {note}")
        print()


def verdict(tested: list[Cell]) -> int:
    """P1 and P2, read off the cells rather than argued."""
    if not tested:
        # WHY there is nothing to test is not one state. A cell below the
        # floor and a cell whose channel gives every drug the same
        # percentile are different findings, and collapsing them into one
        # sentence would report the M1 degeneracy as a sample-size problem.
        print("NOTHING WAS TESTED. P1 is not testable on this release,")
        print("which is not the same as P1 holding. The registration's")
        print("events-per-variable limit applies a fortiori.")
        return 1

    sig = [c for c in tested if c.p_adj is not None and c.p_adj < ALPHA]
    by_cell: dict[tuple[str, str], list[Cell]] = {}
    for c in tested:
        by_cell.setdefault((c.virus, c.chan), []).append(c)

    consistent = []
    for (virus, chan), group in sorted(by_cell.items()):
        if len(group) != 2:
            continue
        a, b = group
        if (a.p_adj < ALPHA and b.p_adj < ALPHA
                and a.b_c is not None and b.b_c is not None
                and (a.b_c > 0) == (b.b_c > 0)):
            consistent.append((virus, chan, a.b_c, b.b_c))

    print(f"P1  no cell clears corrected alpha under BOTH scorings with the "
          f"same sign")
    if consistent:
        print(f"    REJECTED. {len(consistent)} cell(s) do:")
        for virus, chan, ba, bb in consistent:
            print(f"      {virus} {chan}: b_c {ba:+.3f} / {bb:+.3f}")
    else:
        print(f"    NOT REJECTED. {len(sig)} of {len(tested)} cells are "
              f"significant after correction, none of them under both "
              f"scorings with the same sign.")
        print("    At 2 features and a floor of 10 positives the")
        print("    events-per-variable is 5, so this is WEAK EVIDENCE OF")
        print("    ABSENCE, as the registration states. It is not a finding")
        print("    that the channels are empty.")

    print()
    print("P2  a significant cell will not be sign-consistent across scorings")
    if not sig:
        print("    UNTESTED. No cell is significant after correction, so the")
        print("    falsification test has nothing to run on.")
    elif consistent:
        print(f"    REJECTED for {len(consistent)} cell(s), which is the one")
        print("    outcome that would justify the full regularized model.")
    else:
        print(f"    HELD. {len(sig)} significant cell(s), none sign-consistent")
        print("    across the two scorings -- the signature of degree wearing")
        print("    a sign rather than a channel carrying mechanism.")
    return 0 if consistent else 1


if __name__ == "__main__":
    raise SystemExit(main())
