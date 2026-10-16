#!/usr/bin/env python
"""How many more labels would it take, and are there compounds to label.

docs/prereg_partial_degree.md returned NO RESULT and its verdict names one
actionable cell -- SARS-CoV-2 M4, b_c +0.401 / +0.403 under the two
scorings, significant uncorrected and not after Bonferroni -- and says the
action is more labels rather than more modelling. That is a quantity, so
this measures it instead of arguing about it.

Two numbers per cell:

  HEADROOM     compounds the channel reaches that carry no label at all.
               An upper bound on what any second screening panel could add
               WITHOUT extending the graph: a panel can only help where a
               metapath already arrives.

  LABELS NEEDED  the labelled-set size at which the registered test would
               reach 80% power, by simulation from that cell's own fitted
               coefficients and its own feature correlation.

THE SIMULATION USES THE REGISTERED STATISTIC, imported from
partial_degree rather than reimplemented -- including its defect, that the
penalised ratio can come out negative because the two models carry the
same lambda over different numbers of coefficients. A power estimate from
a cleaner statistic would promise a sensitivity the registered test does
not deliver, which is the same error the Phase 4 power table made.

The curve itself is analytic, because simulating it cost hours; the
simulation is kept beside it as a check. Measured over six settings the
two agree within 0.071, and the SIMULATION IS HIGHER in five of six -- the
closed form's n*p(1-p) information term is conservative, so 'needed' errs
toward asking for more labels rather than fewer. At b_c = 0 both sit at
alpha, which is the calibration that matters most here.

Features are simulated as MIDRANK PERCENTILES, not as Gaussians: the real
pipeline percentiles both columns, and a normal-features power curve would
be answering about a different model.

NOT A REGISTRATION. This plans a measurement; it does not perform one and
its output is not evidence about any channel. If a panel is acquired, the
analysis of it needs its own prereg.

USAGE

  python scripts/label_headroom.py
  python scripts/label_headroom.py --virus NCBITaxon:2697049 --channel M4

Read-only. Writes nothing.
"""
from __future__ import annotations

import argparse
import math
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.baselines import (  # noqa: E402
    SCORING_LEGACY,
    SCORING_PDP,
    Graph,
    channel_percentiles,
    degree_ranking,
    dwpc_scores,
)
from kgav.labels import build_labels  # noqa: E402
from kgav.schema import load_schema  # noqa: E402

# THE REGISTERED ARITHMETIC, imported. A second copy would drift, and the
# whole point is to estimate the power of the test that actually ran.
from diagnose_scorer import ranks, spearman  # noqa: E402
from partial_degree import (  # noqa: E402
    ALPHA,
    LAMBDA,
    chi2_sf_1df,
    fit,
    loglik,
    standardise,
)

TARGET_POWER = 0.80
TRIALS = 200

# WHY THERE IS A CLOSED FORM AT ALL. The first version found the required
# label count by simulation alone, doubling n and running 200 fits at each
# step. One point at n = 3,000 costs 19 seconds in pure Python and the
# search reached 64,000, so a single cell ran for hours. The curve is now
# analytic and the simulation is kept as its CHECK -- the better
# arrangement anyway, because the closed form cannot see the registered
# statistic's negativity defect, and comparing the two is how that cost
# gets measured rather than assumed.


def normal_cdf(x: float) -> float:
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def chi2_crit_1df(alpha: float) -> float:
    """The c with P(chi2_1 > c) = alpha, by bisection on the exact tail
    already in partial_degree. No table and no approximation."""
    lo, hi = 0.0, 1000.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if chi2_sf_1df(mid) > alpha:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def power_closed_form(n: int, pos_rate: float, b_c: float, rho: float,
                      correction: int) -> float:
    """Power of the one-df test on b_c, adjusted for a covariate at rho.

    On standardised features the information about b_c after adjusting for
    degree is n * p(1-p) * (1 - rho^2), so the non-centrality is b_c^2
    times that; a one-df chi-square with non-centrality lam is
    (Z + sqrt(lam))^2.
    """
    alpha_c = ALPHA / max(correction, 1)
    crit = chi2_crit_1df(alpha_c)
    lam = (b_c ** 2) * n * pos_rate * (1.0 - pos_rate) * (1.0 - rho * rho)
    root = math.sqrt(max(lam, 0.0))
    s = math.sqrt(crit)
    return normal_cdf(root - s) + normal_cdf(-root - s)


def headroom(pool: set[str], labelled: set[str]) -> tuple[int, int, int]:
    """(reached, labelled within reach, unlabelled within reach)."""
    inside = pool & labelled
    return len(pool), len(inside), len(pool) - len(inside)


def simulate_once(n: int, pos_rate: float, b_c: float, b_d: float,
                  rho: float, rng: random.Random) -> float | None:
    """One synthetic cell, scored by the registered statistic.

    Percentiles, not Gaussians: correlated normals are drawn, then replaced
    by their midranks, which is what channel_percentiles produces.
    """
    z = [rng.gauss(0, 1) for _ in range(n)]
    s = max(0.0, 1.0 - rho * rho) ** 0.5
    w = [rho * zi + s * rng.gauss(0, 1) for zi in z]
    deg = [r / n for r in ranks(z)]
    chan = [r / n for r in ranks(w)]
    c, d = standardise(chan), standardise(deg)
    if c is None or d is None:
        return None
    base = math.log(pos_rate / (1 - pos_rate))
    y = [1 if rng.random() < 1 / (1 + math.exp(-(base + b_c * ci + b_d * di)))
         else 0 for ci, di in zip(c, d)]
    if len(set(y)) < 2:
        return None
    full = [[1.0, ci, di] for ci, di in zip(c, d)]
    red = [[1.0, di] for di in d]
    bf, okf = fit(full, y, LAMBDA)
    br, okr = fit(red, y, LAMBDA)
    if not (okf and okr):
        return None
    stat = 2.0 * (loglik(full, y, bf) - loglik(red, y, br))
    return chi2_sf_1df(max(stat, 0.0))


def power_at(n: int, pos_rate: float, b_c: float, b_d: float, rho: float,
             correction: int, trials: int = TRIALS, seed: int = 0) -> float:
    """Fraction of synthetic cells that clear the CORRECTED alpha.

    The correction is passed in rather than assumed: power against an
    uncorrected alpha would overstate what the registered analysis, which
    Bonferroni-corrects across every cell it runs, can actually detect.
    """
    rng = random.Random(seed)
    hits = ran = 0
    for _ in range(trials):
        p = simulate_once(n, pos_rate, b_c, b_d, rho, rng)
        if p is None:
            continue
        ran += 1
        if min(1.0, p * correction) < ALPHA:
            hits += 1
    return hits / ran if ran else 0.0


def labels_needed(pos_rate: float, b_c: float, rho: float, correction: int,
                  cap: int = 1_000_000) -> tuple[int | None, float]:
    """Smallest labelled-set size reaching TARGET_POWER, from the closed
    form. (None, power_at_cap) when the cap is reached, which is itself the
    answer: the effect is too small to chase with labels."""
    at_cap = power_closed_form(cap, pos_rate, b_c, rho, correction)
    if at_cap < TARGET_POWER:
        return None, at_cap
    lo, hi = 1, cap
    while lo < hi:
        mid = (lo + hi) // 2
        if power_closed_form(mid, pos_rate, b_c, rho, correction) \
                >= TARGET_POWER:
            hi = mid
        else:
            lo = mid + 1
    return hi, power_closed_form(hi, pos_rate, b_c, rho, correction)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path, default=ROOT / "data/releases/v0.1")
    ap.add_argument("--virus", default=None)
    ap.add_argument("--channel", default=None)
    ap.add_argument("--trials", type=int, default=TRIALS)
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    schema = load_schema()
    g = Graph.load(args.release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    labels = build_labels(args.release)
    degree = degree_ranking(g)
    both = {s: dwpc_scores(g, schema.metapaths, schema.retrieval_limits,
                           scoring=s)
            for s in (SCORING_LEGACY, SCORING_PDP)}

    print("PLANNING ONLY. This sizes a future measurement; it is not one,")
    print("and nothing here is evidence about any channel. A panel, if")
    print("acquired, needs its own registration.\n")
    print(f"release {args.release.name} | power target {TARGET_POWER:.0%} | "
          f"{args.trials} trials per point")
    print(f"statistic: the registered one, imported from partial_degree, "
          f"negativity defect included\n")

    viruses = [args.virus] if args.virus else sorted(labels.positives)
    rows = []
    for virus in viruses:
        pos = labels.positives.get(virus, set())
        neg = labels.negatives.get(virus, set())
        chans = [args.channel] if args.channel else sorted(schema.metapaths)
        for chan in chans:
            for scoring in (SCORING_LEGACY, SCORING_PDP):
                raw = {d: hits[virus]
                       for d, hits in both[scoring].get(chan, {}).items()
                       if virus in hits}
                if not raw:
                    continue
                reached, lab, free = headroom(set(raw), pos | neg)
                rows.append((virus, chan, scoring, reached, lab, free,
                             raw, pos, neg, degree))

    print(f"  {'virus':<18} {'chan':<4} {'scoring':<20} {'reached':>8} "
          f"{'labelled':>9} {'headroom':>9}")
    seen = set()
    for virus, chan, scoring, reached, lab, free, *_ in rows:
        if (virus, chan) in seen:
            continue
        seen.add((virus, chan))
        print(f"  {virus:<18} {chan:<4} {'(reach is scoring-independent)':<20} "
              f"{reached:8,} {lab:9,} {free:9,}")
    print()

    fitted: list[tuple[str, str, int, int | None]] = []
    print("  power, for cells the registered test could fit")
    print(f"  {'virus':<18} {'chan':<4} {'scoring':<20} {'b_c':>7} "
          f"{'rho':>7} {'rate':>7} {'now':>7} {'sim':>7} {'needed':>9}")
    correction = sum(1 for *_x, raw, pos, neg, _d in rows
                     if len((pos | neg) & set(raw)) >= 20)
    correction = max(correction, 1)
    for virus, chan, scoring, reached, lab, free, raw, pos, neg, deg in rows:
        inside = sorted((pos | neg) & set(raw))
        n_pos = sum(1 for d in inside if d in pos)
        if n_pos < 10 or len(inside) - n_pos < 10:
            continue
        pct = channel_percentiles(raw)
        dpct = channel_percentiles({d: deg.get(d, 0.0) for d in raw})
        c = standardise([pct[d] for d in inside])
        dd = standardise([dpct[d] for d in inside])
        if c is None or dd is None:
            continue
        y = [1 if d in pos else 0 for d in inside]
        bf, ok = fit([[1.0, ci, di] for ci, di in zip(c, dd)], y, LAMBDA)
        if not ok:
            continue
        rho = spearman([pct[d] for d in inside], [dpct[d] for d in inside])
        rho = 0.0 if rho is None else rho
        rate = n_pos / len(inside)
        now = power_closed_form(len(inside), rate, bf[1], rho, correction)
        sim = power_at(len(inside), rate, bf[1], bf[2], rho, correction,
                       args.trials)
        need, pw = labels_needed(rate, bf[1], rho, correction)
        need_s = f"{need:,}" if need else f">1M ({pw:.0%})"
        fitted.append((virus, chan, reached, need))
        print(f"  {virus:<18} {chan:<4} {scoring:<20} {bf[1]:+7.3f} "
              f"{rho:+7.3f} {rate:7.1%} {now:7.0%} {sim:7.0%} {need_s:>9}")

    print()
    print("  'now' is the closed-form power at the current labelled size and")
    print("  'sim' the same by simulation under the REGISTERED statistic.")
    print("  Checked against each other over six settings: agreement within")
    print("  0.071, with the simulation HIGHER in five of six -- the closed")
    print("  form's n*p(1-p) information term is conservative. At b_c = 0")
    print("  both sit at alpha. Where the two diverge, 'sim' is the one that")
    print("  used the test that actually ran.")
    print()
    print("  'needed' is the LABELLED-set size for 80% power at the")
    print("  corrected alpha, from that cell's own fitted coefficients. It")
    print("  assumes a panel labels compounds the channel already reaches,")
    print("  at the same positive rate, and that the effect is real -- an")
    print("  effect estimated from a non-significant cell is the weakest")
    print("  input here, so read the column as an order of magnitude.")
    print()
    ceiling_verdict(fitted)
    return 0


def ceiling_verdict(fitted: list[tuple[str, str, int, int | None]]) -> None:
    """Compare what each cell needs against what its channel can supply.

    THE SCRIPT USED TO ASK THE READER TO DO THIS. It printed reach in one
    table and the required label count in another, then a paragraph saying
    "where needed exceeds the headroom beside it" -- except the two were
    not beside each other, and nothing performed the comparison. Telling a
    reader to cross-reference two tables is how a conclusion goes unread.

    The ceiling is the channel's REACH, not its headroom: a compound the
    metapath never arrives at cannot enter the regression however well it
    is screened. Labelling every reached compound is the best any panel
    can do over a fixed graph.
    """
    if not fitted:
        print("  No cell was fitted, so there is nothing to size.")
        return
    best: dict[tuple[str, str], tuple[int, int | None]] = {}
    for virus, chan, reached, need in fitted:
        key = (virus, chan)
        prev = best.get(key)
        if prev is None:
            best[key] = (reached, need)
        elif need is not None and (prev[1] is None or need < prev[1]):
            best[key] = (reached, need)

    print("  CAN LABELS SETTLE IT? ceiling is the channel's reach -- a")
    print("  compound the metapath never arrives at cannot enter the")
    print("  regression however well it is screened.")
    print()
    print(f"  {'virus':<18} {'chan':<4} {'ceiling':>8} {'needed':>10} "
          f"{'ratio':>7}  verdict")
    settleable = 0
    for (virus, chan), (reached, need) in sorted(best.items()):
        if need is None:
            print(f"  {virus:<18} {chan:<4} {reached:8,} {'>1M':>10} "
                  f"{'--':>7}  NO -- effect too small to chase with labels")
            continue
        ok = need <= reached
        settleable += ok
        print(f"  {virus:<18} {chan:<4} {reached:8,} {need:10,} "
              f"{need / reached:6.2f}x  "
              f"{'YES' if ok else 'NO -- more graph, not more labels'}")
    print()
    print(f"  {settleable} of {len(best)} cells could be settled by labelling")
    print("  every compound its channel already reaches.")
    if settleable == 0:
        print()
        print("  NO PANEL OVER THIS GRAPH SETTLES ANY CELL. The binding")
        print("  constraint is reach, not labels: the compounds a decisive")
        print("  test would need are ones no metapath arrives at. A second")
        print("  screening panel helps only insofar as it brings compounds")
        print("  the graph can already route to, and the measurement above")
        print("  says there are not enough of those.")


if __name__ == "__main__":
    raise SystemExit(main())
