#!/usr/bin/env python
"""How much of each channel's score is the compound's own degree.

WHY THIS EXISTS. Every table in docs/evaluation.md compares a channel against
a degree baseline, and the verdict of the whole document is that degree wins.
That comparison assumed the two were different quantities. They are not
entirely: `walk_metapath`'s default scoring never divides out the START
node's degree, so a compound with more edges accumulates more score, and the
channel partly IS the baseline it is being compared against.

docs/prereg_m1_potency.md already records the half of this that was noticed:
"nothing drug-specific ever enters", hence 12 distinct M1 values across 3,020
drugs. Two consequences were not recorded, and both are measured here:

  - the score grows with the compound's own edge count, so `rho(score,
    degree)` is high by construction rather than by biology
  - `g.degree` counts every predicate, so a layer the metapath never walks
    still moves the score -- adding the fold layer shifted M1 without
    anything in M1 changing

baselines.SCORING_PDP is the published path-degree product: every node on the
path, once per incident metaedge, at that metaedge's own degree. This script
runs both and prints the difference. It decides nothing. The default scoring
is unchanged and every number in data/results still comes from the legacy
form; switching it is a separate decision that needs this measurement first.

WHAT TO READ

  rho(score, degree)   Spearman between a channel's score and the compound's
                       global node degree, per scoring. The closer to 1, the
                       more the channel is restating the baseline. A channel
                       whose rho FALLS under pdp was carrying degree.

  distinct / tie       the degeneracy table from the prereg, recomputed under
                       both. Expect it NOT to improve: two compounds with one
                       edge each to the same protein have equal degrees on
                       both sides of every metaedge, so they stay tied. If
                       this column moves much, something is wrong.

  rho(ranks) / top-k   whether the change reaches the ranking at all. A
                       channel that reorders nothing cannot change a result.

USAGE

  python scripts/diagnose_scorer.py
  python scripts/diagnose_scorer.py --virus NCBITaxon:64320
  python scripts/diagnose_scorer.py --all-viruses

Read-only. Writes nothing, recommends nothing.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import (  # noqa: E402
    SCORING_LEGACY,
    SCORING_PDP,
    Graph,
    degree_ranking,
    dwpc_scores,
)
from kgav.schema import load_schema  # noqa: E402

TOP_K = 100


# ----------------------------------------------------------------- statistics
def ranks(xs: list[float]) -> list[float]:
    """Competition-free average ranks, which ties REQUIRE.

    Assigning 1..n in sort order would make two equal scores differ by a
    whole rank, and this repository's central channel puts 80.8% of its pool
    in one tie block -- the correction is not a detail here, it is most of
    the data.
    """
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    out = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def spearman(xs: list[float], ys: list[float]) -> float | None:
    """Pearson on average ranks. None when either side is constant.

    None rather than 0.0: a constant vector has no correlation defined, and
    printing 0.00 would read as "measured, and unrelated".
    """
    if len(xs) != len(ys):
        raise ValueError("spearman needs paired vectors")
    n = len(xs)
    if n < 2:
        return None
    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def degeneracy(scores: list[float]) -> tuple[int, int, float]:
    """(distinct values, largest tie block, % of pool inside any tie)."""
    counts: dict[float, int] = {}
    for s in scores:
        counts[s] = counts.get(s, 0) + 1
    tied = sum(c for c in counts.values() if c > 1)
    largest = max(counts.values()) if counts else 0
    pct = 100.0 * tied / len(scores) if scores else 0.0
    return len(counts), largest, pct


def topk(scores: dict[str, float], k: int) -> set[str]:
    """The k highest-scoring drugs.

    Ties at the cut make this one draw from a larger set -- the same defect
    baselines.hits_at_k_range exists for. The overlap figure below is
    therefore read as an upper bound on agreement, and the caller says so.
    """
    return {d for d, _ in sorted(scores.items(), key=lambda kv: -kv[1])[:k]}


# --------------------------------------------------------------------- report
def fmt_rho(r: float | None) -> str:
    return "  n/d " if r is None else f"{r:+.3f}"


def channel_table(name: str, legacy: dict[str, float], pdp: dict[str, float],
                  degree: dict[str, float]) -> str:
    shared = sorted(set(legacy) & set(pdp))
    if not shared:
        return f"  {name:<5} no drug reached"
    lv = [legacy[d] for d in shared]
    pv = [pdp[d] for d in shared]
    dv = [degree.get(d, 0) for d in shared]
    ld, llt, lpct = degeneracy(lv)
    pd_, plt, ppct = degeneracy(pv)
    overlap = len(topk(legacy, TOP_K) & topk(pdp, TOP_K))
    return (f"  {name:<5} {len(shared):6,} "
            f"{fmt_rho(spearman(lv, dv))} {fmt_rho(spearman(pv, dv))}  "
            f"{ld:7,} {pd_:7,}  {llt:6,} {plt:6,}  "
            f"{fmt_rho(spearman(lv, pv))}  {overlap:3d}/{TOP_K}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path, default=ROOT / "data/releases/v0.1")
    ap.add_argument("--virus", default="NCBITaxon:2697049")
    ap.add_argument("--all-viruses", action="store_true")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    schema = load_schema()
    held = schema.held_out_predicates()
    g = Graph.load(args.release, skip_predicates=held,
                   symmetry=schema.symmetry())
    degree = degree_ranking(g)

    print(f"release {args.release.name} | {len(g.nodes):,} nodes | "
          f"{len(g.drugs()):,} drugs | held out: {sorted(held)}")
    print(f"legacy = {SCORING_LEGACY}  (what every results file used)")
    print(f"pdp    = {SCORING_PDP}  (Himmelstein, w unchanged)\n")

    both = {}
    for scoring in (SCORING_LEGACY, SCORING_PDP):
        both[scoring] = dwpc_scores(g, schema.metapaths,
                                    schema.retrieval_limits, scoring=scoring)

    reached = sorted({v for per in both[SCORING_LEGACY].values()
                      for hits in per.values() for v in hits})
    viruses = reached if args.all_viruses else [args.virus]

    # A REPORT THAT PRINTS NOTHING MUST SAY SO. --all-viruses over a graph no
    # metapath completes produced a header, two trailing paragraphs and not
    # one table, which reads as "nothing to flag" rather than "nothing was
    # measured". Same for a --virus no channel reaches.
    if not viruses:
        print("NO VIRUS REACHED by any metapath, so there is nothing to "
              "compare.")
        print(f"  {len(schema.metapaths)} metapaths declared, "
              f"{len(g.drugs()):,} drugs, 0 completed walks.")
        return 1
    if not args.all_viruses and args.virus not in reached:
        print(f"NOTE: no metapath reaches {args.virus} on this release. "
              f"Reached: {reached or 'none'}")
        print("The table below is therefore empty by construction, not by "
              "a null result.\n")

    for virus in viruses:
        print(f"=== {virus}")
        print(f"  {'chan':<5} {'drugs':>6} {'rho(score,deg)':>14}  "
              f"{'distinct':>15}  {'largest tie':>13}  {'rho(L,P)':>8} "
              f"{'top':>7}")
        print(f"  {'':<5} {'':>6} {'legacy':>6} {'pdp':>7}  "
              f"{'legacy':>7} {'pdp':>7}  {'legacy':>6} {'pdp':>6}  "
              f"{'':>8} {'overlap':>7}")
        for name in sorted(schema.metapaths):
            legacy = {d: hits[virus]
                      for d, hits in both[SCORING_LEGACY].get(name, {}).items()
                      if virus in hits}
            pdp = {d: hits[virus]
                   for d, hits in both[SCORING_PDP].get(name, {}).items()
                   if virus in hits}
            print(channel_table(name, legacy, pdp, degree))
        print()

    print("Top-k overlap counts the drugs in both top-100 sets. Ties at the")
    print("cut make each set one draw from a larger one, so the figure is an")
    print("upper bound on agreement, not a measurement of it.")
    print()
    print("Nothing here changes a default. data/results is still the legacy")
    print("scoring; switching it is a separate registered decision.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
