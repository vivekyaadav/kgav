#!/usr/bin/env python
"""Why each biology control ranks where it does.

run_baselines.py states the expectation at its CONTROLS block:

    nirmatrelvir  nsp5 inhibitor  -> must rank high via M1
    remdesivir    nsp12 inhibitor -> must rank high via M1
    hydroxychloroquine / chloroquine: ... If these rank high the graph is
    reproducing a cell-line artifact

On release v0.1 the cross-sectional DWPC ranking returns, out of 6,515:

    nirmatrelvir        6,432     positive   <- near the BOTTOM
    remdesivir            396     positive
    hydroxychloroquine  1,146     artifact
    chloroquine         1,328     artifact

Both artifacts outrank the one approved direct-acting antiviral. The controls
exist to catch exactly this and they caught it.

WHAT THIS SCRIPT IS FOR. The explanation is arithmetic, not mystery, and it
should be read off the numbers rather than reasoned about -- baselines.combine
already carries three recorded rank histories for nirmatrelvir (2,973 under
raw summation, 1,428 under sum-of-percentiles, and a comment predicting that
max-of-percentiles "is biased upward by the number of channels ... harmless
while every channel carries signal and corrosive once some do not"). Since
calibration now reports that NO channel's interval clears chance, that
prediction is a description of the present. This prints the pieces so the
claim is checkable:

  - which channels reach the virus for this drug at all
  - its raw DWPC in each, and its PERCENTILE in that channel's pool
  - the combined score, which is the max of those percentiles
  - its own node degree, and the degree of the proteins its paths traverse

DWPC divides by node degree raised to w, so a route through a heavily studied
protein scores LOW by construction. nsp5 carries the largest number of
resolved activities in the graph, which is the same reason nirmatrelvir is the
best-validated compound and the reason its one route is penalised hardest.
Whether that is the whole story is what the output decides.

The side-by-side that matters is nirmatrelvir against pomotrelvir: both are
Mpro inhibitors, so both route through nsp5, and on v0.1 pomotrelvir ranks 17
while nirmatrelvir ranks 6,432. Same target, same channel, opposite ends.

USAGE

  python scripts/diagnose_controls.py
  python scripts/diagnose_controls.py --virus NCBITaxon:64320
  python scripts/diagnose_controls.py --also pomotrelvir ebselen

Read-only. Writes nothing, recommends nothing.
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import Graph, channel_percentiles, combine, degree_ranking, dwpc_scores, rank
from kgav.labels import build_labels
from kgav.schema import load_schema


def _run_baselines():
    """run_baselines' own CONTROLS and name_index, loaded by path.

    Imported rather than restated so the controls and the name resolution are
    the SAME ones the reported ranking used. A second copy of either would let
    this script explain a ranking that was never computed.
    """
    spec = importlib.util.spec_from_file_location(
        "run_baselines", ROOT / "scripts" / "run_baselines.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


def hits_range(scores: dict[str, float], positives: set[str], k: int):
    """The range Hits@k can take over tie orderings, and the tie at the cut.

    rank() is `sorted(scores.items(), key=lambda kv: -kv[1])` with NO
    tiebreak, and Python's sort is stable, so equal scores come out in dict
    insertion order. hits_at_k then reads ranked[:k]. When a tie block
    STRADDLES position k, which of its members fall inside the cut is decided
    by insertion order rather than by any score -- so Hits@k, MRR and lift
    are reproducible for one dict and not a property of the graph.

    Returned as (low, high, contested, contested_positives): low is the count
    if every tied non-positive is placed first, high if every tied positive
    is. low == high means the cut lands in no tie and the metric is a
    measurement. A wide interval means the published number is one draw from
    it.
    """
    ordered = sorted(scores.items(), key=lambda kv: -kv[1])
    if len(ordered) <= k:
        n = sum(1 for d, _v in ordered if d in positives)
        return n, n, 0, 0
    cut = ordered[k - 1][1]
    inside = [d for d, v in ordered if v > cut]
    contested = [d for d, v in ordered if v == cut]
    c_pos = sum(1 for d in contested if d in positives)
    fixed = sum(1 for d in inside if d in positives)
    slots = k - len(inside)
    low = fixed + max(0, slots - (len(contested) - c_pos))
    high = fixed + min(slots, c_pos)
    # A block of ONE is not a tie: the drug at the cut is simply the kth. It
    # is reported as 0 contested so the column does not imply an ambiguity
    # that is not there.
    if len(contested) == 1:
        return low, high, 0, 0
    return low, high, len(contested), c_pos


def main(argv=None) -> int:
    rb = _run_baselines()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path, default=ROOT / "data/releases/v0.1")
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--virus", default="NCBITaxon:2697049")
    ap.add_argument("--also", nargs="*", default=["pomotrelvir"],
                    help="extra compounds to show alongside the controls")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        sys.exit(f"no release at {args.release}")

    schema = load_schema()
    held = schema.held_out_predicates()
    g = Graph.load(args.release, skip_predicates=held,
                   symmetry=schema.symmetry())
    print(f"release {args.release.name} | held out of traversal: "
          f"{sorted(held)}")
    print(f"{len(g.nodes):,} nodes | {len(g.drugs()):,} drugs\n")

    dwpc = dwpc_scores(g, schema.metapaths, schema.retrieval_limits)
    scores = combine(dwpc, args.virus)
    ranked = rank(scores)
    pos_of = {d: i for i, (d, _v) in enumerate(ranked, 1)}
    deg = degree_ranking(g)
    deg_rank = {d: i for i, (d, _v) in enumerate(rank(deg), 1)}

    # baselines.channel_percentiles, NOT a second copy of its arithmetic.
    # This script used to reimplement it, and when combine moved to midranks
    # the copy did not: it printed a percentile of 0.0182 beside a combined
    # score of 0.4041 for a drug whose only route is one channel -- two
    # numbers that cannot both be right, from the script written to make the
    # scoring checkable. There is now one implementation and nothing to sync.
    pools: dict[str, dict[str, float]] = {}
    pcts: dict[str, dict[str, float]] = {}
    for name, per_drug in dwpc.items():
        pool = {d: by[args.virus] for d, by in per_drug.items()
                if args.virus in by}
        if not pool:
            continue
        pools[name] = pool
        pcts[name] = channel_percentiles(pool)

    # HOW MUCH OF EACH CHANNEL IS ONE NUMBER. DWPC is topology only, so every
    # drug whose sole route is the same single node earns the same product of
    # the same degrees. If most of a pool is one value the channel is not
    # ranking, it is grouping -- and the midrank every member then shares is
    # the honest score, not a defect. This is the measurement that decides
    # whether potency has to enter the direct-acting channel.
    print(f"channel pools for {args.virus}")
    print(f"  {'chan':<6}{'drugs':>8}{'distinct':>10}{'largest tie':>13}"
          f"{'% of pool':>11}{'% in any tie':>14}")
    for name in sorted(pools):
        vals = Counter(pools[name].values())
        _big_val, big_n = vals.most_common(1)[0]
        tied = sum(c for c in vals.values() if c > 1)
        size = len(pools[name])
        print(f"  {name:<6}{size:>8,}{len(vals):>10,}{big_n:>13,}"
              f"{100 * big_n / size:>10.1f}%{100 * tied / size:>13.1f}%")

    names = rb.name_index(g, args.db)
    if not names:
        sys.exit("no name index (needs the ChEMBL db) -- cannot resolve "
                 "control compounds by name")

    wanted = list(rb.CONTROLS.items()) + [(c, "comparison") for c in args.also]
    print()
    for cname, kind in wanted:
        nid = names.get(cname)
        print("=" * 78)
        print(f"{cname}  [{kind}]")
        print("=" * 78)
        if not nid:
            print("  not in the graph under this name\n")
            continue

        print(f"  node {nid}")
        print(f"  DWPC rank {pos_of.get(nid, '-'):>8} of {len(ranked):,}"
              f"   combined score {scores.get(nid, 0.0):.4f}")
        print(f"  degree rank {deg_rank.get(nid, '-'):>6} of "
              f"{len(deg):,}   raw degree {g.degree.get(nid, 0):,}")

        reached = [n for n in sorted(pools) if nid in pools[n]]
        if not reached:
            print("  reaches this virus by NO channel, so it has no DWPC "
                  "route at all\n")
            continue
        print(f"\n  {'channel':<8}{'raw DWPC':>12}{'percentile':>12}"
              f"{'pool':>9}   (combined = max percentile)")
        for n in reached:
            print(f"  {n:<8}{pools[n][nid]:>12.6f}{pcts[n][nid]:>12.4f}"
                  f"{len(pools[n]):>9,}")
        print(f"  channels reached: {len(reached)} of {len(pools)}"
              f"   -> max percentile {max(pcts[n][nid] for n in reached):.4f}")

        # WHAT THE PATHS RUN THROUGH, because DWPC's denominator is node
        # degree: a route through the most-studied protein in the graph is
        # penalised hardest, and that protein is the reason the compound is
        # well validated in the first place.
        # g.neighbours is the PUBLIC walk, which is what walk_metapath uses
        # and which honours declared symmetry. g.out holds (neighbour, edge)
        # TUPLES rather than edge dicts -- reading e["object"] off it throws,
        # and the first version of this script did exactly that.
        #
        # g.degree, not degree_ranking: degree_ranking is
        # `{d: degree for d in g.drugs()}`, so a PROTEIN looked up in it
        # returns 0 and this column would have reported every target as
        # degree 0 -- the silent zero that makes a diagnostic lie.
        #
        # target_resolution is printed because it now DECIDES traversal. M1,
        # M7 and M8 require chain_from_assay_description, so an edge marked
        # parent_unresolved is no longer walked -- and the degree beside it
        # explains why that mattered: the parent accession carries fewer
        # edges than the mature chain, so the unresolved edge used to score
        # HIGHER. The parent's degree is low precisely BECAUSE resolution
        # succeeded elsewhere, which made it a denominator rewarding its own
        # gaps. Shown per edge rather than per target: one drug can hold both
        # kinds.
        inter: Counter = Counter()
        res: dict[str, Counter] = {}
        for pred in ("INHIBITS", "TARGETS"):
            for nbr, edge in g.neighbours(pred, nid):
                inter[nbr] += 1
                r = (edge.get("qualifiers") or {}).get("target_resolution")
                res.setdefault(nbr, Counter())[r or "-"] += 1
        if inter:
            print("\n  first-hop targets, with the degree DWPC divides by:")
            for tgt, n_edges in inter.most_common(6):
                props = (g.nodes.get(tgt) or {}).get("properties") or {}
                lbl = (props.get("protein_family") or props.get("gene_symbol")
                       or tgt)
                how = ", ".join(f"{k}x{v}" if v > 1 else k
                                for k, v in res[tgt].most_common())
                walks = "parent_unresolved" not in res[tgt]
                print(f"    {lbl!s:<22} degree {g.degree.get(tgt, 0):>8,}"
                      f"  ({n_edges} edge(s))  {'M1 walks' if walks else 'M1 REFUSES'}"
                      f"  [{how}]  {tgt}")
        print()

    # ------------------------------------------------- the aggregate check
    print("=" * 78)
    print("DOES RANK TRACK THE NUMBER OF CHANNELS REACHED?")
    print("=" * 78)
    print("""  combine takes the MAX of a drug's percentiles across channels, and its
  own comment predicts the consequence: "the maximum of k independent uniform
  percentiles has expectation k/(k+1): over eight channels that is 0.89 for
  any drug present in all of them, so the ordering reverts to counting pools".
  If that is what is happening, mean rank should fall steeply as the channel
  count rises, and a one-channel drug cannot compete however good its route.
""")
    by_count: dict[int, list[int]] = {}
    for d, _s in ranked:
        c = sum(1 for n in pools if d in pools[n])
        by_count.setdefault(c, []).append(pos_of[d])
    print(f"  {'channels':>9}{'drugs':>9}{'best rank':>11}{'median rank':>13}")
    for c in sorted(by_count):
        rs = sorted(by_count[c])
        print(f"  {c:>9}{len(rs):>9,}{rs[0]:>11,}{rs[len(rs) // 2]:>13,}")

    # ------------------------------------- is the metric a measurement?
    print()
    print("=" * 78)
    print("IS Hits@k A MEASUREMENT, OR AN ARTIFACT OF TIE ORDER?")
    print("=" * 78)
    print("""  rank() sorts by score with NO tiebreak and Python's sort is stable, so
  equal scores emerge in dict insertion order; hits_at_k then takes
  ranked[:k]. Where a tie block straddles the cut at k, which members fall
  inside it is decided by insertion order and not by any score. The range
  below is what Hits@k could be over tie orderings: low places every tied
  NON-positive first, high every tied positive. low == high means the cut
  lands in no tie and the number is real.
""")
    labels = build_labels(args.release)
    positives = labels.positives.get(args.virus, set())
    print(f"  {len(positives):,} compounds measured active against this virus")
    print(f"\n  {'ranking':<10}{'k':>5}{'Hits@k range':>16}{'tied at cut':>13}"
          f"{'of which pos':>14}")
    rankings = {"COMBINED": scores}
    for name in sorted(pools):
        rankings[name] = pools[name]
    for rname, sc in rankings.items():
        for k in (10, 50, 100):
            lo, hi, nc, npos = hits_range(sc, positives, k)
            flag = "" if lo == hi else "   <- NOT DETERMINED"
            print(f"  {rname:<10}{k:>5}{f'{lo} - {hi}':>16}{nc:>13,}"
                  f"{npos:>14,}{flag}")

    print("""
  HOW TO READ THIS
    A control ranking badly with ONE channel and a good raw DWPC is the
    max-over-channels bias: it is being outranked by drugs that appear in
    more pools, not by drugs with better routes.

    A control ranking badly with a LOW percentile in its own channel is the
    degree correction: DWPC divides by node degree, so the better studied its
    target, the lower its score. nsp5 carries the most resolved activities in
    the graph.

    These are different problems with different fixes, and the numbers above
    say which one this is. calibration.combine(..., calibration=cal) already
    restricts the max to channels whose interval clears chance -- on this
    release none do, so it returns no ranking at all, which calibration.py
    argues is the correct output. run_baselines does not pass it.

  This script measures and recommends nothing.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
