#!/usr/bin/env python3
"""Does DWPC weighting read the answer sheet, and by how much.

    python scripts/diagnose_label_degree_weighting.py
    python scripts/diagnose_label_degree_weighting.py --channel M7

READ-ONLY. Writes no results file and changes no registered protocol.
It measures the SENSITIVITY of every metapath channel to one choice
already documented in baselines.py, so the choice can be made on a
number instead of on an argument.

THE QUESTION. `walk_metapath_detail` weights each node the walk ENTERS
by `g.degree ** -w`, and `Graph.degree` counts edges of EVERY predicate
including the held-out labels. That is deliberate -- baselines.py says
so, for hub exclusion and degree-matched binning -- and it is the same
property that made the degree COMPARATOR in the three scorer studies
defective (docs/prereg_structure_scorer.md, amendment 2026-11-17). The
comparator was a bug because `degree_ranking` promised otherwise. This
is not a bug: nothing promised otherwise. It is an open question with a
measurable answer.

WHY IT COULD MATTER, and only for some channels. Held-out predicates
touch three node classes:

  HAS_ANTIVIRAL_ACTIVITY_AGAINST   SmallMolecule -> OrganismTaxon
  MEASURED_INACTIVE_AGAINST        SmallMolecule -> Protein
  TREATS / IN_TRIAL_FOR            SmallMolecule -> Disease

A virus's own degree is a CONSTANT within its cell, so it scales every
drug's score alike and cannot move a ranking. What can move a ranking is
a varying intermediate: a borrowed SmallMolecule neighbour reached over
CHEMICALLY_SIMILAR_TO, and a viral Protein carrying inactive
measurements. So M6 and M7 (similarity hop) and M1/M7/M8 (viral protein)
are the channels at risk, and the rest are a control that should not
move at all.

**M7 IS WHY THIS IS WORTH THE RUN.** It is the only readable cell that
clears the corrected alpha in docs/evaluation.md -- SARS-CoV-2 M7 at
p=0.000841 against 0.05/13 -- and its path enters a borrowed compound
whose degree is mostly label edges.

THREE ARMS, so the weighting is isolated from the hub cut:

  A  as published        g.degree for the step weight, hub set from it
  B  labels out of the   degree_unlabelled for the step weight, HUB SET
     step weight only    HELD AT A's.  <- the isolated measurement
  C  labels out of both   degree_unlabelled for both

B is the headline. C is reported because `_hub_set` reads the same
counter, so a reader asking "what if we took the labels out everywhere"
should not have to guess.

THE ITERATION OVER DRUGS IS NOT `dwpc_scores`, and that is the one
duplication here. `dwpc_scores` recomputes `blocked` from the graph it
is handed, which is exactly the variable arm B has to hold fixed. The
walk itself, the hub set, the oversized-pathway set and the metric are
all imported, so nothing that produces a number is copied.
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import (  # noqa: E402
    DEFAULT_W,
    Graph,
    _hub_set,
    _oversized_pathways,
    lift_row,
    walk_metapath_detail,
)
from kgav.labels import build_labels  # noqa: E402
from kgav.pipeline import release_dir  # noqa: E402
from kgav.provenance import write_results  # noqa: E402
from kgav.schema import load_schema, print_capped  # noqa: E402

ARMS = ("A_as_published", "B_weight_only", "C_weight_and_hubs")

#: THE CELL docs/evaluation.md's claim section rests on: the only
#: readable cell that clears the corrected alpha. Named here so the
#: artifact carries it under a stable key instead of a reader picking it
#: out of sixty-three rows.
CLAIM_VIRUS = "NCBITaxon:2697049"
CLAIM_CHANNEL = "M7"


def score_all(g: Graph, metapaths: dict, blocked: set[str],
              w: float = DEFAULT_W) -> dict[str, dict[str, dict[str, float]]]:
    """{metapath: {drug: {virus: score}}} with a CALLER-SUPPLIED hub set."""
    out: dict[str, dict[str, dict[str, float]]] = {}
    for name, mp in metapaths.items():
        per_drug: dict[str, dict[str, float]] = {}
        for drug in g.drugs():
            hits, _cut = walk_metapath_detail(
                g, drug, mp["path"], mp.get("constraints"),
                blocked if mp.get("exclude_hubs") else set(), w)
            if hits:
                per_drug[drug] = hits
        out[name] = per_drug
    return out


def rows(dwpc: dict, positives: dict[str, set[str]], k: int) -> dict:
    """{(virus, channel): lift_row} over every evaluable cell."""
    out: dict[tuple[str, str], dict] = {}
    for channel, per_drug in dwpc.items():
        viruses = {v for hits in per_drug.values() for v in hits}
        for virus in viruses:
            pos = positives.get(virus, set())
            if not pos:
                continue
            scores = {d: h[virus] for d, h in per_drug.items() if virus in h}
            out[(virus, channel)] = lift_row(scores, pos, k)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--channels", default=None,
                    help="comma-separated metapaths, e.g. M1,M6,M7 "
                         "(default: all)")
    ap.add_argument("-k", type=int, default=100)
    ap.add_argument("--out", type=Path, default=None,
                    help="write the three arms as a STAMPED results file. "
                         "docs/evaluation.md quotes this run beside its own "
                         "M7 claim, and a figure in that document has to be "
                         "derivable by a test rather than typed, so the "
                         "sensitivity ships as an audited artifact")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    s = load_schema()
    held = s.held_out_predicates()
    g = Graph.load(args.release, skip_predicates=held, symmetry=s.symmetry())
    metapaths = dict(s.metapaths)
    if args.channels:
        want = [c.strip() for c in args.channels.split(",") if c.strip()]
        unknown = [c for c in want if c not in metapaths]
        if unknown:
            print(f"no metapath(s) {unknown}; have "
                  f"{', '.join(sorted(metapaths))}")
            return 2
        metapaths = {c: metapaths[c] for c in want}
    positives = build_labels(args.release).positives
    names = {n: d["properties"]["label"] for n, d in g.nodes.items()
             if d["class"] == "OrganismTaxon"}

    true_degree = Counter(g.degree)
    unlabelled = Counter(g.degree_unlabelled)
    moved = {n for n in true_degree if true_degree[n] != unlabelled[n]}
    by_class: Counter = Counter(g.klass(n) for n in moved)
    print(f"release {args.release.name}   held out: "
          f"{', '.join(sorted(held))}")
    print(f"w = {DEFAULT_W}   channels: {', '.join(sorted(metapaths))}")
    print(f"\n{len(moved):,} node(s) carry held-out degree, by class:")
    for klass, n in by_class.most_common():
        print(f"    {klass:<16} {n:>8,}")

    blocked_a = _hub_set(g, s.retrieval_limits) | \
        _oversized_pathways(g, s.retrieval_limits)

    arms: dict[str, dict] = {}
    for arm in ARMS:
        if arm == "A_as_published":
            g.degree = Counter(true_degree)
            blocked = blocked_a
        elif arm == "B_weight_only":
            g.degree = Counter(unlabelled)
            blocked = blocked_a
        else:
            g.degree = Counter(unlabelled)
            blocked = _hub_set(g, s.retrieval_limits) | \
                _oversized_pathways(g, s.retrieval_limits)
            blocked_c = blocked
        print(f"\n--- {arm}: {len(blocked):,} node(s) excluded as path "
              f"connectors")
        arms[arm] = rows(score_all(g, metapaths, blocked), positives, args.k)
    g.degree = Counter(true_degree)

    if blocked_a != (_hub_set(g, s.retrieval_limits) |
                     _oversized_pathways(g, s.retrieval_limits)):
        print("! the hub set is not reproducible from the restored degree; "
              "arms are not comparable")
        return 1

    a, b, c = (arms[x] for x in ARMS)
    print("\n" + "=" * 78)
    print("SENSITIVITY OF EVERY EVALUABLE CELL. 'p' is the hypergeometric "
          "tail\nlift_row computes; (A) is what data/results holds.")
    print("=" * 78)
    head = (f"{'virus':<14}{'ch':<5}{'pool':>6}{'lift A':>9}{'lift B':>9}"
            f"{'p A':>11}{'p B':>11}{'p C':>11}")
    print(head)
    print("-" * len(head))
    worst = []
    for key in sorted(a, key=lambda k: (names.get(k[0], k[0]), k[1])):
        ra, rb, rc = a[key], b.get(key), c.get(key)
        if rb is None:
            print(f"{names.get(key[0], key[0]):<14}{key[1]:<5}"
                  f"{ra['pool']:>6}  CELL ABSENT IN ARM B")
            continue
        print(f"{names.get(key[0], key[0]):<14}{key[1]:<5}{ra['pool']:>6}"
              f"{ra['lift']:>9.2f}{rb['lift']:>9.2f}"
              f"{_p(ra):>11}{_p(rb):>11}{_p(rc):>11}")
        if ra.get("p_value") is not None and rb.get("p_value") is not None:
            worst.append((abs(rb["p_value"] - ra["p_value"]), key, ra, rb))

    # print_capped, NOT a bare slice: the remainder is counted. Here the
    # tail is the cells that did NOT move, so "and N more" is the control
    # this whole measurement rests on -- a head-and-stop would print the
    # four movers and hide the fifty-seven that are the evidence the
    # finding is small.
    print("\nTHE CELLS THAT MOVE MOST (|p_B - p_A|):")
    print_capped(
        sorted(worst, reverse=True), 6,
        lambda row: (f"{row[0]:.6f}   {names.get(row[1][0], row[1][0])} "
                     f"{row[1][1]}   p {row[2]['p_value']:.6f} -> "
                     f"{row[3]['p_value']:.6f}   lift {row[2]['lift']:.2f} "
                     f"-> {row[3]['lift']:.2f}"),
        unit="more cell(s), in descending order of movement")
    crossed = [(k, ra, rb) for _d, k, ra, rb in worst
               if (ra["p_value"] < 0.05) != (rb["p_value"] < 0.05)]
    print(f"\n{len(crossed)} cell(s) cross p=0.05 between arm A and arm B"
          + (":" if crossed else "."))
    for key, ra, rb in crossed:
        print(f"    {names.get(key[0], key[0])} {key[1]}: "
              f"{ra['p_value']:.6f} -> {rb['p_value']:.6f}")

    if args.out:
        # THE HEADLINE IS NAMED, not left to a reader to find among 63
        # cells: docs/evaluation.md's claim section rests on one cell and
        # tests/test_degree_sensitivity.py reads it from here.
        head = a.get((CLAIM_VIRUS, CLAIM_CHANNEL))
        head_b = b.get((CLAIM_VIRUS, CLAIM_CHANNEL))
        payload = {
            "measured_by": "scripts/diagnose_label_degree_weighting.py",
            "registers_nothing": (
                "A SENSITIVITY, NOT A PROTOCOL. Nothing promised that DWPC "
                "step weights exclude held-out edges and baselines.py "
                "argues the opposite on purpose, so switching the "
                "weighting is a registered decision. This file says what "
                "such a registration would be deciding about."),
            "w": DEFAULT_W,
            "k": args.k,
            "held_out": sorted(held),
            "channels": sorted(metapaths),
            "arms": {
                "A_as_published": "g.degree for the step weight; hub set "
                                  "from it. What data/results holds.",
                "B_weight_only": "degree_unlabelled for the step weight, "
                                 "hub set HELD AT A's. The isolated change.",
                "C_weight_and_hubs": "degree_unlabelled for both.",
            },
            "held_out_degree_by_class": dict(by_class),
            "hub_set_size": {"A": len(blocked_a), "B": len(blocked_a),
                             "C": len(blocked_c)},
            "cells": len(a),
            "cells_moved": sum(
                1 for k in a if b.get(k) and
                (a[k]["lift"], a[k].get("p_value")) !=
                (b[k]["lift"], b[k].get("p_value"))),
            "cells_crossing_alpha": len(crossed),
            "claim_cell": {
                "virus": CLAIM_VIRUS, "channel": CLAIM_CHANNEL,
                "p_a": head and head.get("p_value"),
                "p_b": head_b and head_b.get("p_value"),
                "lift_a": head and head.get("lift"),
                "lift_b": head_b and head_b.get("lift"),
            },
            "per_cell": {f"{v}|{ch}": {
                "A": a[(v, ch)], "B": b.get((v, ch)), "C": c.get((v, ch))}
                for (v, ch) in a},
        }
        write_results(args.out, payload, args.release, s.version)
        print(f"\nwrote {args.out}")
    print("\nNOT A PROTOCOL CHANGE. Switching the weighting would restate "
          "every\nfigure in data/results, so it needs its own registration; "
          "this says\nwhat such a registration would be deciding about.")
    return 0


def _p(r: dict) -> str:
    v = r.get("p_value")
    return "--" if v is None else f"{v:.6f}"


if __name__ == "__main__":
    raise SystemExit(main())
