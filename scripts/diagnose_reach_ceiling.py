#!/usr/bin/env python
"""Is a channel's reach limited by missing data or by a retrieval parameter.

label_headroom.py measured that no panel over this graph settles any cell,
and that SARS-CoV-2 M4 is the near miss: it needs 3,869 labelled compounds
against a ceiling of 2,982 reached -- short by 887, a factor of 1.30. So
the question is what sets 2,982, and there are two very different answers.

MISSING DATA. A hop starves because the layer under it is thin, and the fix
is an ingest.

A RETRIEVAL PARAMETER. M4, M5 and M8 declare `exclude_hubs`, so their
INTERMEDIATE nodes are filtered by two schema limits: `hub_degree_cut`
(460, Protein class) and `max_pathway_members` (200). Both exist for good
reasons -- a GO term with 1,389 members makes every drug adjacent to every
host factor -- but both are numbers someone chose, and a reach ceiling set
by a chosen number is a different problem from one set by absent biology.

This prints, per hop, how many nodes survive; and then the reach under four
exclusion settings, so the cost of each limit is a measurement rather than
an assumption.

READ THE COUNTERFACTUALS AS A COST, NOT A RECOMMENDATION. Relaxing the hub
cut is exactly what the cut was added to prevent, and raising either limit
changes every number in data/results, so it would need a registration.
What this script decides is only whether that conversation is worth having
for a given channel, or whether the answer is an ingest instead.

USAGE

  python scripts/diagnose_reach_ceiling.py
  python scripts/diagnose_reach_ceiling.py --channel M5 --virus NCBITaxon:64320

Read-only. Writes nothing, changes no limit.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import (  # noqa: E402
    Graph,
    _hub_set,
    _oversized_pathways,
    walk_metapath,
)
from kgav.schema import load_schema, parse_hop  # noqa: E402

# From the label_headroom run on v0.1, so the gap this is explaining is
# stated rather than remembered. Keyed (virus, channel) -> (reach, needed).
MEASURED_ON = "v0.1"
MEASURED_GAP = {
    ("NCBITaxon:2697049", "M4"): (2982, 3869),
    ("NCBITaxon:2697049", "M5"): (2192, 6504),
    ("NCBITaxon:2697049", "M1"): (3021, 17748),
}


def survivors(g: Graph, path: list, constraints: dict | None,
              blocked: set[str], starts: list[str]) -> list[tuple[str, int]]:
    """Distinct nodes alive after each hop, walking the whole drug set.

    Counts NODES, not paths: the reach ceiling is how many compounds can
    arrive, and a hop that multiplies paths without widening the frontier
    does not raise it.
    """
    steps = [(parse_hop(path[i]), path[i + 1])
             for i in range(1, len(path) - 1, 2)]
    frontier = {s: {s} for s in starts}     # node -> originating compounds
    out = [(path[0], len(frontier))]
    for step_i, ((predicate, reverse), want) in enumerate(steps):
        last = step_i == len(steps) - 1
        nxt: dict[str, set[str]] = defaultdict(set)
        for node, origins in frontier.items():
            for nbr, edge in g.neighbours(predicate, node, reverse):
                if g.klass(nbr) != want:
                    continue
                if constraints:
                    q = edge.get("qualifiers") or {}
                    bad = False
                    for key, wanted in constraints.items():
                        pred, _, qual = key.partition(".")
                        if pred == predicate and q.get(qual) != wanted:
                            bad = True
                            break
                    if bad:
                        continue
                if not last and nbr in blocked:
                    continue
                nxt[nbr] |= origins
        frontier = nxt
        out.append((f"{predicate} -> {want}", len(frontier)))
        if not frontier:
            break
    reached = set()
    for origins in frontier.values():
        reached |= origins
    out.append(("compounds arriving", len(reached)))
    return out


def reach_under(g: Graph, mp: dict, blocked: set[str], virus: str) -> int:
    """Compounds reaching the virus with this exclusion set."""
    n = 0
    for drug in g.drugs():
        hits = walk_metapath(g, drug, mp["path"], mp.get("constraints"),
                             blocked)
        if virus in hits:
            n += 1
    return n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path, default=ROOT / "data/releases/v0.1")
    ap.add_argument("--channel", default="M4")
    ap.add_argument("--virus", default="NCBITaxon:2697049")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    schema = load_schema()
    if args.channel not in schema.metapaths:
        print(f"no metapath {args.channel}; have "
              f"{sorted(schema.metapaths)}")
        return 2
    mp = schema.metapaths[args.channel]
    g = Graph.load(args.release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())

    hubs = _hub_set(g, schema.retrieval_limits)
    pathways = _oversized_pathways(g, schema.retrieval_limits)
    applies = bool(mp.get("exclude_hubs"))

    print(f"release {args.release.name} | {args.channel} | {args.virus}")
    print(f"path: {' '.join(str(p) for p in mp['path'])}")
    print(f"exclude_hubs: {applies}  "
          f"({'limits apply to intermediates' if applies else 'LIMITS DO NOT APPLY -- the counterfactuals below are all equal and that is the finding'})")
    print(f"hub set: {len(hubs):,} proteins over degree "
          f"{schema.retrieval_limits.get('hub_degree_cut')}   "
          f"oversized pathways: {len(pathways):,} over "
          f"{schema.retrieval_limits.get('max_pathway_members')} members\n")

    blocked = (hubs | pathways) if applies else set()
    print("  survivors per hop, with the channel's own exclusions")
    for label, n in survivors(g, mp["path"], mp.get("constraints"), blocked,
                              list(g.drugs())):
        print(f"    {label:<42} {n:8,}")
    print()

    print("  reach under four exclusion settings")
    settings = [("as configured", blocked),
                ("hub cut relaxed", pathways if applies else set()),
                ("pathway cap relaxed", hubs if applies else set()),
                ("both relaxed", set())]
    base = None
    for label, b in settings:
        n = reach_under(g, mp, b, args.virus)
        if base is None:
            base = n
        delta = "" if n == base else f"  {n - base:+,} ({n / base - 1:+.1%})"
        print(f"    {label:<24} {n:8,}{delta}")
        if label == "both relaxed":
            key = (args.virus, args.channel)
            # THE PINNED FIGURES ARE v0.1's. Comparing a reach computed on
            # another release against them would be the stale-before-column
            # problem this repository keeps finding: a verdict that reads as
            # measured while one of its two numbers came from a different
            # graph.
            if args.release.name != MEASURED_ON:
                print(f"\n    The label_headroom figures pinned here were "
                      f"measured on {MEASURED_ON}; this is "
                      f"{args.release.name}, so no gap verdict is given.")
            elif key in MEASURED_GAP:
                _reach, needed = MEASURED_GAP[key]
                verdict = ("CLOSES the gap" if n >= needed
                           else "does NOT close the gap")
                print(f"\n    label_headroom needs {needed:,} labelled "
                      f"compounds for 80% power.")
                print(f"    Relaxing every retrieval limit gives {n:,}, "
                      f"which {verdict}.")
                if n < needed:
                    print(f"    Short by {needed - n:,}. The ceiling is not "
                          f"the parameters; it is the data.")
    print()
    print("  The counterfactuals are a COST, not a recommendation. The hub")
    print("  cut exists because a 1,389-member pathway makes every drug")
    print("  adjacent to every host factor, and raising either limit moves")
    print("  every number in data/results, so it needs a registration.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
