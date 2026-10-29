#!/usr/bin/env python3
"""How large does the frontier actually get, per metapath?

The cap WAS max_paths = 20000, undocumented, and on v0.1 it truncated
114 of M4's 2,983 reaching drugs and nothing else. It is now
MAX_PATHS = 1,000,000, chosen from the output of this script. Those scores
therefore depend on the order edges.jsonl was written in.

diagnose_path_cap.py already showed the UNCAPPED walk completes, so
the cap can simply be raised. This measures by how much: the largest
frontier any single drug produces on each metapath, which is the
number a cap has to clear for the traversal to be complete.

A cap chosen from this is a documented engineering limit. The current
one is a number with no recorded reason that silently changes results.

Read-only. Writes nothing.

Usage:  python scripts/diagnose_frontier_size.py [--release DIR] [--metapath M4]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav import baselines as B            # noqa: E402
from kgav.schema import load_schema        # noqa: E402


def max_frontier(g, start, path, constraints, blocked, w=B.DEFAULT_W):
    """The largest `nxt` the walk builds, with no cap applied."""
    frontier = [(start, 1.0, frozenset({start}))]
    steps = [(B.parse_hop(path[i]), path[i + 1])
             for i in range(1, len(path) - 1, 2)]
    peak = 1
    for step_i, ((predicate, reverse), want_class) in enumerate(steps):
        nxt = []
        last = step_i == len(steps) - 1
        for node, weight, visited in frontier:
            for nbr, edge in g.neighbours(predicate, node, reverse):
                if nbr in visited or g.klass(nbr) != want_class:
                    continue
                if not B.qualifier_ok(edge, predicate, constraints):
                    continue
                if not last and nbr in blocked:
                    continue
                step_w = (g.degree.get(nbr, 1) or 1) ** -w
                nxt.append((nbr, weight * step_w, visited | {nbr}))
        peak = max(peak, len(nxt))
        frontier = nxt
        if not frontier:
            break
    return peak


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path,
                    default=ROOT / "data" / "releases" / "v0.1")
    ap.add_argument("--metapath", default=None,
                    help="measure only this one (M4 is the slow one)")
    args = ap.parse_args()

    s = load_schema()
    g = B.Graph.load(args.release, skip_predicates=s.held_out_predicates(),
                     symmetry=s.symmetry())
    drugs = list(g.drugs())
    print(f"{len(drugs):,} drugs, cap currently {B.MAX_PATHS:,}\n")
    print(f"{'metapath':<10} {'peak frontier':>14} {'over cap':>9} "
          f"{'seconds':>8}")
    worst = 0
    for name, mp in sorted(s.metapaths.items()):
        if args.metapath and name != args.metapath:
            continue
        blocked = (B._hub_set(g, s.retrieval_limits)
                   | B._oversized_pathways(g, s.retrieval_limits)) \
            if mp.get("exclude_hubs") else set()
        t0 = time.time()
        peak = 0
        over = 0
        for d in drugs:
            m = max_frontier(g, d, mp["path"], mp.get("constraints"), blocked)
            peak = max(peak, m)
            if m > B.MAX_PATHS:
                over += 1
        worst = max(worst, peak)
        print(f"{name:<10} {peak:>14,} {over:>9,} {time.time() - t0:>8.1f}")
    print(f"\nLargest frontier seen: {worst:,}")
    print("A cap above that completes every walk on this release. Pick it "
          "with headroom,\nstate the number in the docstring, and the "
          "figures stop depending on edge order.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
