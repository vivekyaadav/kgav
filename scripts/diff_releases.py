#!/usr/bin/env python3
"""Are two releases the same GRAPH, or the same graph in a different ORDER?

Used on 2026-10-29 to establish that four builds of v0.1 differed by
8 to 10 CHEMICALLY_SIMILAR_TO edges rather than by byte order, which
is what edges_sha alone could not say. Kept because the question
recurs after every change to a layer.

release_fingerprint hashes the raw bytes of edges.jsonl, so edges_sha
is order-sensitive. Four builds of this project produced four
edges_sha values -- a9a68fe6..., 2f51efb4..., 3152d60b..., 7587762b...
-- and nothing so far distinguishes "the graph changed" from "the same
edges were written in a different order". Those need completely
different fixes, and the sha cannot tell them apart.

This compares the edge MULTISETS, canonicalised with sorted keys, and
reports which of the two it is.

Usage:
    python scripts/diff_releases.py <release_a> <release_b> [--ignore-source]

--ignore-source drops primary_knowledge_source before comparing, which
isolates the merge-source dedup from everything else: two releases that
differ ONLY by that field are the same graph with different labelling.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path


def edges(p: Path, ignore_source: bool) -> Counter:
    c: Counter = Counter()
    for line in (p / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if ignore_source:
            e.pop("primary_knowledge_source", None)
        c[json.dumps(e, sort_keys=True)] += 1
    return c


def main(argv: list[str]) -> int:
    ignore = "--ignore-source" in argv
    args = [a for a in argv if not a.startswith("--")]
    if len(args) != 2:
        print(__doc__)
        return 2
    a, b = Path(args[0]), Path(args[1])
    for d in (a, b):
        if not (d / "edges.jsonl").exists():
            print(f"no edges.jsonl in {d}")
            return 2

    ba = (a / "edges.jsonl").read_bytes()
    bb = (b / "edges.jsonl").read_bytes()
    print(f"bytes            {'IDENTICAL' if ba == bb else 'differ'}"
          f"   ({len(ba):,} vs {len(bb):,})")
    if ignore:
        print("comparing with primary_knowledge_source REMOVED")

    ca, cb = edges(a, ignore), edges(b, ignore)
    na, nb = sum(ca.values()), sum(cb.values())
    print(f"edges            {na:,} vs {nb:,}")

    if ca == cb:
        if ba == bb:
            print("\nIDENTICAL FILES.")
            print("Same edges, same bytes. Two builds that agree here are"
                  "\nreproducible: the release has one byte form for one"
                  "\ngraph, so edges_sha means what it is read as meaning.")
            return 0
        # NOT reached once the release is written canonically. Kept
        # because a future layer can reintroduce it, and because this
        # paragraph claiming a byte-order defect while the files were
        # IDENTICAL was printed unconditionally in the first version
        # -- a sentence that read as a finding when nothing was found.
        print("\nSAME GRAPH, DIFFERENT ORDER.")
        print("Every edge appears the same number of times in both. The"
              "\nedges_sha difference is byte order, not content: the"
              "\nfingerprint hashes the file, and the file is not written in"
              "\na defined order. That is a FINGERPRINTING defect, not a"
              "\ngraph defect, and the fix is to sort before hashing or to"
              "\nsort before writing.")
        return 0

    only_a = ca - cb
    only_b = cb - ca
    print(f"\nDIFFERENT GRAPHS.")
    print(f"  {sum(only_a.values()):,} edge(s) only in {a.name}")
    print(f"  {sum(only_b.values()):,} edge(s) only in {b.name}")

    def preds(c: Counter) -> Counter:
        out: Counter = Counter()
        for s, n in c.items():
            out[json.loads(s).get("predicate", "?")] += n
        return out

    for label, c in ((f"only in {a.name}", only_a), (f"only in {b.name}",
                                                     only_b)):
        if not c:
            continue
        print(f"\n  BY PREDICATE, {label}")
        for p, n in preds(c).most_common(10):
            print(f"      {n:>8,}  {p}")
    for label, c in ((a.name, only_a), (b.name, only_b)):
        if c:
            print(f"\n  one example only in {label}:")
            print("      " + next(iter(c))[:300])
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
