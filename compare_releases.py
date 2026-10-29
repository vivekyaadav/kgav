#!/usr/bin/env python3
"""Two releases, edge by edge: what the dedup actually changed.

Answers the three questions left open when the pre-dedup release was
overwritten:

  1. WHICH compound left SARS-CoV's pool. The pool is built from
     HAS_ANTIVIRAL_ACTIVITY_AGAINST, the one predicate whose fact
     identity includes the raw source string, so a (compound, virus)
     pair holding two label edges before and one after is the answer.
  2. WHETHER the 542 repeated-curator edges were label edges, which is
     the mechanism claim.
  3. WHETHER temporal.py's `in COMPUTED_SOURCES` was actually firing,
     i.e. whether any edge carried "kgav-computed" twice in one value.

Usage:  python compare_releases.py <before_release_dir> <after_release_dir>
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

LABEL = "HAS_ANTIVIRAL_ACTIVITY_AGAINST"


def load(p: Path):
    for line in (p / "edges.jsonl").read_text().splitlines():
        if line.strip():
            yield json.loads(line)


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    before, after = Path(argv[0]), Path(argv[1])
    for d in (before, after):
        if not (d / "edges.jsonl").exists():
            print(f"no edges.jsonl in {d}")
            return 2

    # THE FINGERPRINTS FIRST. If the rebuilt "before" is not the
    # release that was overwritten, every difference below is a
    # difference between two graphs neither of which is the one in
    # question. The expected pre-dedup edges_sha is a9a68fe6a41380c6
    # and nodes_sha 73638934bacd7827; nodes_sha must match the current
    # release too, since the dedup touched edges only.
    try:
        sys.path.insert(0, "src")
        from kgav.provenance import release_fingerprint
        for tag, d in (("before", before), ("after", after)):
            fp = release_fingerprint(d)
            print(f"{tag:<7} nodes_sha {fp['nodes_sha']}  "
                  f"edges_sha {fp['edges_sha']}")
        print()
    except Exception as e:                       # noqa: BLE001
        print(f"(could not fingerprint: {e})\n")

    counts = {}
    repeated_by_pred: dict[str, Counter] = {}
    computed_twice = {}
    pairs = {}

    for tag, d in (("before", before), ("after", after)):
        n = 0
        rep: Counter = Counter()
        ctwice = 0
        pair: dict[tuple, int] = defaultdict(int)
        for e in load(d):
            n += 1
            s = str(e.get("primary_knowledge_source") or "")
            parts = [q for q in s.split("|") if q]
            if len(parts) != len(set(parts)):
                rep[e["predicate"]] += 1
                if parts.count("infores:kgav-computed") > 1:
                    ctwice += 1
            if e["predicate"] == LABEL:
                pair[(e["subject"], e["object"])] += 1
        counts[tag] = n
        repeated_by_pred[tag] = rep
        computed_twice[tag] = ctwice
        pairs[tag] = pair

    print(f"edges   before {counts['before']:,}   after {counts['after']:,}")

    print("\n--- Q2: repeated-curator edges, BY PREDICATE ---")
    for tag in ("before", "after"):
        r = repeated_by_pred[tag]
        print(f"  {tag:<7} {sum(r.values()):,} total")
        for p, c in r.most_common():
            print(f"            {c:>7,}  {p}")
    if not repeated_by_pred["before"]:
        print("  !! none BEFORE either -- this is not the pre-dedup release")

    print("\n--- Q3: edges naming kgav-computed TWICE in one value ---")
    print(f"  before {computed_twice['before']:,}   "
          f"after {computed_twice['after']:,}")
    print("  before > 0 means temporal.py's `in COMPUTED_SOURCES` WAS "
          "misclassifying real edges;\n  before == 0 means the bug was "
          "latent all along and the v2 message is right.")

    print(f"\n--- Q1: (compound, virus) pairs losing a {LABEL} edge ---")
    b, a = pairs["before"], pairs["after"]
    lost = [(k, b[k], a.get(k, 0)) for k in b if a.get(k, 0) < b[k]]
    gained = [(k, b.get(k, 0), a[k]) for k in a if b.get(k, 0) < a[k]]
    print(f"  {len(lost)} pair(s) lost an edge, {len(gained)} gained one")
    for (subj, obj), nb, na in sorted(lost)[:40]:
        print(f"      {subj}  ->  {obj}   {nb} -> {na}")
    if len(lost) > 40:
        print(f"      ... and {len(lost) - 40} more")
    if not lost and not gained:
        print("      NONE. The pool change is not explained by label-edge "
              "merging, and the source_identity mechanism is wrong.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
