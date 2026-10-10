#!/usr/bin/env python3
"""Which compounds sit on both sides of the temporal cutoff, and what put them there.

WHY THIS EXISTS. temporal_split.py now refuses with

    FAIL  L5 label overlap: 5 positive compounds for NCBITaxon:2697049
          appear in BOTH train and test labels
    FAIL  L5 label overlap: 14 negative compounds for NCBITaxon:2697049
          appear in BOTH train and test labels

and the refusal is right: a compound measured against a virus both before and
after the cutoff is already known, so it is not a prospective test case.

But the last clean temporal run (2026-10-06, edges_sha 73ebb09d) reported NO
L5 failure -- only ceiling notes. So the overlap is new, introduced somewhere
between that graph and the current one (b234bc07).

The candidate is the NCATS CPE layer: one screen, one publication
(PMID:33708112, 2021), so every label it contributes carries the same year.
With cutoff 2021 and `year <= cutoff` routing to TRAIN, any compound that also
carries a post-2021 ChEMBL label lands in TEST, and L5 fires. That is a
hypothesis. This script prints the evidence for or against it instead of
acting on it, because the last two times a layer was changed on a guess about
its contents the guess was wrong.

WHAT TO READ. The per-compound rows, and the source/side table at the end:

  * train side one source at one year, test side another source later
    -> a screen's single publication date meeting the cutoff. The decision is
       which side a different-screen measurement of a known compound belongs
       on, and it is a protocol decision, not a bug.
  * both sides the same source at different years
    -> ordinary re-measurement. L5 is doing exactly its job and the compounds
       should leave the test set.

    python scripts/diagnose_temporal_leak.py
    python scripts/diagnose_temporal_leak.py --cutoff 2021 --undated computed_only

Read-only. Writes nothing.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.schema import load_schema
from kgav.temporal import build_split, edge_year

VALUE_KEYS = ("standard_type", "standard_relation", "standard_value",
              "standard_units")


def _labels(release: Path) -> dict[str, str]:
    """node id -> display label, read straight from nodes.jsonl.

    Deliberately not via baselines.Graph: that applies schema symmetry and
    drops predicates, and this script must see every label edge as written.
    """
    out: dict[str, str] = {}
    with open(release / "nodes.jsonl") as fh:
        for line in fh:
            if not line.strip():
                continue
            n = json.loads(line)
            out[n["id"]] = (n.get("properties") or {}).get("label", n["id"])
    return out


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--cutoff", type=int, default=2021)
    ap.add_argument("--undated", default="include",
                    choices=("include", "exclude", "computed_only"))
    ap.add_argument("--max-compounds", type=int, default=40,
                    help="shown per class per virus, so one pathological "
                         "virus cannot bury the rest")
    args = ap.parse_args()

    if not (args.release / "edges.jsonl").exists():
        sys.exit(f"no edges.jsonl in {args.release}")

    s = load_schema()
    label_preds = s.evaluation_label_predicates()
    split = build_split(args.release, args.cutoff, s.held_out_predicates(),
                        args.undated, label_predicates=label_preds)

    # Index label edges by (subject, object). The point of the script is the
    # source and year columns on each row, so the rows themselves are needed,
    # not just the compound names.
    by_pair: dict[tuple[str, str], list[dict]] = defaultdict(list)
    with open(args.release / "edges.jsonl") as fh:
        for line in fh:
            if not line.strip():
                continue
            e = json.loads(line)
            if e["predicate"] in label_preds:
                by_pair[(e["subject"], e["object"])].append(e)

    name = _labels(args.release)

    print(f"release {args.release}  cutoff {args.cutoff}  "
          f"undated={args.undated}")
    print(f"label predicates: {sorted(label_preds)}\n")

    pairs = 0
    src_side: dict[str, Counter] = {"train": Counter(), "test": Counter(),
                                    "undated": Counter()}

    for kind, test_side, train_side in (
            ("POSITIVE", split.test_labels, split.train_labels),
            ("NEGATIVE", split.test_negatives, split.train_negatives)):
        for virus in sorted(set(test_side) | set(train_side)):
            overlap = test_side.get(virus, set()) & train_side.get(virus, set())
            if not overlap:
                continue
            pairs += len(overlap)
            print("=" * 78)
            print(f"{kind} overlap: {len(overlap):,} compounds for "
                  f"{name.get(virus, virus)} ({virus})")
            print("=" * 78)
            for drug in sorted(overlap)[:args.max_compounds]:
                print(f"\n  {name.get(drug, drug)}  [{drug}]")
                rows = by_pair.get((drug, virus), [])
                if not rows:
                    print("    (no label edge for this pair -- the overlap "
                          "came in under a different object key)")
                for e in sorted(rows, key=lambda r: (edge_year(r) or 0,
                                                     r["predicate"])):
                    year = edge_year(e)
                    side = ("undated" if year is None
                            else "train" if year <= args.cutoff else "test")
                    src = e.get("primary_knowledge_source") or "?"
                    src_side[side][src] += 1
                    print(f"    {side:<7} {year or '-'!s:<6} "
                          f"{e['predicate']:<32} {src:<26} "
                          f"{e.get('publications') or ''}")
                    props = e.get("properties") or {}
                    detail = " ".join(f"{k}={props[k]}" for k in VALUE_KEYS
                                      if props.get(k) is not None)
                    if detail:
                        print(f"            {detail}")
            if len(overlap) > args.max_compounds:
                print(f"\n  ... {len(overlap) - args.max_compounds:,} more "
                      f"(raise --max-compounds)")
            print()

    if not pairs:
        print("no overlap: L5 would not fire on this release and policy.")
        return 0

    print("=" * 78)
    print("WHICH SOURCE PUT A LABEL ON WHICH SIDE  (overlapping compounds only)")
    print("=" * 78)
    for side in ("train", "test", "undated"):
        bound = ("<=" if side == "train" else ">" if side == "test" else "no")
        print(f"\n  {side} side (year {bound} {args.cutoff}):")
        if not src_side[side]:
            print("    none")
        for src, n in src_side[side].most_common():
            print(f"    {n:>6,}  {src}")

    print(f"\n{pairs:,} overlapping compound/virus pairs.")
    print("\nOne source at one year on the train side against a different "
          "source later on\nthe test side is a screen's publication date "
          "meeting the cutoff, not\nre-measurement. Those are different "
          "problems with different fixes, and this\nscript does not choose "
          "between them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
