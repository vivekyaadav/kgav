#!/usr/bin/env python3
"""Which figures in docs/evaluation.md still appear in the results files.

    python scripts/diagnose_doc_figures.py
    python scripts/diagnose_doc_figures.py --doc docs/evaluation.md

ADVISORY, NOT A GATE, and deliberately not a test. The doc records what was
measured AND what it used to say -- "the superseded run read M1 at lift
2.22", "AUC 0.826" as the figure that was wrong. A check demanding every
quoted number be current would fail on the history the doc keeps ON PURPOSE,
which is the most valuable thing in it.

So this reports, and the reading is yours:

    a figure found in a results file   describes the graph on disk
    a figure found nowhere             is EITHER deliberately historical --
                                       and the surrounding text almost
                                       always says so -- OR stale, and
                                       nothing distinguishes them but
                                       reading the sentence

WHY BOTHER. web/client/src/pages/Home.tsx quoted five "real figures from
release v0.1" and four were wrong, in a file whose comment already recorded
being wrong once. evaluation.md is the same kind of artifact with far more
numbers in it, and nothing has ever checked them.

WHAT IT LOOKS FOR. Interval-shaped figures -- 0.512 [0.49, 0.54] -- because
an AUC with its confidence bounds is distinctive enough to match without
guessing, and it is the shape the project's claims are made in. Bare lifts
and p-values are not matched: too many collide by coincidence to mean
anything.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

#: 0.512 [0.49, 0.54] and 0.516 [0.455, 0.577]. Bounds may be 2 or 3 places.
INTERVAL = re.compile(r"([0-9]\.[0-9]{3}) \[([0-9]\.[0-9]{2,3}), ?([0-9]\.[0-9]{2,3})\]")


def _numbers(obj, out: set[str]) -> None:
    """Every float in a results file, at both printed precisions."""
    if isinstance(obj, dict):
        for v in obj.values():
            _numbers(v, out)
    elif isinstance(obj, list):
        for v in obj:
            _numbers(v, out)
    elif isinstance(obj, float):
        out.add(f"{obj:.3f}")
        out.add(f"{obj:.2f}")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--doc", type=Path, default=root / "docs" / "evaluation.md")
    ap.add_argument("--results", type=Path, default=root / "data" / "results")
    args = ap.parse_args()

    if not args.doc.exists():
        print(f"no document at {args.doc}", file=sys.stderr)
        return 2
    files = sorted(args.results.glob("*.json"))
    if not files:
        print(f"no results in {args.results}", file=sys.stderr)
        return 2

    per_file: dict[str, set[str]] = {}
    for p in files:
        try:
            doc = json.loads(p.read_text())
        except json.JSONDecodeError:                 # pragma: no cover
            continue
        seen: set[str] = set()
        _numbers(doc, seen)
        per_file[p.name] = seen

    text = args.doc.read_text()
    lines = text.splitlines()
    found, missing = [], []
    for m in INTERVAL.finditer(text):
        point, lo, hi = m.groups()
        line_no = text[:m.start()].count("\n") + 1
        # All three parts must appear in ONE file. A point estimate alone
        # collides constantly -- 0.504 is in half the files -- and matching
        # it would report everything as current.
        where = [n for n, nums in per_file.items()
                 if {point, lo, hi} <= nums]
        row = (m.group(0), line_no, lines[line_no - 1].strip()[:72], where)
        (found if where else missing).append(row)

    total = len(found) + len(missing)
    if not total:
        # A REGEX THAT STOPS MATCHING PRINTS A CLEAN BILL OF HEALTH. "0
        # figures, 0 present, 0 not" reads like success and means the
        # opposite -- the same shape as a counter nobody prints and a
        # breakdown with no print site, both of which I shipped tonight.
        print(f"NO interval-shaped figures matched in {args.doc.name}.\n"
              f"  Either the document genuinely quotes none, or the pattern\n"
              f"  {INTERVAL.pattern}\n  no longer matches how they are "
              f"written. Check before reading this as clean.", file=sys.stderr)
        return 2
    print(f"{args.doc.name}: {total} interval figures, "
          f"{len(found)} present in a results file, {len(missing)} not\n")

    if found:
        print("PRESENT")
        for fig, n, _ctx, where in found:
            print(f"  L{n:<5} {fig:<24} {', '.join(where)}")

    if missing:
        print("\nNOT FOUND IN ANY RESULTS FILE")
        for fig, n, ctx, _w in missing:
            print(f"  L{n:<5} {fig:<24} {ctx}")
        print("""
  Read each line before changing anything. This document records
  superseded figures deliberately -- the value that was wrong, and why --
  and those SHOULD be absent from the current results. A figure here is a
  finding only if the sentence around it presents it as what the graph
  says now.

  KNOWN AND EXPECTED on evaluation.md: the four protocol rows in "The
  sequence" (selective-only, all, verified-only, +publication-disjoint).
  hard_negatives.json holds ONE protocol at a time and currently holds
  row 5, so rows 1-4 cannot be present however current the document is.
  That table is a record of five separate runs, not of one file.

  A PRE-REGISTRATION IS HISTORICAL BY DEFINITION. Run this against
  docs/prereg_*.md and its figures will read absent, correctly: a prereg
  states what was known when the prediction was declared, and editing
  one to match a later release destroys the only thing it is for. Fix
  the live document; leave the prereg alone.""")

    print("""
WHAT THIS DOES NOT DECIDE
  Nothing here edits the document or judges a sentence. It matches
  interval-shaped figures against the numbers in data/results and reports
  which are which; whether an absent figure is history or rot is a reading
  question, and the doc's own prose is the evidence for it.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
