#!/usr/bin/env python3
"""Flag results files that no longer describe the release on disk.

    python scripts/check_results.py data/results
    python scripts/check_results.py data/results --release data/releases/v0.1
    python scripts/check_results.py data/results --strict   # unstamped also fails

WHY. data/results/metapath_comparison.json was generated one day before the
selectivity layer added 1,568 activity edges, and was still being read as a
baseline afterwards. The numbers looked plausible, the file looked current,
and the drift was indistinguishable from a real change: a positive count that
had actually FALLEN 1,011 -> 889 read as a RISE from 736. Nothing in the file
said which graph it described, so nothing could catch it.

Exit codes: 0 all matched, 1 at least one stale (or unstamped under --strict).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.provenance import check_dir

MARK = {"ok": "ok      ", "stale": "STALE   ", "unstamped": "unstamped"}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--strict", action="store_true",
                    help="treat an unstamped file as a failure; a file with no "
                         "provenance is the case that went undetected")
    args = ap.parse_args()

    if not args.results_dir.exists():
        print(f"no such directory: {args.results_dir}", file=sys.stderr)
        return 2
    if not (args.release / "edges.jsonl").exists():
        print(f"no release at {args.release}", file=sys.stderr)
        return 2

    rows = check_dir(args.results_dir, args.release)
    if not rows:
        print(f"no results files in {args.results_dir}")
        return 0

    print(f"{args.results_dir} against {args.release}\n")
    for name, status, msg in rows:
        print(f"  {MARK[status]}  {name}")
        if status != "ok":
            print(f"            {msg}")

    stale = [r for r in rows if r[1] == "stale"]
    unstamped = [r for r in rows if r[1] == "unstamped"]
    print(f"\n{len(rows) - len(stale) - len(unstamped)} current, "
          f"{len(stale)} stale, {len(unstamped)} unstamped")
    if stale:
        print("\nA stale file describes a graph that no longer exists. Comparing "
              "against it\nattributes release drift to whatever change is under "
              "test. Regenerate it.")
    if unstamped and not args.strict:
        print("\nUnstamped files predate provenance recording. They cannot be "
              "verified either\nway; regenerate them to make them checkable. "
              "Pass --strict to fail on these.")
    return 1 if (stale or (unstamped and args.strict)) else 0


if __name__ == "__main__":
    raise SystemExit(main())
