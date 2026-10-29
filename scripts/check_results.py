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

from kgav.provenance import check_derived, check_dir, code_commit

MARK = {"ok": "ok      ", "stale": "STALE   ",
        "unstamped": "unstamped", "refused": "REFUSED ",
        "absent": "absent  "}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    # AUDIT THE DERIVED RELEASE TOO. temporal_split.py cuts v0.1-train2021
    # from v0.1 and hard_negatives.py reads it, so it sits in the dependency
    # chain of a results file while being invisible to an audit that only
    # walked data/results. A stale derived release produces a fresh-looking
    # result: see provenance.write_derived.
    ap.add_argument("--derived", type=Path, action="append", default=None,
                    help="a release cut from --release (repeatable). "
                         "Defaults to data/releases/v0.1-train2021.")
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

    derived = args.derived if args.derived is not None \
        else [root / "data/releases/v0.1-train2021"]

    rows = check_dir(args.results_dir, args.release)
    if not rows:
        print(f"no results files in {args.results_dir}")
        return 0

    print(f"{args.results_dir} against {args.release}\n")
    for name, status, msg in rows:
        print(f"  {MARK[status]}  {name}")
        if status != "ok":
            print(f"            {msg}")

    # Reported separately and counted separately: a derived release is not a
    # result, and folding it into the results tally would change a number
    # that docs/evaluation.md quotes.
    derived_rows = [(Path(d).name, *check_derived(d, args.release))
                    for d in derived]
    if derived_rows:
        print(f"\nderived releases cut from {args.release.name}:")
        for name, status, msg in derived_rows:
            print(f"  {MARK[status]}  {name}")
            if status != "ok":
                print(f"            {msg}")
    stale_derived = [r for r in derived_rows if r[1] in ("stale", "unstamped")]

    stale = [r for r in rows if r[1] == "stale"]
    unstamped = [r for r in rows if r[1] == "unstamped"]
    refused = [r for r in rows if r[1] == "refused"]
    print(f"\n{len(rows) - len(stale) - len(unstamped) - len(refused)} current, "
          f"{len(stale)} stale, {len(unstamped)} unstamped, "
          f"{len(refused)} refused")
    if stale:
        print("\nA stale file describes a graph that no longer exists. Comparing "
              "against it\nattributes release drift to whatever change is under "
              "test. Regenerate it.")
    if refused:
        # NOT a failure. The file is current and honest: it records that the
        # protocol declined to run on this graph, which is the audit working.
        print("\nA refused protocol is not a stale file -- it is a current, "
              "attributable record\nthat the protocol declined to run on this "
              "graph. Resolve the audit failure\nor the protocol stays "
              "refused; either way nothing here is mistakable for a result.")
    if unstamped and not args.strict:
        print("\nUnstamped files predate provenance recording. They cannot be "
              "verified either\nway; regenerate them to make them checkable. "
              "Pass --strict to fail on these.")
    if stale_derived:
        # A FAILURE even when every results file matched, because this is the
        # state that produces a fresh-looking stale number next time
        # hard_negatives.py runs. It is caught here instead of after the
        # number is written.
        print("\nA derived release that cannot be shown to match this graph "
              "is read by\nhard_negatives.py, whose output is then stamped "
              "against the CURRENT release --\nso the unverifiable input "
              "would be certified current. Re-cut it with\ntemporal_split.py. "
              "(unstamped means it predates the stamp, not that it is known "
              "wrong;\neither way nothing can attribute it to a graph.)")
    _report_release_drift(args.release)
    return 1 if (stale or stale_derived or (unstamped and args.strict)) else 0


def _report_release_drift(release) -> None:
    """Was this release assembled by the code that is checked out now?

    NOT A FAILURE, and deliberately so: working ahead of the release
    is ordinary, and most commits do not change what assembly
    produces. What is not ordinary is not KNOWING, which is the state
    this project was in until 2026-10-29.

    v0.1 was assembled on 2026-10-28 from commit 502d81e. A rebuild
    the next day, from the same raw inputs, produced a different
    edges_sha -- so the shipped graph was not reproducible from the
    code in the tree, every evaluation rested on it, and nothing could
    report that. check_results compares results against the release ON
    DISK and never against a rebuild, and before this the release
    named no commit to compare HEAD against.

    A line of output is not a guarantee of reproducibility. It is the
    one cheap thing that would have made a year-old question visible
    the day it appeared.
    """
    import json as _json

    # ONE GUARD, NOT TWO. An explicit exists() check here was dead
    # code: read_text() on a missing file raises FileNotFoundError,
    # which is an OSError and is caught below. A mutant that deleted
    # the exists() check passed, which is what a redundant branch
    # always looks like -- and an untestable branch is worse than no
    # branch, because it reads as protection.
    try:
        built = (_json.loads((release / "MANIFEST.json").read_text())
                 .get("code") or {})
    except (OSError, ValueError):
        return
    commit = built.get("commit")
    if not commit:
        print("\nThis release records no build commit, so whether the code "
              "in the tree\nstill produces it cannot be checked. "
              "Reassemble to stamp it.")
        return
    now = (code_commit() or {}).get("commit")
    if not now or commit == now:
        return
    print(f"\nThis release was assembled at {commit[:12]} and HEAD is "
          f"{now[:12]}.\nThat is normal and is not a failure. It does mean "
          f"nothing here has shown\nthat the current code still produces "
          f"this graph -- only a rebuild shows that,\nand the one time it "
          f"was tried the graph had moved.")
    if built.get("dirty"):
        print("The tree was DIRTY when it was assembled, so the commit "
              "above does not\nfully identify what built it.")


if __name__ == "__main__":
    raise SystemExit(main())
