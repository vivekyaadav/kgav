#!/usr/bin/env python3
"""Run the whole ingest chain in dependency order and stop on a real failure.

    make rebuild              # the chain
    make rebuild-plan         # what it would run, no execution

WHY NOT A LIST OF COMMANDS IN THE MAKEFILE. Three layers read the ASSEMBLED
release rather than the per-layer ones, so assemble_release has to run in the
middle and again at the end. src/kgav/pipeline.py holds that order as data and
this script executes it; neither restates the other.

WHAT COUNTS AS A FAILURE. A layer whose raw inputs are absent is SKIPPED and
named -- data/raw is gitignored and populated by hand, so an absent input is
an ordinary state of a checkout. A layer whose inputs are present and which
then exits non-zero STOPS the run, because the next assemble would otherwise
write a release missing that layer while every count still looked right.

spine and host are required: without them there is no graph, and the run
refuses to start rather than assembling an empty one.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.pipeline import blocking, plan


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--plan", action="store_true",
                    help="print the steps and exit, running nothing")
    ap.add_argument("--root", type=Path, default=ROOT)
    args = ap.parse_args(argv)

    steps = plan(args.root)
    layers = [s for s in steps if s["kind"] == "layer"]
    runnable = [s for s in layers if s["run"]]
    skipped = [s for s in layers if not s["run"]]

    print(f"{len(runnable)} layer(s) to run, {len(skipped)} skipped, "
          f"2 assemble passes\n")
    for s in steps:
        if s["kind"] == "assemble":
            print(f"  phase {s['phase']}  ASSEMBLE   {s['script']}")
        elif s["run"]:
            print(f"  phase {s['phase']}  run        {s['name']}")
        elif s["required"]:
            print(f"  phase {s['phase']}  REQUIRED   {s['name']}"
                  f"  -- {s['reason']}")
        elif s["existing"]:
            when = time.strftime("%Y-%m-%d %H:%M",
                                 time.localtime(s["written"]))
            print(f"  phase {s['phase']}  REUSE      {s['name']}"
                  f"  -- {s['reason']}; keeping the layer written {when}")
        else:
            print(f"  phase {s['phase']}  absent     {s['name']}"
                  f"  -- {s['reason']}")
    print()

    stop = blocking(args.root)
    if stop:
        print("REFUSING TO START. Required layer(s) cannot run:")
        for line in stop:
            print(f"  {line}")
        print("\nWithout these there is no graph to assemble, and "
              "assemble_release\nwould write a release anyway.")
        return 2

    # TWO KINDS OF SKIP, AND THIS SAID ONLY ONE OF THEM. assemble_release
    # builds its layer set from whatever v0.1-* directories exist, so a
    # skipped ingest whose directory is already on disk does not drop the
    # layer -- it carries the old one forward. The first version of this
    # message claimed every skipped layer was absent from the release, which
    # was wrong for exactly the case that matters most.
    reused = [s for s in skipped if s["existing"]]
    absent = [s for s in skipped if not s["existing"]]
    if absent:
        print("ABSENT from the release this run produces: "
              f"{', '.join(s['name'] for s in absent)}")
        print("  No directory on disk, so assemble_release skips them. Any "
              "metapath that\n  needs one reaches nothing, and the counts "
              "still look right.\n")
    if reused:
        print("REUSED, not rebuilt: "
              f"{', '.join(s['name'] for s in reused)}")
        print("  assemble_release reads every v0.1-* directory that exists, "
              "so these are\n  carried forward from when they were last "
              "written. That is correct only if\n  nothing they depend on "
              "has changed. If you have edited config/viruses.yaml --\n"
              "  promoted a virus, say -- a reused layer has no data for it "
              "and the build\n  will look complete anyway.\n")

    if args.plan:
        print("--plan: nothing was run.")
        return 0

    t0 = time.time()
    for s in steps:
        if s["kind"] == "layer" and not s["run"]:
            continue
        label = s["script"]
        print(f"\n{'=' * 70}\nphase {s['phase']}  {label}\n{'=' * 70}",
              flush=True)
        started = time.time()
        rc = subprocess.call([sys.executable,
                              str(args.root / "scripts" / label)])
        took = time.time() - started
        if rc != 0:
            print(f"\nFAILED after {took:.0f}s: {label} exited {rc}.")
            print("Stopping here. The release on disk is from before this "
                  "run, so it is\nconsistent -- rerunning after the fix "
                  "picks up where this stopped,\nsince every layer already "
                  "written is reused.")
            return 1
        print(f"  ok  ({took:.0f}s)")

    print(f"\n{'=' * 70}")
    print(f"rebuild complete in {time.time() - t0:.0f}s")
    if absent:
        print(f"WITHOUT: {', '.join(s['name'] for s in absent)}")
    if reused:
        print(f"REUSED (not rebuilt this run): "
              f"{', '.join(s['name'] for s in reused)}")
    print("Results in data/results were computed from the previous release "
          "and are\nnow stale. `make check-results` reports which.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
