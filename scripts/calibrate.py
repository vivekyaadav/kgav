#!/usr/bin/env python3
"""What each reasoning channel is measured to be worth, and whether it earns
a place in a ranking.

    python scripts/calibrate.py
    python scripts/calibrate.py --virus NCBITaxon:694009

Reads data/results/hard_negatives.json -- so it describes the LAST evaluation
run, and the protocol it names is that run's, not a default. Re-run
hard_negatives.py before trusting it.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav import calibration as C

# hard_negatives keys its results by DISPLAY LABEL, not CURIE. load()
# accepts either; this default matches what the file holds.
DEFAULT_VIRUS = "SARS-CoV-2"


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path,
                    default=root / "data/results/hard_negatives.json")
    ap.add_argument("--virus", default=DEFAULT_VIRUS,
                    help="display label or taxon CURIE")
    ap.add_argument("--section", default="cross_sectional",
                    choices=["cross_sectional", "temporal",
                             "cross-sectional", "post-2021"])
    args = ap.parse_args()

    if not args.results.exists():
        print(f"not found: {args.results}\n"
              f"Run scripts/hard_negatives.py first -- this file describes a "
              f"specific\nevaluation run and there is no default to fall back "
              f"on.")
        return 1

    cal = C.load(args.results, args.virus, section=args.section)
    if not cal.channels and not cal.refused_reason:
        # SAY WHAT IS THERE. The first version printed only that it found
        # nothing, which is indistinguishable from an empty evaluation and
        # sent the reader looking at the wrong thing.
        have = C.available_viruses(args.results, args.section)
        print(f"no channels recorded for {args.virus!r} in "
              f"{args.section} of {args.results.name}")
        print(f"  that file holds: {', '.join(have) if have else '(nothing)'}")
        print(f"  protocol it records: {cal.protocol}")
        return 1

    for line in cal.report():
        print(line)

    print()
    if cal.can_rank():
        print("A calibrated ranking is possible and will use only the channels "
              "above.")
        print("  python scripts/run_baselines.py   # with --calibrated once wired")
    else:
        print("RANKING IS REFUSED.")
        print(f"  {cal.why_not()}")
        print()
        print("  baselines.combine(..., calibration=cal) returns {} rather "
              "than an order,\n  because max-over-channels is biased upward "
              "by the channel count: a drug\n  at the top of a chance-level "
              "pool scores 1.0 and outranks one with a\n  real route. "
              "Reporting evidence for a named compound still works.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
