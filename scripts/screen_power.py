#!/usr/bin/env python3
"""How many compounds a prospective wet-lab screen would need, per cell.

    make screen-power

THE QUESTION. "Register the prospective screen for every virus, not just
DENV-2" is the obvious generalisation, and it is answerable rather than
arguable: a two-arm screen either has the power to distinguish the
scorer's hit rate from the pool's, at a feasible number of compounds, or
it does not.

THE DESIGN BEING SIZED. Two arms of k compounds each, from the same
candidate pool, on the same plates, read blind:

    arm A   the top k by max-Tanimoto to pre-cutoff known actives
    arm B   k drawn at random from the same pool

One-sided Fisher exact on the 2x2 of hits against misses. The null is
that the two arms hit at the same rate.

WHY TWO ARMS AND NOT A THRESHOLD ON ARM A ALONE. A one-arm screen has to
be compared against an assumed base rate, and the base rate of an
unscreened pool is exactly what nobody knows -- the figure
scripts/hit_rate_at_k.py has to call an upper bound. With a random arm on
the same plates the comparator is measured in the same experiment, under
the same assay, the same cell line and the same batch of cells. It costs
twice the compounds and it is the only version whose answer does not rest
on a number supplied from outside.

WHAT IS ASSUMED, AND IT IS THE ONE SOFT INPUT. Sizing needs an expected
rate for each arm before the experiment exists. Arm A is given the
measured precision@10 under the registered k10 control -- the sterner of
the two prospective controls, so the sizing is not done on the
uncontrolled figure. Arm B is given the cell's own measured base rate.
Both come from compounds that WERE screened, and a real screen draws
from compounds that were not, so these are estimates and the power
column is optimistic. A cell that cannot be powered on optimistic inputs
cannot be powered.

THE FAMILY SIZE IS THE POINT OF THE SECOND COLUMN. Registering ten cells
makes the per-cell alpha 0.05/10, and the table shows what that costs:
DENV-2 goes from 10 compounds per arm to 16, DENV-4 from 26 to 46.
Registering cells that cannot be powered does not merely waste
compounds, it RAISES THE BAR on the one cell that can.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from temporal_structure_scorer import ALPHA as TEMPORAL_ALPHA
from temporal_structure_scorer import FAMILY as TEMPORAL_FAMILY

#: The conventional uncorrected level, named so the pair below reads as
#: a comparison rather than as two magic numbers.
UNCORRECTED_ALPHA = 0.05

#: Power the sizing targets. 0.80 is the conventional floor and is stated
#: here rather than assumed, because a screen sized at 0.5 power is a
#: coin flip dressed as an experiment.
TARGET_POWER = 0.80

#: Arm sizes considered. Beyond this a screen stops being a screen and
#: becomes the library, so a cell needing more is reported as out of
#: reach rather than given a number nobody would fund.
K_MAX = 120
K_STEP = 2

#: The alphas reported: uncorrected, and the registered family's
#: correction -- what registering every evaluable prospective cell would
#: cost.
#:
#: IMPORTED, NOT RESTATED. This read `(0.05, 0.005)` under a comment
#: saying "corrected for a family of ten", beside a summary line that
#: printed "once the family is {len(rows)}". If the number of evaluable
#: cells ever changed, the table would have sized every arm at 0.05/10
#: while the sentence beneath it named a different family -- and the arm
#: size is the one figure a signed registration cannot revise.
#: scripts/hit_rate_at_k.py already imports from this module.
ALPHAS = (UNCORRECTED_ALPHA, TEMPORAL_ALPHA)

RULE = "=" * 78


def fisher_one_sided(a: int, k_a: int, b: int, k_b: int) -> float:
    """P(arm A hits >= a) with both margins fixed. Exact, no scipy."""
    n, tot = a + b, k_a + k_b
    p = 0.0
    for x in range(a, min(k_a, n) + 1):
        if n - x > k_b:
            continue
        p += comb(k_a, x) * comb(k_b, n - x) / comb(tot, n)
    return p


def power(k: int, p_a: float, p_b: float, alpha: float) -> float:
    """P(the test rejects) with k per arm, by enumerating both arms.

    Exact rather than normal-approximate: at these k the approximation
    is wrong in the direction that flatters the design, and the whole
    point of the table is to say which cells cannot be done.
    """
    total = 0.0
    for a in range(k + 1):
        w_a = comb(k, a) * p_a ** a * (1 - p_a) ** (k - a)
        if w_a < 1e-15:
            continue
        for b in range(k + 1):
            w_b = comb(k, b) * p_b ** b * (1 - p_b) ** (k - b)
            if w_b < 1e-15:
                continue
            if fisher_one_sided(a, k, b, k) <= alpha:
                total += w_a * w_b
    return total


def smallest_arm(p_a: float, p_b: float, alpha: float) -> int | None:
    """Smallest k per arm reaching TARGET_POWER, or None past K_MAX."""
    for k in range(4, K_MAX + 1, K_STEP):
        if power(k, p_a, p_b, alpha) >= TARGET_POWER:
            return k
    return None


def measured(k: int = 10) -> dict:
    """Run scripts/hit_rate_at_k.py and read its prospective cells.

    THROUGH THE SCRIPT, not a second implementation. precision@k has one
    home in this repository and the reason to keep it that way is that
    this file would otherwise be a second opinion about what the scorer
    scored.
    """
    script = ROOT / "scripts" / "hit_rate_at_k.py"
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "hit_rate.json"
        r = subprocess.run(
            [sys.executable, str(script), "--out", str(out), "--k", str(k)],
            capture_output=True, text=True, cwd=ROOT, timeout=1800)
        if r.returncode != 0:
            # BOTH STREAMS. hit_rate_at_k.py prints its refusals --
            # "the split was built from a different release than the one
            # on disk" -- to STDOUT and returns 1, so showing stderr
            # alone produced "hit_rate_at_k.py failed:" and nothing
            # after it, in exactly the case an operator needs the reason
            # (a stale split after a rebuild).
            raise SystemExit("hit_rate_at_k.py failed:\n"
                             + (r.stdout or "") + (r.stderr or ""))
        return json.loads(out.read_text())


def _per_k(block: dict, k: int) -> dict | None:
    """per_k survives a JSON round trip with string keys.

    None when the cut depth is absent. hit_rate_at_k.precision_rows
    skips any k larger than the cell's pool, so `--k 100` against a cell
    with 81 candidates has no entry -- and this used to raise a bare
    KeyError that aborted the whole table rather than naming the cell.
    """
    pk = block["per_k"]
    if str(k) in pk:
        return pk[str(k)]
    return pk.get(k)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--k", type=int, default=10,
                    help="the cut depth whose precision sizes arm A")
    ap.add_argument("--max-arm", type=int, default=40,
                    help="the largest arm treated as fundable, for the "
                         "summary line only")
    ap.add_argument("--out", type=Path, default=None,
                    help="also write the table as JSON. Like "
                         "scripts/hit_rate_at_k.py this is NOT written to "
                         "data/results by default: it is a design "
                         "calculation, not a measurement of the graph.")
    args = ap.parse_args(argv)

    doc = measured(args.k)
    cells = doc["cells"]

    print(RULE)
    print("PROSPECTIVE SCREEN — COMPOUNDS PER ARM FOR "
          f"{TARGET_POWER:.0%} POWER")
    print(RULE)
    print()
    print(f"  split {doc['split']}, cutoff {doc['cutoff']}, "
          f"undated {doc['undated_policy']}")
    print(f"  arm A sized on precision@{args.k} under the registered k10 "
          f"control;")
    print("  arm B on the cell's own base rate. One-sided Fisher exact, "
          "k per arm.")
    print()
    print("  'out of reach' means no arm size up to "
          f"{K_MAX} reaches {TARGET_POWER:.0%}.")
    print()

    hdr = (f"  {'virus':<14}{'base':>7}{'prec':>7}{'lift':>7}"
           + "".join(f"{'k @a=' + f'{a:g}':>14}" for a in ALPHAS))
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    rows = []
    skipped = []
    for cell in cells:
        p_b = cell["base_rate"]
        sized = _per_k(cell["controls"]["symmetric_k10"], args.k)
        if sized is None:
            skipped.append(f"{cell['name']} (pool {cell['pool']:,})")
            continue
        p_a = sized["precision_expected"]
        needs = [smallest_arm(p_a, p_b, a) for a in ALPHAS]
        rows.append({"name": cell["name"], "base_rate": p_b,
                     "precision": p_a,
                     "lift": (p_a / p_b) if p_b else None,
                     "arm_per_alpha": dict(zip(map(str, ALPHAS), needs))})
        cols = "".join(f"{(str(n) if n else 'out of reach'):>14}"
                       for n in needs)
        print(f"  {cell['name']:<14}{p_b:>6.1%}{p_a:>7.0%}"
              f"{(p_a / p_b if p_b else 0):>6.1f}x{cols}")
    print("  " + "-" * (len(hdr) - 2))
    if skipped:
        print(f"  not sized at k={args.k}, pool smaller than the cut: "
              + ", ".join(skipped))
    print()

    corrected = str(ALPHAS[-1])
    fundable = [r for r in rows
                if r["arm_per_alpha"][corrected]
                and r["arm_per_alpha"][corrected] <= args.max_arm]
    reach = [r for r in rows if r["arm_per_alpha"][corrected]]
    print(RULE)
    print("WHAT THAT MEANS FOR REGISTERING EVERY VIRUS")
    print(RULE)
    print()
    print(f"  {len(reach)} of {len(rows)} cells can be powered at all, at the "
          f"corrected alpha.")
    print(f"  {len(fundable)} of {len(rows)} need {args.max_arm} compounds "
          f"per arm or fewer: "
          + (", ".join(r["name"] for r in fundable) or "none"))
    total = sum(r["arm_per_alpha"][corrected] * 2 for r in reach)
    print(f"  Registering every cell that CAN be powered: "
          f"{total:,} compounds across {len(reach)} cells.")
    print()
    print("  THE LIMIT IS NOT THE AUC, IT IS ENRICHMENT HEADROOM. A cell")
    print("  whose pool is already half actives cannot show a scorer")
    print("  working: 60% against a 47% base is a lift of 1.3, and")
    print("  separating 1.3 from 1.0 takes hundreds of compounds whatever")
    print("  the AUC says. The cells that are out of reach are out of reach")
    print("  for two different reasons -- no enrichment, or no headroom --")
    print("  and the lift column tells them apart.")
    print()
    print("  AND REGISTERING THE UNPOWERABLE CELLS COSTS THE POWERABLE ONE.")
    for r in rows:
        a_un = r["arm_per_alpha"][str(ALPHAS[0])]
        a_co = r["arm_per_alpha"][corrected]
        if a_un and a_co and a_co > a_un:
            print(f"    {r['name']:<14}{a_un} -> {a_co} compounds per arm "
                  f"at alpha {ALPHAS[-1]:g} (a family of "
                  f"{TEMPORAL_FAMILY})")
    print()

    if args.out:
        args.out.write_text(json.dumps({
            "computed_by": "scripts/screen_power.py",
            "design": "two arms of k, one-sided Fisher exact, hits vs misses",
            "target_power": TARGET_POWER,
            "k_max": K_MAX,
            "k_step": K_STEP,
            "alphas": list(ALPHAS),
            "sizing_cut": args.k,
            "split": doc["split"],
            "cutoff": doc["cutoff"],
            "cells": rows,
        }, indent=2) + "\n")
        print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
