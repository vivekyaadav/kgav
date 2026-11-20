#!/usr/bin/env python3
"""Is the prospective result prediction, or retrieval of a known series?

    make prospective-analogy

THE QUESTION, AND WHY IT IS THE ONE THAT MATTERS. Six of ten prospective
cells clear: trained on activity published at or before the cutoff, they
rank post-cutoff actives above post-cutoff inactives. That is the
strongest positive result in this project.

It admits two readings, and they license different things:

  ANALOGY     the post-cutoff actives are the next compounds in a series
              the pre-cutoff set already contains, and similarity finds
              them because they are near-identical to something known.
              Useful -- it is how medicinal chemistry proceeds -- and a
              far narrower claim than the AUC sounds.
  CHEMOTYPE   the signal survives when the nearest known actives are
              withheld, so it is about a broader resemblance than one
              series.

docs/prereg_temporal_structure_scorer.md left this open in as many
words: H3's analogue question is "a pointer, not a finding". The DENV-4
amendment of 2026-11-21 answered it for ONE cell, where median 100% of
the neighbours its cross-sectional control strips are its own era-mates.
This asks the same question of every cell.

THE MEASUREMENT. A ladder: the registered prospective AUC recomputed
with the k nearest pre-cutoff references withheld from EVERY candidate
alike, for k = 0, 5, 10, 20, 40. k = 0 is the registration's primary and
k = 10 its secondary, so two rungs are the published figures and the
rest extend them.

  - a cell whose AUC collapses by k = 10 was retrieving a near neighbour
  - a cell still clearing at k = 40 is not ranking by one series

Beside it, how close the nearest pre-cutoff active actually is: the
median max-Tanimoto for positives and for negatives. A cell where the
positives sit at 0.9 is reading identity; one where they sit at 0.4 is
not.

NOTHING HERE IS A NEW CLAIM. Every rung is the registered statistic on
the registered split, recomputed at a different exclusion strength, and
the verdict word is the registration's own corrected alpha. It is a
description of a result that already exists. Not written to
data/results.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.labels import evaluate_against_negatives
from kgav.pipeline import release_dir, train_dir
from kgav.similarity import fingerprints

from structure_scorer import corrected, score_cell, small_molecule_props
from temporal_structure_scorer import CUTOFF, Z_CORRECTED, cells, evaluable

#: The ladder. 0 is the registration's primary and 10 its secondary, so
#: two rungs reproduce published figures and the rest extend them.
KS = (0, 5, 10, 20, 40)

RULE = "=" * 84


def nearest_similarity(cand: list[str], refs: list[str], fps: dict) -> dict:
    """max-Tanimoto of each candidate to the reference set, no exclusion."""
    from rdkit import DataStructs

    ref_fps = [fps[r] for r in refs if r in fps]
    out = {}
    for c in cand:
        if c not in fps or not ref_fps:
            continue
        out[c] = max(DataStructs.BulkTanimotoSimilarity(fps[c], ref_fps))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--split", type=Path, default=train_dir(ROOT, CUTOFF))
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    labels = json.loads((args.split / "TEST_LABELS.json").read_text())
    ev = evaluable(cells(labels))
    props = small_molecule_props(args.release)
    everyone = sorted({d for v in ev.values() for part in v for d in part})
    fps, _ = fingerprints(
        {i: (props.get(i, {}).get("smiles") or "") for i in everyone})

    temporal = json.loads(
        (ROOT / "data" / "results" / "temporal_structure_scorer.json")
        .read_text())
    names = {t: c["name"] for t, c in (temporal.get("cells") or {}).items()}
    registered_clears = set(temporal.get("clears") or [])

    print(RULE)
    print("IS THE PROSPECTIVE SIGNAL ANALOGY? — AUC as nearest references "
          "are withheld")
    print(RULE)
    print()
    print(f"  split {args.split.name}, cutoff {temporal.get('cutoff')}, "
          f"corrected z {Z_CORRECTED} (family of {temporal.get('family_size')})")
    print("  k = how many of each candidate's nearest PRE-cutoff actives are "
          "withheld.")
    print("  k=0 is the registration's primary, k=10 its secondary. "
          "* = still clears.")
    print()
    head = (f"  {'virus':<14}" + "".join(f"{'k=' + str(k):>10}" for k in KS)
            + f"{'median Tc':>12}{'  survives':>12}")
    print(head)
    print("  " + "-" * (len(head) - 2))

    rows = []
    for taxon, (refs, pos, neg) in sorted(ev.items()):
        if taxon not in registered_clears:
            continue
        cand = list(pos) + list(neg)
        ladder, clears_at = {}, []
        for k in KS:
            sc, _ = score_cell(cand, refs, fps, props, {},
                               publication_disjoint=False, knn=k)
            m = corrected(evaluate_against_negatives(sc, set(pos), set(neg)),
                          z=Z_CORRECTED)
            ok = bool(m.get("clears_corrected_alpha")
                      and m.get("direction") != "below")
            ladder[k] = {"auc": m["auc"], "clears": ok,
                         "ci_lo": m.get("auc_ci_lo_corrected"),
                         "ci_hi": m.get("auc_ci_hi_corrected")}
            if ok:
                clears_at.append(k)
        sim = nearest_similarity(cand, refs, fps)
        tc_pos = statistics.median([sim[c] for c in pos if c in sim] or [0])
        tc_neg = statistics.median([sim[c] for c in neg if c in sim] or [0])
        last = max(clears_at) if clears_at else None
        rows.append({"taxon": taxon, "name": names.get(taxon, taxon),
                     "ladder": {str(k): v for k, v in ladder.items()},
                     "median_tc_positive": tc_pos,
                     "median_tc_negative": tc_neg,
                     "clears_up_to_k": last})
        cellstr = "".join(
            f"{ladder[k]['auc']:>9.3f}" + ("*" if ladder[k]["clears"] else " ")
            for k in KS)
        survives = "k<=" + str(last) if last is not None else "nowhere"
        print(f"  {names.get(taxon, taxon):<14}{cellstr}"
              f"{tc_pos:>12.3f}{survives:>12}")
    print("  " + "-" * (len(head) - 2))
    print("  median Tc = a POSITIVE's closest pre-cutoff active. The "
          "negatives' median is printed below.")
    print()
    for r in rows:
        print(f"    {r['name']:<14}positives {r['median_tc_positive']:.3f}   "
              f"negatives {r['median_tc_negative']:.3f}   "
              f"gap {r['median_tc_positive'] - r['median_tc_negative']:+.3f}")
    print()

    print(RULE)
    print("READING IT")
    print(RULE)
    print()
    deep = [r for r in rows if (r["clears_up_to_k"] or -1) >= 20]
    shallow = [r for r in rows if (r["clears_up_to_k"] or -1) < 10]
    print(f"  {len(rows)} cells clear in the registered run.")
    print(f"  {len(deep)} still clear with 20 or more nearest references "
          f"withheld: "
          + (", ".join(r["name"] for r in deep) or "none"))
    print(f"  {len(shallow)} stop clearing before k = 10: "
          + (", ".join(r["name"] for r in shallow) or "none"))
    print()
    print("  A cell that stops clearing at small k was ranking by one near")
    print("  neighbour. A cell still clearing at k = 40 is not: forty of its")
    print("  closest known actives were removed and the order survived.")
    print()
    print("  THIS DOES NOT RETRACT ANYTHING. Every rung is the registered")
    print("  statistic on the registered split at a different exclusion")
    print("  strength, and k=0 and k=10 reproduce the published columns.")
    print()

    if args.out:
        args.out.write_text(json.dumps(
            {"computed_by": "scripts/diagnose_prospective_analogy.py",
             "split": args.split.name, "cutoff": temporal.get("cutoff"),
             "z_corrected": Z_CORRECTED, "ks": list(KS),
             "cells": rows}, indent=2) + "\n")
        print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
