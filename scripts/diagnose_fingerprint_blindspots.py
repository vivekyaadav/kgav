#!/usr/bin/env python3
"""What the ECFP4 fingerprint cannot tell apart, and what that is worth.

FOUND 2026-11-12 BY READING STRUCTURES, not by a test -- the same way L27 was
found by reading names. `scripts/repurposing_candidates.py` puts 14 of its 210
suggestions at Tanimoto **1.0**, and all 14 have a DIFFERENT InChIKey
skeleton from the active they resemble. So they are not salts (L27 collapsed
those) and not stereoisomers (the skeleton exclusion removes those). They are
distinct molecules with an identical ECFP4 bit vector:

    CHEMBL4639263 / CHEMBL4526128   cycloheptyl vs cyclohexyl -- one CH2
                                    in a SATURATED RING, where every atom's
                                    radius-2 environment is CH2 flanked by
                                    CH2, so a 6- and a 7-ring fold to the
                                    same bits
    CHEMBL4857345 / CHEMBL1821856   one extra leucine in a peptide
    CHEMBL468049  / CHEMBL522355    Arg-Arg-Lys vs Arg-Lys-Lys
    CHEMBL5283054 / CHEMBL5403511   728- vs 747-character peptides

WHY IT MATTERS FOR THE AUC AND NOT ONLY FOR THE READ-OUT. The registered
control removes a candidate's self, its skeleton group, and either its
publication-sharing or its k nearest references. NONE of those catches a ring
homolog in a different paper: different skeleton, different PMID, Tanimoto
1.0. So a positive can take a perfect score from a compound that is one
methylene away, and the scorer is credited with a retrieval that is closer to
an identity.

WHAT THIS SCRIPT DOES AND DOES NOT DO. It reuses `structure_scorer.score_cell`
unchanged -- the registered function, called with the registered arguments --
and reports the score distribution, how much of each class sits at or above a
near-identity threshold, and what the AUC would be with those candidates
withheld. **That last figure is an EXPOSURE BOUND AND NOT A RESULT.** Removing
candidates after seeing the scores is a protocol change, and this project's
rule is that a protocol change is registered before it runs. The number is
here so the owner can see what is at stake before deciding whether to register
it; quoting it as the scorer's performance would be the thing three
registrations exist to prevent.

    python scripts/diagnose_fingerprint_blindspots.py
    python scripts/diagnose_fingerprint_blindspots.py --virus NCBITaxon:11060
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

#: At or above this, two fingerprints say "the same molecule" and the
#: structures above show they need not be. Reported at several thresholds
#: rather than one, because the choice of threshold is the whole question and
#: a single number would hide it.
THRESHOLDS = (1.0, 0.99, 0.95, 0.90)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=None)
    ap.add_argument("--virus", default="NCBITaxon:2697049")
    ap.add_argument("--control", choices=("publication", "knn", "strength"),
                    default="strength",
                    help="which registered control to score under; the "
                         "default is the strongest one that ran")
    args = ap.parse_args(argv)

    import structure_scorer as S
    from kgav.baselines import Graph
    from kgav.labels import build_labels, evaluate_against_negatives
    from kgav.pipeline import release_dir
    from kgav.schema import load_schema
    from kgav.baselines import label_publication_index
    from kgav.similarity import fingerprints

    release = args.release or release_dir(ROOT)
    if not (release / "nodes.jsonl").exists():
        print(f"no release at {release}")
        return 2

    prior_path = ROOT / "data/results/structure_scorer.json"
    if not prior_path.exists():
        print("REFUSED: the publication run is needed for k_cell")
        return 2
    prior = json.loads(prior_path.read_text())
    k_map = S.k_cells_from_prior(prior)

    schema = load_schema()
    g = Graph.load(release, symmetry=schema.symmetry())
    props = S.small_molecule_props(release)
    pubs = {k: frozenset(v) for k, v in label_publication_index(
        release, {S.ACTIVITY_PREDICATE}).items()}
    labels = build_labels(release, sources={S.PRIMARY_SOURCE})
    pos = labels.positives.get(args.virus, set())
    neg = labels.negatives.get(args.virus, set())
    if not pos or not neg:
        print(f"{args.virus} has no two-class population in "
              f"{S.PRIMARY_SOURCE}")
        return 2
    labelled = pos | neg
    fps, _ = fingerprints(
        {i: (props.get(i, {}).get("smiles") or "") for i in sorted(labelled)})

    knn = {"publication": 0, "knn": S.K_PRIMARY,
           "strength": k_map.get(args.virus, S.K_PRIMARY)}[args.control]
    disjoint = args.control == "publication"
    scores, counts = S.score_cell(
        sorted(labelled), sorted(pos), fps, props, pubs, top_k=1,
        publication_disjoint=disjoint, positives=frozenset(pos), knn=knn)

    name = (g.nodes.get(args.virus, {}).get("properties", {})
            or {}).get("label", args.virus)
    base = S.corrected(evaluate_against_negatives(scores, pos, neg))
    print(f"{name} under the {args.control} control"
          f"{f' (k_cell={knn})' if knn else ''}")
    print(f"  {len(pos):,} positives, {len(neg):,} negatives, "
          f"AUC {base['auc']:.4f} "
          f"[{base['auc_ci_lo_corrected']:.4f}, "
          f"{base['auc_ci_hi_corrected']:.4f}]")
    print("  this reproduces the registered run: it calls score_cell with "
          "the registered arguments\n")

    print(f"  {'threshold':>9}  {'pos at/above':>13}  {'neg at/above':>13}  "
          f"{'AUC withheld':>13}  {'delta':>8}")
    for t in THRESHOLDS:
        hot_pos = {c for c in pos if scores.get(c, 0.0) >= t}
        hot_neg = {c for c in neg if scores.get(c, 0.0) >= t}
        keep_pos, keep_neg = pos - hot_pos, neg - hot_neg
        if len(keep_pos) < S.MIN_CLASS or len(keep_neg) < S.MIN_CLASS:
            line = "    (a class empties)"
            print(f"  {t:>9.2f}  {len(hot_pos):>6,} "
                  f"({len(hot_pos)/len(pos):>5.1%})  {len(hot_neg):>6,} "
                  f"({len(hot_neg)/len(neg):>5.1%}){line}")
            continue
        sub = {c: s for c, s in scores.items()
               if c in keep_pos or c in keep_neg}
        m = S.corrected(evaluate_against_negatives(sub, keep_pos, keep_neg))
        print(f"  {t:>9.2f}  {len(hot_pos):>6,} "
              f"({len(hot_pos)/len(pos):>5.1%})  {len(hot_neg):>6,} "
              f"({len(hot_neg)/len(neg):>5.1%})  {m['auc']:>13.4f}  "
              f"{m['auc'] - base['auc']:>+8.4f}")

    print("\n  THE 'AUC withheld' COLUMN IS AN EXPOSURE BOUND, NOT A RESULT.")
    print("  Withholding candidates chosen by their scores is a protocol")
    print("  change, and this project registers a protocol before it runs.")
    print("  It says how much of the margin rests on near-identities; it")
    print("  does not say what the scorer scores.")
    print(f"\n  for reference, the registered control already removed: "
          f"self {counts['excluded_self']:,}, "
          f"skeleton {counts['excluded_skeleton']:,}, "
          f"knn {counts['excluded_knn']:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
