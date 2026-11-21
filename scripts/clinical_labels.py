#!/usr/bin/env python3
"""The registered clinical-label analysis: doc2/registrations/prereg_clinical_labels.md.

ONE CELL. COVID-19 against SARS-CoV-2, family of one, alpha = 0.05
two-sided, z = 1.9600. The registration fixed all of that on 2026-11-17
before `scripts/ingest_indications.py` existed and before any AUC, and
nothing here may narrow the population, move the threshold, split SARS
from COVID-19, or promote P3 or P4 to a test.

    python scripts/clinical_labels.py

WHAT IT ASKS. Every AUC in this repository is measured against
`HAS_ANTIVIRAL_ACTIVITY_AGAINST` -- in-vitro concentration-response, one
threshold. This asks the same scorer a question that label cannot answer:
does chemical similarity to the compounds measured ACTIVE in vitro
separate the compounds that reached a COVID-19 trial from the compounds
that were screened and did not?

THE STATISTIC IS NOT REIMPLEMENTED. `score_cell` and `corrected` are
imported from `scripts/structure_scorer.py`, so this run uses the exact
exclusion ladder registered by doc2/registrations/prereg_structure_scorer_symmetric.md
-- self, skeleton, then the k = 10 nearest references from every
candidate alike. A second copy of a scoring loop is how two protocols
come to differ without anyone deciding they should, which this
repository has paid for five times.

THE PUBLICATION EXCLUSION IS OFF, by hazard H3: the label does not come
from a publication, so excluding references that share a candidate's
PMID would handicap one class for no reason that applies here.

P2 IS THE LOAD-BEARING TEST, NOT A FORMALITY. A compound in COVID-19
trials is a compound that was noticed, published and edged. If `degree`
matches or beats the structure statistic, P1 is reported as a fame
result. The degree baseline is `baselines.degree_ranking` over a graph
loaded with the held-out predicates skipped -- including `TREATS` and
`IN_TRIAL_FOR` themselves, so the control cannot read its own answer.
That is not the default: `Graph.load` defaults `skip_predicates` to the
empty set, which is how the three scorer studies came to compare against
a null that counted their labels.

WHAT A POSITIVE RESULT WOULD NOT LICENSE, by hazard H2: `drug_indication`
carries no usable per-row date, so there is no temporal split and no way
to separate "the scorer predicts trials" from "trials followed the
assays the scorer's references are drawn from". The verdict says so in
those words.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.baselines import Graph, degree_ranking  # noqa: E402
from kgav.labels import build_labels, evaluate_against_negatives  # noqa: E402
from kgav.pipeline import release_dir  # noqa: E402
from kgav.provenance import write_results  # noqa: E402
from kgav.schema import load_schema  # noqa: E402
from kgav.similarity import fingerprints  # noqa: E402

from structure_scorer import (  # noqa: E402
    K_PRIMARY,
    corrected,
    score_cell,
    small_molecule_props,
)

#: The registered cell. One disease, one virus, and the family is one
#: because SARS and COVID-19 share 354 of 364 indicated compounds.
DISEASE = "MONDO:0100096"
VIRUS = "NCBITaxon:2697049"
CLINICAL_PREDICATES = ("TREATS", "IN_TRIAL_FOR")

#: Family of ONE. Not a Bonferroni correction -- there is one test, and
#: that is registered so a second cell found later cannot be added
#: without restating the alpha.
FAMILY = 1
ALPHA = 0.05 / FAMILY
Z = 1.9600

#: The registration's refusal conditions, as numbers rather than prose.
MIN_POSITIVES = 150
MIN_NEGATIVES = 1000

OUT = ROOT / "data" / "results" / "clinical_labels.json"


def clinical_subjects(release: Path, disease: str) -> dict[str, set[str]]:
    """{predicate: compounds} with a clinical edge to this disease."""
    out: dict[str, set[str]] = {p: set() for p in CLINICAL_PREDICATES}
    for line in (release / "edges.jsonl").open():
        if not any(f'"{p}"' in line for p in CLINICAL_PREDICATES):
            continue
        e = json.loads(line)
        if e["predicate"] in out and e["object"] == disease:
            out[e["predicate"]].add(e["subject"])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    schema = load_schema()
    print(f"release {args.release.name}   cell {DISEASE} / {VIRUS}")
    print(f"family of {FAMILY}, alpha = {ALPHA}, z = {Z}\n")

    clinical = clinical_subjects(args.release, DISEASE)
    in_trial = set().union(*clinical.values()) if clinical else set()
    for pred in CLINICAL_PREDICATES:
        print(f"  {pred:13s} {len(clinical[pred]):5,} compound(s)")
    print(f"  either        {len(in_trial):5,}\n")

    # THE POPULATION IS THE SCREENED ONE. A compound here was measured
    # against SARS-CoV-2, so the absence of a trial is a measured
    # absence rather than an absence of attention.
    labels = build_labels(args.release)
    invitro_pos = set(labels.positives.get(VIRUS, ()))
    invitro_neg = set(labels.negatives.get(VIRUS, ()))
    population = invitro_pos | invitro_neg
    pos = in_trial & population
    neg = population - in_trial

    print(f"  in-vitro labelled population  {len(population):6,}"
          f"   ({len(invitro_pos):,} active / {len(invitro_neg):,} inactive)")
    print(f"  POSITIVES  in a COVID-19 trial or approved  {len(pos):6,}")
    print(f"  NEGATIVES  screened, no curated indication  {len(neg):6,}")
    # THE ORTHOGONALITY THE REGISTRATION EXPECTED, recounted on the run.
    print(f"\n  of the {len(pos):,} positives, {len(pos & invitro_pos):,} are "
          f"in-vitro ACTIVE and {len(pos & invitro_neg):,} INACTIVE")

    refused = []
    if len(pos) < MIN_POSITIVES:
        refused.append(f"{len(pos)} positives, fewer than the registered "
                       f"{MIN_POSITIVES}")
    if len(neg) < MIN_NEGATIVES:
        refused.append(f"{len(neg)} negatives, fewer than the registered "
                       f"{MIN_NEGATIVES}")
    stray = in_trial - population
    if refused:
        print("\nREFUSED, on the registration's own conditions:")
        for r in refused:
            print(f"  - {r}")
        write_results(args.out, {
            "measured_by": "scripts/clinical_labels.py",
            "registered_by": "doc2/registrations/prereg_clinical_labels.md",
            "refused": refused,
        }, args.release, schema.version)
        return 1

    # REFERENCES ARE THE IN-VITRO ACTIVES, which is what makes this a
    # question about the in-vitro axis predicting the clinical one. A
    # candidate that is itself a reference loses itself (exclusion 1).
    props = small_molecule_props(args.release)
    candidates = sorted(population)
    references = sorted(invitro_pos)
    # EVERY COMPOUND THE RUN WILL TOUCH, fingerprinted once: the
    # population, the references, and P3's clinical-label-only set.
    touched = sorted(set(candidates) | set(references) | set(stray))
    fps, fp_stats = fingerprints(
        {i: ((props.get(i) or {}).get("smiles") or "") for i in touched})
    print(f"\n  fingerprinted {fp_stats['fingerprinted']:,} of "
          f"{len(touched):,}  (unparseable {fp_stats['unparseable']:,}, "
          f"no SMILES {fp_stats['no_smiles']:,})")

    pubs: dict[str, frozenset] = {}
    scores, counts = score_cell(
        candidates, references, fps, props, pubs,
        top_k=1, publication_disjoint=False,
        positives=frozenset(pos), knn=K_PRIMARY)
    structure = corrected(evaluate_against_negatives(scores, pos, neg), z=Z)
    structure["counts"] = counts

    g = Graph.load(args.release,
                   skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    degrees = degree_ranking(g)
    degree = corrected(evaluate_against_negatives(
        {d: v for d, v in degrees.items() if d in population}, pos, neg), z=Z)

    # P3: AN OBSERVATION, NOT A TEST. These compounds carry a clinical
    # label and no in-vitro one, so they cannot enter the population
    # above and there is no screened negative class for them.
    unlabelled = sorted(stray & set(props))
    p3_scores, _ = score_cell(
        unlabelled, references, fps, props, pubs, top_k=1,
        publication_disjoint=False, positives=frozenset(unlabelled),
        knn=K_PRIMARY) if unlabelled else ({}, {})

    def pct(xs, q):
        xs = sorted(xs)
        return xs[min(len(xs) - 1, int(q * len(xs)))] if xs else None

    pos_scores = [scores[c] for c in pos if c in scores]
    neg_scores = [scores[c] for c in neg if c in scores]

    print("\n" + "=" * 70)
    print(f"{'channel':12s} {'n+':>5s} {'n-':>7s} {'AUC':>8s}"
          f" {'95% interval':>18s}  clears alpha")
    for name, m in (("structure", structure), ("degree", degree)):
        lo, hi = m.get("auc_ci_lo_corrected"), m.get("auc_ci_hi_corrected")
        print(f"{name:12s} {m['n_pos']:5,} {m['n_neg']:7,} {m['auc']:8.4f}"
              f"   [{lo:.4f}, {hi:.4f}]  "
              f"{'YES' if m.get('clears_corrected_alpha') else 'no'}")
    print("=" * 70)

    verdict = "degree matches or beats the structure statistic" \
        if degree["auc"] >= structure["auc"] else \
        "the structure statistic is above its degree baseline"
    print(f"\nP2: {verdict}.")
    if degree["auc"] >= structure["auc"]:
        print("    So P1 is a FAME result, not a chemistry one -- the "
              "registration\n    fixed that reading in advance.")

    payload = {
        "measured_by": "scripts/clinical_labels.py",
        "registered_by": "doc2/registrations/prereg_clinical_labels.md",
        "disease": DISEASE,
        "virus": VIRUS,
        "family_size": FAMILY,
        "alpha": ALPHA,
        "z": Z,
        "statistic": ("max Tanimoto to in-vitro actives, symmetric k=10 "
                      "exclusion, no publication exclusion (hazard H3)"),
        "clinical_edges": {p: len(v) for p, v in clinical.items()},
        "population": {
            "labelled": len(population),
            "invitro_active": len(invitro_pos),
            "invitro_inactive": len(invitro_neg),
            "positives": len(pos),
            "negatives": len(neg),
            "positives_invitro_active": len(pos & invitro_pos),
            "positives_invitro_inactive": len(pos & invitro_neg),
        },
        "fingerprinted": dict(fp_stats),
        "structure_max_k10": structure,
        "degree": degree,
        "p2_degree_at_or_above_structure": bool(
            degree["auc"] >= structure["auc"]),
        "p3_observation": {
            "clinical_label_only_in_cohort": len(stray),
            "scorable_here": len(unlabelled),
            "scored": len(p3_scores),
            "median_score": pct(list(p3_scores.values()), 0.5),
            "median_score_positives": pct(pos_scores, 0.5),
            "median_score_negatives": pct(neg_scores, 0.5),
            "note": ("a separate population with no screened negative "
                     "class; reported, not tested"),
        },
        "p4_untestable": {
            "note": ("the viral-target stratum is 13 of 222 and is "
                     "registered as not a cell"),
        },
        "no_temporal_split": (
            "drug_indication carries no usable per-row date, so a positive "
            "result does not license a predictive claim (hazard H2)"),
    }
    write_results(args.out, payload, args.release, schema.version)
    print(f"\nwrote {args.out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
