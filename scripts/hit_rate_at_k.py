#!/usr/bin/env python3
"""The wet-lab hit rate: of the top-ranked compounds, how many were active.

    make hit-rate

THE QUESTION. "What is the probability that a suggested compound works in
the lab?" AUC does not answer it. AUC is the probability that a random
active outranks a random inactive -- a statement about ORDERING, with no
units a bench scientist can act on. The quantity asked for is the
positive predictive value of the top of the list: order the top k, screen
them, what fraction come back active.

WHY THE SHORTLIST CANNOT ALREADY ANSWER IT.
scripts/repurposing_candidates.py prints, above every list it produces:

    "So every row below is a chemical-analogy suggestion with a measured
     hit rate."

There was no measured hit rate. Nothing in data/results recorded a
precision@k, a PPV or a hits@k for the structure scorer -- checked across
structure_scorer.json, structure_scorer_strength_matched.json,
temporal_structure_scorer.json and repurposing_candidates.json. The
banner quoted AUC and called it a hit rate, which is the one translation
AUC does not support.

THIS IS A POST-HOC STATISTIC AND IS LABELLED AS ONE. The registration
docs/prereg_temporal_structure_scorer.md fixed AUC as the statistic
before any AUC existed. Precision@k is computed here AFTER those AUCs
were read, on the same held-out labels, so it is a DESCRIPTIVE re-reading
of an already-measured result and not a new test. The hypergeometric
p-value is printed because a hit rate with no null is unreadable, and it
carries no correction for having been chosen after the fact. Nothing here
may be claimed as significant, and nothing here is written into
data/results as a claim-carrying result.

MEASURED IN THE PROSPECTIVE SETTING, which is the only one that answers
the question as asked. The references are activity published at or before
the cutoff; the candidates are compounds whose activity was measured
after it. That is the situation of someone standing in 2021 deciding what
to screen. The cross-sectional cells cannot answer it at any base rate,
because their references and candidates come from the same period.

TIES ARE BOUNDED, NOT RESOLVED. `rank()` breaks ties by dict insertion
order, so a tie block straddling position k makes hits@k a draw rather
than a measurement. baselines.hits_at_k_range returns [low, high] and the
hypergeometric centre; all three are printed, because a point estimate at
the centre of a wide interval reads as a measurement.

THE BASE RATE IS THE WHOLE ANSWER AND IT IS NOT A PROPERTY OF THE SCORER.
A cell where half the screened compounds are active makes any ranker look
excellent in precision and tells you nothing. Enrichment -- precision@k
over the cell's own base rate -- is the part attributable to ranking, and
it is printed beside every precision.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.baselines import hits_at_k_range, lift_pvalue, rank
from kgav.pipeline import release_dir, train_dir
from kgav.provenance import release_fingerprint, write_results
from kgav.similarity import fingerprints

from structure_scorer import POPULATIONS, score_cell, small_molecule_props
from temporal_structure_scorer import CUTOFF, MIN_CLASS, UNDATED, cells, evaluable

#: The cut depths reported. 10 is the length of a tier in
#: repurposing_candidates.py's output, so precision@10 is the figure that
#: describes what somebody is actually handed; the rest bracket it.
KS = (5, 10, 20, 50)

#: The registered secondary control: every candidate loses its 10 nearest
#: references, positives and negatives alike. structure_scorer.K_PRIMARY is
#: the same number; named here so this file reads without it.
K10 = 10

#: Score given to a candidate score_cell could not score. Tanimoto is in
#: [0, 1], so this ranks them below everything real -- which is what the
#: AUC evaluator does when it imputes them to the floor.
UNSCORED = -1.0

RULE = "=" * 100


def precision_rows(cand: list[str], pos: list[str], scores: dict,
                   ks: tuple[int, ...]) -> dict:
    """precision@k over one cell's scored candidates, ties bounded."""
    full = {d: scores.get(d, UNSCORED) for d in cand}
    ranked = rank(full)
    pos_set = set(pos)
    pool, n_pos = len(cand), len(pos)
    base = n_pos / pool if pool else 0.0
    per_k = {}
    for k in ks:
        if k > pool:
            continue
        lo, hi, tied, tied_pos, exp = hits_at_k_range(ranked, pos_set, k)
        per_k[k] = {
            "hits_lo": lo, "hits_hi": hi, "hits_expected": exp,
            "tied_at_cut": tied, "tied_at_cut_positives": tied_pos,
            "determined": lo == hi,
            "precision_lo": lo / k, "precision_hi": hi / k,
            "precision_expected": exp / k,
            "enrichment_expected": (exp / k) / base if base else None,
            "p_value": lift_pvalue(pool, n_pos, k, lo),
        }
    return {"pool": pool, "n_pos": n_pos, "base_rate": base, "per_k": per_k}


def plate_block(release: Path, props: dict, ks: tuple[int, ...]) -> list[dict]:
    """The assay-matched plate cells, where the base rate is realistic.

    WHY THIS BLOCK EXISTS AND WHY IT IS SECOND. The prospective cells
    above answer "would this have worked", and their base rates are 6% to
    70% because a test set is made of compounds somebody chose to screen.
    Nobody screening a library faces a 50% prior. The NCATS plate cells
    do have a library-like base rate -- one protocol, one concentration
    range, both classes off the same plates -- and that is the population
    whose precision translates into "order ten compounds, how many hit".

    IT IS CROSS-SECTIONAL, which is the cost. References and candidates
    come from the same period, so a candidate's nearest known active may
    be its own analogue from the same paper. The strength-matched control
    is this project's answer to that and is what is used here: k_cell
    nearest references removed from every candidate alike, with k_cell
    read from the registered run rather than recomputed.
    """
    from kgav.labels import build_labels

    doc = json.loads((ROOT / "data" / "results"
                      / "structure_scorer_strength_matched.json").read_text())
    k_map = doc["k_cells"]

    vy = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    names = {f"NCBITaxon:{e['taxon']}": e["name"]
             for lst in ("viruses", "surveyed") for e in (vy.get(lst) or [])}

    print(RULE)
    print("THE SAME QUESTION ON A LIBRARY-LIKE BASE RATE (assay-matched "
          "plate cells)")
    print(RULE)
    print()
    print("  Cross-sectional. These are the cells whose base rate "
          "resembles a screening")
    print("  library rather than a curated test set, which is what makes "
          "their precision")
    print("  translate into 'order ten compounds, how many hit'.")
    print()
    print("  Reported at THREE exclusion strengths. 'no exclusion' lets a "
          "candidate match")
    print("  its own analogues among the known actives; 'k10' withholds "
          "each candidate's")
    print("  ten nearest; 'k_cell' withholds this cell's own measured "
          "asymmetry and is")
    print("  the registered control. A precision that collapses down the "
          "ladder was")
    print("  retrieving analogues, not predicting activity.")
    print()
    out = []
    for key, source, assay_types in POPULATIONS:
        if assay_types is None:
            continue                       # ChEMBL: not plate-matched
        labels = build_labels(release, sources={source},
                              assay_types=assay_types)
        for virus in sorted(set(labels.positives) | set(labels.negatives)):
            pos = sorted(labels.positives.get(virus, set()))
            neg = sorted(labels.negatives.get(virus, set()))
            if len(pos) < MIN_CLASS or len(neg) < MIN_CLASS:
                continue
            k_cell = k_map.get(virus)
            if k_cell is None:
                continue
            cand = sorted(set(pos) | set(neg))
            refs = sorted(pos)
            everyone = sorted(set(cand) | set(refs))
            fps, _ = fingerprints(
                {i: (props.get(i, {}).get("smiles") or "") for i in everyone})
            # THE LADDER, NOT ONE RUNG. This is the measurement that
            # decides the question, so it is reported at all three
            # exclusion strengths rather than only at the registered one.
            # Measured on the SARS-CoV-2 cytopathic cell (base 1.20%):
            # precision@10 is 50.0% with no exclusion, 0.0% at k10, 0.0%
            # at this cell's k_cell. An apparent hit rate that vanishes
            # when a candidate's nearest known actives are withheld was
            # never a hit rate -- and the AUC does not show it, because
            # AUC averages over the whole ranking while a screen only
            # ever orders the top of it.
            ladder = {}
            for cname, knn in (("no exclusion", 0), ("k10", K10),
                               (f"k_cell={k_cell}", k_cell)):
                sc, _ = score_cell(cand, refs, fps, props, {},
                                   top_k=1, publication_disjoint=False,
                                   positives=frozenset(pos), knn=knn)
                ladder[cname] = precision_rows(cand, pos, sc, ks)
            r = dict(ladder[f"k_cell={k_cell}"])
            r.update({"taxon": virus, "name": names.get(virus, virus),
                      "population": key, "k_cell": k_cell,
                      "controls": ladder})
            out.append(r)

            print(f"  {r['name']}  —  {'/'.join(sorted(assay_types))}  —  "
                  f"pool {r['pool']:,}, {r['n_pos']} active, "
                  f"base rate {r['base_rate']:.2%}")
            sub = (f"    {'control':<16}" +
                   "".join(f"{'p@' + str(k):>9}" for k in ks) +
                   f"{'enrich@' + str(ks[-1]):>12}")
            print(sub)
            print("    " + "-" * (len(sub) - 4))
            for cname, m in ladder.items():
                line = f"    {cname:<16}"
                for k in ks:
                    c = m["per_k"].get(k)
                    line += (f"{c['precision_expected']:>8.1%}" if c
                             else f"{'--':>9}")
                last = m["per_k"].get(ks[-1])
                line += (f"{last['enrichment_expected']:>11.2f}x" if last
                         else f"{'--':>12}")
                print(line)
            print()
    if not out:
        print("  no plate cell is readable on this release.")
        print()
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--split", type=Path, default=train_dir(ROOT, CUTOFF))
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--k", type=int, action="append", default=None,
                    help="cut depth; repeatable. Default 5/10/20/50.")
    ap.add_argument("--out", type=Path, default=None,
                    help="write the table as JSON. NOT data/results by "
                         "default: this is a post-hoc statistic and is not "
                         "a claim-carrying result.")
    args = ap.parse_args(argv)
    ks = tuple(sorted(set(args.k))) if args.k else KS

    derived = json.loads((args.split / "DERIVED_FROM.json").read_text())
    split_sha = ((derived.get("provenance") or {})
                 .get("release") or {}).get("edges_sha")
    release_sha = release_fingerprint(args.release).get("edges_sha")

    print(RULE)
    print("PROSPECTIVE HIT RATE — of the top-ranked compounds, how many came "
          "back active")
    print(RULE)
    print()
    print(f"  split    {args.split.name}   cutoff {derived.get('cutoff')}   "
          f"undated {derived.get('undated')}")
    print(f"           built from edges_sha {split_sha}")
    print(f"  release  {args.release.name}   edges_sha {release_sha}")
    if split_sha != release_sha:
        print()
        print("  REFUSED: the split was built from a different release than "
              "the one on disk.")
        return 1
    if derived.get("cutoff") != CUTOFF or derived.get("undated") != UNDATED:
        print()
        print(f"  REFUSED: this reads the {CUTOFF}/undated={UNDATED} split "
              f"and the one on disk is "
              f"{derived.get('cutoff')}/{derived.get('undated')}.")
        return 1

    labels = json.loads((args.split / "TEST_LABELS.json").read_text())
    ev = evaluable(cells(labels))
    props = small_molecule_props(args.release)
    everyone = sorted({d for v in ev.values() for part in v for d in part})
    fps, _ = fingerprints(
        {i: (props.get(i, {}).get("smiles") or "") for i in everyone})

    vy = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    names = {f"NCBITaxon:{e['taxon']}": e["name"]
             for lst in ("viruses", "surveyed") for e in (vy.get(lst) or [])}

    print()
    print("  'base' is the fraction of this cell's screened compounds that "
          "are active. It is")
    print("  a property of what was screened, NOT of the scorer, and it is "
          "the number any")
    print("  precision below has to beat. 'enrich' is precision@k / base: "
          "the part that is")
    print("  attributable to ranking. 'p' is the hypergeometric tail, "
          "uncorrected and")
    print("  post-hoc — descriptive only.")
    print()

    rows = []
    for taxon, (refs, pos, neg) in sorted(ev.items()):
        name = names.get(taxon, taxon)
        cand = list(pos) + list(neg)
        # IDENTICAL SETTINGS TO THE REGISTERED PRIMARY in
        # temporal_structure_scorer.py: no exclusion beyond the split,
        # because the reference set is already strictly earlier than the
        # candidate and removing a near neighbour would remove the signal
        # rather than a leak.
        # BOTH REGISTERED CONTROLS, because the primary alone cannot tell
        # a prospective hit from an analogue lookup. The registration's
        # primary applies no exclusion -- the reference set is already
        # strictly earlier, and removing a near neighbour would remove
        # signal rather than a leak. That argument is about TIME and says
        # nothing about CHEMISTRY: a 2022 compound can be an analogue of a
        # 2019 one, and then "it resembles a known active" is the whole of
        # what the top of the list means. k10 is the registered secondary
        # and is what separates the two. Measured: Zika's precision@10
        # falls 70% -> 10% under it, and DENV-2's holds at 70%.
        per_control = {}
        for cname, knn in (("primary", 0), ("symmetric_k10", K10)):
            sc, _ = score_cell(cand, refs, fps, props, {},
                               publication_disjoint=False, knn=knn)
            per_control[cname] = precision_rows(cand, pos, sc, ks)
        r = dict(per_control["primary"])
        r.update({"taxon": taxon, "name": name, "n_refs": len(refs),
                  "controls": per_control})
        rows.append(r)

    hdr = (f"  {'virus':<13}{'pool':>7}{'pos':>7}{'base':>7}{'k':>4}"
           f"{'hits':>9}{'prec@k':>9}{'enrich':>8}"
           f"{'  |':>3}{'prec k10':>10}{'enrich':>8}{'p':>9}  ties")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        first = True
        k10 = r["controls"]["symmetric_k10"]["per_k"]
        for k, m in r["per_k"].items():
            label = f"{r['name']:<13}" if first else " " * 13
            pool = f"{r['pool']:>7,}" if first else " " * 7
            npos = f"{r['n_pos']:>7,}" if first else " " * 7
            baset = f"{r['base_rate']:>6.1%}" if first else " " * 6
            first = False
            hits = (f"{m['hits_lo']}" if m["determined"]
                    else f"{m['hits_lo']}-{m['hits_hi']}")
            p = m["p_value"]
            s2 = k10.get(k) or {}
            # PRECOMPUTED, NOT NESTED. An f-string reusing the outer quote
            # is Python 3.12 syntax and pyproject pins target-version
            # py311, so ruff refused it -- correctly: it ran here only
            # because this interpreter happens to be 3.12.
            p2 = f"{s2['precision_expected']:.1%}" if s2 else "--"
            e2 = f"{s2['enrichment_expected']:.2f}" if s2 else "--"
            pv = f"{p:.2g}" if p is not None else "--"
            print(f"  {label}{pool}{npos}{baset} {k:>4}{hits:>9}"
                  f"{m['precision_expected']:>8.1%}"
                  f"{m['enrichment_expected']:>8.2f}"
                  f"{'  |':>3}{p2:>10}{e2:>8}{pv:>9}"
                  f"  {m['tied_at_cut']}")
        print()
    print("  Left of the bar: the registered primary (no exclusion). Right: "
          "the registered")
    print("  secondary, which removes each candidate's 10 nearest known "
          "actives. A column")
    print("  that collapses across the bar was ranking analogues, not "
          "predicting activity.")
    print()

    plate_rows = plate_block(args.release, props, ks)

    # ------------------------------------------------------------ summary
    print(RULE)
    print("WHAT THIS MEANS FOR A COMPOUND ON THE SHORTLIST")
    print(RULE)
    print()
    k_hdr = 10 if 10 in ks else ks[0]
    usable = [r for r in rows if k_hdr in r["per_k"]]
    if usable:
        precs = [r["per_k"][k_hdr]["precision_expected"] for r in usable]
        bases = [r["base_rate"] for r in usable]
        enr = [r["per_k"][k_hdr]["enrichment_expected"] for r in usable]
        precs_s, bases_s, enr_s = sorted(precs), sorted(bases), sorted(enr)

        def med(xs):
            n = len(xs)
            return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2

        print(f"  Across the {len(usable)} prospective cells, at k={k_hdr}:")
        print(f"    median precision@{k_hdr}        {med(precs_s):.1%}"
              f"   (range {min(precs):.1%} to {max(precs):.1%})")
        print(f"    median base rate          {med(bases_s):.1%}"
              f"   (range {min(bases):.1%} to {max(bases):.1%})")
        print(f"    median enrichment         {med(enr_s):.2f}x"
              f"  (range {min(enr):.2f} to {max(enr):.2f})")
        print()
        print("  THE BASE RATE IS DOING MOST OF THE WORK. These test sets are")
        print("  compounds somebody chose to screen, and people screen things")
        print("  they expect to work, so the base rate here is far above what")
        print("  an arbitrary compound would have. Read the enrichment column")
        print("  for what the ranking adds; read the base column for what the")
        print("  population was already worth.")
        print()
    print("  AND THE SHORTLIST IS NOT THIS POPULATION. Every figure above is")
    print("  measured on compounds that WERE screened, which is how their")
    print("  activity is known. The shortlist in data/results/"
          "repurposing_candidates.json")
    print("  is drawn from compounds that were NOT -- that is what makes them")
    print("  suggestions. Carrying these rates across assumes the unscreened")
    print("  pool has the same base rate as the screened one, and screening is")
    print("  not a random sample: it is enriched for whatever looked promising")
    print("  at the time. The number above is therefore an UPPER BOUND on what")
    print("  a shortlist row should be expected to return.")
    print()
    print("  AND 'ACTIVE' HERE IS AN ASSAY READING, not a drug. The label is")
    print("  HAS_ANTIVIRAL_ACTIVITY_AGAINST at the project's potency")
    print("  threshold, measured in whatever cell line the depositor used --")
    print("  for SARS-CoV-2 very often Vero E6, which lacks meaningful")
    print("  TMPRSS2 and scores endosomal-entry compounds that then fail in")
    print("  airway cells. Selectivity, cytotoxicity and exposure are not in")
    print("  it.")
    print()

    if args.out:
        payload = {
            "measured_by": "scripts/hit_rate_at_k.py",
            "post_hoc": True,
            "not_registered": (
                "precision@k was computed after the registered AUCs were "
                "read, on the same held-out labels. Descriptive only; the "
                "p-values carry no correction for selection."),
            "split": args.split.name,
            "cutoff": derived.get("cutoff"),
            "undated_policy": derived.get("undated"),
            "ks": list(ks),
            "cells": rows,
            "plate_cells": plate_rows,
        }
        # THROUGH write_results, NOT json.dumps. Caught by
        # tests/test_provenance.py::test_every_script_writing_a_result_
        # uses_write_results, and the gate was right: the hand-rolled
        # version recorded `"split": <directory name>`, and a split
        # rebuilt in place keeps its name while being a different graph.
        # That is the exact failure pathway_landscape.json had. The name
        # stays because it is useful; the fingerprint is now beside it.
        write_results(args.out, payload, args.release)
        print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
