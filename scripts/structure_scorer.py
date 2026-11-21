#!/usr/bin/env python3
"""The structure-derived scorer, exactly as registered.

PROTOCOL: doc2/registrations/prereg_structure_scorer.md, signed 2026-11-11 before this
file existed. Nothing here is a choice -- every parameter was fixed in that
document, and `tests/test_prereg_structure_scorer.py` re-derives each one
from the release so this script cannot quietly use a different value.

WHAT IT MEASURES. For each (virus, compound): the maximum Tanimoto
similarity from the compound's ECFP4 fingerprint to the fingerprints of
that virus's known actives, over a reference set from which three things
have been removed --

    1. the compound itself
    2. every active sharing its `inchikey_skel`
    3. every active citing a publication the compound also cites

-- and the ranking of those scores against the measured negatives.

WHY THE THIRD EXCLUSION IS THE WHOLE EXERCISE. A max-similarity scorer
ranks a compound by its nearest known active, and in this label set the
nearest known active is very often the analog published in the same table:
`baselines.label_publication_index`'s docstring records that 1,015 of 1,070
compounds -- 94.9% -- share a PMID between their label and the edge scored
against it. So the run computes the scorer BOTH WAYS, and the difference is
reported as a leak estimate rather than as a result. Registered prediction
P1 is that the publication-disjoint variant reads lower in at least 8 of the
14 cells.

WHY THE POPULATION IS ChEMBL-ONLY. `labels.build_labels`' own docstring:
pooling ChEMBL with the NCATS panel mixes antiviral research compounds
against library compounds, "and a scorer can separate the two classes by
recognising which collection a compound came from". This is the scorer that
would do exactly that, so the pooled task is not computed here at all. The
NCATS panel alone -- one screen, one protocol, both classes off the same
plates -- is the matched secondary cell and carries P5.

A CANDIDATE WHOSE REFERENCE SET EMPTIES IS NOT SCORED, and that is the
registered control biting rather than a bug. A positive whose entire
congeneric series sits in one paper loses every admissible reference, has no
evidence, and takes the imputed floor -- which is the doctrine
`labels.evaluate_against_negatives` already states for a compound with no
path. It is counted as `empty_reference` and it moves `coverage_pos`, so it
shows up in `auc_floor` instead of disappearing.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import (                              # noqa: E402
    Graph,
    degree_ranking,
    label_publication_index,
)
from kgav.labels import (                                 # noqa: E402
    ACTIVITY_PREDICATE,
    auc_interval,
    build_labels,
    evaluate_against_negatives,
)
from kgav.pipeline import release_dir                     # noqa: E402
from kgav.provenance import write_results                 # noqa: E402
from kgav.schema import load_schema                       # noqa: E402
from kgav.similarity import fingerprints                  # noqa: E402

# REGISTERED, NOT CHOSEN HERE. Each of these is pinned to the document by
# tests/test_prereg_structure_scorer.py.
PRIMARY_SOURCE = "infores:chembl"
SECONDARY_SOURCE = "infores:ncats-opendata"

#: (output key, knowledge source, assay_type restriction). The plate-matched
#: populations are split BY PROTOCOL, not just by source, which is what
#: makes each one single-screen -- see run_population. The CPE population
#: yields one cell (SARS-CoV-2) and the entry population three (SARS-CoV-2,
#: SARS-CoV, MERS-CoV), which are the four cells
#: doc2/registrations/prereg_ncats_entry_assays.md fixes at alpha = 0.05/4.
#:
#: `secondary_ncats` KEEPS ITS NAME for the CPE population, because that is
#: the key prereg_structure_scorer.md's registered cell is recorded under
#: and renaming it would silently orphan every comparison against it.
POPULATIONS = (
    ("primary_chembl", PRIMARY_SOURCE, None),
    ("secondary_ncats", SECONDARY_SOURCE, {"cell_based_antiviral"}),
    ("secondary_ncats_entry", SECONDARY_SOURCE,
     {"cell_based_pseudotype_entry"}),
)
MIN_CLASS = 10              # calibration.MIN_REACHED_POSITIVES
FAMILY = 14                 # the readable ChEMBL cells
ALPHA = 0.05 / FAMILY
Z_CORRECTED = 2.9137        # two-sided, for ALPHA

# doc2/registrations/prereg_ncats_entry_assays.md, second amendment 2026-11-17. The
# registration was signed for FOUR plate cells; one of them cannot exist,
# because edge identity excludes `assay_type` and the SARS-CoV-2 entry
# measurements merged into the CPE edges. THREE is the family that was
# achievable, and it is written here so a run that produces a different
# number says so instead of dividing by whatever it found.
PLATE_FAMILY = 3
PLATE_ALPHA = 0.05 / PLATE_FAMILY
PLATE_Z_CORRECTED = 2.3940  # two-sided, for PLATE_ALPHA
TOP_K = 3                   # the one registered alternative statistic

# doc2/registrations/prereg_structure_scorer_symmetric.md, signed 2026-11-11.
# K_PRIMARY is derived, not chosen: JEV holds 18 references after
# self-exclusion, the smallest in the family, and k must leave every cell
# non-empty or the study repeats the prior run's P5 -- a cell that printed
# 0.5000 with an interval because nothing had been ranked. k=30, which the
# prior verdict suggested, leaves -12.
K_CURVE = (1, 3, 10)
K_PRIMARY = 10
CONTROLS = ("publication", "knn", "strength")


def k_cells_from_prior(prior: dict) -> dict[str, int]:
    """virus -> k_cell, from the publication run's own counters.

    doc2/registrations/prereg_structure_scorer_strength_matched.md:

        k_cell = ceil(max(pub_per_positive, pub_per_negative))

    READ, NOT RECOMPUTED, and the registration says why: the k values are
    what make this study auditable, so they come off the artifact a reader
    can inspect rather than out of a fresh traversal that might differ. A
    test re-derives all fourteen from the same file.

    Returns one entry per READABLE cell of the publication run. A cell
    missing from this map gets no strength-matched score, which is correct
    -- there is nothing to match it to.
    """
    out: dict[str, int] = {}
    for virus, cell in prior.get("primary_chembl", {}).items():
        if not cell.get("readable"):
            continue
        n = cell["structure_max"]["counts"]
        ep = n["excluded_publication_pos"] / max(n["candidates_pos"], 1)
        en = n["excluded_publication_neg"] / max(n["candidates_neg"], 1)
        out[virus] = math.ceil(max(ep, en))
    return out


def small_molecule_props(release: Path) -> dict[str, dict]:
    import json
    out = {}
    for line in (release / "nodes.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        n = json.loads(line)
        if n["class"] == "SmallMolecule":
            out[n["id"]] = n.get("properties") or {}
    return out


def score_cell(candidates: list[str], references: list[str],
               fps: dict, props: dict[str, dict],
               pubs: dict[str, set[str]], top_k: int = 1,
               publication_disjoint: bool = True,
               positives: frozenset[str] = frozenset(),
               knn: int = 0) -> tuple[dict, dict]:
    """Max (or mean-of-top-k) Tanimoto to the admissible references.

    Returns (scores by candidate id, counters). A candidate absent from
    `scores` was not scorable and will be imputed to the floor by the
    evaluator -- either it has no fingerprint or its reference set emptied.
    """
    from rdkit import DataStructs

    ref_ids = [r for r in references if r in fps]
    ref_fps = [fps[r] for r in ref_ids]
    ref_skel = [props.get(r, {}).get("inchikey_skel", "") for r in ref_ids]
    ref_pubs = [pubs.get(r, frozenset()) for r in ref_ids]
    index = {r: i for i, r in enumerate(ref_ids)}

    scores: dict[str, float] = {}
    counts = {
        "candidates": len(candidates),
        "references": len(ref_ids),
        "references_without_fingerprint": len(references) - len(ref_ids),
        "no_fingerprint": 0,
        "empty_reference": 0,
        "excluded_self": 0,
        "excluded_skeleton": 0,
        "excluded_publication": 0,
        "without_publication": 0,
        # PER CLASS, BECAUSE THE CONTROL IS NOT SYMMETRIC AND THE FIRST RUN
        # COULD NOT SHOW IT. A positive's nearest known actives are its own
        # paper-mates, which exclusion (3) removes; a negative sitting in a
        # paper that contains no active loses nothing. So the exclusion can
        # push an AUC down by handicapping one class, and four cells read
        # below chance. These two counters are what tells that apart from
        # the scorer genuinely anti-ranking.
        "excluded_publication_pos": 0,
        "excluded_publication_neg": 0,
        "candidates_pos": 0,
        "candidates_neg": 0,
        # EXCLUSION (3'), doc2/registrations/prereg_structure_scorer_symmetric.md. P2 is
        # an invariant over these: every candidate must lose exactly `knn`
        # references, so excluded_knn == knn * candidates in every cell. If
        # it does not, the control is not symmetric and nothing else in the
        # cell means anything.
        "excluded_knn": 0,
        "excluded_knn_pos": 0,
        "excluded_knn_neg": 0,
        "knn_short": 0,
        # THE RESIDUAL. Exclusion (2) still removes a variable number --
        # the registration bounds it at 2 extra for a minority of
        # candidates -- so it is counted per class rather than argued
        # about.
        "excluded_skeleton_pos": 0,
        "excluded_skeleton_neg": 0,
    }
    for c in candidates:
        fp = fps.get(c)
        if fp is None:
            counts["no_fingerprint"] += 1
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
        c_skel = props.get(c, {}).get("inchikey_skel", "")
        c_pubs = pubs.get(c, frozenset())
        is_pos = c in positives
        counts["candidates_pos" if is_pos else "candidates_neg"] += 1
        if not c_pubs:
            counts["without_publication"] += 1
        drop = set()
        if c in index:
            drop.add(index[c])
            counts["excluded_self"] += 1
        for i, s in enumerate(ref_skel):
            if i not in drop and c_skel and s == c_skel:
                drop.add(i)
                counts["excluded_skeleton"] += 1
                counts["excluded_skeleton_pos" if is_pos
                       else "excluded_skeleton_neg"] += 1
        if publication_disjoint and c_pubs:
            for i, pp in enumerate(ref_pubs):
                if i not in drop and (pp & c_pubs):
                    drop.add(i)
                    counts["excluded_publication"] += 1
                    counts["excluded_publication_pos" if is_pos
                           else "excluded_publication_neg"] += 1
        # EXCLUSION (3'): the k NEAREST remaining references, applied after
        # (1) self and (2) skeleton so it removes exactly k beyond them.
        #
        # TIES BROKEN BY REFERENCE INDEX, ASCENDING, and `ref_ids` is
        # sorted, so the choice among equal similarities is deterministic
        # across processes. It has to be: str hashing is randomised per
        # interpreter, and set iteration order reaching a results file is a
        # defect this repository has already paid for once in build_labels.
        if knn:
            order = sorted(((-s, i) for i, s in enumerate(sims)
                            if i not in drop))
            if len(order) <= knn:
                # Would empty the set. The registration's k was derived so
                # this cannot happen in the family -- JEV holds 18
                # references and k is 10 -- so it is counted loudly rather
                # than silently clipped.
                counts["knn_short"] += 1
            for _neg_s, i in order[:knn]:
                drop.add(i)
                counts["excluded_knn"] += 1
                counts["excluded_knn_pos" if is_pos
                       else "excluded_knn_neg"] += 1
        kept = [s for i, s in enumerate(sims) if i not in drop]
        if not kept:
            counts["empty_reference"] += 1
            continue
        if top_k <= 1:
            scores[c] = max(kept)
        else:
            kept.sort(reverse=True)
            scores[c] = sum(kept[:top_k]) / len(kept[:top_k])
    return scores, counts


def corrected(metrics: dict, z: float = Z_CORRECTED,
              also: float | None = None) -> dict:
    """The same AUC at the registered corrected level, and a refusal.

    THE LEVEL IS AN ARGUMENT BECAUSE TWO FAMILIES ARE REGISTERED, and
    until 2026-11-17 this function applied the ChEMBL family's z to every
    cell it was handed. The plate-matched cells are registered as their
    own family of three (doc2/registrations/prereg_ncats_entry_assays.md, second
    amendment): alpha 0.05/3, z 2.3940, not 0.05/14 and 2.9137. Declaring
    PLATE_Z_CORRECTED and then never passing it is the same shape as a
    counter that cannot be printed -- the constant documents an intention
    the code does not carry out. Measured cost of the bug: the MERS-CoV
    entry cell read 0.5873 with corrected bounds [0.4822, 0.6923] and was
    marked as not clearing; at its own registered level the bounds are
    [0.5010, 0.6736] and it clears.

    `also` records a SECOND level without choosing between them. The
    SARS-CoV-2 CPE cell is claimed by two registrations --
    prereg_structure_scorer.md reports it "at the same corrected level" as
    the 14-cell family, and prereg_ncats_entry_assays.md makes it one of
    three plate cells -- so both are computed and both are printed, and
    neither document gets to be the one silently applied.

    A CELL NOTHING WAS SCORED IN IS NOT A NULL AND MUST NOT PRINT AS ONE.
    The first run of this script reported the matched NCATS cell as `0.5000
    [0.4103, 0.5897]` -- an AUC with an interval, reading exactly like a
    scorer that found nothing. It was not: all 7,432 candidates had lost
    their entire reference set to the publication exclusion, so nothing was
    ranked, every compound tied at the imputed floor, and a rank-sum AUC
    over one tie block is 0.5 by arithmetic. `coverage_pos` was 0.0 and
    `auc_skill` was None, and the printed line showed neither.

    That is this repository's recurring defect in a new costume -- a
    quantity that cannot fail reading like one that was tested -- so it is
    refused here rather than annotated downstream.
    """
    if not metrics.get("evaluable"):
        return metrics
    out = dict(metrics)
    if not metrics.get("coverage_pos") or not metrics.get("coverage_neg"):
        out["measurable"] = False
        out["refusal"] = (
            "nothing was scored in this cell: every candidate's reference "
            "set emptied, so the AUC is 0.5 by arithmetic and not by "
            "measurement")
        out.pop("auc", None)
        return out
    out["measurable"] = True
    out["z_corrected"] = z
    lo, hi = auc_interval(metrics["auc"], metrics["n_pos"], metrics["n_neg"],
                          z=z)
    out["auc_ci_lo_corrected"] = lo
    out["auc_ci_hi_corrected"] = hi
    if also is not None:
        lo2, hi2 = auc_interval(metrics["auc"], metrics["n_pos"],
                                metrics["n_neg"], z=also)
        out["z_alternate"] = also
        out["auc_ci_lo_alternate"] = lo2
        out["auc_ci_hi_alternate"] = hi2
        out["clears_alternate_alpha"] = bool(lo2 > 0.5 or hi2 < 0.5)
    # DIRECTION IS PART OF THE VERDICT. A two-sided test is registered, so
    # a cell can clear the corrected alpha by ranking actives BELOW
    # inactives -- and four do. One boolean cannot say which happened.
    out["clears_corrected_alpha"] = bool(lo > 0.5 or hi < 0.5)
    out["direction"] = ("above" if metrics["auc"] > 0.5
                        else "below" if metrics["auc"] < 0.5 else "at")
    return out


def complement_cells(release: Path, source: str,
                     assay_types: set[str] | None,
                     labels) -> dict[str, tuple[int, int]]:
    """Cells whose membership is another protocol's leftovers.

    WHY THIS IS DERIVED AND NOT A LIST OF VIRUS IDS. `assemble.edge_key`
    does not include `assay_type` -- the schema says so and says that
    making it one "would split every merged edge" -- so when one source
    screens one virus under two protocols, the two measurements are one
    fact and they merge. The protocol that loses the merge does not come
    out small and representative. It comes out as exactly the compounds
    the winner never screened, and an AUC on that set measures the
    structural character of what the other library missed.

    The signature is disjointness: SHARING NO COMPOUNDS AT ALL with the
    sibling cell is not what a thin screen looks like, it is what a
    complement looks like. Measured on the 2026-11-17 release, the
    SARS-CoV-2 entry cell holds 196 compounds and the CPE cell 7,427,
    and the overlap is 0 of 196.

    A hard-coded `if virus == "NCBITaxon:2697049"` would be a check
    written with the same literal it is checking, and would say nothing
    the next time a second protocol is added. This asks the release.

    The smaller side is the one refused: the larger cell's population was
    not defined by the leftovers of anything. Ties refuse both, which is
    the conservative direction.
    """
    siblings = {k.split(":", 1)[1] for k in labels.stats
                if k.startswith("other_assay_type:")} - {"none"}
    siblings -= assay_types or set()
    if not siblings:
        return {}
    other = build_labels(release, sources={source}, assay_types=siblings)
    out: dict[str, tuple[int, int]] = {}
    for virus in set(labels.positives) | set(labels.negatives):
        mine = (set(labels.positives.get(virus, ()))
                | set(labels.negatives.get(virus, ())))
        theirs = (set(other.positives.get(virus, ()))
                  | set(other.negatives.get(virus, ())))
        if not mine or not theirs or (mine & theirs):
            continue
        if len(theirs) >= len(mine):
            out[virus] = (len(mine), len(theirs))
    return out


def run_population(release: Path, g: Graph, props: dict, source: str,
                   fps: dict, pubs: dict, control: str = "publication",
                   k_map: dict[str, int] | None = None,
                   strength_matched: bool = True,
                   assay_types: set[str] | None = None,
                   z: float = Z_CORRECTED,
                   z_alternate: float | None = None) -> dict:
    # SOURCE IS NOT ENOUGH TO MAKE A CELL SINGLE-SCREEN, and until
    # 2026-11-17 this function assumed it was. `infores:ncats-opendata` now
    # publishes two protocols against SARS-CoV-2 -- a cytopathic-effect
    # screen and a pseudotyped-particle entry screen -- and grouped by
    # virus alone they merged into ONE cell of 149 positives. The merge is
    # exactly what doc2/registrations/prereg_ncats_entry_assays.md registers against: a
    # plate-matched cell is worth having because one protocol over one
    # concentration range produced both classes, and two protocols pooled
    # destroy that while every count still looks reasonable.
    labels = build_labels(release, sources={source}, assay_types=assay_types)
    complement = complement_cells(release, source, assay_types, labels)
    degrees = degree_ranking(g)
    k_map = k_map or {}
    out: dict[str, dict] = {}
    for virus in sorted(set(labels.positives) | set(labels.negatives)):
        pos = labels.positives.get(virus, set())
        neg = labels.negatives.get(virus, set())
        readable = len(pos) >= MIN_CLASS and len(neg) >= MIN_CLASS
        cell: dict = {"n_positives": len(pos), "n_negatives": len(neg),
                      "readable": readable}
        if virus in complement:
            # REFUSED FOR ITS MEMBERSHIP, NOT ITS SIZE, and refused here
            # rather than in the reader's head -- a cell this shape clears
            # MIN_CLASS and would print an AUC with a plausible interval.
            mine, theirs = complement[virus]
            cell["readable"] = False
            cell["excluded_reason"] = (
                f"complement cell: all {mine} compounds are absent from "
                f"this source's other protocol for this virus ({theirs} "
                f"compounds), so this is what that screen did not cover "
                f"and not a sample of this one")
            out[virus] = cell
            continue
        if not readable:
            # NAMED, NOT DROPPED. A powerless cell removed silently is the
            # row this project has had to go back and label twice.
            cell["excluded_reason"] = (
                f"fewer than {MIN_CLASS} in a class")
            out[virus] = cell
            continue
        candidates = sorted(pos | neg)
        references = sorted(pos)
        if control == "publication":
            variants = (
                ("structure_max", 1, True, 0),
                ("structure_max_no_pub_exclusion", 1, False, 0),
                (f"structure_top{TOP_K}_mean", TOP_K, True, 0),
            )
        elif control == "knn":
            # The registered curve, plus the "both controls" secondary that
            # keeps the two runs comparable. `structure_max` is the primary
            # and is k = K_PRIMARY, named so the comparator code and the
            # report read the same key under either control.
            variants = tuple(
                (f"structure_max_k{k}", 1, False, k) for k in K_CURVE
            ) + (
                ("structure_max", 1, False, K_PRIMARY),
                (f"structure_max_k{K_PRIMARY}_both", 1, True, K_PRIMARY),
            )
        else:
            # STRENGTH-MATCHED. k is this cell's own, and K_PRIMARY rides
            # along as the companion column so every cell's movement is
            # readable in one table -- which the registration requires,
            # because P3 and P5 cannot be read without knowing which cells
            # got a sterner control and which a weaker one.
            k_cell = k_map.get(virus)
            if k_cell is None:
                cell["excluded_reason"] = (
                    "no k_cell: this cell is not readable in the "
                    "publication run, so there is no strength to match")
                cell["readable"] = False
                out[virus] = cell
                continue
            cell["k_cell"] = k_cell
            cell["strength_matched"] = strength_matched
            variants = (
                ("structure_max", 1, False, k_cell),
                (f"structure_max_k{K_PRIMARY}", 1, False, K_PRIMARY),
            )
        for name, top_k, disjoint, knn in variants:
            scores, counts = score_cell(candidates, references, fps, props,
                                        pubs, top_k=top_k,
                                        publication_disjoint=disjoint,
                                        positives=frozenset(pos),
                                        knn=knn)
            m = corrected(evaluate_against_negatives(scores, pos, neg),
                          z=z, also=z_alternate)
            m["counts"] = counts
            cell[name] = m
        cell["degree"] = corrected(
            evaluate_against_negatives(
                {d: v for d, v in degrees.items() if d in pos or d in neg},
                pos, neg), z=z, also=z_alternate)
        out[virus] = cell
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    # TWO REGISTRATIONS, TWO FILES, AND THE DEFAULT IS THE OLDER ONE.
    # doc2/registrations/prereg_structure_scorer.md's run is the comparator for
    # doc2/registrations/prereg_structure_scorer_symmetric.md and its figures are quoted
    # in a verdict and pinned by tests, so `--control knn` writes its own
    # file and the default branch must stay byte-identical.
    ap.add_argument("--control", choices=CONTROLS, default="publication",
                    help="publication: the first registration's exclusion "
                         "(3). knn: the symmetric replacement (3').")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)
    if args.out is None:
        args.out = ROOT / {
            "publication": "data/results/structure_scorer.json",
            "knn": "data/results/structure_scorer_symmetric.json",
            "strength": "data/results/structure_scorer_strength_matched.json",
        }[args.control]

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    doc = {
        "publication": "doc2/registrations/prereg_structure_scorer.md",
        "knn": "doc2/registrations/prereg_structure_scorer_symmetric.md",
        "strength": "doc2/registrations/prereg_structure_scorer_strength_matched.md",
    }[args.control]
    # The publication-disjoint run, for the side-by-side column P1 and P3
    # are about. Absent is not fatal: the column reads `--`, which is what
    # that mark is for.
    prior = None
    if args.control in ("knn", "strength"):
        prior_path = ROOT / "data/results/structure_scorer.json"
        if prior_path.exists():
            prior = json.loads(prior_path.read_text())
        else:
            print("! no publication-disjoint run on disk; the 'prior' "
                  "column will be empty and P1/P3 cannot be read from "
                  "this output alone")
    k_map: dict[str, int] = {}
    if args.control == "strength":
        if prior is None:
            print("REFUSED: the strength-matched control derives every k "
                  "from data/results/structure_scorer.json and that file "
                  "is absent. There is no strength to match.")
            return 2
        k_map = k_cells_from_prior(prior)
        print(f"k_cell per cell, from the publication run's counters: "
              f"{min(k_map.values())} to {max(k_map.values())} over "
              f"{len(k_map)} cells")
    print("STRUCTURE-DERIVED SCORER. Protocol fixed in advance by")
    print(f"{doc}, signed before the code that runs it.")
    print(f"release {args.release.name}   control {args.control}")
    if args.control == "knn":
        print(f"exclusion (3'): the k nearest references, k in {K_CURVE}, "
              f"primary k={K_PRIMARY}")
    elif args.control == "strength":
        print("exclusion (3'): the k_cell nearest references, k_cell = "
              "ceil(max(pub/pos, pub/neg))")
    print()

    schema = load_schema()
    # SKIP THE HELD-OUT PREDICATES, and this script was the only
    # Graph.load caller in the repository that did not.
    #
    # `baselines.degree_ranking` says of itself "OVER THE TRAVERSABLE
    # GRAPH, NOT THE WHOLE ONE (fix-list L7) ... a null hypothesis that
    # reads the held-out answers is not a null hypothesis". It delivers
    # that only if its caller asks for it: `Graph.load` takes
    # `skip_predicates` and defaults it to the EMPTY SET, so without the
    # argument `degree_unlabelled` is `degree` and the scorer's degree
    # comparator counted the very label edges it was being scored
    # against. Twelve other scripts pass `schema.held_out_predicates()`;
    # this one never did, so L7 was fixed everywhere except in the three
    # studies whose entire argument is "a structure column that does not
    # beat its own degree column is not evidence about chemistry".
    #
    # A CONTRACT A CALLER HAS TO OPT INTO is the same defect shape as a
    # constant declared and never passed: the docstring describes the
    # behaviour, the default silently does not provide it, and nothing
    # fails. tests/test_degree_baseline_skips_labels.py now checks every
    # caller rather than this one.
    g = Graph.load(args.release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    props = small_molecule_props(args.release)
    pubs = {k: frozenset(v) for k, v in
            label_publication_index(args.release, {ACTIVITY_PREDICATE}).items()}
    print(f"{len(pubs):,} compounds carry a PMID on a label edge")

    # EVERY LABELLED COMPOUND IN EITHER POPULATION, FINGERPRINTED ONCE.
    # The parse rate is printed because registered prediction P3 -- floor
    # exactly 0.5000 everywhere -- depends on it, and a SMILES is not a
    # fingerprint.
    labelled: set[str] = set()
    for source in (PRIMARY_SOURCE, SECONDARY_SOURCE):
        lab = build_labels(args.release, sources={source})
        for v in set(lab.positives) | set(lab.negatives):
            labelled |= lab.positives.get(v, set())
            labelled |= lab.negatives.get(v, set())
    print(f"{len(labelled):,} labelled compounds across both populations")
    fps, fp_stats = fingerprints(
        {i: (props.get(i, {}).get("smiles") or "") for i in sorted(labelled)})
    print(f"fingerprinted {fp_stats['fingerprinted']:,}, "
          f"unparseable {fp_stats['unparseable']:,}, "
          f"no SMILES {fp_stats['no_smiles']:,}")
    if fp_stats["unparseable"] or fp_stats["no_smiles"]:
        print("  ! P3 PREDICTED A FLOOR OF EXACTLY 0.5000 IN EVERY CELL, and")
        print("    an unfingerprintable compound takes the imputed floor, so")
        print("    that prediction is already false. Reported, not patched.")

    payload = {
        "measured_by": "scripts/structure_scorer.py",
        "registered_by": doc,
        "control": args.control,
        "fingerprint": "ECFP4 radius 2, 2048 bits (kgav.similarity)",
        "statistic": "max Tanimoto to admissible known actives",
        **({"k_curve": list(K_CURVE), "k_primary": K_PRIMARY}
           if args.control == "knn" else {}),
        **({"k_cells": k_map, "k_companion": K_PRIMARY}
           if args.control == "strength" else {}),
        "family_size": FAMILY,
        "alpha": ALPHA,
        "z_corrected": Z_CORRECTED,
        "plate_family_size": PLATE_FAMILY,
        "plate_alpha": PLATE_ALPHA,
        "plate_z_corrected": PLATE_Z_CORRECTED,
        "fingerprint_stats": dict(fp_stats),
        "n_labelled": len(labelled),
    }
    for key, source, assay_types in POPULATIONS:
        print(f"\n{'=' * 68}\n{key}  ({source}"
              f"{', ' + '/'.join(sorted(assay_types)) if assay_types else ''})"
              f"\n{'=' * 68}")
        # THE PLATE-MATCHED CELL IS NOT STRENGTH-MATCHED AND THE OUTPUT
        # MUST SAY SO. Publication-disjointness removed every admissible
        # reference there -- 88.00 per positive and 89.00 per negative out
        # of 88 -- so the matching k is 89 and leaves -1. It runs at the
        # ChEMBL cell's k_cell instead, which makes the two populations
        # comparable at one handicap without pretending to match a
        # strength that was total.
        matched = args.control != "strength" or key == "primary_chembl"
        # A PLATE POPULATION IS TESTED AT THE PLATE FAMILY'S ALPHA. The
        # ChEMBL family's z is kept as the `also` level for the plate
        # cells, because prereg_structure_scorer.md quotes its secondary
        # cell at that level and this script must not pick one
        # registration over the other without printing both.
        plate = assay_types is not None
        cells = run_population(args.release, g, props, source, fps, pubs,
                               control=args.control, k_map=k_map,
                               strength_matched=matched,
                               assay_types=assay_types,
                               z=PLATE_Z_CORRECTED if plate else Z_CORRECTED,
                               z_alternate=Z_CORRECTED if plate else None)
        payload[key] = cells
        if args.control == "strength" and not matched:
            print("  NOT STRENGTH-MATCHED: the publication control removed "
                  "every reference in this")
            print("  cell, so no symmetric control can match it. Run at the "
                  "ChEMBL cell's k_cell")
            print("  for comparability; read it as that and not as a "
                  "strength match.")
        report(cells, control=args.control, prior=prior, population=key,
               plate=plate)
        if args.control in ("knn", "strength"):
            knn_invariant(cells, k_map if args.control == "strength" else None)

    # THE PLATE FAMILY IS REGISTERED TOO, and at three rather than the
    # four the registration was signed for. Reported, not recomputed.
    n_plate = sum(1 for key, _, at in POPULATIONS if at is not None
                  for c in payload[key].values() if c.get("readable"))
    print(f"\n{'=' * 68}\nplate-matched family: {n_plate} readable cells, "
          f"registered at {PLATE_FAMILY}")
    print(f"  alpha = 0.05/{PLATE_FAMILY} = {PLATE_ALPHA:.6f}, "
          f"z = {PLATE_Z_CORRECTED}")
    if n_plate != PLATE_FAMILY:
        print(f"! THE REGISTERED PLATE FAMILY IS {PLATE_FAMILY} AND THIS "
              f"RELEASE GIVES {n_plate}.")
        print("  doc2/registrations/prereg_ncats_entry_assays.md, second amendment. The "
              "family was fixed before")
        print("  the run; this line is the finding, not a new alpha.")

    # THE FAMILY SIZE IS REGISTERED, SO A MISMATCH IS A FINDING.
    n_readable = sum(1 for c in payload["primary_chembl"].values()
                     if c.get("readable"))
    if n_readable != FAMILY:
        print(f"\n! THE REGISTERED FAMILY IS {FAMILY} CELLS AND THIS RELEASE "
              f"GIVES {n_readable}.")
        print("  The corrected alpha in the registration is 0.05/"
              f"{FAMILY}; on this release it would be 0.05/{n_readable}.")
        print("  Reported, not silently recomputed: the family was fixed "
              "before the run.")

    write_results(args.out, payload, args.release, schema.version)
    print(f"\nwrote {args.out}")
    return 0


def knn_invariant(cells: dict, k_map: dict[str, int] | None = None) -> None:
    """P2: every candidate loses exactly k references to (3').

    REGISTERED AS AN INVARIANT, NOT AN OUTCOME, so it is checked where it
    can be observed rather than inferred from an AUC. The prior run's
    per-class counters had to be added after the fact to test a hypothesis
    its own output could not see; this prints the thing the prediction is
    about.
    """
    k_map = k_map or {}
    label = "k_cell" if k_map else f"k={K_PRIMARY}"
    print(f"\n  P2 -- exclusion (3') removes exactly k, per class ({label}):")
    print(f"    {'virus':<20} {'k':>4} {'knn/pos':>9} {'knn/neg':>9} "
          f"{'short':>6} {'skel/pos':>9} {'skel/neg':>9}  invariant")
    bad = []
    for virus, c in sorted(cells.items()):
        if not c.get("readable"):
            continue
        k = c["structure_max"]["counts"]
        # PER CELL WHERE THERE IS A PER-CELL k. Checking against K_PRIMARY
        # under the strength control would compare 10 against 33 and report
        # the invariant violated in thirteen of fourteen cells -- a check
        # that fails for its own reason says nothing about the thing it
        # guards.
        kk = k_map.get(virus, K_PRIMARY)
        exp = kk * k["candidates_pos"], kk * k["candidates_neg"]
        got = k["excluded_knn_pos"], k["excluded_knn_neg"]
        ok = got == exp and not k["knn_short"]
        if not ok:
            bad.append(virus)
        per_pos = got[0] / max(k["candidates_pos"], 1)
        per_neg = got[1] / max(k["candidates_neg"], 1)
        print(f"    {virus:<20} {kk:>4} {per_pos:>9.2f} {per_neg:>9.2f} "
              f"{k['knn_short']:>6,} "
              f"{k['excluded_skeleton_pos'] / max(k['candidates_pos'], 1):>9.2f} "
              f"{k['excluded_skeleton_neg'] / max(k['candidates_neg'], 1):>9.2f}"
              f"  {'HELD' if ok else 'VIOLATED'}")
    if bad:
        print(f"    ! P2 VIOLATED in {len(bad)} cell(s): {bad}. The control "
              f"is not symmetric and no AUC above is interpretable.")
    else:
        print("    P2 held in every cell: the handicap is k for both "
              "classes, by construction.")


def num(x: float | None, width: int = 7) -> str:
    """A number, or a mark that says there is none. NEVER 0.0 for None."""
    return f"{x:.4f}" if x is not None else "--"


def fmt(later: float | None, base: float) -> str:
    return f"{later - base:+.4f}" if later is not None else "--"


def report(cells: dict, control: str = "publication",
           prior: dict | None = None, population: str = "primary_chembl",
           plate: bool = False) -> None:
    """One table per control, because the columns are not the same question.

    The publication control's interesting companion is its own leak -- the
    same scorer with exclusion (3) removed. The knn control's is the CURVE
    and the prior run's figure, because the registered predictions are about
    decay in k (P4) and about movement relative to the publication-disjoint
    result (P1, P3). A shared formatter printing one of those under the
    other's heading is how `structure_max_no_pub_exclusion` came to be read
    for a run that never computed it.
    """
    if control == "publication":
        print(f"  {'virus':<20} {'pos':>6} {'neg':>7} {'AUC':>7} "
              f"{'[corrected 95%]':>18} {'skill':>7} {'leak':>7} "
              f"{'degree':>7}")
    elif control == "strength":
        # k_cell IS A COLUMN, which the registration requires: P3 and P5
        # cannot be read without knowing which cells got a sterner control
        # than k=10 and which got a weaker one.
        print(f"  {'virus':<20} {'pos':>6} {'neg':>7} {'k':>4} "
              f"{'AUC':>7} {'[corrected 95%]':>18} {'skill':>7} "
              f"{'k=10':>7} {'prior':>7} {'degree':>7}")
    else:
        print(f"  {'virus':<20} {'pos':>6} {'neg':>7} "
              f"{'k=1':>7} {'k=3':>7} {'k=10':>7} {'[corrected 95%]':>18} "
              f"{'skill':>7} {'both':>7} {'prior':>7} {'degree':>7}")
    for virus, c in sorted(cells.items()):
        if not c.get("readable"):
            print(f"  {virus:<20} {c['n_positives']:>6,} "
                  f"{c['n_negatives']:>7,}   NOT READABLE "
                  f"({c['excluded_reason']})")
            continue
        m = c["structure_max"]
        d = c["degree"]
        if not m.get("measurable"):
            print(f"  {virus:<20} {m['n_pos']:>6,} {m['n_neg']:>7,} "
                  f"  NOT MEASURABLE — {m['refusal']}")
            continue
        ci = f"[{m['auc_ci_lo_corrected']:.4f}, {m['auc_ci_hi_corrected']:.4f}]"
        # DIRECTION IN THE MARK, because a two-sided test is registered and
        # cells clear it by ranking actives BELOW inactives.
        mark = ("^" if m["clears_corrected_alpha"] and m["direction"] == "above"
                else "v" if m["clears_corrected_alpha"] else " ")
        # NO `or 0.0` ANYWHERE BELOW. That idiom printed a vacuous cell's
        # auc_skill of None as a confident 0.0000.
        if control == "strength":
            at10 = c.get(f"structure_max_k{K_PRIMARY}", {}).get("auc")
            # THE PRIOR IS THIS POPULATION'S OWN PRIOR. Until 2026-11-17
            # this read primary_chembl first and fell back to
            # secondary_ncats only when ChEMBL had no cell for the virus
            # -- so every plate cell printed the CHEMBL figure under a
            # column headed "the publication-disjoint figure this
            # replaces". Measured: the CPE cell showed prior 0.7718 and
            # the two entry cells 0.6851 and 0.5098, which are the ChEMBL
            # SARS-CoV-2, MERS-CoV and SARS-CoV AUCs. Those cells have no
            # prior at all -- their publication-control run REFUSED -- and
            # a reader comparing 0.6062 against "0.7718" would be reading
            # a different population's result as this one's history.
            was = (prior.get(population, {}).get(virus, {})
                   .get("structure_max", {}).get("auc") if prior else None)
            print(f"  {virus:<20} {m['n_pos']:>6,} {m['n_neg']:>7,} "
                  f"{c['k_cell']:>4} {m['auc']:>7.4f}{mark}{ci:>18} "
                  f"{num(m.get('auc_skill')):>7} {num(at10):>7} "
                  f"{num(was):>7} {num(d.get('auc')):>7}")
        elif control == "publication":
            leak = c["structure_max_no_pub_exclusion"]
            print(f"  {virus:<20} {m['n_pos']:>6,} {m['n_neg']:>7,} "
                  f"{m['auc']:>7.4f}{mark}{ci:>18} "
                  f"{num(m.get('auc_skill')):>7} "
                  f"{fmt(leak.get('auc'), m['auc']):>7} "
                  f"{num(d.get('auc')):>7}")
        else:
            curve = [c.get(f"structure_max_k{k}", {}).get("auc")
                     for k in K_CURVE]
            both = c.get(f"structure_max_k{K_PRIMARY}_both", {}).get("auc")
            # THE PRIOR IS THIS POPULATION'S OWN PRIOR. Until 2026-11-17
            # this read primary_chembl first and fell back to
            # secondary_ncats only when ChEMBL had no cell for the virus
            # -- so every plate cell printed the CHEMBL figure under a
            # column headed "the publication-disjoint figure this
            # replaces". Measured: the CPE cell showed prior 0.7718 and
            # the two entry cells 0.6851 and 0.5098, which are the ChEMBL
            # SARS-CoV-2, MERS-CoV and SARS-CoV AUCs. Those cells have no
            # prior at all -- their publication-control run REFUSED -- and
            # a reader comparing 0.6062 against "0.7718" would be reading
            # a different population's result as this one's history.
            was = (prior.get(population, {}).get(virus, {})
                   .get("structure_max", {}).get("auc") if prior else None)
            print(f"  {virus:<20} {m['n_pos']:>6,} {m['n_neg']:>7,} "
                  f"{num(curve[0]):>7} {num(curve[1]):>7} "
                  f"{m['auc']:>7.4f}{mark}{ci:>18} "
                  f"{num(m.get('auc_skill')):>7} {num(both):>7} "
                  f"{num(was):>7} {num(d.get('auc')):>7}")
    if plate:
        # BOTH REGISTERED LEVELS, SIDE BY SIDE. The table above is at the
        # plate family's alpha, which is what this population is
        # registered at. The ChEMBL family's level is printed here rather
        # than left to a reader who may be holding the other registration.
        print(f"    at the plate family's alpha (0.05/{PLATE_FAMILY}, "
              f"z = {PLATE_Z_CORRECTED}) above; at the ChEMBL family's "
              f"(0.05/{FAMILY}, z = {Z_CORRECTED}):")
        for virus, c in sorted(cells.items()):
            m = c.get("structure_max", {}) if c.get("readable") else {}
            if "auc_ci_lo_alternate" not in m:
                continue
            verdict = ("clears" if m["clears_alternate_alpha"]
                       else "does NOT clear")
            agree = "" if (m["clears_alternate_alpha"]
                           == m["clears_corrected_alpha"]) else \
                "   <== THE TWO LEVELS DISAGREE FOR THIS CELL"
            print(f"      {virus:<20} {m['auc']:.4f}  "
                  f"[{m['auc_ci_lo_alternate']:.4f}, "
                  f"{m['auc_ci_hi_alternate']:.4f}]  {verdict}{agree}")
    print("  ^ clears the corrected alpha ABOVE chance, v clears it BELOW.")
    if control == "strength":
        print("    k is this cell's own k_cell = ceil(max(pub/pos, "
              "pub/neg)); 'k=10' is the")
        print("    symmetric study's companion figure and 'prior' the "
              "publication-disjoint one.")
    elif control == "publication":
        print("    'leak' is the no-publication-exclusion variant minus the "
              "primary: positive")
        print("    means same-paper actives were inflating it.")
    else:
        print("    k=10 is the primary. 'both' keeps the publication "
              "exclusion as well;")
        print("    'prior' is the publication-disjoint figure this replaces.")
    print("    '--' is a quantity that does not exist, not a zero.")


if __name__ == "__main__":
    raise SystemExit(main())
