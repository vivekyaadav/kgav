#!/usr/bin/env python3
"""The structure-derived scorer, exactly as registered.

PROTOCOL: docs/prereg_structure_scorer.md, signed 2026-11-11 before this
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
MIN_CLASS = 10              # calibration.MIN_REACHED_POSITIVES
FAMILY = 14                 # the readable ChEMBL cells
ALPHA = 0.05 / FAMILY
Z_CORRECTED = 2.9137        # two-sided, for ALPHA
TOP_K = 3                   # the one registered alternative statistic

# docs/prereg_structure_scorer_symmetric.md, signed 2026-11-11.
# K_PRIMARY is derived, not chosen: JEV holds 18 references after
# self-exclusion, the smallest in the family, and k must leave every cell
# non-empty or the study repeats the prior run's P5 -- a cell that printed
# 0.5000 with an interval because nothing had been ranked. k=30, which the
# prior verdict suggested, leaves -12.
K_CURVE = (1, 3, 10)
K_PRIMARY = 10
CONTROLS = ("publication", "knn")


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
        # EXCLUSION (3'), docs/prereg_structure_scorer_symmetric.md. P2 is
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


def corrected(metrics: dict) -> dict:
    """The same AUC at the registered corrected level, and a refusal.

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
    lo, hi = auc_interval(metrics["auc"], metrics["n_pos"], metrics["n_neg"],
                          z=Z_CORRECTED)
    out["auc_ci_lo_corrected"] = lo
    out["auc_ci_hi_corrected"] = hi
    # DIRECTION IS PART OF THE VERDICT. A two-sided test is registered, so
    # a cell can clear the corrected alpha by ranking actives BELOW
    # inactives -- and four do. One boolean cannot say which happened.
    out["clears_corrected_alpha"] = bool(lo > 0.5 or hi < 0.5)
    out["direction"] = ("above" if metrics["auc"] > 0.5
                        else "below" if metrics["auc"] < 0.5 else "at")
    return out


def run_population(release: Path, g: Graph, props: dict, source: str,
                   fps: dict, pubs: dict, control: str = "publication") -> dict:
    labels = build_labels(release, sources={source})
    degrees = degree_ranking(g)
    out: dict[str, dict] = {}
    for virus in sorted(set(labels.positives) | set(labels.negatives)):
        pos = labels.positives.get(virus, set())
        neg = labels.negatives.get(virus, set())
        readable = len(pos) >= MIN_CLASS and len(neg) >= MIN_CLASS
        cell: dict = {"n_positives": len(pos), "n_negatives": len(neg),
                      "readable": readable}
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
        else:
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
        for name, top_k, disjoint, knn in variants:
            scores, counts = score_cell(candidates, references, fps, props,
                                        pubs, top_k=top_k,
                                        publication_disjoint=disjoint,
                                        positives=frozenset(pos),
                                        knn=knn)
            m = corrected(evaluate_against_negatives(scores, pos, neg))
            m["counts"] = counts
            cell[name] = m
        cell["degree"] = corrected(
            evaluate_against_negatives(
                {d: v for d, v in degrees.items() if d in pos or d in neg},
                pos, neg))
        out[virus] = cell
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    # TWO REGISTRATIONS, TWO FILES, AND THE DEFAULT IS THE OLDER ONE.
    # docs/prereg_structure_scorer.md's run is the comparator for
    # docs/prereg_structure_scorer_symmetric.md and its figures are quoted
    # in a verdict and pinned by tests, so `--control knn` writes its own
    # file and the default branch must stay byte-identical.
    ap.add_argument("--control", choices=CONTROLS, default="publication",
                    help="publication: the first registration's exclusion "
                         "(3). knn: the symmetric replacement (3').")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)
    if args.out is None:
        args.out = ROOT / ("data/results/structure_scorer.json"
                           if args.control == "publication"
                           else "data/results/structure_scorer_symmetric.json")

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    doc = ("docs/prereg_structure_scorer.md" if args.control == "publication"
           else "docs/prereg_structure_scorer_symmetric.md")
    # The publication-disjoint run, for the side-by-side column P1 and P3
    # are about. Absent is not fatal: the column reads `--`, which is what
    # that mark is for.
    prior = None
    if args.control == "knn":
        prior_path = ROOT / "data/results/structure_scorer.json"
        if prior_path.exists():
            prior = json.loads(prior_path.read_text())
        else:
            print("! no publication-disjoint run on disk; the 'prior' "
                  "column will be empty and P1/P3 cannot be read from "
                  "this output alone")
    print("STRUCTURE-DERIVED SCORER. Protocol fixed in advance by")
    print(f"{doc}, signed before the code that runs it.")
    print(f"release {args.release.name}   control {args.control}")
    if args.control == "knn":
        print(f"exclusion (3'): the k nearest references, k in {K_CURVE}, "
              f"primary k={K_PRIMARY}")
    print()

    schema = load_schema()
    g = Graph.load(args.release, symmetry=schema.symmetry())
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
        "family_size": FAMILY,
        "alpha": ALPHA,
        "z_corrected": Z_CORRECTED,
        "fingerprint_stats": dict(fp_stats),
        "n_labelled": len(labelled),
    }
    for key, source in (("primary_chembl", PRIMARY_SOURCE),
                        ("secondary_ncats", SECONDARY_SOURCE)):
        print(f"\n{'=' * 68}\n{key}  ({source})\n{'=' * 68}")
        cells = run_population(args.release, g, props, source, fps, pubs,
                               control=args.control)
        payload[key] = cells
        report(cells, control=args.control, prior=prior)
        if args.control == "knn":
            knn_invariant(cells)

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


def knn_invariant(cells: dict) -> None:
    """P2: every candidate loses exactly k references to (3').

    REGISTERED AS AN INVARIANT, NOT AN OUTCOME, so it is checked where it
    can be observed rather than inferred from an AUC. The prior run's
    per-class counters had to be added after the fact to test a hypothesis
    its own output could not see; this prints the thing the prediction is
    about.
    """
    print(f"\n  P2 -- exclusion (3') removes exactly k, per class "
          f"(k={K_PRIMARY}):")
    print(f"    {'virus':<20} {'knn/pos':>9} {'knn/neg':>9} {'short':>6} "
          f"{'skel/pos':>9} {'skel/neg':>9}  invariant")
    bad = []
    for virus, c in sorted(cells.items()):
        if not c.get("readable"):
            continue
        k = c["structure_max"]["counts"]
        exp = K_PRIMARY * k["candidates_pos"], K_PRIMARY * k["candidates_neg"]
        got = k["excluded_knn_pos"], k["excluded_knn_neg"]
        ok = got == exp and not k["knn_short"]
        if not ok:
            bad.append(virus)
        per_pos = got[0] / max(k["candidates_pos"], 1)
        per_neg = got[1] / max(k["candidates_neg"], 1)
        print(f"    {virus:<20} {per_pos:>9.2f} {per_neg:>9.2f} "
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
           prior: dict | None = None) -> None:
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
        if control == "publication":
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
            was = None
            if prior:
                was = (prior.get("primary_chembl", {}).get(virus, {})
                       .get("structure_max", {}).get("auc"))
                if was is None:
                    was = (prior.get("secondary_ncats", {}).get(virus, {})
                           .get("structure_max", {}).get("auc"))
            print(f"  {virus:<20} {m['n_pos']:>6,} {m['n_neg']:>7,} "
                  f"{num(curve[0]):>7} {num(curve[1]):>7} "
                  f"{m['auc']:>7.4f}{mark}{ci:>18} "
                  f"{num(m.get('auc_skill')):>7} {num(both):>7} "
                  f"{num(was):>7} {num(d.get('auc')):>7}")
    print("  ^ clears the corrected alpha ABOVE chance, v clears it BELOW.")
    if control == "publication":
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
