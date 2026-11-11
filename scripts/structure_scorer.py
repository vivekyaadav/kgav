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
               positives: frozenset[str] = frozenset()) -> tuple[dict, dict]:
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
        if publication_disjoint and c_pubs:
            for i, pp in enumerate(ref_pubs):
                if i not in drop and (pp & c_pubs):
                    drop.add(i)
                    counts["excluded_publication"] += 1
                    counts["excluded_publication_pos" if is_pos
                           else "excluded_publication_neg"] += 1
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
                   fps: dict, pubs: dict) -> dict:
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
        for name, top_k, disjoint in (
            ("structure_max", 1, True),
            ("structure_max_no_pub_exclusion", 1, False),
            (f"structure_top{TOP_K}_mean", TOP_K, True),
        ):
            scores, counts = score_cell(candidates, references, fps, props,
                                        pubs, top_k=top_k,
                                        publication_disjoint=disjoint,
                                        positives=frozenset(pos))
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
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data/results/structure_scorer.json")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    print("STRUCTURE-DERIVED SCORER. Protocol fixed in advance by")
    print("docs/prereg_structure_scorer.md, signed before this file existed.")
    print(f"release {args.release.name}\n")

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
        "registered_by": "docs/prereg_structure_scorer.md",
        "fingerprint": "ECFP4 radius 2, 2048 bits (kgav.similarity)",
        "statistic": "max Tanimoto to admissible known actives",
        "family_size": FAMILY,
        "alpha": ALPHA,
        "z_corrected": Z_CORRECTED,
        "fingerprint_stats": dict(fp_stats),
        "n_labelled": len(labelled),
    }
    for key, source in (("primary_chembl", PRIMARY_SOURCE),
                        ("secondary_ncats", SECONDARY_SOURCE)):
        print(f"\n{'=' * 68}\n{key}  ({source})\n{'=' * 68}")
        cells = run_population(args.release, g, props, source, fps, pubs)
        payload[key] = cells
        report(cells)

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


def num(x: float | None, width: int = 7) -> str:
    """A number, or a mark that says there is none. NEVER 0.0 for None."""
    return f"{x:.4f}" if x is not None else "--"


def fmt(later: float | None, base: float) -> str:
    return f"{later - base:+.4f}" if later is not None else "--"


def report(cells: dict) -> None:
    head = (f"  {'virus':<20} {'pos':>6} {'neg':>7} {'AUC':>7} "
            f"{'[corrected 95%]':>18} {'skill':>7} {'leak':>7} {'degree':>7}")
    print(head)
    for virus, c in sorted(cells.items()):
        if not c.get("readable"):
            print(f"  {virus:<20} {c['n_positives']:>6,} "
                  f"{c['n_negatives']:>7,}   NOT READABLE "
                  f"({c['excluded_reason']})")
            continue
        m = c["structure_max"]
        leak = c["structure_max_no_pub_exclusion"]
        d = c["degree"]
        if not m.get("measurable"):
            print(f"  {virus:<20} {m['n_pos']:>6,} {m['n_neg']:>7,} "
                  f"  NOT MEASURABLE — {m['refusal']}")
            continue
        ci = f"[{m['auc_ci_lo_corrected']:.4f}, {m['auc_ci_hi_corrected']:.4f}]"
        # DIRECTION IN THE MARK, because a two-sided test is registered and
        # four cells clear it by ranking actives BELOW inactives.
        mark = ("^" if m["clears_corrected_alpha"] and m["direction"] == "above"
                else "v" if m["clears_corrected_alpha"] else " ")
        # NO `or 0.0` ON EITHER OF THESE. That idiom is what printed the
        # vacuous cell's auc_skill of None as a confident 0.0000.
        delta = fmt(leak.get("auc"), m["auc"])
        print(f"  {virus:<20} {m['n_pos']:>6,} {m['n_neg']:>7,} "
              f"{m['auc']:>7.4f}{mark}{ci:>18} "
              f"{num(m.get('auc_skill')):>7} {delta:>7} "
              f"{num(d.get('auc')):>7}")
    print("  ^ clears the corrected alpha ABOVE chance, v clears it BELOW. "
          "'leak' is the")
    print("    no-publication-exclusion variant minus the primary: positive "
          "means same-paper")
    print("    actives were inflating it. '--' is a quantity that does not "
          "exist, not a zero.")


if __name__ == "__main__":
    raise SystemExit(main())
