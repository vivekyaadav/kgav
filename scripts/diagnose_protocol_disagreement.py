#!/usr/bin/env python3
"""Why one cell can anti-rank cross-sectionally and clear prospectively.

    make protocol-disagreement            # DENV-4, the cell that does
    python scripts/diagnose_protocol_disagreement.py --virus NCBITaxon:11060

THE QUESTION. DENV-4 is refused by scripts/repurposing_candidates.py --
its cross-sectional strength-matched AUC is 0.4792, below chance, which
means the scorer ranks known actives UNDER known inactives there -- and
it CLEARS prospectively at 0.7896. Those read as a contradiction and
are not one.

WHAT THIS MEASURES, in four steps that rule out the easy explanations:

  1. THE POPULATIONS. Are the two protocols even scoring the same
     compounds? For DENV-4 they are nested exactly: every prospective
     test positive is a cross-sectional positive, and every reference is
     too. So the difference is not which compounds, it is which ones act
     as REFERENCES and which are scored.

  2. ERA POOLING. The cross-sectional cell mixes pre- and post-cutoff
     compounds, and a pooled AUC can sit below both of its parts. So
     the published cell is split by era, under its own control, and the
     post-cutoff half is read on its own.

  3. THE SAME COMPOUNDS, THE OTHER PROTOCOL. The post-cutoff positives
     and negatives are exactly the prospective test set, so their
     cross-sectional AUC and their prospective AUC are directly
     comparable: identical compounds, identical labels, different
     reference set.

  4. THE MECHANISM. For each post-cutoff positive, how many of the k
     nearest references the control removes are its OWN era-mates.

WHAT IT IS NOT. It computes no new claim and amends no verdict. Both
registered results stand: they answer different questions, and this
reports which. Nothing here is written to data/results.
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

from kgav.labels import build_labels, evaluate_against_negatives
from kgav.pipeline import release_dir, train_dir
from kgav.similarity import fingerprints
from kgav import virus_register as vr

from structure_scorer import score_cell, small_molecule_props
from temporal_structure_scorer import CUTOFF

RESULTS = ROOT / "data" / "results"
RULE = "=" * 76

#: The cell the disagreement was found on. Any cell may be passed, and
#: one that does not disagree simply reports that it does not.
DEFAULT_VIRUS = "NCBITaxon:11070"


def _auc(scores: dict, pos: set, neg: set) -> float | None:
    if not pos or not neg:
        return None
    return evaluate_against_negatives(scores, set(pos), set(neg))["auc"]


def analyse(virus: str, release: Path, split: Path) -> dict:
    props = small_molecule_props(release)
    labels = json.loads((split / "TEST_LABELS.json").read_text())
    pre_p = set(labels["train"].get(virus, []))
    post_p = set(labels["test"].get(virus, []))
    post_n = set(labels["test_negatives"].get(virus, []))

    xs = build_labels(release)
    all_p = xs.positives.get(virus, set())
    all_n = xs.negatives.get(virus, set())
    pre_n = all_n - post_n

    k_cell = json.loads(
        (RESULTS / "structure_scorer_strength_matched.json").read_text()
    )["k_cells"].get(virus)
    if k_cell is None:
        raise SystemExit(f"{virus} has no k_cell: it is not readable in the "
                         "publication run, so there is no cross-sectional "
                         "verdict to disagree with")

    cand = sorted(all_p | all_n)
    fps, _ = fingerprints(
        {i: (props.get(i, {}).get("smiles") or "") for i in cand})

    # THE CROSS-SECTIONAL SETUP, exactly as structure_scorer runs it:
    # every positive is a reference, each candidate loses its k_cell
    # nearest, and the compound is never scored against itself.
    scores, _ = score_cell(cand, sorted(all_p), fps, props, {}, top_k=1,
                           publication_disjoint=False,
                           positives=frozenset(all_p), knn=k_cell)

    # THE MECHANISM: of the k nearest references the control removes
    # from a post-cutoff positive, how many are its own era-mates.
    from rdkit import DataStructs

    refs = [r for r in sorted(all_p) if r in fps]
    shares = []
    for c in sorted(post_p):
        if c not in fps:
            continue
        near = sorted(((DataStructs.TanimotoSimilarity(fps[c], fps[r]), r)
                       for r in refs if r != c), reverse=True)[:k_cell]
        if near:
            shares.append(sum(1 for _s, r in near if r in post_p) / len(near))

    temporal = json.loads(
        (RESULTS / "temporal_structure_scorer.json").read_text())
    cell = (temporal.get("cells") or {}).get(virus) or {}

    return {
        "virus": virus,
        "k_cell": k_cell,
        "cutoff": temporal.get("cutoff", CUTOFF),
        "counts": {
            "cross_sectional_positives": len(all_p),
            "cross_sectional_negatives": len(all_n),
            "pre_cutoff_positives": len(pre_p),
            "post_cutoff_positives": len(post_p),
            "post_cutoff_negatives": len(post_n),
            "pre_cutoff_negatives": len(pre_n),
        },
        "nested": {
            "test_positives_in_cross_sectional": len(post_p & all_p),
            "test_negatives_in_cross_sectional": len(post_n & all_n),
            "references_in_cross_sectional": len(pre_p & all_p),
            "train_test_overlap": len(pre_p & post_p),
        },
        "cross_sectional_by_era": {
            "all_vs_all": _auc(scores, all_p, all_n),
            "post_vs_post": _auc(scores, post_p, post_n),
            "pre_vs_pre": _auc(scores, pre_p, pre_n),
            "pre_pos_vs_post_neg": _auc(scores, pre_p, post_n),
        },
        "prospective": {
            "primary": (cell.get("primary") or {}).get("auc"),
            "symmetric_k10": (cell.get("symmetric_k10") or {}).get("auc"),
            "clears": virus in set(temporal.get("clears") or []),
        },
        "era_mate_share": {
            "n_compounds": len(shares),
            "median": statistics.median(shares) if shares else None,
            "mean": statistics.mean(shares) if shares else None,
        },
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--virus", default=DEFAULT_VIRUS)
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--split", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)
    split = args.split or train_dir(ROOT, CUTOFF)

    names = {"NCBITaxon:" + v["taxon"]: v["name"] for v in vr.load()["viruses"]}
    r = analyse(args.virus, args.release, split)
    label = names.get(args.virus, args.virus)

    print(RULE)
    print(f"{label} — cross-sectional against prospective")
    print(RULE)
    print()
    c, n = r["counts"], r["nested"]
    print("1. THE POPULATIONS")
    print(f"   cross-sectional : {c['cross_sectional_positives']} positive / "
          f"{c['cross_sectional_negatives']} negative")
    print(f"   prospective     : {c['pre_cutoff_positives']} references, "
          f"{c['post_cutoff_positives']}+/{c['post_cutoff_negatives']}- tested")
    print(f"   test positives that are cross-sectional positives: "
          f"{n['test_positives_in_cross_sectional']} of "
          f"{c['post_cutoff_positives']}")
    print(f"   references that are cross-sectional positives    : "
          f"{n['references_in_cross_sectional']} of "
          f"{c['pre_cutoff_positives']}")
    nested = (n["test_positives_in_cross_sectional"]
              == c["post_cutoff_positives"]
              and n["references_in_cross_sectional"]
              == c["pre_cutoff_positives"])
    print(f"   -> {'NESTED EXACTLY' if nested else 'NOT nested'}: the "
          "difference is not which compounds, it is which are references.")
    print()

    e = r["cross_sectional_by_era"]
    print(f"2. THE CROSS-SECTIONAL CELL, SPLIT BY ERA "
          f"(strength-matched, k_cell = {r['k_cell']})")
    for key, lab in (("all_vs_all", "all vs all  [the published cell]"),
                     ("post_vs_post", f"post-{r['cutoff']} vs post-{r['cutoff']}"),
                     ("pre_vs_pre", f"pre-{r['cutoff']} vs pre-{r['cutoff']}"),
                     ("pre_pos_vs_post_neg",
                      f"pre-{r['cutoff']} positives vs post-{r['cutoff']} negatives")):
        v = e[key]
        print(f"   {lab:<48}{'—' if v is None else f'{v:.4f}'}")
    print()

    p = r["prospective"]
    print("3. THE SAME COMPOUNDS, THE OTHER PROTOCOL")
    print(f"   post-cutoff compounds, cross-sectional references : "
          f"{e['post_vs_post']:.4f}")
    print(f"   post-cutoff compounds, pre-cutoff references only : "
          f"{p['primary']:.4f}  (k10: {p['symmetric_k10']:.4f})")
    print("   identical compounds and labels. The reference set is the "
          "only difference.")
    print()

    s = r["era_mate_share"]
    print("4. THE MECHANISM")
    if s["median"] is None:
        print("   no post-cutoff positive could be fingerprinted.")
    else:
        print(f"   For a post-{r['cutoff']} positive, of the {r['k_cell']} "
              "nearest references")
        print("   the cross-sectional control REMOVES, the share that are "
              "its own era-mates:")
        print(f"       median {s['median']:.0%}   mean {s['mean']:.0%}   "
              f"n = {s['n_compounds']}")
        print()
        print("   Those era-mates are not in the prospective reference set "
              "at all, so")
        print("   the prospective control cannot remove them. The "
              "cross-sectional control")
        print("   is strictly stronger for exactly these compounds.")
    print()
    print(RULE)
    print("READING IT")
    print(RULE)
    print()
    print("   Both registered results stand. They answer different questions:")
    print()
    print("     cross-sectional, strength-matched")
    print("       'setting aside a compound's own analogues, does chemistry")
    print("        separate actives from inactives?'")
    print("     prospective")
    print("       'does pre-cutoff chemistry rank post-cutoff compounds?'")
    print()
    print("   A prospective pass therefore does NOT license a shortlist: the")
    print("   shortlist ranks never-screened compounds against ALL known")
    print("   actives, which is the cross-sectional setting.")
    print()

    if args.out:
        args.out.write_text(json.dumps(r, indent=2) + "\n")
        print(f"   wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
