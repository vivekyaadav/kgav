#!/usr/bin/env python3
"""Does the structure scorer predict forward? Pre-2021 references, post-2021 labels.

    python scripts/temporal_structure_scorer.py
    make temporal-structure-scorer

REGISTERED IN ADVANCE by doc2/registrations/prereg_temporal_structure_scorer.md,
signed 2026-11-19 before any AUC against a post-2021 label existed.
Every threshold, the family size, the alpha, the refusal conditions and
the 20% comparator rule are read from that document's decisions and are
NOT re-chosen here.

WHAT IS DIFFERENT FROM scripts/structure_scorer.py, and it is only this:
the reference set is the pre-2021 actives and the evaluated population
is the post-2021 labels. The statistic, the fingerprint, the skeleton
exclusion, the imputation and the interval all come from that module by
import, because a second copy of a statistic drifts.

WHY IT IS WORTH RUNNING. Everything else in this repository is
retrospective: known actives separated from known inactives measured at
the same time, largely from the same papers. The graph channels cannot
answer the forward question at all -- 95.2% of post-2021 test compounds
have no pre-2021 evidence, so they can rank 4.8% of them -- while a
fingerprint needs only the molecule. That asymmetry is the hazard H1
registered, not the result: a structure AUC that beats a degree baseline
made entirely of ties has beaten nothing, so this reports how much of
each cell the comparator can even see.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from statistics import NormalDist

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.baselines import Graph, degree_ranking            # noqa: E402
from kgav.labels import evaluate_against_negatives          # noqa: E402
from kgav.pipeline import release_dir, train_dir           # noqa: E402
from kgav.provenance import (                               # noqa: E402
    release_fingerprint,
    write_refusal,
    write_results,
)
from kgav.schema import load_schema, print_capped           # noqa: E402
from kgav.similarity import fingerprints                    # noqa: E402

from structure_scorer import (                              # noqa: E402
    K_PRIMARY,
    MIN_CLASS,
    corrected,
    score_cell,
    small_molecule_props,
)

#: The policies the split must have been built with. Any other split is
#: a different study, and the run refuses rather than measuring one.
#: The DIRECTORY is not named here -- pipeline.train_dir(ROOT, CUTOFF)
#: owns it, so the cutoff below is the single place it is chosen.
CUTOFF = 2021
UNDATED = "include"

#: doc2/registrations/prereg_temporal_structure_scorer.md: ten evaluable cells, so
#: alpha = 0.05/10. Written down so a run that finds a different count
#: REFUSES rather than dividing by whatever it found.
FAMILY = 10
ALPHA = 0.05 / FAMILY
Z_CORRECTED = round(NormalDist().inv_cdf(1 - ALPHA / 2), 4)   # 2.8070

#: H1's rule, fixed before the run. Under this fraction of test
#: compounds carrying any training-graph degree, the degree baseline is
#: mostly ties and the cell is reported as having NO USABLE COMPARATOR
#: rather than as one the scorer won.
MIN_COMPARATOR_REACH = 0.20

#: A refusal condition, not a tuning knob.
MIN_FINGERPRINT_COVERAGE = 0.95


def cells(labels: dict) -> dict[str, tuple[list, list, list]]:
    """{taxon: (references, test positives, test negatives)}."""
    out = {}
    for t in sorted(set(labels["test"]) | set(labels["train"])):
        out[t] = (sorted(labels["train"].get(t, [])),
                  sorted(labels["test"].get(t, [])),
                  sorted(labels["test_negatives"].get(t, [])))
    return out


def evaluable(c: dict) -> dict:
    return {t: v for t, v in c.items()
            if len(v[0]) >= MIN_CLASS and len(v[1]) >= MIN_CLASS
            and len(v[2]) >= MIN_CLASS}


def comparator(g: Graph, pos: list[str], neg: list[str]) -> dict:
    """The degree baseline on the TRAINING graph, with its own reach.

    `degree_ranking` refuses a graph that skipped nothing, which is the
    guard added on 2026-11-17 after the cross-sectional comparator was
    found reading the held-out labels. It is relied on here rather than
    re-checked.
    """
    ranks = degree_ranking(g)
    reached = sum(1 for d in list(pos) + list(neg) if ranks.get(d))
    total = len(pos) + len(neg)
    frac = reached / total if total else 0.0
    m = evaluate_against_negatives(ranks, set(pos), set(neg))
    m = corrected(m, z=Z_CORRECTED)
    m["reach"] = frac
    m["reached"] = reached
    # TWO DIFFERENT WAYS TO HAVE NO COMPARATOR, and they are reported
    # apart. `evaluable` false means the evaluator could not produce an
    # AUC at all; reach below the registered floor means it produced one
    # from a column that is almost entirely ties. H1 anticipated the
    # second; the first turned up on the first run, because a training
    # graph that reaches nobody gives degree_ranking nothing to rank.
    m["usable"] = bool(m.get("auc") is not None
                       and frac >= MIN_COMPARATOR_REACH)
    return m


def _auc(m: dict) -> str:
    """An AUC, or a word saying there is none.

    A cell with no comparator must not format as 0.0000: that reads as a
    measured failure of the baseline rather than the absence of one, and
    it is the same mistake doc2/scorer_aucs.txt's `_cell()` was written
    to avoid.
    """
    v = m.get("auc")
    return "--" if v is None else f"{v:.4f}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    # train_dir and release_dir, not literals: pipeline.py owns the
    # release path and tests/test_release_path_has_one_owner.py enforces
    # it, which is how this script's first draft was caught spelling
    # "data/releases/v0.1" itself.
    ap.add_argument("--split", type=Path, default=train_dir(ROOT, CUTOFF))
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    s = load_schema()
    doc = "doc2/registrations/prereg_temporal_structure_scorer.md"
    refusals: list[str] = []

    derived = json.loads((args.split / "DERIVED_FROM.json").read_text())
    split_sha = ((derived.get("provenance") or {})
                 .get("release") or {}).get("edges_sha")
    # release_fingerprint, not MANIFEST.json: the manifest carries the
    # census and the build commit, and the shas live in provenance.
    release_sha = (release_fingerprint(args.release).get("edges_sha")
                   if (args.release / "nodes.jsonl").exists() else None)
    print(f"split   {args.split.name}   cutoff {derived.get('cutoff')}   "
          f"undated {derived.get('undated')}")
    print(f"        built from edges_sha {split_sha}")
    print(f"release {args.release.name}  edges_sha {release_sha}")
    if split_sha and release_sha and split_sha != release_sha:
        refusals.append(
            f"the split was built from edges_sha {split_sha} and the "
            f"release on disk is {release_sha}: the split is stale")
    if derived.get("cutoff") != CUTOFF or derived.get("undated") != UNDATED:
        refusals.append(
            f"this study is registered for cutoff {CUTOFF}/"
            f"undated={UNDATED}, and the split is "
            f"{derived.get('cutoff')}/{derived.get('undated')}")

    labels = json.loads((args.split / "TEST_LABELS.json").read_text())
    allc = cells(labels)
    ev = evaluable(allc)
    print(f"\n{len(ev)} evaluable cell(s) of {len(allc)} "
          f"(>= {MIN_CLASS} references, test positives and test negatives)")
    if len(ev) != FAMILY:
        refusals.append(
            f"the registration fixed the family at {FAMILY} and this split "
            f"yields {len(ev)}; the alpha would not be the registered one")

    props = small_molecule_props(args.release)
    everyone = sorted({d for v in ev.values() for part in v for d in part})
    fps, fp_stats = fingerprints(
        {i: (props.get(i, {}).get("smiles") or "") for i in everyone})
    cov = fp_stats["fingerprinted"] / len(everyone) if everyone else 0.0
    print(f"fingerprinted {fp_stats['fingerprinted']:,} of {len(everyone):,} "
          f"({cov:.1%}); unparseable {fp_stats['unparseable']:,}, "
          f"no SMILES {fp_stats['no_smiles']:,}")
    if cov < MIN_FINGERPRINT_COVERAGE:
        refusals.append(f"fingerprint coverage {cov:.1%} is below the "
                        f"registered {MIN_FINGERPRINT_COVERAGE:.0%}")

    payload = {
        "measured_by": "scripts/temporal_structure_scorer.py",
        "registered_by": doc,
        "split": args.split.name,
        "cutoff": derived.get("cutoff"),
        "undated_policy": derived.get("undated"),
        "overlap_policy": derived.get("overlap_policy"),
        "statistic": "max Tanimoto to PRE-CUTOFF known actives "
                     "(ECFP4 2048/r2), skeleton-excluded",
        "family_size": FAMILY,
        "alpha": ALPHA,
        "z_corrected": Z_CORRECTED,
        "min_comparator_reach": MIN_COMPARATOR_REACH,
        "fingerprint_stats": dict(fp_stats),
        "n_evaluable": len(ev),
        "n_cells_total": len(allc),
    }

    if refusals:
        print("\nREFUSED, for the reasons the registration fixed in advance:")
        for r in refusals:
            print(f"  - {r}")
        if args.out:
            write_refusal(args.out, payload, refusals, args.release, s.version)
            print(f"\nwrote a refusal record to {args.out}")
        return 1

    g = Graph.load(args.split, skip_predicates=s.held_out_predicates(),
                   symmetry=s.symmetry())
    vy = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    names = {f"NCBITaxon:{e['taxon']}": e["name"]
             for lst in ("viruses", "surveyed") for e in (vy.get(lst) or [])}
    nopubs: dict[str, set[str]] = {}

    print(f"\n{'=' * 90}")
    print("PRIMARY: no exclusion beyond the split. The reference set is "
          "strictly EARLIER\nthan the candidate, which is the prospective "
          "situation; removing a near neighbour\nwould remove the signal "
          "rather than a leak. 'k10' is the registered secondary.")
    print(f"{'=' * 90}")
    print(f"  {'virus':14} {'refs':>5} {'pos':>6} {'neg':>6} {'AUC':>7} "
          f"{'95% CI (corrected)':>22} {'clears':>7} {'k10':>7} "
          f"{'degree':>7} {'reach':>6}")
    print("  " + "-" * 100)

    out_cells = {}
    for t, (refs, pos, neg) in sorted(ev.items()):
        # THROUGH THE REGISTER, NOT THE GRAPH. The first run printed
        # bare taxon ids because an OrganismTaxon node carries no
        # `name` property -- and a taxon id read as a virus name is the
        # error that left a landscape cell open for two days
        # (doc2/registrations/prereg_pathway_overlap.md, amendment 2026-11-18).
        name = names.get(t, t)
        cand = list(pos) + list(neg)
        prim, c1 = score_cell(cand, refs, fps, props, nopubs,
                              publication_disjoint=False, knn=0)
        m = corrected(evaluate_against_negatives(prim, set(pos), set(neg)),
                      z=Z_CORRECTED)
        sec, c2 = score_cell(cand, refs, fps, props, nopubs,
                             publication_disjoint=False, knn=K_PRIMARY)
        m2 = corrected(evaluate_against_negatives(sec, set(pos), set(neg)),
                       z=Z_CORRECTED)
        dg = comparator(g, pos, neg)
        out_cells[t] = {
            "name": name, "n_references": len(refs),
            "n_test_positives": len(pos), "n_test_negatives": len(neg),
            "primary": m, "primary_counts": c1,
            "symmetric_k10": m2, "symmetric_k10_counts": c2,
            "degree": dg,
        }
        print(f"  {name:14} {len(refs):5,} {len(pos):6,} {len(neg):6,} "
              f"{m['auc']:7.4f} "
              f"[{m['auc_ci_lo_corrected']:.4f}, "
              f"{m['auc_ci_hi_corrected']:.4f}] "
              f"{'YES' if m['clears_corrected_alpha'] else '-':>7} "
              f"{m2['auc']:7.4f} {_auc(dg):>7} "
              f"{dg['reach']:5.1%}"
              f"{'' if dg['usable'] else '  <- NO USABLE COMPARATOR'}")

    payload["cells"] = out_cells
    clears = [t for t, c in out_cells.items()
              if c["primary"]["clears_corrected_alpha"]
              and c["primary"]["direction"] == "above"]
    usable = [t for t, c in out_cells.items() if c["degree"]["usable"]]
    beat = [t for t in clears
            if out_cells[t]["degree"]["usable"]
            and out_cells[t]["primary"]["auc"]
            > out_cells[t]["degree"]["auc"]]
    payload |= {"clears": sorted(clears), "comparator_usable": sorted(usable),
                "clears_and_beats_usable_comparator": sorted(beat)}

    print(f"\n{len(clears)} of {len(ev)} cell(s) clear alpha={ALPHA:.6f} "
          f"above chance.")
    print_capped([out_cells[t]["name"] for t in sorted(clears)], 20)
    print(f"{len(usable)} of {len(ev)} cell(s) have a usable degree "
          f"comparator (>= {MIN_COMPARATOR_REACH:.0%} of test compounds "
          f"carry training degree).")
    print(f"{len(beat)} cell(s) BOTH clear and beat a usable comparator.")

    print("\nTHE REGISTERED DECISION RULE:")
    if beat:
        print("  P1 holds and at least one clearing cell beats a comparator "
              "it can\n  be compared against. The scorer has a prospective "
              "result.")
    elif clears:
        print("  MEASURABLE, NOT COMPARABLE. Cells clear, but none of them "
              "has a degree\n  baseline that sees enough of its test set to "
              "be a baseline. H1 said to\n  expect this and said reporting "
              "it as a success would be the failure\n  this registration "
              "exists to prevent.")
    else:
        print("  P1 FAILS. The scorer does not predict forward on this "
              "data, and\n  doc2/08-results.md's claim section must say so "
              "beside its 8 cells.")

    if args.out:
        write_results(args.out, payload, args.release, s.version)
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
