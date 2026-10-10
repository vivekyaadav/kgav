#!/usr/bin/env python3
"""Evaluate against MEASURED negatives, cross-sectionally and temporally.

The question every previous run could not answer: does the graph rank actives
above compounds that were assayed and found INACTIVE, rather than merely above
compounds nobody tested?
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.baselines import (
    Graph,
    combine,
    degree_ranking,
    dwpc_scores,
    label_publication_index,
    metapath_reach,
)
from kgav.labels import build_labels, evaluate_against_negatives
from kgav.provenance import write_results
from kgav.schema import load_schema

MIN_CLASS_SIZE = 10


def run(g: Graph, schema, labels, tag: str,
        k_min: int = MIN_CLASS_SIZE) -> tuple[dict, dict]:
    """(per-virus results, per-metapath reach).

    Reach is returned alongside rather than mixed into the virus-keyed dict,
    because a consumer iterating the results should see viruses and nothing
    else.
    """
    dwpc = dwpc_scores(g, schema.metapaths, schema.retrieval_limits)
    deg = degree_ranking(g)
    out: dict[str, dict] = {}

    # Reported before the per-virus tables, because a metapath absent from
    # every table below is otherwise indistinguishable from one that was never
    # declared. `if not sc: continue` hid M8 here with no stdout line at all.
    reach = metapath_reach(g, schema.metapaths, dwpc)
    dead = [n for n, r in reach.items() if not r["drugs_reached"]]
    if dead:
        print(f"\n[{tag}] metapaths reaching NO drug: {', '.join(dead)}")
        for name in dead:
            miss = reach[name]["missing_predicates"]
            if miss:
                print(f"    {name}: no edges for {', '.join(miss)}")


    viruses = sorted(labels.positives, key=lambda v: -len(labels.positives[v]))
    for virus in viruses:
        pos, neg = labels.positives.get(virus, set()), labels.negatives.get(virus, set())
        label = (g.nodes.get(virus, {}).get("properties", {}) or {}).get("label", virus)
        print(f"\n=== {label} [{tag}]")
        print(f"    {len(pos):,} measured active | {len(neg):,} measured INACTIVE")
        if len(pos) < k_min or len(neg) < k_min:
            print(f"    skipped: need at least {k_min} of each class")
            continue

        scorers = {name: {d: by[virus] for d, by in dwpc[name].items() if virus in by}
                   for name in sorted(dwpc)}
        scorers["COMBINED"] = combine(dwpc, virus)
        scorers["degree"] = deg

        print(f"    {'scorer':<10} {'pos':>5} {'neg':>6} {'AUC':>7} "
              f"{'cov_pos':>8} {'cov_neg':>8}")
        per: dict[str, dict] = {}
        for name, sc in scorers.items():
            if not sc:
                continue
            m = evaluate_against_negatives(sc, pos, neg)
            per[name] = m
            if not m["evaluable"]:
                print(f"    {name:<10} {m['n_pos']:>5} {m['n_neg']:>6} "
                      f"{'-':>7}   one class empty in this pool")
                continue
            print(f"    {name:<10} {m['n_pos']:>5} {m['n_neg']:>6} {m['auc']:>7.3f} "
                  f"{m['coverage_pos']:>7.1%} {m['coverage_neg']:>7.1%}")
        out[label] = per
    return out, reach


def _filter_labels(labels, si: dict, args) -> None:
    """Apply the selectivity filter to the sides `--selectivity-applies` names.

    One function for both sides, so the two cannot drift: the asymmetry was
    not a decision recorded anywhere, it was the negatives simply never being
    passed through the same loop.
    """
    sides = [("actives", labels.positives)]
    if args.selectivity_applies == "both":
        sides.append(("inactives", labels.negatives))
    for what, group in sides:
        for virus, members in list(group.items()):
            if args.selectivity == "verified-only":
                keep = {d for d in members if (d, virus) in si}
            else:
                keep = {d for d in members
                        if si.get((d, virus), -1) >= args.si_threshold}
            dropped = len(members) - len(keep)
            group[virus] = keep
            if dropped:
                print(f"  {virus}: {len(members):,} -> {len(keep):,} {what} "
                      f"({dropped:,} excluded)")


def _report_drops(g: Graph, tag: str) -> None:
    """Name the predicates the publication-disjoint filter withheld.

    Printed per graph rather than summed, because an AUC that moves while
    nothing was dropped would mean the filter is not doing what it says.
    """
    n = g.stats["dropped_publication_coincident"]
    print(f"  {tag}: {n:,} of {g.stats['edges_seen']:,} edges withheld as "
          f"publication-coincident with their own label")
    for k in sorted(k for k in g.stats if k.startswith("dropped_")
                    and k != "dropped_publication_coincident"):
        print(f"    {k.replace('dropped_', ''):<34} {g.stats[k]:>7,}")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--train", type=Path, default=root / "data/releases/v0.1-train2021")
    ap.add_argument("--results", type=Path, default=root / "data/results")
    ap.add_argument("--cutoff", type=int, default=2021)
    ap.add_argument("--inactive-above-nm", type=float, default=10_000.0)
    ap.add_argument("--selectivity", default="all",
                    choices=["all", "selective-only", "verified-only"],
                    help="all: every measured active (the original protocol). "
                         "selective-only: actives with a verified SI >= 10. "
                         "verified-only: actives with ANY verified SI, "
                         "selective or not -- isolates the effect of "
                         "restricting to compounds with matched data from the "
                         "effect of the selectivity filter itself.")
    ap.add_argument("--si-threshold", type=float, default=10.0)
    ap.add_argument("--label-source", default=None,
                    help="RESTRICT THE LABELS TO ONE SCREEN, e.g. "
                         "infores:ncats-opendata (comma-separate for several). "
                         "Pooling ChEMBL and panel labels reintroduces the "
                         "confound the panel removes: ChEMBL actives are "
                         "antiviral research compounds carrying viral-target "
                         "annotations and panel inactives are library "
                         "compounds carrying host-target ones, so a scorer can "
                         "separate the classes by recognising which collection "
                         "a compound came from. One screen means both classes "
                         "come off the same plates under one protocol -- "
                         "matched by construction, not by a filter.")
    ap.add_argument("--publication-disjoint", action="store_true",
                    help="WITHHOLD EVIDENCE THAT SHARES A PAPER WITH ITS OWN "
                         "LABEL. 94.9% of compounds carrying both an INHIBITS "
                         "edge and an antiviral-activity label cite the same "
                         "PMID on both, because one paper reports the target "
                         "IC50 and the cell-based EC50 together. Without this "
                         "flag M1 traverses a measurement taken alongside the "
                         "answer, and neither the compound split nor the "
                         "temporal split separates them. Off by default so "
                         "published numbers reproduce.")
    ap.add_argument("--selectivity-applies", default="positives",
                    choices=["positives", "both"],
                    help="WHICH SIDE THE SELECTIVITY FILTER TOUCHES. "
                         "'positives' (default) reproduces every published "
                         "number: actives are restricted to compounds with a "
                         "verified index while negatives are left whole. "
                         "'both' is the symmetric control -- requiring the "
                         "same paired cytotoxicity evidence of a compound "
                         "before it may be a negative.")
    args = ap.parse_args()

    schema = load_schema()
    held = schema.held_out_predicates()

    label_sources = ({s.strip() for s in args.label_source.split(",") if s.strip()}
                     if args.label_source else None)
    if label_sources:
        print(f"labels restricted to {sorted(label_sources)}")
    all_labels = build_labels(args.release, inactive_above_nm=args.inactive_above_nm,
                              sources=label_sources)
    if label_sources and not any(all_labels.positives.values()):
        print(f"\nNO POSITIVES from {sorted(label_sources)}. Check the source "
              f"string against\n  the release: a typo here yields an empty "
              f"evaluation, not an error.")
        return 1

    if args.selectivity != "all":
        from kgav.selectivity import selectivity_index_by_compound
        si = selectivity_index_by_compound(args.release)
        print(f"selectivity filter '{args.selectivity}' applied to "
              f"{args.selectivity_applies}: "
              f"{len(si):,} (compound, virus) pairs have a verified index")
        _filter_labels(all_labels, si, args)
        if args.selectivity_applies == "positives":
            # NOT A LIKE-FOR-LIKE COMPARISON, and the published AUC is computed
            # this way. Requiring a verified index of a POSITIVE means requiring
            # paired same-document CC50/EC50, which selects for
            # well-characterised compounds -- the ones most likely to carry the
            # nsp5/nsp12 target annotations M1 reads. Measured on v0.1: the
            # positive set halves 1,726 -> 808 and M1's AUC rises 0.729 ->
            # 0.823, while cov_pos jumps 51.5% -> 70.3% and cov_neg does not
            # move from 5.3%. That coverage asymmetry is the signature of
            # enrichment for REACHABILITY, not only for genuine activity.
            #
            # --selectivity-applies both is the control. If the AUC holds near
            # 0.82 the effect is real; if it falls toward 0.73 it was annotation
            # bias. Default stays 'positives' so published numbers reproduce.
            print("  NOTE: negatives are NOT filtered. Actives must have paired "
                  "cytotoxicity\n  data and inactives need not, so the two "
                  "classes are drawn from differently\n  characterised "
                  "populations. Run --selectivity-applies both for the control.")
    print(f"labels from the full release (inactive if censored above "
          f"{args.inactive_above_nm:,.0f} nM):")
    for key in ("active", "inactive", "undecidable", "ambiguous_dropped",
                "other_source"):
        if all_labels.stats[key]:
            print(f"  {key:<20} {all_labels.stats[key]:>7,}")

    print("\n" + "=" * 66)
    print("CROSS-SECTIONAL: full graph, all measured compounds")
    print("=" * 66)
    pub_idx: dict[str, set[str]] | None = None
    if args.publication_disjoint:
        pub_idx = label_publication_index(args.release, held)
        print(f"publication-disjoint: {len(pub_idx):,} compounds carry label "
              f"publications; evidence citing the same paper is withheld")

    g_full = Graph.load(args.release, skip_predicates=held,
                        symmetry=schema.symmetry(),
                        drop_pubs_shared_with=pub_idx)
    if pub_idx:
        _report_drops(g_full, "cross-sectional")
    cross, cross_reach = run(g_full, schema, all_labels, "cross-sectional")

    # Declared AFTER the cross-sectional run on purpose -- cross_reach is
    # already set above and must not be re-initialised here, which is what
    # batch 2a did: it reset cross_reach to {} one line after run() returned
    # it, so the cross-sectional metapath reach never reached the results file.
    temporal: dict = {}
    temporal_reach: dict = {}
    if args.train.exists():
        print("\n" + "=" * 66)
        print(f"TEMPORAL: pre-{args.cutoff} graph, post-{args.cutoff} measurements")
        print("=" * 66)
        test_labels = build_labels(args.release, year_from=args.cutoff + 1,
                                   inactive_above_nm=args.inactive_above_nm,
                                   sources=label_sources)
        # The selectivity filter must apply here too. Filtering only the
        # cross-sectional labels would report a temporal number computed on a
        # different positive set from the one it is being compared against.
        if args.selectivity != "all":
            from kgav.selectivity import selectivity_index_by_compound
            # Through the SAME helper as the cross-sectional side. This was a
            # second hand-written copy of the loop, which is how the two sides
            # came to differ without anyone deciding they should.
            _filter_labels(test_labels, selectivity_index_by_compound(args.release),
                           args)
        print(f"post-{args.cutoff} labels: {test_labels.stats['active']:,} active, "
              f"{test_labels.stats['inactive']:,} inactive")
        # Same index, built from the FULL release: a label is a label wherever
        # it sits, and an edge in the pre-cutoff graph citing the paper that
        # reported the post-cutoff label is leakage across the split, not
        # evidence the split legitimately preserved.
        g_train = Graph.load(args.train, skip_predicates=held,
                             symmetry=schema.symmetry(),
                             drop_pubs_shared_with=pub_idx)
        if pub_idx:
            _report_drops(g_train, f"pre-{args.cutoff} train")
        temporal, temporal_reach = run(g_train, schema, test_labels,
                                       f"post-{args.cutoff}")
    else:
        print(f"\n{args.train} not found -- run temporal_split.py first")

    args.results.mkdir(parents=True, exist_ok=True)
    write_results(args.results / "hard_negatives.json",
                  {"inactive_above_nm": args.inactive_above_nm,
                   "selectivity_filter": args.selectivity,
                   "selectivity_applies": args.selectivity_applies,
                   # THE TWO CONTROLS THAT DECIDE WHAT THESE NUMBERS MEAN, and
                   # neither was recorded until now. M1 reads 0.660 without
                   # publication_disjoint and 0.499 with it, on one release
                   # and one commit -- so the release fingerprint and the
                   # commit that provenance.py stamps do NOT distinguish the
                   # two runs. A file holding 0.660 was indistinguishable from
                   # one holding 0.499 by anything except the number, which is
                   # the exact failure that module exists to prevent.
                   "publication_disjoint": bool(args.publication_disjoint),
                   "label_source": sorted(label_sources) if label_sources else None,
                   "label_stats": dict(all_labels.stats),
                   "cross_sectional": cross, "temporal": temporal,
                   "metapath_reach": {"cross_sectional": cross_reach,
                                      "temporal": temporal_reach}},
                  args.release, schema.version)
    print(f"\nwrote {args.results / 'hard_negatives.json'}")
    print("\nAUC 0.5 is chance. Unlike Hits@k this is prevalence-independent, so "
          "it is comparable\nacross viruses whose base rates differ by two orders "
          "of magnitude.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
