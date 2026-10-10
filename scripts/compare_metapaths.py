#!/usr/bin/env python3
"""Compare metapaths across all seven viruses.

Raw Hits@k is not comparable across viruses here: base rates run from 14% for
SARS-CoV-2 to 0.05% for NL63. Everything below is reported as lift over each
pool's own random expectation, with the expectation printed so a lift built on
fewer than one expected hit can be recognised as noise rather than signal.
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
    confound_warning,
    degree_relative,
    dwpc_scores,
    hits_at_k,
    lift_row,
    metapath_reach,
    mrr,
    rank,
)
from kgav.labels import build_labels
from kgav.provenance import write_results
from kgav.schema import load_schema

def positives(release: Path) -> dict[str, set[str]]:
    """Measured-ACTIVE compounds per virus, via labels.classify.

    Selecting on relation == "=" alone put 578 compounds measured above 10 uM
    into the positive set: they were assayed and found inactive, and a lift
    computed against them measures nothing.
    """
    return build_labels(release).positives


def row(name: str, scores: dict[str, float], pos: set[str], k: int = 100) -> dict:
    """Thin wrapper; the metric lives in baselines so it can be tested."""
    return lift_row(scores, pos, k)


def _vs(m: dict) -> str:
    """Lift relative to the degree baseline, which is the comparison that
    carries information; p alone is against a random ranking."""
    v = m.get("lift_vs_degree")
    if v is None:
        return "      -"
    return f"{v:>7.2f}"


def _p(m: dict) -> str:
    """p shown only where it answers something: None means the question was
    empty (no positives, or k covering the whole pool)."""
    v = m.get("p_value")
    return "      -" if v is None else (f"{v:>7.4f}" if v >= 0.0001
                                        else "< .0001")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--out", type=Path, default=root / "data/results")
    ap.add_argument("-k", type=int, default=100)
    args = ap.parse_args()

    s = load_schema()
    g = Graph.load(args.release, skip_predicates=s.held_out_predicates(),
                   symmetry=s.symmetry())
    pos = positives(args.release)
    labels = {n: d["properties"]["label"] for n, d in g.nodes.items()
              if d["class"] == "OrganismTaxon"}

    print(f"schema v{s.version} | {len(g.drugs()):,} drugs | "
          f"metapaths {', '.join(sorted(s.metapaths))}")
    print(f"lift is Hits@{args.k} over random expectation for that pool; "
          f"(!) marks a lift that is noise: expectation < 1, or a pool no "
          f"larger than k, where lift is exactly pool/k\n")

    dwpc = dwpc_scores(g, s.metapaths, s.retrieval_limits)
    deg = degree_ranking(g)
    results: dict[str, dict] = {}

    for virus, label in sorted(labels.items(), key=lambda kv: -len(pos.get(kv[0], ()))):
        p = pos.get(virus, set())
        print(f"=== {label}  ({len(p):,} positives)")
        if not p:
            print("    unevaluable: no measured activity\n")
            continue
        print(f"    {'scorer':<10} {'pool':>7} {'pos':>5} {'hits':>5} "
                f"{'exp':>7} {'lift':>7} {'p':>7} {'vs deg':>7} {'mrr':>7}")
        per: dict[str, dict] = {}
        for name in sorted(dwpc):
            sc = {d: by[virus] for d, by in dwpc[name].items() if virus in by}
            if not sc:
                print(f"    {name:<10} {'no paths':>7}")
                continue
            per[name] = row(name, sc, p, args.k)
        per["COMBINED"] = row("COMBINED", combine(dwpc, virus), p, args.k)
        per["degree"] = row("degree", deg, p, args.k)
        degree_relative(per)
        for name, m in per.items():
            flag = "" if m["trustworthy"] else " (!)"
            print(f"    {name:<10} {m['pool']:>7,} {m['pos']:>5,} {m['hits']:>5} "
                  f"{m['expected']:>7.1f} {m['lift']:>7.2f} {_p(m)} "
                  f"{_vs(m)} {m['mrr']:>7.4f}{flag}")
        if (warn := confound_warning(per)):
            print(f"    ! {warn}")
        results[label] = per
        print()

    # M6 vs M7 head to head: the reason this script exists.
    print("=== M6 vs M7 (cold-start routes)")
    print(f"    {'virus':<14} {'M6 pool':>8} {'M6 lift':>8} {'M7 pool':>8} "
          f"{'M7 lift':>8}  verdict")
    for label, per in results.items():
        m6, m7 = per.get("M6"), per.get("M7")
        if not m6 and not m7:
            continue
        def cell(m, key, fmt):
            return format(m[key], fmt) if m else "-"
        if m6 and m7:
            verdict = ("M7 better" if m7["lift"] > m6["lift"] * 1.2
                       else "M6 better" if m6["lift"] > m7["lift"] * 1.2
                       else "comparable")
            if not (m6["trustworthy"] or m7["trustworthy"]):
                verdict += " (both noise)"
        else:
            verdict = "M7 only" if m7 else "M6 only"
        print(f"    {label:<14} {cell(m6,'pool','>8,'):>8} {cell(m6,'lift','>8.2f'):>8} "
              f"{cell(m7,'pool','>8,'):>8} {cell(m7,'lift','>8.2f'):>8}  {verdict}")

    # A metapath that reaches nothing was previously printed to stdout and
    # then omitted from this file, so the record could not distinguish "scored
    # badly" from "never ran". M8 has been in that state since it was declared.
    reach = metapath_reach(g, s.metapaths, dwpc)
    dead = {k: v for k, v in reach.items() if not v["drugs_reached"]}
    if dead:
        print("\n=== metapaths that reached NO drug at all")
        for name, r in dead.items():
            why = (f"no edges for {', '.join(r['missing_predicates'])}"
                   if r["missing_predicates"] else
                   "every hop's predicate is present, so the break is in the "
                   "edges themselves (orientation, endpoints or constraints)")
            print(f"    {name:<10} {why}")

    args.out.mkdir(parents=True, exist_ok=True)
    write_results(args.out / "metapath_comparison.json",
                  {"per_virus": results, "metapath_reach": reach},
                  args.release, s.version)
    print(f"\nwrote {args.out / 'metapath_comparison.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
