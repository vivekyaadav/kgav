#!/usr/bin/env python3
"""Compare metapaths across all seven viruses.

Raw Hits@k is not comparable across viruses here: base rates run from 14% for
SARS-CoV-2 to 0.05% for NL63. Everything below is reported as lift over each
pool's own random expectation, with the expectation printed so a lift built on
fewer than one expected hit can be recognised as noise rather than signal.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.baselines import (
    Graph,
    combine,
    degree_ranking,
    dwpc_scores,
    hits_at_k,
    mrr,
    rank,
)
from kgav.schema import load_schema

MIN_TRUSTWORTHY_EXPECTATION = 1.0


def positives(release: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if (e["predicate"] == "HAS_ANTIVIRAL_ACTIVITY_AGAINST"
                and (e.get("qualifiers") or {}).get("relation") == "="):
            out.setdefault(e["object"], set()).add(e["subject"])
    return out


def row(name: str, scores: dict[str, float], pos: set[str], k: int = 100) -> dict:
    r = rank(scores)
    p = len(pos & set(scores))
    exp = k * p / len(r) if r else 0.0
    h = hits_at_k(r, pos, k)
    return {"pool": len(r), "pos": p, "hits": h, "expected": exp,
            "lift": (h / exp) if exp else 0.0, "mrr": mrr(r, pos),
            "trustworthy": exp >= MIN_TRUSTWORTHY_EXPECTATION}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--out", type=Path, default=root / "data/results")
    ap.add_argument("-k", type=int, default=100)
    args = ap.parse_args()

    s = load_schema()
    g = Graph.load(args.release, skip_predicates=s.held_out_predicates())
    pos = positives(args.release)
    labels = {n: d["properties"]["label"] for n, d in g.nodes.items()
              if d["class"] == "OrganismTaxon"}

    print(f"schema v{s.version} | {len(g.drugs()):,} drugs | "
          f"metapaths {', '.join(sorted(s.metapaths))}")
    print(f"lift is Hits@{args.k} over random expectation for that pool; "
          f"(!) marks expectation < {MIN_TRUSTWORTHY_EXPECTATION} where lift is noise\n")

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
                f"{'exp':>7} {'lift':>7} {'mrr':>7}")
        per: dict[str, dict] = {}
        for name in sorted(dwpc):
            sc = {d: by[virus] for d, by in dwpc[name].items() if virus in by}
            if not sc:
                print(f"    {name:<10} {'no paths':>7}")
                continue
            per[name] = row(name, sc, p, args.k)
        per["COMBINED"] = row("COMBINED", combine(dwpc, virus), p, args.k)
        per["degree"] = row("degree", deg, p, args.k)
        for name, m in per.items():
            flag = "" if m["trustworthy"] else " (!)"
            print(f"    {name:<10} {m['pool']:>7,} {m['pos']:>5,} {m['hits']:>5} "
                  f"{m['expected']:>7.1f} {m['lift']:>7.2f} {m['mrr']:>7.4f}{flag}")
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

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "metapath_comparison.json").write_text(json.dumps(results, indent=2))
    print(f"\nwrote {args.out / 'metapath_comparison.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
