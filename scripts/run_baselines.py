#!/usr/bin/env python3
"""Run the three non-learned baselines and the biology controls.

The metric is not the point. The point is whether known antivirals surface and
whether known artifacts do not.
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
    network_proximity,
    rank,
)
from kgav.schema import load_schema

# Compounds whose behaviour tells us whether the graph works.
#   nirmatrelvir  nsp5 inhibitor  -> must rank high via M1
#   remdesivir    nsp12 inhibitor -> must rank high via M1
#   hydroxychloroquine / chloroquine: abundant SARS-CoV-2 activity generated
#   largely in Vero E6, where entry is cathepsin-dependent. Both failed in
#   TMPRSS2-expressing airway cells. If these rank high the graph is
#   reproducing a cell-line artifact, and with no CC50 data we cannot filter
#   them on selectivity.
CONTROLS = {
    "nirmatrelvir": "positive",
    "remdesivir": "positive",
    "GS-441524": "positive",
    "hydroxychloroquine": "artifact",
    "chloroquine": "artifact",
}


def name_index(g: Graph, db: Path | None) -> dict[str, str]:
    """chembl_id -> node id, plus a lowercase pref_name lookup if the ChEMBL
    database is available."""
    by_chembl = {n["properties"].get("chembl_id"): nid
                 for nid, n in g.nodes.items()
                 if n["class"] == "SmallMolecule" and n["properties"].get("chembl_id")}
    out: dict[str, str] = {}
    if db and db.exists():
        import sqlite3
        con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        for cid, pref in con.execute(
                "SELECT chembl_id, lower(pref_name) FROM molecule_dictionary "
                "WHERE pref_name IS NOT NULL"):
            if cid in by_chembl:
                out[pref] = by_chembl[cid]
        con.close()
    return out


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--db", type=Path,
                    default=root / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--virus", default="NCBITaxon:2697049")
    ap.add_argument("--out", type=Path, default=root / "data/results")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    schema = load_schema()
    held = schema.held_out_predicates()
    print(f"held out of traversal (evaluation labels): {sorted(held)}")

    g = Graph.load(args.release, skip_predicates=held)
    print(f"{len(g.nodes):,} nodes | {len(g.drugs()):,} drugs | "
          f"{len(g.viruses())} viruses")

    # Ground truth: the cell-based antiviral edges we withheld.
    positives: set[str] = set()
    for line in (args.release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        # Only exact relations are ground truth: "EC50 > 10000 nM" is a
        # measurement of INACTIVITY and would invert the label.
        if (e["predicate"] == "HAS_ANTIVIRAL_ACTIVITY_AGAINST"
                and e["object"] == args.virus
                and (e.get("qualifiers") or {}).get("relation") == "="):
            positives.add(e["subject"])
    print(f"{len(positives):,} compounds with measured activity against "
          f"{args.virus} (ground truth)")

    print("\ncomputing DWPC over M1-M6")
    dwpc = dwpc_scores(g, schema.metapaths, schema.retrieval_limits)
    for name in sorted(dwpc):
        reach = sum(1 for by_virus in dwpc[name].values() if args.virus in by_virus)
        print(f"  {name}: {reach:,} drugs reach this virus")

    scorers = {
        "degree": degree_ranking(g),
        "dwpc": combine(dwpc, args.virus),
    }
    print("\ncomputing network proximity")
    prox = network_proximity(g, args.virus, schema.retrieval_limits)
    if prox:
        scorers["proximity"] = prox
    print(f"  {len(prox):,} drugs scored")

    print(f"\n{'scorer':<12} {'n':>7} {'Hits@10':>8} {'Hits@50':>8} "
          f"{'Hits@100':>9} {'MRR':>7}")
    ranked_all = {}
    for name, scores in scorers.items():
        r = rank(scores)
        ranked_all[name] = r
        print(f"{name:<12} {len(r):>7,} {hits_at_k(r, positives, 10):>8} "
              f"{hits_at_k(r, positives, 50):>8} {hits_at_k(r, positives, 100):>9} "
              f"{mrr(r, positives):>7.4f}")

    names = name_index(g, args.db)
    if names:
        print("\nbiology controls (rank out of n; lower is higher-ranked):")
        print(f"{'compound':<22} {'expect':<10} " +
              " ".join(f"{s:>12}" for s in ranked_all))
        for cname, kind in CONTROLS.items():
            nid = names.get(cname)
            if not nid:
                print(f"{cname:<22} {kind:<10} " +
                      " ".join(f"{'not in graph':>12}" for _ in ranked_all))
                continue
            cells = []
            for r in ranked_all.values():
                pos = next((i for i, (d, _v) in enumerate(r, 1) if d == nid), None)
                cells.append(f"{pos:>12,}" if pos else f"{'unscored':>12}")
            print(f"{cname:<22} {kind:<10} " + " ".join(cells))

    print(f"\ntop {args.top} by DWPC:")
    for i, (drug, score) in enumerate(ranked_all["dwpc"][:args.top], 1):
        contrib = sorted(((n, dwpc[n][drug][args.virus]) for n in dwpc
                          if drug in dwpc[n] and args.virus in dwpc[n][drug]),
                         key=lambda kv: -kv[1])
        routes = ",".join(f"{n}:{v:.3f}" for n, v in contrib[:3])
        label = next((k for k, v in names.items() if v == drug), None)
        mark = " *" if drug in positives else "  "
        print(f"{i:>3}{mark} {(label or drug)[:34]:<36} {score:.4f}  {routes}")
    print("  * = has measured activity (was held out of traversal)")

    args.out.mkdir(parents=True, exist_ok=True)
    payload = {
        "virus": args.virus,
        "held_out": sorted(held),
        "n_positives": len(positives),
        "metrics": {n: {"n": len(r),
                        "hits@10": hits_at_k(r, positives, 10),
                        "hits@50": hits_at_k(r, positives, 50),
                        "hits@100": hits_at_k(r, positives, 100),
                        "mrr": mrr(r, positives)}
                    for n, r in ranked_all.items()},
        "top": {n: [{"drug": d, "score": s} for d, s in r[:200]]
                for n, r in ranked_all.items()},
    }
    (args.out / f"baselines_{args.virus.replace(':', '_')}.json").write_text(
        json.dumps(payload, indent=2))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
