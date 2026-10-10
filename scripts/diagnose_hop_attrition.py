#!/usr/bin/env python3
"""Where do the host-directed metapaths lose their candidates?

M1 reaches 56.6% of positives, M2 reaches 4.0%, and both start from the same
compound pool. One of the two hops M2 needs is starving it. This prints the
survivor count after every hop of M1 and M2 so the responsible layer is named
rather than guessed.

Run from the repo root:
    python diagnose_hop_attrition.py [RELEASE_DIR]
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


def find_release(arg: str | None) -> Path:
    if arg:
        p = Path(arg)
        if not (p / "edges.jsonl").exists():
            sys.exit(f"no edges.jsonl in {p}")
        return p
    cands = [d for d in sorted(Path("data/releases").iterdir())
             if (d / "edges.jsonl").exists() and (d / "nodes.jsonl").exists()
             and d.name != "smoke"]
    if not cands:
        sys.exit("no assembled release found under data/releases/ — pass the path")
    # prefer the one with no layer suffix (the assembled graph, not a layer)
    plain = [d for d in cands if "-" not in d.name.lstrip("v0123456789.")]
    pick = (plain or cands)[-1]
    if len(cands) > 1:
        print(f"candidates: {[d.name for d in cands]}")
    return pick


def main() -> int:
    rel = find_release(sys.argv[1] if len(sys.argv) > 1 else None)
    print(f"release: {rel}\n")

    # ---------------------------------------------------------------- nodes
    cls: dict[str, str] = {}
    viral: dict[str, bool] = {}
    with (rel / "nodes.jsonl").open() as fh:
        for line in fh:
            if not line.strip():
                continue
            n = json.loads(line)
            cls[n["id"]] = n["class"]
            if n["class"] == "Protein":
                viral[n["id"]] = bool((n.get("properties") or {}).get("is_viral"))

    n_cls = Counter(cls.values())
    print("nodes by class")
    for k, v in n_cls.most_common():
        print(f"  {k:<16} {v:>8,}")
    n_viral = sum(1 for v in viral.values() if v)
    print(f"  (Protein: {n_viral:,} viral / {len(viral) - n_viral:,} host)\n")

    # ---------------------------------------------------------------- edges
    pred = Counter()
    # hop-1 maps
    inhibits_viral: dict[str, set] = defaultdict(set)   # compound -> viral proteins
    targets_host: dict[str, set] = defaultdict(set)     # compound -> host proteins
    targets_any: dict[str, set] = defaultdict(set)
    # hop-2 maps
    hff_any: dict[str, set] = defaultdict(set)          # host protein -> virus
    hff_dep: dict[str, set] = defaultdict(set)          # dependency only
    hff_dir = Counter()

    with (rel / "edges.jsonl").open() as fh:
        for line in fh:
            if not line.strip():
                continue
            e = json.loads(line)
            p = e["predicate"]
            pred[p] += 1
            s, o = e["subject"], e["object"]
            if p == "INHIBITS" and cls.get(s) == "SmallMolecule":
                if viral.get(o):
                    inhibits_viral[s].add(o)
            elif p == "TARGETS" and cls.get(s) == "SmallMolecule":
                targets_any[s].add(o)
                if o in viral and not viral[o]:
                    targets_host[s].add(o)
            elif p == "HOST_FACTOR_FOR":
                d = (e.get("qualifiers") or {}).get("direction")
                hff_dir[d] += 1
                hff_any[s].add(o)
                if d == "dependency":
                    hff_dep[s].add(o)

    print("edges by predicate")
    for k, v in pred.most_common():
        print(f"  {k:<34} {v:>8,}")
    print()

    print("HOST_FACTOR_FOR direction qualifier")
    for k, v in hff_dir.most_common():
        print(f"  {str(k):<16} {v:>8,}")
    print()

    # ------------------------------------------------------------ attrition
    n_sm = n_cls.get("SmallMolecule", 0)

    def pct(x):
        return f"{100.0 * x / n_sm:5.1f}%" if n_sm else "   n/a"

    print(f"HOP ATTRITION  (denominator = all {n_sm:,} SmallMolecule nodes)\n")

    print("M1  direct-acting")
    print(f"  hop1  --INHIBITS-> viral Protein      {len(inhibits_viral):>7,}  {pct(len(inhibits_viral))}")
    print(f"        (distinct viral proteins hit:   {len({p for v in inhibits_viral.values() for p in v}):>7,})")
    print()

    print("M2  host-directed, single hop")
    print(f"  hop1  --TARGETS-> any Protein         {len(targets_any):>7,}  {pct(len(targets_any))}")
    print(f"  hop1  --TARGETS-> HOST Protein        {len(targets_host):>7,}  {pct(len(targets_host))}")
    reach_any = {c for c, ps in targets_host.items() if any(p in hff_any for p in ps)}
    reach_dep = {c for c, ps in targets_host.items() if any(p in hff_dep for p in ps)}
    print(f"  hop2  -> HOST_FACTOR_FOR (any dir)    {len(reach_any):>7,}  {pct(len(reach_any))}")
    print(f"  hop2  -> HOST_FACTOR_FOR dependency   {len(reach_dep):>7,}  {pct(len(reach_dep))}   <- M2's real pool")
    print()

    # which hop is the bottleneck
    host_proteins = {p for v in targets_host.values() for p in v}
    annotated_hf = {p for p in host_proteins if p in hff_dep}
    print("the two candidate bottlenecks")
    print(f"  distinct HOST proteins any drug targets          {len(host_proteins):>7,}")
    print(f"  ... of those, that are dependency host factors   {len(annotated_hf):>7,}")
    print(f"  distinct host factors in the graph (dependency)  {len(hff_dep):>7,}")
    overlap = len(annotated_hf)
    print()
    if not host_proteins:
        print("  VERDICT: the TARGETS layer is empty. Every host-directed metapath")
        print("           is starved at hop 1. Nothing downstream can matter.")
    elif overlap == 0:
        print("  VERDICT: druggable host proteins and screened host factors are")
        print("           DISJOINT sets. The two layers do not meet.")
    else:
        drug_side = len(targets_host) / n_sm if n_sm else 0
        print(f"  VERDICT: hop1 keeps {drug_side*100:.1f}% of compounds; hop2 then keeps")
        print(f"           {len(reach_dep)}/{len(targets_host)} of those "
              f"({100.0*len(reach_dep)/max(1,len(targets_host)):.1f}%).")
        print("           The smaller survival rate is the layer to invest in.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
