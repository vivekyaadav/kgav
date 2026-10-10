#!/usr/bin/env python3
"""Two questions the hop-attrition run raised but could not answer.

1. M2's structural pool is 12.1% of all compounds, yet the evaluation found
   only 4.0% of POSITIVES reachable. Either the labelled compounds are a
   biased sample of the graph, or they are not. This counts it directly.

2. INHIBITS reaches only 19 distinct viral proteins, and M1 is the channel
   that scores. If a compound's INHIBITS edge and its antiviral-activity
   label come from the SAME publication, M1 is partly reading the label it
   is being scored against. This measures that overlap.

Read-only. Run from the repo root:
    python diagnose_label_overlap.py [RELEASE_DIR]
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

POS = "HAS_ANTIVIRAL_ACTIVITY_AGAINST"
NEG = "MEASURED_INACTIVE_AGAINST"


def main() -> int:
    rel = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/releases/v0.1")
    if not (rel / "edges.jsonl").exists():
        sys.exit(f"no edges.jsonl in {rel}")
    print(f"release: {rel}\n")

    cls: dict[str, str] = {}
    viral: dict[str, bool] = {}
    label: dict[str, str] = {}
    with (rel / "nodes.jsonl").open() as fh:
        for line in fh:
            if not line.strip():
                continue
            n = json.loads(line)
            cls[n["id"]] = n["class"]
            props = n.get("properties") or {}
            label[n["id"]] = props.get("label") or n["id"]
            if n["class"] == "Protein":
                viral[n["id"]] = bool(props.get("is_viral"))

    inh: dict[str, set] = defaultdict(set)       # compound -> viral proteins
    tgt: dict[str, set] = defaultdict(set)       # compound -> host proteins
    hff_dep: dict[str, set] = defaultdict(set)   # host protein -> virus
    pos_c: set = set()
    neg_c: set = set()
    # publications per (compound, predicate)
    pubs: dict[tuple, set] = defaultdict(set)
    viral_hit = Counter()

    with (rel / "edges.jsonl").open() as fh:
        for line in fh:
            if not line.strip():
                continue
            e = json.loads(line)
            p, s, o = e["predicate"], e["subject"], e["object"]
            pl = set(e.get("publications") or [])
            if p == "INHIBITS" and cls.get(s) == "SmallMolecule" and viral.get(o):
                inh[s].add(o)
                viral_hit[o] += 1
                pubs[(s, "INHIBITS")] |= pl
            elif p == "TARGETS" and cls.get(s) == "SmallMolecule" and o in viral and not viral[o]:
                tgt[s].add(o)
            elif p == "HOST_FACTOR_FOR" and (e.get("qualifiers") or {}).get("direction") == "dependency":
                hff_dep[s].add(o)
            elif p == POS:
                pos_c.add(s)
                pubs[(s, POS)] |= pl
            elif p == NEG:
                neg_c.add(s)

    n_sm = sum(1 for c in cls.values() if c == "SmallMolecule")
    inh_c, tgt_c = set(inh), set(tgt)
    m2_pool = {c for c, ps in tgt.items() if any(p in hff_dep for p in ps)}

    # ------------------------------------------------- Q1: population overlap
    print("=" * 66)
    print("Q1  are the labelled compounds a biased sample of the graph?")
    print("=" * 66)
    print(f"\nall SmallMolecule nodes                  {n_sm:>7,}")
    print(f"  with INHIBITS -> viral protein         {len(inh_c):>7,}   {100*len(inh_c)/n_sm:5.1f}%")
    print(f"  with TARGETS  -> host protein          {len(tgt_c):>7,}   {100*len(tgt_c)/n_sm:5.1f}%")
    print(f"  in M2's pool (reaches a dep. factor)   {len(m2_pool):>7,}   {100*len(m2_pool)/n_sm:5.1f}%")
    print(f"\n  BOTH INHIBITS and TARGETS              {len(inh_c & tgt_c):>7,}")
    print(f"  INHIBITS only                          {len(inh_c - tgt_c):>7,}")
    print(f"  TARGETS only                           {len(tgt_c - inh_c):>7,}")
    print(f"  neither                                {n_sm - len(inh_c | tgt_c):>7,}")

    for nm, s in (("POSITIVES (%s)" % POS, pos_c), ("NEGATIVES (%s)" % NEG, neg_c)):
        if not s:
            continue
        print(f"\n{nm}: {len(s):,} compounds")
        for lbl, sub in (("reach M1 hop1 (INHIBITS viral)", inh_c),
                         ("have TARGETS -> host protein", tgt_c),
                         ("in M2's pool (dep. host factor)", m2_pool)):
            k = len(s & sub)
            print(f"  {lbl:<34} {k:>6,}   {100*k/len(s):5.1f}%")
        exp = 100 * len(m2_pool) / n_sm
        obs = 100 * len(s & m2_pool) / len(s)
        print(f"  -> M2 reach {obs:.1f}% observed vs {exp:.1f}% if unbiased"
              f"   ({'UNDER' if obs < exp else 'OVER'}-represented {exp/max(obs,0.01):.1f}x)")

    # --------------------------------------------- Q2: same-document leakage
    print("\n" + "=" * 66)
    print("Q2  does M1 read the label it is scored against?")
    print("=" * 66)
    print(f"\nINHIBITS reaches {len(viral_hit)} distinct viral proteins. Top 12:")
    for pid, n in viral_hit.most_common(12):
        print(f"  {label.get(pid, pid)[:44]:<46} {n:>6,}")

    both = sorted(inh_c & pos_c)
    shared = same = 0
    for c in both:
        a, b = pubs.get((c, "INHIBITS"), set()), pubs.get((c, POS), set())
        if a and b:
            shared += 1
            if a & b:
                same += 1
    print(f"\ncompounds with BOTH an INHIBITS edge and a {POS} label: {len(both):,}")
    print(f"  both sides carry publications                        {shared:,}")
    print(f"  and they SHARE at least one publication              {same:,}"
          + (f"   ({100*same/shared:.1f}%)" if shared else ""))
    print("\n  A shared PMID means the binding measurement M1 traverses and the")
    print("  cell-based label it is scored against came from one paper. That is")
    print("  not wrong data -- it is the same experiment counted twice, and it")
    print("  inflates M1 in a way no compound-level split removes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
