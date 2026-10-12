#!/usr/bin/env python
"""Is there enough measured potency to break M1's tie, or is there not?

THE PRECONDITION in docs/prereg_m1_potency.md, and nothing else. That file
states the defect -- M1 resolves 12 distinct values across 3,020 drugs, with
80.8% of the pool sharing one, because walk_metapath's weight is a product of
node degrees and nothing drug-specific ever enters it -- and makes letting
potency into the channel conditional on this measurement:

    If fewer than half of the nsp5 block's edges carry an exact value, this
    change is abandoned, and the honest statement is that the graph records
    which protein a compound hits and not how well.

So this counts, over the INHIBITS edges M1 ACTUALLY WALKS for one virus:

  - edges carrying an exact ic50_nm / ki_nm / kd_nm (relation "=" or absent)
  - edges whose relation is CENSORED (">", "<", ">=", "<=", "~")
  - edges carrying no potency field at all

Censored is kept separate and never counted as coverage. chembl.py's own rule
is that "IC50 > 10000 nM is a NEGATIVE result -- the compound was measured and
found inactive -- and coercing it to a potency number inverts its meaning", so
a censored bound cannot become the number that breaks a tie.

The breakdown is reported three ways: over the whole M1 pool, over the
LARGEST TIE BLOCK alone (which is the population the change has to separate),
and for the named controls, since prediction P2 turns on nirmatrelvir and
chloroquine separating.

USAGE

  python scripts/diagnose_m1_potency.py
  python scripts/diagnose_m1_potency.py --virus NCBITaxon:64320

Read-only. Writes nothing, recommends nothing, and decides nothing: the
decision rule is in the pre-registration.
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import Graph, channel_percentiles, dwpc_scores, potency_tiebreak
from kgav.schema import load_schema

POTENCY = ("ic50_nm", "ki_nm", "kd_nm")
CENSORED = (">", "<", ">=", "<=", "~")
# The named compounds prediction P2 turns on, by InChIKey so no ChEMBL
# database is needed to resolve them. Taken from the diagnose_controls run.
CONTROLS = {
    "INCHIKEY:LIENCHBZNNMNKG-OJFNHCPVSA-N": "nirmatrelvir  (positive)",
    "INCHIKEY:WHTVZRBIWZFKQO-UHFFFAOYSA-N": "chloroquine   (artifact)",
    "INCHIKEY:XXSMGPRMXLTPCZ-UHFFFAOYSA-N": "hydroxychloroq(artifact)",
    "INCHIKEY:RWWYLEGWBNMMLJ-YSOARWBDSA-N": "remdesivir    (positive)",
}


def classify(edge: dict) -> str:
    """exact | censored | none, for one INHIBITS edge.

    A potency value with a censored relation is NOT coverage: it is a bound,
    and on this project's own rule a bound in the wrong direction is a
    measurement of inactivity. An absent relation is treated as exact
    because that is how chembl.py emits an equality measurement.
    """
    q = edge.get("qualifiers") or {}
    has = [f for f in POTENCY if q.get(f) is not None]
    if not has:
        return "none"
    return "censored" if q.get("relation") in CENSORED else "exact"


def spearman(xs: list[float], ys: list[float]) -> float:
    """Rank correlation, ties taking their midrank. No scipy dependency.

    Used for prediction P3: potency exists preferentially for well-studied
    compounds, so the new ordering must not have become a restatement of
    node degree -- which is already the only scorer on this release whose
    interval clears chance, and which is the NULL HYPOTHESIS.
    """
    def ranks(v: list[float]) -> list[float]:
        order = sorted(range(len(v)), key=lambda i: v[i])
        out = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            mid = (i + j) / 2 + 1
            for k in order[i:j + 1]:
                out[k] = mid
            i = j + 1
        return out
    rx, ry = ranks(xs), ranks(ys)
    n = len(rx)
    if n < 2:
        return 0.0
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else 0.0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path, default=ROOT / "data/releases/v0.1")
    ap.add_argument("--virus", default="NCBITaxon:2697049")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        sys.exit(f"no release at {args.release}")

    schema = load_schema()
    g = Graph.load(args.release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    m1 = schema.metapaths["M1"]
    constraints = m1.get("constraints") or {}
    print(f"release {args.release.name} | virus {args.virus}")
    print(f"M1 constraints in force: {constraints or '(none)'}\n")

    # The M1 pool and its tie structure, from the same dwpc_scores the
    # ranking uses -- not recomputed here, because a second walk is a second
    # chance to walk it differently.
    dwpc = dwpc_scores(g, {"M1": m1}, schema.retrieval_limits)
    pool = {d: by[args.virus] for d, by in dwpc.get("M1", {}).items()
            if args.virus in by}
    if not pool:
        sys.exit(f"M1 reaches no drug for {args.virus}")
    vals = Counter(pool.values())
    modal_score, modal_n = vals.most_common(1)[0]
    block = {d for d, v in pool.items() if v == modal_score}
    print(f"M1 pool {len(pool):,} drugs | {len(vals):,} distinct scores")
    print(f"largest tie block {modal_n:,} drugs ({100 * modal_n / len(pool):.1f}%"
          f" of pool) at raw DWPC {modal_score:.6f}\n")

    # Which viral proteins route to this virus at all, so only the edges M1
    # can walk are counted. An edge to a protein with no onward route is not
    # part of this channel whatever its qualifiers say.
    reaching = set()
    for prot in g.nodes:
        if g.klass(prot) != "Protein" or not g.is_viral(prot):
            continue
        for gene, _e in g.neighbours("ENCODED_BY", prot):
            for tax, _e2 in g.neighbours("BELONGS_TO", gene):
                if tax == args.virus:
                    reaching.add(prot)
    print(f"{len(reaching):,} viral proteins route to this virus\n")

    want = constraints.get("INHIBITS.target_resolution")

    def tally(drugs) -> tuple[Counter, Counter]:
        """(edge counts by class, drug counts by best class they hold)."""
        edges: Counter = Counter()
        drugwise: Counter = Counter()
        order = {"exact": 0, "censored": 1, "none": 2}
        for d in drugs:
            best = "none"
            for nbr, edge in g.neighbours("INHIBITS", d):
                if nbr not in reaching:
                    continue
                if want is not None and (edge.get("qualifiers") or {}).get(
                        "target_resolution") != want:
                    continue
                c = classify(edge)
                edges[c] += 1
                if order[c] < order[best]:
                    best = c
            drugwise[best] += 1
        return edges, drugwise

    for title, drugs in (("WHOLE M1 POOL", set(pool)),
                         (f"LARGEST TIE BLOCK ({modal_n:,} drugs)", block)):
        edges, drugwise = tally(drugs)
        tot_e = sum(edges.values()) or 1
        tot_d = sum(drugwise.values()) or 1
        print("=" * 70)
        print(title)
        print("=" * 70)
        print(f"  {'':<10}{'edges':>10}{'% edges':>10}"
              f"{'drugs':>10}{'% drugs':>10}")
        for c in ("exact", "censored", "none"):
            print(f"  {c:<10}{edges[c]:>10,}{100 * edges[c] / tot_e:>9.1f}%"
                  f"{drugwise[c]:>10,}{100 * drugwise[c] / tot_d:>9.1f}%")
        print(f"  {'total':<10}{sum(edges.values()):>10,}{'':>10}"
              f"{sum(drugwise.values()):>10,}\n")
        if "TIE BLOCK" in title:
            share = 100 * edges["exact"] / tot_e
            print(f"  THE PRECONDITION: {share:.1f}% of this block's edges "
                  f"carry an exact value.")
            print(f"  The pre-registration abandons the change below 50%. "
                  f"{'MET' if share >= 50 else 'NOT MET'}.\n")

    print("=" * 70)
    print("THE CONTROLS (prediction P2 turns on these separating)")
    print("=" * 70)
    for nid, label in CONTROLS.items():
        if nid not in g.nodes:
            print(f"  {label}  not in the graph")
            continue
        here = "in M1 pool" if nid in pool else "NOT in M1 pool"
        tie = ", in the tie block" if nid in block else ""
        print(f"\n  {label}  [{here}{tie}]")
        for nbr, edge in g.neighbours("INHIBITS", nid):
            if nbr not in reaching:
                continue
            q = edge.get("qualifiers") or {}
            fam = ((g.nodes.get(nbr) or {}).get("properties") or {}).get(
                "protein_family") or nbr
            measured = ", ".join(f"{f}={q[f]:g}" for f in POTENCY
                                 if q.get(f) is not None) or "no value"
            rel = q.get("relation") or "="
            print(f"    {fam!s:<10} {classify(edge):<9} {rel:<3} "
                  f"{measured:<28} {nbr}")

    # ---------------------------------------- the predictions, scored
    print("\n" + "=" * 70)
    print("THE PREDICTIONS (docs/prereg_m1_potency.md, Amendment 1)")
    print("=" * 70)
    keys, fields = potency_tiebreak(g, args.virus, m1["path"], constraints)
    before = channel_percentiles(pool)
    after = channel_percentiles(pool, keys)
    in_pool = {d: k for d, k in keys.items() if d in pool}

    vals_after = Counter(round(v, 12) for v in after.values())
    big_after = vals_after.most_common(1)[0][1]
    print("\n  P1  tie structure")
    print(f"      distinct M1 scores   {len(vals):>7,}  ->  "
          f"{len(vals_after):>7,}   (predicted: past 1,000)")
    print(f"      largest tie block    {100 * modal_n / len(pool):>6.1f}%  ->  "
          f"{100 * big_after / len(pool):>6.1f}%   (predicted: under 5%)")
    p1 = len(vals_after) > 1000 and big_after / len(pool) < 0.05
    print(f"      P1 {'HOLDS' if p1 else 'FAILS'}")
    if len(pool) < 2000:
        # Said rather than silently rescaled: "past 1,000" is an ABSOLUTE
        # threshold written against v0.1's 3,020-drug pool, so on a smaller
        # pool P1 reads FAIL for a reason that is arithmetic rather than
        # about the change. The pre-registration calls P1 "evidence of
        # nothing except that the change took effect", so the threshold is
        # left as written rather than moved to fit.
        print(f"      note: P1's 1,000 is absolute and this pool holds "
              f"{len(pool):,}; read the two columns, not the verdict.")

    print("\n  P2' controls separate within the block")
    nirma = "INCHIKEY:LIENCHBZNNMNKG-OJFNHCPVSA-N"
    chloro = "INCHIKEY:WHTVZRBIWZFKQO-UHFFFAOYSA-N"
    block_keys = sorted(in_pool[d] for d in block if d in in_pool)
    def decile(d):
        if d not in in_pool or not block_keys:
            return None
        below = sum(1 for k in block_keys if k < in_pool[d])
        return below / len(block_keys)
    for nm, d in (("nirmatrelvir", nirma), ("chloroquine", chloro)):
        if d not in pool:
            print(f"      {nm:<14} not in the M1 pool")
            continue
        pc = decile(d)
        print(f"      {nm:<14} pct {before[d]:.4f} -> {after[d]:.4f}"
              + (f"   potency rank within block: {100 * pc:.1f}th pctile"
                 if pc is not None else "   no potency key"))
    gap = (after.get(chloro) is not None and after.get(nirma) is not None
           and after[nirma] - after[chloro])
    p2 = False
    if gap is not False:
        dn, dc = decile(nirma), decile(chloro)
        p2 = (gap >= 0.30 and dn is not None and dn >= 0.90
              and dc is not None and dc <= 0.50)
        print(f"      gap {gap:+.4f}   (predicted: >= 0.30, nirmatrelvir in "
              f"the top decile, chloroquine in the bottom half)")
    print(f"      P2' {'HOLDS' if p2 else 'FAILS'}")
    print(f"      ceiling, stated in advance: the block spans "
          f"{min(after[d] for d in block):.4f}-"
          f"{max(after[d] for d in block):.4f}, so nothing in it can exceed "
          f"that top.")

    print("\n  P3  the degree confound")
    drugs = sorted(pool)
    deg = [float(g.degree.get(d, 0)) for d in drugs]
    r_before = spearman([before[d] for d in drugs], deg)
    r_after = spearman([after[d] for d in drugs], deg)
    print(f"      rank corr with node degree   {r_before:+.4f}  ->  "
          f"{r_after:+.4f}   (delta {r_after - r_before:+.4f})")
    p3 = (r_after - r_before) <= 0.10
    print(f"      P3 {'HOLDS' if p3 else 'TRIPS -- the pre-registration says revert'}"
          f"   (bound: delta <= 0.10)")

    print("\n  P5  mixed quantity types, sized rather than assumed small")
    tot = sum(fields.values()) or 1
    for f, c in fields.most_common():
        print(f"      ordered on {f:<10} {c:>6,}  {100 * c / tot:>5.1f}%")

    print("\n" + "=" * 70)
    print(f"P1 {'HOLDS' if p1 else 'FAILS'} | P2' {'HOLDS' if p2 else 'FAILS'}"
          f" | P3 {'HOLDS' if p3 else 'TRIPS'}")
    print("The decision rule is in docs/prereg_m1_potency.md: keep only if P1")
    print("and P2' both hold and P3 does not trip. This script does not")
    print("decide, and claims nothing about AUC -- M1 reaches 3 of 89")
    print("measured actives before this change and 3 after.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
