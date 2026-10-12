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

from kgav.baselines import Graph, dwpc_scores
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

    print("\n" + "=" * 70)
    print("This script measures the precondition. The decision rule is in")
    print("docs/prereg_m1_potency.md and is not restated here.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
