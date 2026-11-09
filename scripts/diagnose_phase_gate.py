#!/usr/bin/env python3
"""What host_targets' max_phase gate costs, as three options rather than one.

    make diagnose-phase-gate          # or: python scripts/diagnose_phase_gate.py

FIX-LIST B1 PROPOSES DROPPING THE GATE, and the row estimated "<=1,409
compounds" with nothing in the repository able to check it. That figure
turned out to be counting neither of the two quantities that matter.

"DROP THE PHASE FILTER" IS AMBIGUOUS BY THREE ORDERS OF MAGNITUDE, because
the gate is `md.max_phase >= ?` and `NULL >= 0` is NULL in SQL. So
--min-phase 0 does NOT drop the gate: it keeps the exclusion of every
compound ChEMBL has never put in a trial, which is almost all of them.
The three readings are measured separately here.

THE QUANTITY THAT MATTERS IS NOT COMPOUNDS RECOVERED. A recovered compound
changes what the evaluation can measure only if the graph ALREADY has it and
it already carries antiviral labels -- otherwise it enlarges the candidate
pool, which is a different claim and not an AUC. Both are reported.

WHY IT IS A SCRIPT AND NOT A FIGURE IN A DOCUMENT. The B1 row's own
estimate had no way to be re-derived, which is how it stayed wrong; and the
numbers here move when ChEMBL does. Needs the 30 GB chembl_37 sqlite, which
is gitignored, so this cannot be a test.
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.host_targets import DEFAULT_MIN_PCHEMBL, DEFAULT_MIN_PHASE
from kgav.labels import build_labels
from kgav.pipeline import LAYERS, release_dir

#: Every gate in query_host_activities EXCEPT the phase one, which is what
#: this script varies. Kept in one string so the two queries below cannot
#: disagree about what else is being held constant.
_WITHOUT_PHASE = """
    FROM activities act
    JOIN assays a            ON act.assay_id = a.assay_id
    JOIN target_dictionary t ON a.tid = t.tid
    JOIN molecule_dictionary md ON act.molregno = md.molregno
    LEFT JOIN compound_structures st ON md.molregno = st.molregno
    WHERE t.tax_id = 9606
      AND t.target_type = 'SINGLE PROTEIN'
      AND act.pchembl_value >= ?
      AND act.standard_type IN ('IC50','EC50','Ki','Kd')
"""


def _bucket_sql(min_phase: float) -> str:
    return (f"CASE WHEN md.max_phase IS NULL THEN 'null' "
            f"WHEN md.max_phase >= {min_phase:g} THEN 'kept' "
            f"ELSE 'below' END")


def buckets(db: Path, min_pchembl: float, min_phase: int):
    """Compounds and activity rows either side of the gate.

    Keyed by InChIKey, because that is what becomes a node
    (`host_targets.py`: `drug = f"INCHIKEY:{inchikey}"`), and rows with no
    InChIKey are dropped by the ingest, so they are dropped here too and
    counted separately rather than silently.
    """
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    b = _bucket_sql(min_phase)
    keys: dict[str, set[str]] = {"kept": set(), "below": set(), "null": set()}
    rows, keyless = Counter(), Counter()
    for k, bucket, n in con.execute(
            f"SELECT st.standard_inchi_key, {b}, COUNT(*) "
            f"{_WITHOUT_PHASE} GROUP BY st.standard_inchi_key, {b}",
            (min_pchembl,)):
        rows[bucket] += n
        if k:
            keys[bucket].add(f"INCHIKEY:{k}")
        else:
            keyless[bucket] += n
    con.close()
    return keys, rows, keyless


def proteins_behind(db: Path, min_pchembl: float, min_phase: int) -> set[str]:
    """Human proteins the rows BELOW the gate would reach."""
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    q = """SELECT DISTINCT cs.accession
            FROM activities act
            JOIN assays a            ON act.assay_id = a.assay_id
            JOIN target_dictionary t ON a.tid = t.tid
            JOIN target_components tc ON t.tid = tc.tid
            JOIN component_sequences cs ON tc.component_id = cs.component_id
            JOIN molecule_dictionary md ON act.molregno = md.molregno
            WHERE t.tax_id = 9606
              AND t.target_type = 'SINGLE PROTEIN'
              AND act.pchembl_value >= ?
              AND act.standard_type IN ('IC50','EC50','Ki','Kd')
              AND (md.max_phase IS NULL OR md.max_phase < ?)"""
    out = {f"UniProtKB:{a}" for (a,) in con.execute(q, (min_pchembl, min_phase))
           if a}
    con.close()
    return out


def release_view(release: Path):
    """Nodes, TARGETS subjects, host factors and labels, read once."""
    nodes, small = set(), set()
    for line in open(release / "nodes.jsonl"):
        n = json.loads(line)
        nodes.add(n["id"])
        if n["class"] == "SmallMolecule":
            small.add(n["id"])
    targets, host_factors = set(), set()
    for line in open(release / "edges.jsonl"):
        e = json.loads(line)
        if e["predicate"] == "TARGETS":
            targets.add(e["subject"])
        elif e["predicate"] == "HOST_FACTOR_FOR":
            host_factors.add(e["subject"])
    L = build_labels(release)
    pos = {c for cs in L.positives.values() for c in cs}
    neg = {c for cs in L.negatives.values() for c in cs}
    return nodes, small, targets, host_factors, pos, neg


def layer_contributions(release_root: Path, compounds: set[str]):
    """Which layer first knows each compound, and whether it runs in time.

    Option D ("any phase, but only for compounds the graph already has")
    can only see layers built BEFORE hosttargets, so this decides whether
    it is a flag or a pipeline reordering.
    """
    order = [l.name for l in LAYERS]
    cut = order.index("hosttargets")
    seen, out = set(), []
    for name in order:
        f = release_root / f"v0.1-{name}" / "nodes.jsonl"
        if not f.exists():
            continue
        ids = {json.loads(l)["id"] for l in open(f)
               if json.loads(l).get("class") == "SmallMolecule"}
        here = (compounds & ids) - seen
        if here:
            out.append((name, len(here), order.index(name) < cut))
            seen |= here
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data/raw/chembl/chembl_37/"
                                   "chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--min-pchembl", type=float, default=DEFAULT_MIN_PCHEMBL)
    ap.add_argument("--min-phase", type=int, default=DEFAULT_MIN_PHASE)
    args = ap.parse_args(argv)

    if not args.db.exists():
        print(f"no ChEMBL sqlite at {args.db}")
        return 1
    if not args.release.exists():
        print(f"no release at {args.release}")
        return 1

    keys, rows, keyless = buckets(args.db, args.min_pchembl, args.min_phase)
    nodes, small, targets, hf, pos, neg = release_view(args.release)
    labelled = pos | neg

    print("=" * 74)
    print(f"THE GATE  (max_phase >= {args.min_phase}, "
          f"pChEMBL >= {args.min_pchembl:g})")
    print("=" * 74)
    print(f"{'bucket':<10} {'activity rows':>14} {'no InChIKey':>12} "
          f"{'compounds':>11}")
    for b, label in (("kept", "kept"), ("below", f"< {args.min_phase}"),
                     ("null", "NULL")):
        print(f"{label:<10} {rows[b]:>14,} {keyless[b]:>12,} {len(keys[b]):>11,}")
    added = rows["below"] + rows["null"]
    print(f"\n  rows a dropped gate would add: {added:,}")
    print("  host_targets.gate_census records this independently as "
          "`dropped:\n  max_phase below 2`; if the two disagree, one of them "
          "is wrong.")

    recovered = (keys["below"] | keys["null"]) - keys["kept"]
    options = {
        "B  --min-phase 0 (early RECORDED phase only)":
            (keys["below"] - keys["kept"], rows["below"]),
        "C  no phase gate at all (adds max_phase IS NULL)":
            (recovered, added),
        "D  any phase, but only compounds already in the graph":
            (recovered & nodes, None),
    }
    print()
    print("=" * 74)
    print("THREE READINGS OF B1")
    print("=" * 74)
    print(f"{'option':<52} {'rows+':>11} {'cmpds+':>8} {'new':>8} "
          f"{'routed':>7} {'lab':>6} {'pos':>5}")
    print("-" * 102)
    for name, (ks, r) in options.items():
        gain = ks - targets
        print(f"{name:<52} {(f'{r:,}' if r is not None else '~'):>11} "
              f"{len(ks):>8,} {len(ks - nodes):>8,} {len(gain & nodes):>7,} "
              f"{len(gain & labelled):>6,} {len(gain & pos):>5,}")
    print("-" * 102)
    print("  cmpds+  distinct InChIKeys recovered")
    print(f"  new     of those, absent from the release "
          f"(which holds {len(small):,} small molecules)")
    print("  routed  recovered compounds ALREADY in the graph that gain a "
          "first TARGETS edge")
    print("  lab/pos of the routed, how many carry labels / positive labels")
    print("\n  C AND D ROUTE THE SAME COMPOUNDS, because `routed` is by "
          "definition\n  the in-graph subset. C's extra nodes buy candidate "
          "pool, not evaluation.")

    routed = (recovered & nodes) - targets
    print()
    print("=" * 74)
    print("CAN OPTION D SEE THEM IN TIME?")
    print("=" * 74)
    for name, n, before in layer_contributions(args.release.parent, routed):
        when = "before hosttargets" if before else "AFTER hosttargets"
        print(f"  {name:<14} {n:>6,} compounds   ({when})")
    print("\n  A layer built after hosttargets is invisible to a phase-1 "
          "filter, so\n  those compounds need the layer order changed, not "
          "just a flag.")

    acc = proteins_behind(args.db, args.min_pchembl, args.min_phase)
    print()
    print("=" * 74)
    print("WHAT THE RECOVERED ROWS REACH")
    print("=" * 74)
    print(f"  host dependency/restriction proteins : {len(hf):>7,}")
    print(f"  proteins the recovered rows hit      : {len(acc):>7,}")
    print(f"  ...already nodes in the graph        : {len(acc & nodes):>7,}")
    print(f"  ...that are host factors (M2 routes) : {len(acc & hf):>7,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
