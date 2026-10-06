#!/usr/bin/env python3
"""Ingest the NCATS OpenData SARS-CoV-2 CPE screen.

    curl -L -o "data/raw/ncats/cpe.tsv" \\
      "https://opendata.ncats.nih.gov/public/odp/SARS-CoV-2_cytopathic_effect_(CPE).tsv"

    python scripts/ingest_ncats.py --inspect      # show the header mapping only
    python scripts/ingest_ncats.py                # build the layer

--inspect first, always. The mapping from this file's header onto the fields
this layer needs is the one thing that cannot be verified by a schema check: a
mis-mapped potency column produces a release that validates, counts plausibly
and encodes the wrong measurement.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.emit import Emit
from kgav.ncats import (
    COLUMN_CANDIDATES,
    CPE_PUBLICATION,
    DEFAULT_DATE,
    SARS_COV_2,
    SOURCE,
    ColumnError,
    ingest_cpe,
    join_report,
    read_rows,
    unit_for,
)
from kgav.schema import load_schema

# Below this, assume the normalization is wrong rather than the library being
# unusual. The CPE library is approved and investigational compounds; ChEMBL's
# coronavirus set and DrugCentral overlap it heavily, so a few per cent join
# rate means the keys are being built differently, not that the overlap is
# genuinely tiny.
MIN_PLAUSIBLE_JOIN_RATE = 0.05


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", type=Path,
                    default=root / "data/raw/ncats/cpe.tsv")
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1",
                    help="existing graph, read only, to measure the join rate")
    ap.add_argument("--out", type=Path, default=root / "data/releases/v0.1-ncats")
    ap.add_argument("--date", default=DEFAULT_DATE,
                    help=f"first_asserted_date for every edge (default "
                         f"{DEFAULT_DATE}, the bioRxiv preprint). Never "
                         f"inferred from the file's mtime, which is a "
                         f"re-export date.")
    ap.add_argument("--ac50-units", default=None,
                    choices=["M", "uM", "nM", "logM"],
                    help="override the unit inferred from the column name")
    ap.add_argument("--cc50-units", default=None,
                    choices=["M", "uM", "nM", "logM"])
    ap.add_argument("--max-conc-nm", type=float, default=46_000.0,
                    help="highest concentration the screen reached. An "
                         "inactive row is emitted as '> this'. Wrong here and "
                         "every negative carries a wrong censoring point.")
    ap.add_argument("--inspect", action="store_true",
                    help="print the header, the resolved mapping and a few "
                         "parsed rows, then stop without writing anything")
    args = ap.parse_args()

    if not args.tsv.exists():
        print(f"not found: {args.tsv}\n\nDownload it with:\n"
              f'  mkdir -p {args.tsv.parent}\n'
              f'  curl -L -o "{args.tsv}" \\\n'
              f'    "https://opendata.ncats.nih.gov/public/odp/'
              f'SARS-CoV-2_cytopathic_effect_(CPE).tsv"')
        return 1

    try:
        rows, cols = read_rows(args.tsv)
    except ColumnError as e:
        print(f"HEADER NOT UNDERSTOOD\n\n{e}")
        return 1

    ac50_col = cols.get("ac50") or cols.get("lac50")
    ac50_units = args.ac50_units or (unit_for(ac50_col) if ac50_col else "uM")
    cc50_units = args.cc50_units or (unit_for(cols["cc50"]) if "cc50" in cols else "uM")

    print(f"{args.tsv}  ({len(rows):,} rows)\n")
    print("column mapping  (field <- this file's header)")
    for field in COLUMN_CANDIDATES:
        got = cols.get(field)
        mark = " " if got else "!"
        print(f"  {mark} {field:<10} <- {got if got else '(not found)'}")
    print(f"\n  AC50 read as {ac50_units}, CC50 as {cc50_units}")
    print(f"  inactive rows emitted as '> {args.max_conc_nm:,.0f} nM'")
    print(f"  every edge cites {CPE_PUBLICATION}, dated {args.date}")
    if not ac50_col:
        print("\n  ! NO POTENCY COLUMN MATCHED. Every row will fall back to the "
              "outcome\n    column, so no active will carry an EC50 and no "
              "selectivity index can be\n    computed. Add the real spelling "
              "to COLUMN_CANDIDATES before building.")

    if args.inspect:
        from kgav.ncats import classify_row
        print("\nfirst 5 rows as parsed:")
        for r in rows[:5]:
            v, q, why = classify_row(r, cols, ac50_units=ac50_units,
                                     cc50_units=cc50_units,
                                     max_conc_nm=args.max_conc_nm)
            sm = (r.get(cols["smiles"]) or "")[:30]
            print(f"  {str(v):<9} {why:<26} {json.dumps(q, default=str)[:90]}")
            print(f"            smiles={sm}")
        print("\n--inspect: nothing written. Re-run without it to build.")
        return 0

    try:
        from kgav.normalize_chem import normalize_chemical
    except ImportError as e:
        print(f"\nrdkit is required to normalize structures: {e}")
        return 1

    # Read the existing compounds FIRST: ingest_cpe needs them to decide
    # which compounds it must not write a node for, not merely to report on.
    existing: set[str] = set()
    if (args.release / "nodes.jsonl").exists():
        for line in (args.release / "nodes.jsonl").read_text().splitlines():
            if line.strip():
                n = json.loads(line)
                if n["class"] == "SmallMolecule":
                    existing.add(n["id"])
        print(f"  {len(existing):,} compounds already in {args.release.name}")
    else:
        print(f"  ! {args.release} not found -- every compound will be treated "
              f"as new,\n    so this layer will assert is_approved on "
              f"compounds another layer knows better.")

    em = Emit()
    print("\nnormalizing and ingesting (rdkit, this takes a minute)")
    stats = ingest_cpe(em, rows, cols, normalize_chemical,
                       existing_compounds=existing,
                       taxon=SARS_COV_2, source=SOURCE, date=args.date,
                       ac50_units=ac50_units, cc50_units=cc50_units,
                       max_conc_nm=args.max_conc_nm)

    print(f"\n{stats['active']:,} active | {stats['inactive']:,} measured INACTIVE")
    print(f"  {stats['with_same_plate_si']:,} carry a same-plate selectivity index")
    print("\n  per-row outcome:")
    for k in sorted(k for k in stats if k.startswith("reason_")):
        print(f"    {k.replace('reason_', ''):<28} {stats[k]:>7,}")
    for k in ("no_structure", "unnormalizable", "duplicate_structure"):
        if stats[k]:
            print(f"    {k:<28} {stats[k]:>7,}")

    # --------------------------------------------------- the join rate
    labelled = {e["subject"] for e in em.edges}
    jr = join_report(labelled, existing)
    print(f"\njoin against {args.release.name} ({len(existing):,} compounds)")
    print(f"  compounds labelled {jr['emitted']:>7,}")
    print(f"  already in graph   {jr['already_in_graph']:>7,}   "
          f"{jr['join_rate']:.1%}")
    print(f"  new compounds      {jr['new_compounds']:>7,}")
    print("\n  A compound that joins carries the graph's existing TARGETS "
          "edges, so it can\n  be reached by M2-M6. One that does not is a "
          "label with no path to it.")

    if existing and jr["join_rate"] < MIN_PLAUSIBLE_JOIN_RATE:
        print(f"\nREFUSING TO WRITE — join rate {jr['join_rate']:.1%} is below "
              f"{MIN_PLAUSIBLE_JOIN_RATE:.0%}.\n"
              f"  This library is approved and investigational compounds, which "
              f"ChEMBL and\n  DrugCentral cover heavily. A rate this low means "
              f"the keys are being built\n  differently from the rest of the "
              f"graph, not that the overlap is small --\n  which is the silent "
              f"failure this layer exists to avoid. Compare a few of\n  these "
              f"InChIKeys against the release before overriding.")
        return 1

    violations = load_schema().validate_batch(list(em.nodes.values()), em.edges)
    if violations:
        print(f"\nSCHEMA FAIL — {len(violations)} violations")
        for v in violations[:8]:
            print(f"  {v}")
        return 1
    print("\nschema: PASS")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in em.nodes.values()))
    (args.out / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in em.edges))
    print(f"{len(em.nodes):,} compounds, {len(em.edges):,} edges -> {args.out}")
    print("\nNext: add this layer to the assembler's LAYER_PRECEDENCE, "
          "reassemble, then\n  python scripts/hard_negatives.py --selectivity "
          "verified-only \\\n      --selectivity-applies both "
          "--publication-disjoint")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
