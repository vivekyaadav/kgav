#!/usr/bin/env python3
"""How many INHIBITS edges claim `biochemical` on a CELLULAR measurement?

    python scripts/diagnose_assay_type.py

The agent's first route for chloroquine is "inhibits the viral protein Spike
glycoprotein (IC50 160 nM, moderate)". Chloroquine does not bind spike at
160 nM. That record is near-certainly a pseudovirus ENTRY assay -- cell-based,
reading spike-mediated entry -- filed against "Spike glycoprotein" because
that is the assay's nominal target. Chloroquine blocks endosomal
acidification, so it inhibits the assay without touching the protein.

chembl.py SELECTS a.assay_type and a.description and then writes
quals["assay_type"] = "biochemical" on every viral-protein edge regardless.
The information needed to tell the two apart is fetched and discarded, which
is the H2 bug class: a layer asserting a property it did not measure.

ChEMBL's assay_type is B binding, F functional, A ADMET, T toxicity,
P physicochemical, U unclassified. B/F is binding-vs-functional, NOT
purified-vs-cellular -- an enzyme activity assay is F too -- so the code
alone cannot settle it. The description text can: "pseudotyped", "entry",
"infection", "cytopathic", a cell line name. This counts both and crosses
them, so the fix can be chosen against numbers rather than one example.

Read-only. Touches the ChEMBL SQLite file and nothing else.
"""
from __future__ import annotations

import argparse
import re
import sqlite3
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import yaml

ASSAY_TYPE_NAMES = {"B": "binding", "F": "functional", "A": "ADMET",
                    "T": "toxicity", "P": "physicochemical", "U": "unclassified"}

# Words that indicate the readout came from CELLS or whole virus rather than
# a purified protein. Deliberately specific: "assay" or "inhibition" would
# match everything and measure nothing.
CELLULAR_PATTERNS = (
    r"pseudo[- ]?typed", r"pseudovirus", r"pseudo[- ]?particle",
    r"\bentry\b", r"\binfect", r"cytopath", r"\bCPE\b", r"plaque",
    r"viral load", r"virus yield", r"\breplicon\b", r"\bluciferase\b",
    r"\bvero\b", r"calu-?3", r"\bhuh-?7", r"\bcaco-?2", r"\bhek ?293",
    r"\bA549\b", r"\bcells?\b",
)
CELLULAR_RE = re.compile("|".join(CELLULAR_PATTERNS), re.I)

# A purified-protein readout. Checked separately because a description can
# name both ("protease activity in Vero cell lysate") and the distinction
# then needs a human.
BIOCHEMICAL_PATTERNS = (
    r"recombinant", r"purified", r"enzymatic", r"\bFRET\b", r"fluorogenic",
    r"\bSPR\b", r"thermal shift", r"\bITC\b", r"cell[- ]free",
    # A lysate is a cell-derived but protein-level preparation, and its
    # description necessarily says "cell". Without this it reads as purely
    # cellular, when the honest answer is that it needs a human.
    r"lysate", r"\bextract\b", r"\bin vitro\b",
)
BIOCHEMICAL_RE = re.compile("|".join(BIOCHEMICAL_PATTERNS), re.I)


def classify(desc: str) -> str:
    cell = bool(CELLULAR_RE.search(desc or ""))
    bio = bool(BIOCHEMICAL_RE.search(desc or ""))
    if cell and bio:
        return "both_signals"
    if cell:
        return "cellular"
    if bio:
        return "biochemical"
    return "no_signal"


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path,
                    default=root / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--config", type=Path, default=root / "config/viruses.yaml")
    ap.add_argument("--examples", type=int, default=6)
    args = ap.parse_args()

    if not args.db.exists():
        print(f"not found: {args.db}\nPass --db with the path to chembl_37.db")
        return 1

    cfg = yaml.safe_load(args.config.read_text())
    taxa: set[str] = set()
    for v in cfg["viruses"]:
        taxa.add(str(v["taxon"]))
        for iso in v.get("isolate_taxa") or []:
            taxa.add(str(iso))

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    ph = ",".join("?" * len(taxa))
    # The same join chembl.py uses, narrowed to PROTEIN targets -- the rows
    # that become INHIBITS or MEASURED_INACTIVE_AGAINST edges.
    rows = con.execute(f"""
        SELECT md.chembl_id AS compound, md.pref_name AS drug_name,
               t.pref_name AS target_name, t.target_type,
               act.standard_type, act.standard_value, act.standard_units,
               a.assay_type, a.description, d.pubmed_id
        FROM activities act
        JOIN assays a            ON act.assay_id = a.assay_id
        JOIN target_dictionary t ON a.tid = t.tid
        JOIN molecule_dictionary md ON act.molregno = md.molregno
        LEFT JOIN docs d         ON a.doc_id = d.doc_id
        WHERE t.tax_id IN ({ph})
          AND act.standard_value IS NOT NULL
          AND act.standard_type IN ('IC50','EC50','Ki','Kd')
          AND t.target_type != 'ORGANISM'
    """, sorted(taxa)).fetchall()
    con.close()

    print(f"{len(rows):,} protein-target activity rows across "
          f"{len(taxa)} declared taxa")
    print("every one of these is written with assay_type 'biochemical'\n")

    by_code = Counter(r["assay_type"] for r in rows)
    print("ChEMBL assay_type, as the database records it")
    for code, n in by_code.most_common():
        name = ASSAY_TYPE_NAMES.get(code, "?")
        print(f"  {str(code):<4} {name:<16} {n:>7,}  {n/len(rows):6.1%}")

    desc = Counter(classify(r["description"]) for r in rows)
    print("\ndescription text, as a check on the code")
    for k, n in desc.most_common():
        print(f"  {k:<16} {n:>7,}  {n/len(rows):6.1%}")

    print("\ncrossed — the cell that matters is functional x cellular")
    cross: Counter = Counter()
    for r in rows:
        cross[(r["assay_type"], classify(r["description"]))] += 1
    codes = sorted({c for c, _ in cross}, key=lambda c: -by_code[c])
    kinds = ["cellular", "biochemical", "both_signals", "no_signal"]
    print(f"  {'':<5}" + "".join(f"{k:>14}" for k in kinds))
    for c in codes:
        line = f"  {str(c):<5}"
        for k in kinds:
            line += f"{cross[(c, k)]:>14,}"
        print(line)

    suspect = [r for r in rows
               if classify(r["description"]) == "cellular"]
    print(f"\n{len(suspect):,} rows ({len(suspect)/len(rows):.1%}) read as "
          f"CELLULAR and are filed as biochemical against a viral protein.")
    print("  Each becomes an INHIBITS edge the agent presents as a direct "
          "measurement\n  on that protein.")

    print(f"\nexamples ({min(args.examples, len(suspect))} of {len(suspect):,}):")
    for r in suspect[:args.examples]:
        val = r["standard_value"]
        print(f"  {(r['drug_name'] or r['compound'])[:22]:<24} "
              f"{r['target_name'][:26]:<28} "
              f"{r['standard_type']} {val:>9.1f} {r['standard_units'] or ''}"
              f"  [{r['assay_type']}]")
        print(f"      {(r['description'] or '')[:150]}")

    print("\nspike specifically — the agent's first route for chloroquine:")
    spike = [r for r in rows if "spike" in (r["target_name"] or "").lower()]
    sp_cell = [r for r in spike if classify(r["description"]) == "cellular"]
    print(f"  {len(spike):,} spike rows, {len(sp_cell):,} read as cellular")
    for r in spike:
        if "chloroquine" in ((r["drug_name"] or "") + r["compound"]).lower():
            print(f"  -> {r['drug_name'] or r['compound']}: "
                  f"{r['standard_type']} {r['standard_value']} "
                  f"{r['standard_units'] or ''} [{r['assay_type']}] "
                  f"PMID:{r['pubmed_id']}")
            print(f"     {(r['description'] or '')[:200]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
