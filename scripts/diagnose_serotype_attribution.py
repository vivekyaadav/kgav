#!/usr/bin/env python
"""Does the assay description name the virus the row is filed under?

WHY THIS EXISTS

A row's taxon comes entirely from ChEMBL's target_dictionary.tax_id --
src/kgav/chembl.py:253 is `WHERE t.tax_id IN (...)` and the assay description
gets no vote. That is the right default: ChEMBL's target taxonomy is curated
and a description is free text.

But diagnose_flavivirus_targets printed these among DENV-2's protein rows, all
filed under taxa 11060 and 31634, on 2026-10-11:

    Inhibition of dengue virus-4 NS2B-NS3 protease
    Inhibition of dengue virus-3 NS2B-NS3 protease
    Inhibition of dengue virus-1 NS2B-NS3 protease
    Non-competitive inhibition of Dengue virus 4 NS3 helicase
    Non-competitive inhibition of Dengue virus 3 NS3 helicase

In the same run DENV-1, DENV-3 and DENV-4 each reported ZERO rows reaching the
chain lookup. So those serotypes' protein-level evidence may be sitting in
DENV-2's bucket -- which would inflate DENV-2 and empty the others, and put
INHIBITS edges on the wrong species.

That matters more here than it would in most graphs. M1 is the direct-acting
channel: an edge saying "this compound inhibits DENV-2 NS3" when the
measurement was made on DENV-4 NS3 is a wrong-species claim on the one channel
this project has a defensible result for. It is also the activity-side twin of
the trap config/viruses.yaml already guards on the protein side with
chain_source and isolate_taxa.

WHAT IT DOES NOT CLAIM

Five sampled strings show the phenomenon exists. They say nothing about its
size, and a sample read as a rate is a mistake this project has made before --
"DENV-2's single largest assay" was once asserted from two sampled
descriptions and turned out to be twelve rows. So this counts every row and
reports the rate, the affected targets, and the rows it cannot classify.

A mismatch here is NOT automatically an error. Three readings, and this script
does not choose between them:

  cross-reactivity panel   one paper measures a compound against all four
                           serotypes; ChEMBL may file the series against one
                           representative target.
  coarse target curation   the target is really "dengue virus protease" and
                           the tax_id is the genus-ish stand-in.
  genuine misfiling        the measurement belongs to another serotype.

Only the third is a data error, and telling them apart needs the target
record, so the target chembl_id and name are printed for every mismatch.

USAGE

  python scripts/diagnose_serotype_attribution.py
  python scripts/diagnose_serotype_attribution.py --family Coronaviridae
  python scripts/diagnose_serotype_attribution.py --samples 20

Read-only. Writes nothing, recommends nothing.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.chembl import assay_readout, query_activities  # noqa: E402
from kgav.virus_register import declared_taxa  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = (ROOT / "data" / "raw" / "chembl" / "chembl_37"
              / "chembl_37_sqlite" / "chembl_37.db")

# WHICH VIRUS A DESCRIPTION NAMES.
#
# Anchored on the virus name, never on a bare serotype number. The real WNV
# rows include "West Nile virus serotype 2 NS2B-NS3 protease", so a pattern
# like serotype\s*([1-4]) would read a West Nile row as dengue serotype 2 --
# inventing exactly the mismatch this script is looking for. Every dengue
# pattern therefore requires "dengue" or "denv" adjacent to the digit.
#
# Each entry maps to the register's species taxon.
NAME_PATTERNS: list[tuple[re.Pattern, str, str]] = [
    (re.compile(r"dengue\s+virus(?:\s+(?:type|serotype))?[-\s]*1(?!\d)"
                r"|\bdenv[-\s]?1(?!\d)", re.I), "11053", "DENV-1"),
    (re.compile(r"dengue\s+virus(?:\s+(?:type|serotype))?[-\s]*2(?!\d)"
                r"|\bdenv[-\s]?2(?!\d)", re.I), "11060", "DENV-2"),
    (re.compile(r"dengue\s+virus(?:\s+(?:type|serotype))?[-\s]*3(?!\d)"
                r"|\bdenv[-\s]?3(?!\d)", re.I), "11069", "DENV-3"),
    (re.compile(r"dengue\s+virus(?:\s+(?:type|serotype))?[-\s]*4(?!\d)"
                r"|\bdenv[-\s]?4(?!\d)", re.I), "11070", "DENV-4"),
    (re.compile(r"\bzika\b|\bzikv\b", re.I), "64320", "Zika"),
    (re.compile(r"west\s+nile|\bwnv\b", re.I), "11082", "West Nile"),
    (re.compile(r"yellow\s+fever|\byfv\b", re.I), "11089", "Yellow fever"),
    (re.compile(r"japanese\s+encephalitis|\bjev\b", re.I), "11072", "JEV"),
    (re.compile(r"tick[-\s]?borne\s+encephalitis|\btbev\b", re.I),
     "11084", "TBEV"),
]

# Names that identify no single virus. Counted apart from "named nothing",
# because a generic label is a curation choice rather than a silence -- the
# 1,013-row PubChem screen is titled "Flavivirus Genomic Capping Enzyme".
GENERIC = re.compile(r"\bflavivirus(es)?\b"
                     r"|dengue\s+virus(?!\s*(?:type|serotype)?[-\s]*[1-4])"
                     r"|\bdenv\b(?![-\s]?[1-4])", re.I)


def viruses_named(description: str) -> list[str]:
    """Every register virus this description names, by label."""
    return [label for pat, _tax, label in NAME_PATTERNS
            if pat.search(description or "")]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--config", type=Path,
                    default=ROOT / "config" / "viruses.yaml")
    ap.add_argument("--family", default="Flaviviridae")
    ap.add_argument("--samples", type=int, default=10)
    args = ap.parse_args(argv)

    if not args.db.exists():
        sys.exit(f"no ChEMBL db at {args.db}")

    cfg = yaml.safe_load(args.config.read_text())
    default_family = (cfg.get("defaults") or {}).get("family")
    maps = {sec: declared_taxa(cfg, sec) for sec in ("viruses", "surveyed")}

    wanted = [(v, sec) for sec in ("viruses", "surveyed")
              for v in (cfg.get(sec) or [])
              if (v.get("family") or default_family) == args.family]
    if not wanted:
        sys.exit(f"no virus in either section has family {args.family!r}")

    print(f"{len(wanted)} {args.family} virus(es)\n"
          f"taxon comes from ChEMBL target_dictionary.tax_id; this compares it\n"
          f"against what the assay description says\n")

    grand = Counter()
    for v, section in wanted:
        taxa = [t for t, sp in maps[section].items() if sp == str(v["taxon"])]
        if not taxa:
            print(f"BUG: no taxa for {v['name']} in {section!r}; refusing to "
                  f"query (an empty list returns zero rows via `IN ()`)\n")
            continue

        rows = query_activities(args.db, sorted(taxa))
        # The protein path only: an ORGANISM-target row becomes
        # HAS_ANTIVIRAL_ACTIVITY_AGAINST on the taxon and never reaches a
        # chain, and a cellular readout is emitted on the organism too. Those
        # are labels, not traversable evidence, so a mismatch there does not
        # put a wrong-species edge on a protein.
        protein_rows = [
            r for r in rows
            if (r.get("target_type") or "") != "ORGANISM"
            and "/" not in (r.get("target_name") or "")
            and assay_readout(r.get("description") or "",
                              r.get("assay_type"))[0] == "protein"]

        print("=" * 78)
        print(f"{v['name']}  [{section}]  (taxa {', '.join(sorted(taxa))})")
        print("=" * 78)
        print(f"  {len(rows):,} activity rows, {len(protein_rows):,} on the "
              f"protein path")
        if not protein_rows:
            print("  nothing on the protein path to attribute\n")
            continue

        verdict = Counter()
        mismatched: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
        targets: Counter = Counter()
        for r in protein_rows:
            desc = r.get("description") or ""
            named = viruses_named(desc)
            if not named:
                verdict["generic flavivirus name" if GENERIC.search(desc)
                        else "names no virus"] += 1
            elif named == [v["name"]]:
                verdict["names this virus"] += 1
            elif v["name"] in named:
                verdict["names this virus AND others (panel)"] += 1
                mismatched["panel"].append(
                    (desc, r.get("target") or "?", r.get("target_name") or "?"))
            else:
                verdict["names a DIFFERENT virus"] += 1
                key = "+".join(named)
                mismatched[key].append(
                    (desc, r.get("target") or "?", r.get("target_name") or "?"))
                targets[f"{r.get('target')}  {r.get('target_name')}"] += 1

        n = len(protein_rows)
        print(f"\n  WHAT THE DESCRIPTION NAMES (of {n:,} protein rows)")
        for key, c in verdict.most_common():
            flag = "  <-- WRONG-SPECIES RISK" if "DIFFERENT" in key else ""
            print(f"    {key:<36}{c:>7,}  {c / n:>6.1%}{flag}")
            grand[key] += c

        if targets:
            print(f"\n  THE TARGETS THOSE ROWS SIT ON (the record that would "
                  f"have to change):")
            for t, c in targets.most_common(args.samples):
                print(f"    {c:>6,}  {t}")

        for key in sorted(mismatched):
            if key == "panel":
                continue
            rowset = mismatched[key]
            print(f"\n  [names {key}, filed as {v['name']}] "
                  f"{len(rowset):,} row(s), up to {args.samples}:")
            seen = set()
            for desc, tgt, tname in rowset:
                if desc in seen:
                    continue
                seen.add(desc)
                if len(seen) > args.samples:
                    break
                print(f"      {desc[:92]}")
            print(f"      ({len(set(d for d, _, _ in rowset)):,} distinct "
                  f"description(s) in this bucket)")

        if mismatched.get("panel"):
            rowset = mismatched["panel"]
            print(f"\n  [names this virus AND others] {len(rowset):,} row(s) "
                  f"-- cross-reactivity panels, not misfiling:")
            for desc, _t, _tn in rowset[:args.samples]:
                print(f"      {desc[:92]}")
        print()

    print("=" * 78)
    print("HOW TO READ THIS")
    print("=" * 78)
    print("""  A row under "names a DIFFERENT virus" is not proof of an error. It has
  three readings and this script does not choose between them:

    cross-reactivity panel   one paper tests a compound against all four
                             serotypes and ChEMBL files the series against one
                             representative target.
    coarse target curation   the target really is "dengue virus protease" and
                             the tax_id is a stand-in.
    genuine misfiling        the measurement belongs to another serotype.

  Only the third is a data error. The target chembl_id above is what tells
  them apart, which is why it is printed rather than summarised.

  What makes this worth checking at all: M1 is the direct-acting channel, so a
  row filed under the wrong species puts an INHIBITS edge on a protein the
  compound was never measured against. And if DENV-1, DENV-3 and DENV-4
  protein evidence is sitting under DENV-2, then DENV-2's substrate advantage
  is partly theirs -- which bears on which flavivirus gets built first.

  This script measures and recommends nothing.""")
    if grand["names a DIFFERENT virus"]:
        print(f"\n  {grand['names a DIFFERENT virus']:,} row(s) across all "
              f"viruses name a different virus than they are filed under.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
