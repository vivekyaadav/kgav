#!/usr/bin/env python3
"""Would a flavivirus build produce a graph with labels, or an empty shell?

config/viruses.yaml `surveyed:` holds nine flaviviruses with pinned strain
taxa, verified chain counts and proteome ids. The UniProt side is DONE: the
MERS species/strain trap is pre-solved for all nine. What is NOT known is
whether the other four layers have anything to join to those taxa.

THE HYPOTHESIS THIS EXISTS TO TEST. The register solved the strain trap for
UniProt CHAINS. ChEMBL is a different database with its own taxonomy
decisions, so the strain it files dengue ACTIVITIES under need not be the
strain UniProt curates chains under. Querying only the register's taxa would
return zero rows and look like "no dengue data in ChEMBL" -- the same silent
zero, one database over. So this script queries the register's taxa AND
sweeps target_dictionary for every tax_id whose organism name looks
flaviviral, and reports where the rows actually are.

WHAT DECIDES THE BUILD. After the coronavirus work, two numbers matter more
than total volume:

  publications   SARS-CoV-2's M1 collapsed because 94.9% of compounds shared
                 a PMID between their INHIBITS edge and their label -- a
                 2020 single-burst literature. A virus whose activities are
                 spread over decades and many papers cannot fail that way.
                 Reported as the share held by the single largest PMID.

  host screens   M2/M4/M5 need BioGRID-ORCS dependency factors. SARS-CoV-2
                 got dozens of CRISPR screens in 2020; flaviviruses did not.
                 If ORCS has nothing for these taxa, four of seven channels
                 are dead on arrival and only M1/M6/M7 are testable.

    python scripts/diagnose_flavivirus_join.py
    python scripts/diagnose_flavivirus_join.py --cutoff 2016   # Zika outbreak

Read-only. Writes nothing. Degrades per-source: a missing ChEMBL db or raw
directory is reported and the rest still runs.
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.chembl import assay_readout, query_activities

# Names ChEMBL may file these under. Matched case-insensitively against
# target_dictionary.organism to FIND tax_ids, never to resolve one: every
# match is reported with its id for curation, not silently adopted.
FLAVI_NAME_HINTS = ("dengue", "zika", "west nile", "japanese encephalitis",
                    "yellow fever", "tick-borne encephalitis",
                    "tick borne encephalitis", "flavivirus")


def _surveyed(config: Path) -> list[dict]:
    cfg = yaml.safe_load(config.read_text())
    return cfg.get("surveyed") or []


def _register_taxa(entry: dict) -> list[str]:
    """Species taxon plus every declared isolate, as strings."""
    taxa = [str(entry["taxon"])]
    taxa += [str(x) for x in (entry.get("isolate_taxa") or [])]
    cs = (entry.get("chain_source") or {}).get("taxon")
    if cs and str(cs) not in taxa:
        taxa.append(str(cs))
    return taxa


def _sweep_organisms(db: Path) -> dict[str, tuple[str, int]]:
    """tax_id -> (organism name, target rows) for anything flavivirus-named.

    Finds where ChEMBL actually put the data, which the register cannot know:
    it records UniProt's curation decisions, not ChEMBL's.
    """
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    like = " OR ".join("LOWER(t.organism) LIKE ?" for _ in FLAVI_NAME_HINTS)
    rows = con.execute(f"""
        SELECT t.tax_id, t.organism, COUNT(DISTINCT t.tid) AS targets
        FROM target_dictionary t
        WHERE t.tax_id IS NOT NULL AND ({like})
        GROUP BY t.tax_id, t.organism
    """, [f"%{h}%" for h in FLAVI_NAME_HINTS]).fetchall()
    con.close()
    return {str(r["tax_id"]): (r["organism"], r["targets"]) for r in rows}


def _activity_report(rows: list[dict], cutoff: int) -> dict:
    """The numbers that decide whether a channel is testable."""
    compounds, pmids, accessions = set(), Counter(), set()
    years: Counter = Counter()
    where: Counter = Counter()
    per_compound_pmids: dict[str, set] = defaultdict(set)

    for r in rows:
        cid = r.get("inchikey") or r.get("compound")
        if cid:
            compounds.add(cid)
        if r.get("accession"):
            accessions.add(r["accession"])
        if r.get("pubmed_id"):
            pmids[str(r["pubmed_id"])] += 1
            if cid:
                per_compound_pmids[cid].add(str(r["pubmed_id"]))
        if r.get("year"):
            years[int(r["year"])] += 1
        w, _at = assay_readout(r.get("description") or "",
                              r.get("assay_type") or "")
        where[w] += 1

    dated = sum(years.values())
    top_share = (pmids.most_common(1)[0][1] / dated) if (pmids and dated) else None
    return {
        "rows": len(rows),
        "compounds": len(compounds),
        "accessions": len(accessions),
        "publications": len(pmids),
        "top_pmid": pmids.most_common(1)[0] if pmids else None,
        "top_pmid_share": top_share,
        "dated_rows": dated,
        "undated_rows": len(rows) - dated,
        "pre_cutoff": sum(n for y, n in years.items() if y <= cutoff),
        "post_cutoff": sum(n for y, n in years.items() if y > cutoff),
        "year_span": (min(years), max(years)) if years else None,
        "readout": dict(where),
        "multi_pmid_compounds": sum(1 for v in per_compound_pmids.values()
                                    if len(v) > 1),
    }


def _count_raw(directory: Path, taxa: set[str], label: str) -> dict:
    """How many lines in a raw directory mention any of these taxa.

    Deliberately a TEXT scan and not a parse: this answers "is there anything
    here at all", and a parser tuned to the coronavirus files would report
    zero for a flavivirus file whose columns differ -- the silent zero again.
    """
    if not directory.exists():
        return {"status": f"no {label} directory at {directory}"}
    hits: Counter = Counter()
    files = 0
    for p in sorted(directory.rglob("*")):
        if not p.is_file() or p.suffix.lower() not in (
                ".txt", ".tsv", ".csv", ".mitab", ".psi", ".json"):
            continue
        files += 1
        try:
            text = p.read_text(errors="replace")
        except OSError:
            continue
        for tax in taxa:
            n = text.count(tax)
            if n:
                hits[tax] += n
    return {"status": "ok", "files_scanned": files, "taxon_mentions": dict(hits)}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=root / "config/viruses.yaml")
    ap.add_argument("--db", type=Path,
                    default=root / "data/raw/chembl/chembl_37"
                                   "/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--orcs-dir", type=Path, default=root / "data/raw/orcs")
    ap.add_argument("--vhppi-dir", type=Path, default=root / "data/raw/vhppi")
    ap.add_argument("--cutoff", type=int, default=2016,
                    help="temporal split year; 2016 is the Zika outbreak")
    ap.add_argument("--min-compounds", type=int, default=10,
                    help="matches calibration.MIN_REACHED_POSITIVES: below "
                         "this a channel cannot be evaluated at all")
    args = ap.parse_args()

    surveyed = _surveyed(args.config)
    if not surveyed:
        sys.exit(f"no `surveyed:` entries in {args.config}")

    print(f"{len(surveyed)} surveyed flaviviruses in {args.config}")
    print(f"temporal cutoff {args.cutoff} | evaluability floor "
          f"{args.min_compounds} compounds\n")

    # ---------------------------------------------------------------- ChEMBL
    if not args.db.exists():
        print(f"!! no ChEMBL db at {args.db}")
        print("   the activity half of this report cannot run; "
              "ORCS/VirHostNet below still will\n")
        sweep, per_virus = {}, {}
    else:
        sweep = _sweep_organisms(args.db)
        print("=" * 78)
        print("WHERE ChEMBL ACTUALLY FILES FLAVIVIRUS TARGETS")
        print("=" * 78)
        print("The register pins where UniProt curates CHAINS. ChEMBL makes "
              "its own\ntaxonomy decisions, so a tax_id here that the "
              "register does not list is\nactivity data a build would miss "
              "entirely.\n")
        registered = {t for e in surveyed for t in _register_taxa(e)}
        print(f"  {'tax_id':<10} {'targets':>8}  {'in register':<12} organism")
        for tax, (organism, targets) in sorted(
                sweep.items(), key=lambda kv: -kv[1][1]):
            mark = "yes" if tax in registered else "NOT LISTED"
            print(f"  {tax:<10} {targets:>8}  {mark:<12} {organism}")
        unlisted = {t for t in sweep if t not in registered}
        if unlisted:
            print(f"\n  {len(unlisted)} tax_id(s) carry ChEMBL targets and are "
                  f"not in the register.\n  Curation decision, not a lookup: "
                  f"add them to isolate_taxa only after\n  checking the "
                  f"organism name is the virus you mean.")
        print()

        per_virus = {}
        for e in surveyed:
            taxa = _register_taxa(e)
            rows = query_activities(args.db, taxa)
            per_virus[e["name"]] = (taxa, _activity_report(rows, args.cutoff))

    # ------------------------------------------------------------- the table
    if per_virus:
        print("=" * 78)
        print("CHEMBL ACTIVITY PER SURVEYED VIRUS  (register taxa only)")
        print("=" * 78)
        print(f"  {'virus':<14}{'rows':>7}{'cmpds':>7}{'accs':>6}{'pubs':>6}"
              f"{'top%':>7}{'pre':>7}{'post':>7}  verdict")
        for name, (taxa, r) in per_virus.items():
            share = f"{r['top_pmid_share']:.0%}" if r["top_pmid_share"] else "-"
            verdict = ("NO DATA" if r["compounds"] == 0 else
                       "below floor" if r["compounds"] < args.min_compounds
                       else "evaluable")
            print(f"  {name:<14}{r['rows']:>7,}{r['compounds']:>7,}"
                  f"{r['accessions']:>6}{r['publications']:>6}{share:>7}"
                  f"{r['pre_cutoff']:>7,}{r['post_cutoff']:>7,}  {verdict}")
        print("\n  top% is the share of DATED rows held by the single largest "
              "PMID.\n  SARS-CoV-2's confound was 94.9% of compounds sharing "
              "a PMID with their\n  own label; a low top% here means that "
              "failure mode is structurally\n  unavailable.")

        print("\n" + "=" * 78)
        print("READOUT MIX  (assay_readout, the same classifier the ingest uses)")
        print("=" * 78)
        print(f"  {'virus':<14}{'cells':>8}{'protein':>9}{'unclear':>9}"
              f"{'span':>14}{'multi-PMID':>12}")
        for name, (_taxa, r) in per_virus.items():
            w = r["readout"]
            span = (f"{r['year_span'][0]}-{r['year_span'][1]}"
                    if r["year_span"] else "-")
            print(f"  {name:<14}{w.get('cells', 0):>8,}"
                  f"{w.get('protein', 0):>9,}{w.get('unclear', 0):>9,}"
                  f"{span:>14}{r['multi_pmid_compounds']:>12,}")
        print("\n  protein rows are M1's substrate; cells rows become "
              "organism-level\n  activity. A virus with no protein rows has "
              "no direct-acting channel.")

    # -------------------------------------------------------- host-side data
    all_taxa = {t for e in surveyed for t in _register_taxa(e)}
    print("\n" + "=" * 78)
    print("HOST-SIDE LAYERS  (M2/M4/M5 depend entirely on these)")
    print("=" * 78)
    for label, directory in (("ORCS", args.orcs_dir),
                             ("VirHostNet", args.vhppi_dir)):
        res = _count_raw(directory, all_taxa, label)
        print(f"\n  {label}: {res['status']}")
        if res.get("status") != "ok":
            continue
        print(f"    {res['files_scanned']} file(s) scanned")
        if not res["taxon_mentions"]:
            print("    NO mention of any surveyed taxon -- M2, M4 and M5 "
                  "would have no\n    host route for these viruses, leaving "
                  "M1, M6 and M7 as the only\n    testable channels.")
        for tax, n in sorted(res["taxon_mentions"].items(),
                             key=lambda kv: -kv[1]):
            print(f"    {tax:<10} {n:>8,} line mentions")
    print("\n  A text scan, not a parse: it answers 'is anything here at all'. "
          "A parser\n  tuned to the coronavirus files would report zero for a "
          "file whose columns\n  differ, which is the silent zero this whole "
          "register exists to prevent.")

    print("\n" + "=" * 78)
    print("WHAT THIS DOES NOT DECIDE")
    print("=" * 78)
    print("""  Nothing here promotes a virus. `surveyed:` graduates by being moved into
  `viruses:` by hand, and test_surveyed_viruses_are_not_built enforces that
  until it is. This report is the evidence for that decision, and the
  decision is yours: a tax_id ChEMBL uses that the register does not list is
  a CURATION question about which strain is the virus you mean, which is not
  something a script should answer by itself.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
