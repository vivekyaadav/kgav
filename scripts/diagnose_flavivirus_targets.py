#!/usr/bin/env python3
"""Which flavivirus proteins do ChEMBL's assay descriptions actually name?

THE GATE THIS MEASURES. chembl.NSP_PATTERNS and ingest_spine._nsp_number both
name CORONAVIRUS proteins, and flaviviruses carry NS1-NS5. Measured against
DENV-2's real chain labels, the two halves fail differently:

    NS5   label -> nsp12 and "NS5 RdRp" -> nsp12, so the edge lands on the
          right chain through a wrong name. It works BY ACCIDENT.
    NS3   label -> None, so the chain never enters viral_targets. A
          "NS2B-NS3 protease" description -> None and "NS3 helicase" ->
          nsp13, neither of which has an entry, so both fall through to
          UniProtKB:<polyprotein> -- one node standing in for thirteen
          proteins, which is the chloroquine-on-spike failure relocated.

So the fix is family-scoped naming on BOTH sides at once, and a half-fix
breaks the accident that makes NS5 work. Which patterns to write is a
question about the real description text, not about what the proteins are
called in a textbook. This script reads that text.

    python scripts/diagnose_flavivirus_targets.py
    python scripts/diagnose_flavivirus_targets.py --family Flaviviridae

Read-only. Writes nothing. Reuses chembl.query_activities, chembl.resolve_nsp
and chembl.assay_readout, so every count is what the ingest would see.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.chembl import assay_readout, query_activities, resolve_nsp

# The keys a flavivirus CHAIN can supply, i.e. what ingest_spine.chain_family
# returns for the thirteen Chain descriptions on Q32ZE1 and P29990 (UniProt,
# 2026-10-11 -- identical on both). Resolving a description to a key no chain
# carries leaves the row on the polyprotein exactly as an unresolved one does,
# so the two are counted apart below.
CHAIN_VOCABULARY = {"C", "prM", "pr", "M", "E", "NS1", "NS2A", "NS2B", "NS3",
                    "NS4A", "NS4B", "NS5"}
from kgav.virus_register import declared_taxa

# Candidate tokens, matched only to COUNT them. Nothing here is adopted as a
# pattern: the point is to find which ones carry rows, so the eventual
# patterns are derived from this table rather than from the list.
PROBES: list[tuple[str, re.Pattern]] = [
    ("NS1", re.compile(r"\bNS1\b|non[- ]structural protein 1\b", re.IGNORECASE)),
    ("NS2A", re.compile(r"\bNS2A\b", re.IGNORECASE)),
    ("NS2B", re.compile(r"\bNS2B\b", re.IGNORECASE)),
    ("NS2B-NS3", re.compile(r"NS2B[- /]*NS3|NS3[- /]*NS2B|CF40",
                            re.IGNORECASE)),
    # \bNS3\b CANNOT MATCH "NS3pro": there is no word boundary after the 3.
    # DENV-2's single largest assay is "Inhibitory activity against
    # CF40.NS3pro from Dengue virus type 2" -- CF40 is the 40-residue NS2B
    # cofactor construct -- and it landed in "named no known protein"
    # alongside a PubChem capping-enzyme campaign. Two descriptions, 1,025
    # rows, reported as unexplained because the probe was too strict. The
    # bucket was never unexplained; the probes could not see it.
    ("NS3", re.compile(r"\bNS3\b|NS3pro|NS3[- ]?protease", re.IGNORECASE)),
    ("CF40/NS3pro", re.compile(r"CF40|NS3pro", re.IGNORECASE)),
    ("capping enzyme", re.compile(r"capping enzyme|genomic capping",
                                  re.IGNORECASE)),
    ("PubChem HTS", re.compile(r"PUBCHEM_BIOASSAY", re.IGNORECASE)),
    ("NS4A", re.compile(r"\bNS4A\b", re.IGNORECASE)),
    ("NS4B", re.compile(r"\bNS4B\b", re.IGNORECASE)),
    ("NS5", re.compile(r"\bNS5\b", re.IGNORECASE)),
    ("protease", re.compile(r"protease|proteinase", re.IGNORECASE)),
    ("helicase", re.compile(r"helicase", re.IGNORECASE)),
    ("polymerase/RdRp", re.compile(r"rdrp|rna[- ]dependent rna pol|"
                                   r"rna[- ]directed rna pol|polymerase", re.IGNORECASE)),
    ("methyltransferase", re.compile(r"methyltransferase|mtase|2'-o-methyl", re.IGNORECASE)),
    ("envelope E", re.compile(r"envelope|\bprotein E\b", re.IGNORECASE)),
    ("capsid", re.compile(r"capsid", re.IGNORECASE)),
    ("entry/fusion", re.compile(r"\bentry\b|fusion", re.IGNORECASE)),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "config/viruses.yaml")
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data/raw/chembl/chembl_37"
                                   "/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--family", default="Flaviviridae")
    ap.add_argument("--samples", type=int, default=8,
                    help="verbatim descriptions printed per bucket")
    args = ap.parse_args()

    if not args.db.exists():
        sys.exit(f"no ChEMBL db at {args.db}")

    cfg = yaml.safe_load(args.config.read_text())
    default_family = (cfg.get("defaults") or {}).get("family")
    # BOTH SECTIONS, and `surveyed:` is the point rather than a convenience.
    #
    # This read cfg["viruses"] only. The moment DENV-2 and Zika were demoted
    # back to `surveyed:` it matched nothing and exited -- a script whose only
    # job is to measure a virus's substrate BEFORE deciding to build it,
    # unable to look at any virus that has not already been built. Measuring a
    # surveyed virus is the entire use case.
    wanted = []
    for section in ("viruses", "surveyed"):
        for v in (cfg.get(section) or []):
            if (v.get("family") or default_family) == args.family:
                wanted.append((v, section))
    if not wanted:
        sys.exit(f"no virus in either section has family {args.family!r}; "
                 f"pass --family")

    full_map = declared_taxa(cfg)
    # Say which section each one is in: a surveyed virus's numbers describe
    # what a build WOULD see, not what any release contains.
    print(f"{len(wanted)} {args.family} virus(es): "
          + ", ".join(f"{v['name']} ({sec})" for v, sec in wanted) + "\n")

    for v, section in wanted:
        taxa = [t for t, species in full_map.items()
                if species == str(v["taxon"])]
        rows = query_activities(args.db, sorted(taxa))

        # FILTERED THE WAY THE INGEST FILTERS, which the first version of
        # this script did not do and which made its headline number wrong.
        #
        # ingest_activities branches on target_type == "ORGANISM" FIRST: such
        # a row becomes HAS_ANTIVIRAL_ACTIVITY_AGAINST on the taxon and never
        # needs a chain. Only after that does the protein path run, and it
        # also skips a multicomponent target ("/" in target_name).
        #
        # assay_readout alone is not that filter. When no pattern matches it
        # falls through to "trust the code" and a B-coded row returns
        # "protein" -- so "Antiviral activity against Dengue virus 2", which
        # names no cell line, no CPE, no plaque and no titer, was counted as
        # a protein measurement. 1,070 of the 1,531 rows reported as DENV-2's
        # M1 substrate matched no flavivirus-protein probe at all, which is
        # what exposed it.
        protein_rows = []
        where_mix: Counter = Counter()
        target_mix: Counter = Counter()
        skipped: Counter = Counter()
        for r in rows:
            ttype = r.get("target_type") or "?"
            target_mix[ttype] += 1
            if ttype == "ORGANISM":
                skipped["organism_target"] += 1
                continue
            if "/" in (r.get("target_name") or ""):
                skipped["multicomponent_target"] += 1
                continue
            w, _ = assay_readout(r.get("description") or "",
                                 r.get("assay_type") or "")
            where_mix[w] += 1
            if w == "protein":
                protein_rows.append(r)
            else:
                skipped[f"readout_{w}"] += 1

        print("=" * 78)
        print(f"{v['name']}  [{section}]  (taxa {', '.join(sorted(taxa))})")
        print("=" * 78)
        print(f"  {len(rows):,} activity rows")
        print(f"  target_type: {dict(target_mix)}")
        print(f"  reaching the protein path: readout {dict(where_mix)}")
        for k, n in skipped.most_common():
            print(f"    skipped {k:<24}{n:>7,}")
        print(f"\n  {len(protein_rows):,} rows reach the chain lookup -- M1's "
              f"real substrate\n")

        if not protein_rows:
            print("  nothing to resolve\n")
            continue

        # ---- BEFORE and AFTER, on the same rows, measured not asserted
        #
        # The family-aware table was argued for from 16 quoted descriptions in
        # tests/test_flavivirus_naming.py. That is enough to show the patterns
        # read real text and not enough to say what they do to the whole
        # substrate, so both tables are run over every row here and the
        # difference is printed. CHAIN_VOCABULARY is the set of keys a chain
        # can actually supply, read from UniProt, so "resolves" and "joins"
        # stay separate: a key no chain carries is still a polyprotein row.
        verdicts: Counter = Counter()
        before: Counter = Counter()
        by_verdict: dict[str, list[str]] = defaultdict(list)
        joined = unjoinable = 0
        for r in protein_rows:
            desc = r.get("description") or ""
            nsp, domain = resolve_nsp(desc, args.family)
            key = f"{nsp}{'/' + domain if domain else ''}" if nsp else "UNRESOLVED"
            verdicts[key] += 1
            by_verdict[key].append(desc)

            old_nsp, old_dom = resolve_nsp(desc)          # coronavirus table
            before[f"{old_nsp}{'/' + old_dom if old_dom else ''}"
                   if old_nsp else "UNRESOLVED"] += 1

            if nsp and nsp in CHAIN_VOCABULARY:
                joined += 1
            else:
                unjoinable += 1

        n = len(protein_rows)
        print(f"  resolve_nsp BEFORE (coronavirus table on flavivirus rows):")
        for key, c in before.most_common():
            note = ("  -> falls through to the polyprotein" if key == "UNRESOLVED"
                    else "  -> a CORONAVIRUS family name")
            print(f"    {key:<22}{c:>7,}{note}")

        print(f"\n  resolve_nsp NOW (--family {args.family}):")
        for key, c in verdicts.most_common():
            head = key.split("/")[0]
            if key == "UNRESOLVED":
                note = "  -> still falls through to the polyprotein"
            elif head in CHAIN_VOCABULARY:
                note = "  -> a chain the spine emits"
            else:
                note = "  -> NO CHAIN CARRIES THIS KEY, still the polyprotein"
            print(f"    {key:<22}{c:>7,}{note}")

        was = n - before["UNRESOLVED"]
        print(f"\n  WOULD THE JOIN LAND? (of {n:,} rows reaching the lookup)")
        print(f"    resolved to a real chain   {joined:>7,}  "
              f"{joined / n:>6.1%}")
        print(f"    still the polyprotein      {unjoinable:>7,}  "
              f"{unjoinable / n:>6.1%}")
        print(f"    before this change         {was:>7,}  {was / n:>6.1%}"
              f"   (and {before['UNRESOLVED']:,} unresolved)")
        print("    A key no chain carries is not a join. This counts only keys"
              "\n    present in the UniProt chain vocabulary.")

        # ---- which flavivirus protein the text actually names
        print("\n  WHICH PROTEIN THE DESCRIPTION NAMES (probe counts, "
              "overlapping):")
        probe_hits: Counter = Counter()
        unnamed: list[str] = []
        for r in protein_rows:
            d = r.get("description") or ""
            hit = False
            for label, pat in PROBES:
                if pat.search(d):
                    probe_hits[label] += 1
                    hit = True
            if not hit:
                unnamed.append(d)
        named_none = len(unnamed)
        for label, _ in PROBES:
            n = probe_hits.get(label, 0)
            if n:
                print(f"    {label:<20}{n:>7,}")
        print(f"    {'(no probe matched)':<20}{named_none:>7,}")

        # CONCENTRATION OF THE SUBSTRATE ITSELF, which is the statistic the
        # promotion argument needed and did not have.
        #
        # diagnose_flavivirus_join reports a top-PMID share over ALL rows,
        # and DENV-2's 1% across 135 publications is what argued for building
        # it: a dispersed literature cannot reproduce the SARS-CoV-2 failure
        # where 94.9% of compounds shared a PMID between their INHIBITS edge
        # and their own label.
        #
        # But that figure covers 8,748 rows, of which 7,227 are cellular or
        # organism-level and never reach a chain. The subset M1 actually
        # traverses is 1,483, and 1,013 of them are a single PubChem capping
        # enzyme campaign. A statistic computed over the whole file says
        # nothing about the concentration of the part a channel uses.
        #
        # Rows with NO pubmed_id matter separately: --publication-disjoint
        # withholds an edge that shares a PMID with its own label, so a row
        # carrying no PMID is never withheld by that control at all.
        pmids_p: Counter = Counter()
        assays_p: Counter = Counter()
        no_pmid = 0
        assay_desc: dict[str, str] = {}
        for r in protein_rows:
            pm = r.get("pubmed_id")
            if pm:
                pmids_p[str(pm)] += 1
            else:
                no_pmid += 1
            a = r.get("assay") or "?"
            assays_p[a] += 1
            assay_desc.setdefault(a, r.get("description") or "")

        n = len(protein_rows)
        print(f"\n  CONCENTRATION OF M1's SUBSTRATE ({n:,} protein rows):")
        print(f"    {len(pmids_p):,} distinct PMID(s)")
        print(f"    {no_pmid:,} rows carry NO pubmed_id "
              f"({no_pmid / n:.0%}) -- --publication-disjoint cannot "
              f"withhold these")
        if pmids_p:
            pm, c = pmids_p.most_common(1)[0]
            print(f"    top PMID  PMID:{pm}  {c:,} rows ({c / n:.0%})")
        if assays_p:
            a, c = assays_p.most_common(1)[0]
            print(f"    top assay {a}  {c:,} rows ({c / n:.0%})")
            print(f"              {assay_desc[a][:66]}")
        print("    Compare this against the all-rows top-PMID share in "
              "diagnose_flavivirus_join:\n    they answer different "
              "questions and only this one is about M1.")

        # THE BUCKET THAT MATTERS MOST HAD NO SAMPLES. The verbatim section
        # below is keyed on resolve_nsp's verdict, so the rows that match no
        # protein probe at all were counted and never shown -- 1,025 of
        # DENV-2's 1,483, which is 69% of the substrate with nothing said
        # about it. A diagnostic that reports a large unexplained bucket and
        # then does not print it is asking to be guessed at, and I guessed.
        if unnamed:
            print(f"\n  NAMING NO KNOWN PROTEIN -- {named_none:,} rows, the "
                  f"bucket to explain first:")
            # COUNT DISTINCT OVER EVERYTHING, DISPLAY A SAMPLE. The first
            # version counted inside the display loop, which breaks at the
            # sample cap -- so it would have printed "16 distinct" whenever
            # there were 16 or more, and the one number that matters here is
            # how CONCENTRATED the bucket is. DENV-2's 1,025 rows turned out
            # to be 2 descriptions, which is the whole finding.
            #
            # Dedup on the full string, not d[:100]: two descriptions
            # differing past character 100 are different assays. Printing is
            # truncated; counting is not.
            distinct = Counter(unnamed)
            for d, n in distinct.most_common(args.samples * 2):
                print(f"      {n:>6,}  {d[:96]}")
            print(f"      -- {len(distinct):,} distinct description(s) across "
                  f"{named_none:,} rows")
            if len(distinct) > args.samples * 2:
                print(f"      ({args.samples * 2} shown)")

        # ---- verbatim text, because patterns are written against text
        print(f"\n  VERBATIM BY resolve_nsp VERDICT, up to {args.samples} "
              f"per bucket:")
        for key, _ in verdicts.most_common(4):
            print(f"\n    [{key}]")
            seen = set()
            shown = 0
            for d in by_verdict[key]:
                k = d[:90]
                if k in seen:
                    continue
                seen.add(k)
                print(f"      {k}")
                shown += 1
                if shown >= args.samples:
                    break
        print()

    print("=" * 78)
    print("HOW TO READ THIS")
    print("=" * 78)
    print("""  A probe with many rows and an UNRESOLVED verdict is a pattern worth
  writing. A probe with many rows resolving to an nsp name is a row that
  currently lands somewhere by accident -- check whether the chain it reaches
  is the right one before changing anything, because NS5 works today only
  because both sides are wrong in the same direction.

  The chain side must move with the description side. ingest_spine assigns
  protein_family = f"nsp{N}" and nothing else, so until it can emit NS1..NS5
  a flavivirus-aware resolve_nsp would return names that
  build_viral_target_index has no key for, and every row would fall through
  to the polyprotein -- which is what happens to NS3 today.

  This script measures and recommends nothing.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
