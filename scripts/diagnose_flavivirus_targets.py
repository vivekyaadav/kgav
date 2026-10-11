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
from kgav.virus_register import declared_taxa

# Candidate tokens, matched only to COUNT them. Nothing here is adopted as a
# pattern: the point is to find which ones carry rows, so the eventual
# patterns are derived from this table rather than from the list.
PROBES: list[tuple[str, re.Pattern]] = [
    ("NS1", re.compile(r"\bNS1\b|non[- ]structural protein 1\b", re.IGNORECASE)),
    ("NS2A", re.compile(r"\bNS2A\b", re.IGNORECASE)),
    ("NS2B", re.compile(r"\bNS2B\b", re.IGNORECASE)),
    ("NS2B-NS3", re.compile(r"NS2B[- /]*NS3|NS3[- /]*NS2B", re.IGNORECASE)),
    ("NS3", re.compile(r"\bNS3\b", re.IGNORECASE)),
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
    wanted = [v for v in cfg["viruses"]
              if (v.get("family") or default_family) == args.family]
    if not wanted:
        sys.exit(f"no built virus has family {args.family!r}; "
                 f"promote one first or pass --family")

    full_map = declared_taxa(cfg)
    print(f"{len(wanted)} built {args.family} virus(es): "
          f"{', '.join(v['name'] for v in wanted)}\n")

    for v in wanted:
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
        print(f"{v['name']}  (taxa {', '.join(sorted(taxa))})")
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

        # ---- what resolve_nsp does with them today
        verdicts: Counter = Counter()
        by_verdict: dict[str, list[str]] = defaultdict(list)
        for r in protein_rows:
            nsp, domain = resolve_nsp(r.get("description") or "")
            key = f"{nsp}{'/' + domain if domain else ''}" if nsp else "UNRESOLVED"
            verdicts[key] += 1
            by_verdict[key].append(r.get("description") or "")

        print("  resolve_nsp TODAY (coronavirus families on flavivirus rows):")
        for key, n in verdicts.most_common():
            note = ""
            if key == "UNRESOLVED":
                note = "  -> falls through to the polyprotein"
            elif key.startswith(("nsp", "NSP")):
                note = "  -> a CORONAVIRUS family name"
            print(f"    {key:<14}{n:>7,}{note}")

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

        # THE BUCKET THAT MATTERS MOST HAD NO SAMPLES. The verbatim section
        # below is keyed on resolve_nsp's verdict, so the rows that match no
        # protein probe at all were counted and never shown -- 1,025 of
        # DENV-2's 1,483, which is 69% of the substrate with nothing said
        # about it. A diagnostic that reports a large unexplained bucket and
        # then does not print it is asking to be guessed at, and I guessed.
        if unnamed:
            print(f"\n  NAMING NO KNOWN PROTEIN -- {named_none:,} rows, the "
                  f"bucket to explain first:")
            seen_u: set[str] = set()
            shown_u = 0
            for d in unnamed:
                k = d[:100]
                if k in seen_u:
                    continue
                seen_u.add(k)
                print(f"      {k}")
                shown_u += 1
                if shown_u >= args.samples * 2:
                    break
            print(f"      ({len(seen_u):,} distinct of {named_none:,} shown "
                  f"up to {args.samples * 2})")

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
