#!/usr/bin/env python3
"""What the VirHostNet rows counted as `unresolved_id` actually are.

    python scripts/diagnose_unresolved_ids.py
    python scripts/diagnose_unresolved_ids.py --samples 10

WHY. ingest_vhppi.py reports `unresolved_id 513` (4.4% of rows on v0.1) and
stops there. That single number is consistent with two opposite situations:

    a HOST protein the build does not carry   -> interactome that is missing,
                                                 and M3/M4/M5 are thinner for it
    a VIRAL protein the build does not carry  -> a virus outside the register,
                                                 correctly dropped

Both print as `unresolved_id`. Which one it is decides whether there is
anything to fix, and no existing script says: diagnose_flavivirus_join.py
reports the taxon histogram of the FILES and explicitly refuses to claim
which gate a row dies at ("Only the ingest's own dropped counters say which
gate a row dies at; run it and read them").

HOW IT STAYS HONEST. The gates run in a fixed order and a row is counted at
the FIRST one it fails, so re-testing resolution on every row would count
rows that really died at `unparseable`. This reads the order from
vhppi.classify_row -- the same function ingest_mitab consumes -- rather than
repeating it, and test_the_classifier_and_the_ingest_count_the_same_gates
pins the two together.

WHAT IT DOES NOT DECIDE. Nothing here resolves an identifier or promotes a
virus. An accession absent from the loaded releases may be obsolete,
secondary, from an unbuilt proteome, or from a host species this project
does not model, and this script cannot tell those apart -- it reports the
taxid the FILE gives for that identifier, verbatim, so the ids can be looked
up. Deciding what to do about a taxon is a curation question.
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.aliases import AliasIndex
from kgav.virus_register import declared_taxa
from kgav.vhppi import classify_row


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=root / "config" / "viruses.yaml")
    ap.add_argument("--vhppi-dir", type=Path, default=root / "data" / "raw" / "vhppi")
    ap.add_argument("--spine", type=Path, default=root / "data" / "releases" / "v0.1-spine")
    ap.add_argument("--host", type=Path, default=root / "data" / "releases" / "v0.1-host")
    ap.add_argument("--samples", type=int, default=6,
                    help="verbatim identifiers to print per taxon")
    args = ap.parse_args()

    for p in (args.spine, args.host):
        if not (p / "nodes.jsonl").exists():
            print(f"no release at {p}", file=sys.stderr)
            return 2
    files = sorted(args.vhppi_dir.glob("vhn_*.tab27"))
    if not files:
        print(f"no vhn_*.tab27 files in {args.vhppi_dir}", file=sys.stderr)
        return 2

    index = AliasIndex.from_releases(args.spine, args.host)
    taxon_map = declared_taxa(yaml.safe_load(args.config.read_text()))
    print(f"index: {len(index):,} nodes, {len(index.alias):,} aliases")
    print(f"register: {len(taxon_map)} viral taxids -> "
          f"{len(set(taxon_map.values()))} species")
    print(f"parsing {len(files)} MITAB files\n")

    gates: Counter = Counter()
    rows = 0
    # Per unresolved identifier, the taxid the file gives it. Counted by
    # IDENTIFIER as well as by row: one missing hub protein can account for
    # hundreds of rows, and the row count alone would read as hundreds of
    # separate gaps.
    per_tax_rows: Counter = Counter()
    per_tax_ids: dict[str, set[str]] = defaultdict(set)
    id_rows: Counter = Counter()
    both_sides = 0
    chain_ids = 0
    unparseable = 0
    unparseable_tax: Counter = Counter()
    unparseable_ex: dict[str, list[str]] = defaultdict(list)

    for path in files:
        for line in path.read_text(errors="replace").splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            rows += 1
            v = classify_row(line, index, taxon_map)
            gates[v.gate or "EMITTED (before dedupe)"] += 1
            if v.gate == "unparseable":
                unparseable += 1
                f = line.rstrip("\n").split("\t")
                ta = (v.row or {}).get("taxid_a") or "?"
                tb = (v.row or {}).get("taxid_b") or "?"
                pair = f"{ta} / {tb}"
                unparseable_tax[pair] += 1
                if len(unparseable_ex[pair]) < 2 and len(f) >= 2:
                    unparseable_ex[pair].append(f"{f[0]}  {f[1]}")
                continue
            if v.gate != "unresolved_id":
                continue
            if len(v.unresolved) == 2:
                both_sides += 1
            for raw, tax, _side in v.unresolved:
                tax = tax or "(no taxid)"
                per_tax_rows[tax] += 1
                per_tax_ids[tax].add(raw)
                id_rows[raw] += 1
                if "PRO_" in raw:
                    chain_ids += 1

    print(f"{rows:,} rows\n")
    print("gate reached (first failure wins):")
    for gate, n in gates.most_common():
        print(f"  {gate:<28} {n:>8,}")

    if unparseable:
        # THE GATE BEFORE THIS ONE, broken down the same way, because the
        # register asserts a row count at `unresolved_id` that this script
        # can check -- and the two gates are adjacent, so a claim naming the
        # wrong one is exactly the error diagnose_flavivirus_join.py was
        # corrected for.
        #
        # The raw identifier columns are split out here rather than in
        # parse_row: a row counted `unparseable` is one parse_identifier
        # returned None for, so the verbatim column is the only thing left
        # to look at. No gate logic is repeated -- classify_row still
        # decided which rows these are.
        print(f"\n{'=' * 74}\nUNPARSEABLE: {unparseable:,} rows "
              f"({unparseable / rows:.1%})\n{'=' * 74}")
        print("  parse_identifier recognises a UniProt accession or a PRO_ "
              "chain and nothing\n  else, so a row in another namespace dies "
              "HERE, one gate before unresolved_id.")
        print(f"\n    {'taxid pair':<22} {'rows':>7}  example identifier columns")
        for pair, n in unparseable_tax.most_common(12):
            ex = "  |  ".join(s[:34] for s in unparseable_ex[pair][:2])
            print(f"    {pair:<22} {n:>7,}  {ex}")

    unresolved_rows = gates.get("unresolved_id", 0)
    if not unresolved_rows:
        print("\nNo row died at unresolved_id. Nothing to report.")
        return 0

    print(f"\n{'=' * 74}\nUNRESOLVED IDENTIFIERS: {unresolved_rows:,} rows "
          f"({unresolved_rows / rows:.1%}), {len(id_rows):,} distinct ids")
    print(f"{'=' * 74}")
    print(f"  {both_sides:,} rows had BOTH endpoints unresolved")
    print(f"  {chain_ids:,} id-slots were PRO_ chain references")

    # SLOTS, NOT ROWS. A row with both endpoints unresolved contributes an
    # identifier to two taxids, so this column sums to more than the row
    # count above. Calling it `rows` would make the two disagree by exactly
    # the both-sides count and look like an arithmetic error.
    print("\n  by the taxid the FILE gives the unresolved identifier "
          f"({sum(per_tax_rows.values()):,} id-slots across "
          f"{unresolved_rows:,} rows; a row with both endpoints unresolved "
          "appears under two taxids):")
    print(f"    {'taxid':<14} {'slots':>8} {'ids':>7}  in register   examples")
    for tax, n in per_tax_rows.most_common():
        ids = sorted(per_tax_ids[tax])
        known = "species" if tax in set(taxon_map.values()) else \
                "isolate" if tax in taxon_map else "-"
        ex = ", ".join(i.replace("UniProtKB:", "") for i in ids[:args.samples])
        if len(ids) > args.samples:
            ex += f", +{len(ids) - args.samples:,} more"
        print(f"    {tax:<14} {n:>8,} {len(ids):>7,}  {known:<11}  {ex}")

    print(f"\n  the {min(args.samples * 2, len(id_rows))} identifiers costing "
          f"the most rows:")
    for raw, n in id_rows.most_common(args.samples * 2):
        print(f"    {raw:<34} {n:>6,} row" + ("" if n == 1 else "s"))

    print(f"""
{'=' * 74}
HOW TO READ THIS
{'=' * 74}
  A taxid marked `species` or `isolate` is ALREADY in the register's map,
  so the loss is not a missing declaration: the accession itself is absent
  from the spine and host releases. Check whether it belongs to a proteome
  this build loads before concluding anything.

  A taxid marked `-` is outside the register. If it is a host species this
  project does not model, the drop is correct. If it is a virus you build
  under a different taxid, it may be an isolate_taxa entry -- but look the
  id up first; the taxid column above is verbatim from the file.

  A single identifier costing many rows is one hub protein, not many gaps.
  Resolving it would return all of those rows at once.

WHAT THIS DOES NOT DECIDE
  Nothing here is a reason on its own to add a taxon or an alias. This
  script reports what the files say; whether an accession SHOULD resolve is
  a question about which proteomes the build loads, and promoting a taxon
  is a curation decision that config/viruses.yaml records by hand.
""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
