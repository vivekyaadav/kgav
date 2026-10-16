#!/usr/bin/env python3
"""What kind of literature is ChEMBL, and is it the kind this project needs?

THE QUESTION THIS ANSWERS. Phase 4 found that compounds measured on two
viral proteases are almost always measured in ONE paper: 47 of control C1's
49 pairs, and 10 of arm A's 12. The survey called that "who screened what
together" and left it there. Underneath it is a question about the SOURCE,
and it has a measurable answer.

The first evidence came from a one-line query. The 2021 PLoS ONE screen "A
drug repurposing screen identifies hepatitis C antivirals as inhibitors of
the SARS-CoV2 main protease" -- which is arm B's experiment, published -- is
not in chembl_37. Nor is almost anything else from that journal: 13
documents in the whole database, none from 2021, against a journal
publishing on the order of fifteen thousand papers a year.

WHY THAT MATTERS TO A REPURPOSING PROJECT. A medicinal-chemistry paper is a
congeneric series: forty analogues of one scaffold on one target. A
repurposing screen is the opposite shape -- one library across many targets
-- and it is the shape that answers "would this existing drug work on that
virus". If the chemistry source curates the first and not the second, then
the arm structure phase 4 measured is a property of the database rather than
of the biology, and so is some part of why degree beats every mechanistic
channel in docs/evaluation.md.

WHAT IS MEASURED HERE, and what is NOT. This script prints the journal
profile and the year profile. It does not classify journals into kinds:
"medicinal chemistry" is a judgement, and a judgement asserted in code reads
later as a measurement. The named-venue block below is a lookup of specific
titles, each one a place repurposing screens are published, and it reports
their counts beside the profile rather than rendering a verdict.

THE YEAR PROFILE IS THE POINT, and it is printed COMPLETE. The top fifteen
rows of a per-year table showed 2022-2025 and no 2020 or 2021, which would
say our coronavirus protease data post-dates the repurposing wave entirely.
That may be true and it may be an artifact of taking fifteen rows; a capped
list that decides a conclusion is the failure this repository has now fixed
in four separate reports. Every year with a row appears below.

    python scripts/diagnose_source_scope.py
    python scripts/diagnose_source_scope.py --db path/to/chembl_37.db

Read-only. Writes nothing. Recommends nothing.
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.schema import print_capped

# ONE READER, TWO REPORTS. The organism patterns, the role patterns and the
# polyprotein description rule all live in phase4_survey_chemistry, and a
# second copy here would be a second definition of which targets are viral
# proteases -- on the question of whether the SOURCE and the SURVEY agree.
# They cannot disagree if there is one of them.
from phase4_survey_chemistry import (  # noqa: E402
    FAMILIES,
    load_activities,
    load_targets,
    missing_columns,
    role_key,
    row_role,
)

# Venues where drug-repurposing screens appear. NOT a classification of the
# profile printed above them -- a lookup of named titles, so that "ChEMBL
# does not index this kind of journal" is a count rather than an impression.
# Each is here because screening campaigns and repurposing studies publish
# there; the list makes no claim to be complete.
REPURPOSING_VENUES = (
    "PLoS One", "Sci Rep", "Nat Commun", "PLoS Pathog", "eLife",
    "Antiviral Res", "ACS Infect Dis", "Antimicrob Agents Chemother",
    "SLAS Discov", "J Biomol Screen", "Nature", "Science", "Cell",
    "Proc Natl Acad Sci U S A", "Emerg Microbes Infect", "Viruses",
    "bioRxiv", "Front Pharmacol", "Br J Pharmacol",
)

# The window in which the COVID-19 repurposing wave was published. Declared
# as a constant so the share below cannot be read off a boundary chosen
# after seeing the numbers.
REPURPOSING_WINDOW = (2020, 2021)


def journal_profile(con: sqlite3.Connection) -> tuple[Counter, int]:
    """(journal -> documents, total documents with a journal)."""
    rows = con.execute("""SELECT journal, COUNT(*) FROM docs
                          WHERE journal IS NOT NULL AND journal != ''
                          GROUP BY journal""").fetchall()
    return Counter({j: n for j, n in rows}), sum(n for _, n in rows)


def venue_counts(con: sqlite3.Connection) -> list[tuple[str, int, int]]:
    """(venue, documents, documents in the repurposing window), per venue."""
    lo, hi = REPURPOSING_WINDOW
    out = []
    for v in REPURPOSING_VENUES:
        n, w = con.execute(
            "SELECT COUNT(*), SUM(year BETWEEN ? AND ?) FROM docs "
            "WHERE journal = ? COLLATE NOCASE", (lo, hi, v)).fetchone()
        out.append((v, n or 0, w or 0))
    return out


def protease_rows(acts: list[dict], targets: dict[int, dict]
                  ) -> dict[str, list[dict]]:
    """arm role -> the activity rows that resolved to it."""
    out: dict[str, list[dict]] = defaultdict(list)
    for a in acts:
        t = targets.get(a["tid"])
        if t is None:
            continue
        key = role_key(t, row_role(t, a["description"]))
        if key:
            out[key].append(a)
    return out


def by_year_and_journal(rows: list[dict], docs: dict[int, tuple]
                        ) -> tuple[Counter, Counter, int]:
    """(year -> rows, journal -> rows, rows with no document)."""
    years: Counter = Counter()
    journals: Counter = Counter()
    undated = 0
    for a in rows:
        d = docs.get(a["doc_id"])
        if d is None or d[0] is None:
            undated += 1
        else:
            years[d[0]] += 1
        if d is not None and d[1]:
            journals[d[1]] += 1
    return years, journals, undated


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top", type=int, default=15,
                    help="journals to list in the profile (the YEAR table is "
                         "never capped)")
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    args = ap.parse_args()

    if not args.db.exists():
        print(f"no ChEMBL database at {args.db}\nExiting 2: nothing was "
              f"measured.")
        return 2
    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    gone = missing_columns(con)
    if gone:
        print("this database does not have the columns this script reads:")
        print_capped(gone, 20, indent="    ", unit="more column(s)")
        print("\nExiting 2: nothing was measured.")
        return 2
    print(f"{args.db}\n")

    # ------------------------------------------------- what ChEMBL curates
    prof, total = journal_profile(con)
    print(f"JOURNAL PROFILE — {total:,} documents carrying a journal, "
          f"{len(prof):,} distinct titles.")
    print(f"  the {args.top} largest, with cumulative share:")
    run = 0
    for j, n in prof.most_common(args.top):
        run += n
        print(f"    {j[:44]:<46} {n:>7,}  {run / total:>6.1%}")
    print(f"    the other {len(prof) - args.top:,} titles supply "
          f"{total - run:,} ({1 - run / total:.1%})")
    print("\n  No judgement about what kind of journals these are is made "
          "here.\n  The names are printed; the reader classifies them.")

    lo, hi = REPURPOSING_WINDOW
    print(f"\n  named venues where repurposing screens are published, with "
          f"their\n  {lo}-{hi} share — a lookup of titles, not a "
          f"classification of the list above:")
    print(f"    {'venue':<30} {'docs':>8} {f'{lo}-{hi}':>10}")
    for v, n, w in sorted(venue_counts(con), key=lambda r: -r[1]):
        print(f"    {v:<30} {n:>8,} {w:>10,}")

    # --------------------------------- and what our own protease rows are
    targets = load_targets(con)
    by_tid = {t["tid"]: t for t in targets}
    if not targets:
        print("\nno target in any matched organism.\nExiting 2.")
        con.close()
        return 2
    acts = load_activities(con, [t["tid"] for t in targets])
    doc_ids = {a["doc_id"] for a in acts if a["doc_id"] is not None}
    docs: dict[int, tuple] = {}
    for chunk in [sorted(doc_ids)[i:i + 900]
                  for i in range(0, len(doc_ids), 900)]:
        q = ",".join("?" * len(chunk))
        for did, yr, jr in con.execute(
                f"SELECT doc_id, year, journal FROM docs WHERE doc_id IN ({q})",
                chunk):
            docs[did] = (yr, jr)
    con.close()

    roles = protease_rows(acts, by_tid)
    print(f"\n\nWHERE THIS PROJECT'S PROTEASE ROWS COME FROM — "
          f"{sum(len(v) for v in roles.values()):,} rows that resolved to an "
          f"arm role.")
    for key in sorted(roles):
        years, journals, undated = by_year_and_journal(roles[key], docs)
        inwin = sum(n for y, n in years.items() if lo <= y <= hi)
        tot = sum(years.values())
        print(f"\n  {key}  —  {len(roles[key]):,} rows, {undated:,} with no "
              f"dated document")
        print(f"    {lo}-{hi}, the repurposing window: {inwin:,} of {tot:,} "
              f"dated rows ({inwin / tot:.1%})" if tot else
              "    no dated rows")
        # EVERY YEAR, NOT THE TOP FEW. A truncated version of this table is
        # what suggested the coronavirus data post-dates 2021 entirely, and
        # a capped list that decides a conclusion is the failure this
        # repository has now fixed in four separate reports.
        print("    rows by year, complete:")
        for y in sorted(years):
            bar = "#" * max(1, round(40 * years[y] / max(years.values())))
            print(f"      {y}  {years[y]:>6,}  {bar}")
        print("    journals:")
        print_capped(journals.most_common(), 6,
                     lambda kv: f"{kv[0][:40]:<42} {kv[1]:>6,}",
                     indent="      ", unit="more journal(s)")

    print("\nRead-only. This script measures the SOURCE, not the biology, "
          "and\nrecommends nothing.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
