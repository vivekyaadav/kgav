#!/usr/bin/env python3
"""Rewrite every document table that a test already derives.

FIX-LIST P10, THE PART THAT WAS LEFT. Five tables in docs/ carry a
machine-readable marker and a test that derives the table the results file
implies, cell by cell -- and nothing could WRITE that derivation back. So
after every rebuild the figures were updated by hand against the test's
failure messages. On 2026-11-11 that was 111 failures in one afternoon, and
three of the figures I typed by hand in the same session were wrong: a
pooled compound count I invented (18,914 for 18,900), a delta computed from
two displayed four-decimal values instead of the stored ones, and an AUC
literal typed from memory.

EVERY FIGURE I DERIVED WAS RIGHT AND EVERY FIGURE I TYPED WAS NOT. This
makes the derivation the writer.

THE TEST REMAINS THE CHECK. This script does not contain a formatting rule
or a results-file path: it imports the test module, calls its `derive()`,
and writes what that returns. So the test and the document cannot disagree
about a value -- and if the test's derivation is wrong, the test is still
the thing to fix, exactly as before.

ROW ORDER IS THE DOCUMENT'S, NOT THE DERIVATION'S. The tables are ordered
for reading -- by AUC, by taxon, by channel -- and the tests check coverage
and cells, never order. Rewriting in `derive()` order would reshuffle five
tables for no reason and bury the real changes in the diff.
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

#: Each test module that owns a marked table and derives it.
MODULES = ("test_evaluation_table", "test_temporal_table",
           "test_landscape_table", "test_auc_decomposition",
           "test_release_plan_figures", "test_sequence_table",
           "test_scope_table", "test_landscape")


def rewrite(doc: Path, marker: str, rows: dict[str, list[str]],
            keep: int | None = None) -> int:
    """Replace the marked table's body with `rows`. Returns cells changed.

    THE KEY MAY SPAN SEVERAL COLUMNS, and the first version of this assumed
    one. `test_landscape_table.derive()` keys rows as "taxa | DENV-2" --
    channel and virus -- so every lookup missed, every row was skipped, and
    the script reported "0 cell(s) would change" for a table whose cells
    were wrong in five places. A generator that silently rewrites nothing
    is the same shape as a check that cannot fail, which is the defect this
    repository keeps finding, so the key width is read off the derivation
    and an unmatched row is REFUSED rather than passed over.

    `keep` IS FOR A LEADING COLUMN THAT IS PROSE RATHER THAN DATA. The
    sequence table's second column is the protocol's description -- the
    document's own words, which no results file produces -- so its rows are
    keyed on the row number alone while two leading columns must survive
    untouched. A module declares that with `KEEP_COLUMNS`; without one, the
    key columns and the preserved columns are the same, which is the case
    for the other five tables.
    """
    text = doc.read_text()
    assert marker in text, f"{doc.name} has no marker {marker!r}"
    head, body = text.split(marker, 1)
    lines = body.splitlines(keepends=True)
    first = next(i for i, ln in enumerate(lines) if ln.lstrip().startswith("|"))
    last = first
    while last < len(lines) and lines[last].lstrip().startswith("|"):
        last += 1

    # how many leading columns make up a row's key, from the derivation
    nkey = 1 + max(k.count(" | ") for k in rows) if rows else 1
    nkeep = nkey if keep is None else keep
    assert nkeep >= nkey, (
        f"{doc.name}: KEEP_COLUMNS={keep} is narrower than the "
        f"{nkey}-column key the derivation uses")

    # PADDING IS COPIED ONLY FROM A TABLE THAT HAS SOME. The first version
    # padded every cell to the width the OLD value had, so a value that got
    # shorter left a trailing space: "| 9 hits = **3.61x**  |" for a table
    # whose rows were single-spaced, and
    # test_every_row_of_the_published_table_is_re_derived -- which matches a
    # row as a literal substring -- failed on whitespace while every figure
    # in it was right. So a column is padded only if the document already
    # aligns it: every data row the same width in that column.
    data = [[c.strip() for c in ln.strip().strip("|").split("|")]
            for ln in lines[first + 2:last]]
    aligned: dict[int, int] = {}
    for col in range(max((len(r) for r in data), default=0)):
        ws = {len(r[col]) for r in data if col < len(r)}
        if len(ws) == 1 and len(data) > 1:
            aligned[col] = ws.pop()

    changed, missing = 0, []
    out = []
    for i, ln in enumerate(lines[first:last], start=first):
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if i <= first + 1:
            out.append(ln)                      # header and separator
            continue
        # EMPHASIS IS NOT PART OF THE KEY OR THE VALUE. The temporal table
        # bolds its COMBINED row, and `test_temporal_table` strips `*` from
        # every cell before comparing -- so the writer must too, or the
        # bolded row looks like a row with no derivation. The asterisks are
        # preserved on output: they are the document's emphasis, not data.
        key = " | ".join(c.strip("*") for c in cells[:nkey])
        if key not in rows:
            missing.append(key)
            out.append(ln)
            continue
        want = list(rows[key])
        have = [c.strip("*") for c in cells[nkeep:]]
        # re-wrap each value in whatever emphasis its cell already carried
        marks = [c[:len(c) - len(c.lstrip("*"))] for c in cells[nkeep:]]
        # COMPARED WITH EMPHASIS STRIPPED FROM BOTH SIDES, or a derived
        # value that carries its own asterisks counts as changed on every
        # run: the document holds "**not evaluable**", `have` strips it to
        # "not evaluable", and the derivation returns it emphasised. That
        # made --check report two cells forever while rewriting nothing,
        # which is the drift alarm crying wolf at its own output.
        bare = [v.strip("*") for v in want]
        if have != bare:
            changed += sum(1 for a, b in zip(have, bare) if a != b)
            changed += abs(len(have) - len(bare))
        widths = [aligned.get(j, 0) for j in range(len(cells))]
        # AND A DERIVATION MAY CARRY ITS OWN EMPHASIS. The sequence
        # table's M1 cell is "**not evaluable**" when a protocol reaches
        # too few positives -- the asterisks come from the derivation, not
        # the document -- and re-applying the cell's existing marks on top
        # produced "****not evaluable****", which the test then rejected.
        # So a value that is already emphasised is used as it comes.
        want = [v if v.startswith("*") or not m else f"{m}{v}{m}"
                for v, m in zip(want, marks + [""] * 12)]
        new_cells = cells[:nkeep] + want
        indent = ln[:len(ln) - len(ln.lstrip())]
        out.append(f"{indent}| " + " | ".join(
            c.ljust(w) for c, w in zip(new_cells, widths + [0] * 12)).rstrip()
            + " |\n")
    if missing:
        raise SystemExit(
            f"{doc.name}: {len(missing)} table row(s) have no derivation and "
            f"were left untouched: {missing[:6]}. Either the key width is "
            f"wrong ({nkey} leading column(s) assumed) or the document "
            f"carries rows the results file no longer has -- both are "
            f"failures, not things to skip.")
    doc.write_text(head + marker + "".join(lines[:first]) + "".join(out)
                   + "".join(lines[last:]))
    return changed


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="report what would change and exit 1 if anything "
                         "would, writing nothing")
    args = ap.parse_args(argv)

    total = 0
    for name in MODULES:
        mod = importlib.import_module(name)
        marker, doc = mod.MARKER, mod.DOC
        rows = mod.derive()
        keep = getattr(mod, "KEEP_COLUMNS", None)
        if args.check:
            before = doc.read_text()
            n = rewrite(doc, marker, rows, keep)
            doc.write_text(before)
        else:
            n = rewrite(doc, marker, rows, keep)
        total += n
        print(f"  {name:<28} {doc.name:<20} {n:>3} cell(s) "
              f"{'would change' if args.check else 'rewritten'}")
    print(f"\n{total} cell(s) total")
    if args.check and total:
        print("A marked table disagrees with its own derivation. Run "
              "`make rewrite-tables`.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
