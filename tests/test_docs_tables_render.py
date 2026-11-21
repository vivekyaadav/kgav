"""A table cell past the header's column count is not displayed at all.

WHAT WAS WRONG. `doc2/records/audit_2026-11-03_fixes.md` declares its section
tables as three or four columns and its rows carried five and six, because
the Status column was appended to the rows and never to the header. GFM's
rule is explicit -- "if a row has greater [cells than the header], the
excess is ignored" -- so on 2026-11-17 **twenty rows rendered without their
status**, including every L-row closed in the preceding week. Checked with
pandoc rather than inferred from the spec: a document whose job is to say
what is still open was not showing whether anything was.

The same rule hides a literal pipe. `| E4 |` quoted
`seed_scheme: "per-row: Random(f'{seed}|{channel}|{virus}')"` inside
backticks; a code span does NOT protect a pipe in a table, so the sentence
was cut at `{seed}` and the rest of it became two phantom columns.

THIS IS THE REPOSITORY'S OWN DEFECT CLASS -- a counter that cannot be
printed is the same shape as a check that cannot fail, and a status that
cannot be read is the same shape as a row that was never closed. The two
tests below are the sweep that stops a twenty-first.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXES = ROOT / "doc2" / "records" / "audit_2026-11-03_fixes.md"

# A renderer splits on every unescaped pipe. It does not know about code
# spans, which is exactly how E4 lost half a sentence, so neither does this.
SPLIT = re.compile(r"(?<!\\)\|")


def rows(path: Path):
    """(line number, cells, header cell count) per table row in `path`."""
    header = None
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if not line.startswith("|"):
            header = None
            continue
        if set(line.replace("|", "").replace("-", "").strip()) <= set(" :"):
            continue                       # the |---|---| separator
        cells = [c.strip() for c in SPLIT.split(line.strip().strip("|"))]
        if header is None:
            header = len(cells)
            continue
        yield n, cells, header


def test_no_table_row_in_docs_loses_cells_when_rendered():
    """EVERY markdown document, not just the fix list: the defect is a
    property of the format and nothing stops another table acquiring it."""
    lost = []
    for p in sorted(ROOT.glob("docs/**/*.md")):
        for n, cells, header in rows(p):
            if len(cells) > header:
                dropped = " / ".join(cells[header:])[:70]
                lost.append(f"{p.relative_to(ROOT)}:{n} ({cells[0]}) "
                            f"header {header}, row {len(cells)} -- "
                            f"NOT DISPLAYED: {dropped}")
    assert not lost, (
        "these rows carry cells past their header's column count, and a "
        "renderer drops them silently -- widen the header, and escape any "
        "literal pipe as \\| even inside backticks:\n  "
        + "\n  ".join(lost))


def test_no_table_row_holds_an_unescaped_pipe_inside_a_code_span():
    """THE OTHER HALF OF THE SAME BUG, and test 1 only catches it when the
    stray cells happen to push past the header. A code span does NOT
    protect a pipe in a GFM table: `| E4 |` quoted a seed scheme reading
    f'{seed}|{channel}|{virus}' and the sentence was cut at {seed}. The
    escape is \\| and it works inside backticks."""
    import re as _re
    offenders = []
    for p in sorted(ROOT.glob("docs/**/*.md")):
        for n, line in enumerate(p.read_text().splitlines(), 1):
            if not line.startswith("|"):
                continue
            for m in _re.finditer(r"`[^`]*`", line):
                if _re.search(r"(?<!\\)\|", m.group(0)):
                    offenders.append(
                        f"{p.relative_to(ROOT)}:{n}  {m.group(0)[:60]}")
    assert not offenders, (
        "a pipe inside backticks still splits the table cell; write it as "
        "\\| :\n  " + "\n  ".join(offenders))
