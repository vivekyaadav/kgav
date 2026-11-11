"""The sequence table in docs/evaluation.md, re-derived rather than trusted.

WHAT THIS IS FOR. That table is the document's argument: M1 falls as each
control is applied, and every later section reasons from those five rows. On
2026-10-30 four of them were found to report numbers the graph had stopped
producing -- row 1 read 0.823 against a measured 0.778, and the negative set
had more than doubled since the figures were taken. Row 5 was the only one
that still reproduced, because its labels come from one screen the rebuild
did not touch.

WHY NOTHING CAUGHT IT. check-results compares a results FILE's stamp against
the graph, and the file was current. The document is not the file. The
"re-derived by" map names hard_negatives.py beside these figures, but naming
a script is not running it. And every protocol wrote one
hard_negatives.json, so by the end of a rebuild four of the five rows had no
artifact left to check against at all.

So each protocol now writes its own file and this re-derives the table from
all five, hard-coding no figure. The comparison itself is a pure function
over (cells, document) and is exercised on synthetic input below, so the
rules are tested even in a checkout with no release on disk.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "evaluation.md"
RESULTS = ROOT / "data" / "results"
sys.path.insert(0, str(ROOT / "src"))

from kgav.calibration import MIN_REACHED_POSITIVES
from kgav.labels import auc_interval

MARKER = "<!-- TABLE: hard_negatives_sequence"
VIRUS = "SARS-CoV-2"
CHANNEL = "M1"
NOT_EVALUABLE = "**not evaluable**"

#: table row number -> the results file that must reproduce it.
#:
#: Row 5 keeps the canonical name because calibrate.py and the agent read
#: it, and it is the protocol the block runs last. Listed one per row rather
#: than derived from a pattern: a single `hard_negatives_row\d` entry would
#: still be satisfied with one of the four files missing, which is the
#: membership-satisfied-elsewhere defect this suite keeps finding.
SEQUENCE_FILES = {
    "1": "hard_negatives_row1.json",
    "2": "hard_negatives_row2.json",
    "3": "hard_negatives_row3.json",
    "4": "hard_negatives_row4.json",
    "5": "hard_negatives.json",
}

HEADER = ["#", "Protocol", "pos", "neg", "cov_pos", "cov_neg", "reached",
          "**M1**"]


# ----------------------------------------------------------------- the rules
def m1_cell(auc: float, n_pos: int, n_neg: int, reached: int) -> str:
    """The M1 cell a results document implies.

    ONE FORMATTING RULE, owned here. A second copy in the document's prose
    is how the figures drifted the first time.
    """
    if reached < MIN_REACHED_POSITIVES:
        return NOT_EVALUABLE
    lo, hi = auc_interval(auc, n_pos, n_neg)
    return f"{auc:.3f} [{lo:.3f}, {hi:.3f}]"


def derive_row(doc: dict, virus: str = VIRUS, channel: str = CHANNEL
               ) -> list[str]:
    """The pos..M1 cells one results document implies, as rendered text."""
    per = doc["cross_sectional"][virus]
    m = per[channel]
    pos, neg = m["n_pos"], m["n_neg"]
    cov_pos, cov_neg = m["coverage_pos"], m["coverage_neg"]
    reached = round(cov_pos * pos)
    return [f"{pos:,}", f"{neg:,}", f"{cov_pos:.1%}", f"{cov_neg:.1%}",
            f"{reached:,}", m1_cell(m["auc"], pos, neg, reached)]


def mismatches(cells: list[str], doc: dict) -> list[str]:
    """Which of a table row's cells the results document does not produce.

    Returns a list rather than asserting, so one failure can report every
    wrong cell instead of only the leftmost -- when this table drifts it
    drifts wholesale, and fixing it one cell per test run is how a stale
    table survives a correction.
    """
    want = derive_row(doc)
    names = ["pos", "neg", "cov_pos", "cov_neg", "reached", "M1"]
    return [f"{n}: document says {g!r}, results give {w!r}"
            for n, g, w in zip(names, cells, want) if g != w]


KEEP_COLUMNS = 2        #: '#' and the protocol prose are the document's


def derive() -> dict[str, list[str]]:
    """The table the five results files imply: row number -> pos..M1 cells.

    FOR scripts/rewrite_derived_tables.py, so this table is written from
    its own derivation rather than by hand against `mismatches()`. The
    keys are row numbers; the protocol column beside them is prose no
    results file produces, which is why this module declares
    `KEEP_COLUMNS = 2`.

    A row whose results file is absent is LEFT OUT rather than filled with
    a placeholder, and the rewriter refuses a document row it cannot
    derive -- so a missing file stops the rewrite by name instead of
    quietly publishing whatever was there.
    """
    return {n: derive_row(doc) for n, doc in sorted(_available().items())}


def parse() -> dict[str, list[str]]:
    """The table as the document prints it: row number -> pos..M1 cells."""
    text = DOC.read_text()
    assert MARKER in text, "the sequence table lost its machine-readable marker"
    lines = [ln.strip() for ln in text.split(MARKER, 1)[1].splitlines()]
    start = next(i for i, ln in enumerate(lines) if ln.startswith("| # |"))
    header = [c.strip() for c in lines[start].strip("|").split("|")]
    # THE HEADER TOO. A column deleted from the header while the rows stay
    # eight wide would relabel the data with nothing to say so.
    assert header == HEADER, header
    rows: dict[str, list[str]] = {}
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        cells = [c.strip() for c in line.strip("|").split("|")]
        assert len(cells) == 8, f"expected 8 columns, got {len(cells)}: {line}"
        rows[cells[0]] = cells[2:]          # drop '#' and the Protocol prose
    return rows


def _available() -> dict[str, dict]:
    return {n: json.loads((RESULTS / f).read_text())
            for n, f in SEQUENCE_FILES.items() if (RESULTS / f).exists()}


# ------------------------------------------------- against the real results
@pytest.fixture(scope="module")
def documents() -> dict[str, dict]:
    docs = _available()
    if not docs:
        pytest.skip("no hard_negatives results on disk")
    return docs


@pytest.mark.parametrize("row", sorted(SEQUENCE_FILES))
def test_every_cell_of_the_row_matches_its_results_file(row, documents):
    if row not in documents:
        pytest.skip(f"{SEQUENCE_FILES[row]} not built in this checkout")
    table = parse()
    assert row in table, f"the table has no row {row}"
    bad = mismatches(table[row], documents[row])
    assert not bad, (
        f"table row {row} does not match {SEQUENCE_FILES[row]}:\n  "
        + "\n  ".join(bad)
        + "\n\nRe-derive it; do not edit the number to match the other one.")


def test_the_table_has_a_row_for_every_file_and_a_file_for_every_row(documents):
    """Both ways. A protocol added to the block with no row here would go
    unreported, and a row here with no protocol is a figure nothing
    re-derives -- which is what this whole file exists to stop."""
    assert set(parse()) == set(SEQUENCE_FILES), (
        f"table rows {sorted(parse())} vs files {sorted(SEQUENCE_FILES)}")
    missing = sorted(set(SEQUENCE_FILES) - set(documents))
    # NOT A FAILURE WHEN NONE OF THEM EXIST. Before the block is first run
    # with --results-name, rows 1-4 have never been written and failing here
    # would mean a fresh checkout cannot pass its own suite. A PARTIAL set is
    # the dangerous state -- an interrupted run, or a protocol quietly
    # dropped from the block -- and that still fails.
    per_row = sorted(set(SEQUENCE_FILES) - {"5"})
    if sorted(set(missing)) == per_row:
        pytest.skip("rows 1-4 have not been written yet; run the Reproducing "
                    "block, which now writes one file per protocol")
    assert not missing, (
        f"rows {missing} have no results file while others do -- a partial "
        f"run, or a protocol dropped from the block. Re-run the Reproducing "
        f"block, which writes one file per protocol")


def test_each_row_came_from_the_protocol_the_table_says_it_did(documents):
    """A file-to-row mapping is only as good as the protocol inside the
    file. Row 3 reading row 2's document would match plenty of cells."""
    expect = {
        "1": ("selective-only", "positives", False, None),
        "2": ("all", "positives", False, None),
        "3": ("verified-only", "both", False, None),
        "4": ("verified-only", "both", True, None),
        "5": ("all", "positives", True, ["infores:ncats-opendata"]),
    }
    for row, doc in documents.items():
        sel, applies, pub, src = expect[row]
        got = (doc["selectivity_filter"], doc["selectivity_applies"],
               doc["publication_disjoint"], doc["label_source"])
        assert got == (sel, applies, pub, src), f"row {row}: {got}"


def test_the_files_describe_the_graph_on_disk(documents):
    """A guard asserting against a stale artifact guards nothing."""
    ids = {d["provenance"]["release"]["id"] for d in documents.values()}
    shas = {d["provenance"]["release"]["edges_sha"] for d in documents.values()}
    assert len(ids) == 1 and len(shas) == 1, (
        f"the five rows describe different graphs: {ids}, {shas}")


def test_the_table_is_marked_as_machine_checked():
    text = DOC.read_text()
    assert text.count(MARKER) == 1
    head = text.split(MARKER, 1)[1].split("-->", 1)[0]
    assert "tests/test_sequence_table.py" in head
    assert "Do not hand-edit" in head


# ------------------------------------------- the rules, on synthetic input
#
# Everything above skips in a checkout with no release. These do not: the
# formatting and comparison rules are what drifted, so they are tested on
# documents built here rather than only on whatever happens to be on disk.
def _doc(auc=0.7, n_pos=100, n_neg=200, cov_pos=0.5, cov_neg=0.1):
    return {"cross_sectional": {VIRUS: {CHANNEL: {
        "auc": auc, "n_pos": n_pos, "n_neg": n_neg,
        "coverage_pos": cov_pos, "coverage_neg": cov_neg}}}}


def test_a_matching_row_reports_nothing():
    d = _doc()
    assert mismatches(derive_row(d), d) == []


def test_every_wrong_cell_is_reported_not_just_the_first():
    """One run must show the whole drift. Reporting the leftmost cell only
    is how a stale table survives being corrected cell by cell."""
    d = _doc()
    cells = derive_row(d)
    cells[0] = "999"
    cells[5] = "0.999 [0.9, 1.0]"
    bad = mismatches(cells, d)
    assert len(bad) == 2, bad
    assert bad[0].startswith("pos:") and bad[1].startswith("M1:")


def test_a_channel_below_the_floor_is_not_given_an_auc():
    """The row-4 case. 5 reached positives is an AUC of the imputed floor,
    and the document's own MIN_REACHED_POSITIVES says so."""
    below = _doc(auc=0.497, n_pos=974, cov_pos=5 / 974)
    assert derive_row(below)[-1] == NOT_EVALUABLE
    at_floor = _doc(auc=0.497, n_pos=974, cov_pos=MIN_REACHED_POSITIVES / 974)
    assert derive_row(at_floor)[-1] != NOT_EVALUABLE


def test_the_interval_is_the_one_the_codebase_computes():
    """Not a second Hanley-McNeil implementation living in a test."""
    lo, hi = auc_interval(0.7, 100, 200)
    assert derive_row(_doc())[-1] == f"0.700 [{lo:.3f}, {hi:.3f}]"


def test_reached_is_derived_from_the_coverage_and_not_read_separately():
    assert derive_row(_doc(n_pos=974, cov_pos=0.4804928131416838))[4] == "468"


@pytest.mark.parametrize("field,value", [
    ("n_pos", 101), ("n_neg", 201), ("coverage_pos", 0.51),
    ("coverage_neg", 0.11), ("auc", 0.71),
])
def test_a_results_file_that_moves_any_field_fails_the_row(field, value):
    """Each column independently. A comparison that only looked at the AUC
    would have passed the drift found on 2026-10-30, where pos and neg
    moved while row 4's AUC barely did."""
    cells = derive_row(_doc())
    moved = _doc()
    moved["cross_sectional"][VIRUS][CHANNEL][field] = value
    assert mismatches(cells, moved), f"{field} moving to {value} went unnoticed"


def test_the_parser_reads_the_marked_table_and_not_another():
    """There are several tables in the document. Taking the first would
    check whichever happened to come earliest."""
    rows = parse()
    assert sorted(rows) == ["1", "2", "3", "4", "5"], sorted(rows)
    assert all(len(c) == 6 for c in rows.values())
    # The Protocol column is dropped, so no cell should be prose.
    assert not any(re.search(r"[a-z]{4}", c) and "evaluable" not in c
                   for cells in rows.values() for c in cells)
