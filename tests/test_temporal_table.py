"""The temporal AUC table in docs/evaluation.md, re-derived cell by cell.

THIS TABLE HAD NO TEST FOR THREE WEEKS AND EVERY CELL IN IT WENT STALE.
It read M1 2.49/0.0033, M4 4.16/0.081, M7 0.85, COMBINED 3.36/0.0154 with
AUC 0.512 and degree 0.97/0.5909 with AUC 0.504 -- the figures of a release
two rebuilds earlier -- while the claim block at the top of the document
quoted that COMBINED/degree pair as the one-line summary of the whole
temporal protocol. Every other results table on the page is re-derived by a
test file; this one was the summary of them all and nothing could contradict
it, which is also how "degree is at chance on it" survived degree's interval
moving to [0.533, 0.578].

So: the document quotes what this file produces, and a regeneration that
moves any figure fails here until the document is corrected.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "evaluation.md"
RESULTS = ROOT / "data" / "results" / "temporal_2021_include.json"

MARKER = "<!-- TABLE: temporal_include_auc"
VIRUS = "SARS-CoV-2"
DASH = "—"
HEADER = ["scorer", "pool", "pos", "lift", "p", "AUC [95% CI]", ""]


# --------------------------------------------------------------- the rules
def fmt_p(p: float | None) -> str:
    """The document's own rule for a p-value, owned here and nowhere else."""
    return DASH if p is None else f"{p:.4f}"


def marker_of(m: dict) -> str:
    """The last column: why a row cannot be read, derived from the flags.

    `~ties` and `(!)` are the two mechanisms the prose above the table
    names, and a row that is both would be printed with both -- which no
    row currently is, so the order is fixed here rather than guessed.
    """
    marks = []
    if not m.get("hits_determined"):
        marks.append("~ties")
    if m.get("trustworthy") is False:
        marks.append("(!)")
    return " ".join(marks)


def derive() -> dict[str, list[str]]:
    """scorer -> the rendered cells the results file implies."""
    per = json.loads(RESULTS.read_text())["results"][VIRUS]
    rows = {}
    for scorer, m in per.items():
        if not isinstance(m, dict) or "auc" not in m:
            continue
        rows[scorer] = [
            f"{m['pool']:,}",
            f"{m['pos']:,}",
            f"{m['lift']:.2f}",
            fmt_p(m.get("p_value")),
            f"{m['auc']:.3f} [{m['auc_ci_lo']:.3f}, {m['auc_ci_hi']:.3f}]",
            marker_of(m),
        ]
    return rows


def parse() -> dict[str, list[str]]:
    """The table as the document prints it."""
    text = DOC.read_text()
    assert MARKER in text, "the temporal AUC table lost its marker"
    lines = [ln.strip() for ln in text.split(MARKER, 1)[1].splitlines()]
    start = next(i for i, ln in enumerate(lines) if ln.startswith("| scorer |"))
    header = [c.strip() for c in lines[start].strip("|").split("|")]
    # THE HEADER TOO: a column renamed or dropped would relabel the data.
    assert header == HEADER, header
    rows: dict[str, list[str]] = {}
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        cells = [c.strip() for c in line.strip("|").split("|")]
        assert len(cells) == 7, f"expected 7 columns: {line}"
        # Bold is emphasis, not data.
        scorer = cells[0].strip("*")
        assert scorer not in rows, (
            f"{scorer} appears twice in the table; keying rows by scorer "
            f"would hide the duplicate from every check here")
        rows[scorer] = [c.strip("*") for c in cells[1:]]
    return rows


# -------------------------------------------------------------- the checks
@pytest.fixture(scope="module")
def tables():
    if not RESULTS.exists():
        pytest.skip("the temporal split is not built in this checkout")
    return parse(), derive()


def test_the_table_covers_exactly_the_scorers_in_the_results(tables):
    """Set equality both ways. M5 reaches nothing on this split and is
    absent from the results, so it must be absent here; a channel that
    starts reaching must appear rather than be quietly dropped."""
    doc, got = tables
    assert set(doc) == set(got), (
        f"only in the document: {sorted(set(doc) - set(got))}; "
        f"only in the results: {sorted(set(got) - set(doc))}")


@pytest.mark.parametrize("column", range(6))
def test_every_cell_matches_the_results_file(tables, column):
    doc, got = tables
    labels = ["pool", "pos", "lift", "p", "AUC [95% CI]", "marker"]
    bad = {s: (doc[s][column], got[s][column])
           for s in sorted(set(doc) & set(got))
           if doc[s][column] != got[s][column]}
    assert not bad, (
        f"{labels[column]} is stale -- document vs re-derived: " +
        "; ".join(f"{s}: {d!r} != {g!r}" for s, (d, g) in bad.items()))


def test_the_table_is_marked_as_machine_checked():
    text = DOC.read_text()
    assert text.count(MARKER) == 1
    head = text.split(MARKER, 1)[1].split("-->", 1)[0]
    flat = " ".join(head.split())
    assert "tests/test_temporal_table.py" in flat
    assert "Do not hand-edit" in flat


def test_the_auc_sample_size_is_the_test_set_and_the_document_says_so():
    """THE DEFECT THAT PUT A WRONG BAR IN A PRE-REGISTRATION.

    `docs/prereg_temporal_auc.md` computed the AUC bar at "COMBINED's 48
    reached positives", but `evaluate_against_negatives` ranks the whole
    labelled test set and imputes the floor to everything a scorer cannot
    reach, so the n is the same for every scorer and is far larger. The
    document must print the real pair, because the difference between
    them is the difference between COMBINED clearing the bar and missing
    it by 0.074.
    """
    per = json.loads(RESULTS.read_text())["results"][VIRUS]
    ns = {(m["auc_n_pos"], m["auc_n_neg"]) for m in per.values()
          if isinstance(m, dict) and "auc_n_pos" in m}
    assert len(ns) == 1, f"scorers disagree on the AUC sample size: {ns}"
    n_pos, n_neg = ns.pop()
    text = " ".join(DOC.read_text().split())
    assert f"**{n_pos:,} and {n_neg:,}**" in text, (
        f"the AUC is computed over {n_pos:,} positives and {n_neg:,} "
        f"negatives and the document does not state that pair")
    # AND THE BAR THAT FOLLOWS FROM IT, not the one computed for the
    # reach count.
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.labels import auc_interval
    bar = 0.5
    while auc_interval(bar, n_pos, n_neg)[0] <= 0.5:
        bar += 0.0005
    assert f"is **{bar:.3f}**" in text, (
        f"the smallest AUC whose interval excludes 0.5 at n={n_pos:,}/"
        f"{n_neg:,} is {bar:.3f}; the document states a different bar")


def test_the_document_does_not_call_degree_at_chance_on_this_split():
    """It did, for three weeks, on a figure of 0.504. Degree's interval
    here excludes 0.5. The claim is allowed about the LIFT and must name
    it; unqualified, it is the sentence this test exists to forbid."""
    deg = json.loads(RESULTS.read_text())["results"][VIRUS]["degree"]
    if deg["auc_ci_lo"] <= 0.5 <= deg["auc_ci_hi"]:
        pytest.skip("degree's AUC contains chance on this split")
    flat = " ".join(DOC.read_text().split())
    # A WITHDRAWN CLAIM MAY BE QUOTED; IT MAY NOT BE ASSERTED. The
    # correction paragraphs quote these sentences verbatim in order to
    # retract them, so a flat substring search calls the retraction the
    # defect -- which it did on the first run of this test. Quoted spans
    # are stripped before searching, and the retraction is required to
    # survive that stripping.
    unquoted = re.sub(r'"[^"]*"', " ", flat)
    for claim in ("degree is at chance here",
                  "degree is at chance on it",
                  "Degree reads 0.504 on the same split"):
        assert claim not in unquoted, (
            f"degree reads {deg['auc']:.3f} "
            f"[{deg['auc_ci_lo']:.3f}, {deg['auc_ci_hi']:.3f}] on this "
            f"split, which excludes chance, and the document asserts "
            f"{claim!r} outside a quotation")
    assert "degree is at chance here" in flat, (
        "the retraction of that sentence has been deleted rather than "
        "the sentence corrected; the document should record what it "
        "used to claim")
