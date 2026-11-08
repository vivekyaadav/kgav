"""The AUC decomposition table in docs/evaluation.md, re-derived.

WHY THIS TABLE EXISTS AT ALL. The claim block prints degree at 0.5767 and
the COMBINED ensemble at 0.5746 and called them indistinguishable, two
thousandths apart. They are not comparable: every AUC on that page is
imputed, so a scorer reaching a larger share of the positives than of the
negatives scores above chance with no ranking ability, and degree reaches
everything while COMBINED reaches 43.8% of the positives. 86% of COMBINED's
margin above chance is that gap. The document now prints the floor, the
ceiling and the fraction of the span between them each scorer used, and
this file asserts every cell of it against the results file.

ONE IMPLEMENTATION OF THE ARITHMETIC. The floor, ceiling and skill are
`kgav.labels.auc_floor`, `auc_ceiling` and `auc_skill` called here -- not
re-derived in the test -- because a second copy of a formula in a test is
how the document and the code come to disagree while both pass. What this
file owns is the RENDERING: how many decimals, how a percentage is written,
and what `n/a` means.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.labels import auc_ceiling, auc_floor, auc_skill  # noqa: E402

DOC = ROOT / "docs" / "evaluation.md"
RESULTS = ROOT / "data" / "results" / "hard_negatives.json"
MARKER = "<!-- TABLE: auc_decomposition"
VIRUS = "SARS-CoV-2"
HEADER = ["scorer", "cov_pos", "cov_neg", "floor", "AUC", "ceiling", "skill"]


def reached(m: dict) -> tuple[int, int, int, int]:
    """(reached_pos, reached_neg, n_pos, n_neg) for one cell.

    The results file stores coverages rather than counts, so the counts
    are recovered by multiplying out. `round` and not `int`: 0.43820224
    * 89 is 38.9999..., and truncating it would silently move COMBINED's
    floor by a thousandth and make this test's own figures wrong.
    """
    n_pos, n_neg = m["n_pos"], m["n_neg"]
    return (round(m["coverage_pos"] * n_pos), round(m["coverage_neg"] * n_neg),
            n_pos, n_neg)


def derive() -> dict[str, list[str]]:
    """scorer -> rendered cells, from the results file and kgav.labels."""
    per = json.loads(RESULTS.read_text())["cross_sectional"][VIRUS]
    rows = {}
    for scorer, m in per.items():
        if not isinstance(m, dict) or "auc" not in m:
            continue
        rp, rn, n_pos, n_neg = reached(m)
        floor = auc_floor(rp, rn, n_pos, n_neg)
        ceiling = auc_ceiling(rp, rn, n_pos, n_neg)
        skill = auc_skill(m["auc"], floor, ceiling)
        rows[scorer] = [
            f"{m['coverage_pos']:.1%}",
            f"{m['coverage_neg']:.1%}",
            f"{floor:.4f}",
            f"{m['auc']:.4f}",
            f"{ceiling:.4f}",
            "n/a" if skill is None else f"{skill:.1%}",
        ]
    return rows


def parse() -> dict[str, list[str]]:
    text = DOC.read_text()
    assert MARKER in text, "the decomposition table lost its marker"
    lines = [ln.strip() for ln in text.split(MARKER, 1)[1].splitlines()]
    start = next(i for i, ln in enumerate(lines) if ln.startswith("| scorer |"))
    header = [c.strip() for c in lines[start].strip("|").split("|")]
    assert header == HEADER, header
    rows: dict[str, list[str]] = {}
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        cells = [c.strip().strip("*") for c in line.strip("|").split("|")]
        assert len(cells) == 7, f"expected 7 columns: {line}"
        assert cells[0] not in rows, f"{cells[0]} appears twice"
        rows[cells[0]] = cells[1:]
    return rows


@pytest.fixture(scope="module")
def tables():
    if not RESULTS.exists():
        pytest.skip("hard_negatives.json is not built in this checkout")
    return parse(), derive()


def test_the_table_covers_exactly_the_scorers_in_the_results(tables):
    doc, got = tables
    assert set(doc) == set(got), (
        f"only in the document: {sorted(set(doc) - set(got))}; "
        f"only in the results: {sorted(set(got) - set(doc))}")


@pytest.mark.parametrize("column", range(6))
def test_every_cell_matches_the_results_file(tables, column):
    doc, got = tables
    labels = ["cov_pos", "cov_neg", "floor", "AUC", "ceiling", "skill"]
    bad = {s: (doc[s][column], got[s][column])
           for s in sorted(set(doc) & set(got))
           if doc[s][column] != got[s][column]}
    assert not bad, (
        f"{labels[column]} is stale -- document vs re-derived: " +
        "; ".join(f"{s}: {d!r} != {g!r}" for s, (d, g) in bad.items()))


def test_the_floor_reduces_to_the_coverage_gap():
    """The long form and the short form agree.

    `auc_floor` counts pairs; the document explains it as
    `0.5 + (cov_pos - cov_neg) / 2`. If those ever disagree the prose is
    wrong, the code is wrong, or both -- and a reader checking the table
    by hand would use the short form.
    """
    per = json.loads(RESULTS.read_text())["cross_sectional"][VIRUS]
    for scorer, m in per.items():
        if not isinstance(m, dict) or "auc" not in m:
            continue
        rp, rn, n_pos, n_neg = reached(m)
        short = 0.5 + (rp / n_pos - rn / n_neg) / 2
        assert abs(auc_floor(rp, rn, n_pos, n_neg) - short) < 1e-12, scorer


def test_degree_has_no_coverage_advantage_and_the_document_says_so():
    """THE POINT OF THE WHOLE SECTION, asserted rather than written.

    Degree ranks every compound, so its floor must be exactly 0.5 and
    all of its margin is ranking. If a future release gives degree
    partial coverage, the comparison this section draws changes shape
    and the prose has to be rewritten rather than left standing.
    """
    per = json.loads(RESULTS.read_text())["cross_sectional"][VIRUS]
    deg = per["degree"]
    assert deg["coverage_pos"] == 1.0 and deg["coverage_neg"] == 1.0
    rp, rn, n_pos, n_neg = reached(deg)
    assert auc_floor(rp, rn, n_pos, n_neg) == 0.5
    flat = " ".join(DOC.read_text().split())
    assert "Degree's floor is exactly 0.5" in flat


def test_the_share_of_combineds_margin_that_is_coverage_is_the_files_own():
    """The 86% figure, derived. It is the number a reader takes away from
    this section, so it cannot be a figure that was typed once."""
    per = json.loads(RESULTS.read_text())["cross_sectional"][VIRUS]
    m = per["COMBINED"]
    rp, rn, n_pos, n_neg = reached(m)
    floor = auc_floor(rp, rn, n_pos, n_neg)
    share = (floor - 0.5) / (m["auc"] - 0.5)
    flat = " ".join(DOC.read_text().split())
    assert f"{share:.0%}" in flat, (
        f"{share:.1%} of COMBINED's margin above chance is its coverage "
        f"gap and the document does not state that share")
    # AND THE TWO COMPONENTS, so the percentage can be checked.
    assert f"{floor - 0.5:.4f}" in flat
    assert f"{m['auc'] - floor:.4f}" in flat


def test_the_temporal_protocol_still_cannot_be_decomposed():
    """The half of B2 that is not done, asserted so it cannot be
    forgotten or quietly claimed. When temporal_split.py starts recording
    coverages this fails, and the document must then either decompose
    that table too or say why not."""
    path = ROOT / "data" / "results" / "temporal_2021_include.json"
    if not path.exists():
        pytest.skip("the temporal split is not built in this checkout")
    rows = json.loads(path.read_text())["results"]["SARS-CoV-2"]
    missing = [s for s, m in rows.items()
               if isinstance(m, dict) and "auc" in m
               and "coverage_pos" not in m]
    assert missing, (
        "the temporal results now record coverage; the decomposition can "
        "be computed for that table and the document must do it rather "
        "than say it cannot")
    flat = " ".join(DOC.read_text().split())
    assert "does not record them" in flat
