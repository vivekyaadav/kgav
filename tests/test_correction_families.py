"""The multiple-testing sensitivity table in docs/evaluation.md (fix-list E1).

FOUR NUMBERS TYPED BY HAND. The table read `10 / 61 / 120 / 150` with an
alpha beside each, and nothing in the repository could recompute one of
them. Three of the four turned out right when they were finally derived on
2026-11-09; the point is that no release between them could have said so,
and the paragraph above the table has quoted 28, 51 and 23 rows for the
temporal protocol on three successive releases.

`kgav.baselines.correction_families` now derives all four for any run, both
scripts emit it, and this file asserts the document prints what the results
file holds.

WHAT IS NOT TESTED HERE, DELIBERATELY: which family is correct. That is a
pre-registration decision and this project's rule is that the family is
fixed before the p-values are in hand. The test asserts the sensitivity is
reported and accurate, and that the document says the choice has not been
made.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import (  # noqa: E402
    CORRECTION_FAMILIES,
    SIGNIFICANCE,
    correction_families,
)

DOC = ROOT / "docs" / "evaluation.md"
CROSS = ROOT / "data" / "results" / "metapath_comparison.json"
TEMPORAL = ROOT / "data" / "results" / "temporal_2021_include.json"

#: The document's label for each family, in the order the table prints them.
LABELS = {
    "readable_cells": "readable cells only",
    "cells_computed": "every cell computed",
    "cells_planned": "every cell planned",
    "planned_plus_baselines": "planned, plus COMBINED and degree",
}


@pytest.fixture(scope="module")
def cross() -> dict:
    if not CROSS.exists():
        pytest.skip("metapath_comparison.json is not built in this checkout")
    return json.loads(CROSS.read_text())


def parse_table() -> dict[str, tuple[int, str]]:
    """The sensitivity table as the document prints it: label -> (n, alpha)."""
    text = DOC.read_text()
    lines = [ln.strip() for ln in text.splitlines()]
    start = next(i for i, ln in enumerate(lines) if ln.startswith("| family |"))
    rows: dict[str, tuple[int, str]] = {}
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        cells = [c.strip() for c in line.strip("|").split("|")]
        assert len(cells) == 4, f"expected 4 columns: {line}"
        n = int(cells[1].replace(",", ""))
        rows[cells[0]] = (n, cells[2])
    return rows


# -------------------------------------------------------------- the helper
def test_the_helper_reproduces_the_four_family_sizes(cross):
    """The sizes the document has carried since 2026-10-28, derived.

    Three of these were right when typed and one -- `cells_computed` at 61
    -- had a definition nobody had written down: it is every (virus,
    metapath) cell the run PRODUCED, which is not the same as every cell
    with a p-value (34 of them) and not the same as every cell planned
    (120). An undocumented definition that happens to give the right
    number is the thing this file exists to retire.
    """
    got = correction_families(cross["per_virus"], 8)
    sizes = {k: v["n"] for k, v in got["families"].items()}
    assert sizes == {"readable_cells": 10, "cells_computed": 61,
                     "cells_planned": 120, "planned_plus_baselines": 150}, sizes


def test_a_placeholder_row_is_not_a_test(cross):
    """THE TEMPORAL OUTPUT'S TWO PLACEHOLDER SHAPES.

    A virus refused for too few reachable positives is recorded as
    `{"not_evaluable": {...}}` and a channel that could not run as
    `{"not_evaluated": {...}}` -- both dicts sitting exactly where a
    channel's metrics would be. Counting them would inflate every family
    with tests nobody performed, and on the temporal split twelve of
    fifteen viruses are refused, so the error would be large rather than
    marginal.
    """
    if not TEMPORAL.exists():
        pytest.skip("the temporal split is not built in this checkout")
    per = json.loads(TEMPORAL.read_text())["results"]
    refused = [v for v, ch in per.items()
               if isinstance(ch, dict) and "not_evaluable" in ch]
    assert refused, "no virus is refused, so this test proves nothing"
    got = correction_families(per, 8)
    scored = len(per) - len(refused)
    assert got["families"]["cells_planned"]["n"] == scored * 8
    # And no placeholder key reached the computed family.
    assert got["families"]["cells_computed"]["n"] <= scored * 8


def test_alpha_is_the_uncorrected_alpha_over_the_family_size(cross):
    got = correction_families(cross["per_virus"], 8)
    assert got["uncorrected_alpha"] == SIGNIFICANCE
    for name, f in got["families"].items():
        assert f["alpha"] == pytest.approx(SIGNIFICANCE / max(1, f["n"])), name


def test_survivors_are_only_ever_readable_cells(cross):
    """An unreadable cell surviving a correction is still unreadable.

    Listing one would invite the reader to count it, which is the defect
    the `trustworthy` and `hits_determined` gates exist to prevent -- and
    the temporal run has four rows that survive Bonferroni and fail both
    gates, so this is a live case rather than a hypothetical.
    """
    per = cross["per_virus"]
    got = correction_families(per, 8)
    readable = {f"{v} {n}" for v, ch in per.items() for n, m in ch.items()
                if isinstance(m, dict) and m.get("hits_determined")
                and m.get("trustworthy")}
    for name, f in got["families"].items():
        assert set(f["survivors"]) <= readable, name


# ------------------------------------------------------------ the document
def test_the_table_matches_the_results_file(cross):
    """Every n and alpha in the document, re-derived.

    Reads the committed `correction` block rather than recomputing from
    `per_virus`, so a script that stopped emitting it fails here instead
    of leaving the table checked against a function nothing calls.
    """
    assert "correction" in cross, (
        "metapath_comparison.json carries no `correction` block; re-run "
        "scripts/compare_metapaths.py")
    doc = parse_table()
    for key in CORRECTION_FAMILIES:
        label = LABELS[key]
        row = next((v for k, v in doc.items() if k.startswith(label)), None)
        assert row is not None, (
            f"the sensitivity table has no row for {label!r}; it has "
            f"{sorted(doc)}")
        n, alpha = row
        want = cross["correction"]["families"][key]
        assert n == want["n"], f"{label}: document says n={n}, run gives {want['n']}"
        assert alpha == f"{want['alpha']:.5f}", (
            f"{label}: document says alpha={alpha}, run gives "
            f"{want['alpha']:.5f}")


def test_the_verdict_column_follows_from_the_alpha(cross):
    """"clears" and "fails" are not editorial.

    The table's fourth column judges one cell -- the project's single
    nominally-significant readable result -- against each family's
    threshold, and the cell's own p-value is in the header. Both are
    derived so a p-value that moves cannot leave a verdict standing.
    """
    doc = parse_table()
    text = DOC.read_text()
    lines = [ln.strip() for ln in text.splitlines()]
    start = next(i for i, ln in enumerate(lines) if ln.startswith("| family |"))
    m = re.search(r"at p=([0-9.]+)", lines[start])
    assert m, f"the table header names no p-value: {lines[start]}"
    p_value = float(m.group(1))
    # The cell that p belongs to, found rather than assumed.
    cells = [(v, n, mm) for v, ch in cross["per_virus"].items()
             for n, mm in ch.items()
             if isinstance(mm, dict) and mm.get("p_value") is not None
             and f"{mm['p_value']:.6f}" == f"{p_value:.6f}"]
    assert cells, (
        f"the table judges a cell at p={p_value} and no cell in the results "
        f"file reads that p-value")
    for key in CORRECTION_FAMILIES:
        label = LABELS[key]
        row_key = next(k for k in doc if k.startswith(label))
        alpha = float(doc[row_key][1])
        verdict = next(
            ln for ln in lines[start + 2:] if ln.startswith(f"| {row_key}"))
        clears = p_value < alpha
        assert ("clears" in verdict) == clears, (
            f"{label}: alpha={alpha}, p={p_value}, so the cell "
            f"{'clears' if clears else 'fails'} -- the table says "
            f"otherwise: {verdict}")


def test_the_document_says_the_family_is_not_chosen_here():
    """The stance, asserted. Reporting a sensitivity and then quietly
    reading the smallest family as the answer is the ordering this
    project's rule forbids, and the table is one sentence away from
    inviting it."""
    flat = " ".join(DOC.read_text().split())
    assert "belongs in the prereg" in flat, (
        "the sensitivity table does not say that the choice of family "
        "belongs in a pre-registration")
