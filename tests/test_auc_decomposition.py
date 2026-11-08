"""The AUC decomposition table in docs/evaluation.md, re-derived.

WHY THIS TABLE EXISTS AT ALL. The claim block prints degree at 0.5767 and
the COMBINED ensemble at 0.5725 and called them indistinguishable, four
thousandths apart. They are not comparable: every AUC on that page is
imputed, so a scorer reaching a larger share of the positives than of the
negatives scores above chance with no ranking ability, and degree reaches
everything while COMBINED reaches 43.8% of the positives. 89% of COMBINED's
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


# ------------------------------------------------ the temporal decomposition
TEMPORAL_MARKER = "<!-- TABLE: temporal_auc_decomposition"
TEMPORAL_HEADER = ["policy", "virus", "scorer", "cov_pos", "cov_neg",
                   "floor", "AUC", "ceiling", "skill"]
POLICIES = ("include", "computed_only", "exclude")


def temporal_rows() -> dict[tuple[str, str, str], list[str]]:
    """(policy, virus, scorer) -> rendered cells, for every temporal cell
    whose AUC interval excludes 0.5.

    WHICH ROWS BELONG HERE IS DERIVED, not listed. The table's claim is
    "every temporal cell that clears chance", so a cell that starts
    clearing must appear and one that stops must go; a hand-written row
    set would make the claim unfalsifiable, which is the defect the
    cross-sectional table had for three weeks.

    The coverages come from `auc_coverage_pos` / `auc_coverage_neg`.
    THOSE FIELDS WERE ALWAYS THERE. On 2026-11-08 this document, the
    prereg and the fix list all stated that `temporal_split.py` records
    `auc_n_pos` and `auc_n_neg` "but not the coverages", so the one
    interval in this project excluding chance had "no known floor
    beneath it". It had a floor equal to its AUC to four decimals. The
    claim was wrong by one grep over the writer, not the data.
    """
    out = {}
    for policy in POLICIES:
        path = ROOT / "data" / "results" / f"temporal_2021_{policy}.json"
        if not path.exists():
            continue
        for virus, channels in json.loads(path.read_text())["results"].items():
            for scorer, m in channels.items():
                if not isinstance(m, dict) or m.get("auc_ci_lo") is None:
                    continue
                if m["auc_ci_lo"] <= 0.5:
                    continue
                n_pos, n_neg = m["auc_n_pos"], m["auc_n_neg"]
                cp, cn = m["auc_coverage_pos"], m["auc_coverage_neg"]
                rp, rn = round(cp * n_pos), round(cn * n_neg)
                floor = auc_floor(rp, rn, n_pos, n_neg)
                ceiling = auc_ceiling(rp, rn, n_pos, n_neg)
                skill = auc_skill(m["auc"], floor, ceiling)
                out[(policy, virus, scorer)] = [
                    f"{cp:.1%}", f"{cn:.1%}", f"{floor:.4f}",
                    f"{m['auc']:.4f}", f"{ceiling:.4f}",
                    "n/a" if skill is None else f"{skill:.1%}"]
    return out


def parse_temporal() -> dict[tuple[str, str, str], list[str]]:
    """The table as the document prints it, or {} when it prints none.

    EMPTY IS A STATE THE DOCUMENT MUST ASSERT, not one it may fall into
    by leaving a table out. On 2026-11-09 the C4 fix removed all five
    temporal cells that cleared chance, so there is no table -- and a
    parser that simply found no rows would read the same whether the
    document said so or whether somebody had deleted the section.
    """
    text = DOC.read_text()
    assert TEMPORAL_MARKER in text, "the temporal decomposition lost its marker"
    body = text.split(TEMPORAL_MARKER, 1)[1]
    lines = [ln.strip() for ln in body.splitlines()]
    if not any(ln.startswith("| policy |") for ln in lines):
        flat = " ".join(body.split())
        assert "No temporal cell has an AUC interval excluding 0.5" in flat, (
            "the temporal decomposition prints no table and the document "
            "does not say that no temporal cell clears chance")
        return {}
    start = next(i for i, ln in enumerate(lines) if ln.startswith("| policy |"))
    header = [c.strip() for c in lines[start].strip("|").split("|")]
    assert header == TEMPORAL_HEADER, header
    rows = {}
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        cells = [c.strip().strip("*") for c in line.strip("|").split("|")]
        assert len(cells) == 9, f"expected 9 columns: {line}"
        key = (cells[0], cells[1], cells[2])
        assert key not in rows, f"{key} appears twice"
        rows[key] = cells[3:]
    return rows


@pytest.fixture(scope="module")
def temporal_tables():
    if not (ROOT / "data/results/temporal_2021_include.json").exists():
        pytest.skip("the temporal split is not built in this checkout")
    return parse_temporal(), temporal_rows()


def test_the_temporal_table_lists_exactly_the_cells_that_clear_chance(
        temporal_tables):
    doc, got = temporal_tables
    assert set(doc) == set(got), (
        f"only in the document: {sorted(set(doc) - set(got))}; "
        f"only in the results: {sorted(set(got) - set(doc))}")


@pytest.mark.parametrize("column", range(6))
def test_every_temporal_cell_matches_the_results(temporal_tables, column):
    doc, got = temporal_tables
    labels = ["cov_pos", "cov_neg", "floor", "AUC", "ceiling", "skill"]
    bad = {k: (doc[k][column], got[k][column])
           for k in sorted(set(doc) & set(got))
           if doc[k][column] != got[k][column]}
    assert not bad, (
        f"{labels[column]} is stale -- document vs re-derived: " +
        "; ".join(f"{k}: {d!r} != {g!r}" for k, (d, g) in bad.items()))


def test_the_temporal_auc_is_a_coverage_statistic():
    """THE FINDING THAT REPLACED THE ONE BELOW, on 2026-11-09.

    Most temporal AUC cells equal their own imputation floor exactly, and
    none of them has more than a few hundredths of AUC available between
    floor and ceiling. The document's claim -- that the tie-safe AUC on
    this protocol CANNOT separate, rather than merely does not -- rests
    on those two numbers, so both are derived here.
    """
    cells = []
    for policy in POLICIES:
        path = ROOT / "data" / "results" / f"temporal_2021_{policy}.json"
        if not path.exists():
            continue
        for channels in json.loads(path.read_text())["results"].values():
            if not isinstance(channels, dict):
                continue
            for m in channels.values():
                if not isinstance(m, dict) or m.get("auc") is None:
                    continue
                if m.get("auc_coverage_pos") is None:
                    continue
                n_pos, n_neg = m["auc_n_pos"], m["auc_n_neg"]
                rp = round(m["auc_coverage_pos"] * n_pos)
                rn = round(m["auc_coverage_neg"] * n_neg)
                floor = auc_floor(rp, rn, n_pos, n_neg)
                cells.append((m["auc"], floor,
                              auc_ceiling(rp, rn, n_pos, n_neg)))
    if not cells:
        pytest.skip("the temporal split is not built in this checkout")
    flat = " ".join(DOC.read_text().split())
    exact = sum(1 for auc, fl, _ in cells if f"{fl:.4f}" == f"{auc:.4f}")
    assert f"{len(cells)} TEMPORAL AUC CELLS, {exact} EQUAL" in flat.upper(), (
        f"{exact} of {len(cells)} temporal AUC cells equal their floor to "
        f"four decimals; the document states a different pair")
    worst = max(abs(auc - fl) for auc, fl, _ in cells)
    assert f"{worst:.4f}" in flat, (
        f"the largest deviation from the floor is {worst:.4f} and the "
        f"document does not state it")
    span = max(ce - fl for _, fl, ce in cells)
    assert f"{span:.4f}" in flat, (
        f"the largest floor-to-ceiling span is {span:.4f} and the document "
        f"does not state it")


def test_the_one_chance_clearing_ensemble_cell_is_its_own_floor():
    """THE FINDING OF 2026-11-09, asserted rather than written down.

    COMBINED on the temporal include split is the only non-degree cell
    in this project whose AUC interval excludes 0.5, and its AUC equals
    its imputation floor to four decimals with a ceiling 0.0003 above.
    If a rebuild ever separates those three numbers, the document's
    "artifact of reach" conclusion stops being true and must be
    rewritten rather than left standing.
    """
    path = ROOT / "data" / "results" / "temporal_2021_include.json"
    if not path.exists():
        pytest.skip("the temporal split is not built in this checkout")
    m = json.loads(path.read_text())["results"]["SARS-CoV-2"]["COMBINED"]
    n_pos, n_neg = m["auc_n_pos"], m["auc_n_neg"]
    rp = round(m["auc_coverage_pos"] * n_pos)
    rn = round(m["auc_coverage_neg"] * n_neg)
    floor = auc_floor(rp, rn, n_pos, n_neg)
    ceiling = auc_ceiling(rp, rn, n_pos, n_neg)
    flat = " ".join(DOC.read_text().split())
    # STILL ITS FLOOR -- that part never stopped being true, and it is
    # the whole reason the cell was never a result.
    assert f"{floor:.4f}" == f"{m['auc']:.4f}", (
        f"the AUC {m['auc']:.4f} and the floor {floor:.4f} no longer "
        f"agree to four decimals")
    assert ceiling - floor < 0.001
    if m["auc_ci_lo"] > 0.5:
        # The state this test was written for, 2026-11-08: the cell
        # cleared chance ON ITS FLOOR. Kept so that a regression which
        # reintroduces the leak fails here rather than looking like a
        # discovery.
        assert "IS its floor" in flat or "equal to its floor" in flat
        return
    # 2026-11-09: the C4 fix dropped it below the bar. The document must
    # record what it used to read, because "an interval excluding chance
    # that was exactly an imputation floor" is the finding, not a
    # footnote.
    assert "0.526 [0.503, 0.549]" in flat, (
        "this cell no longer clears chance and the document does not "
        "record that it once did, on a figure equal to its floor")


def test_the_withdrawn_claim_about_missing_coverages_is_recorded():
    """The fields are named `auc_coverage_pos` / `auc_coverage_neg`, and
    three documents said they did not exist. The retraction stays."""
    path = ROOT / "data" / "results" / "temporal_2021_include.json"
    if not path.exists():
        pytest.skip("the temporal split is not built in this checkout")
    m = json.loads(path.read_text())["results"]["SARS-CoV-2"]["COMBINED"]
    assert "auc_coverage_pos" in m and "auc_coverage_neg" in m
    flat = " ".join(DOC.read_text().split())
    assert "The fields were there under a different prefix" in flat
