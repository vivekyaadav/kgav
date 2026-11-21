"""The prospective verdict is derived from its results file, not typed.

WHY THIS FILE EXISTS SEPARATELY from the operand tests. Those check a
registration that had not run. This checks one that has, and the thing
most likely to rot is the verdict paragraph: a figure quoted in prose
beside a table that a rebuild can move. Every number asserted below is
read out of data/results/temporal_structure_scorer.json or computed
from it.

THE DECISION RULE IS ALSO TESTED, not just the numbers. The registration
fixed three outcomes in advance and said which one each result shape
produces. A verdict that announced the wrong one while quoting the right
figures would pass a figures-only check.
"""
from __future__ import annotations

import json
import re
import statistics
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DOC = ROOT / "doc2" / "registrations" / "prereg_temporal_structure_scorer.md"
RESULTS = ROOT / "data" / "results" / "temporal_structure_scorer.json"
CROSS = ROOT / "data" / "results" / "structure_scorer_strength_matched.json"


@pytest.fixture(scope="module")
def run() -> dict:
    if not RESULTS.exists():
        pytest.skip("the prospective study has not run on this checkout")
    return json.loads(RESULTS.read_text())


@pytest.fixture(scope="module")
def verdict() -> str:
    doc = DOC.read_text()
    marker = "# Verdict, 2026-11-19"
    assert marker in doc, "the registration carries no verdict"
    return doc[doc.index(marker):]


def test_every_cell_in_the_verdict_table_matches_the_results_file(run,
                                                                  verdict):
    """All ten rows, both directions: a table that omits the cells which
    went backwards is the selection failure this project is built
    against."""
    rows = re.findall(
        r"^    (\S[^\n]*?)\s+([\d,]+)\s+([\d,]+)\s+([\d,]+)\s+"
        r"(\d\.\d{4})\s+\[(\d\.\d{4}), (\d\.\d{4})\]", verdict, re.M)
    assert len(rows) == len(run["cells"]), (
        f"{len(rows)} rows in the verdict, {len(run['cells'])} cells "
        f"measured")
    byname = {c["name"]: c for c in run["cells"].values()}
    for name, refs, pos, neg, auc, lo, hi in rows:
        c = byname[name.strip()]
        assert int(refs.replace(",", "")) == c["n_references"]
        assert int(pos.replace(",", "")) == c["n_test_positives"]
        assert int(neg.replace(",", "")) == c["n_test_negatives"]
        p = c["primary"]
        assert float(auc) == round(p["auc"], 4)
        assert float(lo) == round(p["auc_ci_lo_corrected"], 4)
        assert float(hi) == round(p["auc_ci_hi_corrected"], 4)


def test_p1_holds_and_the_count_is_derived(run, verdict):
    above = [c for c in run["cells"].values()
             if c["primary"]["clears_corrected_alpha"]
             and c["primary"]["direction"] == "above"]
    assert above, "P1 failed; the verdict says it holds"
    low = verdict.lower()
    assert (f"{_word(len(above))} cells clear" in low
            or f"{len(above)} cells clear" in low
            or f"{len(above)} of {len(run['cells'])} cell" in low), (
        f"{len(above)} cells clear above chance and the verdict does not "
        f"say so")


def test_p2_is_the_derived_comparison_not_the_predicted_one(run, verdict):
    """The registration predicted against 0.6756, which was typed. The
    verdict must quote the DERIVED cross-sectional median and say the
    earlier figure was wrong rather than quietly using the right one."""
    med = statistics.median(c["primary"]["auc"]
                            for c in run["cells"].values())
    assert f"{med:.4f}" in verdict, f"temporal median {med:.4f} not quoted"
    if CROSS.exists():
        x = json.loads(CROSS.read_text()).get("primary_chembl") or {}
        xs = [x[k]["structure_max"]["auc"] for k in run["cells"] if k in x]
        xm = statistics.median(xs)
        assert f"{xm:.4f}" in verdict, f"{xm:.4f} not quoted"
        assert med < xm, "P2 failed; the verdict says it holds"
        assert "0.6756" in verdict, (
            "the registration's typed figure was wrong and the verdict "
            "must say so, not silently replace it")


def test_p3_is_reported_as_the_count_of_usable_comparators(run, verdict):
    usable = run["comparator_usable"]
    assert len(usable) == 0, usable
    assert "ZERO of ten" in verdict or "0 of 10" in verdict
    reaches = [c["degree"]["reach"] for c in run["cells"].values()]
    assert f"{min(reaches):.1%}" in verdict
    assert f"{max(reaches):.1%}" in verdict


def test_the_decision_rule_announced_matches_the_result_shape(run, verdict):
    """The registration fixed three outcomes. Exactly the one the data
    produces may be announced."""
    clears = run["clears"]
    beat = run["clears_and_beats_usable_comparator"]
    if beat:
        want, forbid = "prospective result", "MEASURABLE, NOT COMPARABLE"
    elif clears:
        want, forbid = "MEASURABLE, NOT COMPARABLE", "P1 FAILS"
    else:
        want, forbid = "P1 FAILS", "MEASURABLE, NOT COMPARABLE"
    assert want in verdict, f"the verdict does not announce '{want}'"
    assert forbid not in verdict.replace(want, ""), (
        f"the verdict also announces '{forbid}'")


def test_the_cells_that_rank_backwards_are_named(run, verdict):
    """Two cells clear the alpha BELOW chance. A verdict that reported
    only the six above would be true and misleading."""
    below = [c["name"] for c in run["cells"].values()
             if c["primary"]["clears_corrected_alpha"]
             and c["primary"]["direction"] == "below"]
    for name in below:
        assert name in verdict, f"{name} ranks backwards and is not named"
    assert below, "no cell ranks backwards; the verdict claims two do"


def test_the_graph_reach_claim_is_the_measured_one(run, verdict):
    """100% fingerprint coverage against the graph's few percent is the
    verdict's central contrast and must come from the run."""
    fp = run["fingerprint_stats"]
    n = fp["fingerprinted"]
    assert f"{n:,} of {n:,}" in verdict


def test_no_cell_reports_a_comparator_auc_of_zero(run):
    """A cell whose comparator could not be evaluated must carry None,
    not 0.0, which would read as a measured failure of the baseline."""
    for t, c in run["cells"].items():
        d = c["degree"]
        if not d.get("evaluable", True):
            assert d.get("auc") is None, (t, d.get("auc"))


def _word(n: int) -> str:
    return {0: "zero", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
            6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}.get(n,
                                                                        str(n))
