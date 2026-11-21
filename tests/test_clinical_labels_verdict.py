"""Every figure in the clinical-label verdict is the run's.

The verdict in doc2/registrations/prereg_clinical_labels.md turns on three numbers --
0.4087, 0.7751 and the direction of the first -- and on four control
counters that are the only evidence the exclusion behaved as registered.
All of them are derived here from data/results/clinical_labels.json, so
a re-run that moves any of them fails rather than leaving the prose
standing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

PREREG = ROOT / "doc2" / "registrations" / "prereg_clinical_labels.md"
RESULTS = ROOT / "data" / "results" / "clinical_labels.json"


@pytest.fixture(scope="module")
def run() -> dict:
    if not RESULTS.exists():
        pytest.skip("the clinical-label protocol has not been run")
    return json.loads(RESULTS.read_text())


@pytest.fixture(scope="module")
def verdict() -> str:
    text = PREREG.read_text()
    assert "# Verdict, 2026-11-17" in text, "the verdict section is gone"
    return text.split("# Verdict, 2026-11-17")[1]


def test_the_two_aucs_and_their_intervals_are_the_runs(run, verdict):
    for key in ("structure_max_k10", "degree"):
        m = run[key]
        auc = f"{m['auc']:.4f}"
        lo = f"{m['auc_ci_lo_corrected']:.4f}"
        hi = f"{m['auc_ci_hi_corrected']:.4f}"
        assert auc in verdict, f"{key} AUC {auc} is not in the verdict"
        assert f"[{lo}, {hi}]" in verdict, (
            f"{key}'s interval [{lo}, {hi}] is not in the verdict")


def test_p1_is_reported_as_below_chance_while_it_is_below_chance(run, verdict):
    """THE DIRECTION IS THE WHOLE RESULT.

    The registered prediction is satisfied by an interval that excludes
    0.5, and this one excludes it from BELOW. A verdict that reported
    "P1 TRUE" without the direction would be literally correct and
    completely misleading, so the direction is asserted against the
    interval rather than trusted to the prose.
    """
    m = run["structure_max_k10"]
    assert m["clears_corrected_alpha"] is True
    below = m["auc_ci_hi_corrected"] < 0.5
    assert below, (
        "the structure statistic no longer sits wholly below chance; the "
        "verdict's 'TRUE, AND IN THE WRONG DIRECTION' is then wrong")
    assert "WRONG DIRECTION" in verdict
    assert "anti-ranking" in verdict


def test_p2_is_decided_by_the_margin_the_run_measured(run, verdict):
    gap = run["degree"]["auc"] - run["structure_max_k10"]["auc"]
    assert run["p2_degree_at_or_above_structure"] is True, (
        "degree no longer beats the structure statistic; the verdict's "
        "fame reading has to be rewritten, not reinterpreted")
    assert f"{gap:.2f} AUC" in verdict, (
        f"the verdict does not state the measured margin {gap:.2f}")


def test_the_class_sizes_and_their_in_vitro_split_are_the_runs(run, verdict):
    pop = run["population"]
    assert f"{pop['positives']} positives" in verdict
    assert f"{pop['negatives']:,} negatives" in verdict
    flat = " ".join(verdict.split())
    assert (f"**{pop['positives_invitro_inactive']} of the "
            f"{pop['positives']} positives were measured INACTIVE") in flat, (
        "the verdict does not state the inactive share of the positive "
        "class with its denominator")
    assert pop["positives_invitro_inactive"] > pop["positives_invitro_active"], (
        "the positive class is no longer mostly in-vitro inactive, which is "
        "the premise the verdict explains the anti-ranking with")


def test_the_symmetric_exclusion_invariant_holds_and_is_quoted(run, verdict):
    """P2 of the symmetric registration, on this cell.

    If every candidate did not lose exactly k references, the control is
    not symmetric and no comparison in this cell means anything.
    """
    c = run["structure_max_k10"]["counts"]
    k = 10
    assert c["excluded_knn"] == k * c["candidates"], (
        f"excluded_knn {c['excluded_knn']} against {k} x "
        f"{c['candidates']}; the exclusion is not symmetric")
    assert c["knn_short"] == 0
    assert f"{c['excluded_knn']:,} = {k} x {c['candidates']:,}" in verdict
    assert f"`excluded_self` = {c['references']:,}" in verdict


def test_nothing_was_imputed_and_the_verdict_says_so(run, verdict):
    c = run["structure_max_k10"]["counts"]
    fp = run["fingerprinted"]
    assert c["empty_reference"] == 0 and c["no_fingerprint"] == 0
    assert fp.get("unparseable", 0) == 0 and fp.get("no_smiles", 0) == 0
    assert f"{fp['fingerprinted']:,} of {fp['fingerprinted']:,}" in verdict
    assert c["excluded_publication"] == 0, (
        "the publication exclusion fired; hazard H3 registered it as off")


def test_p3_is_reported_as_an_observation_with_its_medians(run, verdict):
    p3 = run["p3_observation"]
    assert f"{p3['median_score']:.4f}" in verdict
    assert f"{p3['median_score_negatives']:.4f}" in verdict
    assert "reported not tested" in verdict or "reported, not tested" in verdict


def test_the_registered_operand_gap_is_named_rather_than_absorbed(run, verdict):
    """222 registered, 221 measured. A protocol that quietly reported the
    number it got is the failure mode a pre-registration exists to stop."""
    pop = run["population"]
    registered = 222
    assert pop["positives"] != registered, (
        "the cell now holds the registered 222 positives, so this test is "
        "about a gap that no longer exists and should be removed")
    flat_doc = " ".join(PREREG.read_text().split())
    assert (f"**{pop['positives']} positives against {pop['negatives']:,} "
            f"negatives, where this document registered {registered}") in \
        flat_doc, "the verdict does not name the operand gap"
    assert "no ClinicalTrials reference" in flat_doc


def test_the_run_is_the_registered_protocol_and_not_a_variant(run):
    assert run["family_size"] == 1
    assert run["alpha"] == 0.05
    assert run["z"] == 1.96
    assert run["disease"] == "MONDO:0100096"
    assert run["virus"] == "NCBITaxon:2697049"
    assert "k=10" in run["statistic"]
    assert "no publication exclusion" in run["statistic"]


def test_the_verdict_refuses_the_predictive_claim(verdict):
    """Hazard H2 required this sentence, in those words."""
    flat = " ".join(verdict.split())
    assert "no temporal split" in flat
    assert "does not license a predictive claim" in flat or \
        "would not have licensed a predictive claim" in flat
