"""The highest AUC a cell could have scored, given what it reaches.

WHAT THIS IS FOR. An AUC near 0.5 means one of two completely different
things -- the scorer ranked badly, or it could not have ranked well --
and nothing in a results cell distinguished them. Unreached compounds
take the imputed floor, so a scorer reaching few positives cannot
separate the classes however good its ranking is.

Measured on release v0.1, SARS-CoV-2, every channel sits within 0.04 of
its ceiling and three sit ON it: M1 measured 0.5135 against a ceiling of
0.5136, M8 0.506 against 0.506, M7 0.4973 against 0.4973. The channels
are not failing to rank. They are failing to reach.

AND THE EXISTING FLOOR DOES NOT CATCH IT. MIN_REACHED_POSITIVES gates on
the COUNT of positives reached, and the count does not determine the
ceiling -- both coverages do. M2 reaches 11 of 89, clears the floor of
10, and has a ceiling of 0.5125.

The ceiling is REPORTED, not enforced: whether it replaces the count as
the readability gate changes what the evaluation claims, and that
belongs in the prereg.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.labels import auc_ceiling, evaluate_against_negatives


def cell(k: int, j: int, n_pos: int = 89, n_neg: int = 7343) -> dict:
    """A best-case cell: k positives and j negatives reached, ranked
    perfectly. Reached positives outrank everything; reached negatives
    sit above the floor and below them."""
    pos = {f"P{i}" for i in range(n_pos)}
    neg = {f"N{i}" for i in range(n_neg)}
    scores = {f"P{i}": 10.0 for i in range(k)}
    scores.update({f"N{i}": 1.0 for i in range(j)})
    return evaluate_against_negatives(scores, pos, neg)


# ------------------------------------ the formula is the thing it claims
@pytest.mark.parametrize("k,j", [(3, 500), (0, 300), (11, 400), (32, 1200),
                                 (89, 2000), (1, 50), (10, 1000), (45, 0)])
def test_the_ceiling_is_what_a_perfect_ranking_actually_scores(k, j):
    """DERIVED AND THEN MEASURED. The formula is four terms of algebra
    and could be wrong in a way no reading catches, so it is checked
    against the evaluator it describes -- build the best-case score map,
    score it, and the two must agree exactly.
    """
    got = cell(k, j)
    assert got["auc"] == pytest.approx(auc_ceiling(k, j, 89, 7343), abs=1e-12)
    assert got["auc_ceiling"] == pytest.approx(got["auc"], abs=1e-12)


def test_full_coverage_can_reach_one():
    """The sanity bound at the top. A scorer reaching every positive is
    limited by nothing, so a ceiling below 1.0 there would mean the
    formula penalises coverage it should not."""
    assert auc_ceiling(89, 7343, 89, 7343) == pytest.approx(1.0)
    assert auc_ceiling(89, 0, 89, 7343) == pytest.approx(1.0)


def test_reaching_no_positive_but_some_negatives_caps_below_chance():
    """M7 ON THE REAL RELEASE. It reaches 39 negatives and no positives,
    so those negatives necessarily outrank the floor-tied positives and
    the ceiling falls under 0.5. Its measured 0.4973 is not a weak
    result; it is the only value it could have taken."""
    c = auc_ceiling(0, 39, 89, 7343)
    assert c < 0.5, c
    assert c == pytest.approx(cell(0, 39)["auc"], abs=1e-12)


def test_reaching_nothing_at_all_is_exactly_chance():
    """The one case that must sit ON 0.5: every compound tied at the
    floor is a coin, not a bias in either direction."""
    assert auc_ceiling(0, 0, 89, 7343) == pytest.approx(0.5)


def test_the_ceiling_rises_with_positives_reached():
    """Monotone in the quantity the floor is supposed to gate on. If it
    were not, 'reaches more' would not mean 'could show more' and the
    whole argument for a reach floor would be unfounded."""
    seq = [auc_ceiling(k, 800, 89, 7343) for k in range(0, 90, 5)]
    assert seq == sorted(seq)
    assert seq[-1] > seq[0]


def test_reaching_more_negatives_lowers_the_ceiling():
    """The half of it a COUNT of positives cannot see, and the reason
    MIN_REACHED_POSITIVES is set on the wrong quantity: two cells
    reaching the same positives have different ceilings."""
    high = auc_ceiling(11, 100, 89, 7343)
    low = auc_ceiling(11, 5000, 89, 7343)
    assert low < high, (low, high)


def test_a_cell_can_clear_the_reach_floor_and_still_be_incapable():
    """M2 ON THE REAL RELEASE: 11 positives reached against a floor of
    10, and a ceiling of 0.5125. This is the finding that makes the
    ceiling worth recording at all, so it is pinned as a test rather
    than left in a comment.
    """
    from kgav.calibration import MIN_REACHED_POSITIVES

    reached = 11
    assert reached >= MIN_REACHED_POSITIVES, (
        "M2's reach no longer clears the floor, so this test is no longer "
        "about a cell the floor admits")
    c = auc_ceiling(reached, 826, 89, 7343)
    assert c == pytest.approx(0.5125, abs=5e-4), c
    assert c - 0.5 < 0.02, (
        "a cell the reach floor admits has less than 0.02 of headroom; "
        "the floor gates on the count and the count does not determine "
        "what a cell can show")


# ----------------------------------------------- it reaches every cell
def test_every_evaluated_cell_carries_a_ceiling():
    """A field computed and not reported is the silence it was meant to
    end. It rides evaluate_against_negatives, so every results file in
    the project gets it without a caller opting in."""
    got = cell(32, 1200)
    assert "auc_ceiling" in got and "auc_headroom" in got
    assert got["auc_headroom"] == pytest.approx(got["auc_ceiling"] - 0.5)


def test_an_unscorable_cell_reports_no_ceiling_rather_than_a_wrong_one():
    """One class empty: there is no AUC, so there is no ceiling either.
    Returning 0.5 here would read as 'could have shown nothing', which
    is a claim about a measurement that was never made."""
    got = evaluate_against_negatives({"P0": 1.0}, {"P0"}, set())
    assert got["evaluable"] is False
    assert "auc_ceiling" not in got


@pytest.mark.parametrize("n_pos,n_neg", [(0, 10), (10, 0), (0, 0)])
def test_the_helper_does_not_divide_by_zero(n_pos, n_neg):
    assert auc_ceiling(0, 0, n_pos, n_neg) == 0.5


def test_counts_beyond_the_class_are_clamped_not_trusted():
    """reached > total is a caller bug, and a ceiling above 1.0 would
    propagate it into a results file as a plausible-looking number."""
    assert auc_ceiling(500, 99_999, 89, 7343) == pytest.approx(1.0)
    assert 0.0 <= auc_ceiling(-5, -5, 89, 7343) <= 1.0


def test_without_imputation_there_is_no_ceiling_to_impose():
    """THE SURVIVING MUTANT THAT WAS RIGHT. The ceiling exists because
    unreached compounds take the floor and drag the AUC toward 0.5. With
    impute_missing=False they are dropped from the evaluation instead,
    so nothing holds the scorer down and a perfect ranking over what it
    reaches scores 1.0.

    The first version computed the ceiling against the full class sizes,
    which imposed an imputation ceiling on a cell with no imputation.
    Every test passed, because every other test uses the imputed path
    where the two are identical.
    """
    pos = {f"P{i}" for i in range(89)}
    neg = {f"N{i}" for i in range(7343)}
    scores = {f"P{i}": 10.0 for i in range(3)}
    scores.update({f"N{i}": 1.0 for i in range(49)})
    got = evaluate_against_negatives(scores, pos, neg, impute_missing=False)
    assert got["evaluable"] is True
    assert got["n_pos"] == 3 and got["n_neg"] == 49, (
        "the non-imputed path should rank only what the scorer reached")
    assert got["auc"] == pytest.approx(1.0), "the fixture ranks perfectly"
    assert got["auc_ceiling"] == pytest.approx(1.0), (
        "a cell with no imputation has no ceiling below 1.0, and this one "
        "is reporting one")


def test_with_imputation_the_same_cell_is_capped_near_chance():
    """The contrast, on identical input. The ceiling is a property of
    the imputation, not of the scores, and a test that only ever ran one
    path could not show that."""
    pos = {f"P{i}" for i in range(89)}
    neg = {f"N{i}" for i in range(7343)}
    scores = {f"P{i}": 10.0 for i in range(3)}
    scores.update({f"N{i}": 1.0 for i in range(49)})
    got = evaluate_against_negatives(scores, pos, neg, impute_missing=True)
    assert got["auc_ceiling"] == pytest.approx(0.5136, abs=5e-4)
