"""Gate for the label-headroom planner.

It sizes a future measurement, so the way it goes wrong is by promising a
sensitivity the registered test does not have. Two guards matter: it must
use the REGISTERED statistic rather than a cleaner one, and its power
curve must behave like a power curve -- monotone in n, and at the null
equal to the corrected alpha rather than to anything flattering.
"""
from __future__ import annotations

import json
import random
import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

H = import_module("label_headroom")


def test_headroom_is_reach_minus_what_is_already_labelled():
    pool = {f"d{i}" for i in range(100)}
    labelled = {f"d{i}" for i in range(30)} | {"outside"}
    assert H.headroom(pool, labelled) == (100, 30, 70)
    assert H.headroom(set(), labelled) == (0, 0, 0)


def test_the_statistic_is_imported_not_reimplemented():
    """A second copy would drift from the registered test, and the point of
    the power curve is to describe the test that actually ran -- including
    the defect where its ratio can come out negative."""
    src = (ROOT / "scripts" / "label_headroom.py").read_text()
    assert "from partial_degree import" in src
    assert "def fit(" not in src
    assert "def chi2_sf_1df(" not in src
    assert "negativity defect included" in src


@pytest.mark.parametrize("alpha,want", [
    (0.05, 3.8414588206941245),
    (0.01, 6.634896601021213),
    (0.001, 10.827566170662733),
])
def test_the_critical_value_is_found_not_tabulated(alpha, want):
    assert H.chi2_crit_1df(alpha) == pytest.approx(want, abs=1e-6)


def test_the_normal_cdf_is_the_standard_one():
    assert H.normal_cdf(0.0) == pytest.approx(0.5, abs=1e-12)
    assert H.normal_cdf(1.959963985) == pytest.approx(0.975, abs=1e-9)
    assert H.normal_cdf(-1.959963985) == pytest.approx(0.025, abs=1e-9)


def test_power_rises_with_n_when_an_effect_is_real():
    """The defining property. If it did not, the 'needed' column would be
    noise dressed as a plan."""
    small = H.power_closed_form(200, 0.05, 0.6, 0.5, 1)
    large = H.power_closed_form(3000, 0.05, 0.6, 0.5, 1)
    assert large > small
    assert large > 0.8


def test_power_at_the_null_IS_alpha():
    """With b_c = 0 the non-centrality is 0 and the rejection rate is
    exactly the corrected alpha. Anything else and every 'needed' figure
    is wrong."""
    assert H.power_closed_form(1500, 0.05, 0.0, 0.5, 1) == \
        pytest.approx(0.05, abs=1e-6)
    assert H.power_closed_form(1500, 0.05, 0.0, 0.5, 20) == \
        pytest.approx(0.05 / 20, abs=1e-6)


def test_the_correction_makes_the_test_harder_not_easier():
    """Power is computed against the Bonferroni-corrected alpha, because
    that is what the registered analysis applies."""
    assert H.power_closed_form(800, 0.05, 0.5, 0.5, 26) < \
        H.power_closed_form(800, 0.05, 0.5, 0.5, 1)


def test_collinearity_costs_power():
    """The (1 - rho^2) term. A channel that IS degree carries no
    independent information, and the curve has to say so -- four cells in
    diagnose_scorer.py read rho = -1.000 exactly."""
    assert H.power_closed_form(2000, 0.05, 0.5, 0.0, 1) > \
        H.power_closed_form(2000, 0.05, 0.5, 0.9, 1)
    assert H.power_closed_form(2000, 0.05, 0.5, 1.0, 1) == \
        pytest.approx(0.05, abs=1e-6)


def test_the_closed_form_agrees_with_the_registered_statistic():
    """THE CHECK THAT MATTERS. The curve is analytic for speed; the
    simulation is what uses the test that actually ran.

    The full comparison -- six settings, agreement within 0.071, the
    simulation higher in five of six, both at alpha when b_c = 0 -- was
    run once and is recorded in the script's docstring. This is a cheap
    regression guard on one point, not that comparison: the simulation
    costs seconds per point and `pytest -q` runs everything, including
    what `make test` excludes.

    The `slow` marker is deliberately NOT used. Its documented meaning in
    pyproject is "integration check against data/raw"; this is a
    self-contained simulation, and borrowing the marker would make the
    exclusion list mean two things.
    """
    cf = H.power_closed_form(400, 0.08, 0.7, 0.4, 1)
    sim = H.power_at(400, 0.08, 0.7, 0.3, 0.4, 1, trials=60, seed=5)
    assert abs(cf - sim) < 0.15, f"closed {cf:.3f} vs simulated {sim:.3f}"
    assert H.power_closed_form(1500, 0.05, 0.0, 0.5, 1) == \
        pytest.approx(0.05, abs=1e-6)


def test_features_are_simulated_as_percentiles_not_gaussians():
    """The real pipeline percentiles both columns. A Gaussian-feature power
    curve answers about a different model."""
    src = (ROOT / "scripts" / "label_headroom.py").read_text()
    assert "ranks(" in src
    assert "MIDRANK PERCENTILES" in src


def test_an_undetectable_effect_reports_the_cap_rather_than_a_number():
    """A tiny effect must not come back as a confident label count. Saying
    ">cap" is the useful answer: the thing cannot be chased with labels."""
    need, power = H.labels_needed(0.04, 0.005, 0.9, 26, cap=800)
    assert need is None
    assert power < H.TARGET_POWER


def test_labels_needed_returns_a_size_reaching_the_target():
    need, power = H.labels_needed(0.10, 1.0, 0.3, 1, cap=16000)
    assert need is not None
    assert power >= H.TARGET_POWER
    assert H.power_closed_form(need - 1, 0.10, 1.0, 0.3, 1) < H.TARGET_POWER


def test_labels_needed_is_fast_enough_to_run_in_a_suite():
    """The first version searched by simulation, doubling to 64,000 at 19
    seconds a point. A planning tool nobody waits for is not used."""
    import time
    start = time.time()
    for _ in range(20):
        H.labels_needed(0.037, 0.403, 0.5, 16)
    assert time.time() - start < 5.0


# ------------------------------------------------------------------ framing
def test_it_says_it_is_planning_and_not_evidence(tmp_path, capsys):
    """The failure mode for a planning script is being read as a result."""
    H.main(["--release", str(_release(tmp_path))])
    out = capsys.readouterr().out
    assert "PLANNING ONLY" in out
    assert "is evidence about any channel" in out
    assert "needs its own registration" in out


def test_it_warns_where_the_need_exceeds_the_reach(tmp_path, capsys):
    """If a cell needs more labelled compounds than the channel reaches, no
    panel over this graph can settle it, and the output has to say so
    rather than printing a number that cannot be acted on."""
    H.main(["--release", str(_release(tmp_path))])
    out = capsys.readouterr().out
    assert "no panel over" in out
    assert "more graph, not more labels" in out


def test_it_writes_nothing():
    src = (ROOT / "scripts" / "label_headroom.py").read_text()
    assert "write_results" not in src
    assert "Read-only. Writes nothing." in src


def test_a_missing_release_exits_two(tmp_path, capsys):
    assert H.main(["--release", str(tmp_path / "nope")]) == 2
    assert "no release at" in capsys.readouterr().out


def _release(tmp_path: Path) -> Path:
    sys.path.insert(0, str(ROOT / "tests"))
    return import_module("test_partial_degree")._release(tmp_path)
