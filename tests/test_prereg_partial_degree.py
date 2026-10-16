"""The registration is binding only if something reads it.

docs/prereg_partial_degree.md fixes six parameters before the code exists.
A registration nobody checks is a description written early, so these pin
the numbers now and bind the implementation to them when it lands -- the
cross-check runs against the script the moment the file appears, and
asserts nothing about it before.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PREREG = ROOT / "docs" / "prereg_partial_degree.md"
SCRIPT = ROOT / "scripts" / "partial_degree.py"


@pytest.fixture(scope="module")
def text() -> str:
    return PREREG.read_text()


def test_it_is_a_preregistration_in_the_house_shape(text):
    """Body first, verdict appended beneath -- as prereg_temporal_auc.md and
    prereg_phase4_hcv.md do. Until the measurement runs there is no verdict,
    and the absence is what says it has not run."""
    heads = [l for l in text.splitlines() if l.startswith("# ")]
    assert heads[0].startswith("# Pre-registration")
    assert "Written **before** the code" in text
    if any(h.startswith("# Verdict") for h in heads):
        assert heads[-1].startswith("# Verdict"), \
            "the verdict goes at the end, beneath the body it judges"


@pytest.mark.parametrize("fixed", [
    "lambda = 1.0",
    "alpha = 0.05",
    "10 positives AND 10 negatives",
    "Bonferroni",
    "likelihood-ratio",
    "midrank percentiles",
])
def test_every_fixed_parameter_is_stated(text, fixed):
    """Each of these decides a number in the output. One left unstated is one
    the implementation gets to choose after seeing the data."""
    assert fixed in text


def test_all_four_predictions_are_registered(text):
    for p in ("**P1.", "**P2.", "**P3.", "**P4."):
        assert p in text, f"{p} missing"
    assert "What would change the conclusion" in text


def test_it_states_what_it_cannot_do(text):
    """A registration that only lists what it will show is an advertisement.
    The events-per-variable limit in particular must be stated BEFORE the
    result, or a null reads as evidence of absence."""
    assert "## What this cannot do" in text
    assert "events-per-variable" in text
    assert "weak evidence of absence" in text


def test_it_does_not_let_the_implementation_copy_the_percentile_arithmetic(text):
    """diagnose_controls.py once reimplemented channel_percentiles and the
    copy went stale, printing a percentile of 0.0182 beside a combined score
    of 0.4041. One implementation, named in the registration."""
    assert "baselines.channel_percentiles" in text
    assert "rather than a second copy" in text


def test_it_names_the_absent_dependencies_it_must_hand_roll(text):
    """pyproject declares pyyaml, pandas, pyarrow. A registration that
    assumed scipy would be rewritten at implementation time, which is the
    thing a registration exists to prevent."""
    assert "no scipy and no" in text
    deps = (ROOT / "pyproject.toml").read_text()
    body = deps.split("dependencies = ", 1)[1].split("]", 1)[0]
    for absent in ("scipy", "sklearn", "scikit-learn", "numpy"):
        assert absent not in body, f"{absent} is now a dependency; " \
            "the registration's hand-rolling constraint is stale"


# ------------------------------------------- binding, once the code lands
@pytest.mark.skipif(not SCRIPT.exists(), reason="implementation not written yet")
def test_the_implementation_uses_the_registered_constants():
    """THE POINT OF THE FILE ABOVE. Runs only once scripts/partial_degree.py
    exists, and then refuses any drift from what was registered."""
    src = SCRIPT.read_text()
    assert re.search(r"\bLAMBDA\s*=\s*1\.0\b", src), "registered lambda = 1.0"
    assert re.search(r"\bALPHA\s*=\s*0\.05\b", src), "registered alpha = 0.05"
    assert "MIN_REACHED_POSITIVES" in src, \
        "the floor is calibration's constant, not a literal 10"
    assert "channel_percentiles" in src, \
        "the registration forbids a second copy of the percentile arithmetic"
    assert "docs/prereg_partial_degree.md" in src, \
        "the script must name the protocol it is executing"


@pytest.mark.skipif(not SCRIPT.exists(), reason="implementation not written yet")
def test_the_implementation_runs_both_scorings():
    """P2 is the falsification test and needs both. A script that ran one
    would make P2 unfalsifiable while appearing to execute the protocol."""
    src = SCRIPT.read_text()
    assert "SCORING_LEGACY" in src and "SCORING_PDP" in src
