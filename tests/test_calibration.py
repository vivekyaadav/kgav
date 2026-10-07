"""A channel has to earn its place in a ranking.

The measured position on this release is that none of the eight does. These
tests fix the behaviour that follows from that: refuse to rank, say why, and
never let a point AUC above 0.5 stand in for discrimination.
"""
from __future__ import annotations

import json

import pytest

from kgav import calibration as C
from kgav.baselines import combine
from kgav.labels import auc_interval


def _row(auc, lo, hi, n_pos=89, n_neg=7343, cov_pos=0.1, cov_neg=0.1):
    return {"auc": auc, "auc_ci_lo": lo, "auc_ci_hi": hi, "n_pos": n_pos,
            "n_neg": n_neg, "evaluable": True,
            "coverage_pos": cov_pos, "coverage_neg": cov_neg}


# The real single-screen numbers, with their Hanley-McNeil intervals.
MEASURED = {
    "M1": _row(0.513, 0.452, 0.574), "M2": _row(0.492, 0.432, 0.552),
    "M3": _row(0.523, 0.462, 0.584), "M4": _row(0.507, 0.446, 0.568),
    "M5": _row(0.516, 0.455, 0.577), "M6": _row(0.502, 0.442, 0.562),
    "M7": _row(0.499, 0.439, 0.559), "COMBINED": _row(0.516, 0.455, 0.577),
    "degree": _row(0.582, 0.520, 0.644, cov_pos=1.0, cov_neg=1.0),
}
V = "NCBITaxon:2697049"


def _results(tmp_path, rows, args=None, name="hard_negatives.json"):
    doc = {"results": {"cross-sectional": {V: rows}},
           "args": args if args is not None else {
               "label_source": "infores:ncats-opendata",
               "selectivity": "all", "publication_disjoint": True}}
    p = tmp_path / name
    p.write_text(json.dumps(doc))
    return p


# ------------------------------------------------------------- the interval
def test_the_interval_is_what_decides_not_the_point_estimate():
    """0.523 looks like signal and is not. Half of pure noise posts a point
    AUC above 0.5; the interval is the only thing that separates them."""
    good = C.ChannelCalibration("X", 0.523, 0.462, 0.584, 89, 7343)
    assert good.auc > C.CHANCE
    assert not good.discriminates
    assert good.weight == 0.0


def test_a_channel_whose_interval_clears_chance_discriminates():
    c = C.ChannelCalibration("degree", 0.582, 0.520, 0.644, 89, 7343)
    assert c.discriminates
    assert c.weight == pytest.approx(0.082)


def test_weight_is_the_margin_over_chance_not_the_auc():
    """A channel at 0.52 is worth a fifth of one at 0.60, not 87% of it."""
    a = C.ChannelCalibration("a", 0.52, 0.51, 0.53, 500, 500)
    b = C.ChannelCalibration("b", 0.60, 0.58, 0.62, 500, 500)
    assert b.weight / a.weight == pytest.approx(5.0)


def test_an_unevaluable_channel_has_no_weight():
    c = C.ChannelCalibration("M8", None, None, None, 0, 0)
    assert not c.discriminates and c.weight == 0.0
    assert "not evaluable" in c.describe()


@pytest.mark.parametrize("auc,n_pos,n_neg", [
    (0.523, 89, 7343), (0.726, 1730, 4874), (0.99, 5, 5), (0.01, 5, 5),
])
def test_the_interval_stays_inside_zero_to_one(auc, n_pos, n_neg):
    """The normal approximation runs outside [0,1] at small n, and an AUC of
    1.04 shown to a reader is worse than a wide interval."""
    lo, hi = auc_interval(auc, n_pos, n_neg)
    assert 0.0 <= lo <= auc <= hi <= 1.0


def test_a_larger_sample_narrows_the_interval():
    narrow = auc_interval(0.52, 2000, 2000)
    wide = auc_interval(0.52, 50, 50)
    assert (narrow[1] - narrow[0]) < (wide[1] - wide[0])


# --------------------------------------------------------------- loading
def test_load_reads_channels_and_excludes_combined_and_degree(tmp_path):
    cal = C.load(_results(tmp_path, MEASURED), V)
    assert set(cal.channels) == {"M1", "M2", "M3", "M4", "M5", "M6", "M7"}
    assert "COMBINED" not in cal.channels      # circular
    assert "degree" not in cal.channels        # the null hypothesis


def test_the_protocol_is_recorded_from_the_run_not_assumed(tmp_path):
    cal = C.load(_results(tmp_path, MEASURED), V)
    assert "infores:ncats-opendata" in cal.protocol
    assert "publication-disjoint" in cal.protocol


def test_a_results_file_with_no_recorded_args_is_protocol_unknown(tmp_path):
    """An unnamed protocol is what REFUSED_PROTOCOLS exists to catch, so it
    must not be guessed into something that looks trustworthy."""
    cal = C.load(_results(tmp_path, MEASURED, args={}), V)
    assert cal.protocol.endswith("unknown")


def test_the_confounded_protocol_is_refused_by_name(tmp_path):
    """The run that produced 0.823 filtered positives and left negatives
    whole. Its numbers must never become ranking weights."""
    rows = dict(MEASURED, M1=_row(0.823, 0.805, 0.841, 816, 4874))
    cal = C.load(_results(tmp_path, rows, args={"selectivity": "selective-only",
                                                "selectivity_applies": "positives"}), V)
    assert not cal.usable
    assert cal.weights() == {} and cal.discriminating() == []
    assert "17x coverage asymmetry" in cal.why_not()


# -------------------------------------------------- the measured position
def test_no_channel_on_this_release_earns_a_place(tmp_path):
    """THE RESULT. Every metapath interval spans 0.5, so a calibrated ranking
    has nothing to rank with."""
    cal = C.load(_results(tmp_path, MEASURED), V)
    assert cal.discriminating() == []
    assert not cal.can_rank()


def test_why_not_names_the_strongest_channel_and_its_interval(tmp_path):
    cal = C.load(_results(tmp_path, MEASURED), V)
    why = cal.why_not()
    assert "M3" in why and "0.523" in why
    assert "no reasoning channel beats chance" in why
    # and it says what the graph CAN still do
    assert "report evidence for a named compound" in why


def test_the_report_is_readable_and_states_every_channel(tmp_path):
    lines = C.load(_results(tmp_path, MEASURED), V).report()
    joined = "\n".join(lines)
    for m in ("M1", "M2", "M3", "M4", "M5", "M6", "M7"):
        assert m in joined
    assert "indistinguishable from chance" in joined


# ------------------------------------------------- combine, calibrated
def _dwpc(channels: dict[str, dict[str, float]]):
    return {name: {d: {V: v} for d, v in scores.items()}
            for name, scores in channels.items()}


def test_uncalibrated_combine_is_unchanged(tmp_path):
    """Every published number was computed this way and must reproduce."""
    d = _dwpc({"M1": {"a": 1.0, "b": 2.0}, "M7": {"a": 5.0, "c": 1.0}})
    assert combine(d, V) == combine(d, V, calibration=None)
    assert set(combine(d, V)) == {"a", "b", "c"}


def test_calibrated_combine_refuses_when_nothing_discriminates(tmp_path):
    """An order built from eight chance-level channels is not a result."""
    cal = C.load(_results(tmp_path, MEASURED), V)
    d = _dwpc({"M1": {"a": 1.0, "b": 2.0}, "M7": {"a": 5.0, "c": 1.0}})
    assert combine(d, V, calibration=cal) == {}
    assert combine(d, V) != {}            # and the uncalibrated path still ranks


def test_a_noise_channel_cannot_win_once_calibration_applies(tmp_path):
    """M7 is at 0.499. Under max-over-channels a drug at the top of its pool
    scores 1.0 and outranks a drug with a real route; dropping the channel is
    the fix, and this is the behaviour that was wrong."""
    rows = dict(MEASURED, M1=_row(0.70, 0.65, 0.75))     # M1 earns its place
    cal = C.load(_results(tmp_path, rows), V)
    assert cal.discriminating() == ["M1"]

    d = _dwpc({"M1": {"real": 2.0, "other": 1.0},
               "M7": {"noise_king": 9.0, "other": 1.0}})
    uncal = combine(d, V)
    calibrated = combine(d, V, calibration=cal)
    assert uncal["noise_king"] == 1.0          # tops the noise pool, wins
    assert uncal["noise_king"] >= uncal["real"]
    assert "noise_king" not in calibrated      # gone
    assert calibrated["real"] > calibrated.get("other", 0.0)


def test_weights_scale_the_percentile_so_a_weak_channel_counts_less(tmp_path):
    rows = dict(MEASURED, M1=_row(0.70, 0.65, 0.75), M4=_row(0.55, 0.52, 0.58))
    cal = C.load(_results(tmp_path, rows), V)
    assert set(cal.discriminating()) == {"M1", "M4"}
    d = _dwpc({"M1": {"x": 1.0, "top_of_m1": 2.0},
               "M4": {"x": 1.0, "top_of_m4": 2.0}})
    out = combine(d, V, calibration=cal)
    # both sit at percentile 1.0 in their own pool, so only the weight differs
    assert out["top_of_m1"] > out["top_of_m4"]
    assert out["top_of_m1"] == pytest.approx(0.20)   # 0.70 - 0.5
    assert out["top_of_m4"] == pytest.approx(0.05)   # 0.55 - 0.5
