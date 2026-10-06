"""Gate for the hard-negative driver: the two label sides are filtered by one
function, and the metapath reach survives into the results.
"""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _mod():
    spec = importlib.util.spec_from_file_location(
        "hard_negatives", ROOT / "scripts" / "hard_negatives.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class _Labels:
    def __init__(self, pos, neg):
        self.positives, self.negatives = pos, neg


V = "NCBITaxon:2697049"


def _labels():
    return _Labels({V: {"A", "B", "C", "D"}}, {V: {"W", "X", "Y", "Z"}})


# si: A and W are selective, B and X have an index but are cytotoxic,
# C, D, Y, Z have no paired cytotoxicity data at all.
SI = {("A", V): 50.0, ("B", V): 2.0, ("W", V): 80.0, ("X", V): 1.5}


def test_default_filters_positives_only_which_is_not_like_for_like():
    """The published 0.826 is computed this way, and the asymmetry is real.

    Requiring a verified index of a POSITIVE requires paired same-document
    CC50/EC50, which selects for well-characterised compounds -- the ones most
    likely to carry the target annotations M1 reads. Measured on v0.1 the
    positive set halves 1,726 -> 808 and M1 rises 0.729 -> 0.823, while
    cov_pos goes 51.5% -> 70.3% and cov_neg stays at 5.3%.
    """
    m = _mod()
    lab = _labels()
    args = SimpleNamespace(selectivity="selective-only", si_threshold=10.0,
                           selectivity_applies="positives")
    m._filter_labels(lab, SI, args)
    assert lab.positives[V] == {"A"}          # only the selective active
    assert lab.negatives[V] == {"W", "X", "Y", "Z"}   # untouched


def test_both_applies_the_same_requirement_to_negatives():
    """The control. A compound must clear the same evidence bar to be a
    negative as to be a positive."""
    m = _mod()
    lab = _labels()
    args = SimpleNamespace(selectivity="selective-only", si_threshold=10.0,
                           selectivity_applies="both")
    m._filter_labels(lab, SI, args)
    assert lab.positives[V] == {"A"}
    assert lab.negatives[V] == {"W"}


def test_verified_only_keeps_cytotoxic_compounds_on_both_sides():
    """verified-only isolates 'has matched data' from 'passed the threshold',
    so a cytotoxic compound with an index stays in."""
    m = _mod()
    lab = _labels()
    args = SimpleNamespace(selectivity="verified-only", si_threshold=10.0,
                           selectivity_applies="both")
    m._filter_labels(lab, SI, args)
    assert lab.positives[V] == {"A", "B"}
    assert lab.negatives[V] == {"W", "X"}


def test_one_helper_serves_both_the_cross_sectional_and_temporal_sides():
    """The temporal block held a second hand-written copy of the loop, which
    is how the two sides came to differ without anyone deciding they should."""
    src = (ROOT / "scripts" / "hard_negatives.py").read_text()
    assert src.count("_filter_labels(") >= 3      # definition + both call sites
    assert "si_t.get(" not in src                 # the duplicated loop is gone


def test_cross_sectional_reach_is_not_discarded_before_it_is_written():
    """Batch 2a reset cross_reach to {} one line AFTER run() returned it, so
    the cross-sectional metapath reach never reached the results file -- the
    exact thing 2a was written to record."""
    src = (ROOT / "scripts" / "hard_negatives.py").read_text()
    assigned = src.index("cross, cross_reach = run(")
    reinit = src.find("cross_reach: dict = {}")
    assert reinit == -1 or reinit < assigned, \
        "cross_reach is re-initialised after run() set it"
