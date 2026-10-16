"""Gate for the registered partial-degree measurement.

Two things must not be wrong. The arithmetic -- a ridge logistic fit and a
chi-square tail, both hand-rolled because pyproject declares pyyaml, pandas
and pyarrow -- and the POWER: a measurement whose registered prediction is a
null has to be shown capable of rejecting that null when signal is present,
or "not rejected" means only that the instrument is blunt.
"""
from __future__ import annotations

import json
import math
import random
import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

P = import_module("partial_degree")


# ---------------------------------------------------------------- arithmetic
def test_the_fit_matches_sklearn_on_a_pinned_dataset():
    """COMPUTED, not recalled. sklearn is NOT a dependency, which is why the
    reference is pinned here rather than imported; the agreement was also
    checked over 200 random datasets of 30 to 300 rows, zero
    non-convergences, max |beta - sklearn| = 5.38e-08."""
    x = [[-1.0, -1.0], [-1.0, 1.0], [1.0, -1.0], [1.0, 1.0], [0.5, -0.5],
         [-0.5, 0.5], [2.0, 0.0], [0.0, 2.0], [-2.0, 0.0], [0.0, -2.0]]
    y = [0, 0, 1, 1, 1, 0, 1, 0, 0, 0]
    beta, ok = P.fit([[1.0, a, b] for a, b in x], y, 1.0)
    assert ok
    for got, want in zip(beta, (-0.609784083927121, 1.4282109116739075,
                                -0.09627971243036944)):
        assert got == pytest.approx(want, abs=1e-8)


def test_the_intercept_is_not_penalised():
    """Penalising it shrinks the fitted base rate toward 0.5, which on an
    imbalanced cell -- every cell here -- biases every other coefficient."""
    y = [1] * 90 + [0] * 10
    x = [[1.0, 0.0]] * 100
    beta, ok = P.fit(x, y, 1.0)
    assert ok
    assert beta[0] == pytest.approx(math.log(90 / 10), abs=1e-6)


@pytest.mark.parametrize("stat,want", [
    (0.0, 1.0),
    (3.841458820694124, 0.05),      # the 95th percentile of chi2(1)
    (6.634896601021213, 0.01),
    (10.827566170662733, 0.001),
])
def test_the_chi_square_tail_is_exact(stat, want):
    """chi2 with one df is a squared standard normal, so the upper tail is
    erfc(sqrt(x/2)) from the standard library -- not a series approximation.
    Checked against scipy.stats.chi2.sf at 400 random points, max error
    4.30e-15."""
    assert P.chi2_sf_1df(stat) == pytest.approx(want, abs=1e-12)


def test_a_singular_system_returns_none_rather_than_a_number():
    assert P.solve([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0]) is None
    assert P.solve([[2.0, 0.0], [0.0, 4.0]], [2.0, 8.0]) == [1.0, 2.0]


def test_the_log_likelihood_does_not_overflow_on_a_separated_cell():
    """Separation is expected on small cells. math.exp(eta) for a large
    positive eta raises OverflowError, which would crash the run rather than
    report NOT FITTED."""
    big = P.loglik([[1.0, 50.0], [1.0, -50.0]], [1, 0], [0.0, 40.0])
    assert math.isfinite(big)
    # A perfectly separated fit assigns probability 1 to every observation,
    # so its log-likelihood is exactly 0. The first version of this test
    # asserted < 0 and failed on correct arithmetic.
    assert big == pytest.approx(0.0, abs=1e-12)
    mixed = P.loglik([[1.0, 1.0], [1.0, -1.0]], [1, 0], [0.0, 0.5])
    assert mixed < 0


def test_a_constant_column_is_not_a_feature():
    """A channel giving every drug the same percentile cannot be told from
    the intercept, and fitting it returns a coefficient whose size is set by
    the penalty alone. M8 is exactly this on three viruses."""
    assert P.standardise([0.5] * 20) is None
    got = P.standardise([0.0, 1.0])
    assert got == pytest.approx([-1.0, 1.0])


# --------------------------------------------------------------------- power
def _cell(n=400, pos_rate=0.25, channel_effect=0.0, degree_effect=1.5,
          seed=0):
    """A labelled pool where the channel's true effect is set by the caller.

    Both features are correlated, as they are in the graph -- the whole
    reason the measurement exists is that they are not independent.
    """
    rng = random.Random(seed)
    chan, deg, y = [], [], []
    for _ in range(n):
        d = rng.gauss(0, 1)
        c = 0.7 * d + 0.7 * rng.gauss(0, 1)       # correlated with degree
        eta = math.log(pos_rate / (1 - pos_rate)) \
            + channel_effect * c + degree_effect * d
        chan.append(c)
        deg.append(d)
        y.append(1 if rng.random() < 1 / (1 + math.exp(-eta)) else 0)
    return chan, deg, y


def test_it_REJECTS_the_null_when_the_channel_really_does_add_signal():
    """THE POWER CHECK. The registration predicts a null. A measurement that
    cannot reject when signal is present makes that prediction unfalsifiable,
    and "not rejected" would mean only that the instrument is blunt."""
    chan, deg, y = _cell(channel_effect=1.2, seed=1)
    cell = P.Cell("v", "M1", "s", "fitted", len(y), sum(y))
    got = P.evaluate(cell, chan, deg, y)
    assert got.status == "fitted"
    assert got.p < 0.001, f"no power: p = {got.p}"
    assert got.b_c > 0


def test_it_does_NOT_reject_when_the_channel_is_only_degree():
    """The complement. With no independent contribution the test must sit at
    chance, or every channel would look significant because both features
    track the same thing."""
    chan, deg, y = _cell(channel_effect=0.0, seed=2)
    cell = P.Cell("v", "M1", "s", "fitted", len(y), sum(y))
    got = P.evaluate(cell, chan, deg, y)
    assert got.p > 0.05, f"false positive: p = {got.p}"


def test_the_false_positive_rate_is_near_alpha_over_many_null_cells():
    """A single null cell passing proves little. Over 200 null cells the
    rejection rate must sit near alpha, or the p-values are miscalibrated
    and every verdict built on them is too."""
    rejected = 0
    for seed in range(200):
        chan, deg, y = _cell(channel_effect=0.0, seed=1000 + seed)
        cell = P.Cell("v", "M1", "s", "fitted", len(y), sum(y))
        got = P.evaluate(cell, chan, deg, y)
        if got.p is not None and got.p < 0.05:
            rejected += 1
    assert rejected <= 20, f"{rejected}/200 null cells rejected at alpha=0.05"


# ---------------------------------------------------- the registered floor
def test_the_floor_is_calibrations_constant_on_both_classes():
    """"10 positives AND 10 negatives". Nine of either is NOT RUN, and NOT
    RUN is not a null."""
    from kgav.calibration import MIN_REACHED_POSITIVES as F
    pool = {f"d{i}": 1.0 for i in range(100)}
    pos = {f"d{i}" for i in range(F)}
    neg = {f"d{i}" for i in range(F, 2 * F)}
    assert P.run_cell("v", "M1", "s", pool, pos, neg).status == "fitted"
    assert P.run_cell("v", "M1", "s", pool, set(list(pos)[:F - 1]),
                      neg).status == "NOT RUN"
    assert P.run_cell("v", "M1", "s", pool, pos,
                      set(list(neg)[:F - 1])).status == "NOT RUN"


def test_the_registered_constants_are_the_documents():
    assert P.LAMBDA == 1.0
    assert P.ALPHA == 0.05
    assert P.PREREG == "docs/prereg_partial_degree.md"


# ------------------------------------------------------------------ the gap
def test_the_undefined_negative_is_declared_a_gap_not_quietly_resolved(
        tmp_path, capsys):
    """The registration says "10 negatives" and never says what a negative
    is -- the same shape of hole prereg_phase4_hcv.md had around "inactive".
    Resolving it silently in the better direction is improving a protocol
    after the fact."""
    P.main(["--release", str(_release(tmp_path))])
    out = capsys.readouterr().out
    assert "GAP:" in out
    assert "never defines a" in out
    assert "--negatives all is the" in out
    src = (ROOT / "scripts" / "partial_degree.py").read_text()
    assert "except where marked GAP" in src


# ------------------------------------------------------------------ end to end
def _release(tmp_path: Path, n_drugs: int = 60, n_proteins: int = 6) -> Path:
    """A release completing M1, with measured-active and measured-inactive
    labels so the floor can be met AND varying reach so the features have
    spread.

    The first version gave every drug one edge to one protein. That is the
    M1 degeneracy in miniature -- identical topology, identical DWPC,
    identical degree -- so every cell came back NO SPREAD and the end-to-end
    test asserted against a run that measured nothing.
    """
    d = tmp_path / "rel"
    d.mkdir(parents=True, exist_ok=True)
    nodes = [{"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
              "properties": {}}]
    edges = []
    for j in range(n_proteins):
        nodes.append({"id": f"UniProtKB:P{j}", "class": "Protein",
                      "properties": {"taxon_id": "NCBITaxon:2697049",
                                     "is_viral": True}})
        nodes.append({"id": f"GENE:{j}", "class": "Gene", "properties": {}})
        edges.append((f"UniProtKB:P{j}", "ENCODED_BY", f"GENE:{j}", {}))
        edges.append((f"GENE:{j}", "BELONGS_TO", "NCBITaxon:2697049", {}))
    q = {"target_resolution": "chain_from_assay_description"}
    for i in range(n_drugs):
        c = f"INCHIKEY:D{i:03d}"
        nodes.append({"id": c, "class": "SmallMolecule", "properties": {}})
        for j in range(1 + i % n_proteins):        # varying reach and degree
            edges.append((c, "INHIBITS", f"UniProtKB:P{j}", q))
        nm = 10.0 if i % 2 == 0 else 100000.0
        edges.append((c, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", "NCBITaxon:2697049",
                      {"ic50_nm": nm}))
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (d / "edges.jsonl").write_text("\n".join(json.dumps(
        {"subject": s, "predicate": p, "object": o, "qualifiers": qu,
         "primary_knowledge_source": "infores:chembl", "evidence_tier": 1,
         "first_asserted_date": "2021-01-01"})
        for s, p, o, qu in edges))
    return d


def test_it_names_the_protocol_and_reports_every_cell(tmp_path, capsys):
    rel = _release(tmp_path)
    code = P.main(["--release", str(rel)])
    out = capsys.readouterr().out
    assert code in (0, 1)
    assert "docs/prereg_partial_degree.md" in out
    assert "Nothing below departs from it except where marked GAP." in out
    assert "lambda = 1.0" in out and "alpha = 0.05" in out
    assert "P1 " in out and "P2 " in out
    assert "NOTHING WAS TESTED" not in out, \
        "the fixture must produce at least one fitted cell"


def test_a_missing_release_exits_two(tmp_path, capsys):
    assert P.main(["--release", str(tmp_path / "nope")]) == 2
    assert "no release at" in capsys.readouterr().out


def test_it_writes_nothing_and_ranks_nothing():
    src = (ROOT / "scripts" / "partial_degree.py").read_text()
    assert "write_results" not in src
    assert "Read-only. Writes nothing, recommends nothing, ranks no drug." \
        in src


def test_a_run_with_no_fitted_cell_does_not_claim_P1_held(tmp_path, capsys):
    """NOT TESTABLE and NOT REJECTED are different states, and the registered
    prediction is a null -- so the empty run is exactly where a report would
    quietly turn "we measured nothing" into "the prediction held"."""
    assert P.verdict([]) == 1
    out = capsys.readouterr().out
    assert "NOTHING WAS TESTED" in out
    assert "not the same as P1 holding" in out


def test_an_untested_cell_is_reported_with_its_reason(capsys):
    """Below the floor and no-spread are different findings about the graph.
    The first version collapsed both into "NO CELL CLEARED THE FLOOR", which
    would have reported the M1 degeneracy as a sample-size problem."""
    P.why_not_tested([
        P.Cell("v", "M1", "s", "NO SPREAD"),
        P.Cell("v", "M2", "s", "NOT RUN"),
        P.Cell("v", "M3", "s", "NOT RUN"),
    ])
    out = capsys.readouterr().out
    assert "NOT RUN" in out and "2" in out
    assert "NO SPREAD" in out
    assert "cannot be told from the intercept" in out
    assert "below the registered floor" in out
