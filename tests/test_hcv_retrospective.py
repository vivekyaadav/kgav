"""The registered test, and the ways a registered test gets quietly bent.

The prereg is the specification. So the failures worth catching are not
arithmetic ones -- they are the protocol drifting from what was written:
a threshold restated instead of deferred to, a case the wording did not
consider silently resolved the better way, or a power table in the document
that the shipped code does not actually reproduce.
"""
from __future__ import annotations

import re
import sqlite3
import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PREREG = ROOT / "docs" / "prereg_phase4_hcv.md"


def _mod(name):
    sys.path.insert(0, str(ROOT / "scripts"))
    return import_module(name)


R = _mod("phase4_hcv_retrospective")


def _fixture(tmp_path: Path, **kw) -> Path:
    sys.path.insert(0, str(ROOT / "tests"))
    return import_module("test_phase4_survey")._db(
        tmp_path / "retro.db", **kw)


# ------------------------------------------------------- the exact test
@pytest.mark.parametrize("table,greater,expected", [
    # COMPUTED, not recalled. The first draft of this list carried five
    # values written from memory and four of them were wrong -- which is
    # the whole failure this repository spent a session removing, in the
    # test file for a pre-registered measurement.
    #
    # Each value below is scipy.stats.fisher_exact's answer, printed and
    # pasted. scipy is NOT a dependency (pyproject: pyyaml, pandas,
    # pyarrow), which is why it is pinned here rather than imported: the
    # agreement was also checked over 400 random tables, both tails, zero
    # mismatches.
    ((3, 9, 1, 11), True, 0.2950310559006211),
    ((3, 9, 1, 11), False, 0.953416149068323),
    ((8, 2, 1, 9), True, 0.002738747320790664),
    ((8, 2, 1, 9), False, 0.9999404620147654),
    ((0, 12, 6, 6), True, 1.0),
    ((0, 12, 6, 6), False, 0.006864988558352404),
])
def test_fisher_matches_the_exact_distribution(table, greater, expected):
    a, b, c, d = table
    assert R.fisher_one_sided(a, b, c, d, greater=greater) == \
        pytest.approx(expected, abs=1e-12)


def test_a_degenerate_table_is_one_not_a_crash():
    """An empty row or an all-active column has no tail to sum."""
    assert R.fisher_one_sided(0, 0, 5, 5) == 1.0
    assert R.fisher_one_sided(5, 0, 5, 0) == 1.0
    assert R.fisher_one_sided(0, 5, 0, 5) == 1.0


def test_the_shipped_test_reproduces_the_preregs_power_table():
    """THE CROSS-CHECK THAT MATTERS. The prereg's table was computed from a
    binomial before this script existed; the shipped code uses Fisher
    against a large background. If they disagree, either the document
    promises a sensitivity the code does not deliver or the code is not the
    test that was registered."""
    rows = re.findall(r"^\s{10,}(\d+)%\s+(\d+)\s+\d+%\s+\d+\s+\d+%\s*$",
                      PREREG.read_text(), re.M)
    assert len(rows) == 6, "the registered power table is not 6 rows"
    big = 10_000
    for pct, needed in rows:
        p = int(pct) / 100
        actives = round(p * big)
        got = next(k for k in range(13)
                   if R.fisher_one_sided(k, 12 - k, actives, big - actives,
                                         greater=True) < R.ALPHA)
        assert got == int(needed), f"background {pct}%: doc {needed}, code {got}"


# --------------------------------------------- the rule is not restated
def test_the_threshold_lives_in_labels_not_in_this_script():
    """The prereg: "The implementation CALLS that function; it does not
    restate the rule." A 10,000 written here would be a second definition
    of activity on the measurement the whole verdict turns on."""
    src = (ROOT / "scripts" / "phase4_hcv_retrospective.py").read_text()
    body = src.split('"""', 2)[2]
    assert "10_000" not in body and "10000" not in body
    assert "from kgav.labels import classify" in src


def test_the_edge_is_what_classify_reads():
    from kgav.labels import DEFAULT_ACTIVE_BELOW_NM, classify
    potent = {"standard_type": "IC50", "standard_value": 50.0,
              "standard_units": "nM", "standard_relation": "="}
    weak = {**potent, "standard_value": DEFAULT_ACTIVE_BELOW_NM * 10}
    assert classify(R.coronavirus_edge(potent)) == "active"
    assert classify(R.coronavirus_edge(weak)) == "inactive"


def test_a_row_classify_cannot_use_is_none():
    assert R.coronavirus_edge({"standard_type": "Inhibition",
                               "standard_value": 1.0,
                               "standard_units": "%",
                               "standard_relation": "="}) is None
    assert R.coronavirus_edge({"standard_type": "IC50",
                               "standard_value": 1.0,
                               "standard_units": "ug.mL-1",
                               "standard_relation": "="}) is None


def test_censored_rows_are_classifys_business_not_ours():
    """A ">" past the threshold proves inactivity and a "<" under it proves
    activity; both rules are classify's, and this script passes the
    relation through rather than deciding."""
    from kgav.labels import classify
    over = {"standard_type": "IC50", "standard_value": 100000.0,
            "standard_units": "nM", "standard_relation": ">"}
    under = {"standard_type": "IC50", "standard_value": 10.0,
             "standard_units": "nM", "standard_relation": "<"}
    vague = {"standard_type": "IC50", "standard_value": 10.0,
             "standard_units": "nM", "standard_relation": ">"}
    assert classify(R.coronavirus_edge(over)) == "inactive"
    assert classify(R.coronavirus_edge(under)) == "active"
    assert classify(R.coronavirus_edge(vague)) is None


# ------------------------- the case the registered wording did not cover
def test_the_literal_wording_is_the_primary_and_the_gap_is_shown():
    """The prereg defines "active" and not "inactive". Read literally, a
    compound whose every row was undecidable is not-active -- and an
    undecidable row is not evidence of inactivity. The literal reading is
    PRIMARY because it is what was registered; resolving it the better way
    and saying nothing would be improving the protocol after the fact."""
    verdicts = {"a": {"active"}, "b": {"inactive"}, "c": {None},
                "d": {None}, "e": {"active"}}
    hcv = {"a", "c"}
    lit = R.two_by_two(hcv, verdicts, undecidable_is_inactive=True)
    sens = R.two_by_two(hcv, verdicts, undecidable_is_inactive=False)
    assert lit == (1, 1, 1, 2, 0)      # c counted as not-active
    assert sens == (1, 0, 1, 1, 2)     # c and d excluded
    assert "PRIMARY" in (ROOT / "scripts"
                         / "phase4_hcv_retrospective.py").read_text()


def test_a_compound_active_in_one_row_is_active():
    """"active if classify returns 'active' for at least one coronavirus
    row" -- so a compound inactive in four assays and active in a fifth is
    active, and the prereg says so."""
    a, b, c, d, _ = R.two_by_two(
        {"x"}, {"x": {"inactive", "active", None}},
        undecidable_is_inactive=True)
    assert (a, b) == (1, 0)


# ------------------------------------------------------------ end to end
def test_it_runs_and_names_the_protocol(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["x", "--db", str(_fixture(tmp_path))])
    code = R.main()
    out = capsys.readouterr().out
    assert code in (0, 1)
    assert "docs/prereg_phase4_hcv.md" in out
    assert "Nothing below departs from it." in out


def test_below_the_floor_it_does_not_run(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv",
                        ["x", "--db", str(_fixture(tmp_path, arm_b=1))])
    code = R.main()
    out = capsys.readouterr().out
    assert code == 1
    assert "BELOW THE REGISTERED FLOOR" in out
    assert "NOT RUN." in out


def test_the_null_branch_says_nothing_is_confirmed(tmp_path, monkeypatch,
                                                   capsys):
    """The gate predicts a null, so the common branch must not read as
    support for it."""
    monkeypatch.setattr(sys, "argv",
                        ["x", "--db", str(_fixture(tmp_path, broad_docs=14))])
    R.main()
    out = capsys.readouterr().out
    if "NO RESULT" in out:
        assert "remains an UNTESTED safety measure" in out
    assert "three-fold enrichment and nothing finer" in out


def test_a_missing_database_exits_two(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["x", "--db", str(tmp_path / "no.db")])
    assert R.main() == 2
    assert "no ChEMBL database at" in capsys.readouterr().out


def test_it_writes_no_results_file():
    """The verdict belongs beneath the prereg, not in a results file
    stamped against a release this measurement never reads."""
    src = (ROOT / "scripts" / "phase4_hcv_retrospective.py").read_text()
    assert "write_results" not in src
    assert "Read-only. Writes nothing" in src
