"""What the agent says about a HAS_ANTIVIRAL_ACTIVITY_AGAINST edge (L1).

THE PREDICATE NAME IS WRONG FOR MOST OF ITS INSTANCES. On release v0.1,
`labels.classify` over the 26,017 edges carrying it reads 20,239 inactive
(77.8%), 5,232 active (20.1%) and 546 undecidable (2.1%) -- and the agent
rendered every one of them as "{s} shows activity against {o} (EC50 ... nM)",
dropping `relation` entirely. An edge reading "EC50 > 50000 nM" -- a screen
reporting that it tested to 50 uM and saw nothing -- was presented as
evidence of activity with the bound's own number quoted as a potency.

WHY NO TEST CAUGHT IT. Both label edges in `test_agent_tools.py`'s fixture
are active (EC50 100 nM and 200 nM, no relation), so the inactive branch was
never exercised, and the one assertion about the wording --
`all("shows activity against" in ln for ln in marked)` -- was satisfied by
the only case the fixture could produce. 442 agent tests passed over this
bug. A fixture too small to reproduce the condition is the defect this
suite finds most often.

So this file exercises every branch, and the last test asserts the release
actually contains the case that matters rather than trusting the fixture.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.agent.tools import GraphTools  # noqa: E402
from kgav.labels import (  # noqa: E402
    DEFAULT_ACTIVE_BELOW_NM,
    DEFAULT_INACTIVE_ABOVE_NM,
    classify,
)

PRED = "HAS_ANTIVIRAL_ACTIVITY_AGAINST"


class _Renderer:
    """Just the two methods, so the sentence can be tested without a graph.

    `_potency_note` is wrapped in staticmethod deliberately: it IS one on
    GraphTools, and assigning the plain function here instead makes
    `self._potency_note(value)` pass `self` as the value. That mistake
    cost a debugging round on the first run of this file, and it would
    have been read as a bug in the code under test.
    """

    _potency_note = staticmethod(GraphTools._potency_note)
    _activity_sentence = GraphTools._activity_sentence


def say(**quals) -> str:
    e = {"predicate": PRED, "qualifiers": quals}
    return _Renderer()._activity_sentence(e, "remdesivir", "SARS-CoV-2")


# ------------------------------------------------------------- the branches
def test_a_censored_bound_above_the_floor_is_reported_as_inactive():
    """THE BUG, stated as a test. This read "shows activity against
    SARS-CoV-2 (EC50 50000 nM)"."""
    out = say(relation=">", ec50_nm=50_000.0)
    assert "FOUND INACTIVE" in out
    assert "not evidence of activity" in out
    assert "shows activity" not in out


def test_an_exact_value_above_the_floor_is_also_inactive():
    """"= 50000 nM" is a measurement of inactivity too, not a weak hit.
    The relation is not what decides this -- the threshold is."""
    out = say(relation="=", ec50_nm=50_000.0)
    assert "FOUND INACTIVE" in out


def test_the_relation_is_always_shown_with_the_number():
    """A bound quoted without its relation is a different claim.

    `>` must survive into the sentence; `=` is dropped because "(= 75
    nM)" reads worse than "(75 nM)" and loses nothing.
    """
    assert "(> 50 uM)" in say(relation=">", ec50_nm=50_000.0)
    assert "(75 nM)" in say(relation="=", ec50_nm=75.0)
    assert "(= 75" not in say(relation="=", ec50_nm=75.0)


def test_an_active_value_keeps_its_potency_qualifier():
    """The same three bands `_potency_note` already gave INHIBITS, so a
    75 nM hit and a 7.28 uM hit do not read as equivalent claims."""
    assert ", potent" in say(relation="=", ec50_nm=75.0)
    assert ", moderate" in say(relation="=", ec50_nm=400.0)
    assert "WEAK" in say(relation="=", ec50_nm=7_280.0)


def test_a_weak_bound_decides_nothing_and_says_why():
    """"> 100 nM" proves neither activity nor its absence: the true value
    could be anywhere above 100 nM, including potent. Reporting it as
    either is the error, and the floor it fails to clear is named."""
    out = say(relation=">", ec50_nm=100.0)
    assert "decides nothing" in out
    assert "10 uM" in out
    assert "FOUND INACTIVE" not in out and "shows activity" not in out


def test_no_potency_value_is_not_a_bound_below_the_floor():
    """Two ways to decide nothing, and they are not the same to a reader.

    The first version of that branch said "a bound below the 10 uM
    inactivity floor" about an edge carrying no bound at all.
    """
    out = say(relation="=")
    assert "no potency value was recorded" in out
    assert "bound below" not in out


def test_an_ic50_is_read_when_there_is_no_ec50():
    assert "(250 nM)" in say(relation="=", ic50_nm=250.0)


def test_micromolar_is_written_as_micromolar():
    """50,000 nM read as a number is hard to place against a 10 uM cut."""
    assert "50 uM" in say(relation=">", ec50_nm=50_000.0)
    assert "7.28 uM" in say(relation="=", ec50_nm=7_280.0)


# --------------------------------------------------- one implementation only
@pytest.mark.parametrize("quals", [
    {"relation": ">", "ec50_nm": 50_000.0},
    {"relation": "=", "ec50_nm": 75.0},
    {"relation": ">", "ec50_nm": 100.0},
    {"relation": "=", "ec50_nm": 10_000.0},
    {"relation": "<", "ec50_nm": 400.0},
    {"relation": "<", "ec50_nm": 50_000.0},
])
def test_the_sentence_agrees_with_classify(quals):
    """The agent must not reach a different verdict from the evaluation.

    Not a restatement of `classify` -- it CALLS it, and this asserts the
    sentence matches. A second copy of the threshold logic in the agent
    is the two-implementations failure this repository keeps finding, and
    it would show up as the agent describing an edge the evaluation
    counts the other way.
    """
    e = {"predicate": PRED, "qualifiers": quals}
    verdict = classify(e, DEFAULT_INACTIVE_ABOVE_NM, DEFAULT_ACTIVE_BELOW_NM)
    out = _Renderer()._activity_sentence(e, "d", "v")
    if verdict == "inactive":
        assert "FOUND INACTIVE" in out, out
    elif verdict == "active":
        assert "shows activity against" in out, out
    else:
        assert "decides nothing" in out, out


def test_the_release_contains_the_case_this_file_is_about():
    """Not a formatting check: it asserts the population.

    If inactives ever stop carrying this predicate -- because they were
    routed to MEASURED_INACTIVE_AGAINST, or the predicate was renamed --
    then the branch above is dead and the docstring explaining why it
    exists is describing a graph that no longer exists. Either way this
    fails and the decision gets recorded rather than drifting.
    """
    release = ROOT / "data" / "releases" / "v0.1"
    if not (release / "edges.jsonl").exists():
        pytest.skip("release v0.1 is not built in this checkout")
    counts = {"active": 0, "inactive": 0, "undecidable": 0}
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] != PRED:
            continue
        v = classify(e, DEFAULT_INACTIVE_ABOVE_NM, DEFAULT_ACTIVE_BELOW_NM)
        counts[v or "undecidable"] += 1
    total = sum(counts.values())
    assert counts["inactive"] > counts["active"], (
        f"most edges carrying {PRED} are no longer inactive: {counts}. The "
        f"agent's renderer and its docstring are written for a graph where "
        f"they are")
    # THE SHARE IS READ OUT OF THE FILES THAT QUOTE IT, not written here.
    # This assertion was `f"{share:.1%}" == "77.0%"`, which is the defect
    # it was built to catch: a literal in the check, so the rebuild of
    # 2026-11-09 moved the graph to 77.8% and the test failed on its own
    # copy of the number rather than on tools.py's. The figure now has one
    # home per file and this compares them with the release.
    share = counts["inactive"] / total
    rendered = f"{share:.1%}"
    for path in (ROOT / "src/kgav/agent/tools.py",
                 ROOT / "schema/kg_schema.yaml"):
        quoted = set(re.findall(r"\b(\d\d\.\d)% of\b|inactive\s+[\d,]+\s+(\d\d\.\d)%",
                                path.read_text()))
        stated = {g for pair in quoted for g in pair if g}
        assert stated, f"{path.name} no longer quotes an inactive share"
        assert stated == {rendered.rstrip('%')}, (
            f"{counts['inactive']:,} of {total:,} ({rendered}) edges "
            f"carrying {PRED} are inactive and {path.name} states "
            f"{sorted(stated)}")


def test_the_schema_documents_what_the_predicate_means():
    """The name is misleading, so the schema has to say so.

    Any consumer of this graph reads the schema before it reads
    labels.py, and "HAS_ANTIVIRAL_ACTIVITY_AGAINST" invites exactly one
    reading. The `doc` is the only place a downstream user is told that
    three quarters of these edges measured the opposite, and the share
    is asserted against the release so the sentence cannot go stale.
    """
    import yaml
    doc_path = ROOT / "schema" / "kg_schema.yaml"
    spec = yaml.safe_load(doc_path.read_text())
    entry = next(e for e in spec["edge_classes"]
                 if e["predicate"] == PRED)
    assert "doc" in entry, (
        f"{PRED} carries no `doc`, so nothing tells a consumer that the "
        f"predicate name is not the claim")
    text = entry["doc"]
    assert "ACTIVE OR INACTIVE" in text
    assert "classify" in text
    release = ROOT / "data" / "releases" / "v0.1"
    if not (release / "edges.jsonl").exists():
        pytest.skip("release v0.1 is not built in this checkout")
    inactive = total = 0
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] != PRED:
            continue
        total += 1
        if classify(e, DEFAULT_INACTIVE_ABOVE_NM,
                    DEFAULT_ACTIVE_BELOW_NM) == "inactive":
            inactive += 1
    assert f"{inactive / total:.1%}" in text, (
        f"the schema says a share the release does not: measured "
        f"{inactive / total:.1%} ({inactive:,} of {total:,})")
