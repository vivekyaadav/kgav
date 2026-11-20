"""Can I repurpose drugs with this? -- answered, per virus, from the runs.

THE QUESTION THE SYSTEM COULD NOT ANSWER. "can i repurpose drugs using
this" returned "No entity in this question matches a node in the graph"
after silently assuming SARS-CoV-2. It names no compound because it is
not asking about one: it asks what the measurements license, which is a
question with an answer, and the answer differs sharply by virus.

AND THE ONE IT ANSWERED INCONSISTENTLY. scripts/repurposing_candidates.py
names HYDROXYCHLOROQUINE for HCoV-229E and MIGLUSTAT for DENV-2 from
this release under a registered verdict, while the agent said it does
not generate candidates. One of the two was wrong and it was the agent:
its refusal was written before any cell had been measured and never
revisited.

THE GATE IS NOW ONE PREDICATE, calibration.Capability.shortlist_licensed:
the cross-sectional control must clear AND the prospective control must
not contradict it. That second half is new and it matters -- Zika and
West Nile clear cross-sectionally and were measured ranking future
actives LAST, and the shortlist printed for them anyway.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "results"
RELEASE = ROOT / "data" / "releases" / "v0.1"

needs_release = pytest.mark.skipif(
    not (RELEASE / "nodes.jsonl").exists(), reason="no release")
needs_results = pytest.mark.skipif(
    not (RESULTS / "temporal_structure_scorer.json").exists(),
    reason="no results")


@pytest.fixture(scope="module")
def caps():
    import sys
    sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
    from kgav import virus_register as vr
    from kgav.calibration import capabilities
    from repurposing_candidates import tier

    names = {"NCBITaxon:" + v["taxon"]: v["name"]
             for v in vr.load()["viruses"]}
    return capabilities(RESULTS, names=dict(names), cross_tier=tier)


# ------------------------------------------- the predicate is the point
@needs_results
def test_a_cell_that_antiranks_prospectively_is_not_licensed(caps):
    """THE DEFECT THIS CLOSES. Clearing cross-sectionally was the whole
    gate, so Zika (0.703 cross-sectional) got a usable list although
    the prospective run measured it at 0.422 -- clearing BELOW chance,
    meaning it put the post-cutoff actives last. A shortlist ranks
    compounds nobody has screened, which is the prospective act."""
    con = [c for c in caps.values() if c.contradicted]
    assert con, (
        "no cell clears cross-sectionally while anti-ranking "
        "prospectively; this test asserts nothing")
    for c in con:
        assert c.cross_sectional == "clears"
        assert c.prospective == "below"
        assert not c.shortlist_licensed, (
            f"{c.name} anti-ranks prospectively and is still licensed")
        # And it must SAY which measurement refuses it.
        assert f"{c.prospective_auc:.3f}" in c.sentence()
        assert "LAST" in c.sentence()


@needs_results
def test_the_licensed_set_is_neither_everything_nor_nothing(caps):
    lic = [c for c in caps.values() if c.shortlist_licensed]
    assert lic, "nothing is licensed; the agent can answer no question"
    assert len(lic) < len(caps), (
        "every cell is licensed, so the gate is not gating")
    for c in lic:
        assert c.cross_sectional == "clears"
        assert c.prospective != "below"


@needs_results
def test_a_cell_below_chance_cross_sectionally_is_refused(caps):
    """DENV-4. Established by the 2026-11-21 amendment: it clears
    PROSPECTIVELY and is still correctly refused, because a shortlist
    ranks against ALL known actives, the cross-sectional setting."""
    d4 = next((c for c in caps.values() if c.name == "DENV-4"), None)
    if d4 is None:
        pytest.skip("DENV-4 not measured")
    assert d4.cross_sectional == "below"
    assert not d4.shortlist_licensed
    assert d4.prospective == "clears", (
        "DENV-4 no longer clears prospectively, so it is no longer the "
        "cross-protocol case the amendment describes")


@needs_results
def test_every_cell_gets_a_sentence_naming_its_own_numbers(caps):
    """A verdict with no figure is an assertion."""
    for c in caps.values():
        s = c.sentence()
        assert c.name in s
        if c.cross_auc is not None:
            assert f"{c.cross_auc:.3f}" in s, (
                f"{c.name}'s sentence quotes no cross-sectional AUC: {s}")


# --------------------------------------------- the agent uses that gate
@needs_release
@needs_results
def test_the_agent_answers_the_capability_question(caps):
    """Not "no entity matches a node in the graph"."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.agent.guards import QueryKind, classify_intent
    from kgav.agent.orchestrator import Orchestrator
    from kgav.agent.tools import GraphTools

    assert classify_intent("can i repurpose drugs using this") \
        is QueryKind.CAPABILITY
    orch = Orchestrator(GraphTools(RELEASE), results=RESULTS)
    b = orch.build_brief("can i repurpose drugs using this")
    assert b.intent == "capability"
    txt = b.clarification
    assert "No entity in this question matches" not in txt
    lic = sorted(c.name for c in caps.values() if c.shortlist_licensed)
    for name in lic:
        assert name in txt, f"{name} is licensed and unnamed in the answer"
    # The refusals are named too, with their numbers.
    for c in caps.values():
        if c.contradicted:
            assert c.name in txt and f"{c.prospective_auc:.3f}" in txt
    # And it does not claim the graph works.
    assert "chemistry, not the graph" in txt


@needs_release
@needs_results
def test_a_capability_question_does_not_assume_a_virus():
    """It used to assume SARS-CoV-2 and then fail to resolve anything."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.agent.orchestrator import Orchestrator
    from kgav.agent.tools import GraphTools

    orch = Orchestrator(GraphTools(RELEASE), results=RESULTS)
    b = orch.build_brief("can i repurpose drugs using this")
    assert not any("assumed" in t for t in b.trace), b.trace


@needs_release
@needs_results
def test_asking_for_names_is_still_routed_as_a_request(caps):
    """The split is modal. "can i repurpose drugs" asks what is
    possible; "what drugs should i repurpose for X" asks for names and
    must reach the prospective path, which either lists or refuses with
    the measurement -- never the capability summary."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.agent.guards import QueryKind, classify_intent

    for q in ("what drugs should i repurpose for DENV-2",
              "which drug should i test for zika"):
        assert classify_intent(q) is QueryKind.PROSPECT, q


@needs_release
@needs_results
def test_the_agent_without_results_does_not_invent_a_verdict():
    """A contract a caller can forget is not a contract: the results
    directory is optional, and absent it the honest answer is that no
    run is on disk."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.agent.orchestrator import Orchestrator
    from kgav.agent.tools import GraphTools

    orch = Orchestrator(GraphTools(RELEASE), results=None)
    txt = orch.build_brief("can i repurpose drugs using this").clarification
    assert "No registered run" in txt
    assert "licensed" not in txt.lower() or "no registered run" in txt.lower()


@needs_results
def test_the_tool_and_the_agent_share_one_tiering(caps):
    """Passed in rather than reimplemented, so "clears" cannot come to
    mean two things."""
    import inspect
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav import calibration

    src = inspect.getsource(calibration.capabilities)
    assert "cross_tier" in src
    doc = json.loads((RESULTS / "temporal_structure_scorer.json").read_text())
    assert set(doc["cells"]) <= set(caps), (
        "a prospectively measured cell is missing from the capability map")
