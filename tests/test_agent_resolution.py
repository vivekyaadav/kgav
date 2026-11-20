"""Which entity the agent decided the question was about.

FIVE DEFECTS, ALL OF THEM SILENT SUBSTITUTIONS, found on 2026-11-20 when
the frontend put the agent in front of someone who had not written it.
Each answered a different question from the one asked, without saying so:

  1. Two viruses named -> one answered about, the other dropped with no
     caveat, no trace entry and not even the "(assumed)" marker. The
     winner was whichever name sorted first in a dict built by
     DESCENDING NAME LENGTH, so "denv-2" beat "zika".
  2. A compound the graph does not hold -> the retrieval fell through to
     a PROFILE OF THE VIRUS and returned edge counts, with nothing
     anywhere saying the compound had not been found.
  3. "yellow fever" -> MONDO:0020502, the Disease node added by the
     indications layer three days earlier, so no taxon resolved and the
     default virus answered. Also "japanese encephalitis" and
     "tick-borne encephalitis".
  4. "JEV" -> nothing. It is in the index; MIN_SEARCH_NAME = 4 made the
     resolver skip it. It is this project's own canonical name for
     NCBITaxon:11072 and the name every table prints.
  5. A typo -> the default virus under a caveat.

THE COMMON SHAPE IS THE POINT. Each had a plausible local justification,
and together they meant the agent would answer about SARS-CoV-2 whenever
it did not understand the question. The frontend made that visible
because a reader who did not write the code does not read the trace.

AND THE FIX FOR 5 IS NOT AN AUTOCORRECT. Measured on this release: 6,499
pairs of DIFFERENT compounds have names one character apart, and all four
dengue serotypes are one edit from each other -- DENV-2 being the only
cell in this project that survives a power analysis, and DENV-1 one that
cannot be powered at any arm size. Correcting silently would move the
substitution rather than remove it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "data" / "releases" / "v0.1"
sys.path.insert(0, str(ROOT / "src"))

needs_release = pytest.mark.skipif(
    not (RELEASE / "nodes.jsonl").exists(),
    reason="no assembled release on disk")


@pytest.fixture(scope="module")
def orch():
    from kgav.agent.orchestrator import Orchestrator
    from kgav.agent.tools import GraphTools
    return Orchestrator(GraphTools(RELEASE))


def _trace(out) -> str:
    return " ".join(out.get("trace") or [])


# ============================================ 3 and 4: every name resolves
@needs_release
def test_every_register_name_resolves_to_its_own_taxon(orch):
    """THE CHECK WHOSE ABSENCE LET 3 AND 4 THROUGH.

    Not "a name resolves" and not "most of them do": every string the
    register declares for every virus it declares, against the taxon that
    declared it. Three were shadowed by a Disease and one was below a
    length floor, and nothing asked.
    """
    from kgav import virus_register as vr

    doc = vr.load()
    wrong = []
    for entry in doc["viruses"]:
        taxon = "NCBITaxon:" + str(entry["taxon"])
        if taxon not in orch.tools.nodes:
            continue                      # declared, not built: not this test
        for raw in [entry.get("name"), entry.get("full_name"),
                    *(entry.get("aliases") or [])]:
            if not raw:
                continue
            key = " ".join(str(raw).lower().split())
            got = orch.tools.names.get(key)
            if got != taxon:
                cls = orch.tools.nodes.get(got, {}).get("class") if got else None
                wrong.append(f"{key!r} -> {got} ({cls}), want {taxon}")
    assert not wrong, "register names that do not reach their taxon:\n  " \
                      + "\n  ".join(wrong)


@needs_release
@pytest.mark.parametrize("typed,expect", [
    ("what is known about yellow fever?", "NCBITaxon:11089"),
    ("what is known about JEV?", "NCBITaxon:11072"),
    ("what is known about japanese encephalitis?", "NCBITaxon:11072"),
    ("what is known about tick-borne encephalitis?", "NCBITaxon:11084"),
])
def test_the_four_names_that_were_broken(orch, typed, expect):
    out = orch.answer(typed, None)
    assert f"virus={expect}" in _trace(out), _trace(out)
    assert "(assumed)" not in _trace(out)


@needs_release
def test_a_disease_name_reaches_the_virus_that_causes_it(orch):
    """Derived from the graph's own CAUSES edges, so a disease added
    later needs no edit here. "covid-19" is a Disease node and names a
    virus as far as any reader is concerned."""
    out = orch.answer("what is known about covid-19?", None)
    assert "virus=NCBITaxon:2697049" in _trace(out), _trace(out)


# ======================================================== 1: two viruses
@needs_release
def test_two_viruses_are_not_silently_reduced_to_one(orch):
    out = orch.answer("explain remdesivir for Zika and DENV-2", None)
    assert out.get("clarification"), _trace(out)
    assert "ambiguous_virus" in _trace(out)
    for name in ("Zika", "DENV-2"):
        assert name in out["answer"], f"{name} is not named in the reply"


@needs_release
def test_the_dropped_virus_was_chosen_by_name_length(orch):
    """The mechanism, pinned so the fix cannot be undone by restoring a
    'pick the first one' that happens to pick differently."""
    res, _problem, _om = orch._resolve_mentions(
        "explain remdesivir for Zika and DENV-2")
    taxa = orch._viruses_in(res)
    assert len(taxa) == 2, taxa
    assert orch._virus_in(res) is None, (
        "_virus_in returned a single virus for a question naming two")


# =================================================== 2: unknown compound
@needs_release
def test_an_unknown_compound_is_not_answered_as_the_virus(orch):
    out = orch.answer("why might frobnicatinib work against SARS-CoV-2?",
                      None)
    assert out.get("clarification"), _trace(out)
    assert "no_compound_resolved" in _trace(out)
    assert "HOST_FACTOR_FOR" not in " ".join(out.get("facts") or []), (
        "the reply still contains a profile of the virus")


@needs_release
def test_a_known_compound_still_answers(orch):
    """The guard above must not be satisfied by refusing everything."""
    out = orch.answer("why might nirmatrelvir work against SARS-CoV-2?",
                      None)
    assert not out.get("clarification"), _trace(out)
    assert "tool=explain ok=True" in _trace(out)


@needs_release
def test_profile_and_evidence_still_take_any_node(orch):
    """EXPLAIN and TRIAGE are the compound intents; the others are not,
    and the no-compound guard must not reach them."""
    for q in ("what is known about SIGMAR1?",
              "what is known about SARS-CoV-2?"):
        out = orch.answer(q, None)
        assert not out.get("clarification"), (q, _trace(out))


# ======================================================= 5: suggestions
@needs_release
@pytest.mark.parametrize("typed,expect", [
    ("what drug should i test for yellew fever virus?", "Yellow fever"),
    ("what is known about zikaa?", "Zika"),
])
def test_a_near_miss_is_suggested_not_applied(orch, typed, expect):
    out = orch.answer(typed, None)
    assert "near_virus_name" in _trace(out), _trace(out)
    assert expect in out["answer"]
    assert "did you mean" in out["answer"]
    # Applied would mean an answer about that virus. It must not be one.
    assert "virus=" not in _trace(out), (
        "the suggestion was applied rather than offered")


@needs_release
def test_an_ambiguous_near_miss_lists_every_candidate(orch):
    """THE REASON THIS SUGGESTS RATHER THAN CORRECTS.

    The four dengue serotypes are one edit from each other. DENV-2 is the
    cell registered for a wet-lab screen at 10.7x enrichment; DENV-1
    cannot be powered at any arm size. A one-edit autocorrect between
    them inverts the conclusion.
    """
    out = orch.answer("explain remdesivir for denv-5", None)
    assert "near_virus_name" in _trace(out), _trace(out)
    for serotype in ("DENV-1", "DENV-2", "DENV-3", "DENV-4"):
        assert serotype in out["answer"], f"{serotype} not offered"


def test_the_suggester_never_returns_an_exact_match():
    """A name that IS correct must not be offered as a correction of
    itself, or every answerable question grows a clarification."""
    from kgav.agent.orchestrator import near_names

    cands = {"zika": "NCBITaxon:64320", "denv-1": "NCBITaxon:11053"}
    assert near_names("what about zika", cands) == {}


def test_the_suggester_is_bounded():
    """Two edits. Unbounded, "zika" would suggest every short name in the
    register and the reply would be a list of sixteen."""
    from kgav.agent.orchestrator import MAX_SUGGEST_EDITS, near_names

    assert MAX_SUGGEST_EDITS == 2
    cands = {"zika": "NCBITaxon:64320"}
    assert near_names("what about zikaaa", cands)        # 2 edits
    assert not near_names("what about zikaaaa", cands)   # 3


# ============================== 6: the agent stopped assuming a virus
#
# The owner's decision, 2026-11-20. guards.virus_assumed_caveat defended
# the default: "SARS-CoV-2 is the only virus any protocol here evaluated.
# Defaulting without saying so is [wrong]." Fourteen ChEMBL cells, three
# plate cells and ten prospective cells later, the premise is gone.

@needs_release
@pytest.mark.parametrize("typed", [
    "is chloroquine useful?",                 # UNKNOWN intent, a compound
    "could remdesivir be repurposed?",        # PROSPECT, a compound
    "what is the evidence that remdesivir inhibits nsp12?",  # falls to explain
])
def test_a_virus_decided_question_asks_instead_of_assuming(orch, typed):
    out = orch.answer(typed, None)
    assert "virus_not_named" in _trace(out), _trace(out)
    assert out.get("clarification")
    assert "does not name a virus" in out["answer"]


@needs_release
def test_the_condition_is_what_retrieval_uses_not_the_intent(orch):
    """"is chloroquine useful?" classifies as UNKNOWN and still reaches
    explain(compound, virus). Keying the rule on intent alone missed it,
    which is the case virus_assumed_caveat's own docstring cites."""
    from kgav.agent.guards import QueryKind, classify_intent
    assert classify_intent("is chloroquine useful?") is QueryKind.UNKNOWN
    out = orch.answer("is chloroquine useful?", None)
    assert "virus_not_named" in _trace(out)


@needs_release
@pytest.mark.parametrize("typed", [
    "what is known about SIGMAR1?",           # PROFILE of a protein
    "what is known about SARS-CoV-2?",        # PROFILE of a taxon
])
def test_a_question_the_virus_does_not_decide_still_answers(orch, typed):
    """The guard must not be satisfied by refusing everything. A profile
    of a node that is not a virus does not need one."""
    out = orch.answer(typed, None)
    assert not out.get("clarification"), _trace(out)
    assert "tool=profile ok=True" in _trace(out)


@needs_release
def test_the_assumed_caveat_still_fires_where_a_virus_IS_assumed(orch):
    """The caveat was not deleted, only made rarer. PROFILE of a
    non-virus node still carries a default, and still says so."""
    out = orch.answer("what is known about SIGMAR1?", None)
    assert "(assumed)" in _trace(out), _trace(out)
    assert any("VIRUS ASSUMED" in w for w in (out.get("warnings") or [])), (
        "the assumed-virus caveat stopped firing where it still applies")


def test_the_suggester_ignores_short_common_words():
    """"the" is two edits from "tbev", which is two thirds of a
    three-letter word. With a flat budget an ordinary English sentence
    produced a clarification about a virus nobody typed."""
    from kgav.agent.orchestrator import near_names

    cands = {"tbev": "NCBITaxon:11084", "zika": "NCBITaxon:64320"}
    assert near_names("what is the evidence for this", cands) == {}
    # And a long enough phrase still gets the full budget.
    assert near_names("yellew fever", {"yellow fever": "NCBITaxon:11089"})


@needs_release
def test_a_question_resolving_nothing_is_an_ENTITY_problem_not_a_virus_one():
    """"which drugs should we test?" names nothing the graph holds. The
    useful reply says what kinds of thing it does hold -- which the
    existing `problem` path already did -- not "name a virus", which
    answers a question the reader has not reached yet.

    Pinned because the first version of the ask ran before that path and
    took it over, and because the boundary between the two is the only
    thing deciding which of two true messages a reader gets.
    """
    from kgav.agent.orchestrator import Orchestrator
    from kgav.agent.tools import GraphTools
    o = Orchestrator(GraphTools(RELEASE))
    out = o.answer("explain the mechanism of frobnicatinib", None)
    assert "virus_not_named" not in _trace(out), _trace(out)
    assert "No entity" in (out.get("answer") or "")

    # PROSPECT is the exception: it is a complete question missing only a
    # virus, and its refusal quotes that virus's own figures.
    out = o.answer("which drugs should we test?", None)
    assert "virus_not_named" in _trace(out), _trace(out)
