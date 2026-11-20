"""The agent can say what it is, how it works, and which way an edge runs.

"the output seems very hard wired, agent dont know how to talk or
explain the logics behind" -- the owner, 2026-11-21, and it was right.

FOUR THINGS WERE WRONG, all in how the agent was asked rather than in
what it knows:

  "tell me about yourself"   -> PROFILE, the node-profile intent, then
                                "No entity in this question matches a
                                node in the graph". It introduces
                                itself correctly for one exact phrasing
                                and fails on the obvious variants,
                                which is what reads as hardwired: a
                                fixed list with no notion of what the
                                question is ABOUT.
  "how i can use you ?"      -> same.
  "explain the logic behind
   these suggestions"        -> EXPLAIN, the compound-mechanism intent,
                                then the same dead end. It is the most
                                likely question after a ranking, and a
                                user told to name a compound concludes
                                the ranking was unexplainable.
  the prose layer            -> inverted BELONGS_TO. See below.

THE INVERSION IS THE INTERESTING ONE. Direction markers were added to
the evidence rows after a model wrote "DENV-2 has antiviral activity
against 4,162 entities". It then wrote "Zika belongs to the gene
KGAV:GENE_64320_POLY" from a row reading `<- BELONGS_TO` -- the gene
belongs to Zika. The marker fixed the predicate that carries a loud
warning paragraph and taught it nothing about arrows.

So the facts are now finished English, written in code from the edge
tuples, and the model only interprets. A text marker is not an
enforcement mechanism; that took two predicates to establish.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "data" / "releases" / "v0.1"
RESULTS = ROOT / "data" / "results"
sys.path.insert(0, str(ROOT / "src"))

needs_release = pytest.mark.skipif(
    not (RELEASE / "nodes.jsonl").exists(), reason="no release")


@pytest.fixture(scope="module")
def orch():
    from kgav.agent.orchestrator import Orchestrator
    from kgav.agent.tools import GraphTools
    return Orchestrator(GraphTools(RELEASE), results=RESULTS)


# ------------------------------------------------------ it can be asked
@pytest.mark.parametrize("q", [
    "tell me about yourself", "how i can use you ?", "what are your limits",
    "introduce yourself", "who are you", "what can i ask you",
])
@needs_release
def test_a_question_about_the_agent_reaches_the_agent(orch, q):
    """Never entity resolution. No node in the graph is the referent
    of "you", so resolving first can only fail."""
    b = orch.build_brief(q)
    assert b.intent in ("about", "capability"), (q, b.intent, b.trace)
    assert "No entity in this question matches" not in (
        b.clarification or ""), q


@pytest.mark.parametrize("q", [
    "can you explain the logic behind these suggestions ?",
    "how did you rank these", "why these compounds",
    "what does Tc mean", "explain the method", "what does resembles mean",
])
@needs_release
def test_a_question_about_the_method_is_answered_without_an_entity(orch, q):
    b = orch.build_brief(q)
    assert b.intent == "method", (q, b.intent, b.trace)
    txt = b.clarification
    assert "No entity" not in txt
    # The three things that make it an explanation and not a boast.
    assert "ECFP4" in txt and "Tanimoto" in txt
    assert "analogy" in txt.lower()
    for control in ("publication-disjoint", "symmetric k=10",
                    "strength-matched"):
        assert control in txt, f"{control} unmentioned in the method"
    # It must not claim the graph contributes.
    assert "sits at chance" in txt or "sit at chance" in txt


@needs_release
def test_asking_about_a_compound_still_reaches_the_compound(orch):
    """The self-reference and method patterns must not swallow real
    questions. This is the failure mode of adding phrase lists."""
    from kgav.agent.guards import QueryKind, classify_intent
    for q, want in (
            ("why might nirmatrelvir work against SARS-CoV-2?",
             QueryKind.EXPLAIN),
            ("rank these compounds for Zika", QueryKind.TRIAGE),
            ("what is known about remdesivir", QueryKind.PROFILE),
            ("what evidence supports remdesivir against SARS-CoV-2",
             QueryKind.EVIDENCE)):
        assert classify_intent(q) is want, q


# ------------------------------------------- it says which way edges run
@needs_release
def test_an_incoming_edge_never_makes_the_node_the_subject(orch):
    """THE REGRESSION THE REVIEWER ASKED FOR, generalised past the one
    predicate that was caught. For every incoming relation, the focal
    node must not be the grammatical subject of the clause."""
    for node in ("NCBITaxon:64320", "NCBITaxon:11060"):
        r = orch.tools.profile(node)
        assert r.ok
        focal = orch.tools._label(node)
        incoming = {d["predicate"] for d in r.data["directed"]
                    if d["direction"] == "<-"}
        assert incoming, f"{node} has no incoming edge to test"
        for line in r.data["relation_sentences"]:
            pred_of = [p for p in incoming
                       if p.replace("_", " ").lower() in line.lower()
                       or p in line]
            if not pred_of:
                continue
            assert not line.startswith(focal), (
                f"an incoming relation is voiced with {focal} as the "
                f"subject: {line!r}")


@needs_release
def test_the_specific_inversions_that_were_reported_are_gone(orch):
    r = orch.tools.profile("NCBITaxon:64320")
    said = " ".join(r.data["relation_sentences"])
    # 1. the virus does not have antiviral activity; compounds do.
    assert "Zika has measured antiviral activity" not in said
    assert "have measured antiviral activity against Zika" in said
    # 2. the gene belongs to Zika, not Zika to the gene.
    assert "Zika belongs to" not in said
    assert "belongs to Zika" in said
    # 3. and the small-count neighbours are NAMED.
    assert "Zika virus infectious disease" in said


@needs_release
def test_counts_read_as_prose_not_as_schema(orch):
    """"3,401 SmallMolecules" is a class identifier leaking out."""
    r = orch.tools.profile("NCBITaxon:64320")
    said = " ".join(r.data["relation_sentences"])
    for leak in ("SmallMolecule", "OrganismTaxon", "Protein:"):
        assert leak not in said, f"{leak} appears in prose: {said}"
    # And singular agrees.
    assert "1 gene belongs to" in said
    assert "belong to Zika" not in said


@needs_release
def test_the_prompt_gets_the_finished_clauses_not_the_arrow_rows(orch):
    """The model must not be the thing deciding direction. The arrow
    rows stay in verbalised for the audit trail shown in the UI."""
    r = orch.tools.profile("NCBITaxon:64320")
    assert any(ln.strip().startswith("<-") for ln in r.verbalised), (
        "the audit rows lost their direction markers")
    assert r.data["relation_sentences"], "no finished clauses produced"

    # THE BRIEF, not the source. Grepping build_brief for
    # "relation_sentences" passed with the injection replaced by `if
    # False: pass`, because the dead assignment above it still matched
    # -- a check that cannot fail is worse than no check.
    b = orch.build_brief("tell me about zika")
    facts = "\n".join(b.facts)
    assert "have measured antiviral activity against Zika" in facts, (
        "the brief hands the model count rows rather than finished "
        "clauses, so the model is still deciding direction")
    assert "belongs to Zika" in facts
    arrows = [f for f in b.facts if f.strip().startswith(("<-", "->"))]
    assert not arrows, (
        f"the brief still carries bare arrow rows the model must "
        f"interpret: {arrows}")


# --------------------------------------------- one policy, not two
def test_the_refusal_and_the_licensed_path_describe_one_policy():
    """"Which drugs should we test against Zika?" was refused with
    "cannot support prospective candidate generation", full stop, while
    "name a drug for yellow fever" returned twelve ranked compounds. A
    reader of both concludes the gate is arbitrary."""
    from kgav.agent.guards import QueryKind, check
    v = check(QueryKind.PROSPECT)
    assert not v.allowed
    assert "NOVEL chemical matter" in v.reason, (
        "the refusal does not say what it is scoped to")
    assert "licensed" in v.reason, (
        "the refusal does not tell the reader some viruses ARE licensed")


# ------------------------------------- it does not damage the repository
def test_the_agent_never_writes_into_data_results():
    """A DEFECT INTRODUCED BY THE FIX ABOVE IT.

    _shortlist invokes scripts/repurposing_candidates.py, which writes
    data/results/repurposing_candidates.json by default -- a REGISTERED
    run covering every virus, stamped and verified by `make
    check-results`. Called per-virus with no --out, it replaced that
    file with single-virus output on EVERY prospective question the
    agent answered. It was caught only because the file turned up
    modified at staging; nothing in the agent reported it, and a user
    running the frontend would have corrupted their own results without
    a word.
    """
    import inspect
    from kgav.agent import orchestrator as o

    src = inspect.getsource(o.Orchestrator._shortlist)
    assert "--out" in src, (
        "_shortlist invokes the shortlist script without --out, so it "
        "overwrites the registered results file")
    assert "tempfile" in src, (
        "_shortlist passes --out but not to a throwaway path")
    # The --out argument is the temp path and not anything else.
    assert "str(tmp)" in src, (
        "--out is passed but not the throwaway path")
