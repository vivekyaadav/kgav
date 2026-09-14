"""Gate for the orchestrator: the pipeline is deterministic, the refusal comes
before retrieval, and the model is given no discretion.
"""
import json

import pytest

from kgav.agent.orchestrator import Orchestrator
from kgav.agent.tools import GraphTools

V = "NCBITaxon:2697049"
DRUG = "INCHIKEY:NIRMA"
CYTO = "INCHIKEY:CYTOX"


@pytest.fixture
def orch(tmp_path):
    def prot(pid, viral, sym=None, fam=None):
        p = {"taxon_id": V if viral else "NCBITaxon:9606", "is_viral": viral,
             "sequence_hash": "h", "reviewed": True}
        if sym:
            p["gene_symbol"] = sym
        if fam:
            p["protein_family"] = fam
        return {"id": pid, "class": "Protein", "properties": p}

    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        prot("UniProtKB:NSP5", True, fam="nsp5"),
        {"id": "KGAV:G", "class": "Gene",
         "properties": {"symbol": "rep", "taxon_id": V, "is_viral": True}},
        prot("UniProtKB:SIGMAR1", False, sym="SIGMAR1"),
        {"id": DRUG, "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "NIRMA", "label": "nirmatrelvir",
                        "chembl_id": "CHEMBL1", "is_approved": True,
                        "salt_collapsed": False, "stereo_collapsed": False}},
        {"id": CYTO, "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "CYTOX", "label": "chloroquine",
                        "chembl_id": "CHEMBL2", "is_approved": True,
                        "salt_collapsed": False, "stereo_collapsed": False}},
    ]

    def e(s, p, o, q=None):
        return {"subject": s, "predicate": p, "object": o, "qualifiers": q or {},
                "primary_knowledge_source": "infores:chembl", "evidence_tier": 1,
                "first_asserted_date": "2021-03-04"}

    edges = [
        e(DRUG, "INHIBITS", "UniProtKB:NSP5", {"assay_type": "biochemical", "ic50_nm": 25.0}),
        e("UniProtKB:NSP5", "ENCODED_BY", "KGAV:G"),
        e("KGAV:G", "BELONGS_TO", V),
        e(CYTO, "TARGETS", "UniProtKB:SIGMAR1",
          {"direction": "inhibitor", "assay_type": "binding"}),
        e("UniProtKB:SIGMAR1", "HOST_FACTOR_FOR", V,
          {"direction": "dependency", "screen_type": "CRISPRko", "cell_line": "CVCL_0574"}),
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))
    return Orchestrator(GraphTools(tmp_path))


# ------------------------------------------------------------------ refusal
def test_prospective_question_is_refused_before_retrieval(orch):
    """No partial answer may leak out alongside a refusal."""
    out = orch.answer("what should we test against SARS-CoV-2?")
    assert out["allowed"] is False
    assert out["facts"] == [] and out["citable_edge_ids"] == []
    assert "0.500" in out["answer"]


def test_refusal_happens_even_when_entities_resolve(orch):
    out = orch.answer("recommend a drug for SARS-CoV-2 like nirmatrelvir")
    assert out["allowed"] is False


# ------------------------------------------------------------- determinism
def test_same_question_gives_the_same_brief(orch):
    """A researcher who runs a query twice and gets two answers cannot cite
    the tool."""
    a = orch.answer("why might nirmatrelvir work against SARS-CoV-2?")
    b = orch.answer("why might nirmatrelvir work against SARS-CoV-2?")
    assert a["facts"] == b["facts"]
    assert a["citable_edge_ids"] == b["citable_edge_ids"]


# --------------------------------------------------------------- resolution
def test_longest_match_wins(orch):
    """hydroxychloroquine must not resolve as chloroquine -- the same
    substring hazard as Vero E6 inside Vero."""
    orch.tools.names["hydroxychloroquine"] = "INCHIKEY:HCQ"
    orch.tools.nodes["INCHIKEY:HCQ"] = {
        "id": "INCHIKEY:HCQ", "class": "SmallMolecule",
        "properties": {"label": "hydroxychloroquine"}}
    brief = orch.build_brief("explain hydroxychloroquine")
    assert "INCHIKEY:HCQ" in brief.resolved.values()
    assert CYTO not in brief.resolved.values()


def test_unresolvable_question_asks_rather_than_guesses(orch):
    out = orch.answer("explain the mechanism of aspirin")
    assert out["allowed"] is True
    assert "does not match" in out["answer"] or "No entity" in out["answer"]


# ------------------------------------------------------------------ routing
def test_explain_intent_retrieves_a_mechanistic_path(orch):
    out = orch.answer("why might nirmatrelvir work against SARS-CoV-2?")
    assert out["intent"] == "explain"
    assert any("M1 direct-acting" in f for f in out["facts"])
    assert out["citable_edge_ids"]


def test_triage_intent_ranks_the_named_compounds(orch):
    out = orch.answer("rank nirmatrelvir and chloroquine for SARS-CoV-2")
    assert out["intent"] == "triage"
    assert any("nirmatrelvir" in f for f in out["facts"])


def test_profile_intent_summarises(orch):
    out = orch.answer("what is known about SIGMAR1")
    assert out["intent"] == "profile"
    assert out["facts"]


# ----------------------------------------------------------------- warnings
def test_warnings_survive_into_the_brief(orch):
    out = orch.answer("why might chloroquine work against SARS-CoV-2?")
    joined = " ".join(out["facts"]) + " ".join(out["warnings"])
    assert "CELL CONTEXT" in joined or "SELECTIVITY UNKNOWN" in joined


def test_prompt_requires_every_warning_to_appear(orch):
    brief = orch.build_brief("why might chloroquine work against SARS-CoV-2?")
    p = brief.to_prompt()
    assert "MUST appear in your answer" in p
    assert "Use ONLY the facts listed below" in p


def test_prompt_lists_only_retrieved_edge_ids(orch):
    brief = orch.build_brief("why might nirmatrelvir work against SARS-CoV-2?")
    p = brief.to_prompt()
    for eid in brief.citable_edge_ids:
        assert eid in p


# -------------------------------------------------------------------- model
def test_model_is_optional_and_receives_only_the_brief(orch):
    seen = {}

    def fake_model(prompt):
        seen["prompt"] = prompt
        return "answer text"

    out = orch.answer("why might nirmatrelvir work against SARS-CoV-2?", fake_model)
    assert out["answer"] == "answer text"
    assert "FACTS:" in seen["prompt"]
    # The model must not be handed a query language or the graph itself.
    # Checking for "SELECT" alone false-positives on "SELECTIVITY UNKNOWN" --
    # the fourth substring collision in this project, after SARS-CoV inside
    # SARS-CoV-2, Vero inside Vero E6, and chloroquine inside
    # hydroxychloroquine.
    p = seen["prompt"]
    assert "SELECT " not in p.upper().replace("SELECTIVITY", "")
    assert "jsonl" not in p and "MATCH (" not in p


def test_no_model_returns_the_brief_unsynthesised(orch):
    out = orch.answer("why might nirmatrelvir work against SARS-CoV-2?")
    assert out["answer"] is None and "prompt" in out


def test_trace_records_each_stage(orch):
    out = orch.answer("why might nirmatrelvir work against SARS-CoV-2?")
    joined = " ".join(out["trace"])
    assert "intent=" in joined and "resolved=" in joined and "allowed=" in joined
