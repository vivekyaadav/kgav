"""Gate for the agent's guardrails: the refusal is deterministic, and every
caveat the evaluation established is attached where a user will see it.
"""
import pytest

from kgav.agent.guards import (
    CELL_CONTEXT_SENSITIVE,
    QueryKind,
    check,
    classify_intent,
    evidence_note,
    path_warnings,
    selectivity_note,
)

SARS2 = "NCBITaxon:2697049"


# ------------------------------------------------------------------ intent
@pytest.mark.parametrize("q", [
    "what should we test against Nipah?",
    "suggest a drug for MERS",
    "recommend drugs for SARS-CoV-2",
    "predict which compounds will work",
    "find a treatment for COVID",
    "propose candidates for HCoV-229E",
])
def test_prospective_questions_are_recognised(q):
    assert classify_intent(q) is QueryKind.PROSPECT


@pytest.mark.parametrize("q,kind", [
    ("why might remdesivir work against MERS?", QueryKind.EXPLAIN),
    ("explain the mechanism of nirmatrelvir", QueryKind.EXPLAIN),
    ("rank these 20 compounds", QueryKind.TRIAGE),
    ("which of these should I follow up", QueryKind.TRIAGE),
    ("what evidence supports that interaction", QueryKind.EVIDENCE),
    ("what is known about ACE2", QueryKind.PROFILE),
])
def test_supported_questions_are_routed(q, kind):
    assert classify_intent(q) is kind


# ----------------------------------------------------------------- refusal
def test_prospective_queries_are_refused():
    """The question the graph cannot answer. Answering anyway produces
    confident, mechanistically coherent, wrong suggestions."""
    v = check(QueryKind.PROSPECT)
    assert v.allowed is False
    assert "0.50" in v.reason and "77%" in v.reason


def test_refusal_explains_and_offers_what_is_supported():
    """A bare refusal is unhelpful; the user needs the reason and the
    alternative."""
    reason = check(QueryKind.PROSPECT).reason
    assert "temporal" in reason.lower()
    assert "explain" in reason.lower() and "rank" in reason.lower()


def test_supported_queries_are_allowed():
    assert check(QueryKind.EXPLAIN, SARS2).allowed is True
    assert check(QueryKind.TRIAGE, SARS2).warnings == []


def test_unsupported_virus_carries_a_warning():
    """Cross-sectional AUC is ~0.50 for every metapath on every coronavirus
    except SARS-CoV-2."""
    v = check(QueryKind.EXPLAIN, "NCBITaxon:11137")
    assert v.allowed is True
    assert any("not validated" in w for w in v.warnings)


# ---------------------------------------------------------- path warnings
def _prot(sym):
    return {"properties": {"gene_symbol": sym}}


def test_cell_context_warning_fires_on_the_chloroquine_route():
    w = path_warnings([_prot("SIGMAR1"), _prot("TMEM97")])
    assert w and "CELL CONTEXT" in w[0]
    assert "Vero E6" in w[0] and "chloroquine" in w[0].lower()


def test_no_cell_context_warning_on_a_direct_acting_route():
    assert path_warnings([_prot("nsp5")], "M1 direct-acting") == []


def test_host_directed_routes_carry_a_performance_warning():
    w = path_warnings([_prot("SOMEGENE")], "M4")
    assert any("below chance" in x for x in w)


def test_endosomal_machinery_is_flagged():
    assert {"CTSL", "RAB7A", "ATP6V1A", "NPC1"} <= CELL_CONTEXT_SENSITIVE


# ------------------------------------------------------------ selectivity
def test_unknown_selectivity_does_not_read_like_a_pass():
    """No data and passed are different states. 38% of checkable compounds
    were cytotoxic."""
    note = selectivity_note(None, None)
    assert "UNKNOWN" in note and "38%" in note


def test_cytotoxic_compound_is_named_as_such():
    assert "CYTOTOXIC" in selectivity_note(2.0, True)


def test_selective_compound_reports_its_index():
    note = selectivity_note(50.0, True)
    assert "SELECTIVE" in note and "50.0" in note


def test_verified_flag_is_required_for_a_pass():
    """An index without the same-document flag is not a verified index."""
    assert "UNKNOWN" in selectivity_note(50.0, False)


# -------------------------------------------------------------- provenance
def test_proximity_methods_are_flagged_as_not_binding():
    note = evidence_note(1, "MI:1314")
    assert "PROXIMITY" in note and "1.3%" in note


def test_evidence_tiers_are_described():
    assert "experimental" in evidence_note(1)
    assert "Text-mined" in evidence_note(4)
