"""Gate for verification: hallucination is detected deterministically, not by
asking a model to be honest.
"""

from kgav.agent.verify import verify, verify_response

FACTS = [
    ("nirmatrelvir inhibits the viral protein nsp5 (IC50 75 nM, potent) "
    "[chembl, PMID:34798775, 2021; E335139]"),
    "nsp5 is encoded by rep [uniprot, 2020; E7]",
    "SARS-CoV-2 carries the gene rep [uniprot, 2020; E1]",
]
IDS = ["E335139", "E7", "E1"]
WARN = [("SELECTIVE: selectivity index 63291.1 (CC50/EC50, both measured in the "
        "same study).")]


# ------------------------------------------------------------- happy path
def test_a_grounded_cited_answer_passes():
    ans = ("Nirmatrelvir inhibits nsp5 with an IC50 of 75 nM [E335139]. "
           "nsp5 is encoded by the rep gene [E7], which SARS-CoV-2 carries [E1]. "
           "SELECTIVE: the selectivity index is 63291.1.")
    v = verify(ans, FACTS, IDS, WARN)
    assert v.passed, v.report()
    assert set(v.cited) == {"E335139", "E7", "E1"}


# ---------------------------------------------------------- hallucination
def test_invented_citation_is_caught():
    """The cheapest hallucination to detect, and the most damaging to trust."""
    ans = "Nirmatrelvir also inhibits nsp12 [E999999]. SELECTIVE: index 63291.1."
    v = verify(ans, FACTS, IDS, WARN)
    assert not v.passed
    assert any(f.code == "invented_citation" for f in v.failures)


def test_ungrounded_number_is_caught():
    """A transposed digit in a potency value is harder to notice than a wrong
    sentence."""
    ans = "Nirmatrelvir inhibits nsp5 with an IC50 of 57 nM [E335139]. SELECTIVE: index 63291.1."
    v = verify(ans, FACTS, IDS, WARN)
    assert not v.passed
    assert any(f.code == "ungrounded_number" for f in v.failures)


def test_years_and_ordinals_are_not_treated_as_claims():
    ans = ("Route 1: nirmatrelvir inhibits nsp5 at 75 nM [E335139], reported in "
           "2021 [E7]. SELECTIVE: index 63291.1.")
    v = verify(ans, FACTS, IDS, WARN)
    assert v.passed, v.report()


def test_thousands_separators_compare_equal():
    facts = ["remdesivir inhibits nsp12 (IC50 1560 nM) [E5]"]
    v = verify("Remdesivir inhibits nsp12 at 1,560 nM [E5].", facts, ["E5"])
    assert v.passed, v.report()


# --------------------------------------------------------------- warnings
def test_dropped_warning_fails_verification():
    """The check that stops a cell-context caveat being cut for brevity."""
    warn = [("CELL CONTEXT: this path runs through SIGMAR1, whose contribution "
            "depends on the cell type used.")]
    ans = "Chloroquine engages SIGMAR1 [E7]."
    v = verify(ans, FACTS, IDS, warn)
    assert not v.passed
    assert any(f.code == "warning_dropped" for f in v.failures)


def test_warning_present_passes():
    warn = ["CELL CONTEXT: measured in Vero E6."]
    ans = "Chloroquine engages SIGMAR1 [E7]. CELL CONTEXT: measured in Vero E6."
    v = verify(ans, FACTS, IDS, warn)
    assert not any(f.code == "warning_dropped" for f in v.findings)


# -------------------------------------------------------------- citations
def test_uncited_claim_is_flagged_as_a_warning_not_a_failure():
    ans = ("Nirmatrelvir inhibits nsp5 [E335139]. It is widely used as a "
           "treatment for coronavirus disease in many countries.")
    v = verify(ans, FACTS, IDS)
    assert any(f.code == "uncited_claim" for f in v.findings)
    assert v.uncited_sentences


def test_hedging_and_short_sentences_need_no_citation():
    ans = ("Nirmatrelvir inhibits nsp5 [E335139]. "
           "No other route was found. "
           "This is the only evidence available.")
    v = verify(ans, FACTS, IDS)
    assert not v.uncited_sentences


# --------------------------------------------------------------- refusals
def test_refusal_alone_passes():
    r = {"allowed": False, "answer": "This graph cannot support that.",
         "facts": [], "citable_edge_ids": [], "warnings": []}
    assert verify_response(r).passed


def test_refusal_with_an_answer_attached_fails():
    """A refusal accompanied by content is the pipeline producing exactly what
    the guard exists to prevent."""
    r = {"allowed": False, "answer": "This graph cannot support that.",
         "facts": [], "citable_edge_ids": [], "warnings": []}
    v = verify_response(r, answer="But you could try remdesivir.")
    assert not v.passed
    assert any(f.code == "refusal_contaminated" for f in v.failures)


# ------------------------------------------------------------- model check
def test_model_dispute_is_advisory_not_decisive():
    """The deterministic checks decide. A second model is a second opinion."""
    ans = "Nirmatrelvir inhibits nsp5 at 75 nM [E335139]. SELECTIVE: index 63291.1."
    v = verify(ans, FACTS, IDS, WARN,
               model_check=lambda a, f: {"supported": False, "reason": "unsure"})
    assert v.passed          # still passes: no deterministic failure
    assert any(f.code == "model_dispute" for f in v.findings)


def test_report_is_readable():
    v = verify("Something [E999].", FACTS, IDS)
    assert "verification failed" in v.report()
    assert "invented_citation" in v.report()
