"""Gate for verification: hallucination is detected deterministically, not by
asking a model to be honest.
"""

import pytest

from kgav.agent import guards
from kgav.agent import verify as verify_mod
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


# ------------------------------------------- no citable edges (triage/profile)
# triage() and profile() return summary counts and never set edge_ids, so an
# empty citable list is their NORMAL shape, not an edge case. An earlier version
# returned early here and ran neither of the two checks below, so a fabricated
# potency and a dropped caveat both printed as "verified".
TRIAGE_FACTS = [
    "1. nirmatrelvir: best direct IC50 3 nM, SI 63291 (2 direct-acting, 0 host-directed routes)",
    "2. chloroquine: best direct IC50 160 nM, SI 15 (2 direct-acting, 1 host-directed routes)",
]
TRIAGE_WARN = [
    ("HOST-DIRECTED ROUTE: measured below chance against known-inactive "
     "compounds, and worse once cytotoxic compounds are excluded."),
]


def test_invented_number_is_caught_with_no_citable_edges():
    """A potency the graph never reported must not pass merely because the
    brief carried no edge ids to cite."""
    ans = ("HOST-DIRECTED ROUTE: shown for context only. Nirmatrelvir is the "
           "clear winner with an IC50 of 250 nM.")
    v = verify(ans, TRIAGE_FACTS, [], TRIAGE_WARN)
    assert not v.passed, v.report()
    assert any(f.code == "ungrounded_number" for f in v.failures)


def test_dropped_warning_is_caught_with_no_citable_edges():
    """The host-directed caveat is what makes a triage ranking interpretable.
    Dropping it must fail whether or not anything was citable."""
    ans = "Nirmatrelvir ranks first with a best direct IC50 of 3 nM."
    v = verify(ans, TRIAGE_FACTS, [], TRIAGE_WARN)
    assert not v.passed, v.report()
    assert any(f.code == "warning_dropped" for f in v.failures)


def test_faithful_triage_answer_still_passes_with_no_citable_edges():
    """The fix must not make an empty citable list unpassable: a faithful
    answer that keeps the caveat and invents no figure still verifies."""
    ans = ("Nirmatrelvir ranks first with a best direct IC50 of 3 nM and SI "
           "63291; chloroquine follows at 160 nM and SI 15. HOST-DIRECTED "
           "ROUTE: measured below chance against known-inactive compounds.")
    v = verify(ans, TRIAGE_FACTS, [], TRIAGE_WARN)
    assert v.passed, v.report()


# ------------------------------------------- caveats are keyed, not quoted
# Enforcement keys on the caveat's LABEL, so the prose may be paraphrased. The
# trap this closes: the key used to be the text before the first colon, or the
# first 40 characters when there was no colon, so a colon-less caveat demanded
# near-verbatim reproduction and editing any wording silently moved the goalposts.
def test_paraphrased_warning_passes_when_the_label_survives():
    w = guards.caveat("host_directed_route",
                      "measured below chance against known-inactive compounds.")
    ans = ("Nirmatrelvir ranks first. HOST-DIRECTED ROUTE: routes of this kind "
           "did worse than chance in evaluation, so read them as hypotheses "
           "rather than evidence.")
    v = verify(ans, TRIAGE_FACTS, [], [w])
    assert v.passed, v.report()
    assert not any(f.code == "warning_dropped" for f in v.findings)


def test_omitted_warning_still_fails_when_paraphrase_is_allowed():
    """Tolerating paraphrase must not tolerate omission."""
    w = guards.caveat("host_directed_route",
                      "measured below chance against known-inactive compounds.")
    ans = "Nirmatrelvir ranks first on direct-acting evidence and selectivity."
    v = verify(ans, TRIAGE_FACTS, [], [w])
    assert not v.passed, v.report()
    assert any(f.code == "warning_dropped" for f in v.failures)


def test_rewording_a_caveat_does_not_change_what_is_enforced():
    """The point of the explicit key: two wordings of the same caveat impose
    the same requirement, because neither the key nor the label came from the
    prose."""
    original = guards.caveat("ranking", "compounds are ordered by X then Y.")
    reworded = guards.caveat("ranking", "a completely different explanation.")
    assert original.key == reworded.key == "ranking"
    assert original.label == reworded.label == "RANKING"
    ans = "RANKING: ordering explained. Nirmatrelvir first."
    assert verify(ans, TRIAGE_FACTS, [], [original]).passed
    assert verify(ans, TRIAGE_FACTS, [], [reworded]).passed


def test_colonless_caveat_no_longer_demands_verbatim_prose():
    """The specific regression: the triage caveat has no colon of its own, so
    its key used to be its first 40 characters."""
    w = guards.caveat("ranking", "compounds are ordered by direct-acting "
                                 "evidence, then selectivity, then potency.")
    verbatim = f"{w} Nirmatrelvir first."
    paraphrase = "RANKING: ordered on direct evidence first. Nirmatrelvir first."
    assert verify(verbatim, TRIAGE_FACTS, [], [w]).passed
    assert verify(paraphrase, TRIAGE_FACTS, [], [w]).passed


def test_a_dropped_caveat_is_named_by_key_in_the_report():
    w = guards.caveat("cell_context", "runs through SIGMAR1.")
    v = verify("Nirmatrelvir ranks first.", TRIAGE_FACTS, [], [w])
    assert "cell_context" in v.report()


def test_unlabelled_string_warnings_are_still_enforced():
    """A plain string -- not built by guards.caveat(), or round-tripped through
    JSON -- must not silently stop being checked."""
    v = verify("Nirmatrelvir ranks first.", TRIAGE_FACTS, [],
               ["CELL CONTEXT: this path runs through SIGMAR1."])
    assert not v.passed
    assert any(f.code == "warning_dropped" for f in v.failures)


def test_citation_coverage_stays_gated_on_having_citable_edges():
    """A claim cannot be checked against citations that do not exist, so the
    uncited-claim check must stay off -- demanding citations where none are
    possible is what drove the model to invent them."""
    ans = ("Nirmatrelvir ranks first with a best direct IC50 of 3 nM, ahead of "
           "chloroquine on selectivity. HOST-DIRECTED ROUTE: below chance.")
    v = verify(ans, TRIAGE_FACTS, [], TRIAGE_WARN)
    assert not any(f.code == "uncited_claim" for f in v.findings)


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


# ------------------- a sentence the brief supplied is not the model's claim
CAVEAT = guards.caveat(
    "cell_context",
    "this path runs through SIGMAR1, whose contribution depends on the cell "
    "type used. The endosomal and sigma-receptor routes are prominent in Vero "
    "E6 cells and largely absent in TMPRSS2-expressing airway cells.")


def _uncited(v):
    return [f for f in v.findings if f.code == "uncited_claim"]


def test_a_caveats_later_sentences_are_not_uncited_claims():
    """The hedge-prefix check only ever exempted a caveat's FIRST sentence.
    The cell-context caveat is three, so its continuations were reported as
    uncited claims on text the brief supplied -- and forcing the caveats in
    verbatim made that happen on every answer."""
    answer = ("CELL CONTEXT: this path runs through SIGMAR1, whose "
              "contribution depends on the cell type used. The endosomal and "
              "sigma-receptor routes are prominent in Vero E6 cells and "
              "largely absent in TMPRSS2-expressing airway cells.")
    v = verify(answer, facts=["a fact [E1]"], citable_edge_ids=["E1"],
               required_warnings=[CAVEAT])
    assert _uncited(v) == [], v.report()


def test_a_sentence_the_model_invented_is_still_reported():
    """The exemption must not become a blanket pass."""
    answer = ("CELL CONTEXT: this path runs through SIGMAR1, whose "
              "contribution depends on the cell type used. Chloroquine is "
              "also an established treatment for COVID-19 in hospital "
              "settings worldwide.")
    v = verify(answer, facts=["a fact [E1]"], citable_edge_ids=["E1"],
               required_warnings=[CAVEAT])
    assert len(_uncited(v)) == 1


def test_a_short_sentence_is_not_exempted_by_accident():
    """"it does not transfer." occurs inside almost any corpus; exempting it
    would let a real uncited claim through."""
    assert not verify_mod._is_supplied("it does not transfer.",
                                       "a corpus where it does not transfer.")


# ---------------------------------------------- truncation is a failure
def test_a_truncated_answer_fails():
    """num_predict cut the chloroquine answer at "measured in Vero E" --
    mid-word, mid-caveat. Every other check passed: the caveat's label had
    already appeared and each citation so far resolved."""
    answer = ("CELL CONTEXT: this path runs through SIGMAR1 and the "
              "selectivity index was measured in Vero E")
    v = verify(answer, facts=["f [E1]"], citable_edge_ids=["E1"],
               required_warnings=[CAVEAT])
    codes = [f.code for f in v.findings]
    assert "truncated_answer" in codes
    assert not v.passed
    assert any(f.severity == "fail" for f in v.findings
               if f.code == "truncated_answer")


@pytest.mark.parametrize("ending", [
    "it is cited [E1].", "is that so [E1]?", "stop [E1]!",
    'he said "yes" [E1]', "(see above) [E1]", "a list [E1]:",
])
def test_a_properly_closed_answer_is_not_called_truncated(ending):
    v = verify(ending, facts=["f [E1]"], citable_edge_ids=["E1"])
    assert "truncated_answer" not in [f.code for f in v.findings], ending


def test_an_empty_answer_is_not_reported_as_truncated():
    """Nothing to truncate. The refusal and no-model paths both arrive here."""
    v = verify("", facts=[], citable_edge_ids=[])
    assert "truncated_answer" not in [f.code for f in v.findings]
