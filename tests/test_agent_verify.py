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
    'he said "yes" [E1]', "(see above) [E1]",
])
def test_a_properly_closed_answer_is_not_called_truncated(ending):
    v = verify(ending, facts=["f [E1]"], citable_edge_ids=["E1"])
    assert "truncated_answer" not in [f.code for f in v.findings], ending


# "a list [E1]:" USED TO BE IN THE LIST ABOVE, asserted as properly
# closed. It is not: an answer ending by promising a list and not
# delivering one is the definition of cut off. The case carried no
# rationale of its own and sat among five genuinely complete endings, so
# it reads as having been included by analogy with the other punctuation
# rather than decided. Replaced by its opposite, with the reason, because
# contradicting an existing assertion quietly is how the next person
# re-adds it.
@pytest.mark.parametrize("ending", [
    "a list [E1]:",
    "The direct-acting routes for this compound are:",
    "VIRUS NOT VALIDATED:",
    "Nirmatrelvir inhibits nsp5 [E1]. VIRUS NOT VALIDATED:",
    "measured in Vero E6 cells...",
    "measured in Vero E6 cells\u2026",
])
def test_an_answer_ending_open_is_truncated(ending):
    v = verify(ending, facts=["nirmatrelvir inhibits nsp5 [chembl, 2021; E1]"],
               citable_edge_ids=["E1"])
    assert "truncated_answer" in [f.code for f in v.findings], ending


def test_a_bare_caveat_label_does_not_verify_clean():
    """THE WORST COMBINATION, and the reason the colon matters. A caveat
    is enforced by its LABEL and a label ends in a colon, so a model cut
    off one token after emitting its first label produced:

        "VIRUS NOT VALIDATED:"

    which passed EVERY check and returned no findings at all -- the
    label is present so no caveat was dropped, the colon ended the
    answer so it was not truncated, the sentence is under six words so
    it is not an uncited claim, and it carries no citation and no
    number. A fully verified answer that says nothing."""
    v = verify("VIRUS NOT VALIDATED:", facts=["f [chembl, 2021; E1]"],
               citable_edge_ids=["E1"], required_warnings=[CAVEAT])
    assert not v.passed, v.report()
    assert "truncated_answer" in [f.code for f in v.findings]


def test_a_colon_inside_a_finished_answer_is_fine():
    """The check is about how the answer ENDS. A caveat renders as
    "LABEL: prose", so a colon mid-answer is the normal case and
    rejecting it would fail every answer that carries a caveat."""
    v = verify("VIRUS NOT VALIDATED: reliability is unknown. "
               "Nirmatrelvir inhibits nsp5 [E1].",
               facts=["nirmatrelvir inhibits nsp5 [chembl, 2021; E1]"],
               citable_edge_ids=["E1"])
    assert "truncated_answer" not in [f.code for f in v.findings]


def test_an_empty_answer_is_not_reported_as_truncated():
    """Nothing to truncate. The refusal and no-model paths both arrive here."""
    v = verify("", facts=[], citable_edge_ids=[])
    assert "truncated_answer" not in [f.code for f in v.findings]


# ------------------------------------------- a number with a unit is a claim
def test_a_fabricated_single_digit_potency_is_caught():
    """THE HOLE THIS CLOSES. _numbers dropped every integer at or under ten
    as "an ordinal or list index", so a claimed "Ki of 3 nM" verified clean
    against facts reading 190 nM -- in the check this module's own docstring
    calls the one that matters most. Single-digit potencies are not
    hypothetical here: nirmatrelvir's own Ki on v0.1 is 3 nM, so the exempted
    range is exactly where the most potent compounds live.
    """
    facts = ["nirmatrelvir INHIBITS nsp5, ic50_nm=75, ki_nm=190 [Eabc123]"]
    for claimed in ("3", "9", "11", "750"):
        v = verify(f"Nirmatrelvir inhibits nsp5 with a Ki of {claimed} nM "
                   f"[Eabc123].", facts, ["Eabc123"])
        assert not v.passed, f"{claimed} nM passed"
        assert any(f.code == "ungrounded_number" for f in v.findings)


def test_a_year_with_a_unit_after_it_is_still_a_claim():
    """The year rule had the same hole: "1960 nM" was dropped as a year."""
    v = verify("The compound was measured at 1960 nM [Eabc123].",
               ["compound INHIBITS target, ic50_nm=75 [Eabc123]"], ["Eabc123"])
    assert any(f.code == "ungrounded_number" for f in v.findings)


def test_ordinals_and_years_without_units_are_still_exempt():
    """So the fix is targeted and does not turn every list marker into a
    finding. A number is a magnitude claim when it carries a unit."""
    v = verify("Route 1 of 3 comes from the 2021 study and is sound here "
               "[Eabc123].",
               ["compound INHIBITS target, ic50_nm=75 [Eabc123]"], ["Eabc123"])
    assert not any(f.code == "ungrounded_number" for f in v.findings)


def test_the_facts_are_taken_whole_not_filtered():
    """The check is "every number in the ANSWER appears in the facts", so
    filtering the FACTS can only invent findings. It was symmetric while both
    sides dropped the same tokens; once the answer side became unit-aware it
    was not, and an answer reading "5 nM" against a fact reading "ic50_nm=5"
    was reported as an ungrounded figure."""
    v = verify("Nirmatrelvir has a best potency of 5 nM across its routes "
               "[Eabc123].",
               ["nirmatrelvir: 5 routes, best ic50_nm=5 [Eabc123]"],
               ["Eabc123"])
    assert v.passed, [f.detail for f in v.findings]


# ------------------------------------ numbers are paired with their units (A8)
#
# WHAT THIS IS FOR. The number check held BARE numeric tokens, so any
# number in the facts could be reattached to any claim with any unit.
# Reproduced against a fact reading `IC50 19 nM`, every one of these
# verified clean before 2026-11-02:
#
#     "IC50 19 uM"              a 1000-fold error
#     "19 mM"                   a million-fold error
#     "inhibits nsp5 by 19%"    a potency reread as an efficacy
#     "inhibits nsp5 19-fold"   the same, as a ratio
#
# The earlier UNIT_RE work made the ordinal exclusions unit-aware -- it
# stopped "3 nM" being dropped as a small number -- and never checked the
# pairing. Found by external audit.
_POTENCY_FACT = "nirmatrelvir INHIBITS nsp5 (IC50 19 nM) [E1]"


def _fails(answer: str, facts: list[str]) -> list[str]:
    v = verify(answer, facts, ["E1"])
    return [f.code for f in v.findings if f.severity == "fail"]


@pytest.mark.parametrize("answer,why", [
    ("nirmatrelvir inhibits nsp5 with an IC50 of 19 uM [E1]", "1000-fold"),
    ("nirmatrelvir inhibits nsp5 with an IC50 of 19 mM [E1]", "million-fold"),
    ("nirmatrelvir inhibits nsp5 with an IC50 of 19 pM [E1]", "1000-fold down"),
    ("nirmatrelvir inhibits nsp5 by 19% [E1]", "potency as efficacy"),
    ("nirmatrelvir inhibits nsp5 19-fold [E1]", "potency as a ratio"),
    ("nirmatrelvir inhibits nsp5 19x [E1]", "the same, abbreviated"),
])
def test_the_right_number_with_the_wrong_unit_is_a_finding(answer, why):
    assert "ungrounded_number" in _fails(answer, [_POTENCY_FACT]), why


def test_the_right_number_with_the_right_unit_passes():
    """The check must not fire on the answer it exists to allow."""
    assert _fails(
        "nirmatrelvir inhibits nsp5 with an IC50 of 19 nM [E1]",
        [_POTENCY_FACT]) == []


@pytest.mark.parametrize("spelling", ["19 nM", "19 nm", "19nM", "19 nmol"])
def test_spellings_of_one_unit_compare_equal(spelling):
    """Pairing must not turn a formatting difference into a fabrication
    report -- that would make the check unusable and it would be turned
    off, which is worse than the hole it closes."""
    assert _fails(f"the IC50 is {spelling} [E1]", [_POTENCY_FACT]) == []


def test_a_unit_carried_in_the_fact_key_still_grounds_the_answer():
    """A fact may write `ic50_nm=19` rather than `IC50 19 nM`. The
    module's own docstring records an earlier change breaking exactly
    this, so the facts side reads the key suffix too."""
    assert _fails("the IC50 is 5 nM [E1]",
                  ["qualifiers ic50_nm=5 for this edge [E1]"]) == []


def test_a_fact_grounds_the_same_number_mentioned_without_a_unit():
    """The facts side emits the bare number as well as the pair. An
    answer that names the figure without repeating the unit is vague,
    not fabricated, and reporting it would bury the real findings."""
    assert _fails("the value is 19 [E1]", [_POTENCY_FACT]) == []


def test_a_number_absent_from_the_facts_is_still_a_finding():
    """The original purpose of the check, which the pairing must not
    weaken."""
    assert "ungrounded_number" in _fails(
        "the IC50 is 750 nM [E1]", [_POTENCY_FACT])


def test_a_fold_change_in_the_facts_grounds_a_fold_change_in_the_answer():
    """Allowing a hyphen as the unit separator must not make every
    hyphenated number ungrounded."""
    assert _fails("a 19-fold increase [E1]",
                  ["measured a 19-fold increase [E1]"]) == []


def test_two_potencies_cannot_be_swapped_between_claims():
    """THE SHAPE THE AUDIT CALLS SET CONTAINMENT. Both numbers and both
    units are present in the facts, so a bare-token check passes a
    claim that attaches each to the wrong edge."""
    facts = ["drugA INHIBITS nsp5 (IC50 19 nM) [E1]",
             "drugB INHIBITS nsp12 (IC50 400 uM) [E2]"]
    # Right pairing.
    assert verify("drugA reads 19 nM and drugB reads 400 uM [E1] [E2]",
                  facts, ["E1", "E2"]).passed
    # Units swapped: 19 uM and 400 nM are each individually "present".
    bad = [f.code for f in verify(
        "drugA reads 19 uM and drugB reads 400 nM [E1] [E2]",
        facts, ["E1", "E2"]).findings if f.severity == "fail"]
    assert "ungrounded_number" in bad


# ---------------------------------- the unit vocabulary itself, not just pairs
def test_no_alias_is_unreachable():
    """DEAD CODE THAT LOOKS LIKE NORMALISING IS WORSE THAN NONE.
    Python's alternation is ordered, and UNIT_RE listed `nm` before
    `nmol`, so "nmol" was captured as "nm" and three entries in
    _UNIT_ALIASES could never fire. Found by a mutant that deleted the
    whole alias map with no test noticing.
    """
    from kgav.agent.verify import UNIT_RE, _UNIT_ALIASES

    unreachable = [
        k for k in _UNIT_ALIASES
        if not UNIT_RE.match(" " + k, 0)
        or UNIT_RE.match(" " + k, 0).group(0).strip(" -").lower() != k]
    assert unreachable == [], (
        f"UNIT_RE never produces {unreachable}, so their alias entries "
        f"are dead -- put the longer alternative first")


@pytest.mark.parametrize("a,b", [
    ("19 nmol", "19 nM"),      # alias, and the one the reorder fixed
    ("19 percent", "19%"),
    ("19 x", "19-fold"),
    ("19 µm", "19 um"),
    ("19 hours", "19 h"),
    ("19-fold", "19 fold"),    # the separator, not the unit
    ("19-day", "19 days"),
])
def test_two_spellings_of_one_unit_are_one_token(a, b):
    """A formatting difference must not read as a fabricated figure, or
    the check becomes unusable and gets switched off."""
    from kgav.agent.verify import UNIT_RE, _unit_token

    ta = _unit_token(UNIT_RE.match(a, a.index(" ") if " " in a else 2).group(0))
    tb = _unit_token(UNIT_RE.match(b, b.index(" ") if " " in b else 2).group(0))
    assert ta == tb, f"{a} -> {ta!r} but {b} -> {tb!r}"


def test_the_separator_is_stripped_from_the_unit():
    """UNIT_RE matches its own separator, so `.strip()` alone left the
    hyphen attached and "19-fold" did not equal "19 fold". Introduced
    by the hyphen change in the commit that needed it."""
    from kgav.agent.verify import _unit_token

    assert _unit_token("-fold") == "fold"
    assert _unit_token(" nM") == "nm"
    assert _unit_token("  -  ") == ""


# ------------------------------- an uncited claim is a failure, not a note (A9)
#
# TWO HOLES IN ONE CHECK, both reproduced before changing anything.
#
#   SEVERITY. `uncited_claim` was "warn", and `passed` is
#   `not any(f.severity == "fail")`. So an answer carrying a fabricated
#   sentence with no citation was reported as verified.
#
#   SCOPE. HEDGE_STARTS exempted the WHOLE SENTENCE by prefix, and every
#   caveat label was on the list. The audit's example:
#
#     "CALIBRATION: scoring for this virus is fully validated and the
#      ranking can be trusted."
#
#   -- a fabrication inverting the one thing that caveat exists to say,
#   waved through because it began with the label.
_FACT = "nirmatrelvir INHIBITS nsp5 (IC50 19 nM) [E1]"


def _cav():
    from kgav.agent.guards import caveat
    return caveat("cell_context",
                  "the chloroquine route was measured in Vero cells, "
                  "which lack TMPRSS2.")


@pytest.mark.parametrize("claim", [
    # The audit's own example, verbatim.
    "CALIBRATION: scoring for this virus is fully validated and the "
    "ranking can be trusted.",
    # Any other caveat label works the same way.
    "THIN ROUTE: this compound is the best available option for every "
    "patient.",
    # And the generic hedges, which were never only about caveats.
    "Note that nirmatrelvir cures every coronavirus infection in humans.",
    "This compound cures every coronavirus infection in humans everywhere.",
    "These results show the drug is safe and effective in all patients.",
    "However, this compound eradicates the virus in every tested host.",
])
def test_a_claim_behind_a_hedge_or_a_label_is_a_failure(claim):
    v = verify(claim, [_FACT], ["E1"])
    assert "uncited_claim" in [f.code for f in v.findings], claim
    assert v.passed is False, (
        "an uncited claim must fail the answer, not annotate it")


def test_a_restated_caveat_is_still_exempt():
    """_is_supplied covers this BY CONTENT -- its own docstring says it
    replaces prefix-matching the label, and the labels stayed on the
    hedge list anyway. Removing them must not start reporting genuine
    caveats as fabrications, or the check gets switched off."""
    cav = _cav()
    assert verify(str(cav), [_FACT], ["E1"], required_warnings=[cav]).passed


def test_half_a_restated_caveat_is_still_exempt():
    """A model legitimately borrows a clause and ends it with a full
    stop where the brief had a comma. Containment then failed on the
    punctuation alone -- harmless while this was a warning, and a
    broken answer once it was not."""
    cav = _cav()
    half = "CELL CONTEXT: the chloroquine route was measured in Vero cells."
    v = verify(half, [_FACT], ["E1"], required_warnings=[cav])
    assert "uncited_claim" not in [f.code for f in v.findings], v.findings


def test_a_short_remark_is_still_exempt():
    """The length test applies to the CLAIM, not to the sentence
    carrying it, so a genuine aside stays exempt while "Note that
    <twelve words of fabrication>" does not."""
    assert verify("However, no.", [_FACT], ["E1"]).passed
    assert verify("This is the only route found [E1].", [_FACT], ["E1"]).passed


def test_the_caveat_labels_are_no_longer_a_hedge_prefix():
    """Pinned directly, because the fix is a DELETION and a deletion is
    the easiest thing to put back by accident."""
    from kgav.agent.guards import CAVEAT_LABELS
    from kgav.agent.verify import HEDGE_STARTS

    overlap = {lbl.lower() for lbl in CAVEAT_LABELS.values()} & set(HEDGE_STARTS)
    assert overlap == set(), (
        f"{sorted(overlap)} exempt a whole sentence by prefix again; "
        f"_is_supplied covers a restated caveat by content")


def test_openings_are_stripped_repeatedly():
    """"CELL CONTEXT: this was ..." carries a label and a hedge, and one
    pass would leave the other doing the exempting."""
    from kgav.agent.verify import _strip_openings

    assert _strip_openings(
        "CELL CONTEXT: this compound binds nsp5").lower().startswith("compound")
    assert _strip_openings("Note that these results hold").lower().startswith(
        "results")
    assert _strip_openings("plain sentence") == "plain sentence"


# ----------------------------- a citation only exempts a SHOWN fact
# A13. `if EDGE_ID_RE.search(s): continue` exempted any sentence
# carrying any citable id, and a tool may return far more ids than it
# shows facts for -- pathways() returns 29 census lines and 1,500
# PARTICIPATES_IN ids, none of which appears in a line. Verified before
# the fix: the sentence below passed CLEAN with every caveat kept.
def test_a_citation_to_an_unshown_fact_does_not_exempt_a_claim():
    """The hole the uncited_claim promotion was meant to close, reached
    through the ids instead of the hedges."""
    facts = ["342 pathway(s) contain at least one dependency factor."]
    v = verify(
        "Promiscuousdrug eradicates SARS-CoV-2 in every patient tested "
        "and is the standard of care worldwide [E1001].",
        facts, ["E1001", "E1002"])
    assert not v.passed
    assert any(f.code == "uncited_claim" for f in v.findings), v.findings
    # NOT invented: the id really is citable. That is what made it work.
    assert not any(f.code == "invented_citation" for f in v.findings)


def test_a_citation_to_a_shown_fact_still_exempts():
    """No regression. An id the model was shown a fact for is exactly
    what a citation is for."""
    facts = ["nirmatrelvir inhibits the viral protein nsp5 "
             "(IC50 19 nM) [chembl, 2021; E1]"]
    v = verify("Nirmatrelvir acts on the nsp5 protease by direct "
               "inhibition of its active site [E1].", facts, ["E1"])
    assert not any(f.code == "uncited_claim" for f in v.findings), v.findings


def test_one_shown_id_among_several_cited_is_enough():
    facts = ["nirmatrelvir inhibits the viral protein nsp5 [chembl, 2021; E1]"]
    v = verify("Nirmatrelvir engages the protease through its catalytic "
               "cysteine residue [E1, E9].", facts, ["E1", "E9"])
    assert not any(f.code == "uncited_claim" for f in v.findings), v.findings


def test_the_shown_set_comes_from_the_facts_not_the_allowed_list():
    """If shown_ids were derived from citable_edge_ids the fix would be a
    no-op, because that is the list the hole was built on."""
    facts = ["342 pathway(s) contain at least one dependency factor."]
    v = verify("This compound is curative in all recorded human cases "
               "of infection [E1].", facts, ["E1"])
    assert any(f.code == "uncited_claim" for f in v.findings), v.findings


# ------------------------------------------- A5: a label, not a word
# `dropped_warnings` was `marker not in low` -- a bare substring search over
# the lowercased answer -- so a caveat was satisfied by its label appearing
# as an ordinary English word. The paragraph below states none of the
# nineteen declared caveats and satisfied SEVEN of them.

_PROSE_USING_THE_WORDS = (
    "The evidence for this compound is a single binding assay. Its "
    "selectivity is unknown and the ranking below is calibration-free. "
    "This cytotoxic effect was contradicted in a second cell line. The "
    "pathway route is thin."
)


def test_no_declared_caveat_is_satisfied_by_prose_using_its_words():
    """EVERY label, derived from CAVEAT_LABELS rather than listed.

    The seven that used to pass -- CELL LINE, CYTOTOXIC, EVIDENCE,
    CALIBRATION, CONTRADICTED, RANKING, PATHWAY ROUTE -- are the caveats
    most likely to be REQUIRED, because they are the ones whose subject
    matter the answer is about. The enforcement failed hardest exactly
    where it mattered.
    """
    from kgav.agent.guards import CAVEAT_LABELS, caveat
    from kgav.agent.verify import dropped_warnings

    satisfied = [CAVEAT_LABELS[k] for k in CAVEAT_LABELS
                 if not dropped_warnings(_PROSE_USING_THE_WORDS,
                                         [caveat(k, "the real warning text")])]
    assert satisfied == [], (
        f"{satisfied} are satisfied by a paragraph that states none of them")


def test_every_declared_caveat_is_satisfied_by_its_labelled_sentence():
    """The other half, and not a formality: a check nothing can satisfy
    withholds every answer, which is how a verifier becomes a switch that
    is always off. Each caveat renders as "LABEL: text", so its own
    rendering must pass."""
    from kgav.agent.guards import CAVEAT_LABELS, caveat
    from kgav.agent.verify import dropped_warnings

    for key in CAVEAT_LABELS:
        c = caveat(key, "the real warning text")
        assert not dropped_warnings(f"Some prose first. {c} And after.", [c]), (
            f"{CAVEAT_LABELS[key]} cannot be satisfied by its own rendering")


def test_a_longer_label_cannot_answer_for_the_one_it_contains():
    """"EVIDENCE" is a substring of "PROXIMITY EVIDENCE", so an answer
    carrying the proximity caveat satisfied the evidence-tier one for free
    -- enforcement that can never bite, which guards.py already calls
    decoration in its comment on `held_out_label`. The pairs are derived,
    so a new label that swallows an existing one is covered here the day it
    is added rather than the day it is noticed."""
    from kgav.agent.guards import CAVEAT_LABELS, caveat
    from kgav.agent.verify import _LABEL_CONTAINERS, dropped_warnings

    by_label = {v.lower(): k for k, v in CAVEAT_LABELS.items()}
    pairs = [(inner, outer) for inner, outers in _LABEL_CONTAINERS.items()
             for outer in outers]
    assert pairs, (
        "no label contains another any more; delete this test rather than "
        "letting it pass vacuously")
    for inner, outer in pairs:
        outer_caveat = caveat(by_label[outer], "the outer warning")
        inner_caveat = caveat(by_label[inner], "the inner warning")
        assert dropped_warnings(str(outer_caveat), [inner_caveat]), (
            f"an answer carrying {outer.upper()!r} satisfies "
            f"{inner.upper()!r}, which can then never be required on its own")


def test_a_word_ending_in_a_label_is_not_that_label():
    from kgav.agent.guards import caveat
    from kgav.agent.verify import dropped_warnings

    c = caveat("evidence_tier", "tier 3 only")
    assert dropped_warnings("nonevidence: nothing here", [c])
    assert dropped_warnings("cross-evidence: nothing here", [c])
    assert not dropped_warnings("EVIDENCE: tier 3 only", [c])
