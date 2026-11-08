"""Gate for the agent's guardrails: the refusal is deterministic, and every
caveat the evaluation established is attached where a user will see it.
"""
from pathlib import Path

import pytest

from kgav.agent import guards
from kgav.agent.guards import (
    CAVEAT_LABELS,
    CELL_CONTEXT_SENSITIVE,
    ROUTED_KINDS,
    QueryKind,
    VirusRoutes,
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


# ------------------------------------------ candidate-seeking, off the list
#
# `prospective` held "which drug should" and not "which compounds
# should", so
#
#     which compounds should be considered for Yellow fever
#
# classified UNKNOWN, fell through the orchestrator's catch-all `else`,
# and came back as a PROFILE of the taxon node -- 1,030 antiviral
# measurements listed, for a virus with no route to any compound, whose
# refusal text exists and is correct and was never reached. The one
# phrasing on the list refused and the obvious paraphrase did not.
#
# Found by running the agent against one virus per tier in section 9
# rather than SARS-CoV-2 again; every earlier hand-check used the one
# tier where everything is present.
@pytest.mark.parametrize("q", [
    "which compounds should be considered for Yellow fever",
    "which drugs should be considered for DENV-4",
    "what compounds could treat Yellow fever",
    "what drugs would work for HCoV-HKU1",
    "any antivirals worth trying for TBEV",
    "list therapeutics for JEV",
    "candidates for Yellow fever",
    "shortlist for DENV-1",
    "which molecules look promising for Zika",
])
def test_candidate_seeking_phrasings_off_the_list_are_still_refused(q):
    """The backstop is structural -- a candidate noun plus a seeking
    word -- because a keyword list leaks by construction and this one
    leaked on the most natural phrasing of the question it exists to
    catch."""
    assert classify_intent(q) is QueryKind.PROSPECT, q


@pytest.mark.parametrize("q,kind", [
    # Each of these contains a candidate NOUN and would be swallowed by
    # the backstop, which since 2026-11-02 runs BEFORE the explicit
    # branches rather than after them.
    #
    # THE ORDERING ARGUMENT WAS THE WRONG WAY ROUND, and this comment
    # used to state it: "that ordering is the whole correctness
    # argument", meaning the backstop had to run last. Running last is
    # what let any permitted keyword appended to a prospective question
    # defeat it -- "which drugs could work against Zika and what is the
    # mechanism?" classified EXPLAIN. The backstop now runs first and
    # these cases are held by explicit carve-outs in
    # _NOT_SEEKING_CANDIDATES, which is why that list is load-bearing
    # rather than decoration.
    ("rank these compounds for SARS-CoV-2", QueryKind.TRIAGE),
    ("which of these compounds should I prioritise", QueryKind.TRIAGE),
    ("which paper shows this compound inhibits nsp5", QueryKind.EVIDENCE),
    ("what evidence links this drug to ACE2", QueryKind.EVIDENCE),
    ("tell me about this compound", QueryKind.PROFILE),
    ("why might this drug work against MERS", QueryKind.EXPLAIN),
    ("which pathways does this compound act through", QueryKind.PATHWAY),
])
def test_the_backstop_does_not_swallow_the_explicit_branches(q, kind):
    assert classify_intent(q) is kind, q


@pytest.mark.parametrize("q", [
    "what is the capital of France",
    "UniProtKB:P0DTD1",
    "how many nodes does this graph have",
])
def test_the_backstop_does_not_fire_on_everything(q):
    """A backstop that returned PROSPECT for any unrecognised text
    would refuse the whole long tail and look identical, in the test
    suite, to one that works."""
    assert classify_intent(q) is QueryKind.UNKNOWN, q


# ----------------------------------------------------------------- refusal
# REAL Calibration objects, not the _Cal stub further down this file. A stub
# whose author also writes the expectation cannot falsify the contract it
# stands for -- which is on record in this repo as the MITAB fixture that
# encoded the same wrong assumption as the regex it was testing.
def _real_cal(auc, lo, hi, *, protocol="single-screen", virus="Zika",
              n_pos=763, coverage_pos=0.025):
    from kgav.calibration import Calibration, ChannelCalibration
    cal = Calibration(virus=virus, protocol=protocol)
    cal.channels["M1"] = ChannelCalibration(
        name="M1", auc=auc, ci_lo=lo, ci_hi=hi, n_pos=n_pos, n_neg=2604,
        coverage_pos=coverage_pos, coverage_neg=0.002)
    return cal


def test_prospective_queries_are_refused():
    """The question the graph cannot answer. Answering anyway produces
    confident, mechanistically coherent, wrong suggestions.

    ASSERTS THE REFUSAL AND ITS MECHANISM, NOT A FIGURE. This read
    `"0.50" in v.reason and "77%" in v.reason`, pinning two hardcoded measured
    numbers into a test -- so the module could not stop quoting a stale figure
    without the suite failing, which is backwards. 77% was also the wrong
    shape: the post-cutoff absent fraction is 85.3% for SARS-CoV-2 and 99.0%
    for Zika, so no single constant could be right.
    """
    v = check(QueryKind.PROSPECT)
    assert v.allowed is False
    assert "cannot support prospective" in v.reason
    assert "novel chemical matter" in v.reason, "the durable reason, not a number"


def test_the_refusal_claims_no_result_without_a_calibration():
    """A VERDICT WITH NO FIGURE IS STILL A VERDICT (fix-list A7).

    This test used to assert the defect. It required that the
    no-calibration branch quote no AUC and DO say "every scorer performs
    at chance", under the docstring "absent a measured run, the limit is
    stated without a number" -- which treats a number as the only thing
    that can be unsupported. No caller passed `temporal` at all, so the
    sentence every user actually got asserted a measurement that had never
    been loaded, and called it measured in the same breath.

    What the branch must now do: describe the split, say that nothing was
    loaded, and state no result either way.
    """
    import re
    reason = check(QueryKind.PROSPECT).reason
    assert not re.search(r"AUC \d", reason), \
        f"quoted an AUC with no calibration supplied:\n{reason}"
    for claim in ("every scorer performs at chance",
                  "no scorer was even evaluable",
                  "which clears chance",
                  "measured rather than assumed"):
        assert claim not in reason, (
            f"the refusal asserts {claim!r} with no prospective measurement "
            f"loaded:\n{reason}")
    assert "temporal split" in reason, "the split's design is still the framing"
    assert "no scorer performance" in reason.lower()
    # The durable reason has to survive the withdrawal, or the refusal is
    # left resting on nothing at all.
    assert "novel chemical matter" in reason


def test_the_refusal_quotes_the_measured_figure_when_given_one():
    """Supplied a temporal calibration, it names the run and the number."""
    cal = _real_cal(0.502, 0.470, 0.534, protocol="post-2021")
    reason = check(QueryKind.PROSPECT, temporal=cal).reason
    assert "0.502" in reason and "post-2021" in reason


def test_a_prospectively_discriminating_channel_is_not_hidden_by_the_refusal():
    """If a scorer ever clears chance prospectively, a refusal built on
    "every scorer is at chance" has stopped being true, and the person is
    entitled to know that before being refused."""
    cal = _real_cal(0.70, 0.62, 0.78, protocol="post-2021")
    reason = check(QueryKind.PROSPECT, temporal=cal).reason
    assert "clears chance" in reason


def test_refusal_explains_and_offers_what_is_supported():
    """A bare refusal is unhelpful; the user needs the reason and the
    alternative."""
    reason = check(QueryKind.PROSPECT).reason
    assert "temporal" in reason.lower()
    assert "explain" in reason.lower() and "rank" in reason.lower()


def test_supported_queries_are_allowed():
    assert check(QueryKind.EXPLAIN, SARS2).allowed is True


def test_the_calibrated_virus_is_not_exempt_from_the_caveat_by_name():
    """THIS TEST USED TO ASSERT THE DEFECT: `check(TRIAGE, SARS2).warnings
    == []`, with no calibration passed at all.

    The gate was `virus != SUPPORTED_VIRUS`, so SARS-CoV-2 skipped the
    not-validated caveat by string comparison rather than by anything
    measured. Run the agent with the results file absent, unreadable, or
    measured on another release and every other virus correctly reports
    that its reliability is unknown, while the flagship virus -- the one
    a reader is most likely to ask about and most likely to believe --
    says nothing. Silence is the strongest claim the output can make,
    and it was being made by a constant.

    Handed nothing, the agent now says it knows nothing, about this
    virus like any other.
    """
    v = check(QueryKind.TRIAGE, SARS2)
    assert v.allowed is True
    assert any(getattr(w, "label", "") == "VIRUS NOT VALIDATED"
               for w in v.warnings), (
        f"no calibration was supplied and nothing was said: {v.warnings}")


def test_a_measured_discriminating_channel_is_reported_not_implied():
    """The other half: the good case gets a label of its own rather than
    being signalled by the absence of a bad one. docs/evaluation.md, on
    the identical mistake in its verdict column: "printed only as an
    absence of warning, the one condition every real result here depends
    on was left for the reader to notice."
    """
    v = check(QueryKind.TRIAGE, SARS2, calibration=_real_cal(0.675, 0.660,
                                                             0.690))
    assert v.allowed is True
    labels = [getattr(w, "label", "") for w in v.warnings]
    assert "CALIBRATION" in labels, labels
    assert "VIRUS NOT VALIDATED" not in labels


def test_unsupported_virus_carries_a_warning():
    """ASSERTS THE LABEL, WHICH IS WHAT IS ENFORCED.

    This checked lowercase "not validated" in the caveat's PROSE, against the
    module's own contract: "Editing a LABEL is now a visible, deliberate act;
    editing prose is free." verify.py enforces the label, so the prose was
    free to change and this test said otherwise -- it broke the moment the
    caveat started rendering from calibration rather than from a fixed
    sentence.
    """
    v = check(QueryKind.EXPLAIN, "NCBITaxon:11137")
    assert v.allowed is True
    assert any(getattr(w, "label", "") == "VIRUS NOT VALIDATED"
               for w in v.warnings)


def test_the_virus_caveat_asserts_nothing_about_other_viruses():
    """It used to read "approximately 0.50 for every metapath on all
    coronaviruses except SARS-CoV-2, and on four of them a degree-only
    baseline outperforms every mechanistic route". It is now five, not four
    (Zika 0.527, HCoV-OC43 0.552, SARS-CoV 0.696, MERS-CoV 0.723,
    HCoV-229E 0.646), and "all coronaviruses" does not describe a graph
    containing a flavivirus."""
    w = " ".join(check(QueryKind.EXPLAIN, "NCBITaxon:64320").warnings)
    for leaked in ("coronavirus", "four of them", "SARS-CoV-2", "0.50"):
        assert leaked not in w, f"still asserts {leaked!r} about other viruses:\n{w}"


def test_the_virus_caveat_renders_this_virus_own_measurement():
    """check() had the virus in hand and never consulted its calibration,
    while _channel_caveat twenty lines below did exactly that."""
    at_chance = _real_cal(0.511, 0.488, 0.534)
    w = " ".join(check(QueryKind.EXPLAIN, "NCBITaxon:64320",
                       calibration=at_chance).warnings)
    assert "single-screen" in w, f"must name the run:\n{w}"
    assert "not validated" in w.lower()

    clears = _real_cal(0.675, 0.660, 0.690)
    w2 = " ".join(check(QueryKind.EXPLAIN, "NCBITaxon:64320",
                        calibration=clears).warnings)
    assert "M1" in w2 and "clear chance" in w2, \
        f"a measured, discriminating channel must be reported as such:\n{w2}"


# ---------------------------------------------------------- path warnings
def _prot(sym):
    return {"properties": {"gene_symbol": sym}}


def _taxon(family, label="X"):
    """A virus node as ingest_spine actually writes it: with its family."""
    return {"id": "NCBITaxon:1", "class": "OrganismTaxon",
            "properties": {"label": label, "family": family}}


def test_cell_context_warning_fires_on_the_chloroquine_route():
    """Now supplies a CORONAVIRUS, because the explanation is family-scoped.

    This passed a path with no OrganismTaxon and asserted the SARS-CoV-2 text,
    so it was implicitly requiring that a path naming no virus be explained by
    coronavirus entry biology -- which is the bug this patch fixes.
    """
    w = path_warnings([_prot("SIGMAR1"), _prot("TMEM97"),
                       _taxon("Coronaviridae", "SARS-CoV-2")])
    assert w and "CELL CONTEXT" in w[0]
    assert "Vero E6" in w[0] and "chloroquine" in w[0].lower()


def test_a_flavivirus_route_is_not_explained_by_coronavirus_entry():
    """THE BUG. The genes are genuinely involved in flavivirus entry too --
    clathrin-mediated, acidification-dependent -- so the warning belongs on
    both. The EXPLANATION does not transfer: flaviviruses do not use TMPRSS2,
    and chloroquine's 2020 SARS-CoV-2 story is not theirs. An explanation
    about the wrong virus is worse than none, because it reads as specific
    knowledge."""
    w = path_warnings([_prot("RAB7A"), _prot("ATP6V1A"),
                       _taxon("Flaviviridae", "Zika")])
    assert w and "CELL CONTEXT" in w[0]
    assert "interferon" in w[0].lower(), f"expected the real reason:\n{w[0]}"
    for leaked in ("TMPRSS2", "chloroquine", "SARS-CoV-2", "airway"):
        assert leaked.lower() not in w[0].lower(), \
            f"coronavirus biology leaked into a flavivirus path: {leaked}\n{w[0]}"


def test_an_unknown_family_gets_no_borrowed_explanation():
    """A path with no virus, or a family this module has no account of, must
    say the dependence is unassessed rather than borrow one."""
    w = path_warnings([_prot("SIGMAR1")])
    assert w and "CELL CONTEXT" in w[0]
    assert "unassessed" in w[0]
    for leaked in ("TMPRSS2", "chloroquine", "interferon"):
        assert leaked.lower() not in w[0].lower(), f"borrowed {leaked}:\n{w[0]}"


# ------------------------------------------------- cell line, per family
def test_a_vero_measurement_is_explained_by_the_right_mechanism():
    """Vero is the standard line for BOTH families and the reason a result
    does not transfer differs: TMPRSS2 absence for a coronavirus, the type I
    interferon deletion for a flavivirus. One text said TMPRSS2 for both."""
    corona = guards.cell_line_note("Vero E6", family="Coronaviridae")
    assert "TMPRSS2" in corona

    flavi = guards.cell_line_note("Vero E6", family="Flaviviridae")
    assert "interferon" in flavi.lower()
    assert "TMPRSS2" not in flavi, f"TMPRSS2 is not a flavivirus route:\n{flavi}"


def test_the_interferon_fact_is_claimed_only_for_vero():
    """Vero's interferon locus deletion is a property of that line. BHK,
    Huh-7, HEL 299 and MRC5 are on the non-transferring list for their own
    reasons, and this must not assert the interferon fact about them."""
    other = guards.cell_line_note("Huh-7", family="Flaviviridae")
    assert "interferon" not in other.lower(), \
        f"asserted Vero's property about another line:\n{other}"
    assert "not to transfer reliably" in other


def test_a_transferring_cell_line_is_reported_plainly():
    note = guards.cell_line_note("Calu-3", family="Coronaviridae")
    assert "Calu-3" in note and "TMPRSS2" not in note


def test_no_cell_context_warning_on_a_direct_acting_route():
    assert path_warnings([_prot("nsp5")], "M1 direct-acting") == []


def test_host_directed_routes_carry_a_performance_warning():
    """It must always be attached. WHAT it says comes from the calibration --
    this test asserted the literal words "below chance", which was the
    hardcoded figure, so it passed while the number went stale."""
    w = path_warnings([_prot("SOMEGENE")], "M4")
    assert any("HOST-DIRECTED ROUTE" in x for x in w)
    assert any("mechanistic hypothesis" in x for x in w)


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


# --------------------------------- the host-directed caveat is measured
class _Ch:
    def __init__(self, auc, lo, hi, n_neg=7343, discriminates=False,
                 n_pos=89, n_reached_pos=15, evaluable=True):
        self.auc, self.ci_lo, self.ci_hi, self.n_neg = auc, lo, hi, n_neg
        self.discriminates = discriminates
        # Coverage decides whether the AUC describes the channel or the
        # imputed floor, so a stub without it reads as never evaluated.
        self.n_pos, self.n_reached_pos = n_pos, n_reached_pos
        self.evaluable = evaluable


class _Cal:
    def __init__(self, channels, protocol="single-screen/publication-disjoint",
                 usable=True):
        self.channels, self.protocol, self.usable = channels, protocol, usable


HOST_PATH = [{"id": "INCHIKEY:X", "class": "SmallMolecule", "properties": {}},
             {"id": "UniProtKB:Q99720", "class": "Protein",
              "properties": {"gene_symbol": "SIGMAR1", "is_viral": False}},
             {"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
              "properties": {"label": "SARS-CoV-2",
                             "family": "Coronaviridae"}}]


def _host_caveat(cal):
    out = guards.path_warnings(HOST_PATH, "M4 pathway-mediated", calibration=cal)
    hits = [w for w in out if "HOST-DIRECTED ROUTE" in w]
    assert len(hits) == 1, out
    return hits[0]


@pytest.mark.parametrize("auc,lo,hi", [
    (0.507, 0.447, 0.568), (0.391, 0.300, 0.482), (0.650, 0.600, 0.700),
])
def test_every_figure_in_the_caveat_comes_from_the_calibration(auc, lo, hi):
    """The caveat used to read "AUC 0.391-0.482", from a run nobody could
    identify; by the time the single-screen protocol existed the real values
    were 0.492-0.523 -- at chance, not below it. Grepping the source for the
    stale constant would also flag the comment explaining it, so this asserts
    the behaviour: the emitted figures track whatever calibration is handed
    in, and nothing else appears."""
    w = _host_caveat(_Cal({"M4": _Ch(auc, lo, hi)}))
    assert f"{auc:.3f}" in w and f"{lo:.3f}" in w and f"{hi:.3f}" in w
    # no OTHER three-decimal figure leaked in from a constant
    import re
    assert set(re.findall(r"0\.\d{3}", w)) == {f"{auc:.3f}", f"{lo:.3f}", f"{hi:.3f}"}


def test_the_caveat_quotes_the_measured_auc_and_names_the_protocol():
    cal = _Cal({"M4": _Ch(0.507, 0.447, 0.568)})
    w = _host_caveat(cal)
    assert "0.507" in w and "0.447" in w and "0.568" in w
    assert "7,343" in w                      # what it was measured against
    assert "single-screen/publication-disjoint" in w
    assert "indistinguishable from chance" in w


def test_without_a_calibration_the_route_is_uncalibrated_not_numbered():
    """A figure that might be stale is worse than no figure."""
    w = _host_caveat(None)
    assert "no measured performance is attached" in w
    assert not any(c.isdigit() for c in w.split("Read")[0])


def test_a_refused_calibration_is_treated_as_absent():
    cal = _Cal({"M4": _Ch(0.823, 0.805, 0.841)}, usable=False)
    w = _host_caveat(cal)
    assert "0.823" not in w
    assert "no measured performance is attached" in w


def test_an_unevaluable_channel_says_so():
    w = _host_caveat(_Cal({"M4": _Ch(None, None, None)}))
    assert "not evaluable" in w and "reliability is unknown" in w


def test_a_discriminating_channel_is_reported_as_above_chance():
    """And still as a hypothesis: beating chance on a ranking task is not a
    measured activity for the compound in front of the reader."""
    w = _host_caveat(_Cal({"M4": _Ch(0.70, 0.65, 0.75, discriminates=True)}))
    assert "0.700" in w and "above chance" in w
    assert "mechanistic hypothesis" in w


def test_a_direct_acting_route_gets_no_host_directed_caveat():
    out = guards.path_warnings(HOST_PATH, "M1 direct-acting",
                               calibration=_Cal({"M1": _Ch(0.513, 0.452, 0.574)}))
    assert not any("HOST-DIRECTED ROUTE" in w for w in out)


# ------------------------ no evaluation figure is written into the source
STALE_FIGURES = ("0.826", "0.823", "0.729", "0.660", "0.391", "0.482")


def _code_constants(path: Path) -> list[str]:
    """Every literal the module can EXECUTE, docstrings excluded.

    Parsed rather than grepped: the docstring explaining this bug necessarily
    names the figures it is about, and a grep cannot tell prose from a value a
    user will be shown.
    """
    import ast
    tree = ast.parse(path.read_text())
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            body = getattr(node, "body", None) or []
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                docstrings.add(id(body[0].value))
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and id(node) not in docstrings:
            if isinstance(node.value, (str, int, float)):
                out.append(str(node.value))
    return out


def test_no_hardcoded_evaluation_figure_remains_in_agent_text():
    """0.826 was BOTH a dead constant in guards and a separate string literal
    in tools' ranking caveat -- two copies with nothing keeping them equal,
    and the figure has since fallen to 0.726, 0.660 and 0.499. A number a
    reader is shown has to come from the run that measured it."""
    import kgav.agent.tools as tools_mod
    for mod in (guards, tools_mod):
        for const in _code_constants(Path(mod.__file__)):
            for stale in STALE_FIGURES:
                assert stale not in const, (
                    f"{stale} is an executable literal in {mod.__name__}: "
                    f"{const[:90]!r}")


def test_the_ranking_caveat_quotes_the_measured_direct_acting_figure():
    cal = _Cal({"M1": _Ch(0.513, 0.452, 0.574)})
    w = guards.direct_acting_caveat(cal)
    assert "0.513" in w and "0.452" in w and "0.574" in w
    assert "indistinguishable from chance" in w
    assert "evidence available rather than a demonstrated ability to rank" in w


def test_a_discriminating_direct_acting_route_is_stated_plainly():
    w = guards.direct_acting_caveat(_Cal({"M1": _Ch(0.78, 0.74, 0.82,
                                                    discriminates=True)}))
    assert "0.780" in w and "separates measured actives" in w


def test_without_a_calibration_the_ranking_claims_nothing():
    w = guards.direct_acting_caveat(None)
    assert "no measured performance is attached" in w
    assert not any(c.isdigit() for c in w)


def test_an_unevaluated_channel_is_not_given_an_auc():
    """M1 reads 0.513 while reaching 3 of 89 positives, so that figure
    describes the imputed floor. Quoting it in a caveat is the hardcoded
    figure's error in a new place."""
    cal = _Cal({"M4": _Ch(0.513, 0.453, 0.574, n_reached_pos=3,
                          evaluable=False)})
    w = _host_caveat(cal)
    assert "0.513" not in w
    assert "NOT been evaluated" in w and "3 of 89" in w


def test_the_ranking_caveat_also_refuses_an_unevaluated_figure():
    w = guards.direct_acting_caveat(
        _Cal({"M1": _Ch(0.513, 0.453, 0.574, n_reached_pos=3,
                        evaluable=False)}))
    assert "0.513" not in w
    assert "NOT been evaluated" in w
    assert "rather than a demonstrated ability to rank" in w


# ------------------------------------------------ the biology carries sources
def test_every_biology_claim_renders_a_source():
    """Every performance figure in this module is traceable to the evaluation
    that produced it, while the biology was asserted in bare prose -- four
    specific, falsifiable claims about deleted gene loci, entry routes and
    interferon restriction, with nothing a reader could check. A graph whose
    whole argument is provenance should not make its most confident
    statements the only unsourced ones.
    """
    rendered = [
        guards._cell_context_text(["NPC1"], "Flaviviridae"),
        guards._cell_context_text(["CTSL"], "Coronaviridae"),
        guards.cell_line_note("Vero E6", family="Flaviviridae"),
        guards.cell_line_note("Vero E6", family="Coronaviridae"),
    ]
    for text in rendered:
        assert "doi:" in text, f"no source in: {text[:70]}"


def test_the_four_claims_cite_the_four_sources():
    """Named individually, so removing a citation fails here rather than
    leaving a claim quietly bare."""
    flavi = guards._cell_context_text(["NPC1"], "Flaviviridae")
    for key in ("flavivirus_entry", "vero_ifn_deletion",
                "flavivirus_ifn_restriction"):
        assert guards.BIOLOGY_SOURCES[key][0] in flavi, key
    corona = guards._cell_context_text(["CTSL"], "Coronaviridae")
    assert guards.BIOLOGY_SOURCES["chloroquine_not_in_airway"][0] in corona


def test_cite_raises_on_an_unknown_key():
    """A claim whose source has been removed must break the render rather
    than quietly lose its citation."""
    with pytest.raises(KeyError):
        guards.cite("no_such_source")


def test_every_source_carries_a_resolvable_identifier():
    """A reference without a DOI is a sentence, not a citation. DOIs
    throughout rather than a mix: a citation style that varies by what
    happened to be findable invites the reader to wonder what the difference
    means."""
    for key, (short, full) in guards.BIOLOGY_SOURCES.items():
        assert short.startswith("doi:"), key
        assert short[4:] in full, f"{key}: short id absent from the reference"
        assert len(full) > 120, f"{key}: reference too thin to check"


def test_a_cited_caveat_still_passes_number_provenance():
    """THE INTERACTION. verify checks that every number in an answer appears
    in the facts, and a DOI is full of digits. Caveats are passed to verify
    as required_warnings, which are part of that corpus -- so a cited caveat
    reproduced in an answer must not read as a fabricated figure."""
    from kgav.agent.verify import verify
    warning = guards.cell_line_note("Vero E6", family="Flaviviridae")
    answer = f"{warning} The route is therefore not transferable [Eabc123]."
    v = verify(answer, ["compound INHIBITS target [Eabc123]"], ["Eabc123"],
               required_warnings=[warning])
    assert not [f for f in v.findings if f.code == "ungrounded_number"], \
        [f.detail for f in v.findings]


def test_a_measured_at_chance_calibration_is_labelled_not_validated():
    """ASSERTS THE LABEL, not the prose.

    Mutating _is_validated to accept `usable` alone survived the suite:
    the prose still read "scoring is not validated for this virus",
    because _virus_caveat picks its sentence from discriminating() while
    check() picked the KEY from the mutated predicate. Prose and label
    disagreed and only the prose was under test -- the same split this
    file's other docstring warns about, from the other side.
    """
    v = check(QueryKind.TRIAGE, SARS2,
              calibration=_real_cal(0.511, 0.488, 0.534))
    labels = [getattr(w, "label", "") for w in v.warnings]
    assert "VIRUS NOT VALIDATED" in labels, labels
    assert "CALIBRATION" not in labels


def test_a_refused_calibration_is_not_validated():
    """A protocol whose numbers must never become weights is refused by
    name in calibration.REFUSED_PROTOCOLS. Refused is not measured, so
    it cannot earn the CALIBRATION label -- and discriminating() is the
    single place that rule lives."""
    cal = _real_cal(0.675, 0.660, 0.690)
    cal.refused_reason = "single-screen weights are not transferable"
    assert not guards._is_validated(cal)
    labels = [getattr(w, "label", "")
              for w in check(QueryKind.TRIAGE, SARS2, calibration=cal).warnings]
    assert "VIRUS NOT VALIDATED" in labels, labels


# --------------------------- the prospective refusal's closing offer
#
# The generic offer promises "explain the mechanism of a compound with
# known activity, rank a set of compounds you already have". For Yellow
# fever and DENV-4 the agent refuses BOTH -- no metapath terminates at
# them. A refusal that declines the question and then offers two things
# it will also decline sends the reader straight into a second refusal.
#
# It became reachable on 2026-10-29: until the candidate-seeking
# backstop landed, "which compounds should be considered for Yellow
# fever" never classified PROSPECT, so nothing prospective had ever
# arrived here carrying a routeless virus.
def _routeless(label="Yellow fever", n_labels=1030):
    return VirusRoutes(label=label, n_labels=n_labels)


def _routed(label="SARS-CoV-2"):
    return VirusRoutes(label=label, n_labels=14996, inhibits=400,
                       crispr=4582, ppi=900)


def test_a_routeless_virus_is_not_offered_what_it_will_also_be_refused():
    v = check(QueryKind.PROSPECT, virus="NCBITaxon:11089",
              routes=_routeless())
    assert not v.allowed
    assert "rank a set of compounds you already have, and retrieve" \
        not in v.reason, "the generic offer was made for a routeless virus"
    assert "Yellow fever" in v.reason
    assert "1,030" in v.reason
    assert "no route to reach a compound at all" in v.reason
    assert "cannot explain a mechanism or rank a set" in v.reason
    # and it still says what DOES work, rather than only what does not
    assert "retrieving evidence for a named claim" in v.reason


def test_a_virus_with_routes_still_gets_the_general_offer():
    """The other direction. A narrower offer printed for every virus
    would understate what the graph does for the eight EVALUABLE ones,
    and would pass any test that only checked Yellow fever."""
    v = check(QueryKind.PROSPECT, virus="NCBITaxon:2697049",
              routes=_routed())
    assert not v.allowed
    assert "What this graph can do" in v.reason
    assert "no route to reach a compound at all" not in v.reason
    assert "Yellow fever" not in v.reason


def test_the_prospective_refusal_still_works_with_no_routes_supplied():
    """`routes` is optional throughout this module -- a caller that
    cannot compute it gets the old behaviour rather than a message
    built from a virus it knows nothing about. `None` must not be read
    as "no routes"."""
    v = check(QueryKind.PROSPECT, virus="NCBITaxon:11089", routes=None)
    assert not v.allowed
    assert "What this graph can do" in v.reason
    assert "no route to reach a compound at all" not in v.reason


def test_the_general_limit_is_still_stated_for_a_routeless_virus():
    """The virus-specific part is ADDED to the categorical refusal, not
    substituted for it. Prospective generation is unsupported for every
    virus in this graph, and a reader told only "this virus is thin"
    could reasonably conclude a richer virus would be answered."""
    v = check(QueryKind.PROSPECT, virus="NCBITaxon:11089",
              routes=_routeless())
    assert "cannot support prospective candidate generation" in v.reason
    assert "temporal split" in v.reason


# ------------------------------------------- the intent gate, audit A4 and A5
#
# Reproduced before any change, all eight exactly as the audit reported.
@pytest.mark.parametrize("q", [
    # A4, FIRST DIRECTION: the backstop ran LAST, so appending any
    # permitted keyword to a prospective question defeated it. One extra
    # word, and the answer came back under full caveats that make it
    # read as vetted.
    "which drugs could work against Zika and what is the mechanism?",
    "which compounds should we test against Yellow fever? explain why",
    "suggest compounds for JEV and rank them",
    "which inhibitors could be repurposed for yellow fever, which process",
    # A4, SECOND DIRECTION: no candidate noun at all, so the
    # noun-plus-seeking-word pair could not see them. All four fell to
    # UNKNOWN, and the orchestrator's catch-all routes UNKNOWN to
    # explain() -- a routed operation. That is the live path by which a
    # refused-tier virus got a mechanism walk attempted.
    "what might be repurposed for West Nile?",
    "which approved medicines could be repurposed for dengue?",
    "what should be screened against MERS-CoV?",
    "top 10 molecules for Yellow fever",
])
def test_a_permitted_keyword_no_longer_defeats_the_backstop(q):
    assert classify_intent(q) is QueryKind.PROSPECT, q


@pytest.mark.parametrize("q", [
    # A5: bare substring matching. "any" inside *company*, "try" inside
    # *country*, "list" inside *listing*. On launch day a refusal here
    # reads as the tool being broken.
    "what is the company behind this compound",
    "which country reported this drug first",
    "show the listing for this compound",
    # And questions about what the graph HOLDS, which are a census and
    # not a prediction. These are the ones the reordering would have
    # broken without carve-outs.
    "which compounds are in this graph",
    "how many drugs does the graph hold for Zika",
    "what antivirals are stored for DENV-2",
])
def test_legitimate_questions_are_not_read_as_requests_for_a_shortlist(q):
    assert classify_intent(q) is not QueryKind.PROSPECT, q


def test_screen_is_matched_as_a_phrase_not_as_a_bare_verb():
    """"screened against X" is a request for candidates. "screen" alone
    is not: this project asks about CRISPR screens constantly, and a
    bare verb on the list would refuse its own vocabulary."""
    assert classify_intent(
        "what should be screened against MERS-CoV") is QueryKind.PROSPECT
    assert classify_intent(
        "which CRISPR screens exist for Zika") is not QueryKind.PROSPECT
    assert classify_intent(
        "how many screens does BioGRID ORCS have for SARS-CoV-2"
    ) is not QueryKind.PROSPECT


def test_a_compound_in_hand_is_not_a_request_for_compounds():
    """"this compound" pairs the seeking word "what" with the candidate
    noun "compound" and asks about one the user already has. The comment
    on _SELF_SEEKING always said so and nothing acted on it."""
    assert classify_intent(
        "what is the company behind this compound") is not QueryKind.PROSPECT
    assert classify_intent(
        "what is the structure of that drug") is not QueryKind.PROSPECT


def test_a_plural_noun_overrides_the_in_hand_carve_out():
    """FAILING OPEN IS THE DIRECTION THAT MATTERS. Without this, "which
    other compounds are like this compound" is carved out by the
    determiner and answered."""
    assert classify_intent(
        "which other compounds are like this compound") is QueryKind.PROSPECT
    assert classify_intent(
        "what drugs resemble this drug") is QueryKind.PROSPECT


def test_a_prospective_verb_beats_the_in_hand_carve_out():
    """Ordering inside _seeks_candidates: "repurposed for" is the
    request whatever noun it points at."""
    assert classify_intent(
        "could this compound be repurposed for dengue") is QueryKind.PROSPECT


def test_the_word_matcher_is_the_one_the_orchestrator_already_had():
    """_mentions handled this correctly for virus and compound names and
    was simply not used here. Both must agree, or a name is a mention in
    one layer and not in the other."""
    from kgav.agent.guards import _word
    from kgav.agent.orchestrator import _mentions

    for text, name in [("the company filed", "any"), ("in this country", "try"),
                       ("a listing of", "list"), ("what is known", "what"),
                       ("hydroxychloroquine", "chloroquine"),
                       ("sars-cov-2.", "sars-cov-2")]:
        assert _word(text, name) == _mentions(text, name), (text, name)


def test_the_plural_noun_list_is_derived_from_the_noun_list():
    """A second hand-written list of the same nouns is the drift this
    project pays for repeatedly."""
    from kgav.agent.guards import _CANDIDATE_NOUNS, _PLURAL_CANDIDATE_NOUNS

    assert set(_PLURAL_CANDIDATE_NOUNS) <= set(_CANDIDATE_NOUNS)
    assert all(n.endswith("s") for n in _PLURAL_CANDIDATE_NOUNS)
    assert "compounds" in _PLURAL_CANDIDATE_NOUNS
    assert "compound" not in _PLURAL_CANDIDATE_NOUNS


@pytest.mark.parametrize("q", [
    # THE SUBSTRINGS, WITH NO OTHER CARVE-OUT IN REACH. The cases above
    # all contain "this compound" or "this drug" and are held by the
    # in-hand rule, so they passed with substring matching restored and
    # two mutants survived. These carry a candidate noun, no determiner,
    # and no real seeking word -- only one hiding inside another word.
    "the company sells antivirals",          # "any" inside company
    "a listing of drugs follows",            # "list" inside listing
    "this country produces generic drugs",   # "try" inside country
    # "shortlist" inside "shortlisted", which is the one self-seeking
    # noun that occurs inside a longer English word -- and it means the
    # opposite: a shortlisted compound is one already in hand.
    "the shortlisted compounds we already have",
])
def test_a_seeking_word_inside_a_longer_word_is_not_a_seeking_word(q):
    from kgav.agent.guards import _seeks_candidates

    assert _seeks_candidates(q) is False, q
    assert classify_intent(q) is not QueryKind.PROSPECT, q


# ------------------------------------- the route that reaches no label
# THE THIRD STATE, added 2026-11-02. docs/kg_audit_report.txt section 9
# prints four per-virus tiers and says "THE AGENT ENFORCES THIS". It
# enforced three: a virus with routes that reach NONE of its held-out
# positives was answered exactly like one whose routes reach plenty,
# with a caveat saying its reliability was "unknown" -- where it is
# known and it is zero. Five viruses on v0.1: DENV-1, DENV-3,
# HCoV-NL63, JEV, TBEV.
#
# DENV-1's figures throughout: 349 held-out positives, 53 virus-host
# interactions, reach 0. It is deliberately NOT thin -- 53 clears
# THIN_ROUTE of 25 -- because the virus that looked most answerable in
# this tier is the one worth testing.
def _unreachable(label="DENV-1", n_labels=349, ppi=53, reached=0):
    return VirusRoutes(label=label, n_labels=n_labels, ppi=ppi,
                       reached_labels=reached)


def _labels_of(verdict):
    return [w.label for w in verdict.warnings]


@pytest.mark.parametrize("kind", sorted(ROUTED_KINDS, key=lambda k: k.value))
def test_a_route_reaching_no_label_is_caveated_on_every_routed_kind(kind):
    v = check(kind, virus="NCBITaxon:11053", routes=_unreachable())
    assert v.allowed, "the ranking exists; it is the evaluation that does not"
    assert "LABELS UNREACHABLE" in _labels_of(v), _labels_of(v)


def test_the_caveat_names_the_labels_and_the_routes():
    v = check(QueryKind.TRIAGE, virus="NCBITaxon:11053",
              routes=_unreachable())
    text = next(w for w in v.warnings if w.label == "LABELS UNREACHABLE")
    assert "349" in text
    assert "53" in text
    assert "DENV-1" in text
    assert "reach NONE of them" in text
    assert "mark it right or wrong" in text


def test_unmeasured_reach_claims_nothing():
    """None is not 0. A missing or stale metapath_comparison.json must
    leave reach UNKNOWN, not assert that no route reaches a label --
    the opposite-claims mistake labels.measured_reach exists to avoid,
    and it would be undone here by one truthiness test."""
    r = _unreachable(reached=None)
    assert not r.reach_measured
    assert not r.unreachable
    v = check(QueryKind.TRIAGE, virus="NCBITaxon:11053", routes=r)
    assert "LABELS UNREACHABLE" not in _labels_of(v), \
        "unmeasured reach was reported as zero reach"


def test_a_reachable_route_is_not_caveated():
    r = _unreachable(label="DENV-2", n_labels=400, ppi=200, reached=117)
    assert r.reach_measured and not r.unreachable
    v = check(QueryKind.TRIAGE, virus="NCBITaxon:11060", routes=r)
    assert "LABELS UNREACHABLE" not in _labels_of(v)


def test_reach_of_one_is_enough_to_clear_the_caveat():
    """The boundary. The tier is reach == 0, not reach below some
    floor: MIN_REACHED_POSITIVES is a readability gate on the
    evaluation and a different question from whether anything at all
    was reached."""
    assert not _unreachable(reached=1).unreachable
    assert _unreachable(reached=0).unreachable


def test_a_virus_with_no_label_is_not_unreachable():
    """HCoV-HKU1: routes exist, no label exists -- FEATURES, NOTHING TO
    SCORE. Reach is 0 because there is nothing to reach, and saying its
    routes "reach none of its 0 measurements" would be true and
    meaningless."""
    r = VirusRoutes(label="HCoV-HKU1", n_labels=0, ppi=90,
                    reached_labels=0)
    assert not r.unreachable
    v = check(QueryKind.TRIAGE, virus="NCBITaxon:290028", routes=r)
    assert "LABELS UNREACHABLE" not in _labels_of(v)


def test_a_routeless_virus_is_refused_rather_than_caveated():
    """DENV-4 holds 312 labels and NO route, so no metapath terminates
    at it: the LABELS, NO FEATURE ROUTE tier, which check() refuses
    outright. The external audit listed DENV-4 as a sixth virus in the
    reaches-no-label tier. It is not one -- a refusal is already the
    strongest possible answer, and adding a caveat about reach to a
    reply that was never given would be the weaker outcome."""
    r = VirusRoutes(label="DENV-4", n_labels=312, reached_labels=0)
    assert not r.unreachable
    v = check(QueryKind.TRIAGE, virus="NCBITaxon:11070", routes=r)
    assert not v.allowed
    assert "LABELS UNREACHABLE" not in _labels_of(v)


def test_thin_and_unreachable_are_both_said_reach_first():
    """JEV: 46 labels through 7 virus-host interactions, reach 0. Both
    are true and they are different claims -- one that the route is
    small, one that it lands nowhere a measurement covers. Reach is
    stated first because it subsumes the reason a reader would care
    that the route is small."""
    r = VirusRoutes(label="JEV", n_labels=46, ppi=7, reached_labels=0)
    assert r.thin and r.unreachable
    said = _labels_of(check(QueryKind.TRIAGE, virus="NCBITaxon:11072",
                            routes=r))
    assert "LABELS UNREACHABLE" in said and "THIN ROUTE" in said
    assert said.index("LABELS UNREACHABLE") < said.index("THIN ROUTE")


def test_unreachable_is_not_the_same_claim_as_thin():
    """DENV-1's 53 edges clear THIN_ROUTE, so before this tier existed
    the virus whose ranking is least checkable carried no caveat about
    its route at all."""
    r = _unreachable()
    assert not r.thin, "DENV-1's 53 PPI edges are above THIN_ROUTE"
    assert r.unreachable


def test_the_caveat_key_is_registered():
    """caveat() raises on an unknown key, so a label added to the
    emitter and not to CAVEAT_LABELS fails loudly. Asserted because
    the verifier enforces the LABEL and nothing else."""
    assert CAVEAT_LABELS["labels_unreachable"] == "LABELS UNREACHABLE"


# ---------------------- the readability threshold the USER is shown
# S3's other half. pathway_caveat quoted the module constant -- 3,
# derived from the uncorrected 0.05 -- in a sentence that names the
# CORRECTED alpha two clauses earlier. So the user was shown both
# thresholds and told the wrong one gated the ratio. A mutation
# restoring the constant survived every test until this one.
def test_the_pathway_caveat_quotes_the_rows_own_threshold():
    from kgav.agent.guards import pathway_caveat
    from kgav.landscape import Landscape, TaxonOverlap

    V = "NCBITaxon:11060"
    # DENV-2's real figures: 50 dependency factors, 6.22 expected, 2
    # observed, lower-tail p 0.045 -- nominally interesting, nowhere
    # near a corrected 0.001136, and its ratio is NOT readable downward
    # at that alpha even though it is at 0.05.
    row = TaxonOverlap(virus=V, n_dependency=50, n_observed=2,
                       expected=6.22, p_naive=0.02, p_matched=0.9,
                       p_matched_less=0.045, alpha=0.05 / 44, n_tests=44)
    land = Landscape(source_file="pathway_landscape.json", release="v0.1",
                     permutations=2000, seed=0, taxa={V: row},
                     other_channel_tests=28)
    said = pathway_caveat(land, V)
    assert not row.ratio_readable
    assert "6.8" in said, said
    assert " 3," not in said and "below 3" not in said, (
        "the caveat still quotes the uncorrected threshold")


# --------------------- A7: the refusal's figure has to reach the refusal
# guards.check has taken a `temporal` argument all along, and nothing passed
# it: the orchestrator supplied `calibration` and `routes` only, so the
# branch every user reached was the one with no measurement in it. Both ends
# are tested here -- the sentence, and the wiring that feeds it.


def test_the_tools_expose_the_prospective_calibration_separately(tmp_path):
    """Keyed by the same rule as the cross-sectional one, through the same
    indexer, so the two cannot be looked up differently."""
    import json
    from kgav.agent.tools import GraphTools
    V = "NCBITaxon:2697049"
    nodes = [{"id": V, "class": "OrganismTaxon",
              "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                             "baltimore_class": "IV", "is_enveloped": True}}]
    (tmp_path / "nodes.jsonl").write_text(json.dumps(nodes[0]))
    (tmp_path / "edges.jsonl").write_text("")
    cross = _real_cal(0.577, 0.514, 0.639, protocol="cross-sectional")
    temp = _real_cal(0.502, 0.470, 0.534, protocol="post-2021")
    for cal in (cross, temp):
        cal.virus = V
    g = GraphTools(tmp_path, calibration={V: cross}, temporal={V: temp})
    assert g.calibration_for(V) is cross
    assert g.temporal_for(V) is temp
    assert g.temporal_for("NCBITaxon:11060") is None, (
        "a virus the prospective run did not cover must resolve to None, "
        "not to another virus's numbers")


def test_a_tools_object_given_no_prospective_run_has_none(tmp_path):
    import json
    from kgav.agent.tools import GraphTools
    (tmp_path / "nodes.jsonl").write_text(json.dumps(
        {"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}}))
    (tmp_path / "edges.jsonl").write_text("")
    g = GraphTools(tmp_path)
    assert g.temporal_for("NCBITaxon:2697049") is None


def test_the_refusal_states_which_scorer_and_which_protocol():
    """The measured branch, which no caller could reach before."""
    cal = _real_cal(0.510, 0.486, 0.533, protocol="post-2021")
    reason = check(QueryKind.PROSPECT, temporal=cal).reason
    assert "0.510" in reason and "post-2021" in reason
    assert "measured rather than assumed" in reason, (
        "the claim and the figure have to travel together; either without "
        "the other is what A7 was about")
