"""Gate for the agent's guardrails: the refusal is deterministic, and every
caveat the evaluation established is attached where a user will see it.
"""
from pathlib import Path

import pytest

from kgav.agent import guards
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


def test_the_refusal_quotes_no_figure_without_a_calibration():
    """Absent a measured run, the limit is stated without a number. A figure
    with no run behind it is what _channel_caveat and direct_acting_caveat
    were each rewritten to remove."""
    import re
    reason = check(QueryKind.PROSPECT).reason
    assert "every scorer performs at chance" in reason
    assert not re.search(r"AUC \d", reason), \
        f"quoted an AUC with no calibration supplied:\n{reason}"


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
    assert check(QueryKind.TRIAGE, SARS2).warnings == []


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
