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
    def __init__(self, auc, lo, hi, n_neg=7343, discriminates=False):
        self.auc, self.ci_lo, self.ci_hi, self.n_neg = auc, lo, hi, n_neg
        self.discriminates = discriminates


class _Cal:
    def __init__(self, channels, protocol="single-screen/publication-disjoint",
                 usable=True):
        self.channels, self.protocol, self.usable = channels, protocol, usable


HOST_PATH = [{"id": "INCHIKEY:X", "class": "SmallMolecule", "properties": {}},
             {"id": "UniProtKB:Q99720", "class": "Protein",
              "properties": {"gene_symbol": "SIGMAR1", "is_viral": False}},
             {"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
              "properties": {"label": "SARS-CoV-2"}}]


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
