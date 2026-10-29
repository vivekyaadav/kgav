"""The virus-centred pathway census.

pathways() asks "does this compound hit anything the virus depends on".
This asks the other direction -- "which pathways does this virus depend
on" -- which is probably the first question a reader puts to a graph
built for pathway analysis, and which the orchestrator answered with a
node profile until 2026-10-29 because its PATHWAY branch required a
compound.

The fixture is small and hand-built so every figure in these tests is
checkable by counting the fixture rather than by running the code and
recording what it said.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from kgav.agent import guards
from kgav.agent.orchestrator import Orchestrator
from kgav.agent.tools import GraphTools

V = "NCBITaxon:64320"          # Zika in the register
OTHER = "NCBITaxon:11060"      # DENV-2
DRUG = "INCHIKEY:AAAAA"


def _nodes(n_members: int) -> list[dict]:
    out = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "Zika", "family": "Flaviviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        {"id": OTHER, "class": "OrganismTaxon",
         "properties": {"label": "DENV-2", "family": "Flaviviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        {"id": DRUG, "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "AAAAA",
                        "label": "adrug", "chembl_id": "CHEMBL9",
                        "is_approved": True, "salt_collapsed": False,
                        "stereo_collapsed": False}},
        {"id": "REACT:BIG", "class": "Pathway",
         "properties": {"label": "Big pathway"}},
        {"id": "REACT:SMALL", "class": "Pathway",
         "properties": {"label": "Small pathway"}},
    ]
    out.append({"id": "UniProtKB:VP", "class": "Protein",
                "properties": {"taxon_id": V, "is_viral": True,
                               "sequence_hash": "h", "reviewed": True}})
    for i in range(n_members):
        out.append({"id": f"UniProtKB:P{i}", "class": "Protein",
                    "properties": {"taxon_id": "NCBITaxon:9606",
                                   "is_viral": False, "sequence_hash": "h",
                                   "reviewed": True,
                                   "gene_symbol": f"G{i}"}})
    return out


def _edge(s, p, o, q=None):
    return {"subject": s, "predicate": p, "object": o, "qualifiers": q or {},
            "primary_knowledge_source": "infores:biogrid-orcs",
            "evidence_tier": 1, "first_asserted_date": "2021-01-01"}


def _release(tmp_path: Path, *, deps: list[int], big: list[int],
             small: list[int], n_members: int = 12,
             druggable: list[int] | None = None,
             direction: str = "dependency", routable: bool = True) -> Path:
    """A release with two pathways of known size and known overlap.

    deps      proteins the virus depends on
    big/small which proteins are members of each pathway
    druggable which proteins carry a TARGETS edge
    routable  give the virus a virus-host PPI edge, so guards.check does
              not refuse the query before the census runs. Without it a
              virus whose only host factors are RESTRICTION factors has
              no feature route at all and is refused for that reason --
              true, and not the thing these tests are about. Real
              viruses in this register are in exactly that position:
              HCoV-NL63 carries restriction and PPI and no dependency.
    """
    d = tmp_path / "rel"
    d.mkdir()
    (d / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in _nodes(n_members)))
    edges = []
    if routable:
        edges.append(_edge("UniProtKB:VP", "PHYSICALLY_INTERACTS_WITH",
                           "UniProtKB:P0"))
    for i in deps:
        edges.append(_edge(f"UniProtKB:P{i}", "HOST_FACTOR_FOR", V,
                           {"direction": direction, "screen_type": "CRISPRko",
                            "cell_line": "CVCL_0574"}))
    for i in big:
        edges.append(_edge(f"UniProtKB:P{i}", "PARTICIPATES_IN", "REACT:BIG"))
    for i in small:
        edges.append(_edge(f"UniProtKB:P{i}", "PARTICIPATES_IN",
                           "REACT:SMALL"))
    for i in (druggable or []):
        edges.append(_edge(DRUG, "TARGETS", f"UniProtKB:P{i}",
                           {"direction": "inhibitor",
                            "assay_type": "binding"}))
    (d / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return d


def _tools(release: Path) -> GraphTools:
    return GraphTools(release)


# ------------------------------------------------------------- the census
def test_it_reports_the_pathways_a_virus_depends_on(tmp_path):
    rel = _release(tmp_path, deps=[0, 1, 2], big=[0, 1, 2, 3, 4, 5],
                   small=[2, 6])
    r = _tools(rel).virus_pathways(V)
    assert r.ok
    by_id = {row["pathway"]: row for row in r.data}
    assert set(by_id) == {"REACT:BIG", "REACT:SMALL"}
    assert by_id["REACT:BIG"]["n_members"] == 6
    assert by_id["REACT:BIG"]["n_dependency_members"] == 3
    assert by_id["REACT:SMALL"]["n_members"] == 2
    assert by_id["REACT:SMALL"]["n_dependency_members"] == 1


def test_a_restriction_factor_is_not_a_dependency(tmp_path):
    """The sign is the difference between a therapeutic target and one
    it would be harmful to inhibit, and this census must count only the
    same set the measured null was computed over."""
    rel = _release(tmp_path, deps=[0, 1], big=[0, 1, 2], small=[],
                   direction="restriction")
    r = _tools(rel).virus_pathways(V)
    assert not r.ok
    assert "none carrying direction=dependency" in r.note
    assert "restriction or undirected" in r.note


def test_a_virus_with_no_screen_at_all_says_which_absence_it_is(tmp_path):
    """Five of the sixteen viruses in this register have no dependency
    factor. "No pathways" would read as a statement about the virus;
    it is a statement about the screen layer."""
    rel = _release(tmp_path, deps=[], big=[0, 1], small=[])
    r = _tools(rel).virus_pathways(V)
    assert not r.ok
    assert "no HOST_FACTOR_FOR edge at all" in r.note
    assert "absence in the screen layer" in r.note
    assert "not a statement that this virus depends on nothing" in r.note


def test_druggable_members_are_counted(tmp_path):
    """The census is for repurposing, so "which of these can be reached
    by a compound in this graph" is the column that decides whether a
    row is actionable."""
    rel = _release(tmp_path, deps=[0, 1, 2], big=[0, 1, 2], small=[],
                   druggable=[1])
    row = _tools(rel).virus_pathways(V).data[0]
    assert row["n_druggable"] == 1
    assert row["druggable_members"] == ["UniProtKB:P1"]


def test_density_is_reported_against_pathway_size(tmp_path):
    """A 90-member pathway holds more of anything than a 4-member one,
    so a census ranked on the raw count ranks on size. Each row carries
    the count its size predicts."""
    rel = _release(tmp_path, deps=[0, 1], big=[0, 1, 2, 3, 4, 5],
                   small=[0, 1], n_members=10)
    rows = {r["pathway"]: r for r in _tools(rel).virus_pathways(V).data}
    # 2 dependency factors among 10 host proteins = 20% density.
    assert rows["REACT:BIG"]["expected_members"] == pytest.approx(6 * 0.2)
    assert rows["REACT:SMALL"]["expected_members"] == pytest.approx(2 * 0.2)
    # Same count, different size -> the smaller pathway is the denser.
    assert (rows["REACT:SMALL"]["density_ratio"]
            > rows["REACT:BIG"]["density_ratio"])


def test_it_claims_no_significance_for_a_row(tmp_path):
    """The set is chosen BY holding a dependency factor, so an
    enrichment tail would be circular: every row contains one with
    probability 1. The text must not imply otherwise."""
    rel = _release(tmp_path, deps=[0, 1], big=[0, 1, 2], small=[])
    text = " ".join(_tools(rel).virus_pathways(V).verbalised)
    assert "No significance is claimed for a row" in text
    assert "nothing survived correction" in text
    assert "lead to examine" in text


def test_dependency_factors_in_no_pathway_are_reported(tmp_path):
    """Coverage, stated rather than implied. A census over 3 of 40
    factors reads identically to one over 40 of 40 unless it says so."""
    rel = _release(tmp_path, deps=[0, 1, 2, 3], big=[0, 1], small=[])
    text = " ".join(_tools(rel).virus_pathways(V).verbalised)
    assert "2 of the dependency factors sit in no pathway" in text
    assert "covers 50%" in text


def test_factors_with_no_pathway_annotation_at_all(tmp_path):
    """Distinct from the virus having no screen: the screen found
    something and the annotation did not place it."""
    rel = _release(tmp_path, deps=[7, 8], big=[0, 1], small=[2])
    r = _tools(rel).virus_pathways(V)
    assert not r.ok
    assert "gap in pathway annotation, not in the screen" in r.note


def test_the_cap_excludes_and_says_so(tmp_path):
    """The same cap the landscape applied, for the reason it applied
    it: a null computed over capped pathways beside an uncapped census
    would be two questions under one answer."""
    rel = _release(tmp_path, deps=[0, 1], big=list(range(12)), small=[0, 1],
                   n_members=12)
    t = _tools(rel)
    t._max_pathway_members = 5
    r = t.virus_pathways(V)
    assert {row["pathway"] for row in r.data} == {"REACT:SMALL"}
    assert "excluded for exceeding" in " ".join(r.verbalised)


# ------------------------------------------------------------- the guard
def test_a_virus_with_no_feature_route_is_refused(tmp_path):
    """The no-route refusal applies here as it does to any routed
    query: a census for a virus nothing reaches is not an empty
    result, it is a question that cannot be asked of this release."""
    rel = _release(tmp_path, deps=[], big=[], small=[], routable=False)
    r = _tools(rel).virus_pathways(V)
    assert not r.ok
    assert "none can be reached" in r.note, (
        "a virus nothing reaches must get the no-route refusal, not an "
        "empty census")


# ------------------------------------------------------- the orchestrator
def test_the_question_reaches_the_tool(tmp_path):
    """THE DEFECT THIS CLOSES. The PATHWAY branch required a compound,
    so a question naming only a virus fell through every branch to a
    node profile -- an answer about the taxon node, to a question about
    pathways, with nothing saying the question was not understood."""
    rel = _release(tmp_path, deps=[0, 1], big=[0, 1, 2], small=[],
                   druggable=[0])
    orch = Orchestrator(GraphTools(rel))
    brief = orch.build_brief("which pathways does Zika depend on?")
    assert brief.allowed
    # THE TRACE, NOT A DISJUNCTION. The first version of this asserted
    # `"virus_pathways" in kinds or "pathway" in text`, reading kinds
    # off a `results` attribute Brief does not have -- so that half was
    # always empty and the test passed entirely on the word "pathway"
    # appearing somewhere in the facts, which it would have done for a
    # node profile too. The trace records which tool ran.
    assert "tool=virus_pathways ok=True" in brief.trace, brief.trace
    assert "Zika depends on 2 host protein(s)" in " ".join(brief.facts)


def test_a_compound_question_still_goes_to_the_compound_tool(tmp_path):
    """The new branch must not capture the question the old one
    answered."""
    rel = _release(tmp_path, deps=[0, 1], big=[0, 1, 2], small=[],
                   druggable=[0])
    orch = Orchestrator(GraphTools(rel))
    brief = orch.build_brief("which pathways does adrug act through "
                             "against Zika?")
    assert "tool=virus_pathways ok=True" not in brief.trace, brief.trace
    assert any(s.startswith("tool=pathways") for s in brief.trace), brief.trace


def test_the_intent_classifier_routes_the_phrasing():
    for q in ("which pathways does Zika depend on?",
              "what reactome pathways does zika rely on",
              "which process does Zika depend on"):
        assert guards.classify_intent(q) is guards.QueryKind.PATHWAY, q


# ------------------------------- the ratio column must select its own list
def test_a_dense_small_pathway_reaches_the_second_list(tmp_path):
    """THE DEFECT THIS CLOSES. The census sorted by COUNT, which sorts
    by SIZE -- a 200-member pathway holds more of everything. On the
    live SARS-CoV-2 run every one of the 25 rows shown was large, and
    the pathways the ratio column would have chosen were not on the
    page. Here BIG holds more factors and SMALL is far denser."""
    rel = _release(tmp_path, deps=[0, 1, 2, 3, 4, 5, 6],
                   big=list(range(12)), small=[0, 1, 2, 3],
                   n_members=12)
    t = _tools(rel)
    t._max_pathway_members = 200
    text = " ".join(t.virus_pathways(V).verbalised)
    assert "BY DENSITY instead" in text
    head, tail = text.split("BY DENSITY instead", 1)
    # SMALL: 4 of 4 members are factors, expects 4 * 7/12 = 2.33 -> 1.71
    # BIG:   7 of 12,            expects 12 * 7/12 = 7.0  -> 1.00
    #
    # SCOPED TO EACH LIST. Both pathways appear in both, so a search
    # over the whole text finds whichever is printed first and says
    # nothing about either ordering.
    assert "Small pathway" in tail and "Big pathway" in tail
    assert tail.index("Small pathway") < tail.index("Big pathway"), (
        "the second list is not ordered by density")
    assert head.index("Big pathway") < head.index("Small pathway"), (
        "the first list is not ordered by count")


def test_a_pathway_expecting_under_one_factor_is_not_ranked_by_ratio(tmp_path):
    """A 2-member pathway holding both factors reads as a huge ratio and
    means nothing. The floor is expectation >= 1, which is this
    project's existing rule for when a lift is readable rather than a
    threshold invented for this tool."""
    rel = _release(tmp_path, deps=[0, 1], big=[0, 1, 2, 3, 4, 5],
                   small=[0, 1], n_members=12)
    t = _tools(rel)
    t._max_pathway_members = 200
    r = t.virus_pathways(V)
    rows = {x["pathway"]: x for x in r.data}
    # SMALL expects 2 * 2/12 = 0.33 -- under one.
    assert rows["REACT:SMALL"]["expected_members"] < 1.0
    text = " ".join(r.verbalised)
    tail = text[text.index("BY DENSITY instead"):] if "BY DENSITY" in text \
        else ""
    assert "Small pathway" not in tail, (
        "a pathway expecting under one factor was ranked by its ratio")
    assert "excluded from this second list" in text


def test_the_second_list_says_how_many_it_excluded(tmp_path):
    rel = _release(tmp_path, deps=[0, 1], big=[0, 1, 2, 3, 4, 5],
                   small=[0, 1], n_members=12)
    t = _tools(rel)
    t._max_pathway_members = 200
    text = " ".join(t.virus_pathways(V).verbalised)
    assert "1 pathway(s) are excluded from this second list" in text


# ------------- a pathway curated from this virus's own literature
def _named_release(tmp_path, pathway_label):
    """A release whose one pathway carries the given label."""
    import json as _json
    rel = _release(tmp_path, deps=[0, 1, 2], big=[0, 1, 2, 3, 4, 5],
                   small=[])
    nodes = [_json.loads(ln) for ln
             in (rel / "nodes.jsonl").read_text().splitlines() if ln.strip()]
    for n in nodes:
        if n["id"] == "REACT:BIG":
            n["properties"]["label"] = pathway_label
    (rel / "nodes.jsonl").write_text(
        "\n".join(_json.dumps(n) for n in nodes))
    return rel


def test_a_pathway_naming_the_virus_is_flagged(tmp_path):
    """THE CONFOUND THIS SURFACES. Reactome carries virus-specific
    pathways, and the live SARS-CoV-2 census put "SARS-CoV-2 modulates
    autophagy" fourth by density. A screen for a virus's dependency
    factors overlapping a pathway curated FROM that virus's literature
    is partly definitional -- the shared-publication confound this
    project chases elsewhere, in the pathway vocabulary this time."""
    rel = _named_release(tmp_path, "Zika modulates autophagy")
    r = _tools(rel).virus_pathways(V)
    row = next(x for x in r.data if x["pathway"] == "REACT:BIG")
    assert row["names_the_virus"] is True
    # SCOPED TO THE ROW. The explanation paragraph also contains the
    # literal "[NAMES THIS VIRUS]", so searching the whole answer passed
    # with the per-row mark deleted -- the mark is what a reader sees
    # against the pathway, and the paragraph only explains it.
    row_line = next(ln for ln in r.verbalised
                    if "Zika modulates autophagy" in ln
                    and "ratio" in ln)
    assert "[NAMES THIS VIRUS]" in row_line, row_line
    text = " ".join(r.verbalised)
    assert "partly definitional" in text
    assert "shown rather than dropped" in text


def test_an_ordinary_pathway_is_not_flagged(tmp_path):
    """The other side, so the mark is not simply always printed."""
    rel = _named_release(tmp_path, "clathrin-mediated endocytosis")
    r = _tools(rel).virus_pathways(V)
    row = next(x for x in r.data if x["pathway"] == "REACT:BIG")
    assert row["names_the_virus"] is False
    assert "[NAMES THIS VIRUS]" not in " ".join(r.verbalised)


def test_the_flag_follows_the_register_not_one_spelling(tmp_path):
    """The register gives each virus aliases, and a pathway label may
    use any of them. Matching only the display label would miss
    "SARS-CoV-2 infection" for a virus the graph labels differently."""
    from kgav.agent.tools import GraphTools
    # THE ALIAS MUST NOT BE A SUPERSTRING OF THE DISPLAY LABEL, or the
    # test passes under a display-label-only implementation: "zika" is
    # inside "zika virus infection". ZIKV is the real abbreviation and
    # shares no substring with "Zika" long enough to match.
    rel = _named_release(tmp_path, "zikv infection")
    t = GraphTools(rel)
    # Read from the graph, not asserted between two literals: the point
    # is that the DISPLAY label does not occur in this pathway's label,
    # so only the alias can produce the flag.
    label = t._label(V).lower()
    pathway = t.nodes["REACT:BIG"]["properties"]["label"].lower()
    assert label not in pathway, (label, pathway)
    t.names["zikv"] = V
    r = t.virus_pathways(V)
    row = next(x for x in r.data if x["pathway"] == "REACT:BIG")
    assert row["names_the_virus"] is True


def test_the_flag_is_marked_in_both_lists(tmp_path):
    """It was only in the by-density list, so a reader of the one
    printed first -- the one most likely to be read at all -- saw the
    pathway without the thing that qualifies it."""
    rel = _named_release(tmp_path, "Zika modulates autophagy")
    r = _tools(rel).virus_pathways(V)
    marked = [ln for ln in r.verbalised
              if "Zika modulates autophagy" in ln and "ratio" in ln]
    assert len(marked) == 2, f"expected the row in both lists: {marked}"
    for ln in marked:
        assert "[NAMES THIS VIRUS]" in ln, ln
