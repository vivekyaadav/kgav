"""A5: which pathways a compound's targets sit in, and what that is not.

The question had no answer before this. _find_paths walks M1 and M2 only, so
the pathway metapath the schema declares as M4 is never traversed, and an
agent asked about pathway membership returned a silence that reads like a
negative finding.

These tests fix what the answer may and may not say. That it is a CENSUS --
every pathway counted, every exclusion reported, every target that
contributed none named -- rather than a sample, because a path sampler tells
a reader nothing about coverage. That it applies the same member cap the
measured null was computed over. And above all that membership is never
presented as benefit: most TARGETS edges in this graph record potency with no
direction, so the compound may activate a protein the virus depends on as
readily as inhibit it.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from kgav.agent.tools import GraphTools

V = "NCBITaxon:2697049"
DRUG = "INCHIKEY:DRUGA"
OTHERDRUG = "INCHIKEY:DRUGB"
DIRECT = "INCHIKEY:DRUGC"
ROGUE = "INCHIKEY:DRUGD"


def _node(nid, klass, **props):
    return {"id": nid, "class": klass, "properties": props}


def _edge(s, p, o, **q):
    return {"subject": s, "predicate": p, "object": o,
            "qualifiers": q, "provenance": {"source": "test"}}


@pytest.fixture
def release(tmp_path) -> Path:
    """A compound whose two targets sit in three pathways.

    REACT:DEP     holds a dependency factor         -> the interesting case
    REACT:PLAIN   holds none                        -> reported, not hidden
    REACT:BIG     exceeds the cap                   -> excluded, and said so
    UniProtKB:P3  a target in no pathway at all     -> named
    """
    nodes = [
        _node(DRUG, "SmallMolecule", smiles="C", inchikey_skel="DRUGA",
              chembl_id="CHEMBL1"),
        _node(OTHERDRUG, "SmallMolecule", smiles="C", inchikey_skel="DRUGB",
              chembl_id="CHEMBL2"),
        # Remdesivir's shape: a well-annotated VIRAL target and no human
        # one. Not an unannotated compound, and must not be described as
        # one.
        _node(DIRECT, "SmallMolecule", smiles="C", inchikey_skel="DRUGC",
              chembl_id="CHEMBL3"),
        # A release that violates its own schema: INHIBITS carries
        # object_constraints {is_viral: true}, and this one does not.
        _node(ROGUE, "SmallMolecule", smiles="C", inchikey_skel="DRUGD",
              chembl_id="CHEMBL4"),
        _node(V, "OrganismTaxon", name="SARS-CoV-2"),
        _node("UniProtKB:P1", "Protein", gene_symbol="P1", is_viral=False,
              taxon_id="NCBITaxon:9606", sequence_hash="h", reviewed=True),
        _node("UniProtKB:P2", "Protein", gene_symbol="P2", is_viral=False,
              taxon_id="NCBITaxon:9606", sequence_hash="h", reviewed=True),
        _node("UniProtKB:P3", "Protein", gene_symbol="P3", is_viral=False,
              taxon_id="NCBITaxon:9606", sequence_hash="h", reviewed=True),
        _node("UniProtKB:DEPF", "Protein", gene_symbol="DEPF", is_viral=False,
              taxon_id="NCBITaxon:9606", sequence_hash="h", reviewed=True),
        _node("UniProtKB:RESF", "Protein", gene_symbol="RESF", is_viral=False,
              taxon_id="NCBITaxon:9606", sequence_hash="h", reviewed=True),
        _node("UniProtKB:VIRAL", "Protein", protein_family="nsp5",
              is_viral=True, taxon_id=V, sequence_hash="h", reviewed=True),
        _node("REACT:DEP", "Pathway", label="Dependency-bearing pathway"),
        _node("REACT:PLAIN", "Pathway", label="Plain pathway"),
        _node("REACT:AAA", "Pathway", label="Busy pathway, no dependency"),
        _node("REACT:BIG", "Pathway", label="Oversized pathway"),
    ]
    edges = [
        # One directed target edge, one undirected: the caveat must count.
        _edge(DRUG, "TARGETS", "UniProtKB:P1", value_nm=50.0,
              direction="antagonist"),
        _edge(DRUG, "TARGETS", "UniProtKB:P2", value_nm=90.0),
        _edge(DRUG, "TARGETS", "UniProtKB:P3", value_nm=120.0),
        # A viral target must not enter a human-pathway census.
        _edge(DRUG, "TARGETS", "UniProtKB:VIRAL", value_nm=10.0),
        _edge("UniProtKB:P1", "PARTICIPATES_IN", "REACT:DEP"),
        _edge("UniProtKB:DEPF", "PARTICIPATES_IN", "REACT:DEP"),
        _edge("UniProtKB:RESF", "PARTICIPATES_IN", "REACT:DEP"),
        _edge("UniProtKB:P2", "PARTICIPATES_IN", "REACT:PLAIN"),
        # Two of the compound's targets and no dependency factor. Sorts
        # FIRST on target count and FIRST alphabetically, so only an
        # ordering that puts dependency hits above both keeps REACT:DEP.
        _edge("UniProtKB:P1", "PARTICIPATES_IN", "REACT:AAA"),
        _edge("UniProtKB:P2", "PARTICIPATES_IN", "REACT:AAA"),
        _edge("UniProtKB:P1", "PARTICIPATES_IN", "REACT:BIG"),
        # DEPF is depended on; RESF is a restriction factor, the opposite
        # claim, and must never be counted as one.
        _edge(DIRECT, "INHIBITS", "UniProtKB:VIRAL", value_nm=3.0),
        _edge(ROGUE, "INHIBITS", "UniProtKB:DEPF", value_nm=3.0),
        _edge("UniProtKB:DEPF", "HOST_FACTOR_FOR", V, direction="dependency"),
        _edge("UniProtKB:RESF", "HOST_FACTOR_FOR", V, direction="restriction"),
    ]
    # Pad REACT:BIG past any sane cap.
    for i in range(400):
        pid = f"UniProtKB:BIG{i}"
        nodes.append(_node(pid, "Protein", gene_symbol=f"BIG{i}",
                           is_viral=False, taxon_id="NCBITaxon:9606",
                           sequence_hash="h", reviewed=True))
        edges.append(_edge(pid, "PARTICIPATES_IN", "REACT:BIG"))
    # The agent resolves typed names through the release's prebuilt index,
    # so a fixture without one cannot exercise the orchestrator at all.
    (tmp_path / "NAMES.json").write_text(json.dumps({
        "druga": DRUG, "drugb": OTHERDRUG, "drugc": DIRECT, "drugd": ROGUE,
        "sars-cov-2": V}))
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in edges))
    return tmp_path


@pytest.fixture
def tools(release):
    return GraphTools(release)


# ------------------------------------------------------------- the census
def test_it_answers_a_question_explain_cannot(tools):
    """_find_paths walks M1 and M2 only; M4 is never traversed.

    The pathway question therefore returned nothing at all before this, and
    nothing is indistinguishable from "no".
    """
    assert tools._find_paths(DRUG, V, 20) == [] or all(
        r.startswith(("M1", "M2")) for r, _ in tools._find_paths(DRUG, V, 20))
    r = tools.pathways(DRUG, V)
    assert r.ok and r.data


def test_it_finds_the_pathway_with_a_dependency_factor(tools):
    r = tools.pathways(DRUG, V)
    dep = [row for row in r.data if row["pathway"] == "REACT:DEP"]
    assert len(dep) == 1
    assert dep[0]["dependency_members"] == ["UniProtKB:DEPF"]


def test_a_restriction_factor_is_not_a_dependency_factor(tools):
    """The opposite claim. Pooling the two makes every count meaningless."""
    r = tools.pathways(DRUG, V)
    for row in r.data:
        assert "UniProtKB:RESF" not in row["dependency_members"]


def test_a_pathway_with_no_dependency_factor_is_reported_not_dropped(tools):
    r = tools.pathways(DRUG, V)
    plain = [row for row in r.data if row["pathway"] == "REACT:PLAIN"]
    assert len(plain) == 1 and plain[0]["n_dependency_members"] == 0


def test_a_target_in_no_pathway_is_named(tools):
    """An absence a reader cannot otherwise detect."""
    r = tools.pathways(DRUG, V)
    assert any("P3" in line and "no pathway" in line
               for line in r.verbalised)


def test_a_viral_target_does_not_enter_a_human_pathway_census(tools):
    r = tools.pathways(DRUG, V)
    assert "3 human target(s)" in r.verbalised[0]
    for row in r.data:
        assert "UniProtKB:VIRAL" not in row["targets_in_pathway"]


# ------------------------------------------------------------------- caps
def test_an_oversized_pathway_is_excluded_and_said_so(tools):
    """Silently dropping it would describe a different graph from the one
    the measured null was computed over."""
    r = tools.pathways(DRUG, V)
    assert all(row["pathway"] != "REACT:BIG" for row in r.data)
    assert any("EXCLUDED" in line and "cap" in line for line in r.verbalised)
    assert "over the cap" in r.note


def test_the_cap_comes_from_the_schema(tools):
    from kgav.schema import load_schema
    assert tools._max_pathway_members == \
        load_schema().retrieval_limits.get("max_pathway_members")


def test_the_overflow_says_how_many_hits_it_hides(release):
    """A truncated list that does not say what it truncated is the defect
    this repository keeps finding."""
    t = GraphTools(release)
    r = t.pathways(DRUG, V, limit=1)
    assert any("further pathway(s) not shown" in line
               for line in r.verbalised)


def test_pathways_with_dependency_factors_survive_the_limit(release):
    """Ordering decides what the limit keeps, so the hits must come first.

    Asserted on the VERBALISED lines, which is what the model is handed.
    The first version of this test read r.data -- which is never truncated
    -- so it passed with the ordering removed entirely, and the fixture
    tie-broke alphabetically in the right direction by luck. REACT:AAA now
    carries two of the compound's targets and no dependency factor, so it
    wins every ordering except the correct one.
    """
    t = GraphTools(release)
    r = t.pathways(DRUG, V, limit=1)
    shown = [ln for ln in r.verbalised
             if ln.startswith(("  REACT", "  Dependency", "  Busy",
                               "  Plain"))]
    assert len(shown) == 1
    assert "Dependency-bearing pathway" in shown[0]
    assert "Busy pathway" not in " ".join(shown)


# ------------------------------------------------------ what it may not imply
def test_it_says_the_direction_of_most_target_edges_is_unknown(tools):
    """Membership is not benefit.

    Inhibiting a dependency factor should impede the virus; activating one
    should not, and this graph mostly does not record which.
    """
    r = tools.pathways(DRUG, V)
    labels = [getattr(w, "label", "") for w in r.warnings]
    assert "TARGET DIRECTION" in labels
    text = " ".join(str(w) for w in r.warnings)
    assert "2 of 3 target edge(s) here record potency without a direction" in text
    assert "hypothesis to test, not a reason to expect benefit" in text


def test_the_direction_count_comes_from_the_edges_actually_used(release):
    """Not quoted from the release, so it cannot drift from the answer."""
    t = GraphTools(release)
    text = " ".join(str(w) for w in t.pathways(DRUG, V).warnings)
    assert "2 of 3" in text and "2,226" not in text


def test_it_carries_the_measured_pathway_null(tools):
    """The route this tool reports is the one the landscape measured."""
    labels = [getattr(w, "label", "")
              for w in tools.pathways(DRUG, V).warnings]
    assert "PATHWAY ROUTE" in labels


def test_the_null_is_carried_even_with_no_landscape(tools):
    """Unmeasured and measured-null must not look the same from outside."""
    text = " ".join(str(w) for w in tools.pathways(DRUG, V).warnings)
    assert "has not been measured on this release" in text


# ------------------------------------------------------------- absences
def test_a_compound_not_in_the_graph_is_reported_as_absent(tools):
    r = tools.pathways("INCHIKEY:NOPE", V)
    assert not r.ok and "not in this graph" in r.note


def test_a_compound_with_no_human_target_says_why(tools):
    """A property of the source, not a finding about the compound."""
    r = tools.pathways(OTHERDRUG, V)
    assert not r.ok
    assert "no human target" in r.note
    assert "property of the source" in r.note


def test_an_absent_compound_still_carries_the_caveats(tools):
    """A refusal that drops them teaches the reader the caveats are
    decoration."""
    labels = [getattr(w, "label", "")
              for w in tools.pathways("INCHIKEY:NOPE", V).warnings]
    assert "PATHWAY ROUTE" in labels


def test_no_pathway_hits_is_stated_as_an_absence_in_the_crispr_layer(release):
    """A virus nothing was screened against must not read as a compound
    with no pathway relevance."""
    t = GraphTools(release)
    r = t.pathways(DRUG, "NCBITaxon:11082")
    assert r.ok
    assert any("absence in the CRISPR layer" in line for line in r.verbalised)


def test_every_reported_edge_is_citable(tools):
    """Each id must resolve, or verification rejects the answer built on it."""
    r = tools.pathways(DRUG, V)
    assert r.edge_ids
    known = {e["_id"] for e in tools.edges}
    assert all(i in known for i in r.edge_ids)


# --------------------------------------------------------------- routing
@pytest.mark.parametrize("question", [
    "which pathways does remdesivir's target sit in",
    "what pathway does this drug act through",
    "is nirmatrelvir in the same pathway as a dependency factor",
    "which reactome processes contain its targets",
])
def test_a_pathway_question_routes_to_the_pathway_intent(question):
    from kgav.agent import guards
    assert guards.classify_intent(question) is guards.QueryKind.PATHWAY


def test_it_is_matched_before_explain():
    """"explain which pathway this acts through" contains both triggers.

    EXPLAIN walks M1 and M2 and cannot answer a pathway question, so routing
    there returns a silence that reads as a negative.
    """
    from kgav.agent import guards
    assert guards.classify_intent(
        "explain which pathway this acts through") is guards.QueryKind.PATHWAY


def test_a_prospective_pathway_question_is_still_refused():
    """The refusal must not be routed around by adding the word pathway."""
    from kgav.agent import guards
    q = "what should we test against a pathway SARS-CoV-2 depends on"
    assert guards.classify_intent(q) is guards.QueryKind.PROSPECT
    assert not guards.check(guards.QueryKind.PROSPECT).allowed


def test_the_orchestrator_dispatches_a_pathway_question_to_the_tool(release):
    """Routing the intent nowhere would leave the capability unreachable."""
    from kgav.agent.orchestrator import Orchestrator
    orch = Orchestrator(GraphTools(release))
    brief = orch.build_brief(
        "which pathways do druga targets sit in for sars-cov-2")
    assert brief.intent == "pathway"
    assert brief.allowed
    # The TOOL ran, not just the classifier: its census line is unmistakable
    # and no other tool produces it.
    assert any("member cap" in f for f in brief.facts), brief.facts
    assert "PATHWAY ROUTE" in " ".join(brief.warnings)
    assert "TARGET DIRECTION" in " ".join(brief.warnings)


# ----------------------------------- direct-acting is not "no annotation"
def test_a_direct_acting_compound_is_not_called_unannotated(tools):
    """Remdesivir's shape, and the first thing the live run turned up.

    Its answer read "most small molecules here carry an activity
    measurement and no target annotation" -- true of the population, false
    about the compound, whose target is annotated, well-evidenced and
    viral. A reader is told the graph is missing something when it is not.
    """
    r = tools.pathways(DIRECT, V)
    assert not r.ok
    assert "no HUMAN target" in r.note
    assert "inhibit" in r.note and "viral protein(s)" in r.note
    assert "nsp5" in r.note
    assert "property of the compound, not a gap in the source" in r.note
    assert "no target annotation" not in r.note


def test_it_names_the_viral_targets_rather_than_counting_them_only(tools):
    r = tools.pathways(DIRECT, V)
    assert "UniProtKB:VIRAL" in r.edge_ids[0] or r.edge_ids
    known = {e["_id"] for e in tools.edges}
    assert all(i in known for i in r.edge_ids)


def test_it_points_at_the_question_that_does_have_an_answer(tools):
    """A refusal that does not say what to ask instead is a dead end."""
    assert "ask for its mechanism instead" in tools.pathways(DIRECT, V).note


def test_a_compound_with_no_target_at_all_still_says_so(tools):
    """The other branch must stay distinct, or the fix has just moved the
    conflation to the other side."""
    r = tools.pathways(OTHERDRUG, V)
    assert "no viral one either" in r.note
    assert "no target annotation at all" in r.note
    assert "direct-acting" not in r.note


def test_a_direct_acting_refusal_still_carries_the_caveats(tools):
    labels = [getattr(w, "label", "") for w in tools.pathways(DIRECT, V).warnings]
    assert "PATHWAY ROUTE" in labels


# ---------------------------------- a release that breaks its own schema
def test_a_non_viral_inhibits_edge_is_reported_not_dropped(tools):
    """INHIBITS carries object_constraints {is_viral: true}.

    Filtering those edges out silently would make a release that violates
    its own schema look like one that is merely thin -- the same conflation
    the direct-acting branch exists to end, one level down. The filter was
    there and nothing exercised it: removing it left every test passing.
    """
    r = tools.pathways(ROGUE, V)
    assert not r.ok
    assert "non-viral" in r.note and "object_constraints" in r.note
    assert "does not satisfy its own schema" in r.note
    assert "DEPF" in r.note
    assert "validate_release" in r.note


def test_a_rogue_edge_is_not_counted_as_a_viral_target(tools):
    """The claim the missing filter would have produced."""
    r = tools.pathways(ROGUE, V)
    # "non-viral protein(s)" contains "viral protein(s)", so the assertion
    # has to name the CLAIM rather than a substring of it.
    assert "it does inhibit" not in r.note
    assert "direct-acting" not in r.note


def test_a_rogue_edge_alongside_a_real_one_is_still_named(release):
    """The direct-acting answer stays, with the violation appended rather
    than swallowed by it."""
    import json as _json
    nodes = [_json.loads(x) for x in
             (release / "nodes.jsonl").read_text().splitlines()]
    edges = [_json.loads(x) for x in
             (release / "edges.jsonl").read_text().splitlines()]
    edges.append({"subject": DIRECT, "predicate": "INHIBITS",
                  "object": "UniProtKB:DEPF", "qualifiers": {},
                  "provenance": {"source": "test"}})
    (release / "edges.jsonl").write_text(
        "\n".join(_json.dumps(e) for e in edges))
    assert nodes
    r = GraphTools(release).pathways(DIRECT, V)
    assert "direct-acting" in r.note
    assert "NON-viral" in r.note and "forbid" in r.note
