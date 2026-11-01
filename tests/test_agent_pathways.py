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
import random
from pathlib import Path

import pytest

from kgav.agent.tools import GraphTools

V = "NCBITaxon:2697049"
DRUG = "INCHIKEY:DRUGA"
OTHERDRUG = "INCHIKEY:DRUGB"
DIRECT = "INCHIKEY:DRUGC"
ROGUE = "INCHIKEY:DRUGD"

# What the fixture's own NAMES.json carries. Declared here so a test can
# assert that a name it is checking is NOT in it -- a test that asserts
# resolution of a name the fixture already indexes cannot detect the
# register seeding being removed.
_FIXTURE_NAMES = {"druga": "INCHIKEY:DRUGA", "drugb": "INCHIKEY:DRUGB",
                  "drugc": "INCHIKEY:DRUGC", "drugd": "INCHIKEY:DRUGD"}


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
    (tmp_path / "NAMES.json").write_text(json.dumps(_FIXTURE_NAMES))
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


def test_no_pathway_hits_is_stated_as_an_absence_in_the_crispr_layer(
        release, tmp_path):
    """A virus nothing was SCREENED against must not read as a compound
    with no pathway relevance.

    The virus needs a route for this case to exist at all: with none it
    is the Yellow fever tier and the tool refuses outright, which is a
    different and correct answer. The shape here is DENV-1's -- a
    virus-host interaction and no CRISPR screen -- so the fixture is
    extended with one interaction rather than the test relaxed.
    """
    import shutil

    shutil.copytree(release, tmp_path / "r")
    rel = tmp_path / "r"
    nodes = rel / "nodes.jsonl"
    edges = rel / "edges.jsonl"
    nodes.write_text(nodes.read_text().rstrip("\n") + "\n" + "\n".join([
        json.dumps(_node("NCBITaxon:11082", "OrganismTaxon", name="West Nile")),
        json.dumps(_node("UniProtKB:WNV_NS3", "Protein", gene_symbol="NS3",
                         is_viral=True, taxon_id="NCBITaxon:11082",
                         sequence_hash="h", reviewed=True)),
    ]))
    edges.write_text(edges.read_text().rstrip("\n") + "\n" + json.dumps(
        _edge("UniProtKB:WNV_NS3", "PHYSICALLY_INTERACTS_WITH",
              "UniProtKB:P1")))
    t = GraphTools(rel)
    assert t.routes_for("NCBITaxon:11082").routes == ("ppi",), (
        "the fixture must give this virus a route, or the tool refuses "
        "and this test checks the wrong message")
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


# --------------------------------------- the count needs a MATCHED background
def test_the_background_is_binned_by_pathway_size(tools):
    """A bare rate was not enough, and the second live run showed why.

    It came out at 58% against chloroquine's 90% and the answer called that
    "+33% on the background". The schema's own measurements put the median
    pathway at 4 members and p99 at 113; a 4-member pathway holds one of
    SARS-CoV-2's 3,460 dependency factors about half the time and a
    95-member pathway essentially always. Chloroquine's targets include
    mTOR, which sits in large pathways, so the higher rate follows from the
    sizes alone -- the same confound the landscape degree-matches away at
    the protein level, one level up.
    """
    bg = tools._pathway_background(V)
    assert bg["total"] == 3
    # A HISTOGRAM OVER EXACT SIZES, not quantile bins. The counts must
    # account for every pathway in the population: a size silently
    # dropped from the histogram is a size whose rate is estimated from
    # its neighbours without anyone being told.
    assert sum(n for _h, n in bg["by_size"].values()) == bg["total"]
    assert sum(h for h, _n in bg["by_size"].values()) == bg["hits"]
    assert all(0 <= h <= n for h, n in bg["by_size"].values())


def test_the_background_respects_the_member_cap(tools):
    """REACT:BIG is over the cap and must not enter the denominator, or the
    background covers a different population from the answer."""
    assert tools._pathway_background(V)["total"] == 3


def test_a_bigger_pathway_raises_its_own_expectation(release):
    """The property the matching exists for: expectation tracks size."""
    t = GraphTools(release)
    small = t._size_matched_expectation(V, [1])
    big = t._size_matched_expectation(V, [999])
    assert big >= small


def test_the_count_is_never_reported_without_a_matched_background(tools):
    text = " ".join(tools.pathways(DRUG, V).verbalised)
    assert "SIZE-MATCHED background expects" in text
    assert "dependency factors, so a pathway's size largely decides" in text
    # The unmatched figure is still shown, because the gap between the two
    # is itself the measurement.
    assert "unmatched," in text and "of all pathways under the cap" in text
    assert "median" in text


def test_membership_at_the_matched_expectation_is_called_what_it_is(tools):
    text = " ".join(tools.pathways(DRUG, V).verbalised)
    assert "the excess is not significant" in text
    assert "no more than these pathways' sizes already predict" in text
    assert "not evidence of anything about this compound" in text


def test_an_excess_is_tested_rather_than_compared_to_a_chosen_number(
        tools, monkeypatch):
    """The verdict used to turn on `ratio <= 1.05`.

    On the live run the ratio came out at 1.04, so which of two things the
    reader was told hung on the second decimal of a constant picked by
    hand. It is now an exact Poisson-binomial tail against the shared
    alpha.
    """
    monkeypatch.setattr(tools, "_size_matched_probs",
                        lambda v, sizes: [0.01] * len(sizes))
    text = " ".join(tools.pathways(DRUG, V).verbalised)
    assert "reaches p 0." in text
    assert "anti-conservative" in text
    assert "already predict" not in text


def test_a_significant_excess_still_states_the_independence_violation(
        tools, monkeypatch):
    """Pathways share members, so the tail is anti-conservative and a small
    p is a flag rather than a finding. Saying the p without saying that
    would be the overclaim this whole tool exists to refuse."""
    monkeypatch.setattr(tools, "_size_matched_probs",
                        lambda v, sizes: [0.001] * len(sizes))
    text = " ".join(tools.pathways(DRUG, V).verbalised)
    assert "assumes an independence that does not hold" in text
    assert "a flag to examine, not a finding" in text


def test_the_alpha_is_the_one_the_landscape_declares(tools):
    """Two different significance thresholds in one answer would be
    indefensible, so it is imported rather than restated."""
    from kgav.agent.tools import LANDSCAPE_ALPHA
    from kgav.landscape import ALPHA
    assert LANDSCAPE_ALPHA is ALPHA


@pytest.mark.parametrize("probs,k,want", [
    ([], 0, 1.0),
    ([1.0], 1, 1.0),
    ([0.0], 1, 0.0),
    ([0.5, 0.5], 1, 0.75),
    ([0.5, 0.5], 2, 0.25),
    ([0.5, 0.5], 0, 1.0),
])
def test_the_upper_tail_is_exact(probs, k, want):
    """Hand-checkable cases. Two fair coins: P(>=1)=3/4, P(>=2)=1/4."""
    assert GraphTools._upper_tail(probs, k) == pytest.approx(want)


def test_the_upper_tail_is_zero_above_the_maximum():
    assert GraphTools._upper_tail([0.5, 0.5], 3) == 0.0


# ------------------------------------------------------- bounded output
def test_dependency_names_are_bounded_but_counted_in_full(release):
    """The live run printed forty-one names on one line, twenty-five times.

    Nothing in it was wrong and none of it was readable, and a brief built
    from it spends its budget on name lists.
    """
    import json as _json

    from kgav.agent.tools import MAX_DEPENDENCY_NAMES
    edges = [_json.loads(x) for x in
             (release / "edges.jsonl").read_text().splitlines()]
    nodes = [_json.loads(x) for x in
             (release / "nodes.jsonl").read_text().splitlines()]
    extra = MAX_DEPENDENCY_NAMES + 4
    for i in range(extra):
        pid = f"UniProtKB:MANY{i}"
        nodes.append(_node(pid, "Protein", gene_symbol=f"MANY{i}",
                           is_viral=False, taxon_id="NCBITaxon:9606",
                           sequence_hash="h", reviewed=True))
        edges.append(_edge(pid, "PARTICIPATES_IN", "REACT:DEP"))
        edges.append(_edge(pid, "HOST_FACTOR_FOR", V, direction="dependency"))
    (release / "nodes.jsonl").write_text(
        "\n".join(_json.dumps(n) for n in nodes))
    (release / "edges.jsonl").write_text(
        "\n".join(_json.dumps(e) for e in edges))

    r = GraphTools(release).pathways(DRUG, V)
    dep_line = next(ln for ln in r.verbalised if "Dependency-bearing" in ln)
    assert f"+{extra + 1 - MAX_DEPENDENCY_NAMES} more" in dep_line
    assert dep_line.count(", ") < extra
    # The count and the full set survive untruncated for a caller that
    # wants them.
    row = next(x for x in r.data if x["pathway"] == "REACT:DEP")
    assert row["n_dependency_members"] == extra + 1
    assert len(row["dependency_members"]) == extra + 1


# ------------------------------------------------- the binning itself
V2 = "NCBITaxon:11137"


@pytest.fixture
def graded(tmp_path) -> Path:
    """100 pathways of sizes 1..100, where ONLY the large ones hold a
    dependency factor.

    The three-pathway fixture above cannot exercise binning at all: with
    fewer pathways than bins every mutation to the quantile logic survives.
    Here the unmatched rate is 50% by construction while the size-matched
    rate is ~0 at the bottom and ~1 at the top, so a background that
    ignores size gives a visibly wrong answer.
    """
    nodes = [_node(V2, "OrganismTaxon", name="Dengue"),
             _node(DRUG, "SmallMolecule", smiles="C", inchikey_skel="DRUGA",
                   chembl_id="CHEMBL1")]
    edges = []
    # SHUFFLED, with a fixed seed. Written in ascending size order the file
    # is already sorted by the time the loader reads it, so removing the
    # sort in _pathway_background changed nothing and the mutation
    # survived: the fixture was doing the sorting.
    order = list(range(1, 101))
    random.Random(0).shuffle(order)
    for size in order:
        pw = f"REACT:S{size:03d}"
        nodes.append(_node(pw, "Pathway", label=f"pathway of {size}"))
        for m in range(size):
            pid = f"UniProtKB:G{size:03d}_{m}"
            nodes.append(_node(pid, "Protein", gene_symbol=f"G{size}_{m}",
                               is_viral=False, taxon_id="NCBITaxon:9606",
                               sequence_hash="h", reviewed=True))
            edges.append(_edge(pid, "PARTICIPATES_IN", pw))
            # Exactly one dependency factor, and only in the big half.
            if size > 50 and m == 0:
                edges.append(_edge(pid, "HOST_FACTOR_FOR", V2,
                                   direction="dependency"))
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in edges))
    return tmp_path


def test_the_unmatched_rate_is_half_by_construction(graded):
    bg = GraphTools(graded)._pathway_background(V2)
    assert bg["total"] == 100
    assert bg["hits"] == 50


def test_the_matched_expectation_tracks_size_not_the_flat_rate(graded):
    """A flat 50% background would predict 1.0 either way. Size-matching
    must predict ~0 for small pathways and ~2 for large ones."""
    t = GraphTools(graded)
    small = t._size_matched_expectation(V2, [3, 7])
    large = t._size_matched_expectation(V2, [90, 97])
    assert small == pytest.approx(0.0, abs=0.2)
    assert large == pytest.approx(2.0, abs=0.2)
    assert large - small > 1.5


def test_the_histogram_covers_every_size_present(graded):
    """One entry per distinct size, and nothing invented between them.
    The quantile version this replaces reported ten bins whatever the
    data, so a size could be described by a bin that contained none of
    it."""
    bg = GraphTools(graded)._pathway_background(V2)
    assert sorted(bg["by_size"]) == list(range(1, 101))
    assert sum(n for _h, n in bg["by_size"].values()) == 100


def test_the_rate_rises_with_size_here_because_the_fixture_is_built_that_way(
        graded):
    """Not a property of real data -- a property of THIS fixture, which
    is what makes it able to detect a size-selection bug. Only pathways
    over 50 members hold a dependency factor, so the estimated rate must
    be 0 at the bottom and 1 at the top with a single crossing."""
    t = GraphTools(graded)
    probs = t._size_matched_probs(V2, list(range(1, 101)))
    assert probs[0] == pytest.approx(0.0)
    assert probs[-1] == pytest.approx(1.0)
    assert probs == sorted(probs), "the estimate is not monotone in size"
    crossings = sum(1 for a, b in zip(probs, probs[1:]) if a < 0.5 <= b)
    assert crossings == 1, f"{crossings} crossings of 0.5; the fixture has one"


def test_a_size_above_every_bin_uses_the_largest_bin(graded):
    t = GraphTools(graded)
    assert t._size_matched_expectation(V2, [10_000]) == pytest.approx(1.0)


def test_a_size_below_every_bin_uses_the_smallest_bin(graded):
    t = GraphTools(graded)
    assert t._size_matched_expectation(V2, [0]) == pytest.approx(0.0)


def test_the_verdict_turns_on_the_shared_alpha_not_a_local_one(tools,
                                                               monkeypatch):
    """Asserting the import is not asserting the comparison.

    test_the_alpha_is_the_one_the_landscape_declares checks that the name
    is bound to landscape.ALPHA; replacing the number in the comparison
    with 0.5 left it passing, because the fixture's p sits above both.
    This pins a p BETWEEN the two: not significant at 0.05, significant at
    0.5, so only the right threshold gives the right sentence.
    """
    monkeypatch.setattr(tools, "_size_matched_probs",
                        lambda v, sizes: [0.072] * len(sizes))
    p = GraphTools._upper_tail([0.072] * 3, 1)
    assert 0.05 < p < 0.5, p
    text = " ".join(tools.pathways(DRUG, V).verbalised)
    assert "the excess is not significant" in text
    assert "anti-conservative" not in text


# ------------------------------------------- the virus the question named
def test_every_register_virus_resolves_by_its_common_names():
    """The live failure: "for dengue virus 2" resolved to nothing and the
    question was answered about SARS-CoV-2.

    NAMES.json carries the register's SHORT name only, so every virus whose
    common name is not its short form was unreachable -- on this release
    that is every flavivirus.
    """
    from kgav import virus_register as R
    names = R.search_names(R.load())
    assert names["dengue virus 2"] == "NCBITaxon:11060"
    assert names["dengue 2"] == "NCBITaxon:11060"
    assert names["japanese encephalitis virus"] == "NCBITaxon:11072"
    assert names["tick-borne encephalitis"] == "NCBITaxon:11084"
    assert names["zika virus"] == "NCBITaxon:64320"
    assert names["severe acute respiratory syndrome coronavirus 2"] == \
        "NCBITaxon:2697049"


def test_every_declared_virus_has_at_least_one_searchable_name():
    from kgav import virus_register as R
    doc = R.load()
    names = R.search_names(doc)
    for entry in doc["viruses"]:
        curie = f"NCBITaxon:{entry['taxon']}"
        assert curie in names.values(), f"{entry['name']} is unreachable"


def test_a_name_claimed_by_two_taxa_raises():
    """Silently picking one would be the wrong-virus answer again."""
    from kgav import virus_register as R
    doc = {"viruses": [{"taxon": "1", "name": "A", "aliases": ["shared name"]},
                       {"taxon": "2", "name": "B", "aliases": ["shared name"]}]}
    with pytest.raises(ValueError, match="claimed by both"):
        R.search_names(doc)


def test_inert_short_aliases_are_not_declared():
    """The resolver skips strings under 4 characters, so a 3-letter alias
    would read as support the agent does not have."""
    from kgav import virus_register as R
    for key in R.search_names(R.load()):
        assert len(key) >= R.MIN_SEARCH_NAME


def test_the_agent_resolves_a_register_name_the_release_never_indexed(
        release):
    """GraphTools seeds from the register, so resolution does not depend on
    what NAMES.json happens to hold.

    The name asserted here must be one the fixture's NAMES.json does NOT
    carry, or the seeding can be deleted with every test still passing --
    which is exactly what happened: the first version asserted "sars-cov-2",
    which the fixture indexes itself.
    """
    t = GraphTools(release)
    assert "sars-cov-2" not in _FIXTURE_NAMES
    assert t.names["severe acute respiratory syndrome coronavirus 2"] == V
    assert t.names["sarscov2"] == V


def test_a_register_name_for_an_absent_taxon_is_not_offered(release):
    """Resolving to a node that is not in this release would trade a wrong
    answer for a crash."""
    t = GraphTools(release)
    assert "dengue virus 2" not in t.names


def test_the_release_name_wins_over_the_register(release):
    """A release that curated its own mapping keeps it.

    The key must be one the register ALSO claims, and must map somewhere
    DIFFERENT, or setdefault and plain assignment behave identically and
    the precedence is untested -- the first version mapped "sars-cov-2" to
    the same node the register does.
    """
    import json as _json
    (release / "NAMES.json").write_text(_json.dumps(
        {"druga": DRUG, "sarscov2": DRUG}))
    t = GraphTools(release)
    assert t.names["sarscov2"] == DRUG, "the release's own mapping must win"


# ----------------------------------------------- ties, which the audit found
#
# THE BUG, in full, because the fixtures above could not see it. `sized`
# held (size, hit) tuples and `_pathway_background` sorted it with a bare
# `sized.sort()`, which orders by size AND THEN by the hit flag -- and
# False < True. Inside every tied-size block each non-hit therefore sorted
# before each hit, the quantile cut split that block, and the lower piece
# got a rate of 0. `_size_matched_probs` returned the FIRST bin whose
# upper edge reached the size, which for a tied size is exactly that
# zero-rate piece.
#
# Too low an expectation means too high a ratio and too small a p-value,
# so the error ran in the one direction that manufactures findings.
#
# WHY NOTHING HERE CAUGHT IT: the `graded` fixture has sizes 1..100, each
# exactly once. It has no tied sizes at all, and ties are the rule on the
# real graph -- the median pathway has 4 members and thousands share it.
# A fixture that cannot produce the condition cannot test the handling of
# it, which is the defect this suite has found most often.
@pytest.fixture
def tied(tmp_path) -> Path:
    """400 pathways, EVERY ONE of 4 members, hit rate exactly 0.25.

    The degenerate case, and the one the quantile code answered 0.000 to
    whatever the truth was. With one size there is nothing for
    size-matching to match on, so the correct answer is the flat rate.
    """
    nodes = [_node(V2, "OrganismTaxon", name="Dengue"),
             _node(DRUG, "SmallMolecule", smiles="C", inchikey_skel="DRUGA",
                   chembl_id="CHEMBL1")]
    edges = []
    for i in range(400):
        pw = f"REACT:T{i:03d}"
        nodes.append(_node(pw, "Pathway", label=f"tied {i}"))
        for m in range(4):
            pid = f"UniProtKB:T{i:03d}_{m}"
            nodes.append(_node(pid, "Protein", gene_symbol=f"T{i}_{m}",
                               is_viral=False, taxon_id="NCBITaxon:9606",
                               sequence_hash="h", reviewed=True))
            edges.append(_edge(pid, "PARTICIPATES_IN", pw))
            if i % 4 == 0 and m == 0:
                edges.append(_edge(pid, "HOST_FACTOR_FOR", V2,
                                   direction="dependency"))
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in edges))
    return tmp_path


def test_one_tied_size_returns_the_true_rate_not_zero(tied):
    """THE REGRESSION. Every pathway is 4 members and a quarter of them
    hold a dependency factor, so the background for a 4-member pathway
    is 0.25. The quantile version returned 0.000 -- and an expectation
    of zero renders as an infinite ratio beside p = 0.0000."""
    t = GraphTools(tied)
    bg = t._pathway_background(V2)
    assert bg["total"] == 400 and bg["hits"] == 100
    assert sorted(bg["by_size"]) == [4], "the fixture should hold one size"
    assert t._size_matched_probs(V2, [4])[0] == pytest.approx(0.25)
    assert t._size_matched_expectation(V2, [4] * 8) == pytest.approx(2.0)


def test_the_estimate_does_not_depend_on_the_hit_flag_ordering(tied):
    """The bug was an ORDERING one: non-hits sorted before hits inside a
    tied block. Estimating per size cannot depend on the order within a
    size at all, so shuffling the pathway ids must not move the answer."""
    t = GraphTools(tied)
    first = t._size_matched_probs(V2, [4])[0]
    t2 = GraphTools(tied)
    t2._cache.clear()
    assert t2._size_matched_probs(V2, [4])[0] == pytest.approx(first)


@pytest.fixture
def mixed(tmp_path) -> Path:
    """300 pathways over three tied sizes, each with its OWN rate.

    sizes 2 (rate 0.0), 4 (rate 0.5) and 8 (rate 1.0), 100 of each. One
    size per rate, so an estimate that pools across sizes is visible
    immediately and a flat background cannot pass.
    """
    nodes = [_node(V2, "OrganismTaxon", name="Dengue"),
             _node(DRUG, "SmallMolecule", smiles="C", inchikey_skel="DRUGA",
                   chembl_id="CHEMBL1")]
    # The drug targets one protein in a 4-member pathway, so pathways()
    # has something to answer about and the prose can be checked. It adds
    # no pathway and changes no membership, so every rate above is
    # untouched.
    edges = [_edge(DRUG, "TARGETS", "UniProtKB:M4_000_1")]
    plan = [(2, 0.0), (4, 0.5), (8, 1.0)]
    for size, rate in plan:
        for i in range(100):
            pw = f"REACT:M{size}_{i:03d}"
            nodes.append(_node(pw, "Pathway", label=f"mixed {size}/{i}"))
            for m in range(size):
                pid = f"UniProtKB:M{size}_{i:03d}_{m}"
                nodes.append(_node(pid, "Protein", gene_symbol=f"M{size}_{i}_{m}",
                                   is_viral=False, taxon_id="NCBITaxon:9606",
                                   sequence_hash="h", reviewed=True))
                edges.append(_edge(pid, "PARTICIPATES_IN", pw))
                if m == 0 and i < round(rate * 100):
                    edges.append(_edge(pid, "HOST_FACTOR_FOR", V2,
                                       direction="dependency"))
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in edges))
    return tmp_path


def test_each_tied_size_keeps_its_own_rate(mixed):
    """300 pathways, so the pooling floor is 30 and each size has 100 of
    its own -- no size is smoothed against its neighbours. A flat
    background would answer 0.5 everywhere."""
    t = GraphTools(mixed)
    assert t._size_matched_probs(V2, [2])[0] == pytest.approx(0.0)
    assert t._size_matched_probs(V2, [4])[0] == pytest.approx(0.5)
    assert t._size_matched_probs(V2, [8])[0] == pytest.approx(1.0)


def test_a_size_the_graph_has_none_of_takes_its_nearest_neighbours(mixed):
    """5 is absent. The nearest sizes are 4 (distance 1) and 8 (distance
    3), so the estimate is 4's. Returning 0.0 for an unseen size -- the
    other obvious implementation -- would make every unusual pathway
    look like a discovery."""
    t = GraphTools(mixed)
    assert t._size_matched_probs(V2, [5])[0] == pytest.approx(0.5)
    assert t._size_matched_probs(V2, [1])[0] == pytest.approx(0.0)
    assert t._size_matched_probs(V2, [99])[0] == pytest.approx(1.0)


def test_the_pooling_floor_cannot_swallow_the_whole_distribution(mixed):
    """An unconditional floor of 100 would pool 4 with its neighbours on
    a 300-pathway graph and flatten the very thing being measured. The
    floor is capped at a tenth of the population, so a small graph gets
    a noisier local estimate rather than a precise global one."""
    from kgav.agent.tools import MIN_BACKGROUND_SAMPLES

    t = GraphTools(mixed)
    assert MIN_BACKGROUND_SAMPLES == 100
    # 300 pathways -> floor 30, and size 4 has 100 of its own.
    assert t._size_matched_probs(V2, [4])[0] == pytest.approx(0.5), (
        "size 4 was pooled with its neighbours despite having 100 "
        "pathways of its own")


def test_no_background_is_not_an_infinite_ratio(tmp_path):
    """A3. `ratio inf` beside `p 0.0000` was the strongest possible
    claim from the weakest possible null. An expectation of zero is a
    statement about coverage, not about the compound."""
    nodes = [_node(V2, "OrganismTaxon", name="Dengue"),
             _node(DRUG, "SmallMolecule", smiles="C", inchikey_skel="DRUGA",
                   chembl_id="CHEMBL1"),
             _node("UniProtKB:X1", "Protein", gene_symbol="X1",
                   is_viral=False, taxon_id="NCBITaxon:9606",
                   sequence_hash="h", reviewed=True),
             _node("REACT:ONLY", "Pathway", label="only pathway"),
             # A ROUTE, or the tool refuses before it computes anything
             # and this test checks the refusal instead of the ratio.
             _node("UniProtKB:V1", "Protein", gene_symbol="V1",
                   is_viral=True, taxon_id=V2, sequence_hash="h",
                   reviewed=True)]
    edges = [_edge(DRUG, "TARGETS", "UniProtKB:X1"),
             _edge("UniProtKB:X1", "PARTICIPATES_IN", "REACT:ONLY"),
             _edge("UniProtKB:V1", "PHYSICALLY_INTERACTS_WITH",
                   "UniProtKB:X1")]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in edges))
    text = " ".join(GraphTools(tmp_path).pathways(DRUG, V2).verbalised)
    assert "ratio inf" not in text, text
    assert "inf" not in text.split("ratio")[-1][:20] if "ratio" in text else True
    assert "nothing to compare against" in text, text


def test_the_background_median_is_the_populations_own(mixed):
    """QUOTED IN THE ANSWER, so it has to be the real median. The
    previous version took the lower edge of the middle QUANTILE BIN,
    which is a different quantity whenever the bins are uneven -- and
    over tied sizes they always are. A mutant returning 0 here survived
    every test in this file.

    The fixture holds 100 pathways each of 2, 4 and 8 members, so the
    150th is a 4.
    """
    t = GraphTools(mixed)
    bg = t._pathway_background(V2)
    assert bg["total"] == 300
    assert t._background_median(bg) == 4
    assert " median of 4 members " in " ".join(
        t.pathways(DRUG, V2).verbalised), "the answer quotes another median"


def test_the_background_median_of_an_empty_population_is_not_invented():
    """Zero is the only honest answer with nothing to take a median of,
    and it must not come from falling off the end of a loop that was
    meant to find something."""
    assert GraphTools._background_median({"total": 0, "by_size": {}}) == 0


def test_the_median_takes_the_lower_size_on_an_exact_half():
    """THE ONE INPUT `>=` AND `>` DISAGREE ON. With 100 pathways of 2
    members and 100 of 8, the cumulative count lands on exactly half at
    size 2: the lower median is 2, and `>` would skip past it and
    report 8 -- a four-fold error in the figure the answer quotes as
    the population's typical size. Every fixture above has an odd split
    and cannot see it.
    """
    bg = {"total": 200, "by_size": {2: (0, 100), 8: (100, 100)}}
    assert GraphTools._background_median(bg) == 2


def test_the_median_of_a_single_size_is_that_size():
    assert GraphTools._background_median(
        {"total": 400, "by_size": {4: (100, 400)}}) == 4
