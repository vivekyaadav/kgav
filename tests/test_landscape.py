"""The agent may not imply its ranking is predictive.

The measured position on this release is that the drugged human proteome is no
likelier to hit a virus's dependency factors than a degree-matched random
protein set of the same size. An answer that orders compounds, prints a
host-directed route count beside each, and never states that finding has made
the claim the landscape refuses -- and it has made it in the most persuasive
place, inside a fluent ranking.

These tests fix three things. That the finding is an ARTIFACT, loaded from the
run that measured it and refused when it describes a different release, a
subset of taxa, or too few permutations to quote. That an underpowered zero is
reported as UNTESTED rather than as evidence of absence. And that the
statement reaches the user: triage() attaches it unconditionally and verify()
fails an answer that drops it.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from kgav import landscape as L
from kgav.agent import guards
from kgav.agent.tools import GraphTools
from kgav.agent.verify import verify

V = "NCBITaxon:2697049"
OTHER = "NCBITaxon:64320"
DRUG = "INCHIKEY:NIRMA"


@pytest.fixture
def graph_dir(tmp_path) -> Path:
    """The smallest release triage() will rank over.

    Deliberately minimal and local: these tests are about what the ranking
    SAYS, not what it ranks, and a compound with no route at all still gets
    an order and still must carry the caveat.
    """
    nodes = [
        {"id": DRUG, "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "NIRMA",
                        "chembl_id": "CHEMBL1"}},
        {"id": "INCHIKEY:CYTOX", "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "CYTOX",
                        "chembl_id": "CHEMBL2"}},
        {"id": "UniProtKB:P0DTD1", "class": "Protein",
         "properties": {"taxon_id": V, "is_viral": True,
                        "protein_family": "nsp5", "sequence_hash": "h",
                        "reviewed": True}},
        {"id": V, "class": "OrganismTaxon", "properties": {"name": "SARS-CoV-2"}},
    ]
    edges = [
        {"subject": DRUG, "predicate": "INHIBITS", "object": "UniProtKB:P0DTD1",
         "qualifiers": {"value_nm": 3.0, "measure": "Ki"},
         "provenance": {"source": "chembl"}},
        {"subject": "UniProtKB:P0DTD1", "predicate": "ENCODED_BY_GENOME_OF",
         "object": V, "qualifiers": {}, "provenance": {"source": "spine"}},
    ]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in edges))
    return tmp_path


def _row(n_dep=200, n_obs=40, expected=44.0, p_naive=1e-9, p_matched=0.81):
    return {"n_dependency": n_dep, "n_observed": n_obs, "expected": expected,
            "p_naive": p_naive, "p_matched": p_matched}


def _doc(taxa=None, **over):
    doc = {"measured_on": "2026-10-25",
           "measured_by": "scripts/pathway_landscape.py",
           "release": "v0.1", "permutations": 2000, "seed": 0,
           "partial": False, "partial_reason": "",
           "taxa": taxa if taxa is not None else {V: _row()},
           "uncovered": []}
    doc.update(over)
    return doc


def _write(tmp_path: Path, doc) -> Path:
    p = tmp_path / "pathway_landscape.json"
    p.write_text(json.dumps(doc))
    return p


def _load(tmp_path, taxa=None, release="v0.1", **over) -> L.Landscape:
    return L.load(_write(tmp_path, _doc(taxa, **over)), release=release)


# ------------------------------------------------------------ the threshold
def test_the_power_threshold_is_the_poisson_derivation_it_claims():
    """MIN_EXPECTED_HITS is derived, not chosen for looking round.

    The docstring says it is the smallest matched expectation at which a total
    absence of overlap could reach p < 0.05. If that stops being arithmetic
    the constant has become a fitted parameter with a story attached, which is
    the objection calibration.py raises about hardcoded weights.
    """
    assert math.exp(-L.MIN_EXPECTED_HITS) < L.ALPHA
    assert math.exp(-(L.MIN_EXPECTED_HITS - 0.1)) > L.ALPHA


def test_the_permutation_floor_matches_the_resolution_it_claims():
    """MIN_PERMUTATIONS must actually resolve a p below ALPHA."""
    assert 1 / (L.MIN_PERMUTATIONS + 1) < L.ALPHA


# ----------------------------------------------------------------- refusals
def test_a_missing_file_is_refused_not_treated_as_a_null(tmp_path):
    land = L.load(tmp_path / "nope.json", release="v0.1")
    assert not land.usable
    assert land.for_virus(V) is None
    assert "no results file" in land.refused_reason


def test_a_landscape_from_another_release_is_refused(tmp_path):
    """The mis-keyed calibration bug, in a new place.

    A figure measured on one release stated about another is the most
    confident sentence this module can emit and the easiest to get wrong.
    """
    land = _load(tmp_path, release="v0.2-folds")
    assert not land.usable
    assert "v0.1" in land.refused_reason and "v0.2-folds" in land.refused_reason


def test_a_partial_run_is_refused(tmp_path):
    """A --virus run cannot support a statement about every taxon."""
    land = _load(tmp_path, partial=True,
                 partial_reason="run with --virus NCBITaxon:2697049")
    assert not land.usable
    assert "partial" in land.refused_reason.lower()


def test_too_few_permutations_is_refused(tmp_path):
    land = _load(tmp_path, permutations=50)
    assert not land.usable
    assert "50" in land.refused_reason


def test_no_permutation_count_is_refused_distinctly(tmp_path):
    land = _load(tmp_path, permutations=0)
    assert not land.usable
    assert "no permutation count" in land.refused_reason


def test_a_refused_landscape_yields_no_taxon_rows(tmp_path):
    """Refusal must not leave quotable numbers behind."""
    land = _load(tmp_path, release="other")
    assert land.tested() == [] and land.enriched() == []
    assert land.at_null() == [] and land.unreadable_ratio() == []
    assert land.depleted() == [] and not land.depletion_tested()


# -------------------------------------------------------------- the derivation
def test_the_ratio_is_derived_from_the_two_numbers_it_comes_from(tmp_path):
    land = _load(tmp_path, {V: _row(n_obs=40, expected=50.0)})
    assert land.for_virus(V).ratio == pytest.approx(0.8)


def test_a_stored_ratio_cannot_override_the_derived_one(tmp_path):
    """A third number that can disagree with its own inputs is not stored."""
    land = _load(tmp_path, {V: dict(_row(n_obs=40, expected=50.0), ratio=9.9)})
    assert land.for_virus(V).ratio == pytest.approx(0.8)


def test_an_empty_background_gives_an_infinite_ratio_not_a_crash(tmp_path):
    """JSON cannot write infinity, so the loader must produce it."""
    land = _load(tmp_path, {V: _row(n_obs=2, expected=0.0, p_matched=0.0)})
    assert land.for_virus(V).ratio == math.inf


def test_an_empty_background_with_no_hits_is_not_called_depleted(tmp_path):
    land = _load(tmp_path, {V: _row(n_obs=0, expected=0.0, p_matched=1.0)})
    t = land.for_virus(V)
    assert t.ratio == 0.0
    assert not t.ratio_readable
    assert not t.depleted


# ------------------------------------------- what the threshold actually gates
@pytest.mark.parametrize("expected,readable", [(2.9, False), (3.0, True),
                                               (0.4, False), (44.0, True)])
def test_the_threshold_gates_the_ratio_not_the_test(tmp_path, expected,
                                                    readable):
    """It governs whether the LOW end of the ratio means anything.

    It was introduced as a depletion-power gate, which was answering a
    question the run never asks: empirical_p measures the enrichment tail
    only. On the measured release that mismatch classified all eight taxa as
    "testable" and printed a 0.00 ratio with no qualification at all.
    """
    t = _load(tmp_path, {V: _row(expected=expected)}).for_virus(V)
    assert t.ratio_readable is readable
    # Tested either way: a small taxon can still return a significant excess.
    assert t.at_null is True


def test_a_small_taxon_is_still_tested_for_enrichment(tmp_path):
    """And is counted as tested, rather than dropped from the denominator."""
    land = _load(tmp_path, {V: _row(n_obs=0, expected=0.4, p_matched=1.0)})
    t = land.for_virus(V)
    assert t in land.tested() and t in land.at_null()
    assert t in land.unreadable_ratio()
    assert "not readable downward" in t.describe()


def test_a_small_taxon_can_still_clear_the_null_upward(tmp_path):
    """The reason the threshold must not gate testing.

    8 hits against 1.5 expected is a real excess, and a gate that called the
    taxon untested would have thrown the finding away. The expectation is
    deliberately BELOW the threshold: the point is that a taxon whose ratio
    cannot be read downward is still perfectly testable upward.
    """
    t = _load(tmp_path, {V: _row(n_dep=18, n_obs=8, expected=1.5,
                                 p_matched=0.004)}).for_virus(V)
    assert t.enriched
    assert not t.ratio_readable


def test_the_ratio_range_covers_every_tested_taxon(tmp_path):
    """Including the ones whose low end is soft. Trimming them would report
    a tidier range and hide that those taxa were measured at all."""
    land = _load(tmp_path, {V: _row(n_obs=40, expected=44.0),
                            OTHER: _row(n_obs=0, expected=0.4)})
    lo, hi = land.ratio_range()
    assert lo == 0.0 and hi == pytest.approx(40 / 44.0)


# ------------------------------------------------------- the untested direction
def test_a_sub_one_ratio_is_not_called_a_depletion(tmp_path):
    """The defect the measured release exposed.

    West Nile came back 0 of 18 against 3.1 expected -- ratio 0.00 -- and the
    summary printed it inside "ratios 0.00-0.94" beside "no enrichment". A
    reader takes that as depletion. The run never computed a lower tail.
    """
    land = _load(tmp_path, {V: _row(n_dep=18, n_obs=0, expected=3.1,
                                    p_matched=1.0)})
    t = land.for_virus(V)
    assert not t.depletion_tested and not t.depleted
    assert "not a measured depletion" in t.describe()
    assert "Depletion was not tested at all" in land.summary()


def test_a_lower_tail_in_the_file_is_used_when_present(tmp_path):
    """Forward compatible: the writer can start recording it without a
    change here, and the claim upgrades from untested to measured."""
    land = _load(tmp_path, {V: dict(_row(n_dep=18, n_obs=0, expected=3.1,
                                         p_matched=1.0),
                                    p_matched_less=0.004)})
    t = land.for_virus(V)
    assert t.depletion_tested and t.depleted
    assert "DEPLETED" in t.describe()
    assert "clears the null downward" in land.summary()


def test_a_missing_lower_tail_is_not_defaulted_to_one(tmp_path):
    """1.0 would read as "measured, and nowhere near depleted" -- a claim
    about a test that was never run."""
    t = _load(tmp_path, {V: _row()}).for_virus(V)
    assert t.p_matched_less is None


def test_a_partly_measured_direction_does_not_count_as_measured(tmp_path):
    """All-or-nothing: a mixed file cannot support "none is depleted"."""
    land = _load(tmp_path, {V: dict(_row(), p_matched_less=0.9),
                            OTHER: _row()})
    assert not land.depletion_tested()
    assert "Depletion was not tested at all" in land.summary()


# ----------------------------------------------------------------- summary
def test_the_summary_counts_every_taxon_it_tested(tmp_path):
    land = _load(tmp_path, {V: _row(n_obs=40, expected=44.0, p_matched=0.81),
                            OTHER: _row(n_obs=0, expected=0.4, p_matched=1.0)})
    s = land.summary()
    assert "2 taxa with dependency factors" in s
    # And the soft one is named, so the 0.00 at the bottom of the range can
    # be discounted rather than read as a result.
    assert OTHER in s and "ordinary draw" in s
    # The figures themselves, not just the shape of the sentence. Without
    # this the summary could print "ratios 0.00-0.00" over every taxon and
    # every other assertion here would still hold -- a test passing because
    # it does not exercise the case it is named for.
    assert f"{40 / 44.0:.2f}" in s
    assert "0.81" in s and V in s


def test_the_summary_prints_both_ends_of_the_range(tmp_path):
    """Both ends, from the taxa that were actually measured."""
    land = _load(tmp_path, {V: _row(n_obs=40, expected=44.0, p_matched=0.81),
                            OTHER: _row(n_obs=9, expected=10.0, p_matched=0.73),
                            })
    s = land.summary()
    assert f"{40 / 44.0:.2f}" in s and f"{9 / 10.0:.2f}" in s
    assert "2 taxa with dependency factors" in s
    # The taxon that comes CLOSEST to clearing the null, not the one furthest
    # from it. Only the smallest p supports the "none below 0.05" beside it;
    # quoting the largest would put a number next to a claim it cannot carry.
    assert "0.73" in s and OTHER in s
    assert "0.81" not in s


def test_the_summary_qualifies_a_range_built_from_soft_taxa(tmp_path):
    """Every taxon small: the enrichment claim still stands, and the range
    is explicitly not a finding in the low direction."""
    land = _load(tmp_path, {V: _row(n_obs=0, expected=0.4),
                            OTHER: _row(n_obs=1, expected=1.1)})
    s = land.summary()
    assert "not enriched" in s
    assert "The low end of that range is not a finding" in s
    assert V in s and OTHER in s


def test_the_summary_reports_enrichment_when_it_is_there(tmp_path):
    """The branch that says the route WORKS must be reachable.

    A module whose only output is "no enrichment" keeps saying it after the
    finding changes, which is the stale-figure failure this project keeps
    finding, pointing the other way.
    """
    land = _load(tmp_path, {V: _row(n_obs=90, expected=44.0, p_matched=0.0005)})
    s = land.summary()
    assert "IS enriched" in s
    assert "not enriched" not in s


def test_a_release_with_no_dependency_factors_says_so(tmp_path):
    land = _load(tmp_path, {})
    assert land.usable
    assert "no taxon" in land.summary() and "untested rather than null" in land.summary()


# ------------------------------------------------------------------ caveat
def test_the_caveat_with_no_landscape_does_not_assert_a_null():
    """Unmeasured and measured-null must not read the same."""
    text = guards.pathway_caveat(None, V)
    assert "has not been measured" in text
    assert "no enrichment" not in text
    assert "not enriched" not in text


def test_the_caveat_names_why_a_landscape_was_refused(tmp_path):
    land = _load(tmp_path, release="v0.2-folds")
    text = guards.pathway_caveat(land, V)
    assert "has not been measured" in text
    assert "v0.2-folds" in text


def test_the_caveat_quotes_this_taxons_own_numbers(tmp_path):
    land = _load(tmp_path, {V: _row(n_dep=200, n_obs=40, expected=44.0,
                                    p_matched=0.81),
                            OTHER: _row(n_dep=11, n_obs=9, expected=3.1,
                                        p_matched=0.004)})
    text = guards.pathway_caveat(land, V)
    assert "40 of 200" in text and "44.0" in text
    # The other taxon's enrichment may appear in the release-wide tail, but
    # this taxon's own line must be its own.
    assert text.split("Release-wide")[0].count("9 of 11") == 0


def test_the_caveat_refuses_to_borrow_another_taxons_numbers(tmp_path):
    land = _load(tmp_path, {OTHER: _row(n_obs=40, expected=44.0)})
    text = guards.pathway_caveat(land, V)
    assert "untested for it rather than null" in text
    assert "40 of" not in text.split("Release-wide")[0]


def test_the_caveat_renders_figures_from_the_artifact_not_the_source(tmp_path):
    """Change the file, change the text. A hardcoded figure would not move."""
    a = guards.pathway_caveat(
        _load(tmp_path, {V: _row(n_obs=40, expected=44.0)}), V)
    b = guards.pathway_caveat(
        _load(tmp_path, {V: _row(n_obs=7, expected=12.5)}), V)
    assert a != b
    assert "7 of" in b and "12.5" in b


def test_the_caveat_qualifies_a_ratio_that_cannot_be_read_downward(tmp_path):
    """Not "untested" -- the taxon WAS tested, in the enrichment direction.

    The branch this replaces claimed the taxon was untested on a depletion
    power basis, a direction the run never measures. On the real release
    that branch never fired at all, so a 0.00 ratio reached the reader bare.
    """
    land = _load(tmp_path, {V: _row(n_obs=0, expected=0.4, p_matched=1.0)})
    text = guards.pathway_caveat(land, V)
    own = text.split("Release-wide")[0]
    assert "no enrichment" in own
    assert "not a measured depletion" in own
    assert "not readable downward" in own and "0.4 expected" in own


def test_the_caveat_does_not_double_its_full_stops(tmp_path):
    """summary() ends in a full stop and the tail concatenated another."""
    land = _load(tmp_path, {V: _row(n_obs=0, expected=0.4, p_matched=1.0)})
    assert ".." not in guards.pathway_caveat(land, V)


def test_the_caveat_states_a_measured_depletion_when_there_is_one(tmp_path):
    land = _load(tmp_path, {V: dict(_row(n_dep=18, n_obs=0, expected=3.1,
                                         p_matched=1.0),
                                    p_matched_less=0.004)})
    own = guards.pathway_caveat(land, V).split("Release-wide")[0]
    assert "FEWER dependency factors than chance predicts" in own
    assert "0.0040" in own


def test_the_caveat_never_calls_an_untested_direction_a_depletion(tmp_path):
    """The claim the measured release would otherwise have invited."""
    land = _load(tmp_path, {V: _row(n_dep=18, n_obs=0, expected=3.1,
                                    p_matched=1.0)})
    own = guards.pathway_caveat(land, V).split("Release-wide")[0]
    assert "FEWER" not in own and "depleted" not in own.lower()
    assert "not a measured depletion" in own


def test_the_caveat_states_enrichment_plainly_when_measured(tmp_path):
    land = _load(tmp_path, {V: _row(n_obs=90, expected=44.0, p_matched=0.0005)})
    assert "IS enriched" in guards.pathway_caveat(land, V)


def test_the_caveat_is_a_registered_label():
    """An unregistered key raises rather than passing silently."""
    c = guards.caveat("pathway_route", "x")
    assert c.label == "PATHWAY ROUTE" and c.startswith("PATHWAY ROUTE:")


# ------------------------------------------------- the ranking must say it
def test_triage_states_the_pathway_null_with_a_landscape(graph_dir, tmp_path):
    land = _load(tmp_path, {V: _row(n_obs=40, expected=44.0, p_matched=0.81)})
    g = GraphTools(graph_dir, landscape=land)
    out = g.triage(["INCHIKEY:NIRMA", "INCHIKEY:CYTOX"], V)
    labels = [getattr(w, "label", "") for w in out.warnings]
    assert "PATHWAY ROUTE" in labels
    assert any("no enrichment" in str(w) for w in out.warnings)


def test_triage_states_it_even_with_no_landscape_at_all(graph_dir):
    """Silence is the failure A4 closes.

    An agent built without a landscape must say the route is unmeasured, not
    omit the subject. An unstated null and a measured one must not look the
    same from outside the answer.
    """
    out = GraphTools(graph_dir).triage(["INCHIKEY:NIRMA"], V)
    labels = [getattr(w, "label", "") for w in out.warnings]
    assert "PATHWAY ROUTE" in labels
    assert any("has not been measured" in str(w) for w in out.warnings)


def test_triage_does_not_quote_another_viruses_overlap(graph_dir, tmp_path):
    land = _load(tmp_path, {OTHER: _row(n_obs=40, expected=44.0)})
    out = GraphTools(graph_dir, landscape=land).triage(["INCHIKEY:NIRMA"], V)
    pathway = next(str(w) for w in out.warnings
                   if getattr(w, "label", "") == "PATHWAY ROUTE")
    assert "untested for it rather than null" in pathway


def test_the_ranking_caveat_and_the_pathway_caveat_are_both_present(
        graph_dir, tmp_path):
    """One is about the sort key, the other about the route. Not substitutes.

    A reader can accept "host-directed routes contribute nothing to this
    order" and still read the host-directed count printed beside each compound
    as weak supporting evidence. Only the pathway caveat refuses that.
    """
    land = _load(tmp_path, {V: _row()})
    out = GraphTools(graph_dir, landscape=land).triage(["INCHIKEY:NIRMA"], V)
    labels = [getattr(w, "label", "") for w in out.warnings]
    assert "RANKING" in labels and "PATHWAY ROUTE" in labels


# ----------------------------------------- a retrieved host-directed route
HOST_PATH = [{"id": "INCHIKEY:X", "class": "SmallMolecule", "properties": {}},
             {"id": "UniProtKB:Q99720", "class": "Protein",
              "properties": {"gene_symbol": "SIGMAR1", "is_viral": False}},
             {"id": V, "class": "OrganismTaxon",
              "properties": {"label": "SARS-CoV-2",
                             "family": "Coronaviridae"}}]


def test_a_host_directed_path_carries_the_pathway_caveat(tmp_path):
    """explain() shows routes, and a route shown is the pathway hypothesis.

    This path was added with the ranking and went untested: a mutation
    deleting the whole branch left every other test passing.
    """
    land = _load(tmp_path, {V: _row(n_obs=40, expected=44.0, p_matched=0.81)})
    out = guards.path_warnings(HOST_PATH, "M4 pathway-mediated",
                               landscape=land, virus=V)
    labels = [getattr(w, "label", "") for w in out]
    assert "PATHWAY ROUTE" in labels
    assert any("40 of 200" in str(w) for w in out)


def test_a_direct_acting_path_does_not_carry_it(tmp_path):
    """M1 is not the pathway route, so the caveat would be noise there."""
    land = _load(tmp_path, {V: _row()})
    out = guards.path_warnings(HOST_PATH, "M1 direct-acting",
                               landscape=land, virus=V)
    assert "PATHWAY ROUTE" not in [getattr(w, "label", "") for w in out]


def test_a_host_directed_path_without_a_landscape_adds_nothing(tmp_path):
    """Supplied-only here, unlike the ranking, and the asymmetry is tested.

    A retrieved path already carries HOST-DIRECTED ROUTE with its measured
    channel performance, and showing one orders nothing. triage() is where
    silence would be the failure, and it attaches the caveat unconditionally.
    """
    out = guards.path_warnings(HOST_PATH, "M4 pathway-mediated")
    labels = [getattr(w, "label", "") for w in out]
    assert "PATHWAY ROUTE" not in labels
    assert "HOST-DIRECTED ROUTE" in labels


# ------------------------------------------------ and it must reach the user
def test_an_answer_that_drops_the_pathway_caveat_fails_verification():
    """The test A4 asked for, rather than a caveat printed under the answer.

    verify() enforces warning survival by LABEL, so an answer that renders a
    fluent ranking and omits this one does not pass.
    """
    warning = guards.caveat("pathway_route", "no enrichment in any taxon.")
    answer = "Nirmatrelvir ranks first on direct-acting potency [E1]."
    fact = "nirmatrelvir ranks first on direct-acting potency"
    dropped = verify(answer, facts=[fact], citable_edge_ids=["E1"],
                     required_warnings=[warning])
    assert not dropped.passed
    assert any(f.code == "warning_dropped" and "pathway_route" in f.detail
               for f in dropped.failures)

    # And it passes once the label is carried, so the test above is failing on
    # the dropped caveat rather than on something else in the answer.
    kept = verify(answer + " PATHWAY ROUTE: no enrichment in any taxon.",
                  facts=[fact, str(warning)],
                  citable_edge_ids=["E1"], required_warnings=[warning])
    assert not any(f.code == "warning_dropped" for f in kept.failures)
