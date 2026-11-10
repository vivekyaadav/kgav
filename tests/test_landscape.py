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

import itertools
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
# BOTH TESTS BELOW REPLACE ONES THAT PINNED THE UNCORRECTED ALPHA. They
# were right about the arithmetic and wrong about which alpha it answers
# to, which is why they passed for months:
#
#   assert math.exp(-L.MIN_EXPECTED_HITS) < L.ALPHA
#   assert 1 / (L.MIN_PERMUTATIONS + 1) < L.ALPHA
#
# L.ALPHA is 0.05. Every verdict in that module is judged against
# ALPHA / n_tests, which on v0.1 is 0.05/44 = 0.001136. So the two
# constants were checked against a threshold the module does not use,
# and the checks held while the constants were wrong by a factor of two
# (3.0 against 6.78) and by a factor of forty-four respectively. Stated
# here rather than quietly rewritten, because a test that was passing is
# the thing a reader will assume was right.
def test_the_readability_threshold_is_derived_from_THIS_alpha():
    """The smallest matched expectation at which a total absence of
    overlap could reach the CORRECTED alpha. exp(-m) < alpha."""
    for alpha in (0.05, 0.05 / 22, 0.05 / 44, 0.001):
        m = L.min_expected_hits(alpha)
        assert math.exp(-m) <= alpha * (1 + 1e-9)
        assert math.exp(-(m - 0.01)) > alpha


def test_the_readability_threshold_rises_as_the_family_grows():
    """The direction is the point: correcting for more tests makes a
    zero overlap LESS surprising, so more is needed before the low end
    of a ratio means anything."""
    assert L.min_expected_hits(0.05 / 44) > L.min_expected_hits(0.05)
    assert abs(L.min_expected_hits(0.05) - 3.0) < 0.01
    assert abs(L.min_expected_hits(0.05 / 44) - 6.78) < 0.01


def _overlap_at(expected, n_observed, alpha):
    """One TaxonOverlap, built directly so the alpha can be set.

    NAMED _overlap_at BECAUSE _row ALREADY EXISTS in this file and my
    first version shadowed it, breaking 49 tests at once. A module-level
    helper name is a shared namespace.
    """
    return L.TaxonOverlap(virus="NCBITaxon:11060", n_dependency=50,
                          n_observed=n_observed, expected=expected,
                          p_naive=0.5, p_matched=0.5, p_matched_less=0.045,
                          alpha=alpha, n_tests=1)


def test_a_row_reads_its_own_alpha_not_the_module_constant():
    """TaxonOverlap carries `alpha` per row with a comment saying it is
    held there "so a verdict can never be computed against a different
    threshold from the one reported next to it" -- and ratio_readable
    reached past it for the global. An expectation of 6.22 is readable
    at the uncorrected 0.05 and is not at 0.05/44."""
    loose = _overlap_at(expected=6.22, n_observed=2, alpha=0.05)
    tight = _overlap_at(expected=6.22, n_observed=2, alpha=0.05 / 44)
    assert loose.ratio_readable
    assert not tight.ratio_readable, (
        "DENV-2's 6.22 expected passes a threshold its own row's alpha "
        "could never clear")


def test_the_permutation_count_must_resolve_the_corrected_alpha():
    """MIN_PERMUTATIONS is an absolute floor and answers to no alpha.
    The real requirement is that an attainable p lies below the
    threshold the verdicts use."""
    assert L.resolvable_steps(2000, 0.05 / 44) == 2
    assert L.resolvable_steps(200, 0.05 / 44) == 0, (
        "a 200-permutation run resolves nothing below the corrected "
        "alpha, and MIN_PERMUTATIONS accepted it")
    assert L.resolvable_steps(200, 0.05) >= 10
    assert L.resolvable_steps(0, 0.05) == 0


def test_resolvable_steps_is_the_count_below_alpha_not_at_it():
    """Boundary. With N+1 = 2001 and alpha exactly 2/2001, the value
    0.000999... is not BELOW alpha by more than rounding, so the count
    is the floor and the test is one-sided on purpose."""
    assert L.resolvable_steps(2000, 2 / 2001) == 2
    assert L.resolvable_steps(2000, 1 / 2001) == 1


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


# ------------------------------------------------- multiplicity correction
def test_the_threshold_is_corrected_for_every_test_the_run_performed(tmp_path):
    """Eight taxa in two directions is sixteen tests, not one.

    At an uncorrected 0.05 a run this size turns up a finding from noise
    roughly every other time, and this module exists to stop the agent
    claiming more than was measured. Producing the false claim itself would
    be the same failure wearing a lab coat.
    """
    taxa = {f"NCBITaxon:{i}": dict(_row(), p_matched_less=0.9)
            for i in range(8)}
    land = _load(tmp_path, taxa)
    assert land.n_tests == 16
    assert land.alpha == pytest.approx(0.05 / 16)
    assert all(t.alpha == land.alpha for t in land.tested())


def test_one_direction_counts_one_test_per_taxon(tmp_path):
    land = _load(tmp_path, {V: _row(), OTHER: _row()})
    assert land.n_tests == 2 and land.alpha == pytest.approx(0.025)


def test_the_correction_is_named_in_the_protocol(tmp_path):
    land = _load(tmp_path, {V: _row(), OTHER: _row()})
    assert "Bonferroni-corrected for 2 tests" in land.protocol
    assert "0.0250" in land.protocol


def test_a_nominal_finding_that_does_not_survive_correction_is_not_claimed(
        tmp_path):
    """West Nile's real shape: 0 of 18 against 3.1, lower tail about 0.045.

    Uncorrected that reads as depletion. Across sixteen tests it is nothing,
    and the agent must not announce it.
    """
    taxa = {f"NCBITaxon:{i}": dict(_row(), p_matched_less=0.9)
            for i in range(7)}
    taxa[V] = dict(_row(n_dep=18, n_obs=0, expected=3.1, p_matched=1.0),
                   p_matched_less=0.045)
    land = _load(tmp_path, taxa)
    t = land.for_virus(V)
    assert t.depletion_tested
    assert not t.depleted
    assert land.depleted() == []
    assert "FEWER" not in guards.pathway_caveat(land, V)


def test_a_finding_that_does_survive_correction_is_claimed(tmp_path):
    """The branch must stay reachable, or the correction is just a mute."""
    taxa = {f"NCBITaxon:{i}": dict(_row(), p_matched_less=0.9)
            for i in range(7)}
    taxa[V] = dict(_row(n_dep=1800, n_obs=40, expected=310.0, p_matched=1.0),
                   p_matched_less=0.0005)
    land = _load(tmp_path, taxa)
    assert land.for_virus(V).depleted
    assert "FEWER dependency factors than chance predicts" in \
        guards.pathway_caveat(land, V)


def test_enrichment_is_judged_against_the_same_corrected_threshold(tmp_path):
    """Both directions, one rule. Two thresholds would be harder to state
    than either and impossible for a reader to hold."""
    taxa = {f"NCBITaxon:{i}": dict(_row(), p_matched_less=0.9)
            for i in range(7)}
    taxa[V] = dict(_row(n_obs=90, expected=44.0, p_matched=0.02),
                   p_matched_less=0.99)
    land = _load(tmp_path, taxa)
    assert not land.for_virus(V).enriched   # 0.02 > 0.05/16
    assert land.enriched() == []


def test_an_empty_landscape_does_not_divide_by_zero(tmp_path):
    """n_tests floors at 1. Without it alpha is 0.05/0 on a release where
    no taxon has a dependency factor, which is a crash in a caveat."""
    land = _load(tmp_path, {})
    assert land.n_tests == 1
    assert land.alpha == pytest.approx(0.05)
    assert land.protocol and land.summary()
    assert guards.pathway_caveat(land, V)


# ------------------------------- a measured null downward is still a result
def _both_tails(tmp_path, **over):
    row = dict(_row(n_obs=2, n_dep=50, expected=7.5, p_matched=0.998),
               p_matched_less=0.0095)
    row.update(over)
    taxa = {V: row}
    taxa.update({f"NCBITaxon:{i}": dict(_row(), p_matched_less=0.5)
                 for i in range(7)})
    return _load(tmp_path, taxa)


def test_the_summary_states_the_null_in_the_depletion_direction(tmp_path):
    """The whole return on measuring the lower tail.

    Before this branch existed the summary went SILENT about depletion the
    moment it became measurable: the "not tested" clause switched off and
    the "clears downward" clause never fired, so a run that had done more
    work told the reader less.
    """
    s = _both_tails(tmp_path).summary()
    assert "Nor depleted in any" in s
    assert "0.0095" in s and V in s
    assert "not tested at all" not in s


def test_the_summary_names_the_corrected_threshold_not_the_constant(tmp_path):
    """It said "none below 0.05" while the protocol said 0.0031."""
    s = _both_tails(tmp_path).summary()
    assert "none below the corrected 0.0031" in s
    assert "none below 0.05" not in s


def test_a_tested_undepleted_taxon_says_so_in_its_own_line(tmp_path):
    t = _both_tails(tmp_path).for_virus(V)
    assert "tested downward too and not depleted" in t.describe()
    assert "0.0095" in t.describe()


def test_the_caveat_states_the_downward_null_for_this_taxon(tmp_path):
    own = guards.pathway_caveat(_both_tails(tmp_path), V).split("Release-wide")[0]
    assert "tested downward as well and is not depleted" in own
    assert "0.0095" in own and "0.0031" in own


def test_nothing_claims_depletion_at_the_real_denv2_figure(tmp_path):
    """0.0095 is below 0.05 and above 0.05/16. It must not be announced."""
    land = _both_tails(tmp_path)
    assert land.depleted() == []
    assert "FEWER" not in guards.pathway_caveat(land, V)


# ------------- the correction answers to every test the run performed
_FILE_N = itertools.count()


def _file(tmp_path, taxa, **extra):
    import json
    d = tmp_path / f"f{next(_FILE_N)}"
    d.mkdir()
    p = d / "pathway_landscape.json"
    doc = {"release": "v0.1", "permutations": 2000, "seed": 0,
           "measured_on": "2026-10-28", "partial": False,
           "taxa": taxa, "uncovered": []}
    doc.update(extra)
    p.write_text(json.dumps(doc))
    return p


def _rows(n):
    return {f"NCBITaxon:{i}": {"n_dependency": 100, "n_observed": 10,
                               "expected": 10.0, "p_naive": 0.5,
                               "p_matched": 0.5, "p_matched_less": 0.5}
            for i in range(n)}


def test_a_second_channels_tests_tighten_the_threshold(tmp_path):
    """The interactome channel landed 2026-10-28 and measures fourteen
    taxa in two directions. Counting only this module's eight leaves the
    bar at 0.0031 when the family is forty-four tests and the bar is
    0.0011 -- and the agent would have called a p of 0.0035 significant
    on a run where nothing was.
    """
    eight = _rows(8)
    alone = L.load(_file(tmp_path, eight))
    assert alone.n_tests == 16
    both = L.load(_file(tmp_path, eight, n_tests=44))
    assert both.n_tests == 44
    assert both.alpha < alone.alpha
    assert abs(both.alpha - 0.05 / 44) < 1e-12


def test_a_file_cannot_weaken_the_threshold_by_under_reporting(tmp_path):
    """Taken as a FLOOR. A file claiming fewer tests than this module can
    see must not buy a more permissive bar."""
    land = L.load(_file(tmp_path, _rows(8), n_tests=2))
    assert land.n_tests == 16


def test_an_old_file_with_no_count_still_corrects_for_its_own_rows(tmp_path):
    """Files written before the writer recorded a total."""
    assert L.load(_file(tmp_path, _rows(8))).n_tests == 16


def test_the_protocol_sentence_states_the_whole_family(tmp_path):
    """A reader checking the threshold has to see the number it came
    from, or the correction is unauditable."""
    land = L.load(_file(tmp_path, _rows(8), n_tests=44))
    assert "44 tests" in land.protocol
    assert "0.0011" in land.protocol


def test_every_taxon_carries_the_corrected_threshold(tmp_path):
    """Not just the Landscape. Each verdict is judged against it, and a
    row holding the uncorrected 0.05 would decide on its own."""
    land = L.load(_file(tmp_path, _rows(8), n_tests=44))
    for row in land.taxa.values():
        assert row.n_tests == 44
        assert abs(row.alpha - 0.05 / 44) < 1e-12


# -------------------- a virus whose labels cannot be reached from anywhere
def _routes(**kw):
    kw.setdefault("label", "Yellow fever")
    return guards.VirusRoutes(**kw)


def test_a_virus_with_no_feature_route_refuses_a_ranking():
    """Yellow fever holds 1,030 held-out antiviral measurements in this
    release and no INHIBITS edge, no CRISPR host factor and no
    virus-host interaction. No metapath terminates at it.

    The agent would have returned an EMPTY RANKING, and an empty
    ranking reads as "no compound looks promising". The difference
    between finding nothing and being unable to look is the whole
    distinction this module keeps, and it was missing in the one place
    the user reads.
    """
    v = guards.check(guards.QueryKind.TRIAGE, "NCBITaxon:11089",
                     routes=_routes(n_labels=1030))
    assert not v.allowed
    assert "1,030" in v.reason
    assert "none can be reached" in v.reason
    # THE CONTRAST, not a banned phrase. The refusal legitimately says
    # "not because no compound looks promising, but because none can be
    # reached", so banning that wording fails on the sentence drawing
    # the distinction. Assert the distinction is drawn.
    assert "not because" in v.reason


@pytest.mark.parametrize("kind", [guards.QueryKind.EXPLAIN,
                                  guards.QueryKind.TRIAGE,
                                  guards.QueryKind.PATHWAY])
def test_every_routed_question_is_refused(kind):
    assert not guards.check(kind, "NCBITaxon:11089",
                            routes=_routes(n_labels=312)).allowed


@pytest.mark.parametrize("kind", [guards.QueryKind.PROFILE,
                                  guards.QueryKind.EVIDENCE])
def test_a_lookup_about_the_same_virus_still_works(kind):
    """PROFILE and EVIDENCE report what is STORED about a node, which
    exists whether or not a metapath reaches it. Refusing those would
    be refusing to read the graph."""
    assert guards.check(kind, "NCBITaxon:11089",
                        routes=_routes(n_labels=1030)).allowed


def test_one_route_is_enough_to_answer():
    """DENV-1 reaches its 349 labels through virus-host interactions
    alone. Thin, and not nothing."""
    v = guards.check(guards.QueryKind.TRIAGE, "NCBITaxon:11053",
                     routes=_routes(label="DENV-1", n_labels=349, ppi=53))
    assert v.allowed


def test_a_thin_route_is_answered_with_a_caveat():
    """JEV reaches 46 labels through seven interactions. True, and it
    reads exactly like SARS-CoV-2's 8,316 without a word beside it."""
    v = guards.check(guards.QueryKind.TRIAGE, "NCBITaxon:11072",
                     routes=_routes(label="JEV", n_labels=46, ppi=7))
    assert v.allowed
    text = " ".join(v.warnings)   # Caveat subclasses str
    assert "7 edges" in text and "JEV" in text
    assert "should not be read as a comparison" in text


def test_one_large_route_is_not_thin():
    """SARS-CoV has no CRISPR at all and 1,422 interactions. Having
    few CHANNELS is not the same as having few edges."""
    v = guards.check(guards.QueryKind.TRIAGE, "NCBITaxon:694009",
                     routes=_routes(label="SARS-CoV", n_labels=250,
                                    inhibits=240, ppi=1422))
    assert not any(w.key == "thin_route" for w in v.warnings)


def test_a_large_route_beside_two_small_ones_is_not_thin():
    """West Nile: inhibits(106) + crispr(18) + ppi(6). Its SMALLEST
    channel is below the floor and its largest is four times above it.

    SARS-CoV cannot discriminate here -- both of its channels clear the
    floor, so min and max agree and a test using it passes either way.
    This is the shape that separates them: thin means every route is
    small, so the answer rests on nothing, not that one route happens
    to be. Scoring a compound through the 106 inhibits edges is not made
    uncertain by six interactions sitting beside them.
    """
    v = guards.check(guards.QueryKind.TRIAGE, "NCBITaxon:11082",
                     routes=_routes(label="West Nile", n_labels=594,
                                    inhibits=106, crispr=18, ppi=6))
    assert v.allowed
    assert not any(w.key == "thin_route" for w in v.warnings)


def test_a_release_where_every_route_is_small_is_thin():
    """The same three channels, all below the floor. The caveat is on
    the SIZE of the biggest route, not the number of channels."""
    v = guards.check(guards.QueryKind.TRIAGE, "NCBITaxon:11082",
                     routes=_routes(label="West Nile", n_labels=594,
                                    inhibits=11, crispr=18, ppi=6))
    assert any(w.key == "thin_route" for w in v.warnings)


def test_no_routes_object_means_the_old_behaviour():
    """A caller that cannot compute routes must not get a refusal it
    cannot justify."""
    assert guards.check(guards.QueryKind.TRIAGE, "NCBITaxon:11089").allowed


def test_the_refusal_never_quotes_a_label_as_evidence():
    """n_labels is held-out ground truth. It is said to explain what
    exists, and must never become an answer."""
    v = guards.check(guards.QueryKind.TRIAGE, "NCBITaxon:11089",
                     routes=_routes(n_labels=1030))
    assert "held out" in v.reason


# ------------------- the permutation count against this run's own alpha
# S1. Five of these exist because mutation testing found five gaps at
# once: the alpha guard, the zero-step refusal, the underpowered flag,
# the protocol note, and the user-facing sentence in guards. Each was
# code with no test behind it.
def test_an_alpha_outside_the_unit_interval_is_refused():
    """A threshold of 0 or 1 is not a threshold, and -ln of either is
    not a readability bound. Silently returning inf or 0.0 would pass
    a nonsense alpha straight into a verdict."""
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            L.min_expected_hits(bad)


def test_a_run_that_cannot_reject_at_all_is_refused(tmp_path):
    """200 permutations clears MIN_PERMUTATIONS exactly and resolves
    ZERO attainable p-values below a corrected 0.001136. The old floor
    was derived against the uncorrected 0.05 and accepted it, so the
    agent would have quoted a null test that could not have rejected at
    any observation."""
    land = L.load(_file(tmp_path, _rows(8), permutations=200, n_tests=44),
                  release="v0.1")
    assert not land.usable
    assert "could not have rejected" in land.refused_reason
    assert "200" in land.refused_reason
    # and it says how many would do
    assert "8,799" in land.refused_reason, land.refused_reason


def test_the_real_run_is_usable_and_flagged_underpowered(tmp_path):
    """2,000 permutations resolves TWO steps below 0.001136. Not
    refused: failing to reject is still failing to reject, and
    "nothing survived correction" is the finding being published.
    Flagged, because the negative result is weaker than its four
    decimal places suggest."""
    land = L.load(_file(tmp_path, _rows(8), permutations=2000, n_tests=44),
                  release="v0.1")
    assert land.usable
    assert land.resolvable_steps == 2
    assert land.underpowered


def test_a_well_powered_run_is_not_flagged(tmp_path):
    land = L.load(_file(tmp_path, _rows(8), permutations=100_000,
                        n_tests=44), release="v0.1")
    assert land.usable and not land.underpowered
    assert land.resolvable_steps >= L.MIN_RESOLVABLE_STEPS


def test_the_protocol_states_the_resolution_when_it_is_the_limit(tmp_path):
    """The protocol string is what the agent repeats to a user. Quoting
    a four-decimal alpha beside a permutation count that cannot resolve
    it is the false precision this module is written against."""
    land = L.load(_file(tmp_path, _rows(8), permutations=2000, n_tests=44),
                  release="v0.1")
    proto = land.protocol
    assert "2 attainable p-value(s)" in proto, proto
    assert "could not reject, not that the effect is small" in proto
    assert "8,799 permutations would give 10" in proto
    good = L.load(_file(tmp_path, _rows(8), permutations=100_000,
                        n_tests=44), release="v0.1")
    assert "attainable" not in good.protocol, good.protocol


# ----------------------- the minimum-detectable-effect table (S2, S4)
# The table did not exist anywhere in the project, and "nothing survived
# correction" was being published without it. Every row of the one now in
# docs/scope_and_limits.md is re-derived here from the results file, so a
# re-run with more permutations or another taxon fails until the document
# is regenerated.
_ROOT = Path(__file__).resolve().parents[1]
SCOPE = _ROOT / "docs" / "scope_and_limits.md"
RESULTS = _ROOT / "data" / "results" / "pathway_landscape.json"


def test_the_minimum_detectable_ratio_is_the_poisson_tail_it_claims():
    """The smallest k whose upper tail falls below alpha. Checked against
    a direct sum at a size where the factorial still fits."""
    import math as _m
    alpha = 0.05 / 44
    for m in (0.84, 2.53, 6.22, 7.43, 35.33):
        k, ratio, bind = L.minimum_detectable_ratio(m, alpha, 2000)
        assert bind == "poisson"
        assert abs(ratio - k / m) < 1e-12
        # m bound as a default: the lambda is called in this same iteration,
        # so late binding is harmless here, but binding it says so (ruff B023).
        tail = lambda j, m=m: 1.0 - sum(_m.exp(-m) * m ** i / _m.factorial(i)
                                        for i in range(j))
        assert tail(k) < alpha, (m, k)
        assert tail(k - 1) >= alpha, f"{k} is not the SMALLEST for m={m}"


def test_a_large_expectation_does_not_overflow():
    """m = 425 overflows a factorial, which is why the implementation
    iterates the pmf. A table that cannot describe SARS-CoV-2 is not a
    table of this family."""
    k, ratio, bind = L.minimum_detectable_ratio(425.14, 0.05 / 44, 2000)
    assert bind == "poisson" and k == 490
    assert abs(ratio - 1.15) < 0.01


def test_the_permutation_floor_can_bind_instead_of_the_tail():
    """TWO THINGS CAN BIND. Where the grid cannot express a p below
    alpha, no count is detectable however extreme -- an infinite minimum
    detectable effect caused by the run's configuration, not by the
    biology. Reporting only the Poisson tail would hide it."""
    k, ratio, bind = L.minimum_detectable_ratio(10.0, 0.05 / 44,
                                                permutations=200)
    assert bind == "floor" and ratio == float("inf")
    k2, _, bind2 = L.minimum_detectable_ratio(10.0, 0.05 / 44,
                                              permutations=100_000)
    assert bind2 == "poisson"


def test_an_alpha_outside_the_unit_interval_is_refused_here_too():
    for bad in (0.0, 1.0, -1.0):
        with pytest.raises(ValueError):
            L.minimum_detectable_ratio(10.0, bad)


def test_the_family_covers_every_test_the_correction_charges_for():
    """The module corrects alpha for the interactome channel and loaded
    none of it, so fourteen of the twenty-two rows paying the
    multiple-testing price were invisible to the module charging it."""
    if not RESULTS.exists():
        pytest.skip("no results file")
    land = L.load(RESULTS, release="v0.1")
    assert land.usable, land.refused_reason
    assert len(land.family) * 2 == land.n_tests, (
        f"{len(land.family)} family rows against {land.n_tests} tests")
    assert {r.channel for r in land.family} == {"pathway", "interactome"}
    # and they are NOT results: for_virus must still answer from `taxa`
    assert len(land.taxa) < len(land.family)
    for r in land.family:
        assert not hasattr(r, "p_matched"), (
            "a family row carries a p-value and will be read as a result")


def test_the_powerless_depletion_tests_are_named():
    """SIX, AND THE LITERAL IS PINNED DELIBERATELY. The registration lists
    these rows by name under a DECISION, so a count that drifted silently
    would make the document's list wrong without anything failing. The pair
    moves together or this test does its job -- which it did on 2026-11-10,
    when B1's larger drugged proteome lifted `pathway DENV-2` from 6.40
    expected over the 6.78 bar and this read 5 for an afternoon. B1 was
    declined and the sixth row is back; the episode is why the registration
    now fixes the admission policy as its third DECISION."""
    if not RESULTS.exists():
        pytest.skip("no results file")
    land = L.load(RESULTS, release="v0.1")
    dead = land.powerless_depletion()
    assert len(dead) == 6, [(r.channel, r.virus, r.expected) for r in dead]
    for r in dead:
        assert r.expected <= L.min_expected_hits(r.alpha)
    for r in land.family:
        if r not in dead:
            assert r.expected > L.min_expected_hits(r.alpha)


def test_every_row_of_the_published_table_is_re_derived():
    """The document is the deliverable S4 asks for, so it is checked
    against the run rather than trusted."""
    if not RESULTS.exists():
        pytest.skip("no results file")
    import yaml
    reg = {str(v["taxon"]): v["name"]
           for v in yaml.safe_load((_ROOT / "config/viruses.yaml").read_text())
           ["viruses"]}
    land = L.load(RESULTS, release="v0.1")
    text = SCOPE.read_text()
    for r in land.family:
        name = reg.get(r.virus.split(":")[1], r.virus)
        k, ratio, _ = r.minimum_detectable
        row = (f"| {r.channel} | {name} | {r.expected:.2f} | "
               f"{r.n_observed} | {k} hits = **{ratio:.2f}x** |")
        assert row in text, f"missing or stale table row:\n  {row}"
    # DERIVED, NOT TYPED, AND IT CAUGHT ONE IMMEDIATELY. These two counts
    # read "Six" and "Seven" as literals here and in the document, so the
    # test agreed with the sentence it was checking rather than with
    # pathway_landscape.json -- the pairing this file exists to replace.
    # Derived on 2026-11-10: the first count is right (6) and THE SECOND
    # WAS WRONG (8, not 7). The eighth row is pathway HCoV-NL63 at 2.4513
    # and the ninth is 1.5878, so this is not a rounding boundary; the
    # document had simply been left behind by a rebuild that moved an
    # expectation, and the literal here made that invisible.
    words = {5: "Five", 6: "Six", 7: "Seven", 8: "Eight", 9: "Nine",
             10: "Ten", 11: "Eleven", 12: "Twelve"}
    n_dead = len(land.powerless_depletion())
    n_coarse = sum(1 for r in land.family
                   if r.minimum_detectable[1] >= 2.4)
    assert f"{words[n_dead]} of the twenty-two downward tests" in text, (
        f"{n_dead} rows have no rejection region; the document says "
        f"otherwise")
    assert f"{words[n_coarse]} of the twenty-two need a 2.4x excess" in text, (
        f"{n_coarse} rows need a 2.4x excess or more; the document says "
        f"otherwise")
    assert f"{len(land.family) * 2}" in text


def test_the_coarse_count_is_the_same_in_BOTH_documents(prereg):
    """THE SAME FIGURE LIVED IN TWO FILES AND ONLY ONE WAS CHECKED.
    scope_and_limits.md states it as a sentence and the registration states
    it as a limitation -- "eight of 22 rows need a 2.4x excess or more
    before they could" -- and when the derived count moved from 7 to 8 on
    2026-11-10, the test above moved one document and said nothing about
    the other. A second copy of a figure with no check on it is the same
    defect as a literal, one file further away."""
    if not RESULTS.exists():
        pytest.skip("no results file")
    land = L.load(RESULTS, release="v0.1")
    n_coarse = sum(1 for r in land.family if r.minimum_detectable[1] >= 2.4)
    numbers = {6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}
    flat = " ".join(prereg.split())
    assert f"{numbers[n_coarse]} of 22 rows need a 2.4x excess" in flat, (
        f"{n_coarse} rows need a 2.4x excess or more; the registration's "
        f"limitations section states a different count")


# --------------------- the seed scheme a file was written under (S6)
def test_a_file_with_no_seed_scheme_is_reported_as_the_old_one(tmp_path):
    """ABSENCE IS INFORMATIVE AND IS NOT DEFAULTED AWAY. A file written
    under the shared-generator scheme records no scheme at all, and the
    two are otherwise identical while differing in every p-value. The
    committed v0.1 results file is one of these."""
    land = L.load(_file(tmp_path, _rows(8), permutations=2000, n_tests=44),
                  release="v0.1")
    assert land.seed_scheme == ""
    proto = land.protocol
    assert "single generator in row order" in proto, proto
    assert "adding a taxon moves every downstream p-value" in proto


def test_a_file_with_a_seed_scheme_does_not_carry_the_warning(tmp_path):
    land = L.load(_file(tmp_path, _rows(8), permutations=2000, n_tests=44,
                        seed_scheme="per-row: Random(f'{seed}|{channel}|"
                                    "{virus}')"),
                  release="v0.1")
    assert land.seed_scheme
    assert "single generator" not in land.protocol


def test_the_committed_results_file_predates_the_fix():
    """Stated as a test so a regeneration is visible: once `make
    reproduce` is re-run the file WILL carry a scheme and this flips.
    It is not an assertion that the old file is correct -- it is a
    record of which scheme the published figures came from."""
    if not RESULTS.exists():
        pytest.skip("no results file")
    import json
    doc = json.loads(RESULTS.read_text())
    scheme = doc.get("seed_scheme")
    if scheme:
        assert "per-row" in scheme, scheme
    else:
        land = L.load(RESULTS, release="v0.1")
        assert "single generator in row order" in land.protocol


# ------------------------ the registration, checked against the code (S5)
# A pre-registration that drifts from the code it registers is worse than
# none: it reads as a commitment and documents a different study. Every
# constant the document quotes is re-derived here from the thing it
# describes. The DECISION items are deliberately not pinned -- they are
# proposals awaiting a signature, and a test asserting a proposal would
# make the signature meaningless.
PREREG = _ROOT / "docs" / "prereg_pathway_overlap.md"


@pytest.fixture(scope="module")
def prereg() -> str:
    assert PREREG.exists(), "the registration is missing"
    return " ".join(PREREG.read_text().split())


def test_it_says_plainly_that_it_is_not_a_registration_yet(prereg):
    """The other four open by saying they were written BEFORE the code.
    This one cannot, and the honest version of that is a sentence rather
    than a silence."""
    assert "UNSIGNED DRAFT" in prereg
    assert "must not be executed before then" in prereg
    assert "unregistered and its figures are exploratory" in prereg


def test_the_quoted_release_matches_the_census(prereg):
    """DERIVED FROM THE CENSUS, not typed beside it.

    This hard-coded "392,907 edges" and the two content hashes of the
    2026-10-28 release. Its whole purpose is to catch the registration
    quoting a release that no longer exists -- and when the 2026-11-07
    rebuild moved all four figures, the test failed because ITS copy was
    stale, not because it had detected anything. A check written with the
    same literals it is checking cannot survive the event it exists for.

    The census is generated from the release by scripts/write_kg_report.py
    and audited by `make kg-report-check`, so it is the right source.
    """
    import re
    census = _ROOT / "docs" / "kg_release_census.txt"
    if not census.exists():
        pytest.skip("census absent")
    body = census.read_text()

    m = re.search(r"([\d,]+) nodes[^\d]*([\d,]+) edges", body)
    assert m, "the census no longer states a node and edge count"
    for token in (f"{m.group(1)} nodes", f"{m.group(2)} edges"):
        assert token in prereg, (
            f"the census says {token} and the registration does not quote "
            f"it; the registration names a release that no longer exists")

    shas = sorted(set(re.findall(r"\b[0-9a-f]{16}\b", body)))
    assert len(shas) >= 2, f"the census states {len(shas)} content hash(es)"
    for sha in shas:
        assert sha in prereg, (
            f"the census carries {sha} and the registration does not")


def test_the_declared_alpha_is_the_one_the_module_computes(prereg):
    if not RESULTS.exists():
        pytest.skip("no results file")
    land = L.load(RESULTS, release="v0.1")
    assert f"0.05/{land.n_tests}" in prereg
    assert f"{land.alpha:.6f}" in prereg
    assert f"{len(land.family)} rows x 2 = {land.n_tests}" in prereg


def test_the_declared_readability_gate_is_the_derived_one(prereg):
    """It was a constant 3.0 from the uncorrected alpha. The document
    must register the DERIVATION, not a number that goes stale when the
    family changes."""
    assert "expected > -ln(alpha)" in prereg
    assert "exp(-3) = 0.0498 < 0.05" in prereg, (
        "the document should name the wrong derivation it replaces")
    if RESULTS.exists():
        land = L.load(RESULTS, release="v0.1")
        assert f"{land.min_expected_hits:.2f}" in prereg


def test_the_declared_seed_scheme_is_the_implemented_one(prereg):
    import sys
    from importlib import import_module
    sys.path.insert(0, str(_ROOT / "scripts"))
    pl = import_module("pathway_landscape")
    assert 'Random(f"{seed}|{channel}|{virus}")' in prereg
    assert "channel" in pl.SEED_SCHEME and "virus" in pl.SEED_SCHEME


def test_every_powerless_row_is_named_in_the_registration(prereg):
    """The six are the subject of a DECISION, so the document has to
    list them rather than cite a count. Five for one afternoon on
    2026-11-10, under a graph change that was reverted."""
    if not RESULTS.exists():
        pytest.skip("no results file")
    import yaml
    reg = {str(v["taxon"]): v["name"]
           for v in yaml.safe_load((_ROOT / "config/viruses.yaml").read_text())
           ["viruses"]}
    land = L.load(RESULTS, release="v0.1")
    dead = land.powerless_depletion()
    assert len(dead) == 6
    for r in dead:
        name = reg.get(r.virus.split(":")[1], r.virus)
        assert f"{r.channel} {name} {r.expected:.2f}" in prereg, (
            f"{r.channel}/{name} at {r.expected:.2f} is not listed")


def test_the_permutation_arithmetic_is_right(prereg):
    """The document quotes three permutation counts and a step count.
    Each is derived, so each is checked."""
    if not RESULTS.exists():
        pytest.skip("no results file")
    land = L.load(RESULTS, release="v0.1")
    import math as _m
    hard = int(_m.ceil(1 / land.alpha))
    ten = int(_m.ceil(10 / land.alpha)) - 1
    assert f"N+1 > 1/alpha = {hard}" in prereg
    assert f"N >= {ten:,}" in prereg
    assert f"{L.resolvable_steps(100_000, land.alpha)} steps" in prereg


def test_the_smallest_p_and_its_multiple_of_alpha_are_right(prereg):
    """P2 rests on the smallest p being far above the threshold, so the
    prediction is only as good as that number."""
    if not RESULTS.exists():
        pytest.skip("no results file")
    import json
    doc = json.loads(RESULTS.read_text())
    land = L.load(RESULTS, release="v0.1")
    ps = [r[d] for src in (doc["taxa"], doc["interactome"])
          for r in src.values()
          for d in ("p_matched", "p_matched_less") if r.get(d) is not None]
    lo = min(ps)
    assert f"{lo:.6f}" in prereg, f"smallest p is {lo:.6f}"
    assert f"{lo / land.alpha:.1f}x the Bonferroni threshold" in prereg


def test_the_three_decisions_are_marked_as_such(prereg):
    """Three decisions and no more. A fourth appearing silently would mean
    the document started deciding things on its own."""
    # THE HEADER'S COUNT AND THE MARKERS MUST AGREE. The first draft
    # promised three and marked two -- the third belonged to the
    # hard-negatives analysis and had been dropped from the body without
    # the header being corrected. This test caught that.
    #
    # THE THIRD DECISION ARRIVED 2026-11-10: which compounds the drug-target
    # query admits (--phase-gate). It is here because the verdict was shown
    # to depend on it -- interactome/NCBITaxon:2697049 reads p_matched 0.235
    # under `floor` and 0.000500 under `relaxed` -- and a registration that
    # fixes the permutation count to three digits while leaving the operand
    # undefined is incomplete in the way the seeding scheme already was.
    # It stayed after fix-list B1 was declined, which is the point: the
    # registered value is `floor` BECAUSE the alternative was measured, not
    # because nobody tried it.
    assert "Three decisions are flagged **DECISION**" in prereg
    assert prereg.count("**DECISION**") == 4, (
        "one in the header plus one per decision; the count moved")
