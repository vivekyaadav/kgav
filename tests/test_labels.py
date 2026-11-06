"""Gate for measured negatives: a lower bound is only inactivity if it is high
enough, ambiguity is dropped, and no scorer picks its own test set.
"""
import json
from collections import Counter
from pathlib import Path

import pytest

from kgav import labels as L
from kgav.labels import edge_sources, build_labels, classify, evaluate_against_negatives

V = "NCBITaxon:2697049"


def _act(drug, value, relation="=", date="2020-06-01"):
    return {"subject": drug, "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
            "object": V, "qualifiers": {"assay_type": "cell_based_antiviral",
                                        "ec50_nm": value, "relation": relation},
            "primary_knowledge_source": "s", "evidence_tier": 1,
            "first_asserted_date": date}


@pytest.mark.parametrize("value,relation,expect", [
    (200.0, "=", "active"),
    (50_000.0, "=", "inactive"),
    (50_000.0, ">", "inactive"),
    (50_000.0, ">=", "inactive"),
    # A weak lower bound is NOT inactivity: it may sit below the assay's top
    # concentration.
    (100.0, ">", None),
    (100.0, "<", "active"),
    (50_000.0, "<", None),
])
def test_classification_of_censored_values(value, relation, expect):
    assert classify(_act("d", value, relation)) == expect


def test_edge_without_a_potency_value_decides_nothing():
    e = _act("d", 1.0)
    del e["qualifiers"]["ec50_nm"]
    assert classify(e) is None


def _release(tmp_path, edges):
    (tmp_path / "nodes.jsonl").write_text("")
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return tmp_path


def test_negatives_come_from_censored_measurements(tmp_path):
    r = _release(tmp_path, [_act("A", 200.0), _act("N", 50_000.0, ">")])
    lab = build_labels(r)
    assert lab.positives[V] == {"A"}
    assert lab.negatives[V] == {"N"}


def test_compounds_measured_both_ways_are_dropped_not_assigned(tmp_path):
    """Real disagreement -- different cell lines, strains, labs. Forcing a side
    would manufacture a label the data does not support."""
    r = _release(tmp_path, [_act("X", 100.0), _act("X", 50_000.0, ">")])
    lab = build_labels(r)
    assert "X" not in lab.positives.get(V, set())
    assert "X" not in lab.negatives.get(V, set())
    assert lab.stats["ambiguous_dropped"] == 1


def test_year_window_filters(tmp_path):
    r = _release(tmp_path, [_act("OLD", 200.0, date="2020-01-01"),
                            _act("NEW", 200.0, date="2023-01-01")])
    lab = build_labels(r, year_from=2022)
    assert lab.positives[V] == {"NEW"}
    assert lab.stats["out_of_window"] == 1


def test_undated_labels_are_excluded_from_a_window(tmp_path):
    r = _release(tmp_path, [_act("U", 200.0, date="1970-01-01")])
    assert build_labels(r, year_from=2022).positives == {}


# ------------------------------------------------------------------- scoring
def test_auc_is_half_when_scores_are_uninformative():
    scores = {"a": 1.0, "b": 1.0, "n": 1.0, "m": 1.0}
    m = evaluate_against_negatives(scores, {"a", "b"}, {"n", "m"})
    assert m["auc"] == pytest.approx(0.5)


def test_auc_is_one_when_separation_is_perfect():
    scores = {"a": 2.0, "b": 3.0, "n": 0.1, "m": 0.2}
    assert evaluate_against_negatives(scores, {"a", "b"}, {"n", "m"})["auc"] == 1.0


def test_unreached_compounds_take_the_floor_not_an_exclusion():
    """A metapath scorer only reaches compounds with paths, and measured
    inactives frequently have none. Evaluating on just what a scorer reaches
    lets it choose its own test set -- that is how 3% coverage posts a perfect
    AUC."""
    scores = {"a": 5.0}          # only the active is reached
    m = evaluate_against_negatives(scores, {"a"}, {"n"})
    assert m["evaluable"] and m["n_neg"] == 1
    assert m["coverage_neg"] == 0.0
    assert m["auc"] == 1.0


def test_coverage_is_reported_separately_from_class_size():
    """A high AUC on low coverage means the floor is doing the work."""
    scores = {"a": 5.0, "n": 1.0}
    m = evaluate_against_negatives(scores, {"a", "b"}, {"n", "o"})
    assert m["n_pos"] == 2 and m["coverage_pos"] == 0.5
    assert m["n_neg"] == 2 and m["coverage_neg"] == 0.5


def test_ties_count_as_half():
    scores = {"a": 1.0, "n": 1.0}
    assert evaluate_against_negatives(scores, {"a"}, {"n"})["auc"] == pytest.approx(0.5)


def test_empty_class_is_not_evaluable():
    m = evaluate_against_negatives({"a": 1.0}, {"a"}, set())
    assert m["evaluable"] is False and m["auc"] is None


def test_label_sets_are_ordered_deterministically(tmp_path):
    """Virus order came from iterating a SET of ids, and str hashing is
    randomised per interpreter -- so two runs over one release wrote results
    files that differed textually while every number in them was identical.
    That is indistinguishable from real drift to anything reading the file,
    and it would read as spurious staleness against the provenance stamps."""
    viruses = ["NCBITaxon:2697049", "NCBITaxon:694009", "NCBITaxon:11137",
               "NCBITaxon:31631", "NCBITaxon:277944", "NCBITaxon:1335626"]
    edges = []
    for v in viruses:
        for e in (_act("A", 200.0), _act("N", 50_000.0, ">")):
            edges.append({**e, "object": v})
    lab = build_labels(_release(tmp_path, edges))
    assert list(lab.positives) == sorted(lab.positives)
    assert list(lab.negatives) == sorted(lab.negatives)


# ------------------------------------------------- single-screen evaluation
def _av(subject, source, date="2020-08-18", **quals):
    q = {"assay_type": "cell_based_antiviral"}
    q.update(quals)
    return {"subject": subject, "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
            "object": "NCBITaxon:2697049", "qualifiers": q,
            "primary_knowledge_source": source, "evidence_tier": 1,
            "first_asserted_date": date}


@pytest.fixture
def mixed_sources(tmp_path):
    """One compound per source per class, plus a merged edge."""
    edges = [
        _av("INCHIKEY:CA", "infores:chembl", ec50_nm=100.0, relation="="),
        _av("INCHIKEY:CI", "infores:chembl", ec50_nm=50_000.0, relation=">"),
        _av("INCHIKEY:NA", "infores:ncats-opendata", ec50_nm=900.0, relation="="),
        _av("INCHIKEY:NI", "infores:ncats-opendata", ec50_nm=20_000.0, relation=">"),
        # a merged edge: membership must be tested per component
        _av("INCHIKEY:MM", "infores:chembl|infores:ncats-opendata",
            ec50_nm=200.0, relation="="),
    ]
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in edges))
    return tmp_path


def test_edge_sources_splits_a_merged_source_string():
    assert edge_sources(_av("X", "infores:a|infores:b")) == {"infores:a", "infores:b"}
    assert edge_sources(_av("X", "infores:a")) == {"infores:a"}
    assert edge_sources({"primary_knowledge_source": ""}) == set()


def test_unrestricted_labels_pool_every_source(mixed_sources):
    lab = build_labels(mixed_sources)
    v = "NCBITaxon:2697049"
    assert lab.positives[v] == {"INCHIKEY:CA", "INCHIKEY:NA", "INCHIKEY:MM"}
    assert lab.negatives[v] == {"INCHIKEY:CI", "INCHIKEY:NI"}


def test_restricting_to_one_screen_keeps_only_that_screen(mixed_sources):
    """The point: both classes then come off the same plates, so a scorer
    cannot separate them by recognising which collection they came from."""
    lab = build_labels(mixed_sources, sources={"infores:ncats-opendata"})
    v = "NCBITaxon:2697049"
    assert lab.positives[v] == {"INCHIKEY:NA", "INCHIKEY:MM"}   # MM is merged
    assert lab.negatives[v] == {"INCHIKEY:NI"}
    assert lab.stats["other_source"] == 2


def test_a_merged_edge_matches_either_component(mixed_sources):
    for src in ("infores:chembl", "infores:ncats-opendata"):
        lab = build_labels(mixed_sources, sources={src})
        assert "INCHIKEY:MM" in lab.positives["NCBITaxon:2697049"], src


def test_an_unknown_source_yields_no_labels_rather_than_all(mixed_sources):
    """A typo must empty the evaluation, not silently fall back to pooling."""
    lab = build_labels(mixed_sources, sources={"infores:typo"})
    assert not lab.positives and not lab.negatives
    assert lab.stats["other_source"] == 5


# --------------------------------------------------------------- reach
# MEASURED REACH, moved here from scripts/audit_build.py on 2026-11-02.
# It lived in a script, so the agent -- which cannot import a script --
# could not tell a route that reaches labels from one that reaches none,
# and told five viruses their reliability was unknown when it is zero.
# These tests are about the three states the function must keep apart:
# None for "cannot be known", 0 for "reaches nothing", n for reach.
def _comparison(results, release, pos_by_virus, stale=False, raw=None):
    from kgav.provenance import write_results
    results.mkdir(parents=True, exist_ok=True)
    path = results / "metapath_comparison.json"
    per = raw if raw is not None else {
        v: {"COMBINED": {"pos": n, "pool": max(n, 1)}}
        for v, n in pos_by_virus.items()}
    write_results(path, {"per_virus": per}, release=release)
    if stale:
        doc = json.loads(path.read_text())
        doc["provenance"]["release"]["edges_sha"] = "0" * 16
        path.write_text(json.dumps(doc))
    return results


def _tiny_release(tmp_path):
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text('{"id": "NCBITaxon:1"}\n')
    (rel / "edges.jsonl").write_text(
        '{"subject": "a", "predicate": "P", "object": "b"}\n')
    return rel


def test_measured_reach_reads_a_current_file(tmp_path):
    rel = _tiny_release(tmp_path)
    res = _comparison(tmp_path / "r", rel, {"DENV-1": 0, "DENV-2": 117})
    assert L.measured_reach(res, rel) == {"DENV-1": 0, "DENV-2": 117}


def test_an_absent_file_is_unmeasured_not_zero(tmp_path):
    """None and 0 are OPPOSITE claims. A caller that cannot tell them
    apart asserts that nothing is reachable because it failed to look."""
    rel = _tiny_release(tmp_path)
    assert L.measured_reach(tmp_path / "nothing-here", rel) is None


def test_a_stale_file_is_unmeasured_not_zero(tmp_path):
    """The specific way this project has shipped a wrong figure twice:
    a results file measured on a graph that has since moved. Reported as
    unknown rather than quoted."""
    rel = _tiny_release(tmp_path)
    res = _comparison(tmp_path / "r", rel, {"DENV-1": 0}, stale=True)
    assert L.measured_reach(res, rel) is None, \
        "a stale comparison was quoted as a measurement"


def test_a_file_with_no_integer_reach_is_unmeasured(tmp_path):
    rel = _tiny_release(tmp_path)
    res = _comparison(tmp_path / "r", rel, {},
                      raw={"DENV-1": {"COMBINED": {"pos": None}},
                           "JEV": {"COMBINED": {}},
                           "TBEV": {}})
    assert L.measured_reach(res, rel) is None


def test_a_virus_with_no_integer_reach_is_dropped_not_zeroed(tmp_path):
    """One unreadable row must not become a reach of 0 for that virus
    while the readable rows are kept."""
    rel = _tiny_release(tmp_path)
    res = _comparison(tmp_path / "r", rel, {},
                      raw={"DENV-1": {"COMBINED": {"pos": "0"}},
                           "DENV-2": {"COMBINED": {"pos": 117}}})
    assert L.measured_reach(res, rel) == {"DENV-2": 117}


def test_malformed_json_is_unmeasured(tmp_path):
    rel = _tiny_release(tmp_path)
    res = tmp_path / "r"
    res.mkdir()
    (res / "metapath_comparison.json").write_text("{not json")
    assert L.measured_reach(res, rel) is None


# ------------------------------------------- the inactive floor is MEASURED
# The justification for 10 uM was "the concentration above which a cell-based
# antiviral hit is conventionally not pursued" -- an appeal to convention
# with no citation, and the weakest sentence in labels.py. It is not needed:
# a censored row "EC50 > 20000" is a screen stating where it stopped, so the
# distribution of those bounds MEASURES the assay ceilings in the corpus.
#
# These tests re-derive the docstring's figures from the release. A number in
# a comment with no way back to the data is the failure this repository keeps
# finding; this is the way back.
REAL_RELEASE = Path(__file__).resolve().parents[1] / "data" / "releases" / "v0.1"


def test_the_census_reads_only_censored_rows_of_the_labelled_predicate(tmp_path):
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in [
        # counted: censored, right predicate, has a potency
        {"predicate": L.ACTIVITY_PREDICATE,
         "qualifiers": {"ec50_nm": 20000.0, "relation": ">"}},
        {"predicate": L.ACTIVITY_PREDICATE,
         "qualifiers": {"ic50_nm": 10000.0, "relation": ">="}},
        # not counted: exact value, so it states no ceiling
        {"predicate": L.ACTIVITY_PREDICATE,
         "qualifiers": {"ec50_nm": 500.0, "relation": "="}},
        # not counted: a different predicate is a different population
        {"predicate": "MEASURED_INACTIVE_AGAINST",
         "qualifiers": {"ec50_nm": 20000.0, "relation": ">"}},
        # not counted: censored with no value to be a bound
        {"predicate": L.ACTIVITY_PREDICATE, "qualifiers": {"relation": ">"}},
    ]))
    assert L.censored_bound_census(rel) == Counter({20000.0: 1, 10000.0: 1})


def test_the_census_is_scoped_to_the_predicate_build_labels_reads():
    """A census over a wider set of predicates would be a true statement
    about a population the threshold does not decide."""
    assert L.censored_bound_census.__defaults__[0] == L.ACTIVITY_PREDICATE


@pytest.mark.skipif(not REAL_RELEASE.exists(), reason="no release on disk")
def test_the_floor_admits_almost_all_censored_evidence():
    """The claim the docstring now makes instead of citing a convention:
    10 uM admits 98.1% of the censored rows and rejects the 1.9% whose
    bound is too weak to prove anything."""
    c = L.censored_bound_census(REAL_RELEASE)
    total = sum(c.values())
    assert total > 10_000, f"only {total} censored rows; the census is empty"
    above = sum(n for v, n in c.items() if v >= L.DEFAULT_INACTIVE_ABOVE_NM)
    assert above / total > 0.95, (
        f"the floor now admits only {above / total:.1%} of censored evidence; "
        f"the docstring says 98.1% and would need re-deriving")
    below = total - above
    assert below > 0, (
        "no censored row falls below the floor, so the floor is rejecting "
        "nothing and the rule it encodes is untested by this corpus")


@pytest.mark.skipif(not REAL_RELEASE.exists(), reason="no release on disk")
def test_the_floor_sits_at_or_below_the_modal_assay_ceiling():
    """The reason 10 uM rather than 20: the corpus supports 10-20 uM and the
    floor takes the conservative end. A floor ABOVE the modal ceiling would
    discard the evidence of the commonest assay in the corpus."""
    c = L.censored_bound_census(REAL_RELEASE)
    mode_bound, _ = c.most_common(1)[0]
    assert L.DEFAULT_INACTIVE_ABOVE_NM <= mode_bound, (
        f"the floor ({L.DEFAULT_INACTIVE_ABOVE_NM:,.0f} nM) is above the "
        f"modal assay ceiling ({mode_bound:,.0f} nM), so the commonest "
        f"screen in the corpus can no longer prove a compound inactive")


def test_the_two_cuts_cannot_overlap():
    """CONSISTENCY, not convention. A censored row "> X" proves inactivity
    only if X >= active_below_nm; below that the true value could sit
    between X and the active cut, which would make the compound ACTIVE. Set
    the floor lower and one compound can be classified both ways."""
    assert L.DEFAULT_INACTIVE_ABOVE_NM >= L.DEFAULT_ACTIVE_BELOW_NM

    # The overlap, demonstrated rather than asserted.
    censored = {"qualifiers": {"ec50_nm": 2000.0, "relation": ">"}}
    exact = {"qualifiers": {"ec50_nm": 5000.0, "relation": "="}}
    assert classify(censored, 1000.0, 10000.0) == "inactive"
    assert classify(exact, 1000.0, 10000.0) == "active"
    # Tied, the weak bound is correctly undecidable instead.
    assert classify(censored, 10000.0, 10000.0) is None


def test_a_dead_band_is_expressible_and_stays_consistent():
    """The open choice the docstring records: active <= 1 uM, inactive >=
    10 uM. classify takes the two cuts independently and the inequality
    still holds, so a gap is available without a code change -- only the
    CLI ties them, and --active-below-nm does not exist."""
    mid = {"qualifiers": {"ec50_nm": 5000.0, "relation": "="}}
    assert classify(mid, 10_000.0, 1_000.0) == "inactive"
    weak = {"qualifiers": {"ec50_nm": 900.0, "relation": "="}}
    assert classify(weak, 10_000.0, 1_000.0) == "active"
    bound = {"qualifiers": {"ec50_nm": 5000.0, "relation": ">"}}
    assert classify(bound, 10_000.0, 1_000.0) is None


@pytest.mark.skipif(not REAL_RELEASE.exists(), reason="no release on disk")
def test_the_active_cut_is_not_claimed_to_be_derived():
    """The two constants are equal today, and the docstring says that is
    inheritance rather than an argument. If someone later derives the active
    cut, this test should be replaced by one that checks the derivation --
    not deleted."""
    src = (Path(L.__file__)).read_text()
    assert "ONLY HALF DERIVED" in src
    assert "no distribution in this corpus answers it" in src
