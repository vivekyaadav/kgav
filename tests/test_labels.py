"""Gate for measured negatives: a lower bound is only inactivity if it is high
enough, ambiguity is dropped, and no scorer picks its own test set.
"""
import json

import pytest

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
