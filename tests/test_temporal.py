"""Gate for the temporal split: the future is removed from the training graph,
and the audits fail loudly when it is not.
"""
import json

import pytest

from kgav.schema import load_schema
from kgav.temporal import audit, build_split, edge_year, write_split

V = "NCBITaxon:2697049"
OLD = "INCHIKEY:OLD"
NEW = "INCHIKEY:NEW"


def _e(s, p, o, date, q=None, src="s", tier=1):
    return {"subject": s, "predicate": p, "object": o, "qualifiers": q or {},
            "primary_knowledge_source": src, "evidence_tier": tier,
            "first_asserted_date": date}


def _label(drug, date, relation="="):
    return _e(drug, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", V, date,
              {"assay_type": "cell_based_antiviral", "ec50_nm": 100.0,
               "relation": relation, "unquantified": True})


@pytest.fixture
def held():
    return load_schema().held_out_predicates()


@pytest.fixture
def release(tmp_path):
    def prot(pid, viral, fam=None):
        p = {"taxon_id": V if viral else "NCBITaxon:9606", "is_viral": viral,
             "sequence_hash": "h", "reviewed": True}
        if fam:
            p["protein_family"] = fam
        return {"id": pid, "class": "Protein", "properties": p}

    def drug(n):
        return {"id": n, "class": "SmallMolecule",
                "properties": {"smiles": "C", "inchikey_skel": n[-2:],
                               "is_approved": True, "salt_collapsed": False,
                               "stereo_collapsed": False}}

    nodes = [{"id": V, "class": "OrganismTaxon",
              "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                             "baltimore_class": "IV", "is_enveloped": True}},
             prot("UniProtKB:NSP5", True, "nsp5"),
             {"id": "KGAV:G", "class": "Gene",
              "properties": {"symbol": "rep", "taxon_id": V, "is_viral": True}},
             drug(OLD), drug(NEW)]
    edges = [
        _e(OLD, "INHIBITS", "UniProtKB:NSP5", "2020-01-01",
           {"assay_type": "biochemical", "ic50_nm": 25.0, "unquantified": True}),
        _e(NEW, "INHIBITS", "UniProtKB:NSP5", "2023-01-01",
           {"assay_type": "biochemical", "ic50_nm": 10.0, "unquantified": True}),
        _e("UniProtKB:NSP5", "ENCODED_BY", "KGAV:G", "2020-01-01"),
        _e("KGAV:G", "BELONGS_TO", V, "2020-01-01"),
        _e(OLD, "CHEMICALLY_SIMILAR_TO", NEW, "1970-01-01",
           {"tanimoto_ecfp4": 0.9}, "infores:kgav-computed", 3),
        _e(NEW, "TARGETS", "UniProtKB:NSP5", "1970-01-01",
           {"direction": "unknown", "assay_type": "binding", "ic50_nm": 5.0}),
        _label(OLD, "2020-06-01"),
        _label(NEW, "2023-06-01"),
        _label(NEW, "2024-01-01", relation=">"),
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return tmp_path


# ------------------------------------------------------------------- parsing
def test_placeholder_date_is_not_a_year():
    assert edge_year({"first_asserted_date": "2021-03-04"}) == 2021
    assert edge_year({"first_asserted_date": "1970-01-01"}) is None
    assert edge_year({}) is None


# --------------------------------------------------------------------- split
def test_future_evidence_is_removed_from_the_training_graph(release, held):
    """THE trap this module exists for. Filtering only the labels leaves a
    compound whose nsp5 activity was published in 2023 with an M1 path built
    from 2023 evidence predicting a 2023 label."""
    s = build_split(release, 2021, held)
    assert s.stats["edge_dropped_future"] == 1
    for e in s.train_edges:
        y = edge_year(e)
        assert y is None or y <= 2021


def test_labels_are_partitioned_by_year(release, held):
    s = build_split(release, 2021, held)
    assert s.train_labels[V] == {OLD}
    assert s.test_labels[V] == {NEW}


def test_censored_labels_are_excluded(release, held):
    """'EC50 > 9999 nM' is a measurement of INACTIVITY; using it as a positive
    would invert the target."""
    s = build_split(release, 2021, held)
    assert s.stats["label_censored"] == 1


def test_held_out_predicates_never_enter_the_training_graph(release, held):
    s = build_split(release, 2021, held)
    assert not [e for e in s.train_edges if e["predicate"] in held]


# ------------------------------------------------------------ undated policy
def test_undated_include_keeps_computed_and_real(release, held):
    s = build_split(release, 2021, held, "include")
    assert s.stats["undated_kept_computed"] == 1   # similarity
    assert s.stats["undated_kept"] == 1            # the dateless TARGETS edge


def test_undated_exclude_drops_everything_dateless(release, held):
    """Conservative: an unknown date may be post-cutoff."""
    s = build_split(release, 2021, held, "exclude")
    assert s.stats["undated_dropped"] == 2
    assert all(edge_year(e) is not None for e in s.train_edges)


def test_undated_computed_only_keeps_similarity_but_not_unknown_dates(release, held):
    s = build_split(release, 2021, held, "computed_only")
    assert s.stats["undated_kept_computed"] == 1
    assert s.stats["undated_dropped"] == 1


def test_unknown_policy_is_rejected(release, held):
    with pytest.raises(ValueError):
        build_split(release, 2021, held, "nonsense")


# -------------------------------------------------------------------- audits
def test_clean_split_passes_the_audits(release, held):
    s = build_split(release, 2021, held)
    problems = audit(s, release, held)
    assert not [p for p in problems if p.startswith(("L3", "L5", "label leak"))]


def test_l3_detects_planted_future_bleed(release, held):
    s = build_split(release, 2021, held)
    s.train_edges.append(_e("X", "INHIBITS", "Y", "2024-01-01"))
    assert any(p.startswith("L3") for p in audit(s, release, held))


def test_label_leak_is_detected(release, held):
    s = build_split(release, 2021, held)
    s.train_edges.append(_label(NEW, "2020-01-01"))
    assert any(p.startswith("label leak") for p in audit(s, release, held))


def test_l5_detects_a_compound_in_both_label_sets(release, held):
    s = build_split(release, 2021, held)
    s.train_labels[V].add(NEW)
    assert any(p.startswith("L5") for p in audit(s, release, held))


def test_recall_ceiling_is_reported_not_hidden(release, held):
    """A test compound with no pre-cutoff evidence is unpredictable by
    construction. Not a leak, but it caps achievable recall."""
    s = build_split(release, 2021, held)
    s.test_labels[V].add("INCHIKEY:GHOST")
    assert any(p.startswith("ceiling") for p in audit(s, release, held))


# --------------------------------------------------------------------- output
def test_written_graph_keeps_only_reachable_nodes(release, held, tmp_path):
    s = build_split(release, 2021, held)
    out = tmp_path / "train"
    write_split(s, release, out)
    ids = {json.loads(x)["id"] for x in (out / "nodes.jsonl").read_text().splitlines()
           if x.strip()}
    edges = [json.loads(x) for x in (out / "edges.jsonl").read_text().splitlines()
             if x.strip()]
    for e in edges:
        assert e["subject"] in ids and e["object"] in ids
    labels = json.loads((out / "TEST_LABELS.json").read_text())
    assert labels["cutoff"] == 2021
    assert labels["test"][V] == [NEW]
