"""Gate for the cold-start layer: similarity is bipartite, capped, and the
threshold cannot drift below what the schema declares.
"""
import json

import pytest

from kgav.emit import Emit
from kgav.schema import load_schema
from kgav.similarity import (
    DEFAULT_THRESHOLD,
    build,
    fingerprints,
    load_compounds,
    similarity_distribution,
    similarity_edges,
)

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
ASPIRIN_ME = "CC(=O)Oc1ccccc1C(=O)OC"
CYCLOHEXANE = "C1CCCCC1"


def test_fingerprints_count_failures_rather_than_skipping():
    """A compound that cannot be fingerprinted is unreachable by any method,
    so the number belongs in the report."""
    fps, stats = fingerprints({"a": ASPIRIN, "b": "not_a_smiles", "c": ""})
    assert set(fps) == {"a"}
    assert stats["fingerprinted"] == 1
    assert stats["unparseable"] == 1
    assert stats["no_smiles"] == 1


def test_similar_structures_score_high_and_dissimilar_low():
    fps, _ = fingerprints({"asp": ASPIRIN, "me": ASPIRIN_ME, "cyc": CYCLOHEXANE})
    edges = similarity_edges({"me": fps["me"], "cyc": fps["cyc"]},
                             {"asp": fps["asp"]}, threshold=0.0)
    sims = {(q, r): s for q, r, s in edges}
    assert sims[("me", "asp")] > sims[("cyc", "asp")]
    assert sims[("cyc", "asp")] < 0.2


def test_self_pairs_are_excluded():
    fps, _ = fingerprints({"asp": ASPIRIN})
    assert similarity_edges(fps, fps, threshold=0.0) == []


def test_top_k_caps_dense_scaffold_series():
    """Without a cap, a compound in a dense series acquires hundreds of
    near-identical neighbours and DWPC reads one scaffold counted many times
    as overwhelming evidence."""
    refs = {f"r{i}": ASPIRIN for i in range(30)}
    ref_fps, _ = fingerprints(refs)
    q_fps, _ = fingerprints({"q": ASPIRIN_ME})
    assert len(similarity_edges(q_fps, ref_fps, threshold=0.0, top_k=5)) == 5
    assert len(similarity_edges(q_fps, ref_fps, threshold=0.0, top_k=None)) == 30


def test_threshold_filters():
    fps, _ = fingerprints({"asp": ASPIRIN, "cyc": CYCLOHEXANE})
    assert similarity_edges({"cyc": fps["cyc"]}, {"asp": fps["asp"]},
                            threshold=0.9) == []


def test_distribution_reports_reach_per_cut():
    fps, _ = fingerprints({"asp": ASPIRIN, "me": ASPIRIN_ME})
    dist = similarity_distribution({"me": fps["me"]}, {"asp": fps["asp"]})
    assert dist[0.5] == 1
    assert dist[0.95] == 0


def test_edges_are_tier_three_and_undated():
    """A computed edge has no assertion date. Stamping it with the build date
    would make it look newer than every real fact and corrupt a temporal
    split."""
    em = Emit()
    build(em, {"q": ASPIRIN_ME, "r": ASPIRIN}, {"q"}, {"r"},
          "infores:kgav-computed", threshold=0.0)
    assert em.edges
    for e in em.edges:
        assert e["evidence_tier"] == 3
        assert e["first_asserted_date"] == "1970-01-01"
        assert e["primary_knowledge_source"] == "infores:kgav-computed"


def test_default_threshold_matches_the_schema_floor():
    """The schema pins the cut. If these drift apart the ingest emits edges the
    validator rejects."""
    s = load_schema()
    floor = next((ec.qualifiers.get("tanimoto_ecfp4", {}).get("min")
                  for ec in s.edge_classes
                  if ec.predicate == "CHEMICALLY_SIMILAR_TO"), None)
    assert floor is not None
    assert DEFAULT_THRESHOLD >= floor


def test_output_validates_against_schema(tmp_path):
    nodes = [{"id": f"INCHIKEY:{n}", "class": "SmallMolecule",
              "properties": {"smiles": smi, "inchikey_skel": n, "is_approved": True,
                             "salt_collapsed": False, "stereo_collapsed": False}}
             for n, smi in (("Q", ASPIRIN_ME), ("R", ASPIRIN))]
    em = Emit()
    build(em, {n["id"]: n["properties"]["smiles"] for n in nodes},
          {"INCHIKEY:Q"}, {"INCHIKEY:R"}, "infores:kgav-computed", threshold=0.6)
    violations = load_schema().validate_batch(nodes, em.edges)
    # 0.667 sits below the schema floor of 0.7, so this MUST be rejected --
    # the guard is what stops a looser cut slipping in unnoticed.
    assert any(v.code == "RANGE" for v in violations)


def test_load_compounds_returns_only_small_molecules(tmp_path):
    nodes = [
        {"id": "INCHIKEY:A", "class": "SmallMolecule", "properties": {"smiles": ASPIRIN}},
        {"id": "UniProtKB:P1", "class": "Protein", "properties": {}},
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    assert load_compounds(tmp_path) == {"INCHIKEY:A": ASPIRIN}


@pytest.mark.parametrize("q,r,expect_edge", [
    (ASPIRIN_ME, ASPIRIN, True),
    (CYCLOHEXANE, ASPIRIN, False),
])
def test_build_is_bipartite_query_to_reference(q, r, expect_edge):
    """Only (path-less -> annotated) pairs are computed. All-pairs over 12,821
    compounds is 82M comparisons, most producing edges between compounds that
    already have paths -- degree without reach."""
    em = Emit()
    build(em, {"q": q, "r": r}, {"q"}, {"r"}, "s", threshold=0.5)
    subjects = {e["subject"] for e in em.edges}
    assert ("q" in subjects) is expect_edge
    assert "r" not in subjects
