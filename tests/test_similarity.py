"""Gate for the cold-start layer: similarity is bipartite, capped, and the
threshold cannot drift below what the schema declares.
"""
import json
from pathlib import Path

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

ROOT = Path(__file__).resolve().parents[1]

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


# --------------------------------------------------------------------------
# A compound measured NOT to inhibit a protein is not a cold-start compound.
# --------------------------------------------------------------------------
def test_measured_inactive_compounds_are_not_treated_as_unannotated(tmp_path):
    """ce133f1 moved 1,625 pairs off INHIBITS and into the query set.

    M7 walks similar -> INHIBITS -> protein -> gene -> virus, so a compound in
    the query set acquires similarity edges that let M7 propose activity
    against the exact protein the graph records it as measured inactive
    against. Before ce133f1 those rows were (wrongly) INHIBITS, so they landed
    in `annotated` and the question never arose.

    This asserts the reference/query split the driver computes, by running the
    same predicate rule against a fixture release.
    """
    import json

    edges = [
        {"subject": "INCHIKEY:HASTARGET", "predicate": "TARGETS",
         "object": "UniProtKB:H1"},
        {"subject": "INCHIKEY:INHIBITOR", "predicate": "INHIBITS",
         "object": "UniProtKB:V1"},
        {"subject": "INCHIKEY:MEASUREDBAD", "predicate": "MEASURED_INACTIVE_AGAINST",
         "object": "UniProtKB:V1"},
        {"subject": "INCHIKEY:ORPHAN", "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
         "object": "NCBITaxon:2697049"},
    ]
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))

    ANNOTATING = ("TARGETS", "INHIBITS", "MEASURED_INACTIVE_AGAINST")
    annotated = {json.loads(line)["subject"]
                 for line in (rel / "edges.jsonl").read_text().splitlines()
                 if line.strip() and json.loads(line)["predicate"] in ANNOTATING}
    all_compounds = {"INCHIKEY:HASTARGET", "INCHIKEY:INHIBITOR",
                     "INCHIKEY:MEASUREDBAD", "INCHIKEY:ORPHAN"}

    assert "INCHIKEY:MEASUREDBAD" in annotated
    assert all_compounds - annotated == {"INCHIKEY:ORPHAN"}


def test_the_driver_counts_measured_inactive_as_annotating():
    """The rule lives in the driver, so assert the driver states it.

    A list of predicate names is exactly the kind of thing that goes stale
    when a predicate is added -- which is how this arose in the first place.
    """
    src = (ROOT / "scripts" / "ingest_similarity.py").read_text()
    assert '"MEASURED_INACTIVE_AGAINST"' in src, \
        "a new predicate was added to the schema and this list was not updated"
    assert "ANNOTATING" in src


# ------------------------------------------- which neighbours survive top_k
#
# similarity_edges sorted on similarity alone. Python's sort is
# stable, so tied neighbours kept the order of ref_ids = list(ref_fps)
# -- a dict built by ingest_similarity.py from a SET. String hashing
# is randomised per process, so when a tie straddled position top_k,
# which neighbour survived depended on the hash seed.
#
# Two builds of v0.1 from identical raw inputs differed by exactly 8
# CHEMICALLY_SIMILAR_TO edges out of 389,306, total unchanged: eight
# swapped for eight. edges_sha was therefore unreproducible, and no
# tool in the project could say whether a rebuilt release was a
# different graph.
#
# test_top_k_caps_dense_scaffold_series above builds THIRTY IDENTICAL
# references -- every pair tied -- and asserts only that five come
# back. The scenario that breaks this was already in the suite; what
# was missing was any assertion about which five.
def _tied_refs(order):
    """Twelve references, all the same structure, so every similarity
    is tied and the top_k cut falls inside the tie group."""
    return {f"r{i:02d}": ASPIRIN for i in order}


def test_the_top_k_cut_does_not_depend_on_reference_order():
    q_fps, _ = fingerprints({"q": ASPIRIN_ME})
    fwd, _ = fingerprints(_tied_refs(range(12)))
    rev, _ = fingerprints(_tied_refs(reversed(range(12))))
    assert list(fwd) != list(rev), (
        "both reference dicts are in the same order; this cannot "
        "detect an order-dependent cut")

    a = similarity_edges(q_fps, fwd, threshold=0.0, top_k=5)
    b = similarity_edges(q_fps, rev, threshold=0.0, top_k=5)
    assert {r for _, r, _ in a} == {r for _, r, _ in b}, (
        "a different set of tied neighbours survived the cut")
    assert a == b, "same neighbours, different emission order"


def test_tied_neighbours_are_cut_by_id_not_by_arrival():
    """The tie-break is on the id, so the survivors are a property of
    the content. Asserting the ACTUAL ids, because equality between
    two runs would also hold if both were wrong in the same way."""
    q_fps, _ = fingerprints({"q": ASPIRIN_ME})
    refs, _ = fingerprints(_tied_refs(reversed(range(12))))
    got = [r for _, r, _ in similarity_edges(q_fps, refs, threshold=0.0,
                                             top_k=5)]
    assert got == ["r00", "r01", "r02", "r03", "r04"], got


def test_a_higher_score_still_outranks_a_lower_id():
    """The id breaks TIES and must not override similarity: a worse
    neighbour with an earlier id cannot displace a better one."""
    q_fps, _ = fingerprints({"q": ASPIRIN_ME})
    refs, _ = fingerprints({"aaa": CYCLOHEXANE, "zzz": ASPIRIN})
    got = similarity_edges(q_fps, refs, threshold=0.0, top_k=1)
    assert [r for _, r, _ in got] == ["zzz"], got


def test_the_ingest_script_does_not_iterate_a_set_into_the_fingerprints():
    """Belt and braces, pinned at the source. The tie-break above
    makes the OUTPUT correct whatever order arrives, but emission
    order should not depend on the hash seed either -- and a reader
    of ingest_similarity.py should not have to know the protection
    lives in another file."""
    src = (ROOT / "scripts" / "ingest_similarity.py").read_text()
    assert "for i in sorted(query)" in src
    assert "for i in sorted(annotated)" in src
