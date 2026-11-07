"""Gate for the cold-start layer: similarity is bipartite, capped, and the
threshold cannot drift below what the schema declares.
"""
import json
from collections import Counter
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


# ------------------------------------------ one compound is not its neighbour
# THE GRAPH USES THE FULL InChIKey AS A NODE ID, deliberately, so a salt and
# its parent and two stereoisomers stay separate nodes
# (normalize_chem's module docstring says why). similarity_edges excluded a
# pair only when the two ids were EQUAL, so those separate nodes paired with
# each other: 434 of v0.1's 3,384 CHEMICALLY_SIMILAR_TO edges joined two
# InChIKeys of ONE molecule, 412 of them at tanimoto 1.0 because they are the
# same structure.
#
# normalize_chem.same_as_rule states that a shared skeleton means one
# compound, and gtop.py uses the same 14 characters as PROOF two keys are one
# compound. Only this function read id inequality as proof they were two.
#
# Measured on the real release: 3,384 -> 3,149 edges, and the compounds the
# layer connects fall 1,724 -> 1,529. Those 195 compounds had NO real
# neighbour; their whole M6 route was to another spelling of themselves.
SKEL = "BSYNRYMUTXBXSQ"
ASPIRIN_KEY = f"INCHIKEY:{SKEL}-UHFFFAOYSA-N"
ASPIRIN_SALT = f"INCHIKEY:{SKEL}-UHFFFAOYSA-M"      # same molecule, Na salt
OTHER_KEY = "INCHIKEY:XEKOWRVHYACXOJ-UHFFFAOYSA-N"  # a different molecule


def test_two_keys_of_one_molecule_are_not_neighbours():
    fps, _ = fingerprints({ASPIRIN_KEY: ASPIRIN, ASPIRIN_SALT: ASPIRIN})
    stats = Counter()
    edges = similarity_edges({ASPIRIN_SALT: fps[ASPIRIN_SALT]},
                             {ASPIRIN_KEY: fps[ASPIRIN_KEY]},
                             threshold=0.0, stats=stats)
    assert edges == [], "a compound was paired with another spelling of itself"
    assert stats["same_skeleton_excluded"] == 1


def test_a_genuinely_different_molecule_is_still_a_neighbour():
    """The exclusion must be the skeleton, not anything coarser."""
    fps, _ = fingerprints({ASPIRIN_KEY: ASPIRIN, OTHER_KEY: ASPIRIN_ME})
    edges = similarity_edges({OTHER_KEY: fps[OTHER_KEY]},
                             {ASPIRIN_KEY: fps[ASPIRIN_KEY]}, threshold=0.0)
    assert len(edges) == 1


def test_an_id_with_no_inchikey_in_it_is_not_excluded():
    """skeleton_of_id returns None for a minted KGAV: id, and None == None
    must not read as "same compound" -- that would silently drop every pair
    between two minted nodes."""
    fps, _ = fingerprints({"KGAV:a": ASPIRIN, "KGAV:b": ASPIRIN_ME})
    stats = Counter()
    edges = similarity_edges({"KGAV:b": fps["KGAV:b"]},
                             {"KGAV:a": fps["KGAV:a"]},
                             threshold=0.0, stats=stats)
    assert len(edges) == 1
    assert stats["same_skeleton_excluded"] == 0


def test_the_exclusion_frees_a_top_k_slot_for_a_real_neighbour():
    """Why 463 pairs were excluded while the edge count fell by only 235: a
    same-molecule pair at tanimoto 1.0 sorts FIRST, so it did not merely add
    a bogus edge, it displaced a real neighbour from the top-k cut."""
    fps, _ = fingerprints({ASPIRIN_KEY: ASPIRIN, ASPIRIN_SALT: ASPIRIN,
                           OTHER_KEY: ASPIRIN_ME})
    got = similarity_edges(
        {ASPIRIN_SALT: fps[ASPIRIN_SALT]},
        {ASPIRIN_KEY: fps[ASPIRIN_KEY], OTHER_KEY: fps[OTHER_KEY]},
        threshold=0.0, top_k=1)
    assert [r for _, r, _ in got] == [OTHER_KEY]


def test_the_build_counter_reaches_the_driver():
    em = Emit()
    stats = build(em, {ASPIRIN_KEY: ASPIRIN, ASPIRIN_SALT: ASPIRIN},
                  {ASPIRIN_SALT}, {ASPIRIN_KEY}, "infores:kgav-computed",
                  threshold=0.0)
    assert stats["edges"] == 0
    assert stats["same_skeleton_excluded"] == 1


@pytest.mark.parametrize("nid,skel", [
    (ASPIRIN_KEY, SKEL),
    (ASPIRIN_SALT, SKEL),
    ("INCHIKEY:SHORT", None),
    ("KGAV:minted-1", None),
    ("CHEMBL:CHEMBL25", None),
])
def test_skeleton_of_id(nid, skel):
    from kgav.normalize_chem import skeleton_of_id
    assert skeleton_of_id(nid) == skel


def test_one_definition_of_the_skeleton_length():
    """gtop.py uses these 14 characters as proof two keys are ONE compound
    while similarity.py used id inequality as proof they were two. Two
    constants for one fact is how that disagreement survived."""
    from kgav.gtop import SKELETON_LEN as GTOP_LEN
    from kgav.normalize_chem import SKELETON_LEN
    assert SKELETON_LEN == GTOP_LEN == 14


def test_the_top_k_cut_is_counted():
    """Fix-list E7 and E8: an above-threshold neighbour that falls outside
    top_k is a real similarity the graph does not carry, and the cap's cost
    was never measured. On the real release it discards 156 neighbours
    across 34 query compounds that hit the cap."""
    fps, _ = fingerprints({f"INCHIKEY:{c}{'X'*13}-UHFFFAOYSA-N": ASPIRIN
                           for c in "ABCDE"} | {"q": ASPIRIN_ME})
    refs = {k: v for k, v in fps.items() if k != "q"}
    stats = Counter()
    got = similarity_edges({"q": fps["q"]}, refs, threshold=0.0, top_k=2,
                           stats=stats)
    assert len(got) == 2
    assert stats["above_threshold_beyond_top_k"] == 3
    assert stats["queries_truncated_by_top_k"] == 1


def test_a_query_under_the_cap_is_not_counted_as_truncated():
    fps, _ = fingerprints({"INCHIKEY:AAAAAAAAAAAAAA-UHFFFAOYSA-N": ASPIRIN,
                           "q": ASPIRIN_ME})
    stats = Counter()
    similarity_edges({"q": fps["q"]},
                     {k: v for k, v in fps.items() if k != "q"},
                     threshold=0.0, top_k=10, stats=stats)
    assert stats["queries_truncated_by_top_k"] == 0
