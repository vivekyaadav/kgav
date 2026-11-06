"""Gate for the drug -> human protein layer: direction is never invented."""
import inspect
import sqlite3
from pathlib import Path

import pytest

from kgav.aliases import AliasIndex
from kgav.emit import Emit
from kgav.host_targets import (
    ACTION_TO_DIRECTION,
    host_protein_ids,
    ingest_host_targets,
    load_mechanisms,
    query_host_activities,
    to_nm,
)
from kgav.schema import load_schema

ROOT = Path(__file__).resolve().parents[1]
CHEMBL_DB = ROOT / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db"
needs_db = pytest.mark.skipif(not CHEMBL_DB.exists(), reason="ChEMBL not downloaded")

ACE2 = "UniProtKB:Q9BYF1"


@pytest.fixture
def db(tmp_path):
    p = tmp_path / "c.db"
    con = sqlite3.connect(p)
    con.executescript("""
    CREATE TABLE molecule_dictionary(molregno INT, chembl_id TEXT, max_phase INT);
    CREATE TABLE compound_structures(molregno INT, canonical_smiles TEXT, standard_inchi_key TEXT);
    CREATE TABLE target_dictionary(tid INT, chembl_id TEXT, target_type TEXT, tax_id INT);
    CREATE TABLE target_components(tid INT, component_id INT);
    CREATE TABLE component_sequences(component_id INT, accession TEXT);
    CREATE TABLE assays(assay_id INT, tid INT, doc_id INT);
    CREATE TABLE docs(doc_id INT, year INT, pubmed_id INT);
    CREATE TABLE activities(activity_id INT, assay_id INT, molregno INT, standard_type TEXT,
      standard_value REAL, standard_units TEXT, standard_relation TEXT, pchembl_value REAL);
    CREATE TABLE drug_mechanism(molregno INT, tid INT, action_type TEXT);
    """)
    con.executemany("INSERT INTO molecule_dictionary VALUES(?,?,?)",
                    [(1, "CHEMBL1", 4), (2, "CHEMBL2", 2), (3, "CHEMBL3", 4), (4, "CHEMBL4", 1)])
    con.executemany("INSERT INTO compound_structures VALUES(?,?,?)",
                    [(1, "CC1", "A" * 14 + "-" + "B" * 10 + "-C"),
                     (2, "CC2", "D" * 14 + "-" + "E" * 10 + "-F"),
                     (3, "CC3", None),
                     (4, "CC4", "G" * 14 + "-" + "H" * 10 + "-I")])
    con.executemany("INSERT INTO target_dictionary VALUES(?,?,?,?)",
                    [(10, "T1", "SINGLE PROTEIN", 9606), (11, "T2", "SINGLE PROTEIN", 9606)])
    con.executemany("INSERT INTO target_components VALUES(?,?)", [(10, 100), (11, 101)])
    con.executemany("INSERT INTO component_sequences VALUES(?,?)",
                    [(100, "Q9BYF1"), (101, "NOTINGRAPH")])
    # Assays 3 and 4 hit the SAME target as assay 1 from DIFFERENT documents,
    # so activities 8 and 9 collapse into the edge activity 1 creates. Doc 1
    # carries no pubmed_id on purpose: an edge whose FIRST row is uncitable
    # must still gain the citations of the rows that collapse into it.
    con.executemany("INSERT INTO assays VALUES(?,?,?)",
                    [(1, 10, 1), (2, 11, 1), (3, 10, 2), (4, 10, 3)])
    con.executemany("INSERT INTO docs VALUES(?,?,?)",
                    [(1, 2019, None), (2, 2020, 22222222),
                     (3, 2021, 33333333)])
    con.executemany("INSERT INTO activities VALUES(?,?,?,?,?,?,?,?)",
                    [(1, 1, 1, "IC50", 50.0, "nM", "=", 7.3),
                     (2, 1, 2, "Ki", 0.5, "uM", "=", 6.3),
                     (3, 1, 3, "IC50", 10.0, "nM", "=", 8.0),
                     (4, 2, 1, "IC50", 20.0, "nM", "=", 7.7),
                     (5, 1, 1, "IC50", 1.0, "ug.mL-1", "=", 6.5),
                     (6, 1, 4, "IC50", 5.0, "nM", "=", 8.3),
                     (7, 1, 1, "EC50", 9.0, "nM", "=", 5.1),
                     (8, 3, 1, "IC50", 60.0, "nM", "=", 7.2),
                     (9, 4, 1, "IC50", 70.0, "nM", "=", 7.1)])
    con.executemany("INSERT INTO drug_mechanism VALUES(?,?,?)",
                    [(1, 10, "INHIBITOR"), (2, 10, "AGONIST")])
    con.commit()
    con.close()
    return p


@pytest.fixture
def host(tmp_path):
    import json
    d = tmp_path / "host"
    d.mkdir()
    n = [{"id": ACE2, "class": "Protein",
          "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False,
                         "sequence_hash": "h", "reviewed": True, "gene_symbol": "ACE2"}}]
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(x) for x in n))
    (d / "edges.jsonl").write_text("")
    return d


@pytest.fixture
def built(db, host):
    index = AliasIndex.from_releases(host)
    em = Emit()
    stats = ingest_host_targets(em, query_host_activities(db), load_mechanisms(db),
                                host_protein_ids(index), "infores:chembl")
    return em, stats, index


# --------------------------------------------------------------------- units
@pytest.mark.parametrize("v,u,exp", [(50.0, "nM", 50.0), (0.5, "uM", 500.0),
                                     (1000.0, "pM", 1.0), (1.0, "ug.mL-1", None),
                                     (1.0, None, None)])
def test_unit_conversion(v, u, exp):
    assert to_nm(v, u) == exp


# ----------------------------------------------------------------- direction
def test_direction_comes_from_drug_mechanism_when_stated(built):
    em, _stats, _ = built
    dirs = {e["qualifiers"]["direction"] for e in em.edges}
    assert dirs == {"inhibitor", "agonist"}


def test_direction_defaults_to_unknown_not_inhibitor(db, host):
    """Inhibiting a dependency factor should block infection; ACTIVATING one
    should not. ChEMBL's activity table records potency, not mechanism, so
    asserting 'inhibitor' would manufacture a claim the data does not make and
    every M2 path would silently inherit it."""
    index = AliasIndex.from_releases(host)
    em = Emit()
    stats = ingest_host_targets(em, query_host_activities(db), {},  # no mechanisms
                                host_protein_ids(index), "infores:chembl")
    assert stats["direction_unknown"] == stats["edges"] > 0
    assert all(e["qualifiers"]["direction"] == "unknown" for e in em.edges)


def test_antagonist_and_blocker_map_to_inhibitor():
    for a in ("ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR"):
        assert ACTION_TO_DIRECTION[a] == "inhibitor"
    for a in ("AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR"):
        assert ACTION_TO_DIRECTION[a] == "agonist"


# ------------------------------------------------------------------ filtering
def test_targets_absent_from_the_graph_are_dropped(built):
    """Minting a protein node here would give it no proteome data, and it
    contributes to no metapath."""
    _em, stats, _ = built
    assert stats["target_not_in_graph"] == 1


def test_low_pchembl_and_low_phase_are_excluded_by_the_query(db):
    """pChEMBL 5.1 and max_phase 1 rows must never reach the ingest."""
    rows = query_host_activities(db)
    assert all(r["pchembl_value"] >= 6.0 for r in rows)
    assert all(r["max_phase"] >= 2 for r in rows)


def test_rows_without_an_inchikey_are_dropped(built):
    _em, stats, _ = built
    assert stats["no_inchikey"] == 1


def test_unconvertible_units_are_dropped(built):
    _em, stats, _ = built
    assert stats["unconvertible_units"] == 1


# ------------------------------------------------------------------ integrity
def test_compound_properties_match_the_day7_layer(built):
    """A mismatch here produces property conflicts on thousands of shared
    compounds at assembly time."""
    em, _, _ = built
    for n in em.nodes.values():
        p = n["properties"]
        assert p["salt_collapsed"] is False and p["stereo_collapsed"] is False
        assert p["inchikey_skel"] == n["id"].split(":", 1)[1][:14]


def test_edges_carry_a_quantitative_value_and_relation(built):
    em, _, _ = built
    for e in em.edges:
        q = e["qualifiers"]
        assert q["relation"] in {"=", ">", "<", ">=", "<=", "~"}
        assert any(k in q for k in ("ic50_nm", "ec50_nm", "ki_nm", "kd_nm"))


def test_output_validates_against_schema(built):
    em, _, index = built
    nodes = list(index.nodes.values()) + list(em.nodes.values())
    violations = load_schema().validate_batch(nodes, em.edges)
    assert violations == [], [str(v) for v in violations[:6]]


# ----------------------------------------------------------------- real data
@needs_db
def test_real_scale_and_m2_reachability():
    rows = query_host_activities(CHEMBL_DB)
    assert len(rows) > 40_000, f"{len(rows)} activities"
    drugs = {r["compound"] for r in rows}
    assert len(drugs) > 2_000, f"{len(drugs)} drugs"
    targets = {r["accession"] for r in rows if r["accession"]}
    assert len(targets) > 1_000, f"{len(targets)} targets"


# ------------------------------------------------------------- publications
# ZERO OF 18,954 TARGETS EDGES IN v0.1 CARRY A PUBLICATION, the only
# predicate in the graph with no coverage at all, because the query selected
# docs.year and never docs.pubmed_id -- the same row, one column over. An
# agent asked for evidence on a host-directed route could show a potency
# with no source for it. Nothing here failed: no test looked, and the
# driver printed no coverage line to look at.
def test_edges_cite_the_document_their_date_came_from(built):
    em, _stats, _ = built
    cited = [e for e in em.edges if e.get("publications")]
    assert cited, "no TARGETS edge carries a publication"
    for e in cited:
        assert all(p.startswith("PMID:") for p in e["publications"])


def test_a_collapsed_duplicate_contributes_its_paper(built):
    """The key is (compound, target, assay type) and ~3 source rows share
    one. The edge keeps the first row's value; the others are not evidence
    that can be thrown away, because each is a separate measurement in a
    separate paper."""
    em, stats, _ = built
    assert stats["duplicate"] == 2, "the fixture's two extra IC50 rows"
    assert stats["duplicate_added_publication"] == 2
    e = next(e for e in em.edges
             if e["predicate"] == "TARGETS"
             and "ic50_nm" in e["qualifiers"])
    # Document 1, which activity 1 came from, has no pubmed_id. Both of this
    # edge's citations arrive from rows that collapsed into it -- so an edge
    # whose own row is uncitable is not therefore uncitable.
    assert e["publications"] == ["PMID:22222222", "PMID:33333333"]


def test_an_uncitable_row_yields_an_edge_with_no_publications_key(built):
    """Absent, not an empty list. `publications: []` reads as "searched and
    found none"; this project has shipped that confusion before."""
    em, stats, _ = built
    e = next(e for e in em.edges if "ki_nm" in e["qualifiers"])
    assert "publications" not in e
    assert stats["no_publication"] == 1 and stats["with_publication"] == 1


def test_the_counters_add_up_to_the_edges(built):
    """with_publication + no_publication is a partition of the edges, and
    the direction census is too. Both counted ROWS until 2026-11-07, so
    each summed to the source rows while sitting beside an edge count."""
    _em, stats, _ = built
    assert stats["with_publication"] + stats["no_publication"] == stats["edges"]
    assert sum(v for k, v in stats.items()
               if k.startswith("direction_")) == stats["edges"]


def test_the_row_order_the_first_row_wins_by_is_declared(db):
    """Which duplicate sets the potency and the date decided the value of
    thousands of edges and was left to the query planner. activity_id is
    ChEMBL's own order; without ORDER BY two rebuilds of one database can
    disagree, and the disagreement moves edges across a temporal split."""
    rows = query_host_activities(db)
    ids = [r["compound"] for r in rows]
    assert ids == [r["compound"] for r in query_host_activities(db)]
    assert "ORDER BY" in inspect.getsource(query_host_activities)
