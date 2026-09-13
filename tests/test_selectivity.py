"""Gate for the selectivity layer: the ratio is only computed where it means
something, and the worst measurement wins.
"""
import json
import sqlite3

import pytest

from kgav.emit import Emit
from kgav.schema import load_schema
from kgav.selectivity import (
    SI_THRESHOLD,
    TOX_TYPES,
    ingest_selectivity,
    query_pairs,
    selectivity_index,
    selectivity_index_by_compound,
)

V = "NCBITaxon:2697049"
KEY_A = "A" * 14 + "-" + "B" * 10 + "-C"
KEY_B = "D" * 14 + "-" + "E" * 10 + "-F"


def test_selectivity_index_arithmetic():
    assert selectivity_index(100.0, 5000.0) == 50.0
    assert selectivity_index(0, 5000.0) is None
    assert selectivity_index(100.0, 0) is None
    assert selectivity_index(None, 5000.0) is None


def test_gi50_is_not_treated_as_cytotoxicity():
    """GI50 measures growth arrest, not death, and is dominated by oncology
    panels. Mixing it with CC50 would conflate two endpoints."""
    assert "GI50" not in TOX_TYPES
    assert set(TOX_TYPES) == {"CC50", "TC50"}


def _row(inchikey, ec50, cc50, tax=2697049, doc=1, year=2021):
    return {"molregno": 1, "doc_id": doc, "tax_id": tax, "ec50_nm": ec50,
            "n_activity": 1, "cc50_nm": cc50, "n_tox": 1, "compound": "CHEMBL1",
            "inchikey": inchikey, "year": year, "pubmed_id": 33333333}


@pytest.fixture
def tmap():
    return {"2697049": "2697049"}


def test_selective_and_cytotoxic_are_counted(tmap):
    em = Emit()
    stats = ingest_selectivity(em, [_row(KEY_A, 100.0, 5000.0),
                                    _row(KEY_B, 100.0, 200.0)],
                               tmap, "infores:chembl")
    assert stats["selective"] == 1     # SI 50
    assert stats["cytotoxic"] == 1     # SI 2
    assert stats["edges"] == 2


def test_worst_selectivity_wins_across_documents(tmap):
    """A compound shown cytotoxic anywhere is cytotoxic. Taking the flattering
    measurement is how a selectivity filter stops filtering."""
    em = Emit()
    ingest_selectivity(em, [_row(KEY_A, 100.0, 9000.0, doc=1),   # SI 90
                            _row(KEY_A, 100.0, 300.0, doc=2)],   # SI 3
                       tmap, "infores:chembl")
    assert len(em.edges) == 1
    assert em.edges[0]["qualifiers"]["selectivity_index"] == 3.0


def test_edges_are_marked_verified_and_quantified(tmap):
    em = Emit()
    ingest_selectivity(em, [_row(KEY_A, 100.0, 5000.0)], tmap, "infores:chembl")
    q = em.edges[0]["qualifiers"]
    assert q["selectivity_verified"] is True
    assert q["unquantified"] is False
    assert q["cc50_nm"] == 5000.0


def test_rows_without_structure_or_taxon_are_dropped(tmap):
    em = Emit()
    stats = ingest_selectivity(em, [_row(None, 100.0, 5000.0),
                                    _row(KEY_A, 100.0, 5000.0, tax=99999)],
                               tmap, "infores:chembl")
    assert stats["no_inchikey"] == 1 and stats["unmapped_taxon"] == 1
    assert stats["edges"] == 0


def test_threshold_is_configurable(tmap):
    em = Emit()
    stats = ingest_selectivity(em, [_row(KEY_A, 100.0, 500.0)], tmap,
                               "infores:chembl", threshold=3.0)
    assert stats["selective"] == 1     # SI 5 clears a threshold of 3
    em2 = Emit()
    stats2 = ingest_selectivity(em2, [_row(KEY_A, 100.0, 500.0)], tmap,
                                "infores:chembl", threshold=SI_THRESHOLD)
    assert stats2["cytotoxic"] == 1    # but not the default of 10


def test_output_validates_against_schema(tmap):
    em = Emit()
    ingest_selectivity(em, [_row(KEY_A, 100.0, 5000.0)], tmap, "infores:chembl")
    nodes = [
        {"id": f"INCHIKEY:{KEY_A}", "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": KEY_A[:14],
                        "is_approved": True, "salt_collapsed": False,
                        "stereo_collapsed": False}},
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
    ]
    violations = load_schema().validate_batch(nodes, em.edges)
    assert violations == [], [str(v) for v in violations[:5]]


def test_index_lookup_reads_a_release(tmp_path):
    e = {"subject": "INCHIKEY:X", "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
         "object": V, "qualifiers": {"selectivity_index": 42.0, "assay_type": "binding"},
         "primary_knowledge_source": "s", "evidence_tier": 1,
         "first_asserted_date": "2021-01-01"}
    (tmp_path / "edges.jsonl").write_text(json.dumps(e))
    assert selectivity_index_by_compound(tmp_path) == {("INCHIKEY:X", V): 42.0}


def test_query_pairs_requires_same_document(tmp_path):
    """Cross-document pairing yields more pairs and a meaningless ratio: a
    CC50 from one cell system and an EC50 from another."""
    db = tmp_path / "c.db"
    con = sqlite3.connect(db)
    con.executescript("""
    CREATE TABLE molecule_dictionary(molregno INT, chembl_id TEXT);
    CREATE TABLE compound_structures(molregno INT, standard_inchi_key TEXT);
    CREATE TABLE target_dictionary(tid INT, tax_id INT);
    CREATE TABLE assays(assay_id INT, tid INT, doc_id INT);
    CREATE TABLE docs(doc_id INT, year INT, pubmed_id INT);
    CREATE TABLE activities(activity_id INT, assay_id INT, molregno INT,
      standard_type TEXT, standard_value REAL, standard_units TEXT,
      standard_relation TEXT);
    """)
    con.executemany("INSERT INTO molecule_dictionary VALUES(?,?)", [(1, "CHEMBL1")])
    con.executemany("INSERT INTO compound_structures VALUES(?,?)", [(1, KEY_A)])
    con.executemany("INSERT INTO target_dictionary VALUES(?,?)", [(10, 2697049), (11, 9606)])
    con.executemany("INSERT INTO assays VALUES(?,?,?)",
                    [(1, 10, 100), (2, 11, 100), (3, 11, 200)])
    con.executemany("INSERT INTO docs VALUES(?,?,?)", [(100, 2021, 1), (200, 2022, 2)])
    con.executemany("INSERT INTO activities VALUES(?,?,?,?,?,?,?)",
                    [(1, 1, 1, "EC50", 100.0, "nM", "="),
                     (2, 2, 1, "CC50", 5000.0, "nM", "="),    # same doc -> pairs
                     (3, 3, 1, "CC50", 20.0, "nM", "=")])     # other doc -> ignored
    con.commit()
    con.close()

    rows = query_pairs(db, ["2697049"])
    assert len(rows) == 1
    assert rows[0]["cc50_nm"] == 5000.0
