"""Gate for the selectivity layer: the ratio is only computed where it means
something, and the worst measurement wins.
"""
from collections import Counter
from pathlib import Path
import json
import sqlite3

import pytest

from kgav.emit import Emit
from kgav.schema import load_schema
from kgav.selectivity import (
    SI_THRESHOLD,
    TOX_TYPES,
    extract_cell_line,
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


def _row(inchikey, ec50, cc50, tax=2697049, doc=1, year=2021, exact=True,
         tox_text="cytotoxicity in Vero E6 cells"):
    """cc50_exact is a bool: True when the retained CC50 is an exact
    measurement, False when it is a censored lower bound."""
    return {"molregno": 1, "doc_id": doc, "tax_id": tax, "ec50_nm": ec50,
            "n_activity": 1, "cc50_nm": cc50, "n_tox": 1, "compound": "CHEMBL1",
            "inchikey": inchikey, "year": year, "pubmed_id": 33333333,
            "cc50_exact": exact,
            "assay_text": "assay in Vero E6 cells",
            # THE COUNTERSCREEN'S OWN DESCRIPTION. selectivity_verified means
            # EC50 and CC50 come from one experimental context, so the two
            # cell lines have to be comparable; a fixture carrying only one
            # of them cannot express the distinction at all.
            "tox_descr": tox_text}


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


def test_median_of_exact_measurements_is_used(tmap):
    """Neither the worst nor the best measurement: both select a study by its
    answer. Taking the worst of 34 remdesivir measurements picked an outlier
    and reported a licensed antiviral at SI 0.9; taking the best exonerated
    chloroquine at 17.7. The median of exact measurements is the summary that
    survived the named controls."""
    em = Emit()
    ingest_selectivity(em, [_row(KEY_A, 100.0, 9000.0, doc=1),   # SI 90
                            _row(KEY_A, 100.0, 300.0, doc=2),    # SI 3
                            _row(KEY_A, 100.0, 1000.0, doc=3)],  # SI 10
                       tmap, "infores:chembl")
    assert len(em.edges) == 1
    assert em.edges[0]["qualifiers"]["selectivity_index"] == 10.0


def test_censored_measurement_below_threshold_is_uninformative(tmap):
    """"CC50 > 300 nM" means NOT toxic up to 300 nM -- a lower bound on
    selectivity, never evidence of cytotoxicity."""
    em = Emit()
    stats = ingest_selectivity(em, [_row(KEY_A, 100.0, 300.0, exact=False)],
                               tmap, "infores:chembl")
    assert stats["censored_uninformative"] == 1
    assert stats["edges"] == 0


def test_censored_measurement_above_threshold_proves_selectivity(tmap):
    """"CC50 > 50000 nM" against EC50 100 nM establishes SI of at least 500."""
    em = Emit()
    stats = ingest_selectivity(em, [_row(KEY_A, 100.0, 50000.0, exact=False)],
                               tmap, "infores:chembl")
    assert stats["selective"] == 1


def test_row_without_a_relation_flag_is_dropped(tmap):
    """Defaulting a missing direction-bearing field is the error that produced
    SI 1.8 for nirmatrelvir."""
    row = _row(KEY_A, 100.0, 5000.0)
    del row["cc50_exact"]
    em = Emit()
    stats = ingest_selectivity(em, [row], tmap, "infores:chembl")
    assert stats["missing_relation_flag"] == 1 and stats["edges"] == 0


def test_edges_carry_the_evidence_not_a_quantified_flag(tmap):
    """This layer used to write unquantified=false to mean "a selectivity
    index exists" -- a third definition of the field, which then lost the
    merge to chembl's true. It states the measurements and nothing else; the
    flag is derived from them at assembly."""
    em = Emit()
    ingest_selectivity(em, [_row(KEY_A, 100.0, 5000.0)], tmap, "infores:chembl")
    q = em.edges[0]["qualifiers"]
    assert q["selectivity_verified"] is True
    assert "unquantified" not in q
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
    CREATE TABLE assays(assay_id INT, tid INT, doc_id INT, description TEXT);
    CREATE TABLE docs(doc_id INT, year INT, pubmed_id INT);
    CREATE TABLE activities(activity_id INT, assay_id INT, molregno INT,
      standard_type TEXT, standard_value REAL, standard_units TEXT,
      standard_relation TEXT);
    """)
    con.executemany("INSERT INTO molecule_dictionary VALUES(?,?)", [(1, "CHEMBL1")])
    con.executemany("INSERT INTO compound_structures VALUES(?,?)", [(1, KEY_A)])
    con.executemany("INSERT INTO target_dictionary VALUES(?,?)", [(10, 2697049), (11, 9606)])
    con.executemany("INSERT INTO assays VALUES(?,?,?,?)",
                    [(1, 10, 100, "Antiviral activity in Vero E6 cells"),
                     (2, 11, 100, "Cytotoxicity against Vero E6 cells"),
                     (3, 11, 200, "Cytotoxicity against HeLa cells")])
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


@pytest.mark.parametrize("text,expected", [
    ("Cytotoxicity against African green monkey Vero E6 cells", "Vero E6"),
    ("Antiviral activity in Vero cells", "Vero"),
    ("Antiviral activity against SARS-CoV-2 in Calu-3 cells", "Calu-3"),
    ("Inhibition in human Huh-7.5 hepatoma cells", "Huh-7.5"),
    ("Inhibition in human Huh-7 cells", "Huh-7"),
    ("Cytotoxicity against HEK293T cells transfected with ACE2", "HEK293T"),
    ("Antiviral activity with no cell line stated", None),
])
def test_cell_line_extraction_prefers_the_specific_line(text, expected):
    """Vero E6 must be matched before Vero, and Huh-7.5 before Huh-7.
    Substring ordering has caused two prior bugs in this project: SARS-CoV
    matching inside SARS-CoV-2, and nsp1 patterns catching nsp1x."""
    from kgav.selectivity import extract_cell_line
    assert extract_cell_line(text) == expected


def test_vero_lines_are_the_dominant_evidence_base():
    """Not a code test: a recorded measurement. 668 of 1,144 selectivity edges
    with a named cell line come from Vero lines, against 54 from Calu-3 and 43
    from A549 -- roughly twelve times more evidence from the system that
    misled the field in 2020 than from the systems where results transfer."""
    from kgav.selectivity import extract_cell_line
    assert extract_cell_line("Vero E6") == "Vero E6"
    assert extract_cell_line("Calu-3") == "Calu-3"


def test_no_comment_in_this_module_still_claims_the_worst_is_kept():
    """A COMMENT THAT OUTLIVED ITS CODE DOES ITS DAMAGE BY BEING QUOTED.

    This module's selection comment ended "among several exact measurements
    the WORST is kept" long after the emit loop had moved to the median --
    with "collect all; median chosen below" on the very next line and
    test_median_of_exact_measurements_is_used pinning the behaviour.

    It was not harmless. baselines.potency_tiebreak took "keep the worst"
    from it as an inherited rule, citing it verbatim, and fix 3 of
    docs/prereg_m1_potency.md was designed on that citation.
    """
    src = (Path(__file__).resolve().parents[1]
           / "src/kgav/selectivity.py").read_text()
    # The phrase may appear only where it is explicitly marked as superseded.
    for i, line in enumerate(src.splitlines()):
        if "WORST is kept" in line:
            window = "\n".join(src.splitlines()[max(0, i - 8):i + 2])
            assert "USED TO END" in window or "superseded" in window, \
                f"line {i + 1} claims the worst is kept without marking it stale"


# ------------------------------------------- one context, not one document
# selectivity_verified was written True on EVERY edge, and the schema
# declares it means "CC50 and EC50 come from the SAME document, so the
# selectivity index reflects one experimental context". The join guarantees
# the document; a document is not a context. A paper routinely measures
# antiviral activity in one cell line and cytotoxicity in another, and the
# ratio is then a potency in Calu-3 over a toxicity in Vero.
#
# Measured on chembl_37 over 4,167 (compound, document, taxon) pairs:
# 2,244 agree, 337 name DIFFERENT cell lines, 1,586 have a side unnamed.
def test_matching_cell_lines_are_verified(tmap):
    em = Emit()
    stats = ingest_selectivity(
        em, [_row(KEY_A, 100.0, 5000.0,
                  tox_text="cytotoxicity in Vero E6 cells")],
        tmap, "infores:chembl")
    assert em.edges[0]["qualifiers"]["selectivity_verified"] is True
    assert stats["context_matched"] == 1


def test_a_cell_line_mismatch_is_not_verified(tmap):
    """The 337. Potency in one cell line over toxicity in another is not a
    selectivity index, and the flag must not claim it is."""
    em = Emit()
    stats = ingest_selectivity(
        em, [_row(KEY_A, 100.0, 5000.0,
                  tox_text="cytotoxicity in Calu-3 cells")],
        tmap, "infores:chembl")
    q = em.edges[0]["qualifiers"]
    assert q["selectivity_verified"] is False
    assert stats["context_cell_line_mismatch"] == 1
    # The edge is still emitted with its measurements -- the ratio is real,
    # it is the CONTEXT claim that is withdrawn.
    assert q["selectivity_index"] == 50.0
    assert q["cell_line"] == "Vero E6", "cell_line is the EC50's, not the CC50's"


def test_an_unnamed_cell_line_cannot_be_verified(tmap):
    """The 1,586. Not an error in the source -- the limit of what a
    document-level join can establish. False, not True by default."""
    em = Emit()
    stats = ingest_selectivity(
        em, [_row(KEY_A, 100.0, 5000.0, tox_text="cytotoxicity assay")],
        tmap, "infores:chembl")
    assert em.edges[0]["qualifiers"]["selectivity_verified"] is False
    assert stats["context_cell_line_unnamed"] == 1


def _chembl(tmp_path, assays, activities, docs=((100, 2021, 1),)):
    """A minimal ChEMBL with whatever rows a test needs."""
    db = tmp_path / "c.db"
    con = sqlite3.connect(db)
    con.executescript("""
    CREATE TABLE molecule_dictionary(molregno INT, chembl_id TEXT);
    CREATE TABLE compound_structures(molregno INT, standard_inchi_key TEXT);
    CREATE TABLE target_dictionary(tid INT, tax_id INT);
    CREATE TABLE assays(assay_id INT, tid INT, doc_id INT, description TEXT);
    CREATE TABLE docs(doc_id INT, year INT, pubmed_id INT);
    CREATE TABLE activities(activity_id INT, assay_id INT, molregno INT,
      standard_type TEXT, standard_value REAL, standard_units TEXT,
      standard_relation TEXT);
    """)
    con.execute("INSERT INTO molecule_dictionary VALUES(1, 'CHEMBL1')")
    con.execute("INSERT INTO compound_structures VALUES(1, ?)", (KEY_A,))
    con.executemany("INSERT INTO target_dictionary VALUES(?,?)",
                    [(10, 2697049), (11, 9606)])
    con.executemany("INSERT INTO assays VALUES(?,?,?,?)", assays)
    con.executemany("INSERT INTO docs VALUES(?,?,?)", docs)
    con.executemany("INSERT INTO activities VALUES(?,?,?,?,?,?,?)", activities)
    con.commit()
    con.close()
    return db


def test_the_cell_line_is_the_potency_s_own_row(tmp_path):
    """L11. MIN(description) sat beside MIN(standard_value) as two
    INDEPENDENT aggregates, so SQLite minimised each column separately: the
    potency came from whichever row sorted first by VALUE and the
    description from whichever sorted first by TEXT. Measured on chembl_37,
    they are different rows in 1,249 of 2,771 multi-row groups (45.1%), and
    794 of those name a DIFFERENT CELL LINE than the potency.

    The fixture below is that shape exactly: the cheaper EC50 is in Calu-3
    and the dearer one in "A549", which sorts first alphabetically.
    """
    db = _chembl(
        tmp_path,
        assays=[(1, 10, 100, "Antiviral activity in Calu-3 cells"),
                (2, 10, 100, "A549 antiviral activity"),
                (3, 11, 100, "Cytotoxicity against Calu-3 cells")],
        activities=[(1, 1, 1, "EC50", 100.0, "nM", "="),    # the minimum
                    (2, 2, 1, "EC50", 900.0, "nM", "="),    # sorts first by text
                    (3, 3, 1, "CC50", 5000.0, "nM", "=")])
    rows = query_pairs(db, ["2697049"])
    assert len(rows) == 1
    assert rows[0]["ec50_nm"] == 100.0
    assert extract_cell_line(rows[0]["assay_text"]) == "Calu-3", (
        "the cell line came from a row other than the retained potency's")
    assert rows[0]["n_activity"] == 2


def test_a_convertible_unit_is_converted_not_dropped(tmp_path):
    """E7. The SQL demanded units be literally 'nM', so anything else
    vanished with nothing recorded. On the toxicity side that was 4,327 uM
    and 116 mM rows, worth 72 extra (compound, document) pairs."""
    db = _chembl(
        tmp_path,
        assays=[(1, 10, 100, "Antiviral activity in Vero E6 cells"),
                (2, 11, 100, "Cytotoxicity against Vero E6 cells")],
        activities=[(1, 1, 1, "EC50", 100.0, "nM", "="),
                    (2, 2, 1, "CC50", 5.0, "uM", "=")])
    stats = Counter()
    rows = query_pairs(db, ["2697049"], stats=stats)
    assert len(rows) == 1, "a uM counterscreen was dropped"
    assert rows[0]["cc50_nm"] == 5000.0
    assert stats["tox_converted"] == 1


def test_an_unconvertible_unit_is_dropped_and_counted(tmp_path):
    """ug.mL-1 needs a molecular weight. 240 of the 247 non-nM antiviral
    rows are this, so the gate is right -- it just has to be visible."""
    db = _chembl(
        tmp_path,
        assays=[(1, 10, 100, "Antiviral activity in Vero E6 cells"),
                (2, 11, 100, "Cytotoxicity against Vero E6 cells")],
        activities=[(1, 1, 1, "EC50", 1.0, "ug.mL-1", "="),
                    (2, 2, 1, "CC50", 5000.0, "nM", "=")])
    stats = Counter()
    assert query_pairs(db, ["2697049"], stats=stats) == []
    assert stats["av_dropped_unconvertible_units"] == 1
    assert stats["av_unit:ug.mL-1"] == 1


def test_an_exact_cc50_wins_over_a_censored_one(tmp_path):
    """The flag was MIN(CASE relation='=' THEN 0 ELSE 1) -- a minimum, so
    it meant "at least one row is exact" while being named all_exact -- and
    cc50_nm was MIN over every row, censored included. 51 of the 3,388
    SI-contributing groups reported an exact CC50 taken from a censored
    row, and in all 51 the real exact value is HIGHER, so the index was
    understated and the compound looked more cytotoxic than its own exact
    measurement says."""
    db = _chembl(
        tmp_path,
        assays=[(1, 10, 100, "Antiviral activity in Vero E6 cells"),
                (2, 11, 100, "Cytotoxicity against Vero E6 cells")],
        activities=[(1, 1, 1, "EC50", 100.0, "nM", "="),
                    (2, 2, 1, "CC50", 1000.0, "nM", ">"),   # the lower minimum
                    (3, 2, 1, "CC50", 8000.0, "nM", "=")])  # the exact one
    stats = Counter()
    rows = query_pairs(db, ["2697049"], stats=stats)
    assert rows[0]["cc50_nm"] == 8000.0, "the censored bound won the minimum"
    assert rows[0]["cc50_exact"] is True
    assert stats["tox_mixed_group_exact_preferred"] == 1


def test_an_all_censored_group_keeps_its_bound_and_says_so(tmp_path):
    db = _chembl(
        tmp_path,
        assays=[(1, 10, 100, "Antiviral activity in Vero E6 cells"),
                (2, 11, 100, "Cytotoxicity against Vero E6 cells")],
        activities=[(1, 1, 1, "EC50", 100.0, "nM", "="),
                    (2, 2, 1, "CC50", 50000.0, "nM", ">")])
    rows = query_pairs(db, ["2697049"])
    assert rows[0]["cc50_nm"] == 50000.0
    assert rows[0]["cc50_exact"] is False


def test_every_antiviral_row_is_accounted_for(tmp_path):
    """E7's actual demand: emitted + dropped == rows, which only the ChEMBL
    layer had. A gate inside a WHERE clause cannot be given a counter."""
    db = _chembl(
        tmp_path,
        assays=[(1, 10, 100, "Antiviral activity in Vero E6 cells"),
                (2, 11, 100, "Cytotoxicity against Vero E6 cells")],
        activities=[(1, 1, 1, "EC50", 100.0, "nM", "="),
                    (2, 1, 1, "EC50", 900.0, "nM", ">"),        # censored
                    (3, 1, 1, "EC50", 1.0, "ug.mL-1", "="),     # unconvertible
                    (4, 2, 1, "CC50", 5000.0, "nM", "=")])
    stats = Counter()
    query_pairs(db, ["2697049"], stats=stats)
    assert stats["av_rows"] == 3
    assert (stats["av_dropped_censored_or_null_relation"]
            + stats["av_dropped_unconvertible_units"]
            + stats["av_groups"]) == stats["av_rows"]
