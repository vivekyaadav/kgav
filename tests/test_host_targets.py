"""Gate for the drug -> human protein layer: direction is never invented."""
import inspect
import sqlite3
import sys
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
# scripts/ too: the phase-gate tests read ingest_host_targets.build_parser,
# because a default that only a rebuild can observe is a default nothing
# checks. Done here rather than relying on another test file having run
# first -- the defect L16 records.
sys.path.insert(0, str(ROOT / "scripts"))
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
    CREATE TABLE assays(assay_id INT, tid INT, doc_id INT, assay_type TEXT,
                        chembl_id TEXT);
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
    # assay_type is ChEMBL's own classification: B binding, F functional,
    # A ADMET, T toxicity. One of each, because v0.1 wrote "binding" on all
    # 18,954 edges while the rows behind them are 52,485 B, 4,936 F, 1,301 A
    # and 95 T -- and a fixture with only B rows cannot see that.
    con.executemany("INSERT INTO assays VALUES(?,?,?,?,?)",
                    [(1, 10, 1, "B", "CHEMBL_A1"), (2, 11, 1, "B", "CHEMBL_A2"),
                     (3, 10, 2, "F", "CHEMBL_A3"), (4, 10, 3, "A", "CHEMBL_A4")])
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
                     (9, 4, 1, "IC50", 70.0, "nM", "=", 7.1),
                     # New KEYS, from the F and A assays, so the per-edge
                     # assay_type census has more than one value in it.
                     (10, 3, 2, "Kd", 5.0, "nM", "=", 8.1),
                     (11, 4, 2, "EC50", 3.0, "nM", "=", 8.5)])
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
    # The Ki edge's only row is document 1, which has no pubmed_id.
    assert stats["no_publication"] == 1
    assert stats["with_publication"] == stats["edges"] - 1


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


# ----------------------------------------------------------- what assay kind
# ALL 18,954 EDGES IN v0.1 SAY assay_type: "binding", a constant written at
# the emit site, while ChEMBL classifies the 58,817 rows behind them as
# 52,485 B, 4,936 F, 1,301 A and 95 T. src/kgav/chembl.py had found the
# identical bug in the antiviral layer -- "the ingest wrote 'biochemical' on
# every protein row regardless" -- and this layer kept it.
def test_assay_type_comes_from_chembl_not_from_a_constant(built):
    em, _stats, _ = built
    kinds = {e["qualifiers"]["assay_type"] for e in em.edges}
    assert kinds != {"binding"}, "assay_type is still a constant"
    assert kinds == {"binding", "biochemical", "other"}


def test_the_mapping_is_from_chembl_s_code_not_from_standard_type(built):
    """1,868 EC50 rows are B-coded and 1,426 IC50 rows are F-coded, so the
    readout name does not determine the assay kind. A mapping from
    standard_type -- which is what the fix list proposed -- would be wrong
    in both directions at once."""
    em, _stats, _ = built
    # The fixture's EC50 row sits on an ADMET assay and its Kd on a
    # functional one. Derived from standard_type, those would come out
    # "binding" and "binding".
    by_q = {k: e["qualifiers"]["assay_type"]
            for e in em.edges for k in e["qualifiers"]
            if k.endswith("_nm")}
    assert by_q["ec50_nm"] == "other"
    assert by_q["kd_nm"] == "biochemical"


def test_admet_and_toxicity_rows_are_counted_and_not_silent(built):
    """1,301 ADMET and 95 toxicity rows are CYP3A4, CYP2D6, CYP2C9 and hERG
    liability panels. They become TARGETS edges indistinguishable from
    therapeutic target engagement and feed M2-M6. Not dropped -- inhibiting
    CYP3A4 at 100 nM is a true statement -- but a number nobody can see is
    not a disclosure."""
    _em, stats, _ = built
    assert stats["row_chembl_assay_type_A"] == 2  # activities 9 and 11, on assay 4
    # A PARTITION OF EVERY ROW THE QUERY RETURNED, counted before any gate
    # -- the same population the 1,301 figure was measured over.
    assert sum(v for k, v in stats.items()
               if k.startswith("row_chembl_assay_type_")) == stats["rows"]


def test_an_unknown_future_assay_code_becomes_other_and_is_counted(db, host):
    """A code this mapping has never seen must not silently inherit the
    previous row's value or default to 'binding'."""
    con = sqlite3.connect(db)
    con.execute("UPDATE assays SET assay_type = 'Z' WHERE assay_id = 1")
    con.commit()
    con.close()
    index = AliasIndex.from_releases(host)
    em = Emit()
    stats = ingest_host_targets(em, query_host_activities(db), {},
                                host_protein_ids(index), "infores:chembl")
    assert stats["row_chembl_assay_type_Z"] > 0
    assert all(e["qualifiers"]["assay_type"] in
               ("binding", "biochemical", "cytotoxicity", "other")
               for e in em.edges)


def test_the_relation_coercion_is_counted(built):
    """An absent standard_relation became "=" silently and was then
    indistinguishable from a value ChEMBL stated as exact. Zero rows need
    it on chembl_37 -- which is the point: a silent 0 and a silent 40,000
    looked the same from here."""
    _em, stats, _ = built
    assert stats["relation_imputed"] == 0


def test_a_missing_relation_is_imputed_and_says_so(db, host):
    con = sqlite3.connect(db)
    con.execute("UPDATE activities SET standard_relation = NULL "
                "WHERE activity_id = 1")
    con.commit()
    con.close()
    index = AliasIndex.from_releases(host)
    em = Emit()
    stats = ingest_host_targets(em, query_host_activities(db), {},
                                host_protein_ids(index), "infores:chembl")
    assert stats["relation_imputed"] == 1
    assert all(e["qualifiers"]["relation"] in ("=", ">", "<", ">=", "<=", "~")
               for e in em.edges)


def test_every_emitted_assay_type_is_in_the_schema_enum():
    """The schema's assay_type vocabulary was written around antiviral
    assays, so its cellular members all say "antiviral" -- which a CYP3A4
    measurement is not. Checked against the schema rather than a list kept
    here, so a vocabulary change cannot leave this layer writing a value
    the release will reject."""
    from kgav.host_targets import CHEMBL_ASSAY_TYPE
    from kgav.schema import load_schema
    allowed = set(load_schema().enums["assay_type"])
    assert set(CHEMBL_ASSAY_TYPE.values()) <= allowed, (
        f"{set(CHEMBL_ASSAY_TYPE.values()) - allowed} not in the schema enum")


# ------------------------------------------- the four gates inside the WHERE
# Fix-list E7. query_host_activities' WHERE clause applies four gates --
# non-human target, target not a SINGLE PROTEIN, pChEMBL below the floor,
# max_phase below the floor -- and none can count itself from inside a WHERE.
# The last is fix-list B1's subject, a proposal to drop it, and nothing in
# this repository could say what it withholds.
#
# Measured on chembl_37 (a 76s full scan, so the driver runs it behind
# --gate-census): 24,527,044 rows scanned, 10,573,458 not human, 5,654,400
# target not a single protein, 5,530,961 not a potency type, 1,404,641
# pChEMBL below 6, 1,304,767 max_phase below 2, 58,817 kept -- which is
# exactly what the query returns.
def test_the_gate_census_partitions_every_row_scanned(db):
    from kgav.host_targets import gate_census
    c = gate_census(db)
    total = sum(v for k, v in c.items()
                if k.startswith("dropped_") or k == "kept")
    assert total == c["rows_scanned"], (
        f"{c['rows_scanned']} scanned but {total} accounted for; a gate is "
        f"rejecting rows that no counter claims")
    assert c["kept"] == len(query_host_activities(db))


def test_the_census_attributes_each_row_to_one_gate(db):
    """A row failing two gates is attributed to the first, which is what
    makes the six numbers add up rather than double-count."""
    from kgav.host_targets import gate_census
    c = gate_census(db)
    # The fixture's tid 11 is human but NOTINGRAPH; mol 4 is phase 1.
    assert c["dropped_max_phase_below_2"] >= 1
    assert c["kept"] >= 1


def test_the_phase_filter_cost_is_a_number_now(db):
    """B1 asks whether to drop `max_phase >= 2`. This is what it costs."""
    from kgav.host_targets import gate_census
    strict = gate_census(db, min_phase=2)
    loose = gate_census(db, min_phase=0)
    assert loose["kept"] > strict["kept"]
    assert loose["kept"] - strict["kept"] == strict["dropped_max_phase_below_2"]


# ----------------------- P8: "first asserted" has to mean first asserted
#
# About three rows share each (compound, target, assay type) key and the
# first one by activity_id set the edge's value, relation, direction AND
# DATE. ORDER BY activity_id made that choice deterministic across rebuilds
# -- query_host_activities' docstring says why, naming the stake: "the
# second would land in a different half of a temporal split" -- but ChEMBL's
# insertion order is not the same thing as first asserted.
#
# Measured on chembl_37 over 19,020 keys: 766 kept a year LATER than their
# earliest assertion (54 by ten years or more), 143 kept an UNDATED row
# while a dated sibling existed, and 71 changed which side of a 2021 cutoff
# they fall on. On the layer: 453 of 14,381 edges move, 350 earlier and 103
# off the placeholder, and undated falls 2,090 -> 1,987.


def _dup_rows(years, values=None, types="IC50"):
    """Rows sharing one (compound, target, assay type) key.

    In activity_id order, which is the order query_host_activities returns
    and therefore the order the dedup sees.
    """
    values = values or [50.0] * len(years)
    return [{"compound": "CHEMBL1", "max_phase": 4, "smiles": "C",
             "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
             "accession": "Q9BYF1", "standard_type": types,
             "standard_value": v, "standard_units": "nM",
             "standard_relation": "=", "pchembl_value": 7.0,
             "chembl_assay_type": "B", "assay": f"CHEMBL_A{i}",
             "year": y, "pubmed_id": 1000 + i}
            for i, (y, v) in enumerate(zip(years, values))]


def _emit_dups(rows):
    em = Emit()
    stats = ingest_host_targets(em, rows, {}, {ACE2}, "infores:chembl")
    return em, stats


def test_the_date_is_the_earliest_of_the_collapsed_rows():
    """The kept row is still the first by activity_id for value, relation
    and direction; only the date is taken from the earliest assertion."""
    em, stats = _emit_dups(_dup_rows([2019, 2011, 2016]))
    assert stats["edges"] == 1 and stats["duplicate"] == 2
    assert em.edges[0]["first_asserted_date"] == "2011"
    assert stats["date_moved_earlier"] == 1
    # the value still comes from the first row, which is a separate choice
    assert em.edges[0]["qualifiers"]["ic50_nm"] == 50.0


def test_a_dated_duplicate_rescues_an_undated_kept_row():
    """143 keys took the 1970 placeholder while a sibling carried a year.
    An undated edge is excluded or specially handled by every temporal
    protocol, so this is not cosmetic."""
    em, stats = _emit_dups(_dup_rows([None, 2014]))
    assert em.edges[0]["first_asserted_date"] == "2014"
    assert stats["date_recovered_from_a_dated_duplicate"] == 1
    assert stats["date_moved_earlier"] == 0, (
        "recovering from the placeholder is counted separately from moving "
        "an existing year, because they are different events")


def test_an_undated_duplicate_does_not_overwrite_a_year():
    em, stats = _emit_dups(_dup_rows([2014, None]))
    assert em.edges[0]["first_asserted_date"] == "2014"
    assert not stats["date_moved_earlier"]
    assert not stats["date_recovered_from_a_dated_duplicate"]


def test_the_placeholder_is_not_compared_by_value():
    """"2014" < "1970-01-01" is False lexically, so a bare `<` would let
    the placeholder beat every real year. The branch tests it by identity;
    this is the case that catches a regression to `<` alone."""
    em, _ = _emit_dups(_dup_rows([None, 1971]))
    assert em.edges[0]["first_asserted_date"] == "1971", (
        "the 1970 placeholder won a comparison against a real year")


def test_a_real_1970_is_a_year_and_not_the_placeholder():
    em, stats = _emit_dups(_dup_rows([1975, 1970]))
    assert em.edges[0]["first_asserted_date"] == "1970"
    assert stats["date_moved_earlier"] == 1


def test_every_row_undated_still_leaves_the_placeholder():
    em, stats = _emit_dups(_dup_rows([None, None]))
    from kgav.temporal import PLACEHOLDER_DATE
    assert em.edges[0]["first_asserted_date"] == PLACEHOLDER_DATE
    assert not stats["date_recovered_from_a_dated_duplicate"]


# --------------- P13 option (b): the spread travels with the point
#
# The edge keeps the FIRST row's value by activity_id -- deterministic and
# arbitrary -- and nothing on it said so. Measured on chembl_37 over 19,020
# keys: 7,257 have more than one row, 5,517 of those disagree (spread
# median 6.67x, p90 182x, max 55,625x), 4,296 keep a non-median value, and
# 7,116 of the replicate sets come from DIFFERENT assays, so they are not
# repeats of one experiment. The point is unchanged; the range is now
# readable beside it.


def test_a_single_row_still_carries_the_range():
    """potency_n=1 with min==max==the value. Emitting these only for
    replicates would make their absence mean two things -- one
    measurement, or a layer predating this -- and an absent qualifier is
    the ambiguity this project keeps paying for."""
    em, _ = _emit_dups(_dup_rows([2014], values=[50.0]))
    q = em.edges[0]["qualifiers"]
    assert q["potency_n"] == 1
    assert q["potency_min_nm"] == q["potency_max_nm"] == 50.0
    assert q["potency_spread_field"] == "ic50_nm"


def test_the_range_spans_every_collapsed_row_and_the_point_does_not_move():
    em, stats = _emit_dups(
        _dup_rows([2014, 2015, 2016], values=[1.0, 1000.0, 7.0]))
    q = em.edges[0]["qualifiers"]
    assert q["ic50_nm"] == 1.0, "the retained point must not have changed"
    assert q["potency_n"] == 3
    assert (q["potency_min_nm"], q["potency_max_nm"]) == (1.0, 1000.0)
    assert stats["edges_collapsing_several_rows"] == 1
    assert stats["edges_whose_rows_disagree"] == 1


def test_agreeing_rows_are_collapsed_but_not_counted_as_disagreeing():
    _em, stats = _emit_dups(_dup_rows([2014, 2015], values=[7.0, 7.0]))
    assert stats["edges_collapsing_several_rows"] == 1
    assert not stats["edges_whose_rows_disagree"]


def test_the_range_names_the_field_it_describes():
    """WITHOUT THIS THE THREE ABOVE ARE A4 UNDER NEW NAMES. One (drug,
    target) pair can carry an IC50 row and a Ki row -- 3,801 of 14,381 do
    -- and they merge on (subject, predicate, object), so an unnamed range
    would describe one quantity on an edge stating two."""
    em, _ = _emit_dups(_dup_rows([2014], values=[9.0], types="Ki"))
    q = em.edges[0]["qualifiers"]
    assert q["potency_spread_field"] == "ki_nm"
    assert q["ki_nm"] == 9.0 and "ic50_nm" not in q


def test_the_named_field_always_contains_the_point():
    """The invariant, over every assay type this layer emits."""
    for stype in ("IC50", "EC50", "Ki", "Kd"):
        em, _ = _emit_dups(_dup_rows([2014, 2015], values=[2.0, 90.0],
                                     types=stype))
        q = em.edges[0]["qualifiers"]
        named = q["potency_spread_field"]
        assert q["potency_min_nm"] <= q[named] <= q["potency_max_nm"]


def test_an_orphaned_range_is_dropped_by_the_merge():
    """I HAD THE PRECEDENCE BACKWARDS, AND A TEST PASSED ON THE WRONG
    REASON.

    This assertion used to read `index("hosttargets") > index("gtop")`
    with the explanation that a later layer wins. It passed -- the index
    comparison is true -- and the explanation was inverted: add_layer
    merges the INCOMING edge into the existing one and the existing side
    is treated as higher, so the EARLIER layer in LAYER_PRECEDENCE wins
    each qualifier. gtop sits before hosttargets, so gtop's ic50_nm lands
    beside hosttargets' range.

    Measured on the release before the fix: 544 of 14,381 TARGETS edges
    (3.78%) carried a range that did not contain the field it named. The
    case that found it is MTNBRBDFNSGQKB-GXTWGEPZSA-N -> P24941, where
    gtop's 1.159999966621399 sat beside [1.1, 1.1].

    So the invariant is enforced where it can be seen -- at the merge --
    rather than inferred from a list's ordering. The test below checks the
    whole release, which is what caught this; this one checks the rule in
    isolation so a failure says which of the two broke.
    """
    from collections import Counter
    from kgav.assemble import merge_edges

    def edge(quals):
        return {"subject": "INCHIKEY:A", "predicate": "TARGETS",
                "object": "UniProtKB:P1", "qualifiers": dict(quals),
                "primary_knowledge_source": "infores:x", "evidence_tier": 1,
                "first_asserted_date": "2020"}
    held = {"potency_spread_field": "ic50_nm", "potency_n": 1,
            "potency_min_nm": 1.1, "potency_max_nm": 1.1}
    stats = Counter()
    # the winning side states a potency the range does not contain
    merged = merge_edges(edge({"ic50_nm": 1.16}), edge({**held, "ic50_nm": 1.1}),
                         None, stats)
    q = merged["qualifiers"]
    assert q["ic50_nm"] == 1.16
    assert "potency_min_nm" not in q and "potency_spread_field" not in q, (
        "a range survived beside a potency it does not describe")
    assert stats["potency_range_orphaned_by_merge"] == 1

    # and a range that DOES contain the surviving potency is kept
    stats2 = Counter()
    kept = merge_edges(edge({"ic50_nm": 1.1}), edge({**held, "ic50_nm": 1.1}),
                       None, stats2)
    assert kept["qualifiers"]["potency_min_nm"] == 1.1
    assert not stats2["potency_range_orphaned_by_merge"]


def test_the_release_keeps_the_invariant_through_assembly():
    """The pre-assembly tests above cannot see the merge, and the merge is
    where 3,801 pairs collapse several standard_types into one triple."""
    import json
    release = ROOT / "data" / "releases" / "v0.1"
    if not (release / "edges.jsonl").exists():
        pytest.skip("release v0.1 is not built in this checkout")
    checked = 0
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] != "TARGETS":
            continue
        q = e.get("qualifiers") or {}
        if "potency_n" not in q:
            continue
        named = q.get("potency_spread_field")
        assert named, f"a range with no field: {e['subject']} -> {e['object']}"
        assert q.get(named) is not None, (
            f"{e['subject']} -> {e['object']} names {named} as the field its "
            f"range describes and carries no {named}")
        assert q["potency_min_nm"] <= q[named] <= q["potency_max_nm"], (
            f"{e['subject']} -> {e['object']}: {named}={q[named]} outside "
            f"[{q['potency_min_nm']}, {q['potency_max_nm']}]")
        checked += 1
    assert checked, "no TARGETS edge carries a range; the layer is not in"


def test_the_release_lost_no_range_to_a_merge():
    """FIX-LIST P14's CLAIM, where it can be observed.

    host_targets and gtop are the only layers that emit TARGETS and each
    now writes a range over its own rows, so whichever potency wins a merge
    brings the spread belonging to it. The test above says a surviving
    range describes its own number; this one says none had to be dropped to
    make that true -- which is the whole difference between P13's rule and
    P14's fix.

    Before P14, 551 edges arrived here with a potency and no range,
    counted into the MANIFEST as potency_range_orphaned_by_merge.
    """
    import json
    from kgav.assemble import POTENCY_FIELD_ORDER, _RANGE_QUALS
    release = ROOT / "data" / "releases" / "v0.1"
    if not (release / "edges.jsonl").exists():
        pytest.skip("release v0.1 is not built in this checkout")
    fields = set(POTENCY_FIELD_ORDER)
    bare, with_range = [], 0
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] != "TARGETS":
            continue
        q = e.get("qualifiers") or {}
        if not (fields & set(q)):
            continue
        if any(k in q for k in _RANGE_QUALS):
            with_range += 1
            # ALL FOUR OR NONE. A partial group would describe rows no
            # layer collapsed, and the containment check above can pass on
            # a mix by luck.
            assert all(k in q for k in _RANGE_QUALS), (
                f"{e['subject']} -> {e['object']} carries part of the "
                f"group: {sorted(k for k in _RANGE_QUALS if k in q)}")
        else:
            bare.append(f"{e['subject']} -> {e['object']}")
    assert with_range, "no TARGETS edge carries a range; the layers are not in"
    assert not bare, (
        f"{len(bare):,} TARGETS edge(s) state a potency and no range, so a "
        f"spread was dropped at the merge: {bare[:3]}")


# ------------------------------- the analysable-set limit, derived not typed
ANALYSABLE_FIGURE = "<!-- FIGURE: analysable_set -->"


@pytest.mark.skipif(not (ROOT / "data/releases/v0.1/edges.jsonl").exists(),
                    reason="no release on disk")
def test_the_documents_analysable_set_is_the_releases():
    """ALL EIGHT FIGURES RE-DERIVED FROM THE RELEASE.

    The paragraph read 3,002 / 11.8% / 633 / 1,046 / 58 / 108 / 20,536 /
    1,679 for three days after the GtoPdb layer began contributing TARGETS
    edges. Nothing was mismeasured: the numbers were right when written and
    nothing re-derived them when a new source landed. A phase filter, a new
    drug-target source or a changed pChEMBL floor all move these, which is
    exactly when a document stops being checkable by reading it.
    """
    import json
    import re

    rel = ROOT / "data/releases/v0.1"
    doc = (ROOT / "docs" / "scope_and_limits.md").read_text()
    assert ANALYSABLE_FIGURE in doc, "the analysable-set paragraph lost its marker"
    body = doc.split(ANALYSABLE_FIGURE, 1)[1]
    para = " ".join(body.split("\n\n")[0].split())
    ceiling_para = " ".join(body.split("\n\n")[1].split())
    heading = [l for l in doc.splitlines() if "analysable set is" in l][0]

    m = re.search(r"Of ([\d,]+) small molecules, \*\*([\d,]+)\*\* carry a "
                  r"human `TARGETS` edge.*?\*\*([\d,]+)\*\* approved, "
                  r"\*\*([\d,]+)\*\* clinical, \*\*([\d,]+)\*\* phase 1, "
                  r"\*\*([\d,]+)\*\* preclinical and \*\*([\d,]+)\*\* unknown",
                  para)
    assert m, ("the paragraph no longer states the seven figures in the form "
               "this test reads; update both together")
    c = re.search(r"\*\*([\d,]+)\*\* is the ceiling on what one could add — "
                  r"not ([\d,]+)", ceiling_para)
    assert c, "the ceiling sentence no longer states its two figures"
    h = re.search(r"analysable set is ([\d.]+)% of the molecules", heading)
    assert h, "the heading no longer states the percentage"

    def num(x):
        return int(x.replace(",", ""))

    targets = set()
    for line in open(rel / "edges.jsonl"):
        e = json.loads(line)
        if e["predicate"] == "TARGETS":
            targets.add(e["subject"])
    phase = {}
    for line in open(rel / "nodes.jsonl"):
        n = json.loads(line)
        if n["class"] == "SmallMolecule":
            phase[n["id"]] = (n.get("properties") or {}).get("max_phase")

    def bucket(p):
        if p is None:
            return "unknown"
        if p == 4:
            return "approved"
        if p in (2, 3):
            return "clinical"
        if p == 1:
            return "phase 1"
        return "preclinical"

    without = {k: 0 for k in
               ("approved", "clinical", "phase 1", "preclinical", "unknown")}
    for cid, p in phase.items():
        if cid not in targets:
            without[bucket(p)] += 1
    routed = len(set(phase) & targets)

    derived = (len(phase), routed, without["approved"], without["clinical"],
               without["phase 1"], without["preclinical"], without["unknown"])
    assert tuple(num(g) for g in m.groups()) == derived, (
        f"the document says {tuple(m.groups())}; the release says {derived}")

    ceiling = without["approved"] + without["clinical"]
    assert (num(c.group(1)), num(c.group(2))) == (ceiling,
                                                  len(phase) - routed), (
        f"the document's ceiling sentence says {c.groups()}; the release says "
        f"{(ceiling, len(phase) - routed)}")
    assert float(h.group(1)) == round(100 * routed / len(phase), 1), (
        f"the heading says {h.group(1)}%; the release says "
        f"{round(100 * routed / len(phase), 1)}%")


# ======================================================= the phase gate (B1)
@pytest.fixture
def gate_db(tmp_path):
    """A compound at each phase state, so the gate can be tested on all three.

    Its own fixture rather than an addition to `db`: a dozen tests assert
    counts off that one, and widening it would have made them fail for a
    reason unrelated to what they check.

    P4 phase 4 (clears the floor), P1 phase 1 (below it), PN max_phase NULL
    (never in a trial) -- all three with an identical qualifying activity
    against the same human single-protein target, so phase is the ONLY
    difference between them.
    """
    p = tmp_path / "gate.db"
    con = sqlite3.connect(p)
    con.executescript("""
    CREATE TABLE molecule_dictionary(molregno INT, chembl_id TEXT, max_phase INT);
    CREATE TABLE compound_structures(molregno INT, canonical_smiles TEXT, standard_inchi_key TEXT);
    CREATE TABLE target_dictionary(tid INT, chembl_id TEXT, target_type TEXT, tax_id INT);
    CREATE TABLE target_components(tid INT, component_id INT);
    CREATE TABLE component_sequences(component_id INT, accession TEXT);
    CREATE TABLE assays(assay_id INT, tid INT, doc_id INT, assay_type TEXT, chembl_id TEXT);
    CREATE TABLE docs(doc_id INT, year INT, pubmed_id INT);
    CREATE TABLE activities(activity_id INT, assay_id INT, molregno INT, standard_type TEXT,
      standard_value REAL, standard_units TEXT, standard_relation TEXT, pchembl_value REAL);
    CREATE TABLE drug_mechanism(molregno INT, tid INT, action_type TEXT);
    """)
    con.executemany("INSERT INTO molecule_dictionary VALUES(?,?,?)",
                    [(1, "P4", 4), (2, "P1", 1), (3, "PN", None)])
    con.executemany("INSERT INTO compound_structures VALUES(?,?,?)",
                    [(1, "C", "K4" + "X" * 12 + "-" + "Y" * 10 + "-N"),
                     (2, "C", "K1" + "X" * 12 + "-" + "Y" * 10 + "-N"),
                     (3, "C", "KN" + "X" * 12 + "-" + "Y" * 10 + "-N")])
    con.execute("INSERT INTO target_dictionary VALUES(10,'T','SINGLE PROTEIN',9606)")
    con.execute("INSERT INTO target_components VALUES(10,100)")
    con.execute("INSERT INTO component_sequences VALUES(100,'Q9BYF1')")
    con.execute("INSERT INTO assays VALUES(1,10,1,'B','A1')")
    con.execute("INSERT INTO docs VALUES(1,2020,111)")
    con.executemany("INSERT INTO activities VALUES(?,?,?,?,?,?,?,?)",
                    [(1, 1, 1, "IC50", 50.0, "nM", "=", 7.0),
                     (2, 1, 2, "IC50", 50.0, "nM", "=", 7.0),
                     (3, 1, 3, "IC50", 50.0, "nM", "=", 7.0)])
    con.commit()
    con.close()
    return p


def _key(prefix):
    return prefix + "X" * 12 + "-" + "Y" * 10 + "-N"


def _compounds(rows):
    return {r["compound"] for r in rows}


def test_the_floor_alone_keeps_only_the_phased_compound(gate_db):
    """How v0.1 was built, and the baseline the relaxation is measured from.
    `NULL >= 2` is NULL in SQL, so the unphased compound is excluded by the
    floor without anything saying so."""
    rows = query_host_activities(gate_db, 6.0, 2)
    assert _compounds(rows) == {"P4"}


def test_the_relaxed_gate_admits_a_compound_the_graph_already_has(gate_db):
    """B1 OPTION D. An early-phase compound already in the graph gains its
    human targets; this is the 1,327 compounds and 300 positive labels that
    had no route into the host side at all."""
    rows = query_host_activities(gate_db, 6.0, 2, known={_key("K1")})
    assert _compounds(rows) == {"P4", "P1"}


def test_the_relaxed_gate_admits_an_UNPHASED_compound_in_the_graph(gate_db):
    """THE NULL CASE, WHICH IS MOST OF CHEMBL. 678,515 of the 681,988
    compounds behind this gate have no max_phase at all, so a relaxation
    that handled only `< floor` would have missed almost the whole point.
    `NULL OR TRUE` is TRUE, which is why the OR works."""
    rows = query_host_activities(gate_db, 6.0, 2, known={_key("KN")})
    assert _compounds(rows) == {"P4", "PN"}


def test_the_gate_is_RELAXED_AND_NOT_REMOVED(gate_db):
    """The distinction between B1's option D and option C. A compound below
    the floor that the graph does NOT have stays out -- otherwise this
    admits 676,997 new nodes, which buys candidate pool and no AUC."""
    rows = query_host_activities(gate_db, 6.0, 2, known={_key("KZ")})
    assert _compounds(rows) == {"P4"}, (
        "a compound outside the known set came in below the floor; the gate "
        "is being removed rather than relaxed")


def test_an_empty_known_set_is_the_floor_and_not_everything(gate_db):
    """A checkout with no earlier layers on disk yields an empty set. That
    must mean "relax for nobody", not "relax for everybody" -- the failure
    mode would be a 27x graph on a machine that simply had no layers built."""
    for empty in (set(), None):
        assert _compounds(query_host_activities(gate_db, 6.0, 2,
                                                known=empty)) == {"P4"}


def test_known_inchikeys_reads_small_molecules_and_strips_the_prefix(tmp_path):
    """The set is compared against compound_structures.standard_inchi_key in
    SQL, so it must NOT carry the node prefix; and a Protein node in the
    same file must not leak into it."""
    import json

    from kgav.host_targets import known_inchikeys
    d = tmp_path / "v0.1-chembl"
    d.mkdir()
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in [
        {"id": f"INCHIKEY:{_key('K1')}", "class": "SmallMolecule"},
        {"id": "UniProtKB:Q9BYF1", "class": "Protein"},
        {"id": "CHEMBL:CHEMBL25", "class": "SmallMolecule"},
    ]) + "\n")
    got = known_inchikeys(tmp_path, ["chembl"])
    assert got == {_key("K1")}, (
        "either the prefix survived, a Protein leaked in, or a non-INCHIKEY "
        "SmallMolecule id was admitted unprefixed")


def test_known_inchikeys_skips_a_layer_that_is_not_built(tmp_path):
    """data/raw is populated by hand, so an absent layer is an ordinary
    state of a checkout and must not raise."""
    from kgav.host_targets import known_inchikeys
    assert known_inchikeys(tmp_path, ["chembl", "ncats"]) == set()


def test_the_known_set_comes_from_the_layers_before_hosttargets():
    """NOT A HAND-WRITTEN LIST. A literal would go stale the next time a
    layer moves, silently narrowing or widening the gate."""
    from kgav.pipeline import LAYERS, layers_before
    before = layers_before("hosttargets")
    names = [l.name for l in LAYERS]
    assert before == tuple(names[:names.index("hosttargets")])
    assert "chembl" in before and "orcs" in before


def test_the_relaxed_gate_CANNOT_SEE_ncats_FROM_WHERE_hosttargets_SITS():
    """A LIMIT RECORDED, NOT ENFORCED. 660 of the 1,327 compounds option D
    recovers come from ncats, which is built AFTER this layer, so
    `--phase-gate relaxed` reaches the 667 chembl ones. hosttargets was
    moved past ncats to fix that and moved back when B1 was declined.

    This asserts the state the shipped order is in, so that anyone who
    moves the layer again is told what it changes rather than finding out
    from a figure. If ncats ever precedes hosttargets, the relaxed policy
    silently doubles its reach -- which is a decision, not a refactor."""
    from kgav.pipeline import layers_before
    assert "ncats" not in layers_before("hosttargets"), (
        "ncats now precedes hosttargets, so --phase-gate relaxed admits "
        "roughly twice the compounds it did when B1 was measured; re-read "
        "host_targets.PHASE_GATE_POLICIES before accepting this")


def test_the_drivers_default_phase_gate_is_the_decision():
    """The default IS the decision (the rule L10 settled), and an
    unrecognised value is refused rather than read as either policy."""
    import ingest_host_targets as H
    from kgav.host_targets import DEFAULT_PHASE_GATE, PHASE_GATE_POLICIES
    ap = H.build_parser(Path("/nowhere"))
    act = {a.dest: a for a in ap._actions}["phase_gate"]
    assert act.default == DEFAULT_PHASE_GATE == "floor", (
        "the shipped policy is the floor; B1's relaxation was measured and "
        "declined because its admission rule correlates with the labels")
    assert tuple(act.choices) == PHASE_GATE_POLICIES == ("floor", "relaxed")
    with pytest.raises(SystemExit):
        ap.parse_args(["--phase-gate", "relaxed-ish"])
