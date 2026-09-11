"""Day 7 gate: nsp resolution from assay text, censored values preserved,
units never guessed.
"""
import sqlite3
from pathlib import Path

import pytest

from kgav.aliases import AliasIndex
from kgav.chembl import (
    build_viral_target_index,
    ingest_activities,
    query_activities,
    resolve_nsp,
    to_nm,
)
from kgav.emit import Emit
from kgav.schema import load_schema

ROOT = Path(__file__).resolve().parents[1]
CHEMBL_DB = ROOT / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db"
needs_db = pytest.mark.skipif(not CHEMBL_DB.exists(), reason="ChEMBL not downloaded")


# --------------------------------------------------------------- nsp mapping
@pytest.mark.parametrize("desc,nsp,domain", [
    ("Inhibition of SARS-CoV-2 3CLpro using FRET substrate", "nsp5", None),
    ("Inhibition of SARS-CoV-2 main protease by FRET assay", "nsp5", None),
    ("Inhibition of SARS-CoV-2 MPro", "nsp5", None),
    ("Inhibition of SARS-CoV-2 M-protease", "nsp5", None),
    ("Inhibition of 6-his tagged SARS-Cov-2 3-CL protease", "nsp5", None),
    ("Inhibition of SARS-CoV-2 PL protease", "nsp3", "PLpro"),
    ("Inhibition of C-terminal His-tagged PLpro", "nsp3", "PLpro"),
    ("Inhibition of SARS-CoV-2 RdRp", "nsp12", None),
    ("Inhibition of nsp13 helicase", "nsp13", None),
    ("Some assay with no target named", None, None),
])
def test_nsp_resolution_from_description(desc, nsp, domain):
    assert resolve_nsp(desc) == (nsp, domain)


def test_macrodomain_is_distinguished_from_plpro():
    """nsp3 carries both PLpro and the ADP-ribose macrodomain -- separate
    druggable sites, 623 and 1,479 activities. A macrodomain binder says
    nothing about protease inhibition, so they must not collapse."""
    nsp, dom = resolve_nsp("Functional biochemical assay that targets "
                           "SARS-COV-2 nsp3 macrodomain, displacing ADP-Ribose")
    assert (nsp, dom) == ("nsp3", "macrodomain")
    assert resolve_nsp("Inhibition of SARS-CoV-2 PL protease")[1] == "PLpro"


# --------------------------------------------------------------------- units
@pytest.mark.parametrize("value,unit,expected", [
    (25.0, "nM", 25.0), (1.5, "uM", 1500.0), (1000.0, "pM", 1.0), (1.0, "mM", 1_000_000.0),
])
def test_unit_conversion(value, unit, expected):
    assert to_nm(value, unit) == expected


def test_mass_concentration_is_not_guessed():
    """ug.mL-1 needs a molecular weight. Guessing one produces a number that
    looks measured."""
    assert to_nm(1.0, "ug.mL-1") is None
    assert to_nm(1.0, None) is None
    assert to_nm(1.0, "%") is None


# ------------------------------------------------------------------- fixture
@pytest.fixture
def fake_db(tmp_path):
    p = tmp_path / "c.db"
    con = sqlite3.connect(p)
    con.executescript("""
    CREATE TABLE molecule_dictionary(molregno INT, chembl_id TEXT, max_phase INT, pref_name TEXT);
    CREATE TABLE compound_structures(molregno INT, canonical_smiles TEXT, standard_inchi_key TEXT);
    CREATE TABLE target_dictionary(tid INT, chembl_id TEXT, pref_name TEXT, target_type TEXT, tax_id INT);
    CREATE TABLE target_components(tid INT, component_id INT);
    CREATE TABLE component_sequences(component_id INT, accession TEXT);
    CREATE TABLE assays(assay_id INT, tid INT, description TEXT, assay_type TEXT, chembl_id TEXT, doc_id INT);
    CREATE TABLE docs(doc_id INT, year INT, pubmed_id INT);
    CREATE TABLE activities(activity_id INT, assay_id INT, molregno INT, standard_type TEXT,
                            standard_value REAL, standard_units TEXT, standard_relation TEXT);
    """)
    con.executemany("INSERT INTO molecule_dictionary VALUES(?,?,?,?)",
                    [(1, "CHEMBL1", 4, "a"), (2, "CHEMBL2", 4, "b"), (3, "CHEMBL3", None, "c")])
    con.executemany("INSERT INTO compound_structures VALUES(?,?,?)",
                    [(1, "CC1", "AAAAAAAAAAAAAA-BBBBBBBBBB-C"),
                     (2, "CC2", "DDDDDDDDDDDDDD-EEEEEEEEEE-F"), (3, "CC3", None)])
    con.executemany("INSERT INTO target_dictionary VALUES(?,?,?,?,?)",
                    [(10, "CHEMBL4523582", "Replicase polyprotein 1ab", "SINGLE PROTEIN", 2697049),
                     (11, "CHEMBL4303835", "SARS-CoV-2", "ORGANISM", 2697049)])
    con.executemany("INSERT INTO target_components VALUES(?,?)", [(10, 100)])
    con.executemany("INSERT INTO component_sequences VALUES(?,?)", [(100, "P0DTD1")])
    con.executemany("INSERT INTO assays VALUES(?,?,?,?,?,?)",
                    [(1, 10, "Inhibition of SARS-CoV-2 3CLpro FRET", "B", "A1", 1),
                     (2, 10, "assay targeting SARS-COV-2 nsp3 macrodomain ADP-Ribose", "B", "A2", 1),
                     (5, 10, "Inhibition of an unspecified viral enzyme", "B", "A5", 1),
                     (6, 11, "Cell based antiviral assay against SARS-CoV-2", "F", "A6", 1)])
    con.executemany("INSERT INTO docs VALUES(?,?,?)", [(1, 2021, 33333333)])
    con.executemany("INSERT INTO activities VALUES(?,?,?,?,?,?,?)",
                    [(1, 1, 1, "IC50", 25.0, "nM", "="),
                     (2, 2, 2, "IC50", 1.5, "uM", "="),
                     (5, 5, 1, "IC50", 10000.0, "nM", ">"),
                     (6, 6, 2, "EC50", 500.0, "nM", "="),
                     (7, 6, 3, "EC50", 50.0, "nM", "="),
                     (8, 6, 2, "EC50", 1.0, "ug.mL-1", "=")])
    con.commit()
    con.close()
    return p


@pytest.fixture
def spine(tmp_path):
    import json
    d = tmp_path / "spine"
    d.mkdir()

    def prot(pid, fam=None):
        props = {"taxon_id": "NCBITaxon:2697049", "is_viral": True,
                 "sequence_hash": "x", "reviewed": True}
        if fam:
            props |= {"protein_family": fam, "mature_peptide": fam}
        return {"id": pid, "class": "Protein", "properties": props}

    nodes = [{"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
              "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                             "baltimore_class": "IV", "is_enveloped": True}},
             prot("UniProtKB:P0DTD1"), prot("UniProtKB:PRO_1", "nsp5"),
             prot("UniProtKB:PRO_3", "nsp3")]
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (d / "edges.jsonl").write_text("")
    return d


@pytest.fixture
def built(fake_db, spine):
    index = AliasIndex.from_releases(spine)
    tmap = {"2697049": "2697049"}
    targets = build_viral_target_index(index, tmap)
    rows = query_activities(fake_db, ["2697049"])
    em = Emit()
    stats = ingest_activities(em, rows, targets, tmap, "infores:chembl")
    return em, stats, index


def test_activities_route_to_the_right_edge_type(built):
    _em, stats, _ = built
    assert stats["protein_edges"] == 3
    assert stats["organism_edges"] == 1


def test_resolved_activities_attach_to_the_chain_not_the_polyprotein(built):
    em, _, _ = built
    objs = {e["object"] for e in em.edges if e["predicate"] == "INHIBITS"}
    assert "UniProtKB:PRO_1" in objs      # 3CLpro -> nsp5
    assert "UniProtKB:PRO_3" in objs      # macrodomain -> nsp3


def test_unresolved_activities_are_flagged_not_hidden(built):
    """An unresolved edge must never be mistaken for a measured one."""
    em, _, _ = built
    unres = [e for e in em.edges
             if e["qualifiers"].get("target_resolution") == "parent_unresolved"]
    assert len(unres) == 1
    assert unres[0]["object"] == "UniProtKB:P0DTD1"


def test_censored_values_keep_their_relation(built):
    """IC50 > 10000 nM is a NEGATIVE result. Coercing it to a number inverts
    its meaning -- and these are the best hard negatives available."""
    em, stats, _ = built
    assert stats["censored"] == 1
    censored = [e for e in em.edges if e["qualifiers"].get("relation") == ">"]
    assert len(censored) == 1
    assert censored[0]["qualifiers"]["ic50_nm"] == 10000.0


def test_every_activity_edge_is_flagged_unquantified(built):
    """ChEMBL has 2 CC50 records across 21,900 coronavirus activities, so no
    selectivity index is computable for any of them."""
    em, _, _ = built
    for e in em.edges:
        assert e["qualifiers"]["unquantified"] is True
        assert "selectivity_index" not in e["qualifiers"]


def test_rows_without_an_inchikey_are_dropped(built):
    _em, stats, _ = built
    assert stats["no_inchikey"] == 1


def test_unconvertible_units_are_dropped(built):
    _em, stats, _ = built
    assert stats["unconvertible_units"] == 1


def test_domain_qualifier_records_the_binding_site(built):
    em, _, _ = built
    doms = {e["qualifiers"].get("domain") for e in em.edges}
    assert "macrodomain" in doms


def test_output_validates_against_schema(built):
    em, _, index = built
    nodes = list(index.nodes.values()) + list(em.nodes.values())
    violations = load_schema().validate_batch(nodes, em.edges)
    assert violations == [], [str(v) for v in violations[:6]]


# ----------------------------------------------------------------- real data
@needs_db
def test_real_chembl_scale():
    rows = query_activities(CHEMBL_DB, ["2697049", "694009", "1335626", "11137",
                                        "277944", "31631", "290028", "1263720", "443239"])
    assert len(rows) > 15_000, f"{len(rows)} activities"
    resolved = sum(1 for r in rows if resolve_nsp(r.get("description") or "")[0])
    protein_rows = [r for r in rows if r["target_type"] != "ORGANISM"]
    rate = resolved / max(len(protein_rows), 1)
    assert rate > 0.6, f"only {rate:.0%} of protein activities resolve to an nsp"
