"""Day 7 gate: nsp resolution from assay text, censored values preserved,
units never guessed.
"""
import sqlite3
from pathlib import Path

import pytest

from kgav.aliases import AliasIndex
import kgav.chembl as chembl_mod
from kgav.chembl import (
    accounting,
    assay_readout,
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
    """The third viral-protein row is "IC50 > 10000 nM" -- a compound measured
    NOT to inhibit. It used to be counted as an INHIBITS edge like the other
    two; the polarity of the measurement now decides the predicate."""
    _em, stats, _ = built
    assert stats["protein_edges"] == 2
    assert stats["measured_inactive_edges"] == 1
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


def test_this_layer_does_not_state_unquantified(built):
    """It used to write unquantified=true on every edge, because CC50 is
    absent from ChEMBL's coronavirus activities (2 records in 21,900).

    That is a fact about THIS LAYER, not about the assembled edge. The
    selectivity layer supplies same-document CC50 for 1,568 pairs, so the
    claim was false for 1,264 of them once the graph was assembled -- and it
    won the merge, because chembl precedes selectivity. The field is derived
    after merging now, and no ingest may assert it.
    """
    em, _, _ = built
    for e in em.edges:
        assert "unquantified" not in e["qualifiers"]
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


# ------------------------------------------------- polarity decides predicate
def _prow(stype, value, relation, units="nM"):
    """One viral-protein measurement row, as query_activities returns it."""
    return {"compound": "CHEMBL9", "smiles": "CC", "max_phase": 4, "drug_name": "x",
            "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "target": "CHEMBL4523582",
            "target_name": "Replicase polyprotein 1ab", "target_type": "SINGLE PROTEIN",
            "tax_id": 2697049, "accession": "P0DTD1", "standard_type": stype,
            "standard_value": value, "standard_units": units,
            "standard_relation": relation, "description": "Inhibition of SARS-CoV-2 3CLpro",
            "assay_type": "B", "assay": f"A-{stype}-{value}-{relation}",
            "year": 2021, "pubmed_id": 1}


def _one(row, spine):
    index = AliasIndex.from_releases(spine)
    targets = build_viral_target_index(index, {"2697049": "2697049"})
    em = Emit()
    stats = ingest_activities(em, [row], targets, {"2697049": "2697049"}, "infores:chembl")
    return em, stats


def test_a_measured_non_inhibition_is_not_an_INHIBITS_edge(spine):
    """1,433 of 5,486 pairs had NO active measurement and still asserted
    inhibition, so 26% of the channel M1, M7 and M8 walk was evidence that the
    compound does not inhibit."""
    em, stats = _one(_prow("IC50", 50_000.0, "="), spine)
    assert [e["predicate"] for e in em.edges] == ["MEASURED_INACTIVE_AGAINST"]
    assert stats["measured_inactive_edges"] == 1 and stats["protein_edges"] == 0
    assert em.edges[0]["qualifiers"]["ic50_nm"] == 50_000.0, "the measurement is kept"


def test_a_censored_bound_above_the_threshold_is_inactive(spine):
    """"IC50 > 10000 nM" is a compound assayed and found inactive."""
    em, _ = _one(_prow("IC50", 10_000.0, ">"), spine)
    assert [e["predicate"] for e in em.edges] == ["MEASURED_INACTIVE_AGAINST"]


def test_a_weak_lower_bound_decides_nothing_and_is_counted(spine):
    """"IC50 > 100 nM" may sit below the assay's top concentration. It proves
    neither side, so it becomes no edge rather than a guessed one."""
    em, stats = _one(_prow("IC50", 100.0, ">"), spine)
    assert em.edges == []
    assert stats["undecidable_protein_measurement"] == 1


def test_a_binding_constant_is_read_rather_than_discarded(spine):
    """classify() read only EC50/IC50, so all 667 Ki/Kd measurements returned
    None. Routing on that verdict would have deleted them -- 504 are sub-10uM
    binders, the strongest direct-acting evidence in the layer."""
    em, stats = _one(_prow("Ki", 5.0, "="), spine)
    assert [e["predicate"] for e in em.edges] == ["INHIBITS"]
    assert stats["protein_edges"] == 1
    em2, _ = _one(_prow("Kd", 80_000.0, "="), spine)
    assert [e["predicate"] for e in em2.edges] == ["MEASURED_INACTIVE_AGAINST"]


def test_binding_constants_never_reach_the_organism_labels(built):
    """Widening classify() is label-neutral only because Ki/Kd on an ORGANISM
    target are dropped as mis-entered. If that ever changes, the labels move."""
    em, _, _ = built
    for e in em.edges:
        if e["predicate"] == "HAS_ANTIVIRAL_ACTIVITY_AGAINST":
            q = e["qualifiers"]
            assert q.get("ki_nm") is None and q.get("kd_nm") is None


# --------------------------------------------------------------------------
# Row accounting: every input row must be an emitted edge or a counted drop.
# --------------------------------------------------------------------------
def _row(**kw):
    base = {"compound": "CHEMBL1", "smiles": "C", "inchikey": "K" * 27,
            "max_phase": 4, "drug_name": "d", "target": "CHEMBL_T",
            "target_name": "Replicase polyprotein 1ab", "target_type": "PROTEIN",
            "tax_id": 2697049, "accession": "P0DTD1", "standard_type": "IC50",
            "standard_value": 50.0, "standard_units": "nM",
            "standard_relation": "=", "description": "Inhibition of 3CLpro",
            "assay_type": "B", "assay": "CHEMBL_A1", "year": 2021,
            "pubmed_id": 1}
    base.update(kw)
    return base


def test_every_row_is_emitted_or_counted_as_dropped():
    """The invariant the driver asserts, at the level that decides it.

    ce133f1 added measured_inactive_edges and undecidable_protein_measurement
    and the driver printed neither, so INHIBITS fell 26% with nothing on screen
    saying where the rest went. Comparing counts across things that should
    agree is this project's stated primary integrity check.
    """
    rows = [
        _row(assay="A1"),                                   # -> INHIBITS
        _row(assay="A2", standard_value=80000.0),           # -> MEASURED_INACTIVE
        _row(assay="A3", standard_relation=">",
             standard_value=100.0),                         # -> undecidable, dropped
        _row(assay="A4", inchikey=None),                    # -> no_inchikey
        _row(assay="A5", tax_id=9999),                      # -> unmapped_taxon
        _row(assay="A6", standard_units="ug.mL-1"),         # -> unconvertible_units
        _row(assay="A7", target_name="cereblon/Replicase"), # -> multicomponent
        _row(assay="A8", target_type="ORGANISM",
             standard_type="Ki"),                           # -> binding on organism
        _row(assay="A9", target_type="ORGANISM",
             standard_type="EC50"),                         # -> organism edge
        _row(assay="A1"),                                   # -> duplicate
    ]
    em = Emit()
    s = ingest_activities(em, rows, {"2697049": {"nsp5": "UniProtKB:NSP5"}},
                          {"2697049": "2697049"}, "infores:chembl")

    emitted = ("protein_edges", "measured_inactive_edges", "organism_edges")
    dropped = ("no_inchikey", "unmapped_taxon", "unconvertible_units",
               "no_target_node", "duplicate", "multicomponent_target",
               "binding_constant_on_organism", "undecidable_protein_measurement")
    assert sum(s[k] for k in emitted) + sum(s[k] for k in dropped) == s["rows"]

    # And each branch is genuinely reachable, so the invariant is not holding
    # because the fixture only exercises one path.
    assert s["protein_edges"] == 1
    assert s["measured_inactive_edges"] == 1
    assert s["undecidable_protein_measurement"] == 1
    assert s["organism_edges"] == 1
    assert s["binding_constant_on_organism"] == 1
    assert s["multicomponent_target"] == 1
    assert s["duplicate"] == 1


def test_every_counter_the_module_sets_has_a_declared_fate():
    """A new counter in chembl.py must land in EMITTED, DROPPED, or a
    breakdown prefix.

    STATIC, by design: it reads the source rather than running the ingest, so
    it also covers counters on branches no test exercises. That is why it
    caught reattributed_to_organism being unlisted, and why it would have
    failed on ce133f1 rather than four months later.

    It used to grep the DRIVER for the counter names. The lists now live in
    chembl.py beside the code that increments them -- because they drifted
    from the driver, which is the whole reason this test exists -- so it
    checks against the lists themselves.
    """
    import re
    src = (ROOT / "src/kgav/chembl.py").read_text()
    # Plain counters, and the prefix of f-string ones: stats[f"resolved_{nsp}"]
    plain = set(re.findall(r'stats\["([a-z_0-9]+)"\]', src))
    prefixed = set(re.findall(r'stats\[f"([a-z_]+)\{', src))
    declared = set(chembl_mod.EMITTED_COUNTERS) | set(chembl_mod.DROPPED_COUNTERS) \
        | set(chembl_mod.BREAKDOWN_COUNTERS)
    prefixes = chembl_mod.BREAKDOWN_PREFIXES

    missing = {c for c in plain
               if c not in declared and not c.startswith(prefixes)}
    assert not missing, (
        f"chembl.py sets counters with no declared fate: {sorted(missing)}. "
        f"Add each to EMITTED_COUNTERS or DROPPED_COUNTERS, or give it a "
        f"breakdown prefix if it describes rows rather than deciding them.")

    # An f-string counter must match a declared prefix, or it is invisible to
    # the accounting and to this check both.
    stray = {pfx for pfx in prefixed if not pfx.startswith(prefixes)}
    assert not stray, f"f-string counters with no declared prefix: {sorted(stray)}"


# ------------------- what the measurement was made ON, not just what it says
@pytest.mark.parametrize("desc,code,where,atype", [
    # The real chloroquine rows, verbatim from ChEMBL 37.
    ("Inhibition of spike glycoprotein S in SARS-CoV-2 pseudovirus infected in "
     "human 293T/ACE2 cells assessed as inhibition of viral infection",
     "B", "cells", "reporter"),
    ("Inhibition of SARS-CoV-2 spike glycoprotein/ His tagged human ACE2 "
     "interaction incubated for 1 hr followed by substrate addition and "
     "measured after 15 mins by ELISA", "B", "protein", "biochemical"),
    # B-coded and cellular: the reason the code cannot be the discriminator.
    ("Inhibition of eGFP fused SARS-CoV2 main protease transfected in "
     "HEK293T/17 cells incubated for 72 hrs by fluorescence based flow "
     "cytometry analysis", "B", "cells", "reporter"),
    ("Inhibition of recombinant SARS-CoV-2 3CL protease using a fluorogenic "
     "substrate", "B", "protein", "biochemical"),
    ("Antiviral activity against SARS-CoV-2 in Vero E6 cells by plaque "
     "reduction", "F", "cells", "plaque_reduction"),
    ("Antiviral activity against SARS-CoV-2 in Vero E6 cells measured as "
     "cytopathic effect", "F", "cells", "cell_based_antiviral"),
    # A lysate is cell-derived and protein-level. Left for a human.
    ("Inhibition of 3CLpro activity in Calu-3 cell lysate", "B",
     "unclear", "other"),
    # Nothing either way: the code is the only signal left.
    ("IC50 against nsp5", "B", "protein", "binding"),
    ("IC50 against nsp5", "F", "protein", "biochemical"),
])
def test_assay_readout_reads_the_description_not_the_code(desc, code, where, atype):
    assert assay_readout(desc, code) == (where, atype)


def test_the_assay_type_code_alone_would_get_it_wrong():
    """444 B-coded rows describe a cellular readout. Mapping B -> binding and
    F -> functional, as the obvious fix would, mislabels every one."""
    cellular_but_coded_binding = (
        "Inhibition of eGFP fused SARS-CoV2 main protease transfected in "
        "HEK293T/17 cells by flow cytometry")
    assert assay_readout(cellular_but_coded_binding, "B")[0] == "cells"


def _viral_row(desc, value=160.0, stype="IC50", code="B", assay="CHEMBL_A1"):
    return {"compound": "CHEMBL76", "smiles": "C", "inchikey": "K" * 27,
            "max_phase": 4, "drug_name": "chloroquine", "target": "CHEMBL3927",
            "target_name": "Spike glycoprotein", "target_type": "SINGLE PROTEIN",
            "tax_id": "2697049", "accession": "P0DTC2",
            "standard_type": stype, "standard_value": value,
            "standard_units": "nM", "standard_relation": "=",
            "description": desc, "assay_type": code, "assay": assay,
            "year": 2021, "pubmed_id": "33560122"}


PSEUDOVIRUS = ("Inhibition of spike glycoprotein S in SARS-CoV-2 pseudovirus "
               "infected in human 293T/ACE2 cells assessed as inhibition of "
               "viral infection")
ELISA = ("Inhibition of SARS-CoV-2 spike glycoprotein/ His tagged human ACE2 "
         "interaction measured by ELISA")


def test_a_cellular_row_is_emitted_on_the_organism_not_the_protein():
    """THE FIX. A pseudovirus entry assay does not show the compound binds
    spike, it shows the compound blocks spike-mediated entry -- and
    chloroquine blocks endosomal acidification."""
    em = Emit()
    s = ingest_activities(em, [_viral_row(PSEUDOVIRUS)],
                          {"2697049": {}}, {"2697049": "2697049"},
                          "infores:chembl")
    assert s["reattributed_to_organism"] == 1
    assert s["protein_edges"] == 0
    e = em.edges[0]
    assert e["predicate"] == "HAS_ANTIVIRAL_ACTIVITY_AGAINST"
    assert e["object"] == "NCBITaxon:2697049"
    assert e["qualifiers"]["assay_type"] == "reporter"


def test_the_reattribution_records_which_protein_the_source_named():
    """Auditable rather than silent: the nominal target is why the row
    exists at all."""
    em = Emit()
    ingest_activities(em, [_viral_row(PSEUDOVIRUS)], {"2697049": {}},
                      {"2697049": "2697049"}, "infores:chembl")
    assert em.edges[0]["qualifiers"]["nominal_protein_target"] == "UniProtKB:P0DTC2"


def test_a_real_binding_row_still_lands_on_the_protein():
    """And it is the honest number: chloroquine's spike ELISA is 7,000 nM,
    where the pseudovirus surrogate reads 160 nM and won the best-potency
    selection by 44x."""
    em = Emit()
    s = ingest_activities(em, [_viral_row(ELISA, value=7000.0)],
                          {"2697049": {}}, {"2697049": "2697049"},
                          "infores:chembl")
    assert s["reattributed_to_organism"] == 0
    assert em.edges[0]["predicate"] == "INHIBITS"
    assert em.edges[0]["object"] == "UniProtKB:P0DTC2"
    assert em.edges[0]["qualifiers"]["ec50_nm" if False else "ic50_nm"] == 7000.0


def test_an_unclear_row_is_left_on_the_protein():
    """Moving a measurement requires knowing it was cellular. 537 rows say
    both, and those are the ones a human has to read."""
    em = Emit()
    s = ingest_activities(em, [_viral_row("Inhibition of 3CLpro activity in "
                                          "Calu-3 cell lysate")],
                          {"2697049": {}}, {"2697049": "2697049"},
                          "infores:chembl")
    assert s["reattributed_to_organism"] == 0
    assert em.edges[0]["predicate"] == "INHIBITS"
    assert em.edges[0]["qualifiers"]["assay_type"] == "other"


def test_no_protein_edge_still_claims_biochemical_unconditionally():
    """Every protein row used to be written assay_type 'biochemical'
    regardless of what it measured."""
    em = Emit()
    # DISTINCT assay ids: the dedup key is (drug, predicate, target, assay,
    # standard_type), so three rows sharing one assay id collapse to one and
    # the test silently checks a single edge.
    rows = [_viral_row(PSEUDOVIRUS, assay="A1"),
            _viral_row(ELISA, value=7000.0, assay="A2"),
            _viral_row("IC50 against nsp5", value=50.0, assay="A3")]
    ingest_activities(em, rows, {"2697049": {}}, {"2697049": "2697049"},
                      "infores:chembl")
    kinds = {e["qualifiers"]["assay_type"] for e in em.edges}
    assert kinds == {"reporter", "biochemical", "binding"}


# ------------------------------------------- every row's fate is accounted
def _row_covering(kind):
    """One row engineered to take each branch through ingest_activities."""
    base = _viral_row("Inhibition of recombinant 3CL protease", assay=f"A_{kind}")
    if kind == "no_inchikey":
        base["inchikey"] = None
    elif kind == "unconvertible_units":
        base["standard_units"] = "mg/ml"
    elif kind == "multicomponent_target":
        base["target_name"] = "Spike glycoprotein/ACE2"
    elif kind == "no_target_node":
        base["accession"] = None
        base["description"] = "no chain named here"
    elif kind == "undecidable":
        base["standard_value"] = 5000.0
        base["standard_relation"] = ">"
    elif kind == "cellular":
        base["description"] = PSEUDOVIRUS
    elif kind == "organism":
        base["target_type"] = "ORGANISM"
    elif kind == "inactive":
        base["standard_value"] = 50_000.0
        base["standard_relation"] = ">"
    return base


def test_every_row_is_either_emitted_or_dropped():
    """THE INVARIANT. It failed by exactly 1,238 rows when a new counter was
    added and not listed, the layer was never written, and the reassembly
    that followed silently used the previous chembl layer -- so every count
    came back identical and nothing downstream said why.
    """
    kinds = ["active", "no_inchikey", "unconvertible_units",
             "multicomponent_target", "no_target_node", "undecidable",
             "cellular", "organism", "inactive"]
    rows = [_row_covering(k) for k in kinds]
    rows.append(_row_covering("active"))        # a duplicate of the first
    em = Emit()
    s = ingest_activities(em, rows, {"2697049": {}}, {"2697049": "2697049"},
                          "infores:chembl")
    emitted, dropped, unlisted = accounting(s)
    assert emitted + dropped == s["rows"] == len(rows), (
        f"{s['rows']} rows, {emitted} emitted, {dropped} dropped; "
        f"unlisted: {unlisted}")


def test_no_counter_the_module_sets_has_an_unstated_fate():
    """A counter in neither list and matching no breakdown prefix is a row
    whose fate nobody can state."""
    rows = [_row_covering(k) for k in
            ("active", "cellular", "organism", "inactive", "undecidable",
             "no_inchikey", "multicomponent_target")]
    em = Emit()
    s = ingest_activities(em, rows, {"2697049": {}}, {"2697049": "2697049"},
                          "infores:chembl")
    assert accounting(s)[2] == []


def test_the_breakdown_counters_are_excluded_from_the_sum():
    """readout_* and assay_type_* are incremented before dedup and classify,
    so they are a SUPERSET of what gets emitted -- 1,921 rows read as
    cellular while 1,238 were reattributed. Adding them would double-count."""
    em = Emit()
    rows = [_row_covering("cellular"), _row_covering("active")]
    s = ingest_activities(em, rows, {"2697049": {}}, {"2697049": "2697049"},
                          "infores:chembl")
    assert s["readout_cells"] >= 1 and s["readout_protein"] >= 1
    emitted, dropped, _ = accounting(s)
    assert emitted + dropped == len(rows)
