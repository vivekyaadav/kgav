"""The panel layer's contract: a negative must survive the round trip.

This layer exists to supply measured inactives in the approved-drug
population. An inactive that labels.classify() cannot read back as inactive
would be a label the evaluation silently discards, so the round trip is tested
directly rather than inferred from the qualifiers.
"""
from __future__ import annotations

import pytest

from kgav.emit import Emit
from kgav.labels import classify
from kgav.ncats import (
    ANTIVIRAL_PREDICATE,
    CPE_PUBLICATION,
    SARS_COV_2,
    ColumnError,
    classify_row,
    ingest_cpe,
    join_report,
    resolve_columns,
    to_nm,
    unit_for,
)
from kgav.schema import load_schema


class FakeId:
    """Stands in for normalize_chemical's result so these tests need no rdkit."""

    def __init__(self, structure: str) -> None:
        if structure == "BAD":
            raise ValueError("unparseable")
        self.inchikey = "".join(c for c in structure if c.isalnum()).upper().ljust(27, "X")[:27]
        self.curie = f"INCHIKEY:{self.inchikey}"


def norm(s):
    return FakeId(s)


COLS = {"smiles": "SMILES", "activity": "ACTIVITY", "ac50": "AC50 (uM)",
        "cc50": "CC50 (uM)", "sample": "NCGC SID", "efficacy": "EFFICACY"}


def _row(smiles="CCO", activity="", ac50="", cc50="", eff="", sid="NCGC1"):
    return {"SMILES": smiles, "ACTIVITY": activity, "AC50 (uM)": ac50,
            "CC50 (uM)": cc50, "EFFICACY": eff, "NCGC SID": sid}


# ------------------------------------------------------------------- units
@pytest.mark.parametrize("value,unit,nm", [
    (1.0, "uM", 1_000.0),
    (1.0, "nM", 1.0),
    (1e-6, "M", 1_000.0),
    (-5.0, "logM", 10_000.0),     # NCATS LAC50: log10(M); -5 is 10 uM
    (-6.0, "logM", 1_000.0),
    ("", "uM", None),
    (None, "uM", None),
    ("not a number", "uM", None),
])
def test_concentration_conversion(value, unit, nm):
    got = to_nm(value, unit)
    if nm is None:
        assert got is None
    else:
        assert got == pytest.approx(nm)


def test_unknown_unit_raises_rather_than_guessing():
    with pytest.raises(ValueError):
        to_nm(1.0, "molar")


@pytest.mark.parametrize("col,unit", [
    ("LAC50", "logM"), ("log AC50 (M)", "logM"), ("AC50 (uM)", "uM"),
    ("AC50_nM", "nM"), ("AC50 (M)", "M"), ("AC50", "uM"),
])
def test_unit_inferred_from_column_name(col, unit):
    assert unit_for(col) == unit


# ----------------------------------------------------------------- columns
def test_header_variants_resolve():
    cols = resolve_columns(["NCGC SID", "SMILES", "Activity", "LAC50",
                            "Efficacy", "Curve Class", "Cytotox AC50"])
    assert cols["smiles"] == "SMILES"
    assert cols["lac50"] == "LAC50"
    assert cols["cc50"] == "Cytotox AC50"


def test_missing_structure_column_is_fatal_and_shows_the_header():
    """A partial match is a different file, not a degraded mode."""
    with pytest.raises(ColumnError) as e:
        resolve_columns(["compound_id", "AC50", "Activity"])
    assert "compound_id" in str(e.value)      # the real header is in the error
    assert "COLUMN_CANDIDATES" in str(e.value)  # and how to fix it


# ------------------------------------------------------- the round trip
def test_an_active_round_trips_through_labels_classify():
    v, q, why = classify_row(_row(ac50="2.5", cc50="50"), COLS,
                             ac50_units="uM", cc50_units="uM", max_conc_nm=46_000)
    assert (v, why) == ("active", "fitted_ac50")
    assert q["ec50_nm"] == pytest.approx(2_500.0)
    assert q["relation"] == "="
    edge = {"predicate": ANTIVIRAL_PREDICATE, "qualifiers": q}
    assert classify(edge) == "active"


def test_an_inactive_round_trips_through_labels_classify():
    """THE POINT OF THE LAYER. A negative the evaluation cannot read back is
    not a negative."""
    v, q, why = classify_row(_row(activity="inactive"), COLS,
                             ac50_units="uM", cc50_units="uM", max_conc_nm=46_000)
    assert (v, why) == ("inactive", "called_inactive")
    assert q["relation"] == ">" and q["ec50_nm"] == 46_000
    edge = {"predicate": ANTIVIRAL_PREDICATE, "qualifiers": q}
    assert classify(edge) == "inactive"


def test_inactive_below_the_pipeline_threshold_is_refused():
    """Tested to 5 uM and inactive says nothing about 10 uM. labels.classify()
    would return None for such an edge, so counting it as a negative here
    would invent evidence the assay does not contain."""
    v, _, why = classify_row(_row(activity="inactive"), COLS, ac50_units="uM",
                             cc50_units="uM", max_conc_nm=5_000)
    assert v is None and why == "below_threshold"


def test_a_fitted_ac50_outranks_the_outcome_column():
    """The outcome column is the depositor's curve-class heuristic; the AC50
    is the fit. A row called inactive that nonetheless has a potency is
    recorded by its measurement."""
    v, q, _ = classify_row(_row(activity="inactive", ac50="0.8"), COLS,
                           ac50_units="uM", cc50_units="uM", max_conc_nm=46_000)
    assert v == "active" and q["ec50_nm"] == pytest.approx(800.0)


def test_called_active_without_a_value_is_emitted_unquantified():
    """`unquantified` is derived at assembly, so this must NOT invent a
    potency -- the H2 bug was a layer asserting what only the merge can know."""
    v, q, why = classify_row(_row(activity="active"), COLS, ac50_units="uM",
                             cc50_units="uM", max_conc_nm=46_000)
    assert (v, why) == ("active", "called_active_unquantified")
    assert "ec50_nm" not in q and "unquantified" not in q


def test_no_call_and_no_value_yields_nothing():
    v, _, why = classify_row(_row(), COLS, ac50_units="uM", cc50_units="uM",
                             max_conc_nm=46_000)
    assert v is None and why == "no_call"


# ------------------------------------------------------------- selectivity
def test_same_plate_selectivity_index_is_marked_verified():
    _, q, _ = classify_row(_row(ac50="1.0", cc50="30"), COLS, ac50_units="uM",
                           cc50_units="uM", max_conc_nm=46_000)
    assert q["selectivity_index"] == pytest.approx(30.0)
    assert q["selectivity_verified"] is True
    assert q["selectivity_n_studies"] == 1


def test_no_selectivity_index_without_a_cc50():
    _, q, _ = classify_row(_row(ac50="1.0"), COLS, ac50_units="uM",
                           cc50_units="uM", max_conc_nm=46_000)
    assert "selectivity_index" not in q and "selectivity_verified" not in q


# ----------------------------------------------------------------- ingest
@pytest.fixture
def ingested():
    em = Emit()
    rows = [
        _row(smiles="CCO", ac50="1.0", cc50="40", sid="S1"),
        _row(smiles="CCC", activity="inactive", sid="S2"),
        _row(smiles="CCO", ac50="2.0", sid="S3"),        # same structure again
        _row(smiles="", activity="inactive", sid="S4"),  # no structure
        _row(smiles="BAD", activity="inactive", sid="S5"),  # unnormalizable
        _row(smiles="CCCC", sid="S6"),                   # no call
    ]
    # CCC is already in the graph; CCO is not. Both branches of the
    # is_approved decision get exercised by the same fixture.
    stats = ingest_cpe(em, rows, COLS, norm, max_conc_nm=46_000,
                       existing_compounds={FakeId("CCC").curie})
    return em, stats


def test_both_classes_are_emitted(ingested):
    em, stats = ingested
    assert stats["active"] == 1 and stats["inactive"] == 1
    assert len(em.edges) == 2


def test_one_structure_votes_once(ingested):
    """Two plate records for one structure would double its DWPC weight and
    let one compound vote twice in the AUC."""
    _, stats = ingested
    assert stats["duplicate_structure"] == 1


def test_unusable_rows_are_counted_not_silently_dropped(ingested):
    _, stats = ingested
    assert stats["no_structure"] == 1
    assert stats["unnormalizable"] == 1
    assert stats["reason_no_call"] == 1
    assert stats["rows"] == 6


def test_every_edge_cites_the_screen_and_is_dated(ingested):
    em, _ = ingested
    for e in em.edges:
        assert e["publications"] == [CPE_PUBLICATION]
        assert e["object"] == SARS_COV_2
        assert e["predicate"] == ANTIVIRAL_PREDICATE
        assert e["first_asserted_date"].startswith("2020")
        assert e["qualifiers"]["assay_type"] == "cell_based_antiviral"


def test_identifiers_come_from_the_normalizer_not_the_file(ingested):
    """The dominant bug class in this repo is a discarded identifier. A key
    taken from the file would differ on salt, stereo or tautomer and join
    nothing."""
    em, _ = ingested
    for nid in em.nodes:
        assert nid.startswith("INCHIKEY:")
    for e in em.edges:
        assert e["subject"].startswith("INCHIKEY:")
        assert not e["subject"].startswith("INCHIKEY:NCGC")


def test_a_compound_the_graph_already_holds_gets_no_node(ingested):
    """The H2 lesson: this layer cannot know is_approved, so it must not
    overwrite the layer that does. Only the label edge is contributed."""
    em, stats = ingested
    assert stats["joined_existing_compound"] == 1
    assert stats["new_compound_node"] == 1
    assert FakeId("CCC").curie not in em.nodes        # joined -> no node
    assert FakeId("CCO").curie in em.nodes            # new -> node
    # but BOTH got their label edge
    assert {e["subject"] for e in em.edges} == {FakeId("CCO").curie,
                                                FakeId("CCC").curie}


def test_a_new_compound_gets_the_floor_not_a_claim(ingested):
    em, _ = ingested
    n = em.nodes[FakeId("CCO").curie]
    assert n["properties"]["is_approved"] is False
    assert "max_phase" not in n["properties"]   # not invented either


def test_output_validates_against_schema(ingested):
    em, _ = ingested
    nodes = list(em.nodes.values()) + [
        {"id": FakeId("CCC").curie, "class": "SmallMolecule",
         "properties": {"smiles": "CCC", "inchikey_skel": "CCCXXXXXXXXXXX",
                        "is_approved": True, "salt_collapsed": False,
                        "stereo_collapsed": False}},
    ] + [
        {"id": SARS_COV_2, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}}]
    violations = load_schema().validate_batch(nodes, em.edges)
    assert violations == [], [str(v) for v in violations[:6]]


# ------------------------------------------------------------- join report
def test_join_report_counts_what_actually_connects():
    r = join_report({"INCHIKEY:A", "INCHIKEY:B", "INCHIKEY:C"},
                    {"INCHIKEY:A", "INCHIKEY:Z"})
    assert r == {"emitted": 3, "already_in_graph": 1, "new_compounds": 2,
                 "join_rate": pytest.approx(1 / 3)}


def test_join_report_on_an_empty_layer_does_not_divide_by_zero():
    assert join_report(set(), {"INCHIKEY:A"})["join_rate"] == 0.0


# ------------------------------------------------- interaction with the filter
def test_panel_labels_leak_nothing_through_the_publication_filter(tmp_path):
    """The reason this panel was chosen. Its labels all cite one screen, and
    no ChEMBL binding edge cites that screen -- so --publication-disjoint
    withholds nothing when scoring against them. If a future layer did cite
    PMID:33708112 on an INHIBITS edge, this test is where it would surface."""
    import json

    from kgav.baselines import label_publication_index

    em = Emit()
    ingest_cpe(em, [_row(smiles="CCO", ac50="1.0")], COLS, norm,
               max_conc_nm=46_000)
    # a ChEMBL-style binding edge for the same compound, different paper
    other = {"subject": "INCHIKEY:CCOXXXXXXXXXXXXXXXXXXXXXXX",
             "predicate": "INHIBITS", "object": "UniProtKB:X",
             "qualifiers": {}, "primary_knowledge_source": "infores:chembl",
             "evidence_tier": 1, "first_asserted_date": "2020-01-01",
             "publications": ["PMID:99999"]}
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in em.edges + [other]))
    idx = label_publication_index(tmp_path, {ANTIVIRAL_PREDICATE})
    assert list(idx.values()) == [{CPE_PUBLICATION}]
    # the binding edge's paper is NOT in the label index, so it survives
    assert CPE_PUBLICATION != "PMID:99999"
    assert "PMID:99999" not in next(iter(idx.values()))
