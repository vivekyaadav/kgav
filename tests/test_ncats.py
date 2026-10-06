"""The panel layer's two contracts.

1. A negative must survive the round trip. This layer exists to supply
   measured inactives in the approved-drug population; an inactive that
   labels.classify() cannot read back as inactive is a label the evaluation
   silently discards.

2. Direction must come from a measurement, never from a fitted value. CPE is
   gain-of-signal: a protective compound raises the readout and a cytotoxic one
   lowers it, and BOTH produce a clean AC50. The real export's first row is
   ac50 7.08 uM with efficacy -46.9. Filing that as an antiviral is the
   BioGRID ORCS sign error in a new costume, so it gets its own test.
"""
from __future__ import annotations

import json

import pytest

from kgav.baselines import label_publication_index
from kgav.emit import Emit
from kgav.labels import classify
from kgav.ncats import (
    ANTIVIRAL_PREDICATE,
    CANONICAL_ASSAY,
    COUNTERSCREEN_ASSAY,
    CPE_PUBLICATION,
    SARS_COV_2,
    ColumnError,
    classify_row,
    curve_class,
    index_counterscreen,
    ingest_cpe,
    join_report,
    resolve_columns,
    to_nm,
    top_concentration_nm,
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


# The real assay/ export's header, lowercased as NCATS ships it.
HEADER = ["sample_id", "sample_name", "pubchem_sid", "primary_moa", "assay_name",
          "library", "cas", "alias", "ac50", "log_ac50", "auc", "curve_class2",
          "drug_name", "efficacy", "max_response", "p_hill", "r2", "gene_id",
          "gene_symbol", "smiles", "supplier", "supplier_id",
          "conc0", "conc1", "conc2", "conc3"]
COLS = resolve_columns(HEADER)


def _row(smiles="CCO", sid="NCGC1", assay=CANONICAL_ASSAY, ac50="",
         curve="4", eff="0", concs=("1.60E-07", "8.00E-07", "4.00E-06", "2.00E-05")):
    r = {h: "" for h in HEADER}
    r.update({"sample_id": sid, "smiles": smiles, "assay_name": assay,
              "ac50": ac50, "curve_class2": curve, "efficacy": eff,
              "library": "Approved Drugs Collection (NPC)"})
    for i, c in enumerate(concs):
        r[f"conc{i}"] = c
    return r


def _tox(sid="NCGC1", ac50="", curve="4", eff="0"):
    r = _row(sid=sid, assay=COUNTERSCREEN_ASSAY, ac50=ac50, curve=curve, eff=eff)
    return r


# ------------------------------------------------------------------- columns
def test_the_real_header_resolves():
    assert COLS["smiles"] == "smiles"
    assert COLS["ac50"] == "ac50"
    assert COLS["curve"] == "curve_class2"
    assert COLS["efficacy"] == "efficacy"
    assert COLS["sample"] == "sample_id"        # the stable id, not sample_name
    assert COLS["conc_cols"] == ["conc0", "conc1", "conc2", "conc3"]


def test_uppercase_tsv_header_also_resolves():
    cols = resolve_columns(["SAMPLE_ID", "SMILES", "AC50(uM)", "LOG_AC50",
                            "CURVE_CLASS2", "EFFICACY", "CONC0(M)", "CONC1(M)"])
    assert cols["ac50"] == "AC50(uM)" and cols["curve"] == "CURVE_CLASS2"
    assert cols["conc_cols"] == ["CONC0(M)", "CONC1(M)"]


@pytest.mark.parametrize("missing", ["smiles", "curve_class2", "efficacy"])
def test_each_required_column_is_fatal_when_absent(missing):
    """Without curve class there is no activity call; without efficacy there is
    no direction, and a guessed direction files cytotoxic compounds as
    antivirals."""
    with pytest.raises(ColumnError) as e:
        resolve_columns([h for h in HEADER if h != missing])
    assert "COLUMN_CANDIDATES" in str(e.value)


# --------------------------------------------------------------------- units
@pytest.mark.parametrize("value,unit,nm", [
    (1.0, "uM", 1_000.0), (1.0, "nM", 1.0), (1e-6, "M", 1_000.0),
    (-5.0, "logM", 10_000.0), (-4.9, "logM", 12_589.25),
    ("", "uM", None), (None, "uM", None), ("x", "uM", None),
])
def test_concentration_conversion(value, unit, nm):
    got = to_nm(value, unit)
    assert (got is None) if nm is None else got == pytest.approx(nm, rel=1e-4)


def test_log_ac50_and_ac50_agree_on_the_real_row():
    """Protriptyline in the real export: ac50 12.5892541179, log_ac50 -4.9.
    If these disagree the unit mapping is wrong."""
    assert to_nm(12.5892541179, "uM") == pytest.approx(to_nm(-4.9, "logM"), rel=1e-3)


def test_unknown_unit_raises_rather_than_guessing():
    with pytest.raises(ValueError):
        to_nm(1.0, "molar")


@pytest.mark.parametrize("col,unit", [
    ("LOG_AC50", "logM"), ("log_ac50", "logM"), ("AC50(uM)", "uM"),
    ("ac50", "uM"), ("AC50_nM", "nM"),
])
def test_unit_inferred_from_column_name(col, unit):
    assert unit_for(col) == unit


# ------------------------------------------------------------ concentrations
def test_top_concentration_is_read_from_the_row():
    """The real primary screen tops out at 2.0E-5 M = 20 uM. A hardcoded
    46 uM would record a censoring point no compound was tested at."""
    assert top_concentration_nm(_row(), COLS) == pytest.approx(20_000.0)


def test_top_concentration_ignores_blank_and_zero_columns():
    r = _row(concs=("1.60E-07", "8.00E-07"))
    assert top_concentration_nm(r, COLS) == pytest.approx(800.0)


def test_no_concentrations_means_none_not_a_default():
    r = _row(concs=())
    assert top_concentration_nm(r, COLS) is None


# ------------------------------------------------------------- curve classes
@pytest.mark.parametrize("raw,mag,sign", [
    ("4", "4", 1), ("-2.4", "2.4", -1), ("1.1", "1.1", 1),
    ("-1.2", "1.2", -1), ("", "", 0), ("5", "5", 1),
])
def test_curve_class_splits_magnitude_from_sign(raw, mag, sign):
    assert curve_class(_row(curve=raw), COLS) == (mag, sign)


# ----------------------------------------------------- THE ROUND TRIP (1)
def test_an_active_round_trips_through_labels_classify():
    v, q, why = classify_row(_row(ac50="2.5", curve="1.1", eff="80"), COLS,
                             ac50_units="uM")
    assert (v, why) == ("active", "fitted_active")
    assert q["ec50_nm"] == pytest.approx(2_500.0) and q["relation"] == "="
    assert classify({"predicate": ANTIVIRAL_PREDICATE, "qualifiers": q}) == "active"


def test_an_inactive_round_trips_through_labels_classify():
    """THE POINT OF THE LAYER. A negative the evaluation cannot read back is
    not a negative."""
    v, q, why = classify_row(_row(curve="4"), COLS, ac50_units="uM")
    assert (v, why) == ("inactive", "curve_class_4")
    assert q["relation"] == ">" and q["ec50_nm"] == pytest.approx(20_000.0)
    assert classify({"predicate": ANTIVIRAL_PREDICATE, "qualifiers": q}) == "inactive"


def test_inactive_below_the_pipeline_threshold_is_refused():
    """Tested to 4 uM and inactive says nothing about 10 uM. classify() would
    return None for such an edge, so counting it would invent evidence."""
    r = _row(curve="4", concs=("1.60E-07", "8.00E-07", "4.00E-06"))
    v, _, why = classify_row(r, COLS, ac50_units="uM")
    assert v is None and why == "below_threshold"


def test_inactive_without_a_concentration_is_refused():
    v, _, why = classify_row(_row(curve="4", concs=()), COLS, ac50_units="uM")
    assert v is None and why == "inactive_without_a_concentration"


# -------------------------------------------------- THE SIGN PROBLEM (2)
def test_negative_efficacy_is_cytotoxicity_not_antiviral_activity():
    """The real export's first CPE row: ac50 7.079457844, curve_class2 '5',
    efficacy -46.92489053. In a GAIN-of-signal assay a negative efficacy means
    the compound reduced viability. An AC50-only rule files it as an antiviral
    with EC50 7 uM -- the same error as a mis-signed restriction factor
    becoming a recommendation to help the virus."""
    r = _row(ac50="7.079457844", curve="2.1", eff="-46.92489053")
    v, _, why = classify_row(r, COLS, ac50_units="uM")
    assert v is None and why == "negative_efficacy_cytotoxic"


def test_a_fitted_curve_below_the_efficacy_floor_is_not_an_active():
    r = _row(ac50="5.0", curve="1.1", eff="12")
    v, _, why = classify_row(r, COLS, ac50_units="uM")
    assert v is None and why == "below_min_efficacy"


def test_the_efficacy_floor_is_configurable():
    r = _row(ac50="5.0", curve="1.1", eff="60")
    assert classify_row(r, COLS, ac50_units="uM", min_efficacy_pct=50)[0] == "active"
    assert classify_row(r, COLS, ac50_units="uM", min_efficacy_pct=70)[0] is None


def test_curve_class_sign_is_recorded_but_does_not_decide():
    """Efficacy decides direction because it is measured. The class sign is
    kept for provenance only -- one inferred polarity convention per project
    is already one too many."""
    r = _row(ac50="5.0", curve="-1.1", eff="80")
    v, q, _ = classify_row(r, COLS, ac50_units="uM")
    assert v == "active"                          # efficacy is positive
    assert q["source_curve_class"] == "-1.1"      # but the sign is on the edge


# ------------------------------------------------- a fit is not a call
def test_single_point_activity_is_not_an_active():
    """586 class-3 rows carry an AC50. Admitting them on the strength of that
    value alone would multiply the active set with unreplicated single
    points."""
    r = _row(ac50="8.0", curve="3", eff="90")
    v, _, why = classify_row(r, COLS, ac50_units="uM")
    assert v is None and why == "curve_class_3_not_a_fit"


def test_a_poorly_fit_curve_is_not_an_active():
    v, _, why = classify_row(_row(ac50="7.1", curve="5", eff="60"), COLS,
                             ac50_units="uM")
    assert v is None and why == "curve_class_5_not_a_fit"


def test_a_blank_curve_class_is_not_an_active():
    v, _, why = classify_row(_row(ac50="7.1", curve="", eff="60"), COLS,
                             ac50_units="uM")
    assert v is None and why == "curve_class_blank_not_a_fit"


def test_a_real_fit_without_a_value_is_emitted_unquantified():
    """`unquantified` is derived at assembly, so this must NOT invent a
    potency -- the H2 bug was a layer asserting what only the merge can know."""
    v, q, why = classify_row(_row(curve="2.2", eff="70"), COLS, ac50_units="uM")
    assert (v, why) == ("active", "fitted_unquantified")
    assert "ec50_nm" not in q and "unquantified" not in q


# -------------------------------------------------- same-plate selectivity
def test_selectivity_index_comes_off_the_counterscreen_plate():
    v, q, _ = classify_row(_row(ac50="1.0", curve="1.1", eff="80"), COLS,
                           ac50_units="uM", tox_row=_tox(ac50="30", curve="-2.1",
                                                         eff="-70"),
                           tox_cols=COLS, tox_ac50_units="uM")
    assert v == "active"
    assert q["cc50_nm"] == pytest.approx(30_000.0)
    assert q["selectivity_index"] == pytest.approx(30.0)
    assert q["selectivity_verified"] is True and q["selectivity_n_studies"] == 1


def test_an_inactive_counterscreen_means_no_index_not_a_floor():
    """Not cytotoxic up to the top concentration is a real result, but CC50 is
    censored. Writing top_conc/ec50 would UNDERSTATE a compound that is
    cleaner than the assay can measure."""
    v, q, _ = classify_row(_row(ac50="1.0", curve="1.1", eff="80"), COLS,
                           ac50_units="uM", tox_row=_tox(curve="4"),
                           tox_cols=COLS, tox_ac50_units="uM")
    assert v == "active"
    assert "selectivity_index" not in q
    assert q["cytotoxicity_not_detected"] is True


def test_no_counterscreen_row_means_no_index():
    _, q, _ = classify_row(_row(ac50="1.0", curve="1.1", eff="80"), COLS,
                           ac50_units="uM")
    assert "selectivity_index" not in q and "cytotoxicity_not_detected" not in q


def test_counterscreen_index_keys_on_sample_and_filters_by_assay():
    rows = [_tox(sid="A", ac50="5", curve="1.1"),
            _row(sid="B", assay=CANONICAL_ASSAY)]   # wrong assay for this index
    idx = index_counterscreen(rows, COLS, COUNTERSCREEN_ASSAY)
    assert set(idx) == {"A"}


# ------------------------------------------------------------------ ingest
@pytest.fixture
def ingested():
    em = Emit()
    rows = [
        _row(smiles="CCO", sid="S1", ac50="1.0", curve="1.1", eff="80"),
        _row(smiles="CCC", sid="S2", curve="4"),                  # inactive
        _row(smiles="CCCO", sid="S3", ac50="7.1", curve="2.1", eff="-46"),  # tox
        _row(smiles="CCCC", sid="S4", ac50="8.0", curve="3", eff="90"),     # 1pt
        _row(smiles="CCO", sid="S5", ac50="2.0", curve="1.2", eff="75"),    # dup
        _row(smiles="", sid="S6", curve="4"),                     # no structure
        _row(smiles="BAD", sid="S7", curve="4"),                  # unparseable
        _row(smiles="CCCCC", sid="S8", assay="SARS-CoV-2_CPE_SRI", curve="4"),
    ]
    tox = {"S1": _tox(sid="S1", ac50="40", curve="-1.1", eff="-80")}
    stats = ingest_cpe(em, rows, COLS, norm, tox_by_sample=tox, tox_cols=COLS,
                       existing_compounds={FakeId("CCC").curie})
    return em, stats


def test_both_classes_are_emitted(ingested):
    em, stats = ingested
    assert stats["active"] == 1 and stats["inactive"] == 1
    assert len(em.edges) == 2


def test_the_cytotoxic_row_is_not_an_active(ingested):
    _, stats = ingested
    assert stats["reason_negative_efficacy_cytotoxic"] == 1
    assert stats["active"] == 1          # only CCO, not CCCO


def test_a_different_assay_is_not_merged_in(ingested):
    """cpe.tsv pools three protocols with different concentration ranges;
    merging them averages two censoring points into one."""
    _, stats = ingested
    assert stats["other_assay"] == 1
    assert stats["in_assay"] == 7


def test_one_structure_votes_once(ingested):
    _, stats = ingested
    assert stats["duplicate_structure"] == 1


def test_unusable_rows_are_counted_not_silently_dropped(ingested):
    _, stats = ingested
    assert stats["no_structure"] == 1
    assert stats["unnormalizable"] == 1
    assert stats["reason_curve_class_3_not_a_fit"] == 1
    assert stats["rows"] == 8


def test_an_active_on_any_plate_beats_an_inactive_on_another():
    """A compound that worked once is not evidence of inactivity."""
    em = Emit()
    rows = [_row(smiles="CCO", sid="A", curve="4"),
            _row(smiles="CCO", sid="B", ac50="1.0", curve="1.1", eff="80")]
    stats = ingest_cpe(em, rows, COLS, norm)
    assert stats["upgraded_inactive_to_active"] == 1
    assert len(em.edges) == 1
    assert em.edges[0]["qualifiers"]["relation"] == "="


def test_same_plate_selectivity_is_counted(ingested):
    _, stats = ingested
    assert stats["with_same_plate_si"] == 1
    assert stats["paired_with_counterscreen"] == 1


def test_a_compound_the_graph_already_holds_gets_no_node(ingested):
    """The H2 lesson: this layer cannot know is_approved, so it must not
    overwrite the layer that does. Only the label edge is contributed."""
    em, stats = ingested
    assert stats["joined_existing_compound"] == 1
    assert stats["new_compound_node"] == 1
    assert FakeId("CCC").curie not in em.nodes
    assert FakeId("CCO").curie in em.nodes
    assert {e["subject"] for e in em.edges} == {FakeId("CCO").curie,
                                               FakeId("CCC").curie}


def test_a_new_compound_gets_the_floor_not_a_claim(ingested):
    em, _ = ingested
    props = em.nodes[FakeId("CCO").curie]["properties"]
    assert props["is_approved"] is False
    assert "max_phase" not in props


def test_every_edge_cites_the_screen_and_is_dated(ingested):
    em, _ = ingested
    for e in em.edges:
        assert e["publications"] == [CPE_PUBLICATION]
        assert e["object"] == SARS_COV_2
        assert e["predicate"] == ANTIVIRAL_PREDICATE
        assert e["first_asserted_date"].startswith("2020")
        assert e["qualifiers"]["assay_type"] == "cell_based_antiviral"
        assert e["qualifiers"]["cell_line"] == "CVCL:0574"


def test_identifiers_come_from_the_normalizer_not_the_file(ingested):
    """The dominant bug class in this repo is a discarded identifier. A key
    taken from the file would differ on salt, stereo or tautomer and join
    nothing."""
    em, _ = ingested
    for e in em.edges:
        assert e["subject"].startswith("INCHIKEY:")
        assert "NCGC" not in e["subject"]


def test_output_validates_against_schema(ingested):
    em, _ = ingested
    nodes = list(em.nodes.values()) + [
        {"id": FakeId("CCC").curie, "class": "SmallMolecule",
         "properties": {"smiles": "CCC", "inchikey_skel": "CCCXXXXXXXXXXX",
                        "is_approved": True, "salt_collapsed": False,
                        "stereo_collapsed": False}},
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


# ------------------------------------------- interaction with the filter
def test_panel_labels_leak_nothing_through_the_publication_filter(tmp_path):
    """The reason this panel was chosen. Its labels all cite one screen, and no
    ChEMBL binding edge cites that screen -- so --publication-disjoint
    withholds nothing when scoring against them. If a future layer did cite
    PMID:33708112 on an INHIBITS edge, this test is where it would surface."""
    em = Emit()
    ingest_cpe(em, [_row(smiles="CCO", ac50="1.0", curve="1.1", eff="80")],
               COLS, norm)
    binding = {"subject": FakeId("CCO").curie, "predicate": "INHIBITS",
               "object": "UniProtKB:X", "qualifiers": {},
               "primary_knowledge_source": "infores:chembl", "evidence_tier": 1,
               "first_asserted_date": "2020-01-01",
               "publications": ["PMID:99999"]}
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in em.edges + [binding]))
    idx = label_publication_index(tmp_path, {ANTIVIRAL_PREDICATE})
    assert list(idx.values()) == [{CPE_PUBLICATION}]
    assert "PMID:99999" not in next(iter(idx.values()))


def test_the_layer_cannot_self_validate_and_that_is_by_design(ingested):
    """Its own node set omits every compound that joined and the OrganismTaxon
    the spine owns, because it deliberately does not re-emit nodes it does not
    own. Validating it alone reports DANGLING for every edge -- a true
    statement about an incomplete node set, not about the edges. The driver
    must validate against release_nodes + em.nodes, as ingest_chembl.py does.

    If this test ever passes, the layer has started asserting nodes it has no
    authority over, which is the H2 bug.
    """
    em, _ = ingested
    alone = load_schema().validate_batch(list(em.nodes.values()), em.edges)
    assert alone, "the layer validated in isolation -- is it emitting the taxon?"
    assert all(v.code == "DANGLING" for v in alone), \
        sorted({v.code for v in alone})
