"""Day 6 gate: polarity is never guessed, and a flipped sign is rejected.

This is the one layer where an error produces a confidently WRONG
recommendation rather than a missing one, so the tests are about refusal as
much as about correctness.
"""
import pytest

from kgav.emit import Emit
from kgav.orcs import (
    DEPENDENCY_ANCHORS,
    RESTRICTION_SENTINELS,
    calibrate_sign,
    gene_direction,
    ingest_screens,
    is_virus_resistance_screen,
    polarity_mode,
    screen_date,
    screen_type,
    virus_taxon,
)
from kgav.schema import load_schema


def _screen(sid="1", rationale="Increased resistance to virus", score_type="Z-score",
            condition="Virus: SARS-CoV-2 (isolate USA-WA1/2020)", **over):
    s = {"SCREEN_ID": sid, "SOURCE_ID": "32838362", "AUTHOR": "Hoffmann HH (2020)",
         "SCREEN_RATIONALE": rationale, "SCORE.1_TYPE": score_type,
         "CONDITION_NAME": condition, "CELL_LINE": "Huh-7.5",
         "METHODOLOGY": "Knockout", "PHENOTYPE": "", "NOTES": "", "SCREEN_NAME": ""}
    s.update(over)
    return s


def _gene(symbol, gid, score, hit="YES"):
    return {"OFFICIAL_SYMBOL": symbol, "IDENTIFIER_ID": gid,
            "SCORE.1": str(score), "HIT": hit}


# ------------------------------------------------------------ virus matching
@pytest.mark.parametrize("text,expected", [
    ("Virus: SARS-CoV-2 (isolate USA-WA1/2020)", "2697049"),
    ("Virus: SARS-CoV-2/SB3-TYAGNC", "2697049"),
    ("Virus: HCoV-OC43 (Human coronavirus OC43)", "31631"),
    ("Virus: HCoV-229E", "11137"),
    ("Middle East respiratory syndrome coronavirus", "1335626"),
    ("Influenza A virus", None),
])
def test_virus_matching(text, expected):
    assert virus_taxon({"CONDITION_NAME": text}) == expected


def test_sars_cov_2_is_not_filed_under_sars_cov():
    """'SARS-CoV' is a substring of 'SARS-CoV-2'. Pattern order is the only
    thing stopping every SARS-CoV-2 screen landing on the wrong virus."""
    assert virus_taxon({"CONDITION_NAME": "SARS-CoV-2 Wuhan-Hu-1"}) == "2697049"
    assert virus_taxon({"CONDITION_NAME": "SARS-CoV Urbani"}) == "694009"


# ---------------------------------------------------------------- polarity
def test_one_sided_rationale_fixes_direction():
    assert polarity_mode(_screen(rationale="Increased resistance to virus")) == \
        ("screen_level", "dependency")
    assert polarity_mode(_screen(rationale="Decreased resistance to virus")) == \
        ("screen_level", "restriction")


def test_two_sided_rationale_with_signed_score_is_per_gene():
    assert polarity_mode(_screen(rationale="Increased/decreased resistance to virus",
                                 score_type="Z-score")) == ("signed", None)


def test_two_sided_rationale_with_unsigned_score_is_unresolvable():
    """MAGeCK pos score is one-sided; its sign carries no direction."""
    assert polarity_mode(_screen(rationale="Increased/decreased resistance to virus",
                                 score_type="MAGeCK pos score")) == ("unresolvable", None)


def test_sign_convention_verified_against_screen_1379():
    """Hoffmann 2020, Huh-7.5, OC43. Positive Z gave RAB7A/VPS11/ATP6AP1 --
    endosomal acidification, which OC43 entry requires. Negative Z gave TBK1."""
    assert gene_direction("11.31", None) == "dependency"    # RAB7A
    assert gene_direction("-3.33", None) == "restriction"   # TBK1
    # with an inverted screen convention the same values flip
    assert gene_direction("11.31", None, -1) == "restriction"
    assert gene_direction("-3.33", None, -1) == "dependency"
    assert gene_direction("0", None) is None
    assert gene_direction("not a number", None) is None


def test_fixed_direction_overrides_sign():
    assert gene_direction("-5.0", "dependency") == "dependency"


# -------------------------------------------------------------- exclusions
@pytest.mark.parametrize("rationale,keep", [
    ("Increased resistance to virus", True),
    ("Increased/decreased resistance to virus", True),
    ("Cell-essential genes", False),
    ("Identifying cell surface proteins binding the SARS-CoV-2 spike protein", False),
    ("Increased/decreased viral programmed -1 ribosomal frameshifting", False),
])
def test_non_resistance_screens_are_excluded(rationale, keep):
    """A cell-essential gene is not a host factor. Ingesting those would put
    ribosomal proteins on every path."""
    assert is_virus_resistance_screen({"SCREEN_RATIONALE": rationale}) is keep


def test_screen_type_mapping():
    assert screen_type({"METHODOLOGY": "Knockout"}) == "CRISPRko"
    assert screen_type({"METHODOLOGY": "Activation"}) == "CRISPRa"
    assert screen_type({"METHODOLOGY": "something else"}) == "other"


def test_screen_date_from_author_field():
    assert screen_date({"AUTHOR": "Hoffmann HH (2020)"}) == "2020-01-01"
    # No plausible-looking guess: 1970 sorts before any split point and is
    # findable by query.
    assert screen_date({"AUTHOR": "no year here"}) == "1970-01-01"


# ------------------------------------------------------- SIGN CALIBRATION
def test_calibration_detects_positive_means_dependency():
    """Screen 1379 (Hoffmann 2020): RAB7A +11.31, TBK1 -3.33."""
    mult, seen = calibrate_sign([_gene("RAB7A", "7879", 11.31),
                                 _gene("TBK1", "29110", -3.33)])
    assert mult == 1 and set(seen) == {"RAB7A", "TBK1"}


def test_calibration_detects_positive_means_restriction():
    """Screen 2361 (Le Pen 2024) states z-score >= 2 means antiviral."""
    mult, _ = calibrate_sign([_gene("TBK1", "29110", 4.2),
                              _gene("RAB7A", "7879", -8.0)])
    assert mult == -1


def test_calibration_refuses_when_anchors_split_evenly():
    """A coin flip here inverts the biology for every gene in the screen."""
    mult, seen = calibrate_sign([_gene("TBK1", "29110", 4.2),
                                 _gene("STAT1", "6772", -3.0)])
    assert mult == 0 and len(seen) == 2


def test_calibration_refuses_with_no_anchors():
    assert calibrate_sign([_gene("SOMEGENE", "1", 5.0)]) == (0, [])


def test_calibration_tolerates_a_single_dissenter():
    """Real screens are noisy; demanding unanimity would reject good data."""
    mult, _ = calibrate_sign([_gene("RAB7A", "7879", 9.0), _gene("NPC1", "4864", 7.0),
                              _gene("TBK1", "29110", 2.0)])
    assert mult == 1


def test_uncalibratable_screen_is_skipped_not_guessed():
    em = Emit()
    rows = {"1": [_gene("SOMEGENE", "7879", 8.0)]}
    stats = ingest_screens(em, [_screen(rationale="Increased/decreased resistance to virus")],
                           rows, {"NCBIGene:7879": ["UniProtKB:P51149"]},
                           "infores:biogrid-orcs")
    assert stats["uncalibrated_sign"] == 1 and stats["edges"] == 0


def test_screen_level_rationale_contradicted_by_anchors_is_rejected():
    """A one-sided rationale describes the screen's primary selection, not
    every hit's polarity. Screen 1704 says 'increased resistance' yet contains
    IFNAR1, STAT1 and IRF9."""
    em = Emit()
    rows = {"1": [_gene("IFNAR1", "3454", 5.0), _gene("STAT1", "6772", 4.0)]}
    stats = ingest_screens(em, [_screen(rationale="Increased resistance to virus")],
                           rows, {"NCBIGene:3454": ["UniProtKB:P17181"],
                                  "NCBIGene:6772": ["UniProtKB:P42224"]},
                           "infores:biogrid-orcs")
    assert stats["rationale_contradicted_by_anchors"] == 1 and stats["edges"] == 0


def test_anchor_sets_cover_both_poles():
    assert {"TBK1", "IFNAR1", "STAT1", "MAVS", "JAK1"} <= RESTRICTION_SENTINELS
    assert {"ACE2", "TMPRSS2", "CTSL", "RAB7A", "NPC1"} <= DEPENDENCY_ANCHORS


# ------------------------------------------------------------------ ingest
@pytest.fixture
def ingested():
    em = Emit()
    rows = {"1": [_gene("RAB7A", "7879", 11.31), _gene("TBK1", "29110", -3.33),
                  _gene("AAMP", "14", -0.23, hit="NO"),
                  _gene("UNKNOWN", "99999999", 5.0)]}
    stats = ingest_screens(
        em, [_screen(sid="1", rationale="Increased/decreased resistance to virus")],
        rows,
        {"NCBIGene:7879": ["UniProtKB:P51149"], "NCBIGene:29110": ["UniProtKB:Q9UHD2"]},
        "infores:biogrid-orcs")
    return em, stats


def test_both_directions_are_emitted(ingested):
    em, stats = ingested
    assert stats["dependency"] == 1 and stats["restriction"] == 1
    dirs = {e["qualifiers"]["direction"] for e in em.edges}
    assert dirs == {"dependency", "restriction"}


def test_non_hits_and_unmapped_genes_are_skipped(ingested):
    em, stats = ingested
    assert stats["gene_not_in_graph"] == 1
    assert len(em.edges) == 2


def test_required_qualifiers_present(ingested):
    em, _ = ingested
    for e in em.edges:
        q = e["qualifiers"]
        assert q["direction"] and q["screen_type"] and q["cell_line"]
        assert e["publications"][0].startswith("PMID:")
        assert e["first_asserted_date"] == "2020-01-01"


def test_output_validates_against_schema(ingested):
    em, _ = ingested
    nodes = [
        {"id": "UniProtKB:P51149", "class": "Protein",
         "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False,
                        "sequence_hash": "a", "reviewed": True}},
        {"id": "UniProtKB:Q9UHD2", "class": "Protein",
         "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False,
                        "sequence_hash": "b", "reviewed": True}},
        {"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
    ]
    violations = load_schema().validate_batch(nodes, em.edges)
    assert violations == [], [str(v) for v in violations[:6]]
