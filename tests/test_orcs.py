"""Day 6 gate: polarity is never guessed, and a flipped sign is rejected.

This is the one layer where an error produces a confidently WRONG
recommendation rather than a missing one, so the tests are about refusal as
much as about correctness.
"""
from pathlib import Path

import pytest

from kgav.emit import Emit
from kgav.orcs import (
    DEPENDENCY_ANCHORS,
    RESTRICTION_SENTINELS,
    calibrate_sign,
    gene_direction,
    index_coverage,
    ingest_screens,
    is_virus_resistance_screen,
    per_taxon,
    polarity_mode,
    screen_date,
    screen_heading,
    screen_type,
    unbuilt_taxa,
    virus_taxon,
)
from kgav.schema import load_schema

ROOT = Path(__file__).resolve().parents[1]


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
    # The exact CONDITION_NAME strings the real ORCS index carries for its
    # 11 flavivirus screens. Before these patterns existed virus_taxon
    # returned None for every one, so the screens were present and
    # unreachable.
    ("Virus: Dengue virus 2 Thailand/16681/84", "11060"),
    ("Virus: Dengue virus 2 Jamaica/1409/1983", "11060"),
    ("Virus: West Nile virus", "11082"),
    ("Virus: Zika virus", "64320"),
    ("Virus: Yellow Fever Virus 17D", "11089"),
    # No screen exists for these in this ORCS release, so no pattern does
    # either. None is the correct answer, not a gap to paper over.
    ("Virus: Dengue virus 1", None),
    ("Virus: Dengue virus 3", None),
    ("Virus: Japanese encephalitis virus", None),
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


# --------------------------------------------------------------------------
# The polarity override must be able to take effect.
# --------------------------------------------------------------------------
def test_an_overridden_screen_is_still_fetched():
    """The override existed and could never apply.

    ingest_orcs built `fetchable` from polarity_mode alone, so a screen whose
    polarity comes from config/orcs_polarity_overrides.json was never fetched.
    ingest_screens is handed `usable` rather than `fetchable`, so it then
    applied the override -- incrementing override_applied -- and skipped the
    screen on no_gene_data one line later. The counter reported success for a
    screen contributing nothing.

    The only screens an override exists for are exactly the unresolvable ones:
    SCREEN_RATIONALE reporting both directions with an unsigned score, which
    needs curation from the paper.
    """
    src = (ROOT / "scripts" / "ingest_orcs.py").read_text()
    assert "overrides_early" in src, "fetchable is built before overrides load"
    i_fetchable = src.index("fetchable = [")
    i_overrides = src.index("overrides_early =")
    assert i_overrides < i_fetchable, \
        "overrides must be read BEFORE fetchable is built"
    # Span by lines, not by the first "]" -- that one sits inside
    # polarity_mode(s)[0] and ends the slice before the override clause.
    seg = "\n".join(src[i_fetchable:].splitlines()[:4])
    assert "overrides_early" in seg, \
        f"fetchable must keep an unresolvable screen that has an override:\n{seg}"


def test_override_turns_an_unresolvable_screen_into_a_screen_level_one():
    """The module half, so the two cannot drift apart."""
    screen = {"SCREEN_ID": 9001, "CONDITION_NAME": "SARS-CoV-2",
              "SCREEN_RATIONALE": "Increased/Decreased resistance to virus",
              "SCORE.1_TYPE": "MAGeCK pos score",   # unsigned -> unresolvable
              "METHODOLOGY": "Knockout", "CELL_LINE": "Huh-7.5",
              "AUTHOR": "Someone X (2021)", "SOURCE_ID": "12345"}
    assert polarity_mode(screen)[0] == "unresolvable"

    rows = [{"OFFICIAL_SYMBOL": "RAB7A", "IDENTIFIER_ID": "7879",
             "HIT": "YES", "SCORE.1": "3.0"}]
    em = Emit()
    stats = ingest_screens(em, [screen], {"9001": rows},
                           {"NCBIGene:7879": ["UniProtKB:P51149"]},
                           "infores:biogrid-orcs", {"9001": "dependency"})
    assert stats["override_applied"] == 1
    assert stats["edges"] == 1, "an applied override must produce edges"
    assert em.edges[0]["qualifiers"]["direction"] == "dependency"


def test_a_dengue_strain_screen_resolves_to_the_species_taxon():
    """Two ORCS screens name different DENV-2 strains -- Thailand/16681/84
    (31634, the strain the register pins for DENV-2's chains) and
    Jamaica/1409/1983. Both must land on species 11060, because that is the
    node the graph uses and the one M1's OrganismTaxon hop terminates on.
    """
    for strain in ("Thailand/16681/84", "Jamaica/1409/1983", "16681-PDK53"):
        text = f"Virus: Dengue virus 2 {strain}"
        assert virus_taxon({"CONDITION_NAME": text}) == "11060", text


def test_no_generic_dengue_pattern_shadows_the_serotypes():
    """The SARS-CoV/SARS-CoV-2 trap in dengue form. A bare r"dengue" pattern
    would match every serotype and file all four under whichever taxon it
    carried, silently. Serotype IS the distinction this project ingests, so
    each one is matched on its own digit and a generic pattern is banned.

    Asserted on the patterns themselves rather than on behaviour: a generic
    pattern added BELOW the serotypes would pass a behavioural test today and
    swallow DENV-1, DENV-3 and DENV-4 the moment they were added.
    """
    from kgav.orcs import VIRUS_PATTERNS

    for pattern, taxon in VIRUS_PATTERNS:
        src = pattern.pattern.lower()
        if "dengue" in src:
            assert any(d in src for d in "1234"), (
                f"generic dengue pattern {src!r} -> {taxon}: it would shadow "
                f"every serotype matched after it")


def test_the_flavivirus_patterns_do_not_disturb_the_coronaviruses():
    """Seven coronavirus screens were being ingested before these patterns
    were added, and pattern order decides every match."""
    for text, expected in (("Virus: SARS-CoV-2 (isolate USA-WA1/2020)", "2697049"),
                           ("Virus: SARS-CoV Urbani", "694009"),
                           ("Virus: HCoV-OC43", "31631"),
                           ("Virus: HCoV-229E", "11137"),
                           ("Virus: HCoV-NL63", "277944"),
                           ("Virus: HCoV-HKU1", "290028"),
                           ("Middle East respiratory syndrome coronavirus", "1335626")):
        assert virus_taxon({"CONDITION_NAME": text}) == expected, text


# ------------------------------------------- scope, not corruption (DENV-2)
#
# WHAT HAPPENED. The flavivirus patterns matched DENV-2 screens, ingest_screens
# emitted 90 HOST_FACTOR_FOR edges to NCBITaxon:11060, and the spine had never
# built that node because DENV-2 sits in `surveyed:` rather than `viruses:`.
# The schema check reported 90 DANGLING violations and the driver returned
# before writing, so 11,022 good edges for seven other viruses went with them.
# The output read as graph corruption; it was a scope mismatch.

DENV2 = "Virus: Dengue virus 2 (strain 16681)"
G2P = {"NCBIGene:7879": ["UniProtKB:P51149"]}
ROWS = {"1": [_gene("RAB7A", "7879", 3.0)]}


def test_a_screen_for_an_unbuilt_taxon_is_skipped_and_counted():
    em = Emit()
    s = ingest_screens(em, [_screen(condition=DENV2)], ROWS, G2P,
                       "infores:biogrid-orcs", known_taxa={"2697049"})
    assert em.edges == []
    assert s["taxon_not_built"] == 1
    assert s["edges"] == 0


def test_the_unbuilt_taxon_is_named_not_just_counted():
    """"1 screen skipped" does not tell you which virus to promote."""
    em = Emit()
    s = ingest_screens(em, [_screen(condition=DENV2)], ROWS, G2P,
                       "infores:biogrid-orcs", known_taxa={"2697049"})
    assert unbuilt_taxa(s) == [("11060", 1)]


def test_one_unbuilt_taxon_does_not_discard_every_other_virus():
    """THE REGRESSION. This is the whole point of the gate: 90 dangling
    DENV-2 edges cost 11,022 good ones for seven other viruses, because the
    schema check is all-or-nothing and ran over the whole batch."""
    em = Emit()
    s = ingest_screens(
        em,
        [_screen(sid="1", condition=DENV2),
         _screen(sid="2", condition="Virus: SARS-CoV-2 (isolate USA-WA1/2020)")],
        {"1": [_gene("RAB7A", "7879", 3.0)],
         "2": [_gene("RAB7A", "7879", 3.0)]},
        G2P, "infores:biogrid-orcs", known_taxa={"2697049"})
    assert s["edges"] == 1, "the built virus must still produce its edge"
    assert s["taxon_not_built"] == 1
    assert [e["object"] for e in em.edges] == ["NCBITaxon:2697049"]


def test_no_known_taxa_means_no_gate():
    """Callers that have not got a spine to check against are unchanged."""
    em = Emit()
    s = ingest_screens(em, [_screen(condition=DENV2)], ROWS, G2P,
                       "infores:biogrid-orcs")
    assert s["edges"] == 1 and s["taxon_not_built"] == 0


def test_an_unbuilt_taxon_emits_nothing_that_could_dangle():
    """The schema check must have nothing left to fail on."""
    em = Emit()
    ingest_screens(em, [_screen(condition=DENV2)], ROWS, G2P,
                   "infores:biogrid-orcs", known_taxa={"2697049"})
    assert not [e for e in em.edges if e["object"] == "NCBITaxon:11060"]


# ------------------------------------------------------- the per-taxon table
def test_the_per_taxon_rows_partition_the_edge_total():
    """One total for eight viruses cannot answer "which viruses have data",
    which is the question the flavivirus work exists to ask."""
    em = Emit()
    s = ingest_screens(
        em,
        [_screen(sid="1", condition="Virus: SARS-CoV-2 (isolate USA-WA1/2020)"),
         _screen(sid="2", condition="Virus: HCoV-OC43 (Human coronavirus OC43)")],
        {"1": [_gene("RAB7A", "7879", 3.0)],
         "2": [_gene("RAB7A", "7879", 3.0)]},
        G2P, "infores:biogrid-orcs")
    rows = per_taxon(s)
    assert sum(r[4] for r in rows) == s["edges"]
    assert {r[0] for r in rows} == {"2697049", "31631"}


def test_dependency_and_restriction_split_per_taxon():
    """BOTH poles must be present or this test cannot see the split at all.

    The first version of this used a dependency-only fixture, and collapsing
    the per-taxon buckets into one still passed it -- a test that could not
    fail for the reason it was written. RAB7A +11.31 is a dependency hit and
    TBK1 -3.33 a restriction hit in the same signed screen.
    """
    em = Emit()
    rows = {"1": [_gene("RAB7A", "7879", 11.31), _gene("TBK1", "29110", -3.33)]}
    s = ingest_screens(
        em, [_screen(sid="1",
                     rationale="Increased/decreased resistance to virus")],
        rows,
        {"NCBIGene:7879": ["UniProtKB:P51149"],
         "NCBIGene:29110": ["UniProtKB:Q9UHD2"]},
        "infores:biogrid-orcs")
    row = next(r for r in per_taxon(s) if r[0] == "2697049")
    assert row[2] == 1, "one dependency hit"
    assert row[3] == 1, "one restriction hit"
    assert row[2] + row[3] == row[4] == s["edges"] == 2


def test_a_taxon_with_a_screen_but_no_edge_is_distinguishable_from_no_screen():
    """A failed join and an absent source need different work, and a table
    that only lists taxa with edges silently merges them."""
    em = Emit()
    s = ingest_screens(
        em, [_screen(sid="1", condition="Virus: SARS-CoV-2 (isolate USA-WA1/2020)")],
        {"1": [_gene("NOPE", "999999", 3.0)]},      # gene not in the graph
        G2P, "infores:biogrid-orcs")
    assert s["edges"] == 0
    row = next(r for r in per_taxon(s) if r[0] == "2697049")
    assert row[1] == 1 and row[4] == 0, "screen counted, zero edges, still listed"


def test_the_screen_count_label_no_longer_claims_coronavirus():
    """virus_taxon matches every configured pattern, and the flavivirus ones
    landed 2026-10-11. A heading naming one family while counting two is the
    defect this repository keeps finding."""
    line = screen_heading(1952, 100, 8, 86)
    assert "coronavirus" not in line
    assert "100 matched a configured virus (8 taxa)" in line
    assert "1,952 human screens" in line


# --------------------------- what the index holds, before anything is fetched
_REG = {"2697049": "SARS-CoV-2", "64320": "Zika", "11060": "DENV-2",
        "11089": "Yellow fever", "11082": "West Nile", "11053": "DENV-1",
        "11072": "JEV", "11084": "TBEV"}


def _scr(text):
    return {"CONDITION_NAME": text, "PHENOTYPE": "", "NOTES": "",
            "SCREEN_NAME": ""}


def test_a_virus_with_screens_is_counted():
    lines = index_coverage([_scr("SARS-CoV-2 infection")] * 3, _REG)
    assert any("SARS-CoV-2" in ln and "3" in ln for ln in lines)


def test_a_declared_pattern_with_no_screen_says_so():
    """Not a download anyone skipped: the experiment is not in this
    release, and that is a different problem from a missing file."""
    text = "\n".join(index_coverage([_scr("Zika virus")], _REG))
    assert "Yellow fever" in text
    assert "pattern declared, NO screen in this release" in text


def test_a_virus_with_no_pattern_is_not_reported_as_a_zero():
    """The third state, and the one a coverage table silently turns into
    the second. DENV-1, JEV and TBEV cannot match by construction, so a
    zero beside them would be a claim about ORCS that nothing supports.
    """
    lines = index_coverage([_scr("Zika virus")], _REG)
    rows = [ln for ln in lines if ln.startswith("  ") and "<-" in ln]
    assert not any("DENV-1" in r or "JEV" in r or "TBEV" in r for r in rows)
    text = "\n".join(lines)
    assert "cannot match by construction" in text
    assert "DENV-1, JEV, TBEV" in text


def test_it_says_why_an_untested_pattern_is_not_added():
    text = "\n".join(index_coverage([], _REG))
    assert "an untested pattern is worse than an absent one" in text


def test_every_counted_taxon_has_a_declared_pattern():
    """The count comes from virus_taxon, so a taxon reported with screens
    must be one the patterns can produce."""
    from kgav.orcs import VIRUS_PATTERNS
    declared = {t for _p, t in VIRUS_PATTERNS}
    lines = index_coverage([_scr("west nile virus")] * 2, _REG)
    named = [ln for ln in lines if ln.startswith("  ") and "<-" not in ln
             and "register virus" not in ln]
    for ln in named:
        label = ln.strip().rsplit(None, 1)[0]
        taxon = next(t for t, n in _REG.items() if n == label)
        assert taxon in declared


def test_a_screen_is_counted_against_its_own_virus_only():
    """Attributing every screen to one taxon left all tests passing.

    Each fixture above happened to use screens for the virus it asserted
    about, so a count that ignored virus_taxon entirely gave the same
    answer. Here the only screens are West Nile, and SARS-CoV-2 must come
    back flagged as having none.
    """
    lines = index_coverage([_scr("west nile virus")] * 4, _REG)

    def row(label):
        return next(ln for ln in lines
                    if ln.startswith((f"  {label} ", f"  {label:<14}")))

    assert row("West Nile").split()[-1] == "4"
    assert "NO screen in this release" in row("SARS-CoV-2")
    assert "NO screen in this release" in row("Zika")


def test_a_matched_taxon_missing_from_the_register_is_not_dropped():
    """The report loops over the register, so a pattern pointing at a taxon
    the register no longer lists would count screens into a row that is
    never printed -- silent, and exactly the shape of every other defect
    this file has found."""
    partial = {k: v for k, v in _REG.items() if k != "11082"}
    text = "\n".join(index_coverage([_scr("west nile virus")] * 2, partial))
    assert "NOT in the register" in text and "11082" in text
