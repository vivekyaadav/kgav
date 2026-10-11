"""Virus-name detection in assay descriptions, on real strings.

Every description below was printed by diagnose_flavivirus_targets.py against
the real ChEMBL database on 2026-10-11, or is a shortened form of one. None
was composed to make a pattern pass -- the reason being this repo's own
history, where a MITAB regex was verified against a fixture whose author had
written the same wrong assumption into it.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _mod():
    spec = importlib.util.spec_from_file_location(
        "diag_serotype", ROOT / "scripts" / "diagnose_serotype_attribution.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


d = _mod()


# Real DENV-2-filed rows naming each serotype.
@pytest.mark.parametrize("desc,expected", [
    ("Inhibition of dengue virus-4 NS2B-NS3 protease", ["DENV-4"]),
    ("Inhibition of dengue virus-3 NS2B-NS3 protease", ["DENV-3"]),
    ("Inhibition of dengue virus-1 NS2B-NS3 protease", ["DENV-1"]),
    ("Inhibition of dengue virus-2 NS2B-NS3 protease", ["DENV-2"]),
    ("Non-competitive inhibition of Dengue virus 4 NS3 helicase", ["DENV-4"]),
    ("Non-competitive inhibition of Dengue virus 3 NS3 helicase", ["DENV-3"]),
    ("Inhibition of Dengue virus serotype 2 NS2B-NS3 protease by homogeneous "
     "fluorimetric assay", ["DENV-2"]),
    ("Inhibition of Dengue virus type 2 NS5 RNA methyltransferase SAM site",
     ["DENV-2"]),
    ("Inhibition of Dengue virus 2 NS2B/NS3 protease using "
     "Bz-Nle-Lys-Arg-Arg-AMC as substrate", ["DENV-2"]),
    ("Inhibitory activity against CF40.NS3pro from Dengue virus type 2",
     ["DENV-2"]),
    ("Inhibition of RNA unwinding activity of thioredoxin/His6 tagged Dengue "
     "virus 2 TSV01 NS3 helicase", ["DENV-2"]),
])
def test_dengue_serotypes_are_read_from_real_descriptions(desc, expected):
    assert d.viruses_named(desc) == expected


def test_the_west_nile_serotype_trap():
    """THE PATTERN THIS TEST EXISTS FOR.

    A real West Nile row reads "Inhibition of West Nile virus serotype 2
    NS2B-NS3 protease by homogeneous fluorimetric assay". A pattern matching a
    bare serotype digit -- serotype\\s*([1-4]) -- reads that as dengue
    serotype 2 and reports a West Nile row as a dengue mismatch: inventing
    precisely the finding this script is meant to detect, which would then be
    read back as evidence.

    So every dengue pattern requires "dengue" or "denv" adjacent to the digit.
    """
    desc = ("Inhibition of West Nile virus serotype 2 NS2B-NS3 protease by "
            "homogeneous fluorimetric assay")
    assert d.viruses_named(desc) == ["West Nile"]
    assert "DENV-2" not in d.viruses_named(desc)


@pytest.mark.parametrize("desc,expected", [
    ("Binding affinity to Zika virus NS2B-NS3 protease", ["Zika"]),
    ("Binding affinity to ZIKV RdRP by SPR assay", ["Zika"]),
    ("Inhibition of ZIKV 2'-O-MTase assessed as reduction in [3H]AdoMet "
     "transfer on to capped Gpp", ["Zika"]),
    ("Inhibition of WNV recombinant NS2B-NS3 protease", ["West Nile"]),
    ("Inhibition of West Nile virus NS2B-NS3 protease", ["West Nile"]),
    ("Inhibition of WNV 2'-O-MTase assessed as reduction in [3H]AdoMet "
     "transfer on to capped Gppp", ["West Nile"]),
])
def test_the_other_flaviviruses_are_read_from_real_descriptions(desc, expected):
    assert d.viruses_named(desc) == expected


def test_an_unspecified_dengue_row_names_no_serotype():
    """"Inhibition of dengue virus NS2B-NS3 protease" is a real row and names
    no serotype. Reading it as DENV-2 because it is filed there would hide the
    curation gap; reading it as a mismatch would invent one."""
    assert d.viruses_named("Inhibition of dengue virus NS2B-NS3 protease") == []
    assert d.GENERIC.search("Inhibition of dengue virus NS2B-NS3 protease")


def test_the_generic_screen_title_is_generic_not_silent():
    """DENV-2's largest single assay, 1,013 rows: "PUBCHEM_BIOASSAY: Primary
    and Confirmatory Screening for Flavivirus Genomic Capping Enzyme". It
    names a family, not a virus -- counted apart from a description that names
    nothing at all, because a generic label is a curation choice and a silence
    is not."""
    desc = ("PUBCHEM_BIOASSAY: Primary and Confirmatory Screening for "
            "Flavivirus Genomic Capping Enzyme")
    assert d.viruses_named(desc) == []
    assert d.GENERIC.search(desc)


def test_a_panel_description_names_more_than_one():
    """Several viruses in one description is a cross-reactivity panel, not a
    misfiling, and must not be collapsed into either bucket."""
    named = d.viruses_named(
        "Inhibition of dengue virus-2 and West Nile virus NS2B-NS3 protease")
    assert set(named) == {"DENV-2", "West Nile"}


def test_a_two_digit_number_is_not_a_serotype():
    """The (?!\\d) guard. "DENV-24" is not serotype 2, and a residue number or
    concentration must not be read as one."""
    assert d.viruses_named("Inhibition of DENV-24 construct") == []
    assert d.viruses_named("dengue virus 21 isolate") == []


def test_every_pattern_maps_to_a_taxon_in_the_register():
    """A label whose taxon is not declared cannot be compared against the row
    it is filed under, so the mapping would silently do nothing."""
    import yaml

    from kgav.virus_register import declared_taxa

    cfg = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    declared = set(declared_taxa(cfg)) | set(declared_taxa(cfg, "surveyed"))
    for _pat, taxon, label in d.NAME_PATTERNS:
        assert taxon in declared, f"{label} ({taxon}) is not in the register"


def test_the_pattern_labels_match_the_register_names_exactly():
    """The script compares viruses_named() output against v["name"], so a
    label that reads "DENV2" where the register says "DENV-2" would make every
    correctly-filed row look like a mismatch -- a false finding, reported
    confidently, which is the failure mode this whole script is about."""
    import yaml

    cfg = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    names = {v["name"] for sec in ("viruses", "surveyed")
             for v in (cfg.get(sec) or [])}
    for _pat, _taxon, label in d.NAME_PATTERNS:
        assert label in names, (
            f"{label!r} is not a register name; the comparison against "
            f"v['name'] would never match")


# --------------------------------------------------------------------------
# every bucket is printed, and the distinct count is of the bucket
# --------------------------------------------------------------------------
GENERIC_DESC = ("PUBCHEM_BIOASSAY: Primary and Confirmatory Screening for "
                "Flavivirus Genomic Capping Enzyme")


def _rows(n_generic=5, n_named=3, n_other=1, n_none=1):
    def row(desc, pmid):
        return {"description": desc, "target_type": "SINGLE PROTEIN",
                "target_name": "Genome polyprotein", "target": "CHEMBL5980",
                "assay_type": "B", "pubmed_id": pmid}
    return (
        [row(GENERIC_DESC, None)] * n_generic
        + [row("Inhibition of Dengue virus type 2 NS5 RNA methyltransferase "
               "SAM site", 20108931)] * n_named
        + [row("Inhibition of dengue virus-4 NS2B-NS3 protease", 12345)] * n_other
        + [row("Inhibition of NS3 protease construct", 999)] * n_none)


def _run(monkeypatch, capsys, rows):
    monkeypatch.setattr(d, "query_activities", lambda db, taxa: rows)
    d.main(["--db", str(ROOT / "config" / "viruses.yaml"), "--samples", "3"])
    return capsys.readouterr().out.split("=" * 78)[2]      # DENV-2's section


def test_the_largest_bucket_is_printed_verbatim(monkeypatch, capsys):
    """THE OMISSION THIS FIXES. The first version counted "generic flavivirus
    name" and "names no virus" and printed neither -- for DENV-2 that was
    1,045 + 52 of 1,483 rows tallied and never shown, the largest share of the
    substrate left invisible.

    This project already had that exact bug once: a 1,025-row bucket reported
    without being printed, then guessed at, and the guess was wrong. A count
    with no verbatim text behind it is not something a reader can check.
    """
    out = _run(monkeypatch, capsys, _rows())
    assert "[generic flavivirus name]" in out
    assert GENERIC_DESC[:60] in out, "the biggest bucket must show its text"
    assert "[names no virus]" in out
    assert "Inhibition of NS3 protease construct" in out


def test_the_distinct_count_covers_the_bucket_not_the_sample(monkeypatch,
                                                             capsys):
    """Counting distinct values inside a display loop capped at --samples
    reports the sample size as the total. Done once in this repo already, with
    a dedupe on the first 100 characters on top. Here: 8 rows, 4 distinct
    descriptions, --samples 3.
    """
    rows = _rows(n_generic=5, n_named=3, n_other=0, n_none=0)
    out = _run(monkeypatch, capsys, rows)
    assert "5 row(s), 1 distinct" in out, out
    # And with more distinct values than the sample cap:
    varied = [{"description": f"Inhibition of flavivirus assay variant {i}",
               "target_type": "SINGLE PROTEIN", "target_name": "Genome polyprotein",
               "target": "CHEMBL5980", "assay_type": "B", "pubmed_id": None}
              for i in range(7)]
    out2 = _run(monkeypatch, capsys, varied)
    assert "7 row(s), 7 distinct" in out2, out2


def test_the_crosstab_separates_naming_from_publication(monkeypatch, capsys):
    """The cell that decides build order. The two controls are independent:
    naming the virus makes a row attributable, carrying a pubmed_id makes it
    withholdable by --publication-disjoint. A row missing either is not
    defensible evidence, and DENV-2's substrate is mostly rows missing both.
    """
    out = _run(monkeypatch, capsys, _rows(n_generic=5, n_named=3, n_other=1,
                                          n_none=1))
    assert "ATTRIBUTABLE *AND* WITHHOLDABLE" in out
    # 3 named rows all carry a PMID; the 5 generic rows carry none.
    assert "names this virus                   3         0" in out, out
    assert "3 of 10 rows (30.0%)" in out, out


def test_a_virus_whose_rows_all_name_it_reports_cleanly(monkeypatch, capsys):
    """Zika and West Nile both measured 100% on the real database, so the
    clean case must read as clean rather than printing empty buckets."""
    rows = [{"description": "Binding affinity to Zika virus NS2B-NS3 protease",
             "target_type": "SINGLE PROTEIN", "target_name": "Genome polyprotein",
             "target": "CHEMBL123", "assay_type": "B", "pubmed_id": 33596380}] * 4
    monkeypatch.setattr(d, "query_activities", lambda db, taxa: rows)
    d.main(["--db", str(ROOT / "config" / "viruses.yaml"), "--samples", "3"])
    zika = capsys.readouterr().out.split("=" * 78)[4]
    assert "names this virus" in zika
    assert "WRONG-SPECIES RISK" not in zika
    assert "4 of 4 rows (100.0%)" in zika, zika
