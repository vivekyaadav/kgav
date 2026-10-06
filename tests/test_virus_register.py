"""Gate for the virus register: a descendant taxon must be DECLARED.

The trap: Swiss-Prot curates many viruses under a strain taxon while the
species taxon carries a TrEMBL proteome with no Chain features. A
`taxonomy_id:` lookup walks down the lineage and silently succeeds, recording
nothing about which strain it chose. That is how MERS lost Mpro and RdRp
(species 1335626 is TrEMBL-only; curation is under 1263720).
"""
import copy
from pathlib import Path

import pytest

from kgav.virus_register import declared_taxa, load, validate

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def doc():
    return load()


def _entry(**over):
    base = {"name": "Test", "taxon": "111", "mondo": None,
            "chain_source": {"accession": "P00000", "taxon": "111",
                             "chains": 13, "verified": "2026-09-17"}}
    base.update(over)
    return {"viruses": [base]}


# ------------------------------------------------------------ the real file
def test_the_committed_register_is_valid(doc):
    assert validate(doc) == []


def test_every_entry_declares_a_chain_source(doc):
    for section in ("viruses", "surveyed"):
        for v in doc[section]:
            assert v.get("chain_source"), f"{v['name']} has no chain_source"


def test_every_strain_level_virus_declares_its_isolate(doc):
    """The seven flaviviruses and two coronaviruses curated under a strain."""
    strain_level = []
    for section in ("viruses", "surveyed"):
        for v in doc[section]:
            ct = str(v["chain_source"]["taxon"])
            if ct != str(v["taxon"]):
                strain_level.append(v["name"])
                assert ct in [str(x) for x in (v.get("isolate_taxa") or [])], \
                    f"{v['name']}: chains from {ct}, not in isolate_taxa"
    # MERS + HKU1 + DENV 1-4 + JEV + YFV + TBEV
    assert len(strain_level) == 9, strain_level


def test_zika_and_west_nile_are_species_level_and_declare_no_isolate(doc):
    """"Not needed" must be distinguishable from "forgotten"."""
    for name in ("Zika", "West Nile"):
        v = next(x for x in doc["surveyed"] if x["name"] == name)
        assert str(v["chain_source"]["taxon"]) == str(v["taxon"])
        assert not v.get("isolate_taxa")


def test_surveyed_viruses_are_not_built(doc):
    """ingest_spine reads `viruses:` only. Nothing here may leak into v1."""
    built = {v["taxon"] for v in doc["viruses"]}
    assert not built & {v["taxon"] for v in doc["surveyed"]}
    assert len(doc["viruses"]) == 7


# ------------------------------------------------------------- the rule
def test_an_undeclared_descendant_taxon_is_rejected():
    bad = _entry(chain_source={"accession": "P29990", "taxon": "31634",
                               "chains": 13, "verified": "2026-09-17"})
    codes = {v.code for v in validate(bad)}
    assert "UNDECLARED_CHAIN_TAXON" in codes


def test_the_same_descendant_taxon_declared_is_accepted():
    ok = _entry(isolate_taxa=["31634"],
                chain_source={"accession": "P29990", "taxon": "31634",
                              "chains": 13, "verified": "2026-09-17"})
    assert validate(ok) == []


def test_a_missing_chain_source_is_rejected():
    bad = {"viruses": [{"name": "T", "taxon": "111", "mondo": None}]}
    assert "NO_CHAIN_SOURCE" in {v.code for v in validate(bad)}


def test_an_incomplete_chain_source_is_rejected():
    bad = _entry(chain_source={"accession": "P1", "taxon": "111"})
    assert "INCOMPLETE_CHAIN_SOURCE" in {v.code for v in validate(bad)}


def test_zero_chains_is_rejected():
    bad = _entry(chain_source={"accession": "P1", "taxon": "111", "chains": 0,
                               "verified": "2026-09-17"})
    assert "NO_CHAINS" in {v.code for v in validate(bad)}


def test_an_isolate_equal_to_the_species_is_rejected():
    bad = _entry(isolate_taxa=["111"])
    assert "ISOLATE_IS_SPECIES" in {v.code for v in validate(bad)}


def test_a_duplicate_taxon_across_sections_is_rejected(doc):
    bad = copy.deepcopy(doc)
    bad["surveyed"].append(copy.deepcopy(bad["viruses"][0]))
    assert "DUPLICATE_TAXON" in {v.code for v in validate(bad)}


def test_absent_mondo_is_rejected_but_explicit_null_is_fine():
    missing = {"viruses": [{"name": "T", "taxon": "111",
                            "chain_source": {"accession": "P1", "taxon": "111",
                                             "chains": 1, "verified": "x"}}]}
    assert "MONDO_UNSTATED" in {v.code for v in validate(missing)}
    assert "MONDO_UNSTATED" not in {v.code for v in validate(_entry())}


def test_declared_taxa_includes_chain_source_taxa(doc):
    """The ingests' taxon map must cover the strain a chain really came from,
    or the join fails exactly where the species query already failed."""
    m = declared_taxa(doc, "viruses")
    assert m["1263720"] == "1335626"   # MERS isolate -> species
    assert m["443239"] == "290028"     # HKU1 isolate -> species


# --------------------------------------------------------------------------
# Two sources of truth for "where do this virus's chains live". They drift.
# --------------------------------------------------------------------------
def _fold_specs():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "ifolds", ROOT / "scripts" / "ingest_folds.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.SPECS


@pytest.mark.xfail(strict=True, reason=(
    "KNOWN GAP, not yet curated. HCV (11103) and the five picornaviruses are "
    "in ingest_folds.SPECS and absent from the register. Closing it needs a "
    "chain_source verified against UniProt on a real date, which is a curation "
    "decision. strict=True so this turns into a FAILURE the moment the entries "
    "land, prompting removal of this marker rather than letting the check go "
    "quietly green."))
def test_every_fold_spec_taxon_is_declared_in_the_register():
    """ingest_folds.py states that it is, and for HCV that is false.

    Its header reads: "Every strain taxon here is now also registered in
    config/viruses.yaml under `surveyed:`, with chain_source naming the
    accession and the taxon it really lives under, and
    kgav.virus_register.validate enforces that the pairing is declared rather
    than discovered."

    HCV is not in the register. Its SPECS comment records the reassignment --
    P26664 filed under 11104 (genotype 1a isolate 1), declared as 11103
    (species) -- which is precisely the MERS trap the register exists to make
    checkable, and validate() cannot check an entry that does not exist. HCV is
    also the virus the module's own docstring says phase 4 depends on.

    The picornaviruses are undeclared too. They may well be species-curated,
    but nothing states it, so nobody can tell "verified as not needed" from
    "never looked" -- the distinction the Zika and West Nile entries are
    explicitly written to preserve.
    """
    doc = load(ROOT / "config" / "viruses.yaml")
    declared = {str(v["taxon"]) for v in (doc.get("viruses") or [])} | \
               {str(v["taxon"]) for v in (doc.get("surveyed") or [])}
    missing = sorted({
        (s.virus_taxon.split(":")[1], s.label)
        for s in _fold_specs() if s.virus_taxon.split(":")[1] not in declared})
    assert not missing, (
        "fold SPECS name taxa the register does not declare, so the "
        "species/strain pairing for them is unchecked:\n  "
        + "\n  ".join(f"{t}  {lab}" for t, lab in missing))


def test_every_fold_spec_accession_matches_the_registers_chain_source():
    """The accession and the taxon must come from the same declaration.

    A spec pointing at an accession the register does not name for that virus
    means the two were edited independently, which is how the pairing silently
    stops being the one that was verified.
    """
    doc = load(ROOT / "config" / "viruses.yaml")
    by_taxon = {str(v["taxon"]): (v.get("chain_source") or {}).get("accession")
                for v in (doc.get("viruses") or []) + (doc.get("surveyed") or [])}
    mismatched = []
    for s in _fold_specs():
        t = s.virus_taxon.split(":")[1]
        if t in by_taxon and by_taxon[t] and s.accession != by_taxon[t]:
            mismatched.append(f"{s.label}: spec {s.accession}, "
                              f"register {by_taxon[t]}")
    assert not mismatched, "\n  ".join(["spec/register accession mismatch:"]
                                       + mismatched)
