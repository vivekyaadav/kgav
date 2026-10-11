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
    """"Not needed" must be distinguishable from "forgotten".

    Searched across BOTH sections. This looked in `surveyed:` only and broke
    the moment Zika was promoted -- a test that knows which section a virus
    lives in asserts the register's shape as a side effect of asserting a
    fact about the virus, and fails on a move that changes neither.
    """
    for name in ("Zika", "West Nile"):
        v = next(x for s in ("viruses", "surveyed")
                 for x in (doc.get(s) or []) if x["name"] == name)
        assert str(v["chain_source"]["taxon"]) == str(v["taxon"])
        assert not v.get("isolate_taxa")


# Every virus the build is MEANT to contain, by name. An accidental promotion
# fails against this list and names the virus; a deliberate one adds a line
# that says which virus, which `len(...) == 7` could not.
#
# That assertion was the previous guard and it had to be re-derived by hand
# every time: dryrun_promote.py reports "pins 7, would be 9" precisely because
# a count cannot say what changed.
BUILT = (
    "SARS-CoV-2", "SARS-CoV", "MERS-CoV",
    "HCoV-229E", "HCoV-NL63", "HCoV-OC43", "HCoV-HKU1",
    # Zika, promoted 2026-10-12. DENV-2 was promoted alongside it on
    # 2026-10-11 and both were moved back the same day, because neither
    # proteome was cached and both built an OrganismTaxon with zero proteins.
    # This time dryrun_promote reported Section 0 clear first: UP000054557
    # cached, Q32ZE1 carrying the 13 chains the register claims, and M1 and M7
    # reaching the taxon. DENV-2 stays surveyed -- at parity on rows that name
    # the virus they are filed under (392 vs 355) for several times the
    # handling.
    "Zika",
)

# THE FOUR "CANDIDATE" HITS dryrun_promote REPORTS FOR THIS PROMOTION ARE
# FALSE POSITIVES, checked on 2026-10-12 so the next promotion need not
# re-derive it. Its scan flags any assertion comparing against the current
# built count, which is the right heuristic -- it once reported one failure
# where three followed -- but it cannot read intent:
#
#   test_calibration.py:251,261,270  `len(cal.channels) == 7` counts METAPATH
#       channels from the local MEASURED fixture (M1-M7, with COMBINED and
#       degree excluded by load()). Keyed to one virus, independent of how
#       many are built.
#   test_ncats.py:334  `stats["in_assay"] == 7` counts rows in the fixture's
#       cpe.tsv assay.
#
# Only BUILT above had to change. Editing those four would have made them
# assert 8 of things there are 7 of.


def test_surveyed_viruses_are_not_built(doc):
    """ingest_spine reads `viruses:` only. Nothing in `surveyed:` may leak."""
    built = {v["taxon"] for v in doc["viruses"]}
    assert not built & {v["taxon"] for v in doc["surveyed"]}


def test_the_built_set_is_exactly_what_is_declared_built(doc):
    """Named rather than counted, so the failure says which virus."""
    names = [v["name"] for v in doc["viruses"]]
    assert names == list(BUILT), (
        f"built set changed\n  expected: {list(BUILT)}\n  found:    {names}")


def test_every_built_virus_has_what_a_build_needs(doc):
    """A count never checked this. A promoted entry missing proteome_id or an
    incomplete chain_source builds a virus with no targets and no error."""
    for v in doc["viruses"]:
        ref = v.get("name")
        assert v.get("proteome_id"), f"{ref}: no proteome_id"
        cs = v.get("chain_source") or {}
        for field in ("accession", "taxon", "chains", "verified"):
            assert cs.get(field) not in (None, ""), f"{ref}: chain_source.{field}"
        assert cs["chains"] > 0, f"{ref}: {cs['chains']} chains"
        assert "mondo" in v, f"{ref}: mondo must be present, explicitly null"


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


def test_no_script_builds_the_taxon_map_by_hand():
    """declared_taxa exists for this and its own docstring says so -- "the same
    shape the ingests build ad hoc from isolate_taxa, derived here so
    chain_source cannot drift out of the map" -- and for four scripts nothing
    called it. ingest_chembl, ingest_selectivity, ingest_vhppi and
    diagnose_assay_type each rebuilt the mapping inline, so the canonical
    function was reached only by its own test: five implementations of one
    mapping.

    That is the shape of three bugs already paid for in this project -- AUC
    0.826 in guards.py and tools.py, the coverage floor, and lift's
    trustworthy rule, each fixed in one copy and left broken in the other.

    Detected on the AST rather than by grepping for text, so a reformatted
    copy still trips it.
    """
    import ast
    from pathlib import Path

    scripts = Path(__file__).resolve().parents[1] / "scripts"
    offenders = []
    for path in sorted(scripts.glob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            # for iso in v.get("isolate_taxa") or []:  -> the hand-built copy
            if not isinstance(node, ast.For):
                continue
            src = ast.dump(node.iter)
            if "isolate_taxa" in src:
                offenders.append(f"{path.name}:{node.lineno}")
    assert not offenders, (
        "taxon map rebuilt inline instead of via "
        f"virus_register.declared_taxa: {offenders}")


def test_declared_taxa_maps_every_isolate_to_its_species(doc):
    """What the consolidation has to preserve: VirHostNet files MERS under
    the HCoV-EMC/2012 isolate and HKU1 under a strain taxon, and querying
    either species returns zero rows. The map is the only thing standing
    between those interactions and vanishing as an unknown organism.
    """
    m = declared_taxa(doc, "viruses")
    for entry in doc["viruses"]:
        species = str(entry["taxon"])
        assert m[species] == species
        for iso in entry.get("isolate_taxa") or []:
            assert m[str(iso)] == species, f"{iso} must map to {species}"
        chain = str((entry.get("chain_source") or {}).get("taxon") or "")
        if chain:
            assert m[chain] == species, (
                f"chain_source taxon {chain} missing from the map: this is "
                f"the line the hand-built copies omitted")


def test_the_sars_cov_ppi_isolates_stay_declared(doc):
    """VirHostNet files 73 host-viral interactions under two SARS-CoV isolate
    taxa, identified from UniProt taxonomy with 694009 in both lineages:

        227984  SARS coronavirus Tor2      66 rows
        698419  SARS coronavirus wtic-MB    7 rows

    THIS TEST ASSERTS THE MAPPING, NOT AN EFFECT. The first version of this
    docstring claimed removing these entries loses the graph 73 edges. It
    does not: ingest_vhppi produces a byte-identical layer with and without
    them, because `unmapped_viral_taxon` is 0 either way and the rows are
    dropped earlier at `unresolved_id`. The declaration is correct and inert
    until extra_accessions lets the spine resolve the isolates' identifiers.
    """
    m = declared_taxa(doc, "viruses")
    for iso in ("227984", "698419"):
        assert m.get(iso) == "694009", (
            f"{iso} must map to SARS-CoV; undeclared it is discarded as "
            f"unmapped_viral_taxon")


def test_a_ppi_isolate_is_not_mistaken_for_a_chain_source(doc):
    """SARS-CoV's chain_source is the species itself. These isolates exist
    only so the PPI layer can resolve them, and nothing should start reading
    chains from a recombinant clone because a taxon appeared in the list.
    """
    sars = next(v for v in doc["viruses"] if str(v["taxon"]) == "694009")
    assert str(sars["chain_source"]["taxon"]) == "694009"
    assert "227984" in [str(x) for x in sars["isolate_taxa"]]
    assert sars["chain_source"]["accession"] == "P0C6X7"
