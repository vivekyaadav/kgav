"""Three defects in the dry run's own report.

THE SCRIPT EXISTS so a wrong assumption costs a read instead of a rebuild.
That makes its report the product, and a report that mislabels its evidence
is the failure this repository keeps finding -- here in the tool written to
catch it.

Driven end to end rather than through its internals: what is being checked
is the rendered text a person reads before deciding to rebuild.
"""
from __future__ import annotations

import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))


# ---------------------------------------------------------------- fixtures
#
# THESE BUILD THEIR OWN REGISTER. They used to promote DENV-2 and West Nile
# straight out of config/viruses.yaml, which worked only while those two sat
# in `surveyed:`. The 2026-10-23 promotion emptied that section and every one
# of them errored -- not because dryrun_promote changed, but because a test
# that knows which section a virus lives in asserts the register's shape as a
# side effect of asserting something else. test_virus_register.py carries the
# same warning about the same mistake.
#
# So the subject is a register these tests construct: the named viruses moved
# back into `surveyed:`, everything else untouched. What is under test is the
# dry run's behaviour on a promotable entry, which is true wherever the real
# register happens to have put it today.
def _as_surveyed(cfg: dict, *names: str) -> dict:
    import copy
    cfg = copy.deepcopy(cfg)
    moving = [v for v in cfg["viruses"] if v["name"] in names]
    assert len(moving) == len(names), (
        f"not in viruses:: {set(names) - {v['name'] for v in moving}}")
    cfg["viruses"] = [v for v in cfg["viruses"] if v["name"] not in names]
    cfg["surveyed"] = (cfg.get("surveyed") or []) + moving
    return cfg


def _register(tmp: Path, *names: str, edit=None) -> Path:
    """config/viruses.yaml with `names` moved back to `surveyed:`."""
    import yaml
    cfg = _as_surveyed(
        yaml.safe_load((ROOT / "config/viruses.yaml").read_text()), *names)
    if edit is not None:
        edit(cfg)
    tmp.mkdir(parents=True, exist_ok=True)
    path = tmp / "viruses.yaml"
    path.write_text(yaml.safe_dump(cfg))
    import dryrun_promote as dp
    dp._CFG_CACHE.clear()
    return path


@pytest.fixture(scope="module")
def report(tmp_path_factory) -> str:
    """The full dry-run text for promoting DENV-2 and West Nile."""
    import io
    from contextlib import redirect_stdout

    mod = import_module("dryrun_promote")
    tmp = tmp_path_factory.mktemp("reg")
    cfg = _register(tmp, "DENV-2", "West Nile")
    argv = ["dryrun_promote.py", "--promote", "DENV-2", "West Nile",
            "--config", str(cfg),
            "--proteomes", str(tmp_path_factory.mktemp("none"))]
    buf = io.StringIO()
    old = sys.argv
    sys.argv = argv
    try:
        with redirect_stdout(buf):
            mod.main()
    finally:
        sys.argv = old
    return buf.getvalue()


def _labelled(report: str) -> list[tuple[str, str]]:
    """(label line, evidence line) for every CANDIDATE in section 3."""
    lines = report.splitlines()
    return [(lines[i], lines[i + 1]) for i, l in enumerate(lines)
            if "CANDIDATE" in l and i + 1 < len(lines)]


# ------------------------------------------- 1. the label contradicted itself
def test_a_candidate_names_the_count_it_actually_matched():
    """It scanned for TWO numbers and named the first one for every hit.

    A line reading `== 12` printed under "compares against 8, the current
    built count" -- the label contradicting the evidence directly beneath
    it. The scan was right; only the sentence was wrong, and that is the
    kind of wrong that makes a reader discount a real break as noise.
    """
    import dryrun_promote as dp
    import yaml
    from kgav.virus_register import declared_taxa

    cfg = yaml.safe_load((ROOT / "config/viruses.yaml").read_text())
    n_viruses, n_taxa = len(cfg["viruses"]), len(declared_taxa(cfg))
    assert n_viruses != n_taxa, (
        "if these coincide this test cannot see a mislabel; pick another "
        "fixture rather than deleting the test")

    dp._CFG_CACHE.clear()
    dp._CFG_CACHE.update(cfg)
    for ref, n, line, matched in dp._hardcoded_counts():
        if n is None:
            assert f"== {matched}" in line.replace("==", "== ").replace(
                "==  ", "== "), f"{ref}: label says {matched}, line is {line}"


def test_the_two_counts_are_reported_under_different_names(report):
    labelled = _labelled(report)
    assert labelled, "no CANDIDATE rows; this fixture cannot test the labels"
    for label, evidence in labelled:
        num = label.split("compares against ")[1].split(",")[0].strip()
        assert f"== {num}" in evidence, (
            f"label claims {num} but the evidence line is: {evidence}")


def test_the_taxon_map_size_is_not_called_the_built_count(report):
    for label, evidence in _labelled(report):
        if "== 12" in evidence:
            assert "taxon-map size" in label, label
            assert "built count" not in label, label


# ------------------------------------------- 2. a command that cannot be run
def test_the_fetch_command_quotes_a_name_containing_a_space(report):
    """`--virus West Nile` parses as --virus=West plus a stray positional
    and exits 2. The line is printed to be pasted."""
    assert "--virus 'West Nile'" in report
    assert "--virus West Nile" not in report


def test_virus_is_one_string_per_flag_so_the_quoting_matters():
    """Guards the reason: if --virus ever took nargs='+', the quoting
    would be unnecessary and this pair of tests should be revisited rather
    than silently kept."""
    src = (ROOT / "scripts" / "fetch_raw.py").read_text()
    block = src[src.index('"--virus"'):]
    block = block[:block.index(")")]
    assert 'action="append"' in block
    assert "nargs" not in block


# ------------------------------------------- 3. a blocker that was fixed
def test_the_m8_note_reads_the_layer_listing_instead_of_asserting_it(report):
    """It said `folds` was unlisted in LAYER_PRECEDENCE. True when written,
    false since f1b5afe added it -- the build assembles the layer and ships
    FOLD_SIMILAR_TO edges. A dry run stating a stale blocker sends the
    reader to fix something already fixed."""
    from kgav.assemble import LAYER_PRECEDENCE
    assert "folds" in LAYER_PRECEDENCE, (
        "if folds is ever removed, the other branch is correct again and "
        "this test should assert that instead")
    assert "`folds` IS listed" in report
    assert "the only thing keeping M8 empty" not in report


def test_the_m8_note_would_still_report_a_real_blocker(monkeypatch, tmp_path):
    """The fix must not have turned the warning off, only made it read the
    listing. Removing folds must bring the original sentence back."""
    import io
    from contextlib import redirect_stdout

    import dryrun_promote as dp
    monkeypatch.setattr(dp, "LAYER_PRECEDENCE",
                        [x for x in dp.LAYER_PRECEDENCE if x != "folds"])
    cfg = _register(tmp_path / "_r", "DENV-2")
    buf = io.StringIO()
    old = sys.argv
    sys.argv = ["dryrun_promote.py", "--promote", "DENV-2",
                "--config", str(cfg)]
    try:
        with redirect_stdout(buf):
            dp.main()
    finally:
        sys.argv = old
    out = buf.getvalue()
    if "FAMILIES IN THE BUILD" in out:
        assert "the only thing keeping M8 empty" in out


# ----------------------------- 4. a precondition that read the wrong thing
def _section0(proteomes: Path, *promote: str) -> str:
    import io
    from contextlib import redirect_stdout

    mod = import_module("dryrun_promote")
    cfg = _register(proteomes.parent / "_reg", *promote)
    buf = io.StringIO()
    old = sys.argv
    sys.argv = ["dryrun_promote.py", "--promote", *promote,
                "--config", str(cfg), "--proteomes", str(proteomes)]
    try:
        with redirect_stdout(buf):
            mod.main()
    finally:
        sys.argv = old
    out = buf.getvalue()
    return out[out.index("0. PRECONDITIONS"):out.index("WHICH CHANNELS")]


def _fake(d: Path, stem: str, acc: str) -> None:
    import gzip
    import json
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{stem}.json.gz").write_bytes(gzip.compress(json.dumps(
        {"results": [{"primaryAccession": acc,
                      "entryType": "UniProtKB reviewed (Swiss-Prot)",
                      "proteinDescription": {"recommendedName": {
                          "fullName": {"value": "Genome polyprotein"}}},
                      "genes": [{"geneName": {"value": "POLY"}}],
                      "sequence": {"value": "MNNQ"},
                      "entryAudit": {"firstPublicDate": "2008-01-01"},
                      "features": []}]}).encode()))


def test_a_substituted_proteome_no_longer_blocks(tmp_path):
    """DENV-2's UP000095876 returned nothing from UniProt on 2026-10-23
    while West Nile's UP000008600 returned its entry the same minute. The
    register says the file's contents are replaced by EXTRA_11060, so
    requiring it blocked the promotion on a download nothing reads."""
    d = tmp_path / "p"
    _fake(d, "EXTRA_11060", "P29990")
    out = _section0(d, "DENV-2")
    assert "BLOCKED" not in out
    assert "n/a" in out
    assert "the curated accessions below replace it" in out


def test_a_missing_proteome_with_no_substitute_still_blocks(tmp_path):
    """West Nile declares no extra_accessions, so nothing replaces its
    proteome and the precondition must still fire."""
    out = _section0(tmp_path / "empty", "West Nile")
    assert "BLOCKED" in out
    assert "ABSENT" in out


def test_a_missing_extra_file_still_blocks(tmp_path):
    """DENV-2 with neither file: the substitute is declared but absent."""
    out = _section0(tmp_path / "empty", "DENV-2")
    assert "BLOCKED" in out


def test_the_clear_message_does_not_claim_a_proteome_is_cached(tmp_path):
    """It said "every promoted virus has its proteome cached", which is
    false for an entry that clears without one -- and would send someone
    looking for a download the build does not read."""
    d = tmp_path / "p"
    _fake(d, "EXTRA_11060", "P29990")
    out = _section0(d, "DENV-2")
    assert "has what ingest_spine reads" in out
    assert "has its proteome cached" not in out


# ------------------------- 5. one canned reason for five different blocks
def _fake_poly(d: Path, stem: str, acc: str, genes: bool) -> None:
    import gzip
    import json
    d.mkdir(parents=True, exist_ok=True)
    e = {"primaryAccession": acc,
         "entryType": "UniProtKB reviewed (Swiss-Prot)",
         "proteinDescription": {"recommendedName": {
             "fullName": {"value": "Genome polyprotein"}}},
         "sequence": {"value": "MNNQ"},
         "entryAudit": {"firstPublicDate": "2008-01-01"},
         "features": [{"type": "Chain", "featureId": f"PRO_{i}",
                       "description": f"NS{i}",
                       "location": {"start": {"value": 1 + i * 50},
                                    "end": {"value": 50 + i * 50}}}
                      for i in range(13)]}
    if genes:
        e["genes"] = [{"geneName": {"value": "POLY"}}]
    (d / f"{stem}.json.gz").write_bytes(
        gzip.compress(json.dumps({"results": [e]}).encode()))


def test_a_missing_gene_symbol_is_not_described_as_a_missing_proteome(
        tmp_path, monkeypatch):
    """West Nile was blocked on 2026-10-23 with its proteome present and
    validated -- 13 chains on P06935, printed two lines above -- and the
    report then said "ingest_spine would record proteome_missing". The
    cause was stated correctly and contradicted in the same section."""
    import dryrun_promote as dp

    def _strip(cfg):
        for e in cfg["surveyed"]:
            if e["name"] == "West Nile":
                e.pop("gene_symbol", None)
    path = _register(tmp_path / "_r", "West Nile", edit=_strip)

    d = tmp_path / "p"
    _fake_poly(d, "UP000008600", "P06935", genes=False)
    out = _section0_cfg(dp, path, d, "West Nile")
    assert "BLOCKED" in out
    assert "THE PROTEOME IS FINE" in out
    assert "would record proteome_missing" not in out


def test_a_block_no_download_fixes_offers_no_download(tmp_path):
    """It printed a fetch command for every virus being promoted whenever
    anything was blocked, so a missing gene symbol got a command that would
    re-download a file already correct."""
    import dryrun_promote as dp

    def _strip(cfg):
        for e in cfg["surveyed"]:
            if e["name"] == "West Nile":
                e.pop("gene_symbol", None)
    path = _register(tmp_path / "_r", "West Nile", edit=_strip)

    d = tmp_path / "p"
    _fake_poly(d, "UP000008600", "P06935", genes=False)
    out = _section0_cfg(dp, path, d, "West Nile")
    assert "FETCH THEM WITH" not in out


def test_a_missing_proteome_still_offers_the_download(tmp_path):
    """The split must not have removed the fetch hint where it belongs."""
    import dryrun_promote as dp
    dp._CFG_CACHE.clear()
    out = _section0(tmp_path / "empty", "West Nile")
    assert "FETCH THEM WITH" in out
    assert "would record proteome_missing" in out


def _section0_cfg(dp, config: Path, proteomes: Path, *promote: str) -> str:
    import io
    from contextlib import redirect_stdout
    buf = io.StringIO()
    old = sys.argv
    sys.argv = ["dryrun_promote.py", "--promote", *promote,
                "--config", str(config), "--proteomes", str(proteomes)]
    try:
        with redirect_stdout(buf):
            dp.main()
    finally:
        sys.argv = old
        dp._CFG_CACHE.clear()
    out = buf.getvalue()
    return out[out.index("0. PRECONDITIONS"):out.index("WHICH CHANNELS")]


# ------------------------------------------- 6. the declaration that fixes it
def test_west_nile_declares_the_polyprotein_gene_symbol():
    """UniProt carries no `genes` field on the flavivirus polyprotein entry,
    exactly as for Zika's Q32ZE1. POLY names the real ORF all thirteen
    mature peptides are cut from -- not a fiction to fill a path."""
    import yaml
    cfg = yaml.safe_load((ROOT / "config/viruses.yaml").read_text())
    wnv = next(e for e in (cfg["viruses"] + cfg["surveyed"])
               if e["name"] == "West Nile")
    assert wnv["gene_symbol"] == {"P06935": "POLY"}


def test_the_two_flavivirus_poly_genes_are_different_nodes():
    """ingest_spine keys Gene nodes KGAV:GENE_<taxon>_<symbol>, so West
    Nile's POLY and Zika's POLY cannot merge into one gene."""
    import yaml
    cfg = yaml.safe_load((ROOT / "config/viruses.yaml").read_text())
    entries = {e["name"]: e for e in cfg["viruses"] + cfg["surveyed"]}
    wnv, zika = entries["West Nile"], entries["Zika"]
    assert "POLY" in wnv["gene_symbol"].values()
    assert "POLY" in zika["gene_symbol"].values()
    assert wnv["taxon"] != zika["taxon"]
    assert (f"KGAV:GENE_{wnv['taxon']}_POLY"
            != f"KGAV:GENE_{zika['taxon']}_POLY")


# ---------------------- 7. the substituted entry was never actually checked
def test_a_substituted_entry_is_inspected_not_skipped(tmp_path):
    """Section 0's checks were gated on `if ok:` -- the PROTEOME file
    existing. DENV-2's is substituted, so it printed "n/a" and then nothing:
    no chain count, and no gene-symbol check, the one that caught West Nile.
    It read as clear because nothing had looked."""
    d = tmp_path / "p"
    _fake_poly(d, "EXTRA_11060", "P29990", genes=False)
    out = _section0(d, "DENV-2")
    assert "from EXTRA_11060.json.gz" in out
    assert "CHAIN features on P29990" in out
    assert "gene symbols:" in out


def test_a_substituted_entry_with_no_symbol_anywhere_still_blocks(tmp_path):
    """The check must not merely run -- it must reach the same verdict it
    would for a normal proteome."""
    import dryrun_promote as dp

    def _strip(cfg):
        for e in cfg["surveyed"]:
            if e["name"] == "DENV-2":
                e.pop("gene_symbol", None)
    path = _register(tmp_path / "_r", "DENV-2", edit=_strip)
    d = tmp_path / "p"
    _fake_poly(d, "EXTRA_11060", "P29990", genes=False)
    out = _section0_cfg(dp, path, d, "DENV-2")
    assert "NO GENE SYMBOL ANYWHERE" in out
    assert "BLOCKED" in out


# ------------------- 8. the register, after the 2026-10-23 retirement
def _flaviviruses():
    import yaml
    c = yaml.safe_load((ROOT / "config/viruses.yaml").read_text())
    return [v for v in c["viruses"] + c["surveyed"]
            if v.get("family") == "Flaviviridae"]


def test_every_flavivirus_on_a_strain_accession_declares_it():
    """THE DEFECT, SIX TIMES OVER. Seven of nine flavivirus proteome ids
    returned an empty result set from UniProt on 2026-10-23 -- every one of
    the TrEMBL-only species placeholders the register itself describes as
    having "a single protein and zero Chain features". The two that still
    resolve, Zika's UP000054557 and West Nile's UP000008600, are
    species-curated Swiss-Prot proteomes.

    DENV-2 was the only one of the seven that declared extra_accessions, so
    it was the only one whose chains could be reached without the retired
    file. The other six documented the same gap in prose and never applied
    the mechanism -- which is what DENV-2's own comment calls "the MERS
    mechanism, which this entry documented and did not apply".

    Coronaviruses are deliberately out of scope here: SARS-CoV and
    HCoV-HKU1 also pin strain chain sources and need no extra_accessions,
    because their species proteomes contain the strain entries.
    """
    missing = []
    for v in _flaviviruses():
        cs = v.get("chain_source") or {}
        if cs.get("taxon") in (None, v["taxon"]):
            continue                      # species-curated, nothing to reach
        acc = cs["accession"]
        if acc not in (v.get("extra_accessions") or []):
            missing.append(f"{v['name']}: chains on {acc} "
                           f"(strain {cs['taxon']}), extra_accessions="
                           f"{v.get('extra_accessions')}")
    assert not missing, missing


def test_a_declared_extra_accession_is_the_one_carrying_the_chains():
    """extra_accessions REPLACES the proteome, so an extra list that does
    not include chain_source.accession substitutes entries without the
    chains and builds a virus with none."""
    for v in _flaviviruses():
        extras = v.get("extra_accessions") or []
        if not extras:
            continue
        acc = (v.get("chain_source") or {}).get("accession")
        assert acc in extras, f"{v['name']}: {acc} not in {extras}"


def test_a_declared_gene_symbol_names_an_accession_the_entry_uses():
    """A symbol keyed on an accession the build never loads is dead
    config that reads as coverage."""
    for v in _flaviviruses():
        declared = v.get("gene_symbol") or {}
        if not declared:
            continue
        known = set(v.get("extra_accessions") or [])
        known.add((v.get("chain_source") or {}).get("accession"))
        if v["name"] == "Zika":
            known.add("P0DXN9")       # the uORF isoform, same gene
        unknown = set(declared) - known
        assert not unknown, f"{v['name']}: {unknown} not loaded by this entry"
