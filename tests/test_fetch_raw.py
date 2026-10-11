"""scripts/fetch_raw.py -- no network, injected fetcher.

The fixtures here are REAL responses, trimmed: accessions, organism, chain
counts and chain descriptions were read from rest.uniprot.org on 2026-10-11
and are reproduced rather than imagined. That matters because this project has
already shipped a bug verified against a fixture whose author also wrote the
expectation -- the VirHostNet MITAB regex, which assumed
`taxid:9606(Homo sapiens)` when the real files carry bare `taxid:NNNNN`, and
whose fixture encoded the same assumption and so could not falsify it.

A fixture nobody can check is not evidence. Each one below names its query.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


def _mod():
    spec = importlib.util.spec_from_file_location(
        "fetch_raw", ROOT / "scripts" / "fetch_raw.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m          # @dataclass/annotations resolve via this
    spec.loader.exec_module(m)
    return m


fr = _mod()


# --------------------------------------------------------------------------
# fixtures: real, trimmed, each naming its query
# --------------------------------------------------------------------------
# rest.uniprot.org/uniprotkb/stream?query=(proteome:UP000054557) 2026-10-11
# -> {"results": [...]} with 2 entries, Q32ZE1 and P0DXN9.
# The 13 Chain descriptions on Q32ZE1 are verbatim.
ZIKA_CHAINS = [
    "Genome polyprotein", "Capsid protein C", "Protein prM", "Peptide pr",
    "Small envelope protein M", "Envelope protein E",
    "Non-structural protein 1", "Non-structural protein 2A",
    "Serine protease subunit NS2B", "Serine protease NS3",
    "Non-structural protein 4A", "Non-structural protein 4B",
    "RNA-directed RNA polymerase/Methyltransferase NS5",
]


def _entry(acc: str, chain_descriptions: list[str]) -> dict:
    return {
        "primaryAccession": acc,
        "features": [{"type": "Chain", "description": d}
                     for d in chain_descriptions],
    }


def _body(entries: list[dict]) -> bytes:
    return gzip.compress(json.dumps({"results": entries}).encode())


ZIKA_PROTEOME = _body([_entry("Q32ZE1", ZIKA_CHAINS), _entry("P0DXN9", [])])

# P29990 (DENV-2 strain Thailand/16681/1984, taxon 31634) carries the SAME 13
# chain descriptions and has NO proteome cross-reference at all, so a proteome
# download can never produce it. UP000095876 therefore returns other entries.
DENV2_PROTEOME = _body([_entry("A0A2H4YJ12", []), _entry("A0A2H4YJ99", [])])
DENV2_EXTRA = _body([_entry("P29990", ZIKA_CHAINS)])


@pytest.fixture
def register():
    return yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())


def _virus(register, name):
    return fr.find_virus(register, name)[0]


# --------------------------------------------------------------------------
# URLs
# --------------------------------------------------------------------------
def test_proteome_url_asks_the_stream_endpoint_for_compressed_json():
    """`search` paginates and truncates; `stream` returns the whole set.

    A truncated proteome is the quietest possible failure: some proteins, no
    error, and a release that validates.
    """
    url = fr.proteome_url("UP000054557")
    assert url.startswith("https://rest.uniprot.org/uniprotkb/stream?")
    assert "query=%28proteome%3AUP000054557%29" in url
    assert "format=json" in url and "compressed=true" in url


def test_accession_url_ors_every_accession(register):
    url = fr.accession_url(["P29990", "Q32ZE1"])
    assert "accession%3AP29990+OR+accession%3AQ32ZE1" in url


def test_every_built_virus_yields_a_fetchable_target(register):
    """The register must be enough to fetch from. If a virus is declared built
    and this cannot name its files, `make` has no way to get them."""
    for v in register["viruses"]:
        targets = fr.targets_for(v)
        assert targets, v["name"]
        assert targets[0][0].name == f"{v['proteome_id']}.json.gz"
        if v.get("extra_accessions"):
            assert targets[1][0].name == f"EXTRA_{v['taxon']}.json.gz"


# --------------------------------------------------------------------------
# inspect: a 200 is not a proteome
# --------------------------------------------------------------------------
def test_inspect_reads_the_real_zika_shape():
    r = fr.inspect(ZIKA_PROTEOME)
    assert r["n_entries"] == 2
    assert r["accessions"] == ["Q32ZE1", "P0DXN9"]
    assert r["chains"]["Q32ZE1"] == 13
    assert r["total_chains"] == 13


def test_inspect_accepts_an_uncompressed_body():
    """compressed=true is a request, not a guarantee."""
    assert fr.inspect(json.dumps({"results": [_entry("X", [])]}).encode())


@pytest.mark.parametrize("body,because", [
    (b"<html>502 Bad Gateway</html>", "an HTML error page"),
    (gzip.compress(b'{"messages": ["Invalid query"]}'), "a JSON error object"),
    (gzip.compress(b'{"results": []}'), "a query that matched nothing"),
    (gzip.compress(b"not json at all"), "a truncated body"),
])
def test_inspect_refuses_everything_load_proteome_would_read_as_empty(body, because):
    """load_proteome does `.get("results", [])`, so each of these ingests as
    "a virus with no proteins" -- indistinguishable from a missing file, which
    is the confusion that let two orphan taxa into v0.1."""
    with pytest.raises(ValueError):
        fr.inspect(body)


# --------------------------------------------------------------------------
# the register's claims, made checkable
# --------------------------------------------------------------------------
def test_a_clean_zika_proteome_satisfies_the_register(register):
    """Zika's chain_source pins Q32ZE1 at 13 chains, verified 2026-09-17.
    Re-read from UniProt on 2026-10-11: still 13. The claim holds."""
    z = _virus(register, "Zika")
    assert z["chain_source"]["accession"] == "Q32ZE1"
    assert z["chain_source"]["chains"] == 13
    r = fr.inspect(ZIKA_PROTEOME)
    assert fr.check_against_register(r, z, is_extra=False) == []


def test_a_changed_chain_count_is_caught_not_trusted(register):
    """`verified: <date>` was a human claim nothing re-checked. UniProt
    re-annotates; a chain count that moved silently changes which mature
    peptide every ChEMBL row lands on."""
    z = _virus(register, "Zika")
    thinner = _body([_entry("Q32ZE1", ZIKA_CHAINS[:9])])
    problems = fr.check_against_register(fr.inspect(thinner), z, is_extra=False)
    assert len(problems) == 1
    assert "9 Chain features" in problems[0] and "claims 13" in problems[0]
    assert "2026-09-17" in problems[0], "say which claim is being contradicted"


def test_denv2_proteome_alone_is_refused_and_names_the_route(register):
    """The live case. P29990 has NO proteome cross-reference, so UP000095876
    cannot contain it. Writing the proteome alone yields proteins with no
    mature peptides -- which is exactly how DENV-2 would have looked had the
    file been present during the promotion."""
    d = _virus(register, "DENV-2")
    problems = fr.check_against_register(
        fr.inspect(DENV2_PROTEOME), d, is_extra=False)
    assert problems, "a proteome without the pinned accession must not pass"
    joined = " ".join(problems)
    assert "P29990" in joined
    assert "extra_accessions" in joined, "name the documented route"


def test_denv2_extra_file_satisfies_the_extra_check(register):
    d = _virus(register, "DENV-2")
    assert fr.check_against_register(
        fr.inspect(DENV2_EXTRA), d, is_extra=True) == []


def test_an_extra_file_missing_its_accession_is_caught(register):
    d = _virus(register, "DENV-2")
    wrong = _body([_entry("Q32ZE1", ZIKA_CHAINS)])
    problems = fr.check_against_register(fr.inspect(wrong), d, is_extra=True)
    assert problems and "P29990" in problems[0]


def test_the_mers_trap_a_chainless_proteome_is_refused():
    """A TrEMBL-only species proteome downloads cleanly and carries zero Chain
    features. Every ChEMBL row then lands on the polyprotein, and nothing
    reports it. This is the shape MERS is in, which is why MERS has ten
    extra_accessions."""
    virus = {"name": "Hypothetical", "taxon": "1", "proteome_id": "UP0",
             "chain_source": {"accession": "X1", "taxon": "1", "chains": 0}}
    problems = fr.check_against_register(
        fr.inspect(_body([_entry("X1", [])])), virus, is_extra=False)
    assert problems and "zero Chain features" in problems[0]
    assert "polyprotein" in problems[0]


def test_a_virus_declaring_extra_accessions_is_exempt_from_the_chain_floor(register):
    """MERS's own proteome is chainless BY DESIGN -- the chains come from the
    curated accessions, which REPLACE it. Flagging that would be noise on
    every run, so the floor only applies where no curated route is declared."""
    m = _virus(register, "MERS-CoV")
    assert m["extra_accessions"], "precondition: MERS declares curated accessions"
    chainless = _body([_entry(m["chain_source"]["accession"], [])])
    problems = fr.check_against_register(fr.inspect(chainless), m, is_extra=False)
    assert not any("zero Chain" in p for p in problems)


# --------------------------------------------------------------------------
# the behaviour that protects the release
# --------------------------------------------------------------------------
def test_a_refused_file_is_never_written(register, tmp_path, monkeypatch):
    """The whole point. ingest_spine consumes whatever is on disk without
    re-validating it, so a file that disagrees with the register is worse than
    no file -- absence at least records `proteome_missing`."""
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    d = _virus(register, "DENV-2")
    problems = fr.fetch_virus(
        d, force=True, dry_run=False,
        fetch=lambda url: DENV2_EXTRA if "accession" in url else DENV2_PROTEOME)
    assert problems, "the proteome disagreed and must have been refused"
    assert not (tmp_path / "UP000095876.json.gz").exists()
    # The EXTRA file was fine and is written: a refusal is per-file.
    assert (tmp_path / "EXTRA_11060.json.gz").exists()


def test_a_clean_fetch_writes_the_file_and_its_provenance(register, tmp_path,
                                                          monkeypatch):
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    z = _virus(register, "Zika")
    assert fr.fetch_virus(z, force=True, dry_run=False,
                          fetch=lambda url: ZIKA_PROTEOME) == 0
    target = tmp_path / "UP000054557.json.gz"
    assert target.read_bytes() == ZIKA_PROTEOME
    prov = json.loads(fr.sidecar_path(target).read_text())
    assert prov["url"] == fr.proteome_url("UP000054557")
    assert prov["chains"]["Q32ZE1"] == 13
    assert len(prov["sha256"]) == 64
    assert prov["retrieved"].count("-") == 2


def test_dry_run_fetches_nothing(register, tmp_path, monkeypatch):
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)

    def _explode(url):
        raise AssertionError("dry run must not fetch")

    fr.fetch_virus(_virus(register, "Zika"), force=True, dry_run=True,
                   fetch=_explode)
    assert not list(tmp_path.iterdir())


def test_an_existing_file_is_not_clobbered_without_force(register, tmp_path,
                                                         monkeypatch):
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    target = tmp_path / "UP000054557.json.gz"
    target.write_bytes(b"precious hand-placed bytes")
    fr.fetch_virus(_virus(register, "Zika"), force=False, dry_run=False,
                   fetch=lambda url: ZIKA_PROTEOME)
    assert target.read_bytes() == b"precious hand-placed bytes"


def test_check_reports_a_hand_placed_file_as_lacking_provenance(
        register, tmp_path, monkeypatch, capsys):
    """How the current data/raw gets audited: every file there was placed by
    hand and none has a URL, date or checksum."""
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    (tmp_path / "UP000054557.json.gz").write_bytes(ZIKA_PROTEOME)
    fr.check_on_disk(register)
    out = capsys.readouterr().out
    assert "NO PROVENANCE" in out
    assert "UP000054557.json.gz" in out


def test_check_needs_no_network(register, tmp_path, monkeypatch):
    """--check must be safe to run in CI and offline."""
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)

    def _explode(url):
        raise AssertionError("--check must not reach the network")

    monkeypatch.setattr(fr, "fetch_bytes", _explode)
    fr.check_on_disk(register)


def test_check_flags_a_missing_proteome_for_a_built_virus_only(
        register, tmp_path, monkeypatch, capsys):
    """A surveyed virus has no business having its proteome on disk; a built
    one that lacks it is the DENV-2/Zika promotion failure, reported before
    the build rather than after."""
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    problems = fr.check_on_disk(register)
    out = capsys.readouterr().out

    # Derived, not a literal: every FILE a built virus needs, which is one
    # proteome each plus an EXTRA wherever curated accessions are declared.
    # MERS is the reason that is 8 and not 7 -- its own proteome is chainless
    # by design, so a missing EXTRA_1335626.json.gz is a real missing input,
    # not a formality. A hand-written 7 here would have asserted the opposite.
    expected = sum(len(fr.targets_for(v)) for v in register["viruses"])
    assert problems == expected, f"every built input is missing ({expected})"

    for v in register["viruses"]:
        assert v["proteome_id"] in out, f"{v['name']} not reported"

    # Surveyed viruses are absent too, and that is not a problem to report:
    # Zika's proteome is MEANT not to be on disk while Zika is surveyed.
    assert "UP000054557" not in out, \
        "a surveyed virus's missing proteome must not be reported as MISSING"
