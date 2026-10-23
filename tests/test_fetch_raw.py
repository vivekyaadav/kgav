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


def _gz(payload: bytes) -> bytes:
    """gzip with mtime=0, because the default is time.time().

    This broke the suite on 32 workers. gzip stamps an mtime into its header,
    so a fixture built at module import produced different BYTES on workers
    that started either side of a second boundary (0x93 vs 0x94). Those bytes
    were in a parametrize list and therefore in the test IDs, so xdist saw
    different tests collected on gw0 and gw10..gw31 and refused to run at all
    -- "Different tests were collected between gw0 and gw10". It passed on
    fewer workers purely because they all landed in the same second.

    Two independent fixes, because either alone is fragile: mtime=0 here makes
    the bytes deterministic, and the parametrize below passes explicit `ids`
    so no binary ever reaches a test ID again.
    """
    return gzip.compress(payload, mtime=0)


def _body(entries: list[dict]) -> bytes:
    return _gz(json.dumps({"results": entries}).encode())


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
def test_the_fixtures_are_byte_identical_in_every_process():
    """The bug that broke the suite, pinned at its cause.

    gzip writes a 4-byte mtime at offset 4 of its header, defaulting to
    time.time(). Any fixture built at import time therefore differed between
    xdist workers that started in different seconds, and because those bytes
    reached the test IDs, collection itself diverged and the run aborted.

    Asserting the header field rather than comparing two compressions in one
    process: two calls a microsecond apart agree even with the bug present, so
    that comparison would pass and prove nothing.
    """
    for fixture in (ZIKA_PROTEOME, DENV2_PROTEOME, DENV2_EXTRA):
        assert fixture[:2] == b"\x1f\x8b", "gzip magic"
        assert fixture[4:8] == b"\x00\x00\x00\x00", \
            "gzip mtime must be zeroed, or these bytes vary per process"


def test_no_test_id_in_this_module_contains_raw_bytes(request):
    """The second half of the fix. mtime=0 makes the fixtures deterministic,
    but any future parametrize over binary would reintroduce unreadable and
    potentially unstable IDs, so the ids= argument is the durable guard."""
    ids = [item.nodeid for item in request.session.items
           if item.nodeid.startswith("tests/test_fetch_raw.py")]
    assert ids, "precondition: this module's items are collected"
    for nid in ids:
        assert "\\x" not in nid, f"binary leaked into a test id: {nid}"


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
    (_gz(b'{"messages": ["Invalid query"]}'), "a JSON error object"),
    (_gz(b'{"results": []}'), "a query that matched nothing"),
    (_gz(b"not json at all"), "a truncated body"),
], ids=["html-error-page", "json-error-object", "empty-results", "truncated"])
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


def test_mers_proteome_is_not_judged_on_content_it_never_supplies(register):
    """REGRESSION. The first real --check run reported MERS-CoV as DISAGREEING:
    chain_source pins K9N7C7, and UP000112279 holds eleven unrelated
    K0B*/K4L* accessions instead.

    That is how MERS is designed. K9N7C7 is filed under isolate taxon 1263720,
    arrives via EXTRA_1335626.json.gz, and ingest_spine.py:221 does
    `entries = curated` -- the curated accessions REPLACE the proteome, so its
    content is never read. Only its existence matters, via the
    pfile.exists() guard that returns early before the curated branch.

    Checking the proteome asked the wrong file. A check that cries wolf on a
    correctly curated virus is a check that gets switched off.
    """
    m = _virus(register, "MERS-CoV")
    real_mers_proteome = _body([_entry(a, []) for a in (
        "K0BRG7", "K0BWD0", "K4LC41", "K0BU37", "K0BVN3", "K0BRH1",
        "K0BU34", "K0BUR3", "K0BVN0", "K0BWC8", "K4L9K1")])
    assert fr.check_against_register(
        fr.inspect(real_mers_proteome), m, is_extra=False) == [], \
        "a curated virus's proteome must not be judged on content it replaces"


def test_denv2_proteome_is_a_gate_token_not_a_source(register):
    """Same mechanism, the flavivirus case. P29990 has NO proteome
    cross-reference of any kind, so UP000095876 can never contain it -- and
    with extra_accessions declared, it does not need to."""
    d = _virus(register, "DENV-2")
    assert d["extra_accessions"] == ["P29990"]
    assert fr.check_against_register(
        fr.inspect(DENV2_PROTEOME), d, is_extra=False) == []


def test_the_chains_are_verified_in_the_file_that_actually_supplies_them(register):
    """The check did not disappear, it moved. DENV-2's 13 chains are claimed
    on P29990, which lives in the curated file, so that is where the count is
    compared."""
    d = _virus(register, "DENV-2")
    assert fr.check_against_register(
        fr.inspect(DENV2_EXTRA), d, is_extra=True) == []

    thinner = _body([_entry("P29990", ZIKA_CHAINS[:9])])
    problems = fr.check_against_register(fr.inspect(thinner), d, is_extra=True)
    assert len(problems) == 1
    assert "9 Chain features in the curated file" in problems[0]
    assert "claims 13" in problems[0] and "2026-09-17" in problems[0]


def test_a_virus_with_no_curated_route_is_refused_and_told_to_add_one(register):
    """What the DENV-2 case looked like BEFORE extra_accessions was set, and
    what the next promotion of an isolate-curated virus will look like: the
    pinned accession is absent and nothing declares how it arrives."""
    virus = {"name": "Hypothetical", "taxon": "1", "proteome_id": "UP0",
             "chain_source": {"accession": "P29990", "taxon": "31634",
                              "chains": 13}}
    problems = fr.check_against_register(
        fr.inspect(DENV2_PROTEOME), virus, is_extra=False)
    assert problems, "an unexplained absence must not pass"
    joined = " ".join(problems)
    assert "P29990" in joined
    assert "extra_accessions" in joined, "name the documented route"


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
    # Zika declares no extra_accessions, so its proteome IS the source and is
    # checked. A chain count that moved from 13 to 9 silently changes which
    # mature peptide every ChEMBL row lands on.
    z = _virus(register, "Zika")
    thinner = _body([_entry("Q32ZE1", ZIKA_CHAINS[:9])])
    problems = fr.fetch_virus(z, force=True, dry_run=False,
                              fetch=lambda url: thinner)
    assert problems, "the proteome contradicted the register and was refused"
    assert not (tmp_path / "UP000054557.json.gz").exists()
    assert not fr.sidecar_path(tmp_path / "UP000054557.json.gz").exists(), \
        "a refusal must not leave provenance for a file that was not written"


def test_a_refusal_is_per_file_not_per_virus(register, tmp_path, monkeypatch):
    """DENV-2 fetches two files. A bad curated file must not suppress a good
    proteome, or a retry would have to redo both."""
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    d = _virus(register, "DENV-2")
    bad_extra = _body([_entry("P29990", [])])        # zero chains
    problems = fr.fetch_virus(
        d, force=True, dry_run=False,
        fetch=lambda url: bad_extra if "accession" in url else DENV2_PROTEOME)
    assert problems
    assert (tmp_path / "UP000095876.json.gz").exists(), "the good file lands"
    assert not (tmp_path / "EXTRA_11060.json.gz").exists(), "the bad one does not"


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

    # Surveyed viruses are absent too, and that is not a problem to report.
    # This named Zika's UP000054557 until Zika was promoted on 2026-10-12, at
    # which point its proteome became a real required input and the assertion
    # inverted. Taken from the register instead, so it tracks any promotion.
    surveyed_proteomes = {v["proteome_id"] for v in register["surveyed"]}
    built_proteomes = {v["proteome_id"] for v in register["viruses"]}
    for pid in surveyed_proteomes - built_proteomes:
        assert pid not in out, (
            f"{pid} belongs to a surveyed virus; a missing proteome there is "
            f"not a problem to report")


# ------------------- a failed fetch for a file nothing reads is not a problem
#
# UniProt retired the TrEMBL-only flavivirus species proteomes on or before
# 2026-10-23. Six entries then fetched their EXTRA file correctly and their
# proteome not at all -- the intended outcome, because extra_accessions
# REPLACE the proteome in ingest_spine. The run printed:
#
#     6 problem(s)
#     Nothing was written for a refused file.
#
# after writing six files correctly. It reported a working ingest as broken,
# and used one word, "problem", for a refusal (contents disagree with the
# register -- distrust the source) and an empty result (the file is gone,
# and for these entries irrelevant).

def _gone(url):
    raise RuntimeError("'results' is empty -- the query matched nothing")


def test_a_retired_proteome_is_not_a_problem_when_extras_replace_it(
        register, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    d = _virus(register, "DENV-2")
    problems = fr.fetch_virus(
        d, force=True, dry_run=False,
        fetch=lambda url: DENV2_EXTRA if "accession" in url else _gone(url))
    assert problems == 0, "the only file ingest_spine reads was written"
    assert (tmp_path / "EXTRA_11060.json.gz").exists()
    assert not (tmp_path / "UP000095876.json.gz").exists()


def test_the_virus_is_reported_ready_not_left_on_a_failure_line(
        register, tmp_path, monkeypatch, capsys):
    """A FAILED line as the last word on a virus that is in fact ready is
    what sent the reader looking for a broken ingest."""
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    d = _virus(register, "DENV-2")
    fr.fetch_virus(d, force=True, dry_run=False,
                   fetch=lambda url: DENV2_EXTRA if "accession" in url
                   else _gone(url))
    out = capsys.readouterr().out
    assert "DENV-2 is READY" in out
    assert "extra_accessions" in out


def test_both_files_failing_is_still_a_problem(register, tmp_path,
                                               monkeypatch, capsys):
    """The replacement has to actually arrive. Not counting the proteome
    must not mean counting nothing."""
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    d = _virus(register, "DENV-2")
    problems = fr.fetch_virus(d, force=True, dry_run=False, fetch=_gone)
    assert problems == 1, "the EXTRA failure is the real one"
    out = capsys.readouterr().out
    assert "DENV-2 is NOT ready" in out


def test_a_virus_with_no_substitute_still_reports_its_failure(
        register, tmp_path, monkeypatch):
    """West Nile declares no extra_accessions, so nothing replaces its
    proteome and a failed fetch is exactly as serious as before."""
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    w = _virus(register, "West Nile")
    assert fr.fetch_virus(w, force=True, dry_run=False, fetch=_gone) == 1


def test_a_refused_file_is_still_counted_for_a_substituted_entry(
        register, tmp_path, monkeypatch, capsys):
    """REFUSED means the contents disagree with the register. That is a
    reason to distrust the source whatever else the entry has, and must not
    be swept up by the leniency for a replaced proteome."""
    monkeypatch.setattr(fr, "PROTEOMES", tmp_path)
    d = _virus(register, "DENV-2")
    bad_extra = _body([_entry("P29990", [])])        # zero chains
    problems = fr.fetch_virus(
        d, force=True, dry_run=False,
        fetch=lambda url: bad_extra if "accession" in url else _gone(url))
    assert problems >= 1
    out = capsys.readouterr().out
    assert "REFUSED" in out
    assert not (tmp_path / "EXTRA_11060.json.gz").exists()


def test_the_summary_separates_refused_from_failed():
    """Two causes had one word. REFUSED means the contents disagree with
    the register, which is a reason to distrust the source; FAILED means
    the query returned nothing, which for a replaced proteome is expected.
    Asserted against the source because main() reaches the network."""
    src = (ROOT / "scripts" / "fetch_raw.py").read_text()
    assert "A REFUSED file is one whose contents disagree" in src
    assert "because that file is never read" in src
    assert '"Nothing was written for a refused file. A proteome that "' \
        not in src, "the old conflated sentence is gone"
