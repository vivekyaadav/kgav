"""Gate for the virus-host PPI fetcher.

The files it replaces were placed by hand from a host that has since
moved, with nothing recording where they came from and nothing to re-run.
Nine of sixteen built viruses had no virus-host PPI at all because no
flavivirus file was ever fetched.

What these fix: that the query reaches a virus curated under a strain,
that an empty answer and a failed one stay distinguishable, and that a
file already on disk is not silently replaced.
"""
from __future__ import annotations

import json
import re
import sys
import urllib.parse
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import fetch_vhppi as F

BASE = "https://example.invalid/psicquic"
SRC = "infores:virhostnet"


def _entry(**kw):
    e = {"taxon": "11060", "name": "DENV-2", "family": "Flaviviridae"}
    e.update(kw)
    return e


# ------------------------------------------------------------- the query
def test_an_isolate_taxon_is_queried_as_well_as_the_species():
    """Seven of nine flaviviruses are curated under a strain, and DENV-2's
    chains live under 31634. Querying the species alone returns nothing
    for exactly the viruses this script exists to reach, and the emptiness
    reads as an absent interactome rather than a wrong query.
    """
    t = F.taxids_for(_entry(isolate_taxa=["31634"]))
    assert t == ["11060", "31634"]


def test_the_chain_source_taxon_is_queried_too():
    t = F.taxids_for(_entry(chain_source={"taxon": "31634"}))
    assert "31634" in t


def test_taxids_are_deduplicated_keeping_the_species_first():
    t = F.taxids_for(_entry(isolate_taxa=["31634", "11060"],
                            chain_source={"taxon": "31634"}))
    assert t == ["11060", "31634"]


def test_the_query_covers_both_sides_of_the_interaction():
    """A virus protein can be either partner in a MITAB row."""
    q = F.miql(["11060", "31634"])
    for t in ("11060", "31634"):
        assert f"taxidA:{t}" in q and f"taxidB:{t}" in q


def test_each_taxid_clause_is_bracketed():
    """Without brackets the trailing OR binds across taxids and the query
    means something else."""
    q = F.miql(["1", "2"])
    assert q == "(taxidA:1 OR taxidB:1) OR (taxidA:2 OR taxidB:2)"


def test_the_query_is_url_encoded():
    u = F.url_for(BASE, ["11060"], "count")
    assert " " not in u and "(" not in u
    assert u.startswith(BASE + "/webservices/current/search/query/")
    assert u.endswith("?format=count")


# ------------------------------------------------------------- filenames
@pytest.mark.parametrize("name,stem", [
    ("SARS-CoV-2", "sarscov2"), ("HCoV-229E", "hcov229e"),
    ("West Nile", "westnile"), ("DENV-2", "denv2"),
    ("Yellow fever", "yellowfever"),
])
def test_the_filename_stem_is_stable(name, stem):
    assert F.slug(name) == stem


def test_every_register_virus_gets_a_distinct_filename():
    """Two viruses sharing a stem would overwrite each other, and
    ingest_vhppi globs the directory."""
    viruses = yaml.safe_load((ROOT / "config/viruses.yaml").read_text())["viruses"]
    stems = [F.slug(v["name"]) for v in viruses]
    assert len(stems) == len(set(stems)), [s for s in stems
                                           if stems.count(s) > 1]


# --------------------------------------------- empty, failed, and present
def test_an_empty_answer_is_not_a_failure(tmp_path, monkeypatch, capsys):
    """Zero interactions is this service reporting nothing, which is a
    source limit. Exit 0, nothing written, and said out loud."""
    monkeypatch.setattr(F, "fetch", lambda url: b"0")
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    rc = F.main(["--base", BASE, "--source", SRC, "--config", str(cfg),
                 "--out", str(tmp_path / "raw")])
    out = capsys.readouterr().out
    assert rc == 0
    assert "returned NO interaction" in out and "source limit" in out
    assert not list((tmp_path / "raw").glob("*.tab27"))


def test_a_failed_query_is_not_an_empty_result(tmp_path, monkeypatch, capsys):
    """The distinction this whole layer keeps getting wrong. A failure
    leaves the gap open and must not be reported as 'no data'."""
    def boom(url):
        raise OSError("connection reset")
    monkeypatch.setattr(F, "fetch", boom)
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    rc = F.main(["--base", BASE, "--source", SRC, "--config", str(cfg),
                 "--out", str(tmp_path / "raw")])
    out = capsys.readouterr().out
    assert rc == 1
    assert "FAILED to query" in out
    assert "A failure is not an empty result" in out
    assert "returned NO interaction" not in out


def test_a_dry_run_writes_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr(F, "fetch", lambda url: b"42")
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    F.main(["--base", BASE, "--source", SRC, "--config", str(cfg), "--dry-run",
            "--out", str(tmp_path / "raw")])
    assert not (tmp_path / "raw").exists() or not list(
        (tmp_path / "raw").glob("*"))


def test_an_existing_file_is_left_alone(tmp_path, monkeypatch, capsys):
    """A re-run must not silently replace a file whose provenance someone
    may have recorded by hand."""
    monkeypatch.setattr(F, "fetch", lambda url: b"2" if "count" in url
                        else b"a\tb\n")
    raw = tmp_path / "raw"
    raw.mkdir()
    existing = raw / "vhn_denv2.tab27"
    existing.write_bytes(b"ORIGINAL")
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    F.main(["--base", BASE, "--source", SRC, "--config", str(cfg), "--out", str(raw)])
    assert existing.read_bytes() == b"ORIGINAL"
    assert "exists, left alone" in capsys.readouterr().out


def test_overwrite_replaces_it(tmp_path, monkeypatch):
    monkeypatch.setattr(F, "fetch", lambda url: b"1" if "count" in url
                        else b"new\trow\n")
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "vhn_denv2.tab27").write_bytes(b"ORIGINAL")
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    F.main(["--base", BASE, "--source", SRC, "--config", str(cfg), "--out", str(raw),
            "--overwrite"])
    assert (raw / "vhn_denv2.tab27").read_bytes() == b"new\trow\n"


# ------------------------------------------------------------ provenance
def test_the_sidecar_records_what_was_asked_and_what_came_back(
        tmp_path, monkeypatch):
    monkeypatch.setattr(F, "fetch", lambda url: b"2" if "count" in url
                        else b"r1\nr2\n")
    raw = tmp_path / "raw"
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry(isolate_taxa=["31634"])]}))
    F.main(["--base", BASE, "--source", SRC, "--config", str(cfg), "--out", str(raw)])
    side = json.loads((raw / "vhn_denv2.tab27.source.json").read_text())
    assert side["taxids_queried"] == ["11060", "31634"]
    assert side["n_rows"] == 2
    assert side["bytes"] == 6
    assert side["virus"] == "DENV-2" and side["taxon"] == "11060"
    assert side["fetched_by"] == "scripts/fetch_vhppi.py"
    assert side["consumed_by"] == "scripts/ingest_vhppi.py"
    assert side["url"].startswith(BASE)
    assert len(side["sha256"]) == 64


def test_the_sidecar_sits_beside_the_file_like_every_other(
        tmp_path, monkeypatch):
    """fetch_raw uses <file>.source.json and check-raw looks for it."""
    monkeypatch.setattr(F, "fetch", lambda url: b"1" if "count" in url
                        else b"r\n")
    raw = tmp_path / "raw"
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    F.main(["--base", BASE, "--source", SRC, "--config", str(cfg), "--out", str(raw)])
    assert (raw / "vhn_denv2.tab27.source.json").exists()


def test_a_row_count_differing_from_the_interaction_count_is_flagged(
        tmp_path, monkeypatch, capsys):
    """PSICQUIC counts interactions and a MITAB line is a binary pair, so
    the two can legitimately differ -- but a truncated download looks
    exactly like a small one, so the mismatch is said rather than hidden.
    """
    monkeypatch.setattr(F, "fetch", lambda url: b"50" if "count" in url
                        else b"only\tone\n")
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    F.main(["--base", BASE, "--source", SRC, "--config", str(cfg),
            "--out", str(tmp_path / "raw")])
    assert "against a count of" in capsys.readouterr().out


# ----------------------------------------------------------------- scope
def test_family_selects_only_that_family(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(F, "fetch", lambda url: b"0")
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [
        _entry(),
        {"taxon": "2697049", "name": "SARS-CoV-2", "family": "Coronaviridae"},
    ]}))
    F.main(["--base", BASE, "--source", SRC, "--config", str(cfg), "--dry-run",
            "--family", "Flaviviridae", "--out", str(tmp_path / "raw")])
    out = capsys.readouterr().out
    assert "DENV-2" in out and "SARS-CoV-2" not in out


def test_a_family_matching_nothing_exits_two(tmp_path, capsys):
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    rc = F.main(["--base", BASE, "--source", SRC, "--config", str(cfg),
                 "--family", "Nope", "--out", str(tmp_path / "raw")])
    assert rc == 2
    assert "no virus in the register matches" in capsys.readouterr().out


def test_the_base_url_is_required():
    """A URL defaulted from memory is the assumption this project refuses,
    and a stale one fails in a way that looks like an empty database."""
    with pytest.raises(SystemExit):
        F.main(["--config", str(ROOT / "config/viruses.yaml"), "--dry-run"])


# ------------------------------- the instruction the docstring gives out
def test_the_docstring_asks_the_registry_for_the_rest_endpoint():
    """A doc that hands the reader the wrong URL is worse than no doc.

    Without protocol=rest the PSICQUIC registry returns the SOAP endpoint.
    This docstring omitted it, the SOAP endpoint was passed to --base, and
    it answered every REST query with HTTP 500 and a soap:Fault. That
    reads as a broken service rather than a wrong URL, and nearly moved
    this layer to a different data source.
    """
    doc = F.__doc__
    assert "protocol=rest" in doc
    assert "SOAP" in doc, "the trap has to be named, not just avoided"


def test_no_url_in_the_docstring_is_the_soap_shaped_base():
    """The SOAP endpoint ENDS /webservices/psicquic; the REST base does
    not. Showing the former as a --base example is the same mistake again.

    A plain substring test would fail on the EBI registry URL, which
    carries that segment in the middle quite legitimately -- so this
    checks where the segment sits, not whether it appears.
    """
    urls = re.findall(r"https?://[^\s\\]+", F.__doc__)
    assert urls, "the docstring is supposed to show the registry call"
    soap = [u for u in urls if u.rstrip("/").endswith("/webservices/psicquic")]
    assert not soap, soap


def test_the_query_path_is_appended_to_a_base_without_webservices():
    """--base takes everything before /webservices; QUERY_PATH supplies
    it. Both carrying it doubles the segment and the server 500s."""
    assert not F.QUERY_PATH.startswith("/")
    u = F.url_for("http://h:9090/psicquic", ["11060"], "count")
    assert u.count("/webservices/") == 1


# ----------------------------------------- the service is named, not guessed
def test_the_source_is_required():
    """Same reason --base is. The file is labelled with this, and a
    service inferred from a filename is a false provenance claim that
    the tier ablation then reads as fact."""
    with pytest.raises(SystemExit):
        F.main(["--base", BASE, "--config", str(ROOT / "config/viruses.yaml"),
                "--dry-run"])


@pytest.mark.parametrize("source,pfx", [
    ("infores:virhostnet", "vhn"),
    ("infores:intact", "intact"),
])
def test_the_two_services_already_on_disk_keep_their_prefixes(source, pfx):
    """Renaming them would orphan the sidecars beside them, which are
    the provenance those files have."""
    assert F.prefix_for(source) == pfx


@pytest.mark.parametrize("source,pfx", [
    ("infores:mint", "mint"),
    ("infores:biogrid", "biogrid"),
    ("infores:matrix-db", "matrixdb"),
])
def test_a_new_service_derives_a_prefix_of_its_own(source, pfx):
    assert F.prefix_for(source) == pfx


def test_two_services_do_not_write_the_same_filename():
    """A MINT fetch overwriting a VirHostNet file would destroy the
    source it replaced and keep its sidecar."""
    names = {F.prefix_for(s) for s in
             ("infores:virhostnet", "infores:intact", "infores:mint",
              "infores:biogrid", "infores:dip")}
    assert len(names) == 5


def test_the_file_is_named_for_its_service(tmp_path, monkeypatch):
    monkeypatch.setattr(F, "fetch", lambda url: b"1" if "count" in url
                        else b"r\n")
    raw = tmp_path / "raw"
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    F.main(["--base", BASE, "--source", "infores:mint",
            "--config", str(cfg), "--out", str(raw)])
    assert (raw / "mint_denv2.tab27").exists()
    assert not (raw / "vhn_denv2.tab27").exists()


def test_the_sidecar_records_the_service(tmp_path, monkeypatch):
    monkeypatch.setattr(F, "fetch", lambda url: b"1" if "count" in url
                        else b"r\n")
    raw = tmp_path / "raw"
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    F.main(["--base", BASE, "--source", "infores:mint",
            "--config", str(cfg), "--out", str(raw)])
    side = json.loads((raw / "mint_denv2.tab27.source.json").read_text())
    assert side["source"] == "infores:mint"


def test_the_ingest_reads_back_what_the_fetcher_wrote(tmp_path, monkeypatch):
    """The two halves have to agree, or the label is decided twice."""
    from kgav.vhppi import file_source
    monkeypatch.setattr(F, "fetch", lambda url: b"1" if "count" in url
                        else b"r\n")
    raw = tmp_path / "raw"
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    F.main(["--base", BASE, "--source", "infores:mint",
            "--config", str(cfg), "--out", str(raw)])
    assert file_source(raw / "mint_denv2.tab27") == "infores:mint"


def test_family_coronaviridae_selects_the_coronaviruses(tmp_path, monkeypatch,
                                                        capsys):
    """It matched NOTHING and exited 2. Coronaviruses inherit their
    family from the register's defaults block, and the filter read the
    per-entry key -- so every --family run this week covered
    flaviviruses only, and a coronavirus fetch reported the register as
    having none."""
    monkeypatch.setattr(F, "fetch", lambda url: b"0")
    rc = F.main(["--base", BASE, "--source", SRC, "--dry-run",
                 "--config", str(ROOT / "config/viruses.yaml"),
                 "--family", "Coronaviridae", "--out", str(tmp_path)])
    out = capsys.readouterr().out
    assert rc == 0, out
    assert "7 virus(es)" in out
    assert "SARS-CoV-2" in out and "DENV-2" not in out


def test_family_flaviviridae_still_selects_only_those(tmp_path, monkeypatch,
                                                      capsys):
    monkeypatch.setattr(F, "fetch", lambda url: b"0")
    F.main(["--base", BASE, "--source", SRC, "--dry-run",
            "--config", str(ROOT / "config/viruses.yaml"),
            "--family", "Flaviviridae", "--out", str(tmp_path)])
    out = capsys.readouterr().out
    assert "9 virus(es)" in out
    assert "DENV-2" in out and "SARS-CoV-2" not in out


# ------------------------------- the base URL, as the registry prints it
@pytest.mark.parametrize("given", [
    "http://www.ebi.ac.uk/Tools/webservices/psicquic/hpidb/webservices/current/search/",
    "http://www.ebi.ac.uk/Tools/webservices/psicquic/hpidb/webservices/current/search",
    "http://www.ebi.ac.uk/Tools/webservices/psicquic/hpidb",
])
def test_an_ebi_service_normalises_to_one_base(given):
    """`webservices` appears TWICE in an EBI restUrl. The rule the
    docstring gave -- everything before /webservices -- yields
    http://www.ebi.ac.uk/Tools and every query 404s. Eleven of the
    fourteen active services are EBI-hosted."""
    assert F.normalise_base(given) == \
        "http://www.ebi.ac.uk/Tools/webservices/psicquic/hpidb"


def test_virhostnet_still_normalises_to_what_worked():
    """The one service already in use must not move."""
    VHN = "http://virhostnet.prabi.fr:9090/psicquic"
    for given in (VHN, VHN + "/webservices/current/search/"):
        assert F.normalise_base(given) == VHN


def test_a_registry_url_builds_a_query_with_one_query_path():
    """Doubling the segment is what a hand-trimmed base gets wrong."""
    u = F.url_for("https://psicquic.all.innatedb.com/webservices/current/"
                  "search/", ["11053"], "count")
    # COUNTED WITHOUT THE LEADING SLASH. With it, str.count misses the
    # doubling entirely: the repeated segment shares its leading slash
    # with the end of the first match, counting is non-overlapping, and
    # ".../search/webservices/current/search/" scores 1. The assertion
    # read as a duplication check and detected nothing.
    assert u.count("webservices/current/search") == 1, u
    assert u.startswith("https://psicquic.all.innatedb.com/webservices/")


def test_a_base_that_is_already_trimmed_is_untouched():
    assert F.normalise_base("https://example.invalid/psicquic") == \
        "https://example.invalid/psicquic"


def test_the_docstring_does_not_give_the_before_webservices_rule():
    """It is wrong for most of the registry, and a doc that hands the
    reader a broken URL is the SOAP mistake again."""
    doc = F.__doc__
    assert "everything before /webservices" not in doc
    assert "normalise_base" in doc


# --------------------------------------------- paging a capped service
def _pager(total: int, cap: int = 1000):
    """A PSICQUIC service holding `total` rows and capping at `cap`."""
    calls = []

    def fetch(url):
        calls.append(url)
        q = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
        first = int(q.get("firstResult", ["0"])[0])
        want = min(int(q.get("maxResults", [str(cap)])[0]), cap)
        rows = [f"row{i}".encode() for i in range(first, min(first + want,
                                                            total))]
        return b"\n".join(rows) + (b"\n" if rows else b"")
    return fetch, calls


def test_a_capped_service_still_yields_every_row():
    """IMEx caps at 1,000. Asked for SARS-CoV-2's 8,059 it returned
    exactly 1,000, and SARS-CoV's 2,291 likewise: 12,496 interactions
    arrived as 4,146 rows. Every virus under the cap came back whole,
    so the truncation was invisible except in the row-vs-count note."""
    fetch, calls = _pager(8059)
    raw, pages = F.fetch_all(BASE, ["2697049"], "tab25", 8059, fetch=fetch)
    assert len([ln for ln in raw.split(b"\n") if ln.strip()]) == 8059
    assert pages == 9 and len(calls) == 9


def test_the_first_page_asks_from_zero():
    fetch, calls = _pager(10)
    F.fetch_all(BASE, ["1"], "tab25", 10, fetch=fetch)
    assert "firstResult=0" in calls[0] and "maxResults=" in calls[0]


def test_a_short_page_ends_the_loop():
    """A service free to return fewer than asked must not be polled
    forever."""
    fetch, _calls = _pager(1500)
    _raw, pages = F.fetch_all(BASE, ["1"], "tab25", 0, fetch=fetch)
    assert pages == 2


def test_an_under_cap_virus_is_one_request():
    """Most viruses here are well under the cap, and a second request
    per virus against 16 viruses and several services is not free."""
    fetch, _calls = _pager(53)
    _raw, pages = F.fetch_all(BASE, ["11060"], "tab25", 53, fetch=fetch)
    assert pages == 1


def test_a_server_ignoring_first_result_does_not_loop_forever(capsys):
    """It would return page one until the end of time. The run says so
    instead of spinning."""
    def stuck(url):
        return b"\n".join(b"same" for _ in range(1000)) + b"\n"
    _raw, pages = F.fetch_all(BASE, ["1"], "tab25", 10_000_000, fetch=stuck)
    assert pages == F.MAX_PAGES
    assert "ignoring firstResult" in capsys.readouterr().out


def test_the_url_carries_no_paging_when_not_asked():
    """The count query and anything else built by url_for stay as they
    were, so an unpaged caller is unaffected."""
    u = F.url_for(BASE, ["11060"], "count")
    assert "firstResult" not in u and "maxResults" not in u


def test_rows_fetched_across_pages_reach_the_file(tmp_path, monkeypatch):
    """End to end: the written file holds every page, not the first."""
    fetch, _calls = _pager(2500)
    monkeypatch.setattr(F, "fetch",
                        lambda url: b"2500" if "count" in url else fetch(url))
    raw = tmp_path / "raw"
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    F.main(["--base", BASE, "--source", SRC, "--config", str(cfg),
            "--out", str(raw)])
    side = json.loads((raw / "vhn_denv2.tab27.source.json").read_text())
    assert side["n_rows"] == 2500


def test_a_total_that_lands_on_a_page_boundary_costs_no_extra_request():
    """`expected` is an optimisation, not a correctness guard: without
    it the loop still terminates on the empty page after the last full
    one. With it, a virus whose count is an exact multiple of the page
    size costs no spare round trip -- and across 16 viruses and several
    services the spare requests add up."""
    fetch, _calls = _pager(2000)
    raw, pages = F.fetch_all(BASE, ["1"], "tab25", 2000, fetch=fetch)
    assert len([ln for ln in raw.split(b"\n") if ln.strip()]) == 2000
    assert pages == 2, "the empty third page should not be requested"
    # and with no count to go on, it costs exactly one more
    fetch2, _ = _pager(2000)
    _raw2, pages2 = F.fetch_all(BASE, ["1"], "tab25", 0, fetch=fetch2)
    assert pages2 == 3
