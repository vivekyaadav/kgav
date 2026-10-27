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
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import fetch_vhppi as F

BASE = "https://example.invalid/psicquic"


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
    rc = F.main(["--base", BASE, "--config", str(cfg),
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
    rc = F.main(["--base", BASE, "--config", str(cfg),
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
    F.main(["--base", BASE, "--config", str(cfg), "--dry-run",
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
    F.main(["--base", BASE, "--config", str(cfg), "--out", str(raw)])
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
    F.main(["--base", BASE, "--config", str(cfg), "--out", str(raw),
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
    F.main(["--base", BASE, "--config", str(cfg), "--out", str(raw)])
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
    F.main(["--base", BASE, "--config", str(cfg), "--out", str(raw)])
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
    F.main(["--base", BASE, "--config", str(cfg),
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
    F.main(["--base", BASE, "--config", str(cfg), "--dry-run",
            "--family", "Flaviviridae", "--out", str(tmp_path / "raw")])
    out = capsys.readouterr().out
    assert "DENV-2" in out and "SARS-CoV-2" not in out


def test_a_family_matching_nothing_exits_two(tmp_path, capsys):
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [_entry()]}))
    rc = F.main(["--base", BASE, "--config", str(cfg),
                 "--family", "Nope", "--out", str(tmp_path / "raw")])
    assert rc == 2
    assert "no virus in the register matches" in capsys.readouterr().out


def test_the_base_url_is_required():
    """A URL defaulted from memory is the assumption this project refuses,
    and a stale one fails in a way that looks like an empty database."""
    with pytest.raises(SystemExit):
        F.main(["--config", str(ROOT / "config/viruses.yaml"), "--dry-run"])
