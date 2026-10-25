"""Provenance for files that were placed by hand.

`make data-timeline` reports the host layer as five files, 328,846,356
bytes, `mtime only` -- the human proteome and STRING v12.0, feeding over
three-quarters of the graph with no recorded URL, release or checksum.
fetch_raw.py stamps every proteome it fetches; these predate it.
"""
from __future__ import annotations

import io
import json
import sys
from contextlib import redirect_stdout
from importlib import import_module
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

rp = import_module("record_raw_provenance")


def _dir(tmp: Path, **files: bytes) -> Path:
    d = tmp / "raw"
    d.mkdir(parents=True, exist_ok=True)
    for name, body in files.items():
        (d / name.replace("__", ".")).write_bytes(body)
    return d


def _run(d: Path, *extra) -> tuple[int, str]:
    buf = io.StringIO()
    with redirect_stdout(buf):
        code = rp.main(["--dir", str(d), *extra])
    return code, buf.getvalue()


# --------------------------------------------- it does not invent a download
def test_the_url_and_date_are_null_not_guessed(tmp_path):
    """A plausible URL written into a provenance record is worse than
    none: a reader cannot tell it from one that was really fetched."""
    d = _dir(tmp_path, human_proteome__json__gz=b"x")
    _run(d)
    rec = json.loads((d / "human_proteome.json.gz.source.json").read_text())
    assert rec["url"] is None
    assert rec["retrieved"] is None
    assert rec["provenance"] == "reconstructed"


def test_the_record_says_when_it_was_reconstructed(tmp_path):
    d = _dir(tmp_path, human_proteome__json__gz=b"x")
    _run(d)
    rec = json.loads((d / "human_proteome.json.gz.source.json").read_text())
    assert len(rec["recorded"]) == 10 and rec["recorded"][4] == "-"
    assert "not recoverable" in rec["note"]


def test_a_checksum_is_what_it_actually_establishes(tmp_path):
    import hashlib
    d = _dir(tmp_path, human_proteome__json__gz=b"hello")
    _run(d)
    rec = json.loads((d / "human_proteome.json.gz.source.json").read_text())
    assert rec["sha256"] == hashlib.sha256(b"hello").hexdigest()
    assert rec["bytes"] == 5


# ------------------------------------------ a version the file already states
def test_strings_release_is_read_from_its_filename(tmp_path):
    """STRING writes its release into every filename. Reporting that as
    unknown sends someone to establish a fact the file already states."""
    d = _dir(tmp_path, **{"9606.protein.links.detailed.v12.0.txt.gz": b"x"})
    _run(d)
    rec = json.loads(
        (d / "9606.protein.links.detailed.v12.0.txt.gz.source.json").read_text())
    assert rec["version_in_filename"] == "v12.0"


def test_a_filename_with_no_version_says_null_not_a_guess(tmp_path):
    d = _dir(tmp_path, human_proteome__json__gz=b"x")
    _run(d)
    rec = json.loads((d / "human_proteome.json.gz.source.json").read_text())
    assert rec["version_in_filename"] is None


def test_data_timeline_reads_the_same_versions():
    """Two regexes for one question drift. data_timeline reported `-` for
    the host layer while the version sat in the filename."""
    dt = import_module("data_timeline")
    for name, want in (
            ("9606.protein.links.detailed.v12.0.txt.gz", "v12.0"),
            ("data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db",
             "chembl_37")):
        m = dt.VERSION_IN_PATH.search(name)
        assert m and m.group(1) == want, name


# ------------------------------------- a real record is never downgraded
def test_an_existing_sidecar_is_kept(tmp_path):
    """fetch_raw's record carries a real URL and retrieval date.
    Overwriting it with a reconstruction destroys the one thing that
    distinguishes them."""
    d = _dir(tmp_path, human_proteome__json__gz=b"x")
    real = {"url": "https://example.invalid/p.gz", "retrieved": "2026-01-01",
            "sha256": "deadbeef", "fetched_by": "scripts/fetch_raw.py"}
    (d / "human_proteome.json.gz.source.json").write_text(json.dumps(real))
    _, out = _run(d)
    rec = json.loads((d / "human_proteome.json.gz.source.json").read_text())
    assert rec == real
    assert "already stamped" in out


def test_a_sidecar_is_not_itself_stamped(tmp_path):
    d = _dir(tmp_path, human_proteome__json__gz=b"x")
    _run(d)
    _run(d)
    assert not (d / "human_proteome.json.gz.source.json.source.json").exists()


# -------------------------------------------------- drift becomes visible
def test_a_changed_file_is_reported_as_drifted(tmp_path):
    """THE POINT OF THE CHECKSUM. Today one of these could be replaced
    with a different STRING release and nothing in the repo would notice."""
    d = _dir(tmp_path, human_proteome__json__gz=b"original")
    _run(d)
    (d / "human_proteome.json.gz").write_bytes(b"a different release")
    code, out = _run(d, "--check")
    assert code == 1
    assert "DRIFTED" in out


def test_an_unchanged_file_verifies(tmp_path):
    d = _dir(tmp_path, human_proteome__json__gz=b"original")
    _run(d)
    code, out = _run(d, "--check")
    assert code == 0 and "ok" in out and "DRIFTED" not in out


def test_a_file_with_no_sidecar_is_unstamped_not_silently_ok(tmp_path):
    """Returning 0 for a file nothing has ever recorded would report an
    unauditable directory as audited."""
    d = _dir(tmp_path, human_proteome__json__gz=b"x")
    code, out = _run(d, "--check")
    assert code == 1 and "UNSTAMPED" in out


def test_check_writes_nothing(tmp_path):
    d = _dir(tmp_path, human_proteome__json__gz=b"x")
    _run(d, "--check")
    assert not (d / "human_proteome.json.gz.source.json").exists()


def test_an_empty_or_missing_directory_is_an_error(tmp_path):
    assert rp.main(["--dir", str(tmp_path / "nope")]) == 2
    (tmp_path / "empty").mkdir()
    assert rp.main(["--dir", str(tmp_path / "empty")]) == 2


def test_the_consumer_is_recorded(tmp_path):
    """A sidecar that says only what a file hashes to does not say why the
    file is here."""
    d = _dir(tmp_path, human_proteome__json__gz=b"x")
    _run(d)
    rec = json.loads((d / "human_proteome.json.gz.source.json").read_text())
    assert rec["consumed_by"] == "scripts/ingest_host.py"
