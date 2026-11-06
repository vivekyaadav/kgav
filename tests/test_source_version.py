"""Which RELEASE of which database a layer was built from.

BUILT_BY.json answers "what code wrote this layer". Nothing answered "what
release of the source did it read", so the v0.1 MANIFEST records a
checksum per layer and a code commit and names no database version
anywhere. The audit on 2026-11-07 therefore could not say which GtoPdb,
ORCS, NCATS, InterPro or RDKit release any figure in docs/ came from.

GtoPdb is the first one fixed because it is the one where the answer was
already in hand and thrown away: ingest_gtop has parsed

    "GtoPdb Version: 2026.3.1 - published: 2026-10-12"

out of the first line of the file since the quoted-comment fix, printed it
to a terminal, and dropped it.
"""
from __future__ import annotations

import json

import pytest

from kgav.gtop import parse_version
from kgav.provenance import SOURCE_STAMP_NAME, read_source, write_source


# --------------------------------------------------------------- the parser
def test_the_real_declaration_splits_into_release_and_date():
    v = parse_version("GtoPdb Version: 2026.3.1 - published: 2026-10-12")
    assert v == {"raw": "GtoPdb Version: 2026.3.1 - published: 2026-10-12",
                 "release": "2026.3.1", "published": "2026-10-12"}


def test_a_published_date_is_the_source_s_own_claim_not_a_file_mtime():
    """Every other "release date" in this repo is a download timestamp --
    string_release_date is literally a file mtime. This one the source
    states about itself, which is why it gets its own field."""
    assert parse_version(
        "GtoPdb Version: 2026.3.1 - published: 2026-10-12"
    )["published"] == "2026-10-12"


@pytest.mark.parametrize("decl,release,published", [
    ("GtoPdb Version: 2027.1", "2027.1", None),
    ("GtoPdb Version: 2026.3.1 - published: 2026-10-12", "2026.3.1",
     "2026-10-12"),
    ("gtopdb version: 2026.4 - PUBLISHED: 2027-01-05", "2026.4", "2027-01-05"),
    ("some future header with no version at all", None, None),
    ("", None, None),
    (None, None, None),
])
def test_nothing_is_guessed(decl, release, published):
    """A release with no date gives None for the date. A provenance record
    with a plausible invented value in it is worse than one with a gap,
    because a later reader cannot tell which it is looking at."""
    v = parse_version(decl)
    assert (v["release"], v["published"]) == (release, published)


# ---------------------------------------------------------------- the stamp
def test_the_stamp_round_trips(tmp_path):
    write_source(tmp_path, guidetopharmacology={"release": "2026.3.1"})
    assert read_source(tmp_path) == {
        "guidetopharmacology": {"release": "2026.3.1"}}


def test_a_source_with_no_version_is_omitted_not_stamped_unrecorded(tmp_path):
    """`"version": "unrecorded"` reads later as a value. The key is left
    out, and the MANIFEST's absence list is what reports it."""
    write_source(tmp_path, guidetopharmacology=None, orcs={"queried": "x"})
    body = json.loads((tmp_path / SOURCE_STAMP_NAME).read_text())
    assert body["sources"] == {"orcs": {"queried": "x"}}
    assert "unrecorded" not in (tmp_path / SOURCE_STAMP_NAME).read_text()


def test_an_unstamped_layer_reads_as_none_not_empty(tmp_path):
    """None is "this layer does not say"; {} is "it ran the stamp and had
    nothing to record". Every layer but gtop is the first today, and the
    two want different work."""
    assert read_source(tmp_path) is None
    write_source(tmp_path)
    assert read_source(tmp_path) == {}


def test_an_unreadable_stamp_is_not_a_crash(tmp_path):
    (tmp_path / SOURCE_STAMP_NAME).write_text("{not json")
    assert read_source(tmp_path) is None


def test_the_stamp_is_not_written_for_a_layer_that_does_not_exist(tmp_path):
    """A refused ingest writes no layer directory. Stamping a version onto
    nothing would leave a provenance record for a build that did not
    happen."""
    gone = tmp_path / "never-built"
    write_source(gone, guidetopharmacology={"release": "2026.3.1"})
    assert not gone.exists()


# ------------------------------------------------- the ingest refuses a mix
# NOT PART OF V1 AS WRITTEN, found while making the version a checked fact.
# Only interactions.csv was ever read for a version, so two downloads from
# different releases were undetectable -- and GtoPdb reuses ligand ids
# across releases, so the join between the two files SUCCEEDS and produces
# a graph that is neither release. A silent wrong answer, not an error.
def _gtop_raw(tmp_path, inter_version, ligand_version):
    from tests.test_gtop import IHEAD, LHEAD, _inter, _row_text
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "interactions.csv").write_text(
        f'"# {inter_version}"\n{IHEAD}\n{_row_text(IHEAD, **_inter())}\n')
    (raw / "ligands.csv").write_text(f'"# {ligand_version}"\n{LHEAD}\n')
    return raw


def _run_ingest(raw, out, release):
    import subprocess
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    return subprocess.run(
        [sys.executable, str(root / "scripts" / "ingest_gtop.py"),
         "--raw", str(raw), "--out", str(out), "--release", str(release)],
        capture_output=True, text=True, cwd=root)


def test_two_downloads_from_different_releases_are_refused(tmp_path):
    raw = _gtop_raw(tmp_path, "GtoPdb Version: 2026.3.1",
                    "GtoPdb Version: 2025.4")
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text("")
    r = _run_ingest(raw, tmp_path / "out", rel)
    assert r.returncode == 2, r.stdout + r.stderr
    assert "different GtoPdb releases" in r.stdout
    assert "2026.3.1" in r.stdout and "2025.4" in r.stdout
    assert not (tmp_path / "out").exists(), \
        "a refused ingest must not leave a layer directory"


# --------------------------------------------- and into the MANIFEST
def _layer(tmp_path, name, sources=None):
    d = tmp_path / f"v0.1-{name}"
    d.mkdir()
    (d / "edges.jsonl").write_text("")
    if sources is not None:
        write_source(d, **sources)
    return d


def test_the_manifest_block_carries_the_declared_version(tmp_path):
    import sys
    sys.path.insert(0, str(__import__("pathlib").Path(__file__)
                          .resolve().parents[1]))
    from scripts.assemble_release import source_block
    layers = {"gtop": _layer(tmp_path, "gtop", {
        "guidetopharmacology": {"release": "2026.3.1",
                                "published": "2026-10-12"}})}
    block = source_block(layers)
    assert block["declared"]["gtop"]["guidetopharmacology"]["release"] == \
        "2026.3.1"
    assert block["layers_declaring_no_source_version"] == []


def test_the_layers_that_declare_nothing_are_named(tmp_path):
    """The absences are the finding. A `sources` block holding only gtop
    reads as the whole answer, and today gtop is the ONLY layer that can
    say which release it read -- v0.1 names no ORCS, NCATS, InterPro,
    STRING, UniProt or RDKit version at all."""
    from scripts.assemble_release import source_block
    layers = {
        "gtop": _layer(tmp_path, "gtop", {"guidetopharmacology": {"release": "x"}}),
        "orcs": _layer(tmp_path, "orcs"),
        "chembl": _layer(tmp_path, "chembl"),
        "ncats": _layer(tmp_path, "ncats"),
    }
    block = source_block(layers)
    assert set(block["declared"]) == {"gtop"}
    assert block["layers_declaring_no_source_version"] == [
        "chembl", "ncats", "orcs"]


def test_an_empty_declaration_counts_as_no_declaration(tmp_path):
    """A layer that ran the stamp with nothing to record has not answered
    the question, so it belongs in the absence list rather than in
    `declared` with an empty dict beside the layers that did answer."""
    from scripts.assemble_release import source_block
    layers = {"gtop": _layer(tmp_path, "gtop", {})}
    block = source_block(layers)
    assert block["declared"] == {}
    assert block["layers_declaring_no_source_version"] == ["gtop"]
