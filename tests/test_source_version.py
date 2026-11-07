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
from pathlib import Path

import pytest

from kgav.gtop import parse_version
from kgav.provenance import SOURCE_STAMP_NAME, read_source, write_source

ROOT = Path(__file__).resolve().parents[1]


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


# ------------------------------------------ the version is read, not declared
# Fix-list V2/V3/V5. Three sources state their release somewhere a reader
# can find it; the rest genuinely do not, and are left out rather than given
# a plausible value.
#
#   ChEMBL   "chembl_37" from the database PATH. This build has no version
#            table -- SELECT name, label FROM version raises
#            OperationalError on the chembl_37 SQLite distribution -- so
#            the directory name is the only statement that travels with
#            the file.
#   STRING   "v12.0" from the download's FILENAME.
#   RDKit    rdkit.__version__, the library actually imported.
#
#   UniProt  carries the release in a REST response header that
#            fetch_raw.py does not capture. Its 17 proteome sidecars hold
#            a retrieval date and a URL and no version. A retrieval date
#            is not a release, and writing one into a version field is
#            exactly fix-list V4 -- string_release_date was a download
#            mtime called a release date.
#   ORCS     no sidecars at all; screens_all.json is a bare list of 1,952.
#   NCATS    three CSV/TSV files, no sidecars.
#   InterPro the folds cache has no date field.
@pytest.mark.parametrize("path,want", [
    ("data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db", "chembl_37"),
    ("/elsewhere/chembl_35/x.db", "chembl_35"),
    ("/no/version/here.db", None),
])
def test_chembl_release_is_read_from_the_path(path, want):
    from kgav.provenance import chembl_release
    assert chembl_release(Path(path)) == want


@pytest.mark.parametrize("name,want", [
    ("9606.protein.links.detailed.v12.0.txt.gz", "v12.0"),
    ("9606.protein.aliases.v11.5.txt.gz", "v11.5"),
    ("links.txt.gz", None),
])
def test_string_release_is_read_from_the_filename(name, want):
    from kgav.provenance import string_release
    assert string_release(name) == want


def test_rdkit_release_is_the_library_actually_imported():
    from kgav.provenance import rdkit_release
    v = rdkit_release()
    if v is None:
        pytest.skip("rdkit is an optional extra and is not installed")
    import rdkit
    assert v == rdkit.__version__


def test_the_mtime_field_is_no_longer_called_a_release_date():
    """V4. ingest_host's own comment said string_release_date came "from
    the download timestamp", and that value dates all 141,019 STRING
    edges. The VALUE is unchanged on purpose -- the alternatives are
    STRING's real release date, which nothing on disk states, or the
    undated placeholder, which would move every one of those edges across
    the 2021 temporal cut and is the whole of fix-list C2."""
    host = (ROOT / "src" / "kgav" / "host.py").read_text()
    driver = (ROOT / "scripts" / "ingest_host.py").read_text()
    assert 'prov["string_file_mtime"]' in host
    assert '"string_release_date":' not in driver, (
        "the PROV key still claims to be a release date")
    assert '"string_file_mtime": "2023-08-28"' in driver


def test_every_layer_that_can_declare_a_version_does():
    """Checked at the source: these five drivers read a file whose release
    is discoverable, so each must stamp it. A layer that quietly stops
    calling write_source would otherwise reappear in the MANIFEST's
    absence list with nothing to explain why."""
    for script, token in (("ingest_chembl.py", "chembl_release"),
                          ("ingest_host_targets.py", "chembl_release"),
                          ("ingest_selectivity.py", "chembl_release"),
                          ("ingest_similarity.py", "rdkit_release"),
                          ("ingest_host.py", "string_release")):
        src = (ROOT / "scripts" / script).read_text()
        assert "write_source(" in src, f"{script} stamps no source version"
        assert token in src, f"{script} does not read its source's release"


def test_a_source_with_nothing_to_declare_is_not_given_a_date(tmp_path):
    """The V4 error, generalised: a retrieval date or a file mtime must
    never be written where a reader expects a release."""
    for script in ("ingest_orcs.py", "ingest_ncats.py", "ingest_folds.py"):
        src = (ROOT / "scripts" / script).read_text()
        assert "write_source(" not in src, (
            f"{script} stamps a source version, but none of ORCS, NCATS or "
            f"InterPro states a release anywhere on disk -- check it is not "
            f"recording an mtime or a retrieval date as a version")
