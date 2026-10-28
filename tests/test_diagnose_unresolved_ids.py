"""Gate for narrowing the unresolved-id diagnostic to a few viruses.

ingest_vhppi names the viruses whose rows all died before becoming an
edge; this script is the one you run next. If it looked at a different
set of files than the report counted, it would send the reader to the
wrong rows while appearing to have answered -- the failure mode this
whole layer keeps reproducing.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from diagnose_unresolved_ids import select_files

DOC = {"defaults": {"family": "Coronaviridae"}, "viruses": [
    {"name": "DENV-1", "taxon": "11053", "family": "Flaviviridae"},
    {"name": "Yellow fever", "taxon": "11089", "family": "Flaviviridae"},
    # NO per-entry family, like the real register: inherited from
    # defaults. The fixture declared it, so --family Coronaviridae was
    # never exercised against the path that was broken.
    {"name": "SARS-CoV-2", "taxon": "2697049"},
]}
MAP = {"11053": "11053", "11089": "11089", "2697049": "2697049"}


def _files(tmp_path, *stems):
    out = []
    for s in stems:
        f = tmp_path / f"vhn_{s}.tab27"
        f.write_text("")
        out.append(f)
    return out


def test_a_family_keeps_only_that_familys_files(tmp_path):
    files = _files(tmp_path, "denv1", "yellowfever", "sarscov2")
    keep, dropped = select_files(files, DOC, MAP, "Flaviviridae", None)
    assert [f.name for f in keep] == ["vhn_denv1.tab27",
                                      "vhn_yellowfever.tab27"]
    assert dropped == 1


def test_a_named_virus_keeps_only_its_file(tmp_path):
    files = _files(tmp_path, "denv1", "yellowfever", "sarscov2")
    keep, _ = select_files(files, DOC, MAP, None, ["DENV-1"])
    assert [f.name for f in keep] == ["vhn_denv1.tab27"]


def test_virus_is_repeatable_and_unions_with_family(tmp_path):
    files = _files(tmp_path, "denv1", "yellowfever", "sarscov2")
    keep, _ = select_files(files, DOC, MAP, None, ["DENV-1", "SARS-CoV-2"])
    assert len(keep) == 2


def test_the_name_is_matched_case_insensitively(tmp_path):
    files = _files(tmp_path, "denv1")
    assert select_files(files, DOC, MAP, None, ["denv-1"])[0]


def test_an_unknown_virus_is_an_error_not_an_empty_result(tmp_path):
    """"--virus DENV1" quietly returning nothing reads as "that virus has
    no unresolved rows", which is the opposite of the truth."""
    files = _files(tmp_path, "denv1")
    with pytest.raises(SystemExit) as e:
        select_files(files, DOC, MAP, None, ["DENV1"])
    assert "not a register virus" in str(e.value)
    assert "DENV-1" in str(e.value)


def test_the_sidecar_wins_over_the_filename(tmp_path):
    """Attribution is vhppi.file_taxon, the same rule the report used."""
    f = tmp_path / "vhn_sarscov2.tab27"
    f.write_text("")
    (tmp_path / "vhn_sarscov2.tab27.source.json").write_text(
        '{"taxon": "11053"}')
    keep, _ = select_files([f], DOC, MAP, "Flaviviridae", None)
    assert keep == [f], "the sidecar says DENV-1 whatever the name says"


def test_an_unattributable_file_is_left_out_rather_than_kept(tmp_path):
    """An intact_*.tab27 covers many viruses. Keeping it under --virus
    DENV-1 would mix another virus's rows into the answer."""
    f = tmp_path / "intact_viral.tab27"
    f.write_text("")
    keep, dropped = select_files([f], DOC, MAP, None, ["DENV-1"])
    assert keep == [] and dropped == 1


def test_an_isolate_taxon_in_a_sidecar_still_matches_its_species(tmp_path):
    f = tmp_path / "x.tab27"
    f.write_text("")
    (tmp_path / "x.tab27.source.json").write_text('{"taxon": "31634"}')
    keep, _ = select_files([f], DOC, {**MAP, "31634": "11053"},
                           "Flaviviridae", None)
    assert keep == [f]


def test_a_register_virus_declared_under_a_strain_still_matches_its_file():
    """Both SIDES of the comparison need folding, not just the file's.

    Today every register `taxon` is a species, so taxon_map maps it to
    itself and folding the wanted set is a no-op -- which is exactly why
    a mutant removing it survived the other tests. The register's own
    history makes the other shape real: the flaviviruses were surveyed
    under pinned STRAIN taxa, and a virus declared that way would be
    asked for under the strain while its file resolves to the species.
    Matching those would then silently return no files, and an empty
    result here reads as "nothing unresolved for this virus".
    """
    import tempfile
    doc = {"viruses": [{"name": "DENV-2", "taxon": "31634",
                        "family": "Flaviviridae"}]}
    tmap = {"31634": "11060", "11060": "11060"}
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "vhn_denv2.tab27"
        f.write_text("")
        (Path(d) / "vhn_denv2.tab27.source.json").write_text(
            '{"taxon": "11060"}')
        keep, dropped = select_files([f], doc, tmap, "Flaviviridae", None)
    assert keep == [f] and dropped == 0


def test_family_coronaviridae_selects_the_inheriting_entries(tmp_path):
    """They carry no per-entry family. A filter reading the key alone
    matched nothing, and an empty selection here reads as "that virus
    has no unresolved rows"."""
    files = _files(tmp_path, "denv1", "sarscov2")
    keep, dropped = select_files(files, DOC, MAP, "Coronaviridae", None)
    assert [f.name for f in keep] == ["vhn_sarscov2.tab27"]
    assert dropped == 1
