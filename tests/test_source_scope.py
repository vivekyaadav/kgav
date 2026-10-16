"""What kind of literature ChEMBL is, measured rather than asserted.

The claim this script exists to support is strong -- that the arm structure
phase 4 measured is partly a property of the database's curation scope --
so the ways it can be WRONG are the ways that claim becomes an impression
dressed as a measurement: a journal classified in code, a window chosen
after seeing the share it produces, or a year table truncated to the rows
that suit the conclusion.
"""
from __future__ import annotations

import sqlite3
import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _mod(name):
    sys.path.insert(0, str(ROOT / "scripts"))
    return import_module(name)


D = _mod("diagnose_source_scope")
P = _mod("phase4_survey_chemistry")


def _fixture(tmp_path: Path) -> Path:
    """Reuse the phase 4 fixture: one database, two reports, no second
    definition of what a viral protease target is."""
    sys.path.insert(0, str(ROOT / "tests"))
    p4 = import_module("test_phase4_survey")
    return p4._db(tmp_path / "scope.db", arm_a=12, arm_b=12)


# --------------------------------------------------- one reader, two reports
def test_the_target_definitions_are_the_survey_s_own():
    """Not a copy. A second definition of which targets are viral proteases
    would be a second answer to the question of whether the SOURCE and the
    SURVEY agree -- and they cannot disagree if there is one of them."""
    assert D.load_targets is P.load_targets
    assert D.row_role is P.row_role
    assert D.role_key is P.role_key
    assert D.FAMILIES is P.FAMILIES


# -------------------------------------------------------- the journal profile
def test_the_profile_counts_documents_not_rows(tmp_path):
    con = sqlite3.connect(_fixture(tmp_path))
    prof, total = D.journal_profile(con)
    con.close()
    assert prof == {"J Med Chem": 1, "PLoS One": 1}
    assert total == 2


def test_no_journal_is_classified_in_code():
    """"Medicinal chemistry" is a judgement, and a judgement asserted in
    code reads later as a measurement. The script prints names and counts;
    the reader classifies them."""
    src = (ROOT / "scripts" / "diagnose_source_scope.py").read_text()
    body = src.split('"""', 2)[2]          # past the module docstring
    for word in ("med_chem", "MED_CHEM", "is_medchem", "JOURNAL_KIND"):
        assert word not in body


def test_the_venue_lookup_is_case_insensitive(tmp_path):
    """ChEMBL spells the same journal several ways across releases. A
    lookup that misses "PLOS ONE" would report a zero that is a spelling."""
    p = _fixture(tmp_path)
    con = sqlite3.connect(p)
    con.execute("UPDATE docs SET journal = 'PLOS ONE' WHERE doc_id = 102")
    con.commit()
    found = {v: n for v, n, _ in D.venue_counts(con)}
    con.close()
    assert found["PLoS One"] == 1


def test_the_venue_list_is_a_lookup_not_a_verdict():
    """Every name is a title to count, and the list says so rather than
    claiming to be exhaustive."""
    assert "PLoS One" in D.REPURPOSING_VENUES
    assert len(set(D.REPURPOSING_VENUES)) == len(D.REPURPOSING_VENUES)


def test_the_window_is_declared_before_the_share_is_seen():
    """A boundary chosen after seeing the number it produces is not a
    boundary. It is a constant with a comment saying why."""
    assert D.REPURPOSING_WINDOW == (2020, 2021)
    src = (ROOT / "scripts" / "diagnose_source_scope.py").read_text()
    assert "cannot be read off a boundary chosen" in src


# ------------------------------------------------------------- the year table
def test_every_year_is_printed(tmp_path, capsys, monkeypatch):
    """THE REGRESSION THIS SCRIPT EXISTS TO AVOID. A top-fifteen per-year
    table showed 2022-2025 and no 2020 or 2021, which would say the
    coronavirus protease data post-dates the repurposing wave entirely.
    That may be true and it may be an artifact of taking fifteen rows."""
    p = _fixture(tmp_path)
    con = sqlite3.connect(p)
    # spread the 3CLpro rows over many years, one row each for most
    for i, yr in enumerate(range(2005, 2026)):
        con.execute("INSERT INTO docs VALUES (?,?,?,?,?,?)",
                    (300 + i, yr, "J Med Chem", "t", None, None))
        con.execute("INSERT INTO assays VALUES (?,?,?,?,?)",
                    (400 + i, 1, "Inhibition of recombinant SARS-CoV-2 3CL "
                                 "protease by FRET assay", "B", 300 + i))
        con.execute("INSERT INTO activities VALUES (?,?,?,?,?,?,?)",
                    (8000 + i, 400 + i, 1, "IC50", 100.0, "nM", "="))
    con.commit()
    con.close()
    monkeypatch.setattr(sys, "argv",
                        ["diagnose_source_scope.py", "--db", str(p)])
    assert D.main() == 0
    out = capsys.readouterr().out
    assert "rows by year, complete:" in out
    for yr in range(2005, 2026):
        assert f"      {yr}  " in out, f"{yr} missing from the year table"


def test_undated_rows_are_counted_not_dropped(tmp_path):
    """A row whose document has no year is not a row that did not happen."""
    p = _fixture(tmp_path)
    con = sqlite3.connect(p)
    con.execute("INSERT INTO docs VALUES (999, NULL, 'J Med Chem', 't', "
                "NULL, NULL)")
    con.execute("INSERT INTO assays VALUES (995, 1, 'Inhibition of "
                "recombinant SARS-CoV-2 3CL protease by FRET assay', 'B', 999)")
    con.execute("INSERT INTO activities VALUES (8999, 995, 1, 'IC50', 50.0, "
                "'nM', '=')")
    con.commit()
    targets = D.load_targets(con)
    acts = D.load_activities(con, [t["tid"] for t in targets])
    docs = {d: (y, j) for d, y, j in
            con.execute("SELECT doc_id, year, journal FROM docs")}
    con.close()
    rows = D.protease_rows(acts, {t["tid"]: t for t in targets})
    _years, journals, undated = D.by_year_and_journal(rows["3clpro"], docs)
    assert undated >= 1
    # and it still contributes its journal, because that is known
    assert journals["J Med Chem"] >= 1


def test_a_row_with_no_document_at_all_is_undated(tmp_path):
    con = sqlite3.connect(_fixture(tmp_path))
    targets = D.load_targets(con)
    acts = D.load_activities(con, [t["tid"] for t in targets])
    con.close()
    rows = D.protease_rows(acts, {t["tid"]: t for t in targets})
    _y, _j, undated = D.by_year_and_journal(rows["3clpro"], {})
    assert undated == len(rows["3clpro"])


# ------------------------------------------------------------------- the CLI
def test_a_missing_database_exits_two(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["diagnose_source_scope.py", "--db",
                                      str(tmp_path / "absent.db")])
    assert D.main() == 2
    assert "no ChEMBL database at" in capsys.readouterr().out


def test_it_writes_nothing(tmp_path, capsys, monkeypatch):
    p = _fixture(tmp_path)
    before = {f.name: f.stat().st_mtime_ns for f in tmp_path.iterdir()}
    monkeypatch.setattr(sys, "argv",
                        ["diagnose_source_scope.py", "--db", str(p)])
    assert D.main() == 0
    after = {f.name: f.stat().st_mtime_ns for f in tmp_path.iterdir()}
    assert before == after
    assert "recommends nothing" in capsys.readouterr().out


@pytest.mark.parametrize("flag", ["--top"])
def test_the_top_flag_caps_only_the_journal_profile(tmp_path, capsys,
                                                    monkeypatch, flag):
    p = _fixture(tmp_path)
    monkeypatch.setattr(sys, "argv",
                        ["diagnose_source_scope.py", "--db", str(p), flag, "1"])
    assert D.main() == 0
    out = capsys.readouterr().out
    assert "the other 1 titles supply" in out
    assert "rows by year, complete:" in out
