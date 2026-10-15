"""Phase 4's feasibility survey: the classifiers, the rejections, the verdict.

The survey decides whether the HCV retrospective can be measured at all, so
the ways it can be WRONG are the ways phase 4 gets built on sand: a target
claimed for the wrong enzyme, a replicon number counted as an enzyme
measurement, a polyprotein row counted as a protease row, or a rejection
reported as zero because the count was taken over the wrong set.
"""
from __future__ import annotations

import sqlite3
import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _survey():
    sys.path.insert(0, str(ROOT / "scripts"))
    return import_module("phase4_survey_chemistry")


S = _survey()


# ------------------------------------------------------------ the families
def test_the_three_compared_families_are_recognised():
    assert S.classify_family("Severe acute respiratory syndrome "
                             "coronavirus 2") == "Coronaviridae"
    assert S.classify_family("Human rhinovirus 14") == "Picornaviridae"
    assert S.classify_family("Hepatitis C virus") == "Flaviviridae"
    assert S.classify_family("Zika virus") == "Flaviviridae"


def test_an_unmatched_organism_is_out_of_scope_rather_than_guessed():
    assert S.classify_family("Homo sapiens") is None
    assert S.classify_family("Influenza A virus") is None
    assert S.classify_family(None) is None


def test_hepatitis_a_is_a_picornavirus_and_hepatitis_c_is_not():
    """Both start "Hepatitis", and the picornavirus pattern is tried first.
    Getting this backwards puts HAV 3C -- a cysteine protease, one of the
    folds layer's five picornavirus comparators -- into the serine arm."""
    assert S.classify_family("Hepatitis A virus") == "Picornaviridae"
    assert S.classify_family("Hepatitis C virus") == "Flaviviridae"
    assert S.classify_family("Hepacivirus C") == "Flaviviridae"


def test_hcv_is_split_out_of_the_flaviviruses():
    assert S.is_hcv("Hepatitis C virus") and S.is_hcv("Hepacivirus C")
    assert not S.is_hcv("Zika virus")
    assert not S.is_hcv("Dengue virus 2")


# --------------------------------------------------------------- the roles
@pytest.mark.parametrize("family,name,role", [
    ("Coronaviridae", "3C-like proteinase", "3clpro"),
    ("Coronaviridae", "3CL-protease", "3clpro"),
    ("Coronaviridae", "SARS-CoV-2 main protease", "3clpro"),
    ("Coronaviridae", "Papain-like protease nsp3", "plpro"),
    ("Picornaviridae", "Protease 3C", "3cpro"),
    ("Picornaviridae", "Picornain 3C", "3cpro"),
    ("Flaviviridae", "NS2B-NS3 protease", "ns3pro"),
    ("Flaviviridae", "NS3/4A protease", "ns3pro"),
    ("Flaviviridae", "Serine protease NS3", "ns3pro"),
])
def test_each_role_is_claimed_by_its_own_names(family, name, role):
    assert S.classify_role(family, name) == role


def test_3c_does_not_claim_3c_like():
    """"3C" is a substring of "3C-like". ingest_folds calls the substring
    hazard the sixth instance in this project; this is the seventh place it
    could bite. A coronavirus 3CLpro landing in the picornavirus arm would
    make arm A compare the coronaviruses with themselves."""
    assert S.classify_role("Coronaviridae", "3C-like proteinase") == "3clpro"
    # and the picornavirus table never sees a coronavirus name, because the
    # family is decided by the organism before the role is read at all
    assert S.classify_role("Picornaviridae", "3C-like proteinase") == "3cpro"


def test_3cd_is_not_the_3c_protease():
    """3CD is the uncleaved precursor, a different enzyme with different
    specificity. It is rejected by name."""
    assert S.classify_role("Picornaviridae", "Protease 3CD") is None


def test_ns3_helicase_is_the_wrong_enzyme():
    """Flavivirus NS3 is protease AND helicase in one chain -- the split
    ingest_folds makes by domain coordinates. A target named for the
    helicase alone must not supply a protease measurement."""
    assert S.classify_role("Flaviviridae", "NS3 helicase") is None
    assert S.classify_role("Flaviviridae", "NS3 NTPase") is None


def test_the_bifunctional_construct_is_still_the_protease_target():
    """REGRESSION. Written as `helicase(?!.*protease)` the rejection fired on
    "NS3 protease/helicase", because there the word protease comes FIRST and
    the lookahead only sees what follows. That dropped the main HCV NS3
    target -- the one arm B is entirely about -- and the survey would have
    reported arm B as unmeasurable for a reason that was not true."""
    assert S.classify_role("Flaviviridae",
                           "NS3 protease/helicase NS3") == "ns3pro"
    assert S.classify_role("Flaviviridae",
                           "NS3 helicase/protease") == "ns3pro"


def test_a_polyprotein_target_has_no_role():
    """A measurement filed against a polyprotein does not say which of its
    proteins was inhibited. v1 chapter 9, the first phase 0 script and the
    flavivirus NS3 rows are three previous instances."""
    for family in ("Coronaviridae", "Picornaviridae", "Flaviviridae"):
        assert S.classify_role(family, "Genome polyprotein") is None
    assert S.classify_role("Coronaviridae", "Replicase polyprotein 1ab") is None


def test_role_key_separates_hcv_from_the_other_flaviviruses():
    hcv = {"family": "Flaviviridae", "role": "ns3pro", "hcv": True}
    zika = {"family": "Flaviviridae", "role": "ns3pro", "hcv": False}
    cov = {"family": "Coronaviridae", "role": "3clpro", "hcv": False}
    assert S.role_key(hcv) == "hcv_ns3pro"
    assert S.role_key(zika) == "flavi_ns3pro"
    assert S.role_key(cov) == "3clpro"


# ------------------------------------------------------------ the database
def _db(path: Path, *, arm_a: int = 12, arm_b: int = 12,
        picornavirus: bool = True) -> Path:
    """A ChEMBL-shaped file with a known answer, so the SQL really runs."""
    con = sqlite3.connect(path)
    con.executescript("""
    CREATE TABLE target_dictionary (tid INTEGER PRIMARY KEY, chembl_id TEXT,
        pref_name TEXT, organism TEXT, tax_id INTEGER, target_type TEXT);
    CREATE TABLE assays (assay_id INTEGER PRIMARY KEY, tid INTEGER,
        description TEXT, assay_type TEXT, doc_id INTEGER);
    CREATE TABLE activities (activity_id INTEGER PRIMARY KEY, assay_id INTEGER,
        molregno INTEGER, standard_type TEXT, standard_value REAL,
        standard_units TEXT, standard_relation TEXT);
    CREATE TABLE molecule_dictionary (molregno INTEGER PRIMARY KEY,
        chembl_id TEXT, pref_name TEXT, max_phase REAL);
    """)
    targets = [
        (1, "CHEMBL_COV", "3C-like proteinase",
         "Severe acute respiratory syndrome coronavirus 2", 2697049,
         "SINGLE PROTEIN"),
        (3, "CHEMBL_HCV", "NS3 protease/helicase NS3", "Hepatitis C virus",
         11103, "SINGLE PROTEIN"),
        # ORGANISM-level: replicon and cell-culture numbers, no role
        (4, "CHEMBL_HCVORG", "Hepatitis C virus", "Hepatitis C virus", 11103,
         "ORGANISM"),
        (5, "CHEMBL_POLY", "Genome polyprotein", "Zika virus", 64320,
         "SINGLE PROTEIN"),
    ]
    if picornavirus:
        targets.append((2, "CHEMBL_RV", "Protease 3C", "Human rhinovirus 14",
                        12131, "SINGLE PROTEIN"))
    con.executemany("INSERT INTO target_dictionary VALUES (?,?,?,?,?,?)",
                    targets)
    con.executemany("INSERT INTO assays VALUES (?,?,?,?,?)", [
        (10, 1, ("Inhibition of recombinant SARS-CoV-2 3CL protease by FRET "
                 "assay"), "B", 1),
        (20, 2, "Inhibition of purified human rhinovirus 3C protease", "B", 1),
        (30, 3, "Inhibition of recombinant HCV NS3 protease, enzymatic", "B", 1),
        (40, 4, ("Antiviral activity against HCV in Huh-7 cells by replicon "
                 "assay"), "F", 1),
        (50, 1, "Antiviral activity against SARS-CoV-2 in Vero cells", "F", 1),
    ])
    mols = [(i, f"CHEMBL{i:06d}", f"drug{i}", 0.0) for i in range(1, 61)]
    acts, aid = [], 0

    def add(mol, assay, units="nM"):
        nonlocal aid
        aid += 1
        acts.append((aid, assay, mol, "IC50", 100.0, units, "="))

    for i in range(1, arm_a + 1):
        add(i, 10)
        add(i, 20)
    for i in range(21, 21 + arm_b):
        add(i, 10)
        add(i, 30)
    # HCV enzyme + coronavirus CELLULAR only: in the broad set, not in arm B
    for i in range(41, 46):
        add(i, 30)
        add(i, 50)
    # replicon-only, on the ORGANISM target
    add(60, 40)
    con.executemany("INSERT INTO molecule_dictionary VALUES (?,?,?,?)", mols)
    con.executemany("INSERT INTO activities VALUES (?,?,?,?,?,?,?)", acts)
    con.commit()
    con.close()
    return path


def test_missing_columns_names_every_one_it_needs(tmp_path):
    p = tmp_path / "stunted.db"
    con = sqlite3.connect(p)
    con.executescript("CREATE TABLE target_dictionary (tid INTEGER, "
                      "chembl_id TEXT);")
    con.commit()
    gone = S.missing_columns(con)
    con.close()
    assert "target_dictionary.organism" in gone
    assert "activities (no such table)" in gone
    # the columns it DOES have are not reported missing
    assert "target_dictionary.tid" not in gone


def test_a_real_chembl_shaped_db_passes_the_column_check(tmp_path):
    con = sqlite3.connect(_db(tmp_path / "ok.db"))
    assert S.missing_columns(con) == []
    con.close()


def test_load_activities_refuses_an_empty_target_list(tmp_path):
    """`IN ()` is always-false in SQLite, so a caller bug would read as
    these proteases having no chemistry. chembl.query_activities documents
    the same failure costing a 8,748-row virus a reading of zero."""
    con = sqlite3.connect(_db(tmp_path / "ok.db"))
    with pytest.raises(ValueError, match="no targets"):
        S.load_activities(con, [])
    con.close()


def _loaded(tmp_path, **kw):
    con = sqlite3.connect(_db(tmp_path / "ok.db", **kw))
    targets = S.load_targets(con)
    acts = S.load_activities(con, [t["tid"] for t in targets])
    con.close()
    return targets, {t["tid"]: t for t in targets}, acts


def test_only_matched_organisms_are_loaded(tmp_path):
    targets, _, _ = _loaded(tmp_path)
    assert {t["chembl_id"] for t in targets} == {
        "CHEMBL_COV", "CHEMBL_RV", "CHEMBL_HCV", "CHEMBL_HCVORG",
        "CHEMBL_POLY"}


def test_a_cellular_row_is_not_an_enzyme_measurement(tmp_path):
    """Five compounds have an HCV enzyme number and a SARS-CoV-2 Vero
    number. Counting the Vero number as a 3CLpro measurement would put them
    in arm B, which is a claim about the enzyme."""
    _, by_tid, acts = _loaded(tmp_path)
    protein = S.compounds_by_role(acts, by_tid, protein_only=True)
    anywhere = S.compounds_by_role(acts, by_tid, protein_only=False)
    assert len(protein["3clpro"]) == 24          # arm A 12 + arm B 12
    assert len(anywhere["3clpro"]) == 29         # + the 5 cellular
    assert len(protein["hcv_ns3pro"]) == 17      # arm B 12 + the 5


def test_an_organism_level_target_reaches_no_role(tmp_path):
    """The HCV replicon compound must not appear as an NS3 inhibitor."""
    _, by_tid, acts = _loaded(tmp_path)
    protein = S.compounds_by_role(acts, by_tid, protein_only=False)
    assert "CHEMBL000060" not in protein.get("hcv_ns3pro", set())


def test_the_broad_set_keeps_the_cellular_coronavirus_rows(tmp_path):
    """Arm B needs both ends on the enzyme; the broad set is the weaker
    question the 2020 literature actually asked, and it must still see the
    five compounds whose only coronavirus number is a Vero EC50."""
    _, by_tid, acts = _loaded(tmp_path)
    protein = S.compounds_by_role(acts, by_tid, protein_only=True)
    broad = S.broad_coronavirus_compounds(acts, by_tid)
    assert len(protein["hcv_ns3pro"] & broad) == 17


def test_an_unconvertible_unit_is_dropped_rather_than_guessed(tmp_path):
    """ug.mL-1 needs a molecular weight; chembl.to_nm returns None rather
    than inventing one, and a compound measured only that way is not in an
    arm."""
    p = _db(tmp_path / "ok.db")
    con = sqlite3.connect(p)
    con.execute("INSERT INTO activities VALUES (9001, 10, 59, 'IC50', 1.0, "
                "'ug.mL-1', '=')")
    con.execute("INSERT INTO activities VALUES (9002, 20, 59, 'IC50', 1.0, "
                "'ug.mL-1', '=')")
    con.commit()
    targets = S.load_targets(con)
    acts = S.load_activities(con, [t["tid"] for t in targets])
    con.close()
    protein = S.compounds_by_role(acts, {t["tid"]: t for t in targets})
    assert "CHEMBL000059" not in protein["3clpro"]


# ------------------------------------------------------------- the verdict
def _main(monkeypatch, path):
    monkeypatch.setattr(sys, "argv",
                        ["phase4_survey_chemistry.py", "--db", str(path)])
    return S.main()


def test_both_arms_full_is_measurable(tmp_path, monkeypatch, capsys):
    code = _main(monkeypatch, _db(tmp_path / "ok.db", arm_a=12, arm_b=12))
    out = capsys.readouterr().out
    assert code == 0
    assert "MEASURABLE." in out and "NOT MEASURABLE" not in out


def test_thin_arms_refuse_rather_than_report_a_refutation(tmp_path,
                                                          monkeypatch, capsys):
    """An unmeasurable test is not a failed test, and exiting 0 here would
    let phase 4 read "no transfer found" off a sample of four."""
    code = _main(monkeypatch, _db(tmp_path / "thin.db", arm_a=4, arm_b=2))
    out = capsys.readouterr().out
    assert code == 1
    assert "NOT MEASURABLE" in out
    assert "not mistaken for a refutation" in out


def test_a_family_with_no_targets_is_an_input_problem_not_a_finding(
        tmp_path, monkeypatch, capsys):
    code = _main(monkeypatch,
                 _db(tmp_path / "nopic.db", picornavirus=False))
    out = capsys.readouterr().out
    assert code == 2
    assert "No target at all for: Picornaviridae" in out


def test_a_missing_database_exits_two(tmp_path, monkeypatch, capsys):
    code = _main(monkeypatch, tmp_path / "absent.db")
    assert code == 2
    assert "no ChEMBL database at" in capsys.readouterr().out


def test_the_floor_is_the_one_this_project_already_uses():
    """Not a third 10. calibration's own comment says 10 is already
    hard_negatives.MIN_CLASS_SIZE, and a comment keeping two numbers in step
    is the pattern this session replaced nine times."""
    from kgav.calibration import MIN_REACHED_POSITIVES
    assert S.FLOOR is MIN_REACHED_POSITIVES


def test_the_rejections_are_counted_over_every_matched_target(
        tmp_path, monkeypatch, capsys):
    """REGRESSION. The first draft counted non-protein-level rejections only
    among targets that had already taken a role -- and an ORGANISM-level
    target takes no role, so it printed "rejected as NOT PROTEIN-LEVEL: 0"
    on a database holding one. True of the set it counted, and read as
    nothing having been dropped: the same failure as the join list that
    found fourteen dangling chains and named ten."""
    _main(monkeypatch, _db(tmp_path / "ok.db"))
    out = capsys.readouterr().out
    assert "rejected as NOT PROTEIN-LEVEL: 1 target(s), 1 quantified row(s)" in out
    assert "rejected as POLYPROTEIN targets: 1 target(s)" in out
