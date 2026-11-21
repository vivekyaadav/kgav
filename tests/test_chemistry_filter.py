"""A dye matched to a dye, and resemblance weaker than resemblance that failed.

WHAT WAS WRONG. `scripts/repurposing_candidates.py` ranks by max-Tanimoto,
which knows what a compound looks like and nothing about what it is for, and
the committed read-out of 2026-11-19 shows what that costs:

    FLUORESCEIN, 0.5217, under "APPROVED DRUGS -- orderable now" for West
    Nile and for Yellow fever, resembling ERYTHROSINE. An ophthalmic dye
    matched to a food colourant.

    GEMCITABINE, CYTARABINE, CLOFARABINE, NELARABINE, AZACITIDINE,
    FLOXURIDINE, CLADRIBINE -- cytotoxic antimetabolites in the approved
    block, unmarked. They resemble nucleoside antivirals because they ARE
    nucleosides, and in a cell-based assay a compound that kills the host
    cell lowers the viral readout without touching the virus.

    TEMAZEPAM, 0.5167 for Yellow fever: a hypnotic whose only evidence is a
    Tanimoto under the median that Yellow fever's own screened-and-INACTIVE
    compounds reach.

WHAT IS PINNED HERE, AND WHY IT IS NOT A LIST OF NAMES. Three properties,
each derived on both sides:

  THE ATC RULE NEEDS *EVERY* CLASS. TOLBUTAMIDE is a diagnostic (V04CA01)
  AND a sulfonylurea (A10BB03). So are FOLIC ACID, DEXTROSE, SORBITOL and
  GONADORELIN. "Any class is diagnostic" dropped all five.

  WITHDRAWAL IS NOT A FILTER. `withdrawn_flag` reads 1 for CODEINE and
  THALIDOMIDE in this same database. A filter resting on a field that does
  not mean what its name says is the mirror of "a text marker is not an
  enforcement mechanism", and this suite derives both facts from ChEMBL
  rather than trusting the comment that states them.

  THE FLOOR CANNOT BE A CONSTANT. The median resemblance a known inactive
  achieves runs from 0.21 on MERS-CoV to 0.73 on HCoV-OC43. A flat 0.55 --
  which is what the review asked for -- sits above MERS-CoV's entire
  inactive distribution and below HCoV-OC43's median, so it is wrong in
  both directions at once.
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_spec = importlib.util.spec_from_file_location(
    "repurposing_candidates", ROOT / "scripts" / "repurposing_candidates.py")
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)

RESULTS = ROOT / "data" / "results" / "repurposing_candidates.json"
RELEASE = ROOT / "data" / "releases" / "v0.1"

needs_db = pytest.mark.skipif(not R.CHEMBL_DB.exists(),
                              reason="no raw ChEMBL database on this checkout")
needs_release = pytest.mark.skipif(not (RELEASE / "nodes.jsonl").exists(),
                                   reason="no assembled release")
needs_results = pytest.mark.skipif(not RESULTS.exists(),
                                   reason="no read-out on disk")


@pytest.fixture(scope="module")
def artifact() -> dict:
    if not RESULTS.exists():
        pytest.skip("no read-out on disk")
    return json.loads(RESULTS.read_text())


@pytest.fixture(scope="module")
def props() -> dict:
    if not (RELEASE / "nodes.jsonl").exists():
        pytest.skip("no assembled release")
    return R.molecules(RELEASE)


@pytest.fixture(scope="module")
def ann(props) -> dict:
    if not R.CHEMBL_DB.exists():
        pytest.skip("no raw ChEMBL database")
    return R.pool_annotations(props)


# ------------------------------------------------------- the role rule
def _ann(atc=(), withdrawn=False, black_box=False) -> dict:
    return {"name": "X", "atc": list(atc), "withdrawn": withdrawn,
            "black_box": black_box}


def test_a_compound_with_one_therapeutic_class_is_not_a_diagnostic_agent():
    """THE RULE, AS A UNIT. Dual-purpose compounds are the whole reason
    it reads every class instead of any."""
    assert R.ROLE_DIAGNOSTIC in R.compound_flags(_ann(["V08AB02"]))
    assert R.ROLE_DIAGNOSTIC in R.compound_flags(_ann(["S01JA01", "S01JA51"]))
    assert R.ROLE_DIAGNOSTIC not in R.compound_flags(
        _ann(["V04CA01", "A10BB03"])), (
        "a provocation test that is also a sulfonylurea is not a dye")
    assert R.ROLE_DIAGNOSTIC not in R.compound_flags(_ann([])), (
        "no ATC class at all is not evidence of being a diagnostic agent")
    assert R.ROLE_DIAGNOSTIC not in R.compound_flags(None)


def test_an_unannotated_compound_claims_no_role_at_all():
    """Unchecked is not clean. An absent annotation must produce no flags
    rather than a default."""
    assert R.compound_flags(None) == ()
    assert R.compound_flags({}) == ()


def test_the_other_roles_are_reported():
    assert R.ROLE_CYTOTOXIC in R.compound_flags(_ann(["L01BC05"]))
    assert R.ROLE_WITHDRAWN in R.compound_flags(_ann(withdrawn=True))
    assert R.ROLE_BLACK_BOX in R.compound_flags(_ann(black_box=True))
    # therapeutic_flag is only read for an approved compound.
    nt = _ann()
    nt["therapeutic"] = False
    assert R.ROLE_NO_THERAPEUTIC_USE in R.compound_flags(nt, approved=True)
    assert R.ROLE_NO_THERAPEUTIC_USE not in R.compound_flags(nt), (
        "therapeutic_flag=0 reads on 18,251 of 19,884 compounds, "
        "BALAPIRAVIR among them; outside the approved tier it is noise")


@needs_db
def test_the_dual_purpose_compounds_are_real_in_this_database(ann):
    """DERIVED, NOT ASSERTED. The rule's justification is five specific
    compounds; if ChEMBL stopped giving them a therapeutic class the rule
    would be stricter than its reason, and that should fail here rather
    than be discovered in a list."""
    by_name = {a["name"]: a for a in ann.values() if a.get("name")}
    for name in ("TOLBUTAMIDE", "FOLIC ACID", "DEXTROSE", "SORBITOL",
                 "GONADORELIN"):
        a = by_name.get(name)
        if a is None:
            continue
        assert any(c.startswith(p) for c in a["atc"]
                   for p in R.NON_THERAPEUTIC_ATC), (
            f"{name} no longer carries a diagnostic ATC class, so it is no "
            f"longer evidence that the rule must read every class")
        assert R.ROLE_DIAGNOSTIC not in R.compound_flags(a), (
            f"{name} is classified as a diagnostic agent and would be "
            f"dropped; it is a therapy as well")


@needs_db
def test_withdrawal_cannot_be_a_filter_and_the_database_says_why(ann):
    """The check that was worth running BEFORE writing the filter."""
    by_name = {a["name"]: a for a in ann.values() if a.get("name")}
    still_prescribed = [n for n in ("CODEINE", "THALIDOMIDE")
                        if n in by_name]
    assert still_prescribed, (
        "neither CODEINE nor THALIDOMIDE is in the release, so this "
        "suite cannot show why withdrawal is not a filter")
    for n in still_prescribed:
        assert R.ROLE_WITHDRAWN in R.compound_flags(by_name[n]), (
            f"{n} no longer reads withdrawn; the reason withdrawal is only "
            f"reported has gone, and the rule should be revisited rather "
            f"than left resting on a comment")


@needs_db
@needs_results
def test_a_withdrawn_compound_can_still_be_suggested(artifact, ann):
    """THE ENFORCEMENT, not the intention: if the flag had quietly become
    a filter, no row would carry it."""
    flagged = [r for cell in artifact["viruses"].values()
               for r in cell["suggestions"]
               if R.ROLE_WITHDRAWN in (r.get("flags") or [])]
    assert flagged, (
        "no suggestion is flagged withdrawn: either the flag became a "
        "filter, or it stopped being computed")


# --------------------------------------------- the InChIKey half of the pool
@needs_db
@needs_release
def test_the_role_check_covers_the_compounds_with_no_chembl_id(props, ann):
    """THE GAP A FILTER KEYED ON THE CONVENIENT COLUMN WOULD HAVE HAD.
    5,359 of 25,243 molecules carry no chembl_id, and the dyes are
    concentrated there: phenol red, indocyanine green,
    sulfobromophthalein and most of the iodinated contrast series.

    Derived: the InChIKey route must add diagnostic agents that the
    chembl_id route alone cannot see. If it adds none, either the release
    gained ids for all of them -- fine, and this test says so -- or the
    second lookup stopped working."""
    by_id = {n for n, a in ann.items() if a.get("resolved_by") == "chembl_id"}
    by_key = {n for n, a in ann.items() if a.get("resolved_by") == "inchikey"}
    no_id = {n for n, p in props.items() if not p.get("chembl_id")}
    if not no_id:
        pytest.skip("every molecule now carries a chembl_id")
    assert by_key, (
        f"{len(no_id):,} molecules carry no chembl_id and not one was "
        f"resolved by structure, so their clinical role is unchecked")
    extra = [n for n in by_key if R.ROLE_DIAGNOSTIC in R.compound_flags(ann[n])]
    assert extra, (
        "the InChIKey route found no diagnostic agent; the 17 it found on "
        "v0.1 included phenol red and indocyanine green, so this is a "
        "regression unless the release changed")
    assert by_id, "the chembl_id route resolved nothing"


@pytest.mark.slow
@needs_db
@needs_release
def test_unannotated_compounds_are_counted_and_not_treated_as_clean(
        props, ann):
    """A filter whose input covers part of the pool must say which part."""
    unchecked = set(props) - set(ann)
    out = _run(["--virus", "NCBITaxon:11060"])
    assert f"{len(unchecked):,} could NOT" in out.stdout, (
        f"{len(unchecked):,} molecules could not be role-checked and the "
        f"run does not say so")
    assert "unchecked, not clean" in out.stdout


# ------------------------------------------------------------- the floor
@needs_release
@needs_db
def test_the_floor_is_each_cells_own_and_a_constant_would_be_wrong(
        props, artifact):
    """RE-DERIVED FROM THE RELEASE, cell by cell, and then the reason a
    number was refused: the spread."""
    from kgav.labels import build_labels
    from kgav.similarity import fingerprints

    labels = build_labels(RELEASE)
    fps, _ = fingerprints({i: (p.get("smiles") or "")
                           for i, p in sorted(props.items())})
    floors = {}
    for curie, cell in artifact["viruses"].items():
        want = ((cell.get("chemistry_filter") or {})
                .get("resemblance_floor") or {})
        neg = sorted(labels.negatives.get(curie, set()))
        got = R.resemblance_floor(
            R.failed_resemblance(
                sorted(labels.positives.get(curie, set())), neg, fps, props),
            len(neg))
        assert got["floor"] == want.get("floor"), (
            f"{cell['name']}: the file records floor {want.get('floor')} "
            f"and the release gives {got['floor']}")
        assert got["n_inactives_scored"] == want.get("n_inactives_scored")
        if got["floor"] is not None:
            floors[cell["name"]] = got["floor"]
    assert len(floors) > 1
    lo, hi = min(floors.values()), max(floors.values())
    assert hi / lo > 2, (
        f"the floors now span {lo:.3f} to {hi:.3f}, under a factor of two, "
        f"so a single typed threshold would no longer be wrong in both "
        f"directions and the per-cell derivation needs its reason restated")
    # AND THE CONSTANT THE REVIEW ASKED FOR, measured against them.
    above = sorted(n for n, f in floors.items() if f > 0.55)
    below = sorted(n for n, f in floors.items() if f < 0.55)
    assert above and below, (
        f"0.55 is no longer above some cells' floors and below others': "
        f"above {above}, below {below}")


@needs_results
def test_no_suggestion_sits_below_its_own_cells_floor(artifact):
    """THE PROPERTY THE FILTER EXISTS FOR, checked on the file."""
    bad = []
    for cell in artifact["viruses"].values():
        fl = (cell.get("chemistry_filter") or {}).get("resemblance_floor") or {}
        if not fl.get("applied"):
            continue
        for r in cell["suggestions"]:
            if r["tanimoto"] < fl["floor"]:
                bad.append((cell["name"], r.get("name") or r["chembl_id"],
                            r["tanimoto"], fl["floor"]))
    assert not bad, f"rows below their cell's floor: {bad[:6]}"


@needs_results
def test_the_floor_names_what_it_cost_rather_than_counting_it(artifact):
    """"the floor removed 7 rows" cannot be checked against the previous
    list; "TEMAZEPAM 0.517" can."""
    seen = 0
    for cell in artifact["viruses"].values():
        fl = (cell.get("chemistry_filter") or {}).get("resemblance_floor") or {}
        rows = fl.get("rows_displaced")
        if rows is None:
            continue
        assert fl.get("n_rows_displaced") == len(rows), cell["name"]
        for r in rows:
            assert r["tanimoto"] < fl["floor"], cell["name"]
            assert r.get("candidate"), "a displaced row with no identity"
        seen += len(rows)
    assert seen, (
        "the floor displaced no row in any cell, so nothing here shows it "
        "doing anything")


# ------------------------------------------------- the rows the review named
@needs_results
def test_no_diagnostic_agent_survives_into_any_list(artifact):
    bad = [(cell["name"], r.get("name"))
           for cell in artifact["viruses"].values()
           for r in cell["suggestions"]
           if R.ROLE_DIAGNOSTIC in (r.get("flags") or [])]
    assert not bad, f"diagnostic agents in the lists: {bad}"


@needs_results
def test_the_dye_the_review_quoted_is_gone_for_a_recorded_reason(artifact):
    """A row that merely vanished is indistinguishable from one the scorer
    stopped reaching, so the file has to say which mechanism took it."""
    names = {(r.get("name") or "").upper()
             for cell in artifact["viruses"].values()
             for r in cell["suggestions"]}
    dropped = {n.upper() for cell in artifact["viruses"].values()
               for n in ((cell.get("chemistry_filter") or {})
                         .get("diagnostic_agents_dropped") or [])}
    assert "FLUORESCEIN" not in names, "the dye is still being suggested"
    assert "FLUORESCEIN" in dropped, (
        "FLUORESCEIN is absent from the lists and absent from the drop "
        "record, so nothing says why it went")


@needs_results
def test_temazepam_survives_and_the_run_says_where_it_sits(artifact):
    """THE ROW THE REVIEW CALLED TAIL NOISE, AND THE ANSWER THE DATA GAVE.

    TEMAZEPAM sits at 0.5167 for Yellow fever. A flat 0.55 floor would
    have removed it. Yellow fever's OWN failures have a median of 0.3008,
    so by the only yardstick in the same units the row is not tail noise:
    it resembles the known actives better than three quarters of the
    compounds screened against this virus that came back negative.

    Bending the threshold until a disliked row falls out is how a filter
    becomes an opinion. What the row gets instead is its position in that
    distribution, in every cell, so a reader can see the weakness rather
    than trust that it was removed.
    """
    yf = artifact["viruses"].get("NCBITaxon:11089")
    if not yf:
        pytest.skip("Yellow fever is not in this read-out")
    row = next((r for r in yf["suggestions"]
                if (r.get("name") or "").upper() == "TEMAZEPAM"), None)
    if row is None:
        pytest.skip("TEMAZEPAM is no longer ranked for Yellow fever")
    fl = yf["chemistry_filter"]["resemblance_floor"]
    assert row["tanimoto"] > fl["floor"], (
        "TEMAZEPAM is below the floor and still in the list")
    assert row["beats_failures"] is not None, (
        "the row carries no position in the failure distribution, so its "
        "weakness is unstated")
    assert row["beats_failures"] > 0.5, (
        "a row above the median of the failures must beat more than half "
        "of them; the floor and the percentile disagree")


@needs_results
def test_every_row_says_where_it_sits_among_the_failures(artifact):
    """DERIVED BOTH WAYS. The floor is the 50% line on `beats_failures`,
    so every shown row must beat at least half the failures -- and if the
    two were computed from different distributions, this is where it
    shows."""
    for cell in artifact["viruses"].values():
        fl = (cell.get("chemistry_filter") or {}).get("resemblance_floor") or {}
        for r in cell["suggestions"]:
            assert r.get("beats_failures") is not None, (
                f"{cell['name']}/{r.get('name')} carries no position in "
                f"the failure distribution")
            if fl.get("applied"):
                assert r["beats_failures"] >= 0.5, (
                    f"{cell['name']}/{r.get('name')}: Tc {r['tanimoto']} is "
                    f"above the floor {fl['floor']} and yet beats only "
                    f"{r['beats_failures']:.0%} of the failures, so the "
                    f"floor and the percentile are not the same "
                    f"distribution")


@needs_results
def test_an_approved_product_with_no_therapeutic_use_says_so(artifact):
    """PHENOLPHTHALEIN and ERYTHROSINE are what "APPROVED DRUGS --
    orderable now" was hiding. ChEMBL's therapeutic_flag is useless as a
    filter -- it reads 0 for BALAPIRAVIR and OBELDESIVIR, antivirals in
    development -- but on a compound the release flags APPROVED it is a
    precise statement, and it is the label those rows needed."""
    flagged = [(cell["name"], r.get("name"))
               for cell in artifact["viruses"].values()
               for r in cell["suggestions"]
               if R.ROLE_NO_THERAPEUTIC_USE in (r.get("flags") or [])]
    assert flagged, (
        "no approved row is flagged as having no recorded therapeutic use")
    # AND IT IS CONFINED TO THE APPROVED TIER, where it means something.
    wrong = [(cell["name"], r.get("name"), r["phase_tier"])
             for cell in artifact["viruses"].values()
             for r in cell["suggestions"]
             if R.ROLE_NO_THERAPEUTIC_USE in (r.get("flags") or [])
             and r["phase_tier"] != "approved"]
    assert not wrong, (
        f"the flag escaped the approved tier, where therapeutic_flag=0 "
        f"means little: {wrong[:5]}")


@needs_results
def test_the_cytotoxics_are_flagged_and_still_there(artifact):
    """NOT DROPPED, ON PURPOSE: nucleoside analogues are an antiviral
    family and oncology-to-antiviral is a real repurposing direction. The
    requirement is that a reader never meets one unlabelled."""
    cyt = [(cell["name"], r.get("name"))
           for cell in artifact["viruses"].values()
           for r in cell["suggestions"]
           if R.ROLE_CYTOTOXIC in (r.get("flags") or [])]
    assert cyt, "no suggestion is flagged cytotoxic; the label stopped working"


@needs_db
@needs_results
def test_every_antineoplastic_in_the_lists_carries_the_flag(artifact, ann):
    """DERIVED FROM CHEMBL, the other way round: the flag must be on every
    row ChEMBL classes L01, not only on the ones that were noticed."""
    missed = []
    for cell in artifact["viruses"].values():
        for r in cell["suggestions"]:
            a = ann.get(r["candidate"])
            if not a:
                continue
            is_l01 = any(c.startswith(R.CYTOTOXIC_ATC) for c in a["atc"])
            if is_l01 != (R.ROLE_CYTOTOXIC in (r.get("flags") or [])):
                missed.append((cell["name"], a.get("name"), a["atc"]))
    assert not missed, f"the cytotoxic flag disagrees with ChEMBL: {missed[:6]}"


# ------------------------------------------- the run, against its own file
def _run(extra: list[str]) -> subprocess.CompletedProcess:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "r.json"
        r = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "repurposing_candidates.py"),
             "--out", str(out), *extra],
            capture_output=True, text=True, cwd=ROOT, timeout=2400)
        assert r.returncode == 0, r.stdout[-4000:] + r.stderr[-4000:]
        r.payload = json.loads(out.read_text())     # type: ignore[attr-defined]
        return r


@pytest.mark.slow
@needs_release
@needs_db
def test_the_printout_and_the_written_file_agree_about_the_filtering():
    """The defect class this script already names: a file disagreeing with
    the run that wrote it is worse than either being wrong."""
    r = _run(["--virus", "NCBITaxon:11060"])
    cell = r.payload["viruses"]["NCBITaxon:11060"]
    cf = cell["chemistry_filter"]
    assert f"{cf['n_diagnostic_agents_dropped']:,} dropped as diagnostic" \
        in r.stdout
    assert len(cf["diagnostic_agents_dropped"]) == \
        cf["n_diagnostic_agents_dropped"]
    fl = cf["resemblance_floor"]
    assert f"RESEMBLANCE FLOOR {fl['floor']:.4f}" in r.stdout
    assert f"cost the list {fl['n_rows_displaced']} row(s)" in r.stdout
    assert f"{fl['n_candidates_below_floor']:,} of " in r.stdout


@pytest.mark.slow
@needs_release
@needs_db
def test_both_filters_can_be_switched_off_and_the_difference_measured():
    """A filter whose effect cannot be computed is a claim. The off
    switches exist for this measurement and are asserted to change the
    output, so they cannot rot into no-ops."""
    virus = "NCBITaxon:11089"          # Yellow fever: the TEMAZEPAM cell
    on = _run(["--virus", virus]).payload["viruses"][virus]
    off = _run(["--virus", virus, "--keep-diagnostic-agents",
                "--no-resemblance-floor"]).payload["viruses"][virus]
    assert off["chemistry_filter"]["n_diagnostic_agents_dropped"] == 0
    assert off["chemistry_filter"]["resemblance_floor"]["applied"] is False
    assert on["chemistry_filter"]["n_diagnostic_agents_dropped"] > 0
    kept = {r["candidate"] for r in on["suggestions"]}
    unfiltered = {r["candidate"] for r in off["suggestions"]}
    assert unfiltered - kept, (
        "switching both filters off changes nothing, so neither is doing "
        "anything on this cell")
    assert kept - unfiltered, (
        "the filters removed rows and nothing took their slots, so the "
        "floor is being applied after the tier cap rather than before")


# ------------------------------------------- the agent reads the same roles
def test_the_agents_role_markers_are_the_scripts_own_constants():
    """A CONTRACT ACROSS A SUBPROCESS BOUNDARY.

    The orchestrator shells out to this script -- it is not importable
    from src/ -- and then filters the roles it will show by name. Two
    copies of a string, one on each side of a subprocess, is exactly the
    kind of contract a caller can forget: renaming ROLE_CYTOTOXIC would
    leave the agent silently showing no marker at all, and the lists
    would look clean.
    """
    from kgav.agent.orchestrator import Orchestrator

    known = {R.ROLE_DIAGNOSTIC, R.ROLE_CYTOTOXIC, R.ROLE_WITHDRAWN,
             R.ROLE_BLACK_BOX}
    unknown = set(Orchestrator.SHORTLIST_FLAGS) - known
    assert not unknown, (
        f"the agent marks {sorted(unknown)}, which "
        f"scripts/repurposing_candidates.py does not produce")
    assert R.ROLE_CYTOTOXIC in Orchestrator.SHORTLIST_FLAGS, (
        "the agent would hand over a cytotoxic antimetabolite unmarked")
    assert set(Orchestrator.SHORTLIST_TIERS) <= set(R.TIERS)


@needs_db
def test_the_therapeutic_flag_is_unusable_outside_the_approved_tier(ann):
    """DERIVED, because it is the reason the flag is confined. If ChEMBL
    started flagging development antivirals therapeutic, the restriction
    would be stricter than its justification."""
    by_name = {a["name"]: a for a in ann.values() if a.get("name")}
    seen = [n for n in ("BALAPIRAVIR", "OBELDESIVIR", "ALISPORIVIR")
            if n in by_name]
    assert seen, (
        "none of the development antivirals that justify the restriction "
        "is in the release")
    for n in seen:
        assert by_name[n].get("therapeutic") is False, (
            f"{n} now reads therapeutic in ChEMBL, so therapeutic_flag=0 "
            f"may no longer be noise outside the approved tier")


# ---------------------------------------- the documents, against the run
SCOPE = ROOT / "docs" / "scope_and_limits.md"
DENV2_PREREG = ROOT / "docs" / "prereg_denv2_prospective_screen.md"


def _cells_by_name(artifact: dict) -> dict:
    return {c["name"]: c for c in artifact["viruses"].values()}


@needs_results
def test_every_figure_the_scope_document_quotes_is_the_runs_own(artifact):
    """A LIMIT RESTATED WITH ITS OWN COPY OF THE NUMBERS IS A SECOND
    THING TO GO STALE, which is the failure that document has already had
    twice. Each figure below is read out of the artifact and looked for.
    """
    doc = " ".join(SCOPE.read_text().split())
    by = _cells_by_name(artifact)
    pool = artifact["chemistry_filter"]

    mers = by["MERS-CoV"]["chemistry_filter"]["resemblance_floor"]["floor"]
    oc43 = by["HCoV-OC43"]["chemistry_filter"]["resemblance_floor"]["floor"]
    yf = by["Yellow fever"]["chemistry_filter"]["resemblance_floor"]
    assert f"**{mers:.4f}** on MERS-CoV" in doc, (
        f"the scope document does not quote MERS-CoV's floor {mers:.4f}")
    assert f"**{oc43:.4f}** on HCoV-OC43" in doc
    assert f"median of {yf['floor']:.4f}" in doc
    for n, needle in (
            (pool["n_diagnostic_identified"], "compounds in the pool"),
            (pool["n_unchecked"], "resolve to no ChEMBL record"),
            (pool["n_molecules"], "molecules in the release"),
            (pool["n_diagnostic_identified_by_inchikey"],
             "carry no chembl_id in the release")):
        assert f"{n:,}" in doc, (
            f"the scope document does not quote {n:,} ({needle})")
    tz = next((r for r in by["Yellow fever"]["suggestions"]
               if (r.get("name") or "").upper() == "TEMAZEPAM"), None)
    if tz:
        assert f"{tz['beats_failures']:.0%} of this virus's failures" in doc
    # And the ATC classes are the code's, not a second list.
    for code in R.NON_THERAPEUTIC_ATC:
        assert code in doc, f"the scope document omits ATC {code}"


@needs_results
def test_the_denv2_amendments_pool_figures_are_the_runs_own(artifact):
    """The amendment it supersedes published a hand-typed pool size, which
    this filter moved. This one is derived, in both directions: the new
    pool and the gap from the old."""
    doc = " ".join(DENV2_PREREG.read_text().split())
    cell = artifact["viruses"]["NCBITaxon:11060"]
    dropped = cell["chemistry_filter"]["n_diagnostic_agents_dropped"]
    floor = cell["chemistry_filter"]["resemblance_floor"]
    assert f"**{dropped}** fewer" in doc, (
        f"{dropped} compounds left DENV-2's pool and the amendment says "
        f"otherwise")
    assert f"**{floor['floor']:.4f}**" in doc
    assert f"displaces **{floor['n_rows_displaced']} rows**" in doc or \
        f"displaces **{floor['n_rows_displaced']} row" in doc
    # The two cells it quotes as costing rows must still cost them.
    by = _cells_by_name(artifact)
    for name in ("TBEV", "HCoV-229E"):
        n = (by[name]["chemistry_filter"]["resemblance_floor"]
             ["n_rows_displaced"])
        assert f"{n} on {name}" in doc, (
            f"the amendment says the floor costs {name} a different "
            f"number of rows than the {n} it now costs")


@pytest.mark.slow
@needs_release
@needs_db
def test_the_denv2_draw_pool_is_what_the_amendment_publishes():
    """THE FIGURE THE PREVIOUS AMENDMENT HAND-TYPED, now derived from the
    recorded draw command itself.

    It is NOT this cell's `n_candidates` in the committed read-out:
    `--exclude-labelled` removes every compound carrying any DENV-2
    activity edge, which is stricter than the label sets, so the draw's
    pool is smaller than the browsing read-out's and reading one off the
    other would publish the wrong number confidently.
    """
    r = _run(["--virus", "NCBITaxon:11060", "--exclude-labelled",
              "--flat", "--top", "20"])
    pool = r.payload["viruses"]["NCBITaxon:11060"]["n_candidates"]
    doc = " ".join(DENV2_PREREG.read_text().split())
    assert f"**{pool:,}** of 25,243" in doc, (
        f"the recorded draw command leaves {pool:,} compounds in the pool "
        f"and the amendment publishes another figure")


@pytest.mark.slow
@needs_release
@needs_db
def test_the_denv2_arm_a_is_identical_with_the_filters_off():
    """THE CLAIM THE AMENDMENT RESTS ON. If either filter touched the head
    of DENV-2's ranking, the registered arm A would have changed under a
    screen that is already signed off on its method, and "arm A is
    unchanged" would be a sentence rather than a measurement."""
    virus = "NCBITaxon:11060"
    a = _run(["--virus", virus, "--exclude-labelled", "--flat", "--top",
              "20"]).payload["viruses"][virus]["suggestions"]
    b = _run(["--virus", virus, "--exclude-labelled", "--flat", "--top",
              "20", "--keep-diagnostic-agents",
              "--no-resemblance-floor"]).payload["viruses"][virus][
                  "suggestions"]
    assert [r["candidate"] for r in a] == [r["candidate"] for r in b], (
        "the chemistry filter changed the registered arm A; the amendment "
        "says it did not")
    doc = " ".join(DENV2_PREREG.read_text().split())
    assert "the same twenty compounds in the same order" in doc
    assert len(a) == 20
