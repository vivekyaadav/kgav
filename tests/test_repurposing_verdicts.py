"""A candidate list must carry the verdict on the cell that produced it.

WHAT WAS WRONG. Until 2026-11-12 `scripts/repurposing_candidates.py` emitted
an identically-shaped ranked list for every virus with ten known actives, and
the file carried `n_known_actives`, `n_candidates` and NO measure of whether
the scorer works on that virus. So Yellow fever at AUC 0.8117, which clears
its corrected threshold, and DENV-4 at **0.4792**, where the scorer ranks
known actives UNDER known inactives, printed the same way -- fifteen named
drugs each. A reader could not tell a result from an anti-ranking.

The header did say "clears in ten of fourteen", but that is a statement about
the FAMILY and the reader is looking at one row. It also said "none below",
which had stopped being true: two cells are below chance.

WHAT IS PINNED HERE. That a below-chance cell produces NO suggestions rather
than suggestions with a caveat -- a caveat beside fifteen named drugs loses
to the drugs -- that the written file and the printed run agree about which
those are, and that the banner's counts come from the results files rather
than from a literal.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_spec = importlib.util.spec_from_file_location(
    "repurposing_candidates", ROOT / "scripts" / "repurposing_candidates.py")
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)

RESULTS = ROOT / "data" / "results" / "repurposing_candidates.json"


def scorer_file(root: Path, name: str, cells: dict) -> None:
    """A results file shaped like the registered runs'."""
    d = root / "data" / "results"
    d.mkdir(parents=True, exist_ok=True)
    (d / name).write_text(json.dumps({"primary_chembl": {
        v: {"readable": True, "k_cell": 12, "structure_max": m}
        for v, m in cells.items()}}))


def cell(auc, *, lo, hi, clears, direction="above") -> dict:
    return {"auc": auc, "auc_ci_lo_corrected": lo,
            "auc_ci_hi_corrected": hi,
            "clears_corrected_alpha": clears, "direction": direction}


# ------------------------------------------------------------- the tiers
@pytest.mark.parametrize("auc,clears,direction,want", [
    (0.8117, True, "above", "clears"),
    (0.5903, False, "above", "weak"),
    (0.4792, False, "below", "below"),
    (0.4926, False, "above", "below"),      # below 0.5 whatever the label
    (0.4792, True, "below", "below"),       # clearing BELOW chance is not a
])                                           # result, it is an anti-ranking
def test_the_tier_is_the_cells_own_verdict(auc, clears, direction, want):
    assert R.tier(cell(auc, lo=auc - .1, hi=auc + .1, clears=clears,
                       direction=direction)) == want


def test_a_cell_with_no_registered_run_is_unknown_not_usable():
    assert R.tier(None) == "unknown"
    assert R.tier({}) == "unknown"


# ------------------------------------------------------- reading the runs
def test_the_primary_control_is_the_sternest_one_that_ran(tmp_path):
    """The strength-matched control removes each candidate's k_cell nearest
    references, so quoting it is the conservative choice. The other two
    ride along, because a reader deciding whether to chase a row deserves
    the range rather than the one I picked."""
    v = "NCBITaxon:11070"
    scorer_file(tmp_path, "structure_scorer_strength_matched.json",
                {v: cell(0.4792, lo=0.37, hi=0.58, clears=False)})
    scorer_file(tmp_path, "structure_scorer_symmetric.json",
                {v: cell(0.5868, lo=0.48, hi=0.69, clears=False)})
    scorer_file(tmp_path, "structure_scorer.json",
                {v: cell(0.4644, lo=0.36, hi=0.57, clears=False)})
    got = R.verdicts(tmp_path)[v]
    assert got["control"] == R.PRIMARY_RUN
    assert got["auc"] == 0.4792, "the primary must be the strength-matched run"
    assert got["auc_symmetric_k10"] == 0.5868
    assert got["auc_publication"] == 0.4644
    assert R.tier(got) == "below"


def test_a_cell_the_primary_run_lacks_is_not_promoted_by_a_weaker_one(
        tmp_path):
    """A virus readable only under a lenient control must not inherit that
    control's verdict: `auc` is set by the primary alone."""
    v = "NCBITaxon:99999"
    scorer_file(tmp_path, "structure_scorer_strength_matched.json", {})
    scorer_file(tmp_path, "structure_scorer_symmetric.json",
                {v: cell(0.90, lo=0.85, hi=0.95, clears=True)})
    got = R.verdicts(tmp_path).get(v, {})
    assert got.get("auc") is None, (
        "a cell absent from the sternest run must not take its AUC from a "
        "weaker one")
    assert R.tier(got) == "unknown"


# ---------------------------------------------------------- the banner
def test_the_banner_counts_come_from_the_runs_and_not_from_a_literal(
        tmp_path):
    """The hand-typed banner said "none below" while two cells were below
    chance, and quoted 0.7841 for a cell reading 0.7822."""
    scorer_file(tmp_path, "structure_scorer_strength_matched.json", {
        "NCBITaxon:2697049": cell(0.7822, lo=0.76, hi=0.80, clears=True),
        "NCBITaxon:11089": cell(0.8117, lo=0.74, hi=0.88, clears=True),
        "NCBITaxon:31631": cell(0.5903, lo=0.47, hi=0.71, clears=False),
        "NCBITaxon:11070": cell(0.4792, lo=0.37, hi=0.58, clears=False,
                                direction="below"),
    })
    text = " ".join(R.headline(R.verdicts(tmp_path), tmp_path))
    assert "0.7822" in text, "the SARS-CoV-2 figure is not the run's"
    assert "2 of 4 viruses" in text, f"the clearing count is wrong: {text}"
    assert "1 above chance but not clearing" in text
    assert "1 BELOW chance" in text
    assert "none below" not in text


# -------------------------------------------- the committed read-out itself
@pytest.mark.skipif(not RESULTS.exists(), reason="no read-out on disk")
def test_no_below_chance_virus_carries_a_list_in_the_committed_file():
    """THE PROPERTY THAT MATTERS. Anything else about this file is
    presentation; a ranked list from a cell that anti-ranks is a wrong
    answer with names on it."""
    doc = json.loads(RESULTS.read_text())
    bad = []
    for curie, row in doc["viruses"].items():
        v = row.get("scorer_verdict") or {}
        assert v, (
            f"{row.get('name', curie)} carries no scorer_verdict; the list "
            f"and the measure of whether the scorer works must travel "
            f"together")
        if v.get("tier") in ("below", "unknown") and row["suggestions"]:
            bad.append((row.get("name", curie), v.get("auc"),
                        len(row["suggestions"])))
    assert not bad, (
        "these viruses are below chance or unmeasured and still carry "
        f"suggestions: {bad}")


@pytest.mark.skipif(not RESULTS.exists(), reason="no read-out on disk")
def test_every_verdict_in_the_file_is_the_scorer_runs_own():
    """Derived on both sides: the file's verdict against the registered
    run, so the read-out cannot drift from the study behind it."""
    doc = json.loads(RESULTS.read_text())
    live = R.verdicts(ROOT)
    if not live:
        pytest.skip("no registered scorer run on disk")
    for curie, row in doc["viruses"].items():
        v = row.get("scorer_verdict") or {}
        want = live.get(curie)
        if want is None or want.get("auc") is None:
            assert v.get("tier") == "unknown", (
                f"{row.get('name', curie)} has no registered cell and the "
                f"file claims tier {v.get('tier')!r}")
            continue
        assert v.get("auc") == want["auc"], (
            f"{row.get('name', curie)}: the read-out says AUC "
            f"{v.get('auc')} and the registered run says {want['auc']}")
        assert v.get("tier") == R.tier(want)


@pytest.mark.skipif(not RESULTS.exists(), reason="no read-out on disk")
def test_a_flagged_near_identity_says_what_differs():
    """A Tanimoto of 1.0 against a different skeleton is a fingerprint
    limit, not an identity -- a saturated ring's size, a peptide's length,
    a residue's position. Every flagged row must carry the heavy-atom
    difference that says which."""
    doc = json.loads(RESULTS.read_text())
    flagged = [s for row in doc["viruses"].values()
               for s in row["suggestions"] if s.get("fingerprint_identical")]
    if not flagged:
        pytest.skip("no near-identity in the current read-out")
    for s in flagged:
        assert "heavy_atom_delta" in s, (
            f"{s.get('chembl_id')} is flagged identical with no measure of "
            f"how far apart the two structures are")
        assert s["tanimoto"] >= 0.9999


# ------------------------------------------------- the refusal, end to end
@pytest.mark.skipif(
    not (ROOT / "data" / "releases" / "v0.1" / "nodes.jsonl").exists(),
    reason="no release on disk")
def test_with_no_registered_run_it_refuses_rather_than_listing(
        tmp_path, monkeypatch, capsys):
    """THE BRANCH THAT MUST NOT QUIETLY DEGRADE. Without a registered run
    there is no AUC to put beside any row, and the honest output is no
    output -- not a list with an empty column, which reads exactly like a
    list with a full one."""
    monkeypatch.setattr(R, "verdicts", lambda root=ROOT: {})
    code = R.main(["--out", str(tmp_path / "out.json")])
    out = capsys.readouterr().out
    assert code == 2, "it produced a list with no verdict behind it"
    assert "REFUSED" in out
    assert "make structure-scorer-strength" in out, (
        "a refusal that does not say how to satisfy it is a dead end")
    assert not (tmp_path / "out.json").exists(), (
        "it refused and wrote a file anyway")


# ------------------------------------------- the 2026-11-19 tiering and L7
@pytest.fixture(scope="module")
def results() -> dict:
    if not RESULTS.exists():
        pytest.skip("the shortlist has not been generated on this checkout")
    return json.loads(RESULTS.read_text())


def _tiered_rows(res: dict):
    """Every suggestion in the artifact, with its cell."""
    for taxon, cell in res["viruses"].items():
        for row in cell["suggestions"]:
            yield taxon, cell, row


def test_every_suggestion_declares_a_clinical_tier(results):
    """A list that mixes an approved drug with a research analog and
    says which is which answers two questions; one that does not
    answers neither."""
    for _t, _c, row in _tiered_rows(results):
        assert row.get("phase_tier") in R.TIERS, row


def test_the_approved_tier_is_exactly_the_approved_flag(results):
    """phase_tier must not drift from the field the rest of the script
    and its tests read."""
    for _t, _c, row in _tiered_rows(results):
        assert (row["phase_tier"] == "approved") == bool(row["is_approved"])


def test_an_absent_max_phase_is_research_and_not_promoted(results):
    """20,626 of 25,243 compounds carry no max_phase. Unknown maturity
    must read as research; treating it as clinical would put the bulk of
    the pool in the tier a lab acts on first."""
    assert R._phase_tier({}) == "research"
    assert R._phase_tier({"max_phase": None}) == "research"
    assert R._phase_tier({"max_phase": "not a number"}) == "research"
    assert R._phase_tier({"max_phase": 0}) == "research"
    assert R._phase_tier({"max_phase": 2}) == "clinical"
    assert R._phase_tier({"max_phase": 4}) == "clinical"   # flag decides
    assert R._phase_tier({"is_approved": True}) == "approved"


def test_the_pool_is_not_restricted_to_approved_compounds(results):
    """THE POINT OF TIERING RATHER THAN FILTERING. Every small molecule
    is still scored; what changed is that the window no longer belongs
    to whichever tier happens to rank highest."""
    assert results["approved_only"] is False
    tiers = {row["phase_tier"] for _t, _c, row in _tiered_rows(results)}
    assert "research" in tiers, (
        "no research compound survived, so the pool was restricted after "
        "all")


def test_each_tier_is_ranked_on_the_same_statistic(results):
    """Tiering changes WHICH rows are shown, never their score. Within a
    tier the order must still be descending Tanimoto."""
    for taxon, cell in results["viruses"].items():
        for tier in R.TIERS:
            tc = [r["tanimoto"] for r in cell["suggestions"]
                  if r["phase_tier"] == tier]
            assert tc == sorted(tc, reverse=True), f"{taxon}/{tier}"


def test_a_tier_is_capped_and_cannot_crowd_the_others_out(results):
    """Before tiering, SARS-CoV-2's top 60 by Tanimoto held 0 approved
    and 0 named compounds: ~25,000 research analogs outranked every drug
    a lab could order. No tier may exceed its own cap."""
    cap = results["top"]
    for taxon, cell in results["viruses"].items():
        for tier in R.TIERS:
            n = sum(1 for r in cell["suggestions"]
                    if r["phase_tier"] == tier)
            assert n <= cap, f"{taxon}/{tier}: {n} > {cap}"


def test_the_reported_graph_degree_does_not_count_held_out_labels():
    """L7, AT THE POINT THE SHORTLIST IS READ.

    DERIVED FROM THE RELEASE, not from the script. Until 2026-11-19 this
    script loaded the graph without `skip_predicates`, so
    `degree_unlabelled` was `degree` and the column counted the labels
    the shortlist exists to predict: 163 of 195 reported degrees were
    inflated, and 78 of 195 compounds the graph does not reach were
    printed as reached. A comment above the field already claimed the
    labels were excluded, which is why nothing failed.
    """
    from kgav.baselines import Graph
    from kgav.schema import load_schema

    release = ROOT / "data" / "releases" / "v0.1"
    if not (release / "edges.jsonl").exists():
        pytest.skip("no assembled release")
    res = json.loads(RESULTS.read_text())
    schema = load_schema()
    g = Graph.load(release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    wrong = [(c["name"], r["chembl_id"], r["graph_degree"],
              g.degree_unlabelled.get(r["candidate"], 0))
             for _t, c, r in _tiered_rows(res)
             if r["graph_degree"] != g.degree_unlabelled.get(r["candidate"], 0)]
    assert not wrong, (
        "reported graph_degree disagrees with the label-free degree for "
        f"{len(wrong)} suggestion(s); first few: {wrong[:5]}")
