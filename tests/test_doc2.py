"""doc2/ describes the code it is generated from, and nothing in it is typed.

WHAT THIS FOLDER IS. A from-scratch description of the whole repository,
written to one requirement: nothing in it should be wrong. A hand-written
document cannot meet that for longer than one rebuild, and this repository has
the receipts -- a typed 0.7841 beside a cell reading 0.7822, a "none below
chance" in a banner for a family with two below chance, thirteen stale figures
in a hand-drawn metagraph that had also silently lost two predicates, and a
pool size of 21,081 that a filter moved to 21,039 on the day it was written.
Every one of those was true when it was typed.

So the prose lives in doc2/_templates/*.md.in and every figure is interpolated
at render time from the schema, the register, the release MANIFEST, the
artifacts in data/results and the constants imported from src/kgav.

WHAT THESE TESTS PIN, which is not "the documents exist":

  THE GENERATOR CANNOT RENDER A BLANK. A placeholder with no value raises
  rather than producing an empty string, because a figure that renders blank
  reads like prose that was always meant to read that way.

  --check CAN FAIL. The gate is only worth having if it reports a stale
  document, so one is mutated and the check is required to notice -- the
  defect this repository names most often is a check that cannot fail.

  THE FIGURES ARE THE ARTIFACTS'. Specific claims are re-derived here from
  data/results and looked for in the rendered text, so a document cannot
  drift from the run it quotes even while `--check` passes (it would pass:
  --check compares the file against the generator, and both would move
  together if a figure were hardcoded in a TEMPLATE).

  THE SIGNED REGISTRATIONS ARE NOT REPRODUCED. docs/prereg_*.md are amended
  by appending and never rewritten; doc2 describes and links them. A copy
  would destroy the one property that makes a registration one.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC2 = ROOT / "doc2"
TEMPLATES = DOC2 / "_templates"
RESULTS = ROOT / "data" / "results"
SCRIPT = ROOT / "scripts" / "write_doc2.py"
sys.path.insert(0, str(ROOT / "src"))

needs_artifacts = pytest.mark.skipif(
    not (ROOT / "data" / "releases" / "v0.1" / "MANIFEST.json").exists(),
    reason="no assembled release")


def _check(*extra: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), "--check", *extra],
                          capture_output=True, text=True, cwd=ROOT,
                          timeout=900)


def _text(name: str) -> str:
    return (DOC2 / name).read_text()


def _all_text() -> str:
    return "\n".join(sorted(p.read_text() for p in DOC2.glob("*.md")))


# -------------------------------------------------- the folder is generated
def test_every_template_has_a_rendered_document():
    tpls = sorted(p.name[:-3] for p in TEMPLATES.glob("*.md.in"))
    docs = sorted(p.name for p in DOC2.glob("*.md"))
    assert tpls, "no templates"
    assert tpls == docs, (
        f"templates and documents disagree: only a template "
        f"{sorted(set(tpls) - set(docs))}, only a document "
        f"{sorted(set(docs) - set(tpls))}")


@needs_artifacts
def test_the_generated_folder_is_current():
    r = _check()
    assert r.returncode == 0, r.stdout + r.stderr


def test_no_placeholder_survives_into_a_document():
    """A `{{key}}` left in the output is a figure that never rendered."""
    for p in sorted(DOC2.glob("*.md")):
        left = re.findall(r"\{\{[a-z0-9_]+\}\}", p.read_text())
        assert not left, f"{p.name} still contains {left}"


@needs_artifacts
def test_an_unknown_placeholder_is_an_error_and_not_a_blank(tmp_path):
    """THE FIRST OF THE TWO WAYS A GENERATED DOCUMENT GOES WRONG. A
    placeholder with no value must stop the render: an empty string in
    the middle of a sentence reads as prose that was always meant to read
    that way, which is strictly worse than the stale figure it replaced.
    """
    probe = TEMPLATES / "zz-probe.md.in"
    probe.write_text("{{a_value_nothing_derives}}\n")
    try:
        r = _check()
        assert r.returncode != 0
        assert "a_value_nothing_derives" in (r.stdout + r.stderr)
    finally:
        probe.unlink()


@needs_artifacts
def test_the_check_reports_a_stale_document():
    """THE SECOND WAY, and the one this repository names most often: a
    check that cannot fail is worse than no check. One line is mutated
    and restored."""
    target = DOC2 / "01-overview.md"
    original = target.read_text()
    try:
        target.write_text(original.replace("# 01 — Overview",
                                           "# 01 — Overview (edited)", 1))
        r = _check()
        assert r.returncode == 1, (
            "--check passed on a document that had been edited")
        assert "01-overview.md" in r.stdout
        assert "line " in r.stdout, "it does not say which line moved"
    finally:
        target.write_text(original)
    assert _check().returncode == 0, "the mutation was not restored"


# ------------------------------------------- the figures are the artifacts'
@needs_artifacts
def test_the_release_census_is_the_manifests():
    man = json.loads((ROOT / "data" / "releases" / "v0.1"
                      / "MANIFEST.json").read_text())
    doc = " ".join(_text("02-graph-model.md").split())
    for cls, count in man["nodes_by_class"].items():
        assert f"| {cls} | {count:,} |" in doc, (
            f"02-graph-model does not state {cls} = {count:,}")
    for pred, count in man["edges_by_predicate"].items():
        assert f"`{pred}` | {count:,}" in doc, (
            f"02-graph-model does not state {pred} = {count:,}")
    head = " ".join(_text("00-README.md").split())
    assert f"{man['nodes']:,} nodes" in head
    assert f"{man['edges']:,} edges" in head


@needs_artifacts
def test_every_structure_scorer_auc_is_in_the_results_document():
    """THE CLAIM-BEARING TABLE. Each of the three controls publishes one
    AUC per readable cell, and all three tables are in 07-results."""
    doc = _text("07-results.md")
    seen = 0
    for name in ("structure_scorer.json", "structure_scorer_symmetric.json",
                 "structure_scorer_strength_matched.json"):
        d = json.loads((RESULTS / name).read_text())
        for cell in d["primary_chembl"].values():
            m = cell.get("structure_max") or {}
            if m.get("auc") is None:
                continue
            assert f"{m['auc']:.4f}" in doc, (
                f"{name} records {m['auc']:.4f} and 07-results does not "
                f"state it")
            seen += 1
    assert seen >= 40, f"only {seen} AUCs checked; the tables look truncated"


@needs_artifacts
def test_the_prospective_verdicts_are_the_runs_own():
    """SIX CLEAR AND TWO CLEAR BELOW CHANCE, and the second is a measured
    anti-ranking rather than a null. A document that rendered both as
    "below chance" would disagree with its own count."""
    d = json.loads((RESULTS / "temporal_structure_scorer.json").read_text())
    doc = _text("07-results.md")
    clears = set(d["clears"])
    below = [c for c in d["cells"].values()
             if (c.get("primary") or {}).get("direction") == "below"
             and (c.get("primary") or {}).get("clears_corrected_alpha")]
    assert f"**{len(clears)} of {d['n_evaluable']} cells clear above chance**" \
        in doc
    assert f"**{len(below)} clear below chance**" in doc
    for cell in d["cells"].values():
        auc = (cell.get("primary") or {})["auc"]
        assert f"{auc:.4f}" in doc, f"{cell['name']} {auc:.4f} is missing"
    # And the three states are distinguished.
    assert "clears BELOW chance" in doc
    assert "no verdict" in doc


@needs_artifacts
def test_the_capability_verdict_is_the_shared_predicate():
    """doc2 must not re-decide what licenses a shortlist. The predicate
    lives in kgav.calibration and is read by the agent and the read-out
    alike; this asserts the document agrees with it."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from kgav.calibration import capabilities
    from repurposing_candidates import tier

    caps = capabilities(RESULTS, cross_tier=tier)
    licensed = sorted(c.name for c in caps.values() if c.shortlist_licensed)
    contradicted = sorted(c.name for c in caps.values() if c.contradicted)
    for doc in ("01-overview.md", "07-results.md"):
        text = " ".join(_text(doc).split())
        assert f"**{len(licensed)} viruses**" in text \
            or f"**{len(licensed)} viruses are licensed" in text, doc
        for name in contradicted:
            assert name in text, f"{doc} does not name {name}"


@needs_artifacts
def test_the_channel_claim_names_the_two_scorers_that_do_clear():
    """THE SENTENCE THIS TEST EXISTS FOR. The first draft said "none has
    an interval excluding chance" and TWO do -- the ensemble and the
    degree baseline. The true claim is about the eight metapath channels,
    and the distinction is the finding rather than a detail."""
    h = json.loads((RESULTS / "hard_negatives.json").read_text())
    cells = h["cross_sectional"]["SARS-CoV-2"]
    chans = [c for c in ("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
             if cells[c]["auc_ci_lo"] > 0.5]
    text = " ".join(_text("01-overview.md").split())
    assert f"{len(chans)} of the 8 metapath channels" in text
    for name in ("COMBINED", "degree"):
        row = cells[name]
        if row["auc_ci_lo"] > 0.5:
            assert f"{row['auc']:.4f}" in text, (
                f"{name} clears and 01-overview does not state its AUC")


@needs_artifacts
def test_the_refusals_are_quoted_rather_than_summarised():
    """A refused protocol is not a protocol that scored 0.5, and the
    reasons are the artifact."""
    doc = _text("06-protocols.md")
    refused = [p for p in sorted(RESULTS.glob("*.json"))
               if "refused_by" in json.loads(p.read_text())]
    assert refused, "no refusal on disk, so this test proves nothing"
    for p in refused:
        assert p.name in doc, f"{p.name} is a refusal and is not named"
        for reason in json.loads(p.read_text())["refused_by"]:
            assert reason[:40] in doc, f"{p.name}'s reason is not quoted"


# --------------------------------------------- what doc2 must not contain
def test_the_signed_registrations_are_linked_and_not_copied():
    """A pre-registration is amended by appending and never rewritten. A
    doc2 copy would be a second, editable version of a signed record --
    which destroys the only property that makes it one."""
    all_text = _all_text()
    assert "docs/prereg_" in all_text, "doc2 does not link the registrations"
    for prereg in sorted((ROOT / "docs").glob("prereg_*.md")):
        body = prereg.read_text()
        # A long verbatim run from a registration is a copy, not a quote.
        for para in body.split("\n\n"):
            flat = " ".join(para.split())
            if len(flat) < 300 or flat.startswith("#"):
                continue
            assert flat not in " ".join(all_text.split()), (
                f"doc2 reproduces a paragraph of {prereg.name} verbatim")


def test_doc2_does_not_claim_a_figure_it_did_not_derive():
    """FOUR-DECIMAL FIGURES ARE AUCS, and every one in doc2 has to come
    out of an artifact. A figure typed into a TEMPLATE renders and
    --check passes, because the file and the generator move together --
    so this is the check --check cannot be."""
    produced = set()
    for p in sorted(RESULTS.glob("*.json")):
        produced |= {f"{m:.4f}" for m in _floats(json.loads(p.read_text()))}
    stray: dict[str, list[str]] = {}
    for p in sorted(DOC2.glob("*.md")):
        for fig in set(re.findall(r"\b0\.\d{4}\b", p.read_text())):
            if fig not in produced:
                stray.setdefault(p.name, []).append(fig)
    assert not stray, (
        f"doc2 prints four-decimal figures no artifact produces: {stray}")


def _floats(obj, out=None):
    """Every float anywhere in a results file, at any depth."""
    out = [] if out is None else out
    if isinstance(obj, dict):
        for value in obj.values():
            _floats(value, out)
    elif isinstance(obj, list):
        for value in obj:
            _floats(value, out)
    elif isinstance(obj, float):
        out.append(obj)
    return out


def test_the_templates_hold_the_prose_and_the_script_holds_the_figures():
    """THE DIVISION THIS FOLDER RESTS ON. A template containing a long
    run of digits is a figure that stopped being derived."""
    for tpl in sorted(TEMPLATES.glob("*.md.in")):
        text = tpl.read_text()
        # Years, small counts and version-like strings are prose. A
        # four-decimal number is an AUC.
        bad = re.findall(r"\b0\.\d{4}\b", text)
        assert not bad, f"{tpl.name} hardcodes {bad}; derive it instead"


@pytest.mark.slow
def test_the_four_scripts_that_refuse_their_own_help_are_still_four():
    """doc2 states a COUNT of scripts that exit non-zero on `--help` on
    purpose, and the count is the one figure in the appendix that cannot
    be derived cheaply enough to render -- 87 subprocesses is not
    something `make doc2` should pay. So it is checked here instead, by
    running the sweep: a fifth such script, or a fourth that started
    raising, makes this fail rather than making the document quietly
    wrong.

    The criterion is the one tests/test_every_script_help_runs.py uses --
    a traceback, not an exit code -- because each of the four is a
    deliberate message.
    """
    scripts = sorted((ROOT / "scripts").glob("*.py"))
    assert len(scripts) > 40
    deliberate = []
    for script in scripts:
        r = subprocess.run([sys.executable, str(script), "--help"],
                           capture_output=True, text=True, timeout=120,
                           cwd=ROOT)
        out = r.stdout + r.stderr
        assert "Traceback (most recent call last)" not in out, (
            f"{script.name} --help raised")
        if r.returncode != 0:
            deliberate.append(script.name)
    words = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
             6: "Six"}.get(len(deliberate), str(len(deliberate)))
    text = _text("13-appendix.md")
    assert f"{words} of them exit non-zero" in text, (
        f"{len(deliberate)} scripts exit non-zero on --help "
        f"({', '.join(deliberate)}) and 13-appendix says otherwise")
