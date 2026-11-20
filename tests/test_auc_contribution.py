"""scripts/auc_contribution.py reports results; it must not contain them.

THE RULE THIS ENFORCES. Every figure in this project is derived from
data/results, never typed into the thing that displays it. A
decomposition is the most tempting place to break that rule -- the
numbers are the whole output, and typing one is faster than reading it --
and the most damaging, because a share looks like arithmetic whoever
wrote it.

So: the script may contain thresholds it declares (0.5 is chance, 1000 is
a display cap) and no AUC-shaped literal at all. Then the figures it
prints are re-derived here from the same files, independently, and
compared.
"""
from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "auc_contribution.py"
RESULTS = ROOT / "data" / "results"

needs_results = pytest.mark.skipif(
    not (RESULTS / "structure_scorer_strength_matched.json").exists(),
    reason="no results on disk")


def _module():
    spec = importlib.util.spec_from_file_location(
        "auc_contribution_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_the_script_types_no_auc():
    """A literal with three or more decimal places is an AUC or an alpha.

    CHANCE is 0.5 and is declared once as a named constant, which is the
    definition of the statistic and not a measurement.
    """
    src = SCRIPT.read_text()
    bad = [m.group(0) for m in re.finditer(r"\b\d+\.\d{3,}\b", src)]
    assert not bad, (
        f"{SCRIPT.name} contains the literal(s) {bad}. Every figure it "
        f"prints must be read from data/results.")


def test_chance_is_declared_once_and_is_one_half():
    """IN CODE. The first version of this counted 0.5 as text and found
    seven, every one of them prose explaining what chance is -- "an AUC
    of 0.5 is chance", "the distance from 0.5". A test that forbids the
    documentation from naming the constant it documents is a test that
    gets satisfied by deleting the documentation.
    """
    import ast

    mod = _module()
    assert mod.CHANCE == 0.5

    halves = [n for n in ast.walk(ast.parse(SCRIPT.read_text()))
              if isinstance(n, ast.Constant) and n.value == 0.5]
    assert len(halves) == 1, (
        f"0.5 appears {len(halves)} times in code; it is CHANCE and has "
        f"one home, so a second spelling can drift from it")


@needs_results
def test_the_report_runs_and_names_its_release():
    out = subprocess.run([sys.executable, str(SCRIPT)],
                         capture_output=True, text=True, cwd=ROOT, timeout=300)
    assert out.returncode == 0, out.stderr
    from kgav.provenance import release_fingerprint
    from kgav.pipeline import release_dir
    fp = release_fingerprint(release_dir(ROOT))
    assert fp["nodes_sha"] in out.stdout
    assert fp["edges_sha"] in out.stdout


@needs_results
def test_the_prospective_figures_are_the_results_files_figures():
    """Re-derived here from the file, independently of the script."""
    doc = json.loads((RESULTS / "temporal_structure_scorer.json").read_text())
    out = subprocess.run([sys.executable, str(SCRIPT)],
                         capture_output=True, text=True, cwd=ROOT, timeout=300)
    assert out.returncode == 0, out.stderr

    clears = set(doc["clears"])
    shown = 0
    for taxon, cell in doc["cells"].items():
        auc = (cell.get("primary") or {}).get("auc")
        if not isinstance(auc, (int, float)):
            continue
        shown += 1
        assert f"{auc:.4f}" in out.stdout, (
            f"{cell.get('name', taxon)} AUC {auc:.4f} is in the results file "
            f"and not in the report")
    assert shown == doc["n_evaluable"], (
        f"the file says {doc['n_evaluable']} cells are evaluable and "
        f"{shown} carry a primary AUC")
    assert f"{len(clears)} clear above chance" in out.stdout


@needs_results
def test_the_degree_baseline_is_not_counted_as_an_approach():
    """THE DEFECT THE FIRST VERSION OF THIS REPORT HAD.

    With the degree baseline in the denominator, the clinical-labels
    degree row took 36% of "the contribution" -- a control credited with
    carrying the result it exists to be measured against. The baselines
    are reported, below a rule, with no share.
    """
    out = subprocess.run([sys.executable, str(SCRIPT)],
                         capture_output=True, text=True, cwd=ROOT,
                         timeout=300).stdout
    head, _, tail = out.partition("BASELINES")
    assert tail, "the report no longer separates baselines from approaches"
    # The approach table ends before the baselines; no degree row above it.
    approaches = head.split("1. WHAT EACH APPROACH CONTRIBUTES")[-1]
    offenders = [ln for ln in approaches.splitlines()
                 if "degree baseline" in ln and "%" in ln]
    assert not offenders, (
        f"a degree baseline row carries a share: {offenders}")


@needs_results
def test_a_refused_cell_is_never_read_as_an_auc_of_one_half():
    """`measurable: false` means the protocol declined the cell.

    Its arithmetic AUC is 0.5 and that is not a measurement. Reading it
    as one would enter a non-result into the averages as 'exactly
    chance', which is the kind of figure that then gets quoted.
    """
    mod = _module()
    refused = {"structure_max": {"auc": 0.5, "measurable": False,
                                 "refusal": "nothing was scored in this cell"}}
    assert mod.auc(refused, "structure_max") is None

    real = {"structure_max": {"auc": 0.5, "measurable": True}}
    assert mod.auc(real, "structure_max") == 0.5, (
        "a genuine 0.5 must still read as 0.5; the guard is about the "
        "refusal flag, not about the value")


@needs_results
def test_every_cell_the_report_cannot_compare_is_named():
    """A CHECK THAT CANNOT FAIL IS WORSE THAN NO CHECK.

    The first version looked for `measurable: false` in the
    strength-matched file, found none, and skipped -- while the report
    was in fact leaving two cells out of its table. The cells it cannot
    compare are not all refusals: HCoV-NL63 has no positives, so it
    carries no structure channel at all, and a test that only knows
    about the refusal flag cannot see it.

    The condition is therefore asked with the script's own reader: a
    cell whose structure or degree AUC comes back None is a cell the
    table does not have a row for, and every one of them must be named
    under "not compared".
    """
    mod = _module()
    sm = mod.cells(json.loads(
        (RESULTS / "structure_scorer_strength_matched.json").read_text()))
    deg = mod.cells(json.loads((RESULTS / mod.DEGREE[0]).read_text()))

    uncomparable = [
        mod.name_of(label) for label in sm
        if mod.auc(sm[label], "structure_max") is None
        or mod.auc(deg.get(label, {}), mod.DEGREE[1]) is None]
    assert uncomparable, (
        "no cell is uncomparable today, so this test asserts nothing; "
        "drop it or point it at a run where one is")

    out = subprocess.run([sys.executable, str(SCRIPT)],
                         capture_output=True, text=True, cwd=ROOT,
                         timeout=300).stdout
    _, _, tail = out.partition("not compared")
    assert tail, (
        f"{len(uncomparable)} cell(s) have no comparable row and the report "
        f"does not say so: {uncomparable}")
    for nice in uncomparable:
        assert nice in tail.split("\n\n")[0], (
            f"{nice} is left out of the table and not named: {uncomparable}")
