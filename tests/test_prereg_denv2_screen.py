"""docs/prereg_denv2_prospective_screen.md, re-derived.

A PRE-REGISTRATION IS THE ONE DOCUMENT HERE THAT MUST NOT BE EDITED TO
MATCH A LATER RELEASE. Its whole value is that it says what was believed
before the experiment, so `scripts/diagnose_doc_figures.py` is explicit
that prereg figures reading as absent is correct behaviour.

That makes this file narrower than the claim-block tests, and the
distinction is worth stating:

  - WHILE IT IS UNSIGNED -- no arms drawn, no plates read -- every figure
    in it must be the one the code produces today, because it is being
    written now and a typo in a power calculation is how a screen gets
    sized wrong. These tests enforce that.

  - ONCE THE VERDICT IS WRITTEN, the figures become historical and these
    tests must stop re-deriving them. `test_the_registration_is_unsigned`
    is the switch: it fails the moment the verdict line goes, which is
    the prompt to convert this file rather than silently keep checking a
    historical document against a moved release.

THE DESIGN FIGURES ARE THE POINT. The registration fixes 20 compounds per
arm on a stated power, and the arm size is the one number that cannot be
revised afterwards without destroying what the document is for.
"""
from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "prereg_denv2_prospective_screen.md"
SPLIT = ROOT / "data" / "releases" / "v0.1-train2021"

needs_split = pytest.mark.skipif(
    not (SPLIT / "TEST_LABELS.json").exists(),
    reason="no 2021 temporal split on disk")


def _text() -> str:
    return DOC.read_text()


def _flat() -> str:
    return " ".join(_text().split())


def _power_module():
    spec = importlib.util.spec_from_file_location(
        "screen_power_under_test", ROOT / "scripts" / "screen_power.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------ the switch
def test_the_registration_is_unsigned():
    """THE SWITCH DESCRIBED IN THIS FILE'S DOCSTRING.

    While the verdict is unwritten the figures below are checked against
    today's code. When the screen is run and the verdict replaces this
    line, the figures become historical and re-deriving them would be
    the error a prereg exists to prevent -- so this test failing is the
    instruction to convert the rest of this file, not to edit the
    document.
    """
    assert "*Unsigned. Nothing has been drawn and nothing has been "\
           "screened.*" in _text(), (
        "the verdict has been written: STOP re-deriving this document's "
        "figures and convert this test file to assert the recorded "
        "outcome instead")


def test_it_names_no_compound():
    """The list is drawn after signing. A prereg that already contains
    its candidates is a post-hoc rationalisation with a date on it."""
    text = _text()
    # ChEMBL ids and InChIKey-shaped strings are what a drawn list is made
    # of. The virus taxon and the release shas are identifiers, not
    # compounds, and are expected.
    assert not re.search(r"\bCHEMBL\d{3,}\b", text), "names a ChEMBL compound"
    assert not re.search(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b", text), (
        "names an InChIKey")


def test_the_draw_command_is_fixed_and_its_flags_are_named():
    """The command cannot be adjusted after the list is seen, so it is
    written down -- including the flags that do not exist yet, which is
    the point: implementing them cannot change what was promised."""
    flat = _flat()
    assert "scripts/repurposing_candidates.py" in flat
    for flag in ("--exclude-labelled", "--random-arm", "--seed", "--top 20"):
        assert flag in flat, f"the draw command does not fix {flag}"
    assert "do not exist yet" in flat, (
        "the document does not admit that three of its flags are unbuilt; "
        "a reader would assume the command runs today")


# -------------------------------------------------- the design arithmetic
#: The cell this screen is registered for.
DENV2 = "NCBITaxon:11060"
TEMPORAL = ROOT / "data" / "results" / "temporal_structure_scorer.json"


def _denv2_base_rate() -> float:
    """The test-set base rate, from the counts the prospective run recorded.

    DERIVED, NOT TYPED, and the first version of this file typed it. It
    used 0.0651 -- the document's own rounded display -- and computed
    48.5% power where the document says 48.4%. The document was right:
    it was built from the exact 186/2855. A test that rounds its inputs
    and then checks an exact output reports the rounding as a defect in
    the thing it is checking.
    """
    cell = json.loads(TEMPORAL.read_text())["cells"][DENV2]
    pos, neg = cell["n_test_positives"], cell["n_test_negatives"]
    return pos / (pos + neg)


def _stated_precision() -> float:
    """Arm A's assumed hit rate, read from the document's own P1.

    Taken from the document here and checked against the measured run by
    the slow test below, rather than re-measured twice: this test is
    about whether the POWER arithmetic follows from the stated inputs,
    and the slow test is about whether the stated inputs are the
    measured ones.
    """
    m = re.search(r"Arm A hits at \*\*(\d+)% or above\*\*", _text())
    assert m, "P1 does not state arm A's assumed hit rate"
    return int(m.group(1)) / 100


@pytest.mark.skipif(not TEMPORAL.exists(), reason="no prospective results")
def test_the_arm_size_and_its_power_are_what_the_code_computes():
    """20 per arm, at the stated precision and the derived base rate.

    The one number that cannot be revised after the fact, so it is the
    one most worth re-deriving.
    """
    sp = _power_module()
    base, prec = _denv2_base_rate(), _stated_precision()
    flat = _flat()

    for alpha in (0.05, 0.005):
        for k in (10, 16, 20):
            got = sp.power(k, prec, base, alpha)
            assert f"{got:.1%}" in flat, (
                f"power at k={k}, alpha={alpha} is {got:.1%} and the "
                f"document does not print it")


def test_the_rejection_threshold_is_the_one_fisher_gives():
    """'5 or more hits in arm A against 0 in arm B'."""
    sp = _power_module()
    need = next(a for a in range(21)
                if sp.fisher_one_sided(a, 20, 0, 20) <= 0.05)
    assert f"**{need} or more**" in _text(), (
        f"Fisher exact rejects at {need} hits and the document says "
        f"otherwise")


@pytest.mark.skipif(not TEMPORAL.exists(), reason="no prospective results")
def test_the_expected_counts_are_the_products():
    """Arithmetic on the stated precision and the derived base rate, and
    figures of this kind have been wrong in this repository before."""
    flat = _flat()
    base, prec = _denv2_base_rate(), _stated_precision()
    assert f"**{20 * prec:.1f} hits in arm A" in flat
    assert f"{20 * base:.2f} in arm B**" in flat


@pytest.mark.skipif(not TEMPORAL.exists(), reason="no prospective results")
def test_the_base_rate_the_document_prints_is_the_derived_one():
    """Section 2 and section 5 both quote it; neither may round
    differently from the other, and the power figures were computed from
    the exact value."""
    base = _denv2_base_rate()
    flat = _flat()
    assert f"{base:.2%}" in flat, (
        f"the derived base rate is {base:.2%} and the document prints "
        f"something else")


def test_the_family_size_is_one_and_the_document_says_why():
    flat = _flat()
    assert "family size **1**" in flat or "family size 1" in flat
    assert "would cost this one" in flat.lower(), (
        "the document does not give the reason a larger family is refused")
    assert "10 compounds per arm" in flat and "16" in flat


# ------------------------------------------- against the measured figures
@pytest.mark.slow
@needs_split
def test_every_row_of_the_power_table_is_the_scripts():
    """The table that justifies registering one cell and not ten.

    Re-runs scripts/screen_power.py -- one implementation of the power
    calculation, not a second copy here -- and checks every virus row the
    document prints.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "power.json"
        r = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "screen_power.py"),
             "--out", str(out)],
            capture_output=True, text=True, cwd=ROOT, timeout=3600)
        assert r.returncode == 0, r.stderr
        doc = json.loads(out.read_text())

    text = _text()
    block = text.split("every prospective cell", 1)[1].split("\n\n", 2)[1]
    printed = {}
    for line in block.splitlines():
        if not line.startswith("    ") or "virus" in line:
            continue
        name = line[4:18].strip()
        printed[name] = line[18:].split()

    assert set(printed) == {c["name"] for c in doc["cells"]}, (
        f"the table prints {sorted(printed)} and the run holds "
        f"{sorted(c['name'] for c in doc['cells'])}")

    for cell in doc["cells"]:
        row = printed[cell["name"]]
        assert row[0] == f"{cell['base_rate']:.1%}", cell["name"]
        assert row[1] == f"{cell['precision']:.0%}", cell["name"]
        assert row[2] == f"{cell['lift']:.1f}x", cell["name"]
        # The two k columns, each either a number or "out of reach".
        rest = " ".join(row[3:])
        for alpha in doc["alphas"]:
            k = cell["arm_per_alpha"][str(alpha)]
            want = str(k) if k else "out of reach"
            assert want in rest, (
                f"{cell['name']} at alpha {alpha} needs {want} and the "
                f"row reads {rest!r}")


@pytest.mark.slow
@needs_split
def test_the_denv2_precision_and_base_rate_are_the_measured_ones():
    """The two inputs every figure in section 5 rests on."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "hr.json"
        r = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "hit_rate_at_k.py"),
             "--out", str(out), "--k", "10"],
            capture_output=True, text=True, cwd=ROOT, timeout=1800)
        assert r.returncode == 0, r.stderr
        hr = json.loads(out.read_text())

    d2 = next(c for c in hr["cells"] if c["name"] == "DENV-2")
    pk = d2["controls"]["symmetric_k10"]["per_k"]
    prec = (pk["10"] if "10" in pk else pk[10])["precision_expected"]
    flat = _flat()

    assert f"{prec:.0%} against a base rate of {d2['base_rate']:.2%}" in flat, (
        f"section 2 does not state the measured {prec:.0%} / "
        f"{d2['base_rate']:.2%}")
    lift = prec / d2["base_rate"]
    assert f"**{lift:.1f}x**" in flat, f"the lift is {lift:.1f}x"

    # P1 is the CONTROLLED figure and the document must say the
    # uncontrolled one is higher, or P1 reads as the best available guess
    # rather than as a deliberately conservative floor.
    unc = (d2["controls"]["primary"]["per_k"]["10"]
           if "10" in d2["controls"]["primary"]["per_k"]
           else d2["controls"]["primary"]["per_k"][10])["precision_expected"]
    assert f"is {unc:.0%}" in flat, (
        f"the uncontrolled precision is {unc:.0%} and the document does "
        f"not print it beside P1")


# ---------------------------------------------- the protocol's own rules
def test_every_refusal_condition_is_stated_before_the_numbers_exist():
    """Five of them, and the document must not be able to lose one
    quietly: a refusal condition added after a result is not a refusal
    condition."""
    flat = _flat()
    for needle in ("Fewer than 20 eligible orderable compounds",
                   "Any plate control fails",
                   "outside P2's 2%–15%",
                   "The release moves",
                   "Arm membership is disclosed"):
        assert needle in flat, f"refusal condition missing: {needle!r}"


def test_exactly_the_four_named_outcomes_are_declared():
    flat = _flat()
    for outcome in ("**CONFIRMED.**", "**NOT CONFIRMED.**",
                    "**INVERTED.**", "**REFUSED.**"):
        assert outcome in flat, f"{outcome} is not declared"


def test_the_hit_definition_fixes_cytotoxicity_in_advance():
    """A compound that reduces signal by killing cells is a hit by the
    primary definition, and the document says so NOW rather than
    reclassifying it once the plates are read."""
    flat = _flat()
    assert "Cytotoxicity is recorded and is not part of the primary" in flat
    assert "is a hit by this definition" in flat
    assert "Secondary, declared now" in flat


def test_orderability_is_applied_to_the_pool_not_to_arm_a():
    """Otherwise arm A is "the top compounds we could buy" and arm B is
    "anything at all", and the arms differ by purchasability."""
    flat = _flat()
    assert "ORDERABILITY IS A FILTER ON THE POOL, NOT ON ARM A" in flat
    assert "without seeing the ranking" in flat


def test_it_states_what_a_confirmed_result_would_not_license():
    """The section that stops a single-cell result being read as a
    validation of the graph, or of the other fifteen viruses."""
    flat = _flat()
    assert "Would not:" in flat
    assert "other fifteen viruses" in flat
    assert "ranking uses no graph edges" in flat
