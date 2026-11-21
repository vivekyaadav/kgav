"""The headline claims in doc2, re-derived from the artifacts.

WHAT THIS REPLACES. `tests/test_claim_block.py` held four fixed-width
blocks in `docs/evaluation.md` to the artifacts, cell by cell, and it was
the strongest check in the repository: it is the reason the claim
"two scorers clear chance, both cross-sectionally, both on SARS-CoV-2"
was caught describing no run that had ever happened. When `docs/` was
retired on 2026-11-21 those blocks became generated tables, parsed at
fixed column offsets that a generated document has no reason to keep.

So the derivations moved here and the parsing changed. Three of the four
blocks already have a second derivation elsewhere —
`tests/test_auc_decomposition.py` for the channels,
`tests/test_temporal_table.py` for the prospective cells,
`tests/test_evaluation_table.py` for the chance-clearing census — and
this file covers the fourth, the structure scorer against its own degree
comparator, plus the two claims that sit beside it.

WHY A SECOND DERIVATION AT ALL, when `make doc2-check` already fails on a
document that does not match `write_doc2.derive()`. Because `--check`
compares the document against the generator, and a generator with a bug
produces a document that matches it. Everything here derives the same
figures a different way and requires the two to agree.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

RESULTS = ROOT / "data" / "results"
STRENGTH = RESULTS / "structure_scorer_strength_matched.json"
SYMMETRIC = RESULTS / "structure_scorer_symmetric.json"
CLINICAL = RESULTS / "clinical_labels.json"
RESULTS_DOC = ROOT / "doc2" / "08-results.md"
LIMITS_DOC = ROOT / "doc2" / "09-limits.md"
OVERVIEW_DOC = ROOT / "doc2" / "01-overview.md"
MARKER = "<!-- TABLE: scorer_vs_degree"

needs_artifacts = pytest.mark.skipif(
    not STRENGTH.exists(), reason="the structure scorer run is not on disk")


def _names() -> dict[str, str]:
    """CURIE -> display name, from the register.

    Not a second map: a hard-coded one here would be a fifth place the
    virus names live, and the register has already outgrown three.
    """
    from kgav import virus_register as vr
    return {f"NCBITaxon:{e['taxon']}": e["name"] for e in vr.load()["viruses"]}


SECTIONS = (("primary_chembl", "ChEMBL"),
            ("secondary_ncats", "plate CPE"),
            ("secondary_ncats_entry", "plate entry"))


def derive(path: Path = STRENGTH,
           field: str = "structure_max") -> dict[tuple[str, str], dict]:
    """(virus, population) -> the cells the run implies, derived here.

    `readable` is the run's own flag; a cell it marks unreadable is not
    an argument against the scorer, it is the absence of one, and
    counting it either way would be a different claim.
    """
    doc = json.loads(path.read_text())
    label = _names()
    out = {}
    for section, population in SECTIONS:
        for curie, cell in (doc.get(section) or {}).items():
            if not cell.get("readable"):
                continue
            s = cell.get(field) or {}
            if s.get("auc") is None:
                continue
            g = cell.get("degree") or {}
            s_clears = (bool(s.get("clears_corrected_alpha"))
                        and s.get("direction") != "below")
            g_clears = bool(g.get("clears_corrected_alpha"))
            out[(label.get(curie, curie), population)] = {
                "auc": s["auc"], "clears": s_clears,
                "degree_auc": g.get("auc"), "degree_clears": g_clears,
                "group": ("clears, degree does not" if s_clears and not g_clears
                          else "clears, degree clears too" if s_clears
                          else "does not clear"),
            }
    return out


def parse() -> dict[tuple[str, str], dict]:
    """The table as doc2 prints it."""
    text = RESULTS_DOC.read_text()
    assert MARKER in text, "the scorer/degree table lost its marker"
    lines = [ln.strip() for ln in text.split(MARKER, 1)[1].splitlines()]
    start = next(i for i, ln in enumerate(lines) if ln.startswith("| virus |"))
    rows = {}
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        c = [x.strip() for x in line.strip("|").split("|")]
        assert len(c) == 7, f"expected 7 columns: {line}"
        key = (c[0], c[1])
        assert key not in rows, f"{key} appears twice"
        rows[key] = {
            "auc": float(c[2]),
            "clears": c[3].strip("*") == "yes",
            "degree_auc": None if c[4] == "—" else float(c[4]),
            "degree_clears": c[5] == "yes",
            "group": c[6],
        }
    assert rows, "found the marker and parsed no rows"
    return rows


@pytest.fixture(scope="module")
def tables():
    if not STRENGTH.exists():
        pytest.skip("the structure scorer run is not on disk")
    return parse(), derive()


# ------------------------------------------------- the table, cell by cell
def test_the_table_holds_every_readable_cell_and_no_other(tables):
    """A TABLE OF THE FIRST GROUP ALONE IS THE FAILURE THIS PROJECT IS
    BUILT AGAINST, and so is one that is almost complete: the first
    draft of the predecessor printed three of four groups while its own
    comment claimed the table was the whole family."""
    doc, got = tables
    assert set(doc) == set(got), (
        f"only in the document: {sorted(set(doc) - set(got))}; "
        f"only in the run: {sorted(set(got) - set(doc))}")


@pytest.mark.parametrize("field", ["auc", "clears", "degree_auc",
                                   "degree_clears", "group"])
def test_every_cell_matches_the_run(tables, field):
    doc, got = tables
    wrong = []
    for key in sorted(set(doc) & set(got)):
        a, b = doc[key][field], got[key][field]
        if field.endswith("auc"):
            if a is None or b is None:
                same = a is None and b is None
            else:
                same = f"{a:.4f}" == f"{b:.4f}"
        else:
            same = a == b
        if not same:
            wrong.append(f"{key}: document {a!r}, run {b!r}")
    assert not wrong, "\n  " + "\n  ".join(wrong)


def test_the_group_counts_beside_the_table_are_the_tables_own(tables):
    """The three counts in the sentence under it, derived. They are what
    a reader takes away, and the claim is the FIRST of them: the scorer
    clears where degree does not."""
    _doc, got = tables
    counts: dict[str, int] = {}
    for cell in got.values():
        counts[cell["group"]] = counts.get(cell["group"], 0) + 1
    flat = " ".join(RESULTS_DOC.read_text().split())
    assert f"**{counts.get('clears, degree does not', 0)} clear above a " \
        f"degree baseline that does not clear.**" in flat, (
            f"{counts} and the document states different counts")
    assert f"{counts.get('clears, degree clears too', 0)} clear while " \
        f"their baseline clears too" in flat
    assert f"{counts.get('does not clear', 0)} do not clear" in flat
    assert str(len(got)) in flat, (
        f"{len(got)} cells are readable and the document does not say so")


def test_the_two_controls_agree_on_the_cells_the_claim_rests_on():
    """THE RESULT'S MAIN DEFENCE, asserted rather than written down —
    and it is NARROWER than the sentence first written beside it.

    This test was written on 2026-11-21 to check "the membership of the
    clearing set is identical under both registered controls", stated
    over the whole family, and it failed the same day: the two controls
    DISAGREE on the MERS-CoV pseudotyped-entry cell, where symmetric
    k=10 clears and strength-matched does not. They agree on every
    ChEMBL cell and on the cytopathic plate cell.

    The defensible claim is about the group the result rests on — the
    cells that clear above a degree baseline that does NOT clear — and
    that group is identical under both. In the one disagreeing cell the
    baseline clears and sits above the structure statistic under both
    readings, so it is outside the claim in both; the documents say so,
    with the cell named.
    """
    if not SYMMETRIC.exists():
        pytest.skip("the symmetric control is not on disk")
    st = derive(STRENGTH)
    sy = derive(SYMMETRIC, "structure_max_k10")
    g1 = "clears, degree does not"
    st_claim = {k for k, c in st.items() if c["group"] == g1}
    sy_claim = {k for k, c in sy.items() if c["group"] == g1}
    assert st_claim == sy_claim, (
        f"the two controls now disagree about the group the claim rests "
        f"on: only strength-matched {sorted(st_claim - sy_claim)}, only "
        f"symmetric {sorted(sy_claim - st_claim)}. The claim has to be "
        f"rewritten, not this test relaxed")

    # AND WHERE THEY DISAGREE AT ALL, the documents must name it. A
    # defence stated more widely than it holds is the failure this test
    # found in its own first hour.
    differ = sorted({k for k in set(st) | set(sy)
                     if st.get(k, {}).get("clears")
                     != sy.get(k, {}).get("clears")})
    flat = " ".join(RESULTS_DOC.read_text().split())
    assert "identical under both" in flat, (
        "the documents no longer state the property this test checks")
    for virus, population in differ:
        assert f"{virus} ({population})" in flat, (
            f"the controls disagree on {virus} ({population}) and the "
            f"document does not name it")
        cell = st.get((virus, population)) or sy[(virus, population)]
        assert cell["degree_clears"], (
            f"the controls disagree on {virus} ({population}) and its "
            f"degree baseline does NOT clear, so that cell is inside the "
            f"claim and the disagreement is material")


# ----------------------------------------------- the clinical contradiction
@pytest.mark.skipif(not CLINICAL.exists(), reason="no clinical run")
def test_the_clinical_axis_is_reported_as_contradicting():
    """The one independent label axis anti-correlates with the chemistry,
    and the direction is the finding. A document that reported it as
    merely "not significant" would be describing a different result."""
    c = json.loads(CLINICAL.read_text())
    cell = c["cells"][0] if isinstance(c.get("cells"), list) else c
    struct = (cell.get("structure_max") or {}) if "structure_max" in cell \
        else (c.get("structure_max") or {})
    degree = (cell.get("degree") or {}) if "degree" in cell \
        else (c.get("degree") or {})
    if not struct or struct.get("auc") is None:
        pytest.skip("the clinical artifact has a shape this test cannot read")
    flat = " ".join(RESULTS_DOC.read_text().split())
    assert f"{struct['auc']:.4f}" in flat, (
        f"the clinical structure AUC {struct['auc']:.4f} is not in doc2")
    if degree.get("auc") is not None:
        assert f"{degree['auc']:.4f}" in flat
        assert degree["auc"] > struct["auc"], (
            "degree no longer beats the structure scorer on the clinical "
            "axis; the paragraph's argument has changed")
    assert struct.get("direction") == "below", (
        "the clinical cell no longer clears below chance; doc2 calls it a "
        "contradiction and must be rewritten rather than reinterpreted")


# ------------------------------------------------- the hit-rate limitation
def test_the_hit_rate_qualification_says_what_it_is():
    """A POST-HOC STATISTIC PRESENTED BESIDE REGISTERED ONES ACQUIRES
    THEIR STANDING BY ADJACENCY, so it has to say what it is."""
    flat = " ".join(RESULTS_DOC.read_text().split())
    assert "post-hoc" in flat.lower(), (
        "no section of 08-results labels its post-hoc measurements")
    assert "registered" in flat


def test_the_hit_rate_limit_keeps_no_copy_of_the_figures():
    """A LIMIT RESTATED WITH ITS OWN COPY OF THE NUMBERS IS A SECOND
    THING TO GO STALE — the failure the retired scope statement had
    twice.

    So doc2 states the hit-rate finding and deliberately prints no
    percentages for it: `make hit-rate` measures the head of the ranking
    and is the figures' one home. This asserts the limit is there, names
    the mechanism, and carries no ladder of its own.
    """
    text = LIMITS_DOC.read_text()
    assert "An AUC above chance is not an expected hit rate" in text, (
        "the limits document does not carry the hit-rate limit")
    section = text.split("An AUC above chance is not an expected hit rate",
                         1)[1].split("\n## ", 1)[0]
    flat = " ".join(section.split())
    assert "no actives at all" in flat, "it does not state the finding"
    assert "make hit-rate" in flat, (
        "it does not point at the measurement that owns the figures")
    # Whole figures, not substrings: an unrelated "20.0%" elsewhere was
    # how the first version of this check in test_claim_block reported a
    # duplication that was not one, and a test that cries wolf gets
    # relaxed until it reports nothing.
    stray = sorted(set(re.findall(r"\b\d{1,3}\.\d%", flat)))
    assert not stray, (
        f"the hit-rate limit prints {stray}; those figures have one home "
        f"and this is a second place for them to drift")
