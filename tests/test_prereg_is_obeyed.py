"""A registered decision must be what the code actually does.

WHY THIS FILE EXISTS. docs/prereg_pathway_overlap.md declares, under
"Fixed in advance" and flagged **DECISION**:

    Permutations: 100,000. The hard floor is one attainable value below
    alpha, which needs N+1 > 1/alpha = 880. Ten steps needs N >= 8,799.

`pathway_landscape.DEFAULT_PERMUTATIONS` was 2000, and
data/results/pathway_landscape.json records `permutations: 2000`. So the
registration fixed a number and the run used a different one, for as long
as both have existed.

THIS IS WORSE THAN AN UNDECLARED PARAMETER. An undeclared one is visibly
free. A declared one that the code ignores reads, to anyone checking, as a
pre-registration that was honoured -- the document says 100,000 and nothing
pointed at the 2,000 beside it. The project's standing rule is that a
correction family and a permutation count are fixed BEFORE a p-value is in
hand; obeying them afterwards is the other half of that rule and had no
check.

THE FIGURE IS READ FROM THE DOCUMENT, not repeated here. A test that typed
100000 on both sides would pass if someone changed the prereg, which is the
"check written with the literals it is checking" defect this repository
keeps finding.
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

PREREG = ROOT / "docs/prereg_pathway_overlap.md"
RESULTS = ROOT / "data/results/pathway_landscape.json"


def declared_permutations() -> int:
    """The number the pre-registration fixes, from the document."""
    text = PREREG.read_text()
    m = re.search(r"\*\*Permutations:\s*([\d,]+)\.\*\*\s*\*\*DECISION\*\*", text)
    assert m, ("the prereg no longer states its permutation count in the form "
               "this test reads; update both together")
    return int(m.group(1).replace(",", ""))


def test_the_code_default_is_the_registered_count():
    """THE TEST THAT WOULD HAVE CAUGHT IT. Proven to fail against
    DEFAULT_PERMUTATIONS = 2000, which is what it was."""
    import pathway_landscape as PL
    declared = declared_permutations()
    assert PL.DEFAULT_PERMUTATIONS == declared, (
        f"the prereg fixes {declared:,} permutations and the code defaults to "
        f"{PL.DEFAULT_PERMUTATIONS:,}. A registered decision the code ignores "
        f"reads as a registration that was honoured")


def test_the_registered_count_clears_its_own_stated_floor():
    """The prereg derives a floor of 8,799 for ten attainable steps below
    alpha = 0.001136. The declared count must clear it, or the document
    contradicts its own arithmetic."""
    declared = declared_permutations()
    alpha = 0.05 / 44           # the prereg's own family: 22 rows x 2 tails
    one_step = 1 / alpha - 1    # N+1 > 1/alpha
    assert declared > one_step, (declared, one_step)
    assert declared >= 8799, (
        f"{declared:,} gives fewer than ten attainable values below "
        f"{alpha:.6f}; the prereg's floor is 8,799")


@pytest.mark.skipif(not RESULTS.exists(), reason="no landscape result on disk")
def test_the_result_on_disk_was_computed_at_the_registered_count():
    """And the artifact must agree. This is the assertion that fails until
    the landscape is re-run; it is the point of the exercise, not a
    nuisance, because a result computed at 2,000 cannot support the alpha
    the same document fixes.
    """
    got = json.loads(RESULTS.read_text()).get("permutations")
    declared = declared_permutations()
    assert got == declared, (
        f"pathway_landscape.json was computed at {got:,} permutations and the "
        f"prereg fixes {declared:,}. Re-run `make pathway-landscape`")
