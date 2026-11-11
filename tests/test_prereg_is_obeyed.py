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


# ---------------------------------------------------------------------------
# PINNED BY THE SIGNATURE, 2026-11-11. Until the document was signed, only
# the permutation count was checked here, and tests/test_landscape.py said
# in a comment why the rest were not: "a test asserting a proposal would
# make the signature meaningless". That reasoning held for a draft and
# inverts for a registration. The other three registered values follow.


def _normalised() -> str:
    return " ".join(PREREG.read_text().split())


def declared_family_size() -> int:
    """The correction family the document fixes, read from its own alpha."""
    m = re.search(r"alpha\s*=\s*0\.05\s*/\s*(\d+)\s*=", _normalised())
    assert m, ("the prereg no longer states its family size as 0.05/N; "
               "update both together")
    return int(m.group(1))


def declared_drugged_operand() -> int:
    """The first operand the document fixes: the drugged proteome under
    `floor`. Read from the DECISION paragraph, not from the signature's
    restatement of it, because the DECISION paragraph is the registration."""
    m = re.search(r"first operand -- is ([\d,]+)\s*human proteins under "
                  r"`floor`", _normalised())
    assert m, ("the prereg no longer states the drugged-proteome operand in "
               "the form this test reads")
    return int(m.group(1).replace(",", ""))


@pytest.mark.skipif(not RESULTS.exists(), reason="no landscape result on disk")
def test_the_result_used_the_registered_correction_family():
    """`n_tests` is what alpha is divided by, and the document fixes it at
    22 rows x 2 tails. A run that performed a different number of tests is
    corrected against the wrong threshold, and the file records the count
    precisely so this is checkable."""
    got = json.loads(RESULTS.read_text()).get("n_tests")
    declared = declared_family_size()
    assert got == declared, (
        f"the run performed {got} tests and the prereg fixes alpha at "
        f"0.05/{declared}; the Bonferroni threshold on disk is not the "
        f"registered one")


@pytest.mark.slow
def test_the_registered_operand_is_the_one_the_release_gives():
    """THE PARAMETER THE ARTIFACT DOES NOT RECORD.

    The signature names this as its own open item: `pathway_landscape.json`
    stores the permutation count, the seed scheme, the degree basis and the
    family size, and NOT the size of the drugged proteome -- the one operand
    the document says the verdict depends on (0 of 22 rows clearing
    Bonferroni under `floor`, 3 under `relaxed`). So which admission policy
    produced the figures on disk cannot be read off them.

    It can be re-derived, and that is what this does: the drugged proteome
    is a function of the release, so the registered 1,946 is checked against
    the release the registration names. A disagreement means either the
    release moved under the registration or the run was made under a policy
    the document does not register -- and the two are told apart by the
    release identification test in tests/test_landscape.py.

    Marked slow: it loads the assembled release, which is minutes not
    seconds.
    """
    import pathway_landscape as PL
    from kgav.baselines import Graph
    from kgav.pipeline import release_dir
    from kgav.schema import load_schema

    rd = release_dir(ROOT)
    if not (rd / "nodes.jsonl").exists():
        pytest.skip("no assembled release on disk")
    g = Graph.load(rd, symmetry=load_schema().symmetry())
    drugged = len(PL.drug_targets(g))
    declared = declared_drugged_operand()
    assert drugged == declared, (
        f"the registration fixes the drugged proteome at {declared:,} "
        f"proteins under `floor` and the release gives {drugged:,}. The "
        f"registered operand is not the one this release produces")


# ---------------------------------------------------------------------------
# THE RELEASE-IDENTIFICATION COUNT, DERIVED. The document maintained this
# number by hand and it drifted: it read "eight times -- 2026-11-08 and
# seven times on 2026-11-09", which matched neither the total (eleven) nor
# either date (three and five). A prose count nobody re-derives is the same
# shape as the figures this project moved into <!-- FIGURE: --> markers, so
# it is derived here from git and compared against what the document says.

_NUMBER_WORDS = {"once": 1, "one": 1, "twice": 2, "two": 2, "three": 3,
                 "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8,
                 "nine": 9, "ten": 10, "eleven": 11, "twelve": 12,
                 "thirteen": 13, "fourteen": 14, "fifteen": 15}


def _identification_commits() -> list[tuple[str, str]]:
    """(sha, date) for every commit that rewrote the release identification.

    The pickaxe, not `git log -L`: the line-number form needs the number the
    line had at the revision being read, and the `/regex/` form re-resolves
    as it walks and drops d43c316. Both failure modes were observed on this
    file; the pickaxe was checked to return the same twelve commits as the
    correct line-number invocation.
    """
    import subprocess
    out = subprocess.run(
        ["git", "log", "--format=%h %ad", "--date=short",
         "-G", "nodes_sha.*edges_sha", "--", "docs/prereg_pathway_overlap.md"],
        cwd=ROOT, capture_output=True, text=True, timeout=60)
    if out.returncode != 0:
        pytest.skip(f"git unavailable: {out.stderr.strip()[:80]}")
    rows = [ln.split() for ln in out.stdout.splitlines() if ln.strip()]
    return [(r[0], r[1]) for r in rows]


def test_the_identification_update_count_is_what_git_says():
    """The total, derived. The oldest of the commits is the drafting, so the
    number of UPDATES is one less than the number of commits."""
    commits = _identification_commits()
    if len(commits) < 2:
        pytest.skip("history too shallow to count identification updates")
    updates = len(commits) - 1

    text = _normalised()
    m = re.search(r"release identification has now been updated (\w+) times",
                  text)
    assert m, ("the document no longer states how many times the release "
               "identification was updated, in the form this test reads")
    claimed = _NUMBER_WORDS.get(m.group(1).lower())
    assert claimed is not None, f"unrecognised number word {m.group(1)!r}"
    assert claimed == updates, (
        f"git records {updates} update(s) to the release identification and "
        f"the document claims {claimed}. The old hand-maintained count said "
        f"eight when git said eleven")


def test_the_identification_breakdown_by_date_is_what_git_says():
    """AND THE BREAKDOWN, which is the half that was wrong while the total
    was merely stale: the document put one update on 2026-11-08 and seven on
    2026-11-09, where git puts three and five."""
    from collections import Counter
    commits = _identification_commits()
    if len(commits) < 2:
        pytest.skip("history too shallow to count identification updates")
    # the oldest is the drafting, not an update
    by_date = Counter(d for _sha, d in commits[:-1])

    # SCOPED TO THE CURRENT CLAIM, and the first version of this test was
    # not. The document quotes its own superseded breakdown -- "eight times
    # -- 2026-11-08 and seven times on 2026-11-09" -- as the thing it got
    # wrong, and a document-wide regex read that quotation as an assertion
    # and failed on it. A superseded figure recorded on purpose is not a
    # claim about the graph; scripts/diagnose_doc_figures.py says the same
    # about evaluation.md at length. So only the sentence that makes the
    # claim is parsed.
    text = _normalised()
    span = re.search(r"release identification has now been updated .*?"
                     r"-- and the analysis", text)
    assert span, ("the document no longer states the breakdown in the "
                  "sentence this test reads")
    stated = re.findall(r"(\w+) (?:times? )?on (\d{4}-\d\d-\d\d)",
                        span.group(0))
    assert stated, "the document no longer breaks the count down by date"
    for word, day in stated:
        n = _NUMBER_WORDS.get(word.lower())
        if n is None:          # a date mentioned in prose, not a count
            continue
        assert by_date.get(day) == n, (
            f"the document puts {n} update(s) on {day} and git puts "
            f"{by_date.get(day, 0)}")
    # and nothing git knows about may be left out of the breakdown
    named = {day for word, day in stated if word.lower() in _NUMBER_WORDS}
    missing = {d: c for d, c in by_date.items() if d not in named}
    assert not missing, (
        f"the breakdown omits {missing}; a date with updates that the "
        f"document does not mention is how the count drifted before")
