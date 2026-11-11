"""The assembled release's path is written in exactly one place.

FIX-LIST L14. The literal `data/releases/v0.1` was a `default=` in 37
argument definitions across 28 files. The row that filed it said "ten
places", and `pipeline.release_dir`'s own docstring said "nine other
scripts" -- both counted correctly when written, and both were stale,
because every script added afterwards copied the literal from the one
beside it. A defect that grows while its own description apologises for
being small is the argument for converting rather than documenting.

AND 37 WAS NOT THE COUNT EITHER. The migration was scoped to `default=`
expressions, because that is how the row described the defect. This test
was written afterwards, looked for path CONSTRUCTION anywhere, and found
two more: `check_results.py`'s `--derived` fallback and the ORCS edge path
in `ingest_host_targets.py`. **39 sites, in 28 files.** Three successive
counts of one defect -- ten, 37, 39 -- each taken from a description of it
rather than from the tree. The number in this docstring is the only one a
test re-derives.

WHAT THIS TEST DOES NOT CLAIM. Nothing here was broken: the scripts READ
where `release_dir` says WRITE, so they could not disagree with it
silently the way `assemble` and `rebuild` could before L13. The cost was a
release id that could not be changed in one edit. This test keeps that
property rather than asserting a bug is gone.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: The one module allowed to spell the path. It holds `release_dir`,
#: `layer_dir` and `train_dir`, and every other caller goes through them.
OWNER = ROOT / "src/kgav/pipeline.py"

#: The release path being CONSTRUCTED: joined onto a root with `/`, or
#: wrapped in `Path(...)`. Covers the bare directory, a layer suffix and a
#: training cut.
#:
#: MENTIONING THE PATH IS ALLOWED AND MUST BE. `ingest_folds.py` refuses
#: without a spine and tells the operator to "build data/releases/v0.1-spine
#: first"; that sentence is for a human, and a test that forbade it would be
#: asking the error message to be less useful than the error. What this
#: pattern catches is a second place that COMPUTES the path, which is the
#: thing that can drift from pipeline.py without anyone noticing.
LITERAL = re.compile(
    r'(?:[\w.\]\)]\s*/\s*|Path\(\s*)"data/releases/v0\.1[^"]*"'
    #: AND THE COMPONENT-WISE FORM, which escaped both the fix-list L14
    #: migration and the first version of this test. ingest_gtop.py had
    #: `default=ROOT / "data" / "releases" / "v0.1"` -- the same path, built
    #: from four strings instead of one, so a pattern anchored on the joined
    #: literal could not see it. Found on 2026-11-11 by reading that file for
    #: an unrelated reason, which is the 40th site and the third wrong count
    #: of this defect. A test that can only see one spelling of the thing it
    #: forbids is a test that reports the spellings it knows.
    r'|"data"\s*/\s*"releases"\s*/\s*"v0\.1'
    r'|"releases"\s*/\s*"v0\.1')


def _code_lines(path: Path) -> list[tuple[int, str]]:
    """Lines that are not inside a docstring and not a comment.

    Prose MAY name the path -- a help string, a module docstring and a
    usage example all legitimately print it, and a test that forbade that
    would be asking documentation to lie. Only executable spellings count.
    """
    src = path.read_text()
    out, in_doc, delim = [], False, ""
    for i, line in enumerate(src.splitlines(), 1):
        s = line.strip()
        if in_doc:
            if delim in s:
                in_doc = False
            continue
        for d in ('"""', "'''"):
            if s.startswith(d) or re.match(rf'^[a-zA-Z_]*{d}', s):
                if s.count(d) == 1:
                    in_doc, delim = True, d
                break
        if in_doc or s.startswith("#"):
            continue
        out.append((i, line))
    return out


def test_only_pipeline_spells_the_release_path():
    """Proven to fail against the tree of 2026-11-11: run against the
    pre-migration files it reports **39 site(s) spell the release path**."""
    offenders: list[str] = []
    for path in sorted(list((ROOT / "scripts").glob("*.py"))
                       + list((ROOT / "src/kgav").glob("*.py"))):
        if path.resolve() == OWNER.resolve():
            continue
        for lineno, line in _code_lines(path):
            if "help=" in line or "description=" in line:
                continue
            if LITERAL.search(line):
                offenders.append(f"{path.relative_to(ROOT)}:{lineno}: {line.strip()}")
    assert not offenders, (
        f"{len(offenders)} site(s) spell the release path instead of calling "
        f"kgav.pipeline.release_dir / layer_dir / train_dir:\n  "
        + "\n  ".join(offenders))


def test_the_three_helpers_agree_with_each_other():
    """They must share a parent, or `release_dir` and `layer_dir` would
    point at different release families and the rebuild would assemble one
    while fingerprinting the other -- the failure mode L13 fixed."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.pipeline import layer_dir, release_dir, train_dir
    r = Path("/x")
    assert release_dir(r) == r / "data/releases/v0.1"
    assert layer_dir(r, "spine") == r / "data/releases/v0.1-spine"
    assert train_dir(r) == r / "data/releases/v0.1-train2021"
    assert train_dir(r, 2019) == r / "data/releases/v0.1-train2019"
    assert {layer_dir(r, "spine").parent, train_dir(r).parent} == {release_dir(r).parent}
