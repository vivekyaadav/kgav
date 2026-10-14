"""The results directory is the evidence behind docs/evaluation.md, so what
is IN the repository has to be the set that gets audited.

data/results is tracked deliberately -- data/raw and data/releases are
gitignored, this is not -- and its last commit says so: "Diagnostics and
results behind the six causes in docs/evaluation.md".

`make check-results` audits the WORKING COPY. It passed on the workstation
with fourteen files while eleven were in git, so three results existed on
one machine and nowhere else -- among them the Zika baselines, the only
per-virus evidence file for the family the project had just added. A reader
of the repository could not check the claim, and nothing said a file was
missing rather than the virus having no result.

Same shape as every other audit gap here: correct about the thing it looks
at, silent about the thing beside it.
"""
import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "results"


def _tracked() -> set[str] | None:
    """Filenames COMMITTED under data/results, or None if this is not a repo.

    ls-tree HEAD, not ls-files. The first version used ls-files, which
    reports the INDEX -- so `git add data/results` alone satisfied it while
    nothing was published, and that is precisely what happened: the commit
    that was meant to follow aborted on an empty message, the staged files
    made this test pass, and only the provenance test (which already read
    HEAD) still said the repository was wrong.

    A test about what a reader of the repository gets must read what a
    reader of the repository gets.
    """
    try:
        out = subprocess.run(
            ("git", "ls-tree", "-r", "--name-only", "HEAD", "data/results"),
            cwd=ROOT, capture_output=True, text=True, timeout=10, check=False)
    except (OSError, subprocess.SubprocessError):      # pragma: no cover
        return None
    if out.returncode != 0:
        return None
    return {Path(line).name for line in out.stdout.split() if line.strip()}


needs_git = pytest.mark.skipif(_tracked() is None,
                               reason="not a git checkout")
needs_results = pytest.mark.skipif(not RESULTS.exists(),
                                   reason="no results directory")


@needs_git
@needs_results
def test_every_results_file_is_in_the_repository():
    """A result that exists only on the machine that produced it is not
    evidence for anything."""
    on_disk = {p.name for p in RESULTS.glob("*.json")}
    tracked = _tracked()
    untracked = sorted(on_disk - tracked)
    assert not untracked, (
        "these results are audited by check_results.py but are not committed:"
        "\n  " + "\n  ".join(untracked)
        + "\n\ndata/results is tracked on purpose -- it is what docs/"
          "evaluation.md cites.\nStaging is not publishing: commit them, or "
          "they are evidence only you can see.")


@needs_git
@needs_results
def test_no_results_file_is_tracked_but_gone():
    """The reverse: a file committed that nothing regenerates is a figure
    with no way left to check it."""
    on_disk = {p.name for p in RESULTS.glob("*.json")}
    missing = sorted(_tracked() - on_disk)
    assert not missing, (
        "tracked results that no longer exist on disk: " + ", ".join(missing))


def _committed(name: str) -> dict | None:
    """The COMMITTED content of a results file, not the working copy.

    The working copy is whatever the last local run wrote, and on the
    machine that runs the ingests it is current -- which is exactly why
    reading it answers the wrong question. What a reader of the repository
    gets is the blob in HEAD.
    """
    out = subprocess.run(("git", "show", f"HEAD:data/results/{name}"),
                         cwd=ROOT, capture_output=True, text=True,
                         timeout=10, check=False)
    if out.returncode != 0:
        return None
    try:
        return json.loads(out.stdout)
    except json.JSONDecodeError:                       # pragma: no cover
        return None


@needs_git
def test_every_committed_result_carries_provenance():
    """check_results reports an unstamped file and, without --strict, does
    not fail on one -- and `make check-results` does not pass --strict. So
    an unattributable result can sit in the repository behind a green gate.

    Read from HEAD rather than from disk on purpose: the workstation's
    working copy is regenerated and stamped, so checking it would pass
    while the published evidence stayed unstamped.
    """
    unstamped = sorted(n for n in (_tracked() or set())
                       if "provenance" not in ((_committed(n)) or {"provenance": 1}))
    assert not unstamped, (
        f"{len(unstamped)} committed results carry no provenance block:\n  "
        + "\n  ".join(unstamped)
        + "\n\nThey cannot be shown to describe any release, and they are "
          "what docs/evaluation.md\ncites. Regenerate and commit them.")
