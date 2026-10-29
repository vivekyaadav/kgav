"""The Reproducing block runs from the document, or fails loudly.

tests/test_reproduce_block.py checks the block's CONTENT -- that no
results file exists which no step regenerates. This checks that the
block can be extracted and run at all, because until 2026-10-29 it
was transcribed by hand after every rebuild and a dropped flag there
changes a protocol without failing.

Every test here is about refusing to run the wrong thing quietly. The
extractor must not return an empty block, must not silently find the
wrong block, and --print must not execute.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import scripts.run_reproduce_block as RB  # noqa: E402

DOC = """# Evaluation

## Something Else

```bash
echo WRONG BLOCK
```

## Reproducing

Order matters.

```bash
# a comment
python scripts/temporal_split.py --cutoff 2021 --undated include
make check-results
```

## After
"""


def test_it_takes_the_block_under_the_right_heading():
    """A document holds several bash blocks. Taking the first one in
    the file would run whatever happened to appear earliest."""
    body = RB.block(DOC)
    assert "temporal_split.py" in body
    assert "WRONG BLOCK" not in body


def test_the_real_document_still_parses():
    """Pinned against the document that ships, not only a fixture.
    A heading renamed in docs/evaluation.md would break `make
    reproduce` and nothing else would notice."""
    body = RB.block((ROOT / "docs" / "evaluation.md").read_text())
    cmds = [ln for ln in body.splitlines()
            if ln.strip() and not ln.strip().startswith("#")]
    assert len(cmds) > 10, cmds
    assert any("temporal_split.py" in c for c in cmds)
    assert any("check-results" in c for c in cmds)


@pytest.mark.parametrize("text,why", [
    ("# Evaluation\n\nno heading here\n", "heading missing"),
    ("## Reproducing\n\nprose and no block\n", "no bash block"),
    ("## Reproducing\n\n```bash\n# only a comment\n```\n", "no commands"),
])
def test_it_raises_rather_than_returning_an_empty_block(text, why):
    """An empty block would run nothing and exit 0 -- a pass that
    means the opposite of what it looks like, which is the shape of
    every vacuous check removed from this repo."""
    with pytest.raises(LookupError):
        RB.block(text)


def test_a_broken_document_exits_non_zero_without_running(tmp_path, capsys):
    bad = tmp_path / "d.md"
    bad.write_text("nothing here")
    assert RB.main(["--doc", str(bad)]) == 2
    assert "cannot read the Reproducing block" in capsys.readouterr().err


def test_print_shows_the_block_and_runs_nothing(tmp_path, capsys, monkeypatch):
    """--print is what a reader uses before committing 40 minutes to
    it. If it executed, the preview would be the run."""
    doc = tmp_path / "d.md"
    doc.write_text(DOC)
    called = []
    monkeypatch.setattr(RB.subprocess, "call",
                        lambda *a, **k: called.append(a) or 0)
    assert RB.main(["--doc", str(doc), "--print"]) == 0
    assert called == [], "--print executed the block"
    assert "temporal_split.py" in capsys.readouterr().out


def test_the_runner_does_not_use_set_e(tmp_path, monkeypatch):
    """One step is REQUIRED to exit non-zero -- the --overlap-policy
    keep control must refuse. Under `set -e` the run would abort there
    and skip most of the block, and the output would look like a
    completed run that simply ended early."""
    doc = tmp_path / "d.md"
    doc.write_text(DOC)
    seen = {}

    def fake_call(cmd, **kw):
        seen["script"] = Path(cmd[1]).read_text()
        return 0

    monkeypatch.setattr(RB.subprocess, "call", fake_call)
    assert RB.main(["--doc", str(doc)]) == 0
    lines = [ln.strip() for ln in seen["script"].splitlines()]
    # THE DIRECTIVE, NOT THE STRING. `"set -e" not in script` failed on
    # the preamble COMMENT explaining why set -e is absent -- a
    # substring check satisfied by an unrelated occurrence, which is
    # the defect this suite has found most often, caught here by the
    # test itself rather than by review.
    assert "set -u" in lines
    assert "set -e" not in lines
    assert "set -eu" not in lines
    assert any("temporal_split.py" in ln for ln in lines)


def test_the_exit_code_is_the_block_s(tmp_path, monkeypatch):
    """A failing step must fail the target, or `make reproduce` would
    report success over a half-finished regeneration."""
    doc = tmp_path / "d.md"
    doc.write_text(DOC)
    monkeypatch.setattr(RB.subprocess, "call", lambda *a, **k: 3)
    assert RB.main(["--doc", str(doc)]) == 3
