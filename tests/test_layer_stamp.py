"""Which code wrote a layer release.

WHAT THIS IS FOR. On 2026-11-02 src/kgav/gtop.py was patched to emit typed
potency qualifiers in place of `affinity_nm`. `make rebuild` then failed
with 4,270 UNDECLARED_QUALIFIER 'affinity_nm' -- the IDENTICAL count as
before the patch -- and did so twice before anyone looked at the directory.

The cause is structural. gtop is a PHASE 3 layer: it reads the assembled
release, so it runs after assemble_release, and assemble_release runs again
to pick it up. The phase-2 assemble of each run therefore consumes the
phase-3 layer output of the PREVIOUS run. assemble_release validates and
returns 1; rebuild_release stops on a non-zero exit. So the run died on the
old gtop output before the new ingest_gtop.py could replace it, and
`rm -rf data/releases/v0.1-gtop` fixed it -- an absent layer is left out, a
stale one is not.

Nothing anywhere said "this directory was written by code that has since
changed". pipeline.plan's docstring already calls a reused layer the more
dangerous of its two skips; this is the measurement that makes the danger
visible, and it is the same stamp provenance.write_derived gives a derived
release, with the same four statuses.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.provenance import (
    LAYER_STAMP_NAME,
    check_layer,
    source_closure,
    source_sha,
    write_layer,
)

#: The real closure of the layer this was written for, as the parser finds
#: it: the script, its own module, and what that module imports.
GTOP_CLOSURE = ["scripts/ingest_gtop.py", "src/kgav/assemble.py",
                "src/kgav/emit.py", "src/kgav/gtop.py"]


@pytest.fixture
def tree(tmp_path: Path) -> Path:
    """A repository with one real layer's sources copied into it.

    REAL SOURCES, NOT A STUB. The closure is computed by parsing imports,
    so a fixture with hand-written one-line modules would exercise the
    parser on a shape the repository does not contain -- and the first
    bug in this code was a transitive import being missed, which only a
    real multi-hop closure shows.
    """
    (tmp_path / "scripts").mkdir()
    (tmp_path / "src" / "kgav").mkdir(parents=True)
    for rel in GTOP_CLOSURE:
        shutil.copy(ROOT / rel, tmp_path / rel)
    d = tmp_path / "data" / "releases" / "v0.1-gtop"
    d.mkdir(parents=True)
    (d / "edges.jsonl").write_text("{}\n")
    return tmp_path


def _script(tree: Path) -> Path:
    return tree / "scripts" / "ingest_gtop.py"


def _layer(tree: Path) -> Path:
    return tree / "data" / "releases" / "v0.1-gtop"


# ----------------------------------------------------------- the closure
def test_the_closure_is_transitive(tree):
    """ingest_gtop.py imports kgav.gtop, which imports kgav.assemble. A
    change to the merge rules changes that layer's edges while leaving
    both the script and kgav/gtop.py untouched, so a one-hop closure
    would call the layer current."""
    got = sorted(str(p.relative_to(tree)) for p in
                 source_closure(_script(tree), tree))
    assert got == GTOP_CLOSURE


def test_only_kgav_is_followed(tree):
    """site-packages is pinned elsewhere, and hashing it would make every
    layer stale on an unrelated upgrade."""
    for p in source_closure(_script(tree), tree):
        assert p.is_relative_to(tree), p


def test_a_module_that_is_not_imported_is_not_in_the_closure(tree):
    (tree / "src" / "kgav" / "unrelated.py").write_text("X = 1\n")
    before = source_sha(_script(tree), tree)
    (tree / "src" / "kgav" / "unrelated.py").write_text("X = 2\n")
    assert source_sha(_script(tree), tree) == before, (
        "a module outside the closure moved the hash, so every layer "
        "would go stale on any change anywhere")


def test_a_module_becoming_a_package_moves_the_hash(tree):
    """WHAT THE PATH TERM IN THE HASH IS ACTUALLY FOR, and it took a
    surviving mutant to find out.

    The first version of this test swapped two modules' contents and
    asserted the hash moved. It passed with the path term DELETED: the
    content hashes are concatenated in path-sorted order, so a swap
    reorders the sequence and moves the hash on its own. The test named
    the right property and exercised none of it.

    The case the path term does catch is a file whose IMPORTERS do not
    change: `from kgav.gtop import x` resolves to either
    src/kgav/gtop.py or src/kgav/gtop/__init__.py, so converting a
    module to a package moves the file with identical content and no
    edit anywhere else. Without the path in the hash that is one content
    hash in the same position -- indistinguishable from no change at
    all.

    Every other rename inside the closure edits an import statement and
    is caught by content alone, which is why this is the one shape the
    term is load-bearing for.
    """
    k = tree / "src" / "kgav"
    before = source_sha(_script(tree), tree)
    body = (k / "gtop.py").read_text()
    (k / "gtop.py").unlink()
    (k / "gtop").mkdir()
    (k / "gtop" / "__init__.py").write_text(body)
    after = source_sha(_script(tree), tree)
    assert [str(p.relative_to(tree)) for p in source_closure(_script(tree), tree)
            if "gtop" in str(p) and "scripts" not in str(p)] == [
        "src/kgav/gtop/__init__.py"], "the package form did not resolve"
    assert after != before, (
        "kgav.gtop moved from a module to a package with byte-identical "
        "content and the hash did not move; the closure is keyed on "
        "content alone")


# AN EQUIVALENT MUTANT, RECORDED RATHER THAN COVERED. Deleting the
# `a.name.startswith("kgav")` filter on plain `import x` statements does not
# change behaviour: a candidate only enters the closure if src/<name>.py or
# src/<name>/__init__.py exists, so `import json` resolves to nothing either
# way. The filter is there to say what is intended, and if src/yaml.py were
# ever vendored, following it would be correct rather than wrong. There is no
# test for it because there is nothing to observe.


def test_an_unparseable_module_does_not_crash_the_stamp(tree):
    """A syntax error mid-edit must not make the audit unrunnable: the
    thing being audited is a build, and builds are inspected while
    broken."""
    (tree / "src" / "kgav" / "gtop.py").write_text("def (\n")
    assert source_sha(_script(tree), tree)
    write_layer(_layer(tree), _script(tree), tree)
    assert check_layer(_layer(tree), _script(tree), tree)[0] == "ok"


# ------------------------------------------------------- the four statuses
def test_absent_when_the_layer_was_never_built(tree):
    status, why = check_layer(tree / "data" / "releases" / "v0.1-nope",
                              _script(tree), tree)
    assert status == "absent"
    assert "never built" in why


def test_unstamped_when_it_exists_with_no_record(tree):
    """Every layer built before this existed reads this way. It is the
    absence the stamp ends, not an error."""
    status, why = check_layer(_layer(tree), _script(tree), tree)
    assert status == "unstamped"
    assert LAYER_STAMP_NAME in why


def test_ok_after_a_stamp_against_unchanged_code(tree):
    write_layer(_layer(tree), _script(tree), tree)
    assert check_layer(_layer(tree), _script(tree), tree)[0] == "ok"


@pytest.mark.parametrize("rel", GTOP_CLOSURE)
def test_stale_when_any_file_in_the_closure_changes(tree, rel):
    """EVERY file, parametrised. The bug this guards is a closure that
    covers the script and misses what it imports, which a test touching
    only the script would pass."""
    write_layer(_layer(tree), _script(tree), tree)
    p = tree / rel
    p.write_text(p.read_text() + "\n# touched\n")
    status, why = check_layer(_layer(tree), _script(tree), tree)
    assert status == "stale", f"{rel} changed and the layer still reads ok"
    assert rel in why, f"the message does not name {rel}: {why}"


def test_the_stale_message_distinguishes_changed_from_added(tree):
    """A file whose content moved, one added to the closure and one
    dropped from it are three different changes, and "the hash differs"
    distinguishes none of them."""
    write_layer(_layer(tree), _script(tree), tree)
    (tree / "src" / "kgav" / "extra.py").write_text("Y = 1\n")
    s = _script(tree)
    s.write_text("from kgav.extra import Y\n" + s.read_text())
    status, why = check_layer(_layer(tree), s, tree)
    assert status == "stale"
    assert "new in its closure: src/kgav/extra.py" in why, why
    assert "changed: scripts/ingest_gtop.py" in why, why


def test_an_unreadable_stamp_reads_unstamped_not_ok(tree):
    """A truncated write must not pass. `ok` on an unreadable stamp is
    worse than no stamp: it is a positive claim nothing supports."""
    write_layer(_layer(tree), _script(tree), tree)
    (_layer(tree) / LAYER_STAMP_NAME).write_text("{not json")
    assert check_layer(_layer(tree), _script(tree), tree)[0] == "unstamped"


def test_stamping_a_directory_that_does_not_exist_is_a_no_op(tree):
    write_layer(tree / "data" / "releases" / "v0.1-nope", _script(tree), tree)
    assert not (tree / "data" / "releases" / "v0.1-nope").exists()


def test_the_stamp_records_the_commit_and_the_per_file_hashes(tree):
    write_layer(_layer(tree), _script(tree), tree)
    doc = json.loads((_layer(tree) / LAYER_STAMP_NAME).read_text())
    assert doc["script"] == "ingest_gtop.py"
    assert sorted(doc["file_shas"]) == GTOP_CLOSURE
    assert set(doc["code"]) == {"commit", "dirty"}
    assert doc["source_sha"] == source_sha(_script(tree), tree)


# --------------------------------------- through rebuild_release and the audit
#
# A measurement nothing reads is not a check. These drive the two callers.
def _repo(tmp_path: Path, stale_phase3: bool) -> Path:
    """A tree where the plan can be computed, with gtop built and stamped.

    spine and host are `required`, and rebuild_release refuses to start
    without them -- so their raw inputs exist here as empty files.
    `availability` checks for presence, not content.
    """
    for rel in GTOP_CLOSURE:
        (tmp_path / Path(rel).parent).mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, tmp_path / rel)
    for rel in ("data/raw/proteomes/x.json.gz", "data/raw/host/x.tsv",
                "data/raw/gtop/interactions.csv", "data/raw/gtop/ligands.csv"):
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).write_text("")
    d = tmp_path / "data" / "releases" / "v0.1-gtop"
    d.mkdir(parents=True)
    (d / "edges.jsonl").write_text("{}\n")
    write_layer(d, tmp_path / "scripts" / "ingest_gtop.py", tmp_path)
    if stale_phase3:
        m = tmp_path / "src" / "kgav" / "gtop.py"
        m.write_text(m.read_text() + "\n# the patch that was never reached\n")
    return tmp_path


def test_the_plan_carries_the_code_status_for_every_layer(tmp_path):
    from kgav.pipeline import plan
    steps = [s for s in plan(_repo(tmp_path, False)) if s["kind"] == "layer"]
    assert {s["name"]: s["code"] for s in steps}["gtop"] == "ok"
    assert all("code" in s for s in steps)


def test_rebuild_names_a_stale_phase_three_layer_and_its_consequence(
        tmp_path, capsys):
    """THE WHOLE POINT. gtop reused from before its module was patched is
    read by the phase-2 assemble of the next run, which validates it and
    exits 1 -- so the ingest that would fix it never runs. The message
    has to say that, and say to delete the directory."""
    import rebuild_release
    rebuild_release.main(["--plan", "--root", str(_repo(tmp_path, True))])
    out = capsys.readouterr().out
    assert "WRITTEN BY CODE THAT HAS SINCE CHANGED" in out
    assert "src/kgav/gtop.py" in out
    assert "PHASE 3" in out
    # THE CONDITION, NOT JUST THE WORDS. gtop's raw inputs exist in this
    # fixture, so its ingest IS queued to run -- exactly as on 2026-11-02.
    # A warning gated on "will not be re-run" stays silent here, which is
    # silent through the only failure it was written for.
    assert "even though it is queued to run this pass" in out
    assert "the phase-2 assemble runs first" in out
    assert "Delete the directory instead" in out
    assert "absent layer is left out" in out


def test_rebuild_says_nothing_about_code_when_everything_matches(
        tmp_path, capsys):
    """A warning that is always on is a warning nobody reads."""
    import rebuild_release
    rebuild_release.main(["--plan", "--root", str(_repo(tmp_path, False))])
    out = capsys.readouterr().out
    assert "WRITTEN BY CODE THAT HAS SINCE CHANGED" not in out
    assert "UNSTAMPED" not in out


def test_rebuild_reports_an_unstamped_layer_without_calling_it_an_error(
        tmp_path, capsys):
    import rebuild_release
    tree = _repo(tmp_path, False)
    (tree / "data/releases/v0.1-gtop" / LAYER_STAMP_NAME).unlink()
    rebuild_release.main(["--plan", "--root", str(tree)])
    out = capsys.readouterr().out
    assert "UNSTAMPED" in out and "gtop" in out
    assert "Not an error" in out


def test_the_audit_shows_the_code_column_and_explains_a_stale_one(
        tmp_path, capsys):
    import audit_build
    import yaml
    tree = _repo(tmp_path, True)
    cfg = tree / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [], "defaults": {}}))
    audit_build.main(["--config", str(cfg), "--root", str(tree),
                      "--release", str(tree / "data/releases/v0.1")])
    out = capsys.readouterr().out
    assert "code" in out.split("STALE")[0]
    assert "STALE — written by code that has since changed" in out
    assert "src/kgav/gtop.py" in out


def test_the_audit_resolves_the_script_name_from_the_layer_table(tmp_path):
    """hosttargets is ingest_host_targets.py. `f"ingest_{name}.py"` would
    hash a script that does not exist and report the layer unstamped
    forever, which looks exactly like a layer nobody has rebuilt."""
    import audit_build
    assert audit_build._script_for("hosttargets") == "ingest_host_targets.py"
    with pytest.raises(KeyError):
        audit_build._script_for("not-a-layer")
