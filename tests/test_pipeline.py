"""The ingest order, checked against the scripts rather than against itself.

A phase table written by hand is a second copy of a fact that already lives
in the scripts' own defaults. The test that matters here re-derives which
layers read the ASSEMBLED release and compares, so a new layer added in the
wrong phase fails instead of quietly building from the previous release.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.assemble import LAYER_PRECEDENCE
from kgav.pipeline import ASSEMBLE, LAYERS, availability, blocking, plan


# --------------------------------------------------------------- the table
def test_every_layer_names_a_script_that_exists():
    for layer in LAYERS:
        assert (ROOT / "scripts" / layer.script).exists(), layer


def test_the_assemble_script_exists():
    assert (ROOT / "scripts" / ASSEMBLE).exists()


def test_the_layers_are_exactly_the_precedence_entries():
    """A layer here with no precedence entry is never assembled; a
    precedence entry with no layer here is never built."""
    assert {x.name for x in LAYERS} == set(LAYER_PRECEDENCE)


def test_no_layer_appears_twice():
    names = [x.name for x in LAYERS]
    assert len(names) == len(set(names))


# ------------------------------------------- the phase, re-derived from source
def _reads_assembled(script: str) -> bool:
    """Does this script default to data/releases/v0.1, the ASSEMBLED one.

    `v0.1-<name>` is a per-layer release and needs no prior assemble; a bare
    `v0.1` is the output of assemble_release and does.

    TWO SPELLINGS SINCE fix-list L14. The scripts used to name the path as a
    literal and now call `pipeline.release_dir(root)` for the assembled
    release and `layer_dir(root, name)` for a per-layer one. The distinction
    this function draws is exactly the distinction between those two
    helpers, so it reads the calls as well as the literals -- and this test
    failing on the conversion is the test doing its job: putting selectivity
    or similarity in phase 1 still builds, from the PREVIOUS release, and
    nothing in the output says so.
    """
    src = (ROOT / "scripts" / script).read_text()
    if re.search(r"\brelease_dir\s*\(", src):
        return True
    for m in re.finditer(r'releases["\'/\s]*[/"\s]*["\']?v0\.1(-?[a-z]*)', src):
        if not m.group(1):
            return True
    return False


@pytest.mark.parametrize("layer", LAYERS, ids=lambda x: x.name)
def test_NO_LAYER_READS_THE_ASSEMBLED_RELEASE(layer):
    """THE INVARIANT THAT REPLACED THE PHASE MODEL (L20, option (b)).

    This test used to assert the opposite shape: that a layer reading
    data/releases/v0.1 must sit in phase 3, after an assemble, with a
    second assemble to follow. That arrangement is what made L20 and L25
    possible -- the middle assemble consumed the PREVIOUS run's output of
    the very layers that were about to be rebuilt, so a moved subject set
    meant non-convergence and a renamed qualifier meant a hard SCHEMA FAIL.

    Three layers read the assembly and now read `pipeline.peer_nodes` /
    `peer_edges` instead, which serve the layer directories built before
    them. Measured before the change: that union is the assembled node set
    exactly, 78,789 both ways with zero difference, and for similarity's
    partition phase 1 + gtop gives the identical 9,791 compounds.

    So the invariant is stronger and simpler than a phase number: NOTHING
    reads the assembly, therefore ONE assemble suffices and it goes last.
    """
    assert not _reads_assembled(layer.script), (
        f"{layer.name} reads the ASSEMBLED release. That reintroduces the "
        f"mid-pipeline assemble, and with it L20's loop and L25's hard "
        f"failure. Read the layers before it with pipeline.peer_nodes / "
        f"peer_edges instead")


def test_a_peer_reading_layer_declares_it_and_comes_after_its_dependencies():
    """`reads_peers` is what earns a later phase, and the ORDER inside that
    phase is a real dependency: similarity partitions on TARGETS and gtop
    contributes TARGETS, so similarity must follow it. selectivity emits
    HAS_ANTIVIRAL_ACTIVITY_AGAINST, which is in neither set, so nothing
    depends on it."""
    order = [x.name for x in LAYERS]
    peer = {x.name for x in LAYERS if x.reads_peers}
    assert peer == {"selectivity", "gtop", "similarity"}, peer
    assert order.index("similarity") > order.index("gtop"), (
        "similarity partitions compounds by their TARGETS subjects and gtop "
        "contributes TARGETS; built first, similarity sees 1,082 fewer "
        "annotated compounds and gives them cold-start edges they do not need")
    for x in LAYERS:
        assert (x.phase > 1) == x.reads_peers, (
            f"{x.name}: phase {x.phase} and reads_peers={x.reads_peers}")


def test_there_is_exactly_one_assemble_and_it_is_last(tmp_path):
    """The point of (b). Two assembles existed only because three layers
    needed one to have run before them."""
    steps = plan(tmp_path, env={})
    assembles = [i for i, s in enumerate(steps) if s["kind"] == "assemble"]
    assert len(assembles) == 1, (
        f"{len(assembles)} assembles; (b) leaves exactly one")
    assert assembles[0] == len(steps) - 1, "the assemble must be the last step"


def test_hosttargets_follows_orcs():
    """It reads the orcs release, so within phase 1 the order is load-
    bearing even though both are phase 1."""
    order = [x.name for x in LAYERS]
    assert order.index("hosttargets") > order.index("orcs")


def test_spine_is_first_and_host_second():
    assert [x.name for x in LAYERS][:2] == ["spine", "host"]


# ------------------------------------------------------------- availability
def test_a_missing_file_is_named_not_just_counted(tmp_path):
    layer = next(x for x in LAYERS if x.name == "ncats")
    ok, why = availability(layer, tmp_path, env={})
    assert not ok
    assert "cpe_assay.csv" in why


def test_a_missing_glob_names_the_pattern(tmp_path):
    layer = next(x for x in LAYERS if x.name == "vhppi")
    ok, why = availability(layer, tmp_path, env={})
    assert not ok and "vhn_*.tab27" in why


def test_an_unset_key_is_distinguished_from_a_wrong_length_one(tmp_path):
    """ingest_orcs exits before any request when ORCS_KEY is not exactly 32
    characters, and `export ORCS_KEY=...` is a real thing people type."""
    layer = next(x for x in LAYERS if x.name == "orcs")
    (tmp_path / "data/raw/orcs").mkdir(parents=True)
    (tmp_path / "data/raw/orcs/screens_all.json").write_text("[]")

    ok, why = availability(layer, tmp_path, env={})
    assert not ok and "not set" in why

    ok, why = availability(layer, tmp_path, env={"ORCS_KEY": "..."})
    assert not ok and "3 characters" in why

    ok, why = availability(layer, tmp_path, env={"ORCS_KEY": "a" * 32})
    assert ok and why == ""


def test_the_screen_index_is_required_as_well_as_the_key(tmp_path):
    """ingest_orcs reads screens_all.json rather than fetching it."""
    layer = next(x for x in LAYERS if x.name == "orcs")
    ok, why = availability(layer, tmp_path, env={"ORCS_KEY": "a" * 32})
    assert not ok and "screens_all.json" in why


def test_a_layer_with_no_declared_inputs_always_runs(tmp_path):
    layer = next(x for x in LAYERS if x.name == "similarity")
    ok, why = availability(layer, tmp_path, env={})
    assert ok and why == ""


# -------------------------------------------------------------------- plan
def test_every_layer_precedes_the_assemble(tmp_path):
    """There is one assemble and it is last, so EVERY layer precedes it --
    which is the whole of the ordering constraint now. The two tests this
    replaced asserted an assemble between the phases and the phase-1 set
    before it; both described the arrangement (b) removed."""
    steps = plan(tmp_path, env={})
    first = next(i for i, s in enumerate(steps) if s["kind"] == "assemble")
    before = {s["name"] for s in steps[:first] if s["kind"] == "layer"}
    assert before == {x.name for x in LAYERS}
    assert not [s for s in steps[first:] if s["kind"] == "layer"], (
        "no layer may run after the assemble; one that did would not be in "
        "the release it was built from")


def test_a_skipped_layer_still_appears_in_the_plan(tmp_path):
    """Absent from the output is how a missing layer became invisible in the
    first place. It is listed with its reason, not dropped."""
    steps = plan(tmp_path, env={})
    names = {s["name"] for s in steps if s["kind"] == "layer"}
    assert names == {x.name for x in LAYERS}
    assert all(not s["run"] for s in steps
               if s["kind"] == "layer" and s["name"] != "similarity")


def test_blocking_names_only_required_layers(tmp_path):
    stop = blocking(tmp_path, env={})
    assert sorted(s.split(":")[0] for s in stop) == ["host", "spine"]


# ------------------------------------------------------------------ driver
def test_the_driver_refuses_to_start_without_a_graph(tmp_path, capsys):
    from importlib import import_module
    mod = import_module("rebuild_release")
    assert mod.main(["--root", str(tmp_path), "--plan"]) == 2
    out = capsys.readouterr().out
    assert "REFUSING TO START" in out
    assert "spine" in out and "host" in out


def test_the_plan_flag_runs_nothing(tmp_path, monkeypatch):
    """--plan must not shell out, or a dry run rebuilds the release."""
    import subprocess
    from importlib import import_module

    called = []
    monkeypatch.setattr(subprocess, "call",
                        lambda *a, **k: called.append(a) or 0)
    mod = import_module("rebuild_release")
    for rel in ("data/raw/proteomes", "data/raw/host"):
        (tmp_path / rel).mkdir(parents=True)
    (tmp_path / "data/raw/proteomes/UP1.json.gz").write_bytes(b"")
    (tmp_path / "data/raw/host/x").write_bytes(b"")
    assert mod.main(["--root", str(tmp_path), "--plan"]) == 0
    assert called == []


# ------------------------------------- the third state: reused, not absent
#
# assemble_release builds its layer set as
#   layers = {k: v for k, v in declared.items() if v.exists()}
# so a skipped INGEST whose directory is already on disk does not drop the
# layer -- it carries the old one forward. The first version of this module
# reported every skip as "absent from the release", which was wrong for
# exactly the case that matters most.

def _ready(root: Path) -> None:
    """Minimum for the driver to agree to start."""
    (root / "data/raw/proteomes").mkdir(parents=True, exist_ok=True)
    (root / "data/raw/host").mkdir(parents=True, exist_ok=True)
    (root / "data/raw/proteomes/UP1.json.gz").write_bytes(b"")
    (root / "data/raw/host/x").write_bytes(b"")


def test_a_skipped_layer_with_a_directory_on_disk_is_reused(tmp_path):
    _ready(tmp_path)
    (tmp_path / "data/releases/v0.1-orcs").mkdir(parents=True)
    step = next(s for s in plan(tmp_path, env={}) if s["name"] == "orcs")
    assert not step["run"]
    assert step["existing"] is True
    assert step["written"] is not None


def test_a_skipped_layer_with_no_directory_is_absent(tmp_path):
    _ready(tmp_path)
    step = next(s for s in plan(tmp_path, env={}) if s["name"] == "orcs")
    assert not step["run"] and step["existing"] is False


def test_the_report_separates_reuse_from_absence(tmp_path, capsys):
    """THE CORRECTION. One of these produces a zero you can see; the other
    produces plausible numbers from a stale graph."""
    from importlib import import_module
    _ready(tmp_path)
    (tmp_path / "data/releases/v0.1-orcs").mkdir(parents=True)
    mod = import_module("rebuild_release")
    assert mod.main(["--root", str(tmp_path), "--plan"]) == 0
    out = capsys.readouterr().out
    assert "REUSE      orcs" in out
    assert "REUSED, not rebuilt: orcs" in out
    assert "absent     vhppi" in out
    assert "orcs" not in out.split("ABSENT from the release")[1].split("\n")[0]


def test_the_report_no_longer_calls_every_skip_absent(tmp_path, capsys):
    from importlib import import_module
    _ready(tmp_path)
    (tmp_path / "data/releases/v0.1-orcs").mkdir(parents=True)
    mod = import_module("rebuild_release")
    mod.main(["--root", str(tmp_path), "--plan"])
    out = capsys.readouterr().out
    assert "Skipped layers are absent from the release" not in out


def test_the_reuse_warning_names_the_register_as_the_risk(tmp_path, capsys):
    """Promote a virus, skip the orcs ingest, and the new taxon is built
    with no host factors while every count reads correctly. That is the
    failure the warning exists for, so it says it."""
    from importlib import import_module
    _ready(tmp_path)
    (tmp_path / "data/releases/v0.1-orcs").mkdir(parents=True)
    mod = import_module("rebuild_release")
    mod.main(["--root", str(tmp_path), "--plan"])
    out = capsys.readouterr().out
    assert "config/viruses.yaml" in out
    assert "will look complete anyway" in out


def test_layer_dir_matches_what_assemble_reads():
    """Derived from assemble_release rather than restated, so a change to
    its naming fails here instead of making every reuse check a no-op."""
    from kgav.pipeline import layer_dir
    src = (ROOT / "scripts" / "assemble_release.py").read_text()
    assert 'f"v0.1-{name}"' in src, (
        "assemble_release no longer names layer directories this way; "
        "layer_dir is now checking for something that never exists")
    assert layer_dir(Path("/r"), "orcs") == Path("/r/data/releases/v0.1-orcs")


# --------------------------------------- a layer that reads either of two
def _layer(name):
    return next(x for x in LAYERS if x.name == name)


def _touch(root, rel):
    p = root / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("x")
    return p


def test_vhppi_runs_on_intact_alone(tmp_path):
    """The case the old declaration got wrong.

    ingest_vhppi globs vhn_*.tab27 AND intact_*.tab27 -- two sources, one
    layer. The layer declared only the first, and `globs` is an AND, so a
    release holding the filtered IntAct file and no VirHostNet file was
    reported un-runnable while the ingest would have run fine.
    """
    _touch(tmp_path, "data/raw/vhppi/intact_viral.tab27")
    ok, why = availability(_layer("vhppi"), tmp_path, env={})
    assert ok, why


def test_vhppi_runs_on_virhostnet_alone(tmp_path):
    _touch(tmp_path, "data/raw/vhppi/vhn_sars2.tab27")
    ok, why = availability(_layer("vhppi"), tmp_path, env={})
    assert ok, why


def test_vhppi_with_neither_names_both(tmp_path):
    """A reason naming one of two sources sends the reader to fetch the
    wrong thing."""
    ok, why = availability(_layer("vhppi"), tmp_path, env={})
    assert not ok
    assert "vhn_*.tab27" in why and "intact_*.tab27" in why


# ------------------------------------- a layer that needs both, one of which
# may arrive under either of two names
def test_gtop_is_not_ready_on_interactions_alone(tmp_path):
    """It was reported ready with half its inputs, and ingest_gtop then
    exits on the missing ligand file -- a plan that says go, followed by a
    run that will not."""
    _touch(tmp_path, "data/raw/gtop/interactions.csv")
    ok, why = availability(_layer("gtop"), tmp_path, env={})
    assert not ok
    assert "ligands.csv" in why and "ligand_id_mapping.csv" in why


@pytest.mark.parametrize("ligand_file", ["ligands.csv",
                                         "ligand_id_mapping.csv"])
def test_gtop_accepts_either_ligand_file(tmp_path, ligand_file):
    """Both spellings appear across GtP releases and ingest_gtop takes
    either."""
    _touch(tmp_path, "data/raw/gtop/interactions.csv")
    _touch(tmp_path, f"data/raw/gtop/{ligand_file}")
    ok, why = availability(_layer("gtop"), tmp_path, env={})
    assert ok, why


def test_gtop_still_requires_interactions(tmp_path):
    _touch(tmp_path, "data/raw/gtop/ligands.csv")
    ok, why = availability(_layer("gtop"), tmp_path, env={})
    assert not ok and "interactions.csv" in why


def test_the_declared_alternatives_match_what_the_scripts_glob():
    """The declaration is a claim about another file, so it is checked
    against that file rather than against itself."""
    vh = (ROOT / "scripts" / "ingest_vhppi.py").read_text()
    for pat in ("vhn_*.tab27", "intact_*.tab27"):
        assert pat in vh, f"{pat} is declared but ingest_vhppi does not glob it"
    gt = (ROOT / "scripts" / "ingest_gtop.py").read_text()
    for name in ("interactions.csv", "ligands.csv", "ligand_id_mapping.csv"):
        assert name in gt, f"{name} is declared but ingest_gtop never names it"


# ------------------------------------------------- the one layer that fetches
# A REBUILD MUST NOT DEPEND ON A THIRD-PARTY API BEING UP. Every ingest but
# one reads what fetch_raw put under data/raw, which is what makes a rebuild
# reproducible. ingest_folds.py fetches from InterPro's live REST API at
# INGEST time, and on 2026-11-07 a transient "HTTP Error 408: Request
# Timeout" from InterPro ended `make rebuild` 6m31s in, after seven layers
# had been rewritten -- stopping before the assemble, so the release stayed
# consistent, which is the one part that worked.
#
# ingest_folds.py has had an --offline flag replaying its cache since it was
# written, and pipeline.py already declared that cache as the layer's input.
# There was simply no way to pass a flag: Layer had no field for one and
# rebuild_release called every script bare.
def test_only_the_folds_layer_declares_arguments():
    """If a second layer ever needs a flag, this should be read and changed
    deliberately rather than silently widened."""
    withargs = {ly.name: ly.args for ly in LAYERS if ly.args}
    assert withargs == {"folds": ("--offline",)}


def test_the_folds_layer_is_rebuilt_from_its_cache():
    folds = next(ly for ly in LAYERS if ly.name == "folds")
    assert "--offline" in folds.args, (
        "the rebuild would fetch from InterPro's live API; a 408 from it "
        "has already cost one rebuild")
    # The flag's own default cache IS the file this layer declares as input,
    # so the declaration and the behaviour cannot disagree.
    assert "data/raw/folds/phase2_chains.json" in folds.files


def test_the_plan_carries_each_layer_s_arguments(tmp_path):
    """A plan that hides the flags cannot be used to check what will run."""
    steps = plan(ROOT)
    folds = next(s for s in steps
                 if s["kind"] == "layer" and s["name"] == "folds")
    assert folds["args"] == ["--offline"]
    assert all("args" in s for s in steps), \
        "a step with no args key makes rebuild_release's lookup silent"


# --------------------------------------------------- where the release lands
def test_assemble_writes_where_release_dir_says(tmp_path):
    """L13 LEANS ON THIS. rebuild_release.py fingerprints the release
    directory before and after a run to decide whether it changed anything;
    if that path and assemble_release's --out default ever drift, the
    comparison silently measures a directory nobody writes and every run
    reports "no release on disk". Read off the parser rather than grepped
    for the literal, which is why build_parser was extracted."""
    import assemble_release

    from kgav.pipeline import release_dir
    ap = assemble_release.build_parser(tmp_path)
    out = {a.dest: a.default for a in ap._actions}["out"]
    assert out == release_dir(tmp_path)
    # And it must sit beside the per-layer directories, not inside one.
    from kgav.pipeline import layer_dir
    assert layer_dir(tmp_path, "spine").parent == release_dir(tmp_path).parent


# ------------------------------------------- what a rebuild did to data/results
#: L13. The closing line asserted staleness: every run ended "Results in
#: data/results were computed from the previous release and are now stale".
#: The 2026-11-10 re-stamp rebuild came back byte-identical and
#: check-results reported 27 current, 0 stale -- so the sentence was false
#: of the one run a reader most needs it right about.
def _fp(nodes, edges):
    return {"id": "v0.1", "nodes_sha": nodes, "edges_sha": edges}


def test_an_unchanged_release_is_not_reported_as_staling_anything():
    from rebuild_release import staleness_verdict
    said = staleness_verdict(_fp("n", "e"), _fp("n", "e"))
    assert "did not change the release" in said
    assert "now stale" not in said


def test_a_changed_release_is_reported_as_staling_the_results():
    from rebuild_release import staleness_verdict
    said = staleness_verdict(_fp("n1", "e1"), _fp("n1", "e2"))
    assert "now stale" in said
    # Both fingerprints named, so the claim is checkable against MANIFEST.
    assert "e1" in said and "e2" in said


def test_a_changed_NODES_sha_counts_as_changed_too():
    """The edge file is the one that usually moves, which is exactly why a
    comparison on edges alone would pass a node-only change through."""
    from rebuild_release import staleness_verdict
    assert "now stale" in staleness_verdict(_fp("n1", "e"), _fp("n2", "e"))


def test_an_unchanged_release_does_not_claim_the_results_are_CURRENT():
    """THE SCOPE IS NARROWER THAN THE GOOD NEWS. An unchanged release proves
    this run staled nothing; it says nothing about whether the files were
    already stale from an earlier run. Only check_results.py reads each
    file's own stamp, so the pointer to it must survive in every branch."""
    from rebuild_release import staleness_verdict
    for before, after in [(_fp("n", "e"), _fp("n", "e")),
                          (_fp("n1", "e1"), _fp("n1", "e2")),
                          (_fp(None, None), _fp("n", "e")),
                          (_fp("n", "e"), None)]:
        said = staleness_verdict(before, after)
        assert "check-results" in said, said
        assert "are current" not in said
        assert "nothing is stale" not in said


def test_no_release_before_the_run_is_not_reported_as_unchanged():
    """A fresh checkout has no data/releases/v0.1 -- `_sha256` returns None
    for a missing file, so both shas are None and a naive equality check
    would call that "unchanged" and tell the reader nothing went stale."""
    from rebuild_release import staleness_verdict
    for before in (None, _fp(None, None)):
        said = staleness_verdict(before, _fp("n", "e"))
        assert "no release on disk before this run" in said
        assert "did not change the release" not in said
