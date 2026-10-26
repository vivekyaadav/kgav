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
    """
    src = (ROOT / "scripts" / script).read_text()
    for m in re.finditer(r'releases["\'/\s]*[/"\s]*["\']?v0\.1(-?[a-z]*)', src):
        if not m.group(1):
            return True
    return False


@pytest.mark.parametrize("layer", LAYERS, ids=lambda x: x.name)
def test_a_layer_reading_the_assembled_release_runs_after_an_assemble(layer):
    """THE ONE THAT MATTERS. selectivity, similarity and gtop read
    data/releases/v0.1, so they cannot run before assemble_release has
    written it -- and assemble must then run again to include them. A layer
    put in phase 1 by mistake still builds, from the PREVIOUS release, and
    nothing in the output says so."""
    expected = 3 if _reads_assembled(layer.script) else 1
    assert layer.phase == expected, (
        f"{layer.name} reads "
        f"{'the assembled' if expected == 3 else 'only per-layer'} release, "
        f"so it belongs in phase {expected}, not {layer.phase}")


def test_at_least_one_layer_is_in_each_phase():
    """If this ever fails the two-phase structure has gone, and the second
    assemble is doing nothing."""
    phases = {x.phase for x in LAYERS}
    assert phases == {1, 3}


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
def test_assemble_runs_between_the_phases_and_again_at_the_end(tmp_path):
    steps = plan(tmp_path, env={})
    kinds = [(s["kind"], s["phase"]) for s in steps]
    assert kinds.count(("assemble", 2)) == 1
    assert kinds.count(("assemble", 4)) == 1
    first = next(i for i, s in enumerate(steps) if s["kind"] == "assemble")
    after = [s for s in steps[first:] if s["kind"] == "layer"]
    assert after and all(s["phase"] == 3 for s in after)


def test_every_phase_one_layer_precedes_the_first_assemble(tmp_path):
    steps = plan(tmp_path, env={})
    first = next(i for i, s in enumerate(steps) if s["kind"] == "assemble")
    before = [s for s in steps[:first] if s["kind"] == "layer"]
    assert {s["name"] for s in before} == {x.name for x in LAYERS
                                           if x.phase == 1}


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
