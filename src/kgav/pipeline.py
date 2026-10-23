"""The ingest chain: what runs, in what order, and what each layer needs.

WHY THIS EXISTS AS DATA. The order is not obvious and not recoverable from
LAYER_PRECEDENCE, which ranks layers for qualifier conflicts and says nothing
about when to run them. Three layers -- selectivity, similarity and gtop --
read the ASSEMBLED release rather than the per-layer ones, so they cannot run
until assemble_release has written v0.1, and assemble_release must then run
again to pick them up. A hand-typed sequence gets that wrong silently: the
layer still builds, from the PREVIOUS release, and nothing says so.

PHASES

  1  reads raw inputs and per-layer releases only
  2  assemble_release -> data/releases/v0.1
  3  reads the assembled v0.1
  4  assemble_release again, now including the phase-3 layers

One pass of 3 then 4 is enough: no phase-3 layer's output depends on its own
input, so the second assemble is a fixed point rather than the start of an
iteration.

REQUIRED AND OPTIONAL. spine and host are the graph; without them there is
nothing to assemble and the run stops. Every other layer has raw inputs that
may legitimately be absent in a given checkout -- data/raw is gitignored and
was populated by hand -- and a missing input is a SKIP that gets reported by
name, not a failure. The distinction matters because assemble_release writes
a release whatever is present: a layer silently absent produces a release
whose counts all look right and whose metapaths reach nothing.

ONE DEFINITION, TWO CALLERS. scripts/dryrun_promote.py reported some of this
inline and scripts/rebuild_release.py needs the same facts. This repository
has paid five times for the same measurement living in two places, so the
availability check lives here and both read it.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

CHEMBL_DB = Path("data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")


@dataclass(frozen=True)
class Layer:
    name: str                       # matches its LAYER_PRECEDENCE entry
    script: str                     # scripts/<script>
    phase: int                      # 1 or 3; see the module docstring
    required: bool = False
    files: tuple[str, ...] = ()     # paths under the repo root, all needed
    globs: tuple[str, ...] = ()     # at least one match needed
    env: tuple[str, ...] = ()       # environment variables that must be set
    env_len: int = 0                # exact length required, 0 for any


LAYERS: tuple[Layer, ...] = (
    # ---------------------------------------------------------- phase 1
    Layer("spine", "ingest_spine.py", 1, required=True,
          globs=("data/raw/proteomes/*.json.gz",)),
    Layer("host", "ingest_host.py", 1, required=True,
          globs=("data/raw/host/*",)),
    Layer("vhppi", "ingest_vhppi.py", 1,
          globs=("data/raw/vhppi/vhn_*.tab27",)),
    # ORCS needs a key as well as the screen index: ingest_orcs exits before
    # any request when ORCS_KEY is not 32 characters, so a run without it
    # looks like a layer that produced nothing.
    Layer("orcs", "ingest_orcs.py", 1,
          files=("data/raw/orcs/screens_all.json",),
          env=("ORCS_KEY",), env_len=32),
    Layer("curated", "ingest_curated.py", 1,
          globs=("data/raw/curated_screens/*.tsv",)),
    Layer("chembl", "ingest_chembl.py", 1, files=(str(CHEMBL_DB),)),
    # hosttargets reads the orcs release, so it follows orcs within phase 1.
    Layer("hosttargets", "ingest_host_targets.py", 1,
          files=(str(CHEMBL_DB),)),
    Layer("ncats", "ingest_ncats.py", 1,
          files=("data/raw/ncats/cpe_assay.csv", "data/raw/ncats/cpe_tox.csv")),
    Layer("folds", "ingest_folds.py", 1,
          files=("data/raw/folds/phase2_chains.json",)),
    # ---------------------------------------------------------- phase 3
    # These three read data/releases/v0.1, the ASSEMBLED release.
    Layer("selectivity", "ingest_selectivity.py", 3, files=(str(CHEMBL_DB),)),
    Layer("similarity", "ingest_similarity.py", 3),
    Layer("gtop", "ingest_gtop.py", 3,
          globs=("data/raw/gtop/interactions.csv",)),
)

ASSEMBLE = "assemble_release.py"


def availability(layer: Layer, root: Path,
                 env: dict | None = None) -> tuple[bool, str]:
    """(can it run, why not).

    The reason names the thing that is absent, because "skipped" on its own
    sends the reader to look for a bug in a layer that simply had no input.
    """
    env = os.environ if env is None else env
    for var in layer.env:
        value = str(env.get(var, ""))
        if not value:
            return False, f"{var} is not set"
        if layer.env_len and len(value) != layer.env_len:
            return False, (f"{var} is {len(value)} characters, not "
                           f"{layer.env_len}")
    for rel in layer.files:
        if not (root / rel).exists():
            return False, f"no {rel}"
    for pattern in layer.globs:
        if not list(root.glob(pattern)):
            return False, f"nothing matches {pattern}"
    return True, ""


def plan(root: Path, env: dict | None = None) -> list[dict]:
    """The steps a rebuild would run, in order, each with its verdict.

    Returned rather than printed so the caller decides what to do with a
    skip, and so a test can read the order without parsing output.
    """
    steps: list[dict] = []
    for phase in (1, 3):
        for layer in (x for x in LAYERS if x.phase == phase):
            ok, why = availability(layer, root, env)
            steps.append({"kind": "layer", "name": layer.name,
                          "script": layer.script, "phase": phase,
                          "run": ok, "reason": why,
                          "required": layer.required})
        steps.append({"kind": "assemble", "name": "assemble",
                      "script": ASSEMBLE, "phase": phase + 1,
                      "run": True, "reason": "", "required": True})
    return steps


def blocking(root: Path, env: dict | None = None) -> list[str]:
    """Required layers that cannot run. A rebuild should not start."""
    return [f"{s['name']}: {s['reason']}" for s in plan(root, env)
            if s["kind"] == "layer" and s["required"] and not s["run"]]
