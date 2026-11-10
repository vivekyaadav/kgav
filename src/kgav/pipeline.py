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

ONE PASS OF 3 THEN 4 IS NOT ENOUGH, AND THIS SAID IT WAS. The claim here
read "no phase-3 layer's output depends on its own input, so the second
assemble is a fixed point rather than the start of an iteration". The first
half is true and the conclusion does not follow, because phase 2 assembles
the PREVIOUS run's phase-3 layer directories: what similarity reads at
phase 3 contains gtop's output from the last rebuild, and gtop is itself a
phase 3 layer. ingest_similarity splits compounds into query and reference
by their TARGETS / INHIBITS / MEASURED_INACTIVE_AGAINST subjects, and gtop
contributes TARGETS -- so the dependency is mutual between two phase-3
layers through the assembly that sits between them, not self-dependence.

MEASURED 2026-11-10, reverting fix-list B1, a phase-1 change: rebuild one
went edges_sha af01e12ce8a86891 -> 8402bbaf9ed5c531, differing from the
intended release by exactly one CHEMICALLY_SIMILAR_TO edge; rebuild two
went -> f8a1ae20bfd89b37, byte-identical to the release before B1. Two
passes converged, and nothing in the pipeline checks that they did:
rebuild_release.py prints "The release CHANGED this run" either way, which
is true on an unconverged run and on a genuinely new one, and a reader
cannot tell which from the line. Filed as L20. data_audit_2026-11-07's P17
predicted this and was closed with "document the fixed point" -- so the
document asserting one is the fix that was applied.

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

from kgav.provenance import check_layer

CHEMBL_DB = Path("data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")


@dataclass(frozen=True)
class Layer:
    name: str                       # matches its LAYER_PRECEDENCE entry
    script: str                     # scripts/<script>
    phase: int                      # 1 or 3; see the module docstring
    required: bool = False
    files: tuple[str, ...] = ()     # paths under the repo root, all needed
    globs: tuple[str, ...] = ()     # EVERY pattern must match something
    # ALTERNATIVES: each inner tuple is a set of patterns of which at least
    # one must match. `globs` is an AND over its patterns, which cannot
    # express a layer that reads either of two sources -- and two layers
    # here do. Declaring them with `globs` understated one and overstated
    # the other, in opposite directions.
    any_globs: tuple[tuple[str, ...], ...] = ()
    env: tuple[str, ...] = ()       # environment variables that must be set
    env_len: int = 0                # exact length required, 0 for any
    # ARGUMENTS THE REBUILD MUST PASS. Empty for every layer that reads
    # data/raw and writes a layer directory, which is all of them but one.
    #
    # folds is the exception and it cost a rebuild to find out. It fetches
    # from InterPro's live REST API at INGEST time, so `make rebuild` is
    # the only pipeline step that depends on a third-party service being
    # up: on 2026-11-07 a transient `HTTP Error 408: Request Timeout` from
    # InterPro killed the run 6m31s in, after seven layers had already
    # been rewritten. ingest_folds.py has had an --offline flag replaying
    # its cache the whole time, and this file already declares that cache
    # as the layer's input -- there was simply no way to pass the flag.
    #
    # Offline is the RIGHT default for a rebuild, not a workaround. Every
    # other ingest reads what fetch_raw put on disk, which is what makes
    # it reproducible; an ingest that fetches is unversioned either way
    # (fix list V3: no InterPro version is recorded anywhere), so the
    # online path gives an irreproducible result AND can fail. Refreshing
    # the cache stays available as an explicit run without this flag.
    args: tuple[str, ...] = ()


LAYERS: tuple[Layer, ...] = (
    # ---------------------------------------------------------- phase 1
    Layer("spine", "ingest_spine.py", 1, required=True,
          globs=("data/raw/proteomes/*.json.gz",)),
    Layer("host", "ingest_host.py", 1, required=True,
          globs=("data/raw/host/*",)),
    # TWO SOURCES, ONE LAYER, and ingest_vhppi globs both vhn_*.tab27 and
    # intact_*.tab27. Declaring only the first meant a plan never named
    # IntAct among this layer's inputs, so a reader checking whether IntAct
    # was in the build saw nothing -- and a release holding the filtered
    # IntAct file but no VirHostNet file would have been reported as
    # un-runnable while the ingest would in fact have run.
    Layer("vhppi", "ingest_vhppi.py", 1,
          any_globs=(("data/raw/vhppi/vhn_*.tab27",
                      "data/raw/vhppi/intact_*.tab27"),)),
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
    #
    # IT FOLLOWED ncats FOR ONE DAY AND DOES NOT ANY MORE (fix-list B1).
    # Option D relaxed the phase gate for compounds the graph already has,
    # which hosttargets can only know about from layers already built, and
    # the move was worth 660 of the 1,327 recovered compounds -- the ncats
    # half. B1 was then DECLINED, because the admission rule turned out to
    # correlate with the labels, so the order is back where v0.1 was built
    # and `--phase-gate relaxed` sees only the chembl half. That is recorded
    # in host_targets.PHASE_GATE_POLICIES rather than fixed, because the
    # right order depends on a decision nobody should re-make by reading a
    # tuple. Nothing else depends on the order either way: ncats reads
    # neither hosttargets nor the assembled release.
    Layer("hosttargets", "ingest_host_targets.py", 1,
          files=(str(CHEMBL_DB),)),
    Layer("ncats", "ingest_ncats.py", 1,
          files=("data/raw/ncats/cpe_assay.csv", "data/raw/ncats/cpe_tox.csv")),
    Layer("folds", "ingest_folds.py", 1,
          args=("--offline",),
          files=("data/raw/folds/phase2_chains.json",)),
    # ---------------------------------------------------------- phase 3
    # These three read data/releases/v0.1, the ASSEMBLED release.
    Layer("selectivity", "ingest_selectivity.py", 3, files=(str(CHEMBL_DB),)),
    Layer("similarity", "ingest_similarity.py", 3),
    # ingest_gtop needs interactions.csv AND a ligand file, which may be
    # either ligands.csv or ligand_id_mapping.csv. Declaring only
    # interactions.csv reported the layer ready with half its inputs, and
    # the run would then exit on the missing one.
    Layer("gtop", "ingest_gtop.py", 3,
          globs=("data/raw/gtop/interactions.csv",),
          any_globs=(("data/raw/gtop/ligands.csv",
                      "data/raw/gtop/ligand_id_mapping.csv"),)),
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
    for alternatives in layer.any_globs:
        if not any(list(root.glob(p)) for p in alternatives):
            return False, ("nothing matches any of "
                           + " or ".join(alternatives))
    return True, ""


def layer_dir(root: Path, name: str) -> Path:
    """Where assemble_release looks for this layer's output."""
    return root / "data" / "releases" / f"v0.1-{name}"


def layers_before(name: str) -> tuple[str, ...]:
    """Layer names built before `name`, in order.

    EXISTS SO A LAYER CANNOT GUESS WHAT IS ALREADY ON DISK. hosttargets
    relaxes its phase gate for compounds the graph already has (fix-list
    B1), and "already has" means exactly these layers -- a hand-written
    list would be a second copy of LAYERS and would go stale the next time
    a layer moves, silently narrowing or widening the gate.
    """
    names = [l.name for l in LAYERS]
    if name not in names:
        raise ValueError(f"{name!r} is not a layer; have {names}")
    return tuple(names[:names.index(name)])


def release_dir(root: Path) -> Path:
    """Where assemble_release WRITES, as opposed to where it reads.

    assemble_release.py's `--out` default reads this, and rebuild_release.py
    uses it to fingerprint the release before and after a run (fix-list L13).
    Those two must agree or the fingerprint measures a directory nobody
    writes, which is the reason this function exists.

    IT IS NOW THE ONLY COPY (fix-list L14, converted 2026-11-11). This
    docstring used to say it was not, and named "nine other scripts" and
    "ten copies of one path". BOTH FIGURES WERE STALE WHEN COUNTED: there
    were 28 files and 37 `default=` sites -- 39 counting the two
    constructions that are not argument defaults, which the test named
    below found after the migration had already run -- because every script
    added after the row was written copied the literal from the one
    beside it. A
    paragraph that admits to being one of ten is read as bounded; this one
    had been growing for as long as it had been apologising, which is the
    argument for converting rather than documenting.

    tests/test_release_path_has_one_owner.py now fails if the literal
    reappears anywhere but here.
    """
    return root / "data" / "releases" / "v0.1"


def train_dir(root: Path, year: int = 2021) -> Path:
    """Where temporal_split WRITES its training cut, and others read it.

    The third shape of this path, after `release_dir` and `layer_dir`.
    temporal_split.py writes it; hard_negatives.py and check_results.py
    read it, and check_results additionally names it in its `--train` help
    text. Keeping the year a parameter rather than baking 2021 in is not
    speculative generality: `docs/prereg_temporal_auc.md` states the cutoff
    as a protocol choice, and a second cut is the obvious next request.
    """
    return root / "data" / "releases" / f"v0.1-train{year}"


def plan(root: Path, env: dict | None = None) -> list[dict]:
    """The steps a rebuild would run, in order, each with its verdict.

    Returned rather than printed so the caller decides what to do with a
    skip, and so a test can read the order without parsing output.

    `existing` IS THE THIRD STATE, and leaving it out made this module's own
    report wrong on the day it shipped. assemble_release builds its layer set
    as `{k: v for k, v in declared.items() if v.exists()}` -- whatever
    directory is on disk, regardless of whether the ingest ran this time. So
    a skipped layer whose directory already exists is NOT absent from the
    release; it is REUSED, carried over from whenever it was last written.

    That is the more dangerous of the two skips. An absent layer makes a
    metapath reach nothing, which at least shows up as a zero. A reused one
    produces plausible numbers from a graph that no longer matches the
    register -- promote a virus, skip the orcs ingest, and the new taxon is
    built with no host factors while every count reads correctly.
    """
    steps: list[dict] = []
    for phase in (1, 3):
        for layer in (x for x in LAYERS if x.phase == phase):
            ok, why = availability(layer, root, env)
            d = layer_dir(root, layer.name)
            # WHICH CODE WROTE IT, alongside whether it exists. `existing`
            # is what makes a skipped layer get REUSED rather than dropped,
            # and until 2026-11-02 that was the whole story a caller could
            # tell. It is not enough: a phase-3 layer reused from a run
            # before its own module was patched kills the phase-2 assemble
            # of the next run, before the ingest that would fix it. See
            # provenance.check_layer.
            code, code_reason = check_layer(d, root / "scripts" / layer.script,
                                            root)
            steps.append({"kind": "layer", "name": layer.name,
                          "script": layer.script, "phase": phase,
                          "args": list(layer.args),
                          "run": ok, "reason": why,
                          "required": layer.required,
                          "existing": d.exists(),
                          "code": code, "code_reason": code_reason,
                          "written": (d.stat().st_mtime if d.exists()
                                      else None)})
        steps.append({"kind": "assemble", "name": "assemble",
                      "script": ASSEMBLE, "phase": phase + 1,
                      "args": [],
                      "run": True, "reason": "", "required": True,
                      "existing": False, "code": "n/a", "code_reason": "",
                      "written": None})
    return steps


def blocking(root: Path, env: dict | None = None) -> list[str]:
    """Required layers that cannot run. A rebuild should not start."""
    return [f"{s['name']}: {s['reason']}" for s in plan(root, env)
            if s["kind"] == "layer" and s["required"] and not s["run"]]
