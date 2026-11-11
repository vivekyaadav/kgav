"""The ingest chain: what runs, in what order, and what each layer needs.

WHY THIS EXISTS AS DATA. The order is not obvious and not recoverable from
LAYER_PRECEDENCE, which ranks layers for qualifier conflicts and says nothing
about when to run them. A hand-typed sequence gets it wrong silently: a layer
built out of order still builds, from stale inputs, and nothing says so.

PHASES

  1  reads raw inputs only
  2  reads the layer directories built before it (`reads_peers`)
  3  assemble_release -> data/releases/v0.1, ONCE, last

NO LAYER READS THE ASSEMBLED RELEASE, and that is what makes one assemble
enough. Until 2026-11-11 three did -- selectivity, similarity and gtop --
so an assemble had to run in the middle, and THAT ASSEMBLE CONSUMED THE
PREVIOUS RUN'S OUTPUT OF THOSE SAME THREE LAYERS. Two defects followed from
it, both measured rather than argued:

  L20  NON-CONVERGENCE. ingest_similarity partitions compounds by their
       TARGETS / INHIBITS / MEASURED_INACTIVE_AGAINST subjects and gtop
       contributes TARGETS, so the two were mutually dependent through the
       assembly between them. Reverting fix-list B1 went edges_sha
       af01e12ce8a86891 -> 8402bbaf9ed5c531, off by exactly one
       CHEMICALLY_SIMILAR_TO edge, then -> f8a1ae20bfd89b37 and stopped.
       Two passes, and the output said "The release CHANGED this run" both
       times.

  L25  A HARD STOP. Renaming a phase-3 layer's qualifier (L11) made the
       middle assemble validate 4,264 edges carrying the now-undeclared
       name and exit 1 -- before the ingest that would have replaced them.
       Not avoidable by running twice: every pass began by assembling the
       same stale directory. It had killed this script twice on 2026-11-02
       the same way.

FIXED 2026-11-11, L20 option (b), chosen by the owner. The three layers read
`peer_nodes` / `peer_edges`, which serve the layer directories named by
`layers_before(name)` in LAYER_PRECEDENCE order -- the assembler's own
precedence, so property resolution matches. MEASURED BEFORE THE CHANGE, and
both substitutions are exact rather than close: the union of those
directories holds the assembled node set precisely (78,789 both ways, zero
difference, and no peer-reading layer emits nodes), and phase 1 + gtop gives
similarity the identical 9,791-compound partition the assembly gave it.
Reading phase 1 ALONE would NOT have been neutral -- 1,082 compounds are
annotated only by gtop and would have moved into the query set -- which is
why the order inside phase 2 is a declared dependency and not a convention.

CONFIRMED BY REBUILD: byte-identical release, nodes_sha c800b9005c8b074c and
edges_sha ad7121bc2195a93a unmoved, and it converged on the FIRST pass. The
loop is gone rather than reported, which is what distinguishes (b) from
option (a): looping until the fingerprint repeated would not have helped
L25 at all, because every iteration re-read the same stale directory.

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
    phase: int                      # 1 or 2; see the module docstring
    # READS THE LAYERS BUILT BEFORE IT, rather than the assembled release
    # (fix-list L20, option (b)). A layer with this set is given the peer
    # layer DIRECTORIES -- `layers_before(its own name)`, resolved in
    # LAYER_PRECEDENCE order so property resolution matches the
    # assembler's -- and the pipeline needs only ONE assemble, at the end.
    #
    # It replaces the arrangement that made L20 and L25 possible: three
    # layers read the ASSEMBLED release, so an assemble had to run in the
    # middle, and that middle assemble consumed the PREVIOUS run's output
    # of those same three layers. Mutual dependence through the assembly
    # between them -- non-convergence when a subject set moved, and a hard
    # SCHEMA FAIL when one of their qualifiers was renamed.
    reads_peers: bool = False
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
    # ------------------------------------------------------------ phase 2
    # These three read the layers built before them, not the assembled
    # release, so there is ONE assemble and it is last (L20 option (b)).
    #
    # THE ORDER WITHIN THIS PHASE IS A REAL DEPENDENCY, not a convention.
    # similarity partitions compounds by their TARGETS / INHIBITS /
    # MEASURED_INACTIVE_AGAINST subjects, and gtop contributes TARGETS, so
    # similarity must come after it. selectivity emits
    # HAS_ANTIVIRAL_ACTIVITY_AGAINST, which is in neither set, so nothing
    # here depends on selectivity and its position is free.
    #
    # MEASURED BEFORE THE CHANGE, because the naive form of (b) -- read
    # only phase 1 -- would have altered meaning: 1,082 compounds are
    # annotated ONLY by gtop and would have moved into similarity's query
    # set, acquiring cold-start edges they do not need. Reading
    # phase 1 + gtop gives the identical 9,791-compound partition the
    # assembly gave, so this ordering is semantically neutral and L20's
    # "changes what annotated means" caveat does not apply to it.
    Layer("selectivity", "ingest_selectivity.py", 2, reads_peers=True,
          files=(str(CHEMBL_DB),)),
    # ingest_gtop needs interactions.csv AND a ligand file, which may be
    # either ligands.csv or ligand_id_mapping.csv. Declaring only
    # interactions.csv reported the layer ready with half its inputs, and
    # the run would then exit on the missing one.
    Layer("gtop", "ingest_gtop.py", 2, reads_peers=True,
          globs=("data/raw/gtop/interactions.csv",),
          any_globs=(("data/raw/gtop/ligands.csv",
                      "data/raw/gtop/ligand_id_mapping.csv"),)),
    Layer("similarity", "ingest_similarity.py", 2, reads_peers=True),
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


def peer_dirs(root: Path, name: str, releases: Path | None = None) -> list[Path]:
    """The layer directories `name` reads, in the order the assembler merges.

    LAYER_PRECEDENCE ORDER, NOT BUILD ORDER, and the distinction is the
    whole point. A caller that unions these dedupes by node id with the
    first writer winning -- which is what `assemble_release` does -- so it
    reproduces the assembly's property resolution only if it walks them in
    the same order. Build order would give a different `smiles` for any
    compound two layers describe differently, and fix-list L18 found 222 of
    those in one layer alone.

    Restricted to `layers_before(name)`, so a peer-reading layer sees
    exactly what was built before it and nothing that depends on it. That
    is what makes the dependency a DAG and lets the pipeline assemble once.

    Directories that do not exist are dropped: data/raw is hand-populated,
    an absent layer is an ordinary state of a checkout, and the caller
    reports what it found rather than failing on what it did not.
    """
    from .assemble import LAYER_PRECEDENCE
    allowed = set(layers_before(name))
    out = []
    for layer in LAYER_PRECEDENCE:
        if layer in allowed:
            # `releases` lets a caller point at a directory of per-layer
            # releases that is not root/data/releases -- which is what the
            # ingests' --peers flag passes, and what a test needs in a
            # tmp_path. Without it the override would have to be
            # reconstructed from the path, which the first version of this
            # did by inspecting the directory NAME.
            d = (releases / f"v0.1-{layer}") if releases else layer_dir(root, layer)
            if (d / "nodes.jsonl").exists() or (d / "edges.jsonl").exists():
                out.append(d)
    return out


def peer_nodes(root: Path, name: str, releases: Path | None = None) -> list[dict]:
    """Every node the layers before `name` hold, deduped, first writer wins.

    THE SUBSTITUTE FOR READING THE ASSEMBLED RELEASE. Measured on v0.1
    before the change: the union of the phase-1 layer directories holds
    78,789 node ids and the assembled release holds 78,789, with zero in
    either direction and zero nodes emitted by any peer-reading layer. So
    for the question every one of them actually asks -- which compounds and
    proteins exist -- the union is not an approximation.

    First writer wins in LAYER_PRECEDENCE order, which is `assemble_release`
    own rule; `ingest_ncats` has read its peers this way since 2026-10-09
    and documents the same reasoning. Concatenating without the dedupe hands
    a validator the same id repeatedly and it correctly reports
    DUPLICATE_NODE, which says nothing about the layer under test.
    """
    import json
    by_id: dict[str, dict] = {}
    for d in peer_dirs(root, name, releases):
        f = d / "nodes.jsonl"
        if not f.exists():
            continue
        for line in f.read_text().splitlines():
            if line.strip():
                node = json.loads(line)
                by_id.setdefault(node["id"], node)
    return list(by_id.values())


def peer_edges(root: Path, name: str, releases: Path | None = None):
    """Every edge the layers before `name` hold, in precedence order.

    Not deduped: the caller asks a question about predicates and subjects,
    and two layers stating the same edge is a merge concern rather than a
    reading one. `ingest_similarity` is the only caller -- it needs the
    subjects carrying TARGETS / INHIBITS / MEASURED_INACTIVE_AGAINST, and a
    subject counted twice is still one member of a set.
    """
    import json
    for d in peer_dirs(root, name, releases):
        f = d / "edges.jsonl"
        if not f.exists():
            continue
        for line in f.read_text().splitlines():
            if line.strip():
                yield json.loads(line)


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
    for phase in sorted({x.phase for x in LAYERS}):
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
    # ONE ASSEMBLE, AND IT IS LAST (fix-list L20, option (b)). There were
    # two, one after each phase, because three layers read the ASSEMBLED
    # release and so needed one to exist before they ran. They now read the
    # layer directories built before them, so nothing needs an assembly
    # until everything has been built -- which removes the loop L20
    # describes and the hard failure L25 describes, rather than reporting
    # them.
    steps.append({"kind": "assemble", "name": "assemble",
                  "script": ASSEMBLE, "phase": max(x.phase for x in LAYERS) + 1,
                  "args": [],
                  "run": True, "reason": "", "required": True,
                  "existing": False, "code": "n/a", "code_reason": "",
                  "written": None})
    return steps


def blocking(root: Path, env: dict | None = None) -> list[str]:
    """Required layers that cannot run. A rebuild should not start."""
    return [f"{s['name']}: {s['reason']}" for s in plan(root, env)
            if s["kind"] == "layer" and s["required"] and not s["run"]]
