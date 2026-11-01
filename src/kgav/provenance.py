"""What a results file was computed FROM.

WHY THIS EXISTS. data/results/metapath_comparison.json and baselines_*.json
were generated on 2026-09-12, the selectivity layer landed on 2026-09-13
adding 1,568 activity edges all carrying relation "=", and nothing recorded
the mismatch. Five months later those files were still being read as the
"before" column of a comparison, and the drift was indistinguishable from the
change under test: SARS-CoV-2's M1 positive count appeared to RISE 736 -> 889
when the true movement on one release was 1,011 -> 889, a fall.

That error was undetectable by inspection. A results file recorded numbers and
nothing about where they came from, so there was no way to ask whether it
still applied to the graph on disk. The fix is not care -- it is that every
results file now carries the fingerprint of the release and the code that
produced it, and check_results.py rejects one whose release has moved.

WHAT IS FINGERPRINTED. The content of nodes.jsonl and edges.jsonl, not the
directory name: a release rebuilt in place keeps its name and is a different
graph. The code commit is recorded too, because the same release read by
different code gives different answers -- which is exactly what the label
polarity fix changed.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

PROVENANCE_KEY = "provenance"

# A RUN THAT REFUSED STILL HAS TO LEAVE A RECORD. temporal_split.py aborts on
# a fatal leakage audit and returns 1 without writing, which left the PREVIOUS
# run's numbers on disk -- computed from a graph three revisions old -- with
# nothing in the file saying the protocol that produced them now refuses to
# run. check_results.py then reported "STALE: regenerate it", and regenerating
# was impossible: the script will not write until the leak is resolved. The
# advice was wrong and the numbers were still readable as current work.
#
# So a refusal writes a stamped record carrying the audit failures and NO
# results block. The dead numbers go, the reason stays, and the file is
# checkable like any other.
REFUSAL_KEY = "refused_by"


def _sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def release_fingerprint(release: Path) -> dict[str, Any]:
    """Identity of the graph a result was computed from.

    Hashes the content rather than trusting the directory name, because a
    release rebuilt in place keeps its name.
    """
    release = Path(release)
    return {
        "id": release.name,
        "nodes_sha": _sha256(release / "nodes.jsonl"),
        "edges_sha": _sha256(release / "edges.jsonl"),
    }


# THE FLAG WAS REPORTING ON ITS OWN OUTPUT. `dirty` lives under `code` and is
# meant to say whether the CODE that produced a number differed from the
# commit beside it. It was computed from a bare `git status --porcelain` over
# the whole tree -- and data/results/*.json are tracked, so the first results
# file a run writes dirties the tree and every stamp after it records
# dirty: True, including later stamps from the same script.
#
# The 2026-10-16 regeneration shows the shape exactly: 13 of the 14 results
# files say dirty: True, and the single False is simply the first one written
# (temporal_2021_exclude.json, 05:20:02; the next is 05:20:45). Nothing was
# wrong with that run -- the clean first stamp is positive evidence the code
# was clean when it started -- but a flag that is True for 13 of 14 files
# whatever happens cannot distinguish that run from a genuinely dirty one.
#
# So the status check excludes the generated-results directory. Nothing reads
# this field today, which is why it could be wrong for as long as it was.
#
# GIT DOES THE PATH MATCHING, NOT US. The first version of this fix parsed
# `git status --porcelain` and compared line[3:] against these prefixes --
# and did not work, because _git() below returns stdout.strip() and the strip
# ate the leading space of the first entry: " M data/results/x.json" became
# "M data/results/x.json", so line[3:] was "ta/results/x.json" and matched
# nothing. The unit test written alongside it passed, because it fed the
# parser a hand-written unstripped string rather than the command's real
# output. A pathspec has no column to miscount, and it covers untracked
# files under the directory too.
#
# `root` EXISTS SO THIS CAN BE TESTED AT ALL. The first test for it ran
# against the repository itself and skipped when the tree was dirty -- which
# is to say it skipped in exactly the situation where someone is changing
# this function, and both the broken versions above passed it by being
# skipped. A test needs a tree whose state it controls, so the caller may
# supply one.
GENERATED_PREFIXES = ("data/results",)
_EXCLUDE = tuple(f":(exclude){p}" for p in GENERATED_PREFIXES)


def code_commit(root: Path | None = None) -> dict[str, Any]:
    """The commit that produced a result, and whether the CODE was dirty.

    Dirty means a tracked file OUTSIDE the generated-results directory
    differs from HEAD. Results files are excluded because writing them is
    what this function's callers do -- see GENERATED_PREFIXES above.

    `root` defaults to this repository and is overridden only by tests.

    A dirty tree is recorded rather than refused: results are frequently
    generated mid-change, and a result that says so is more useful than one
    that cannot be produced at all.
    """
    where = Path(root) if root is not None \
        else Path(__file__).resolve().parents[2]

    def _git(*args: str) -> str | None:
        try:
            out = subprocess.run(("git", *args), capture_output=True, text=True,
                                 cwd=where, timeout=10, check=False)
            return out.stdout.strip() if out.returncode == 0 else None
        except (OSError, subprocess.SubprocessError):
            return None

    commit = _git("rev-parse", "HEAD")
    status = _git("status", "--porcelain", "--", *_EXCLUDE)
    return {"commit": commit,
            "dirty": bool(status) if status is not None else None}


def stamp(release: Path, schema_version: str | None = None,
          **extra: Any) -> dict[str, Any]:
    """The provenance block to embed in a results file."""
    return {
        "release": release_fingerprint(release),
        "code": code_commit(),
        "schema_version": schema_version,
        "generated_utc": _dt.datetime.now(_dt.UTC).isoformat(timespec="seconds"),
        **extra,
    }


def write_results(path: Path, payload: dict, release: Path,
                  schema_version: str | None = None) -> None:
    """Write a results file with its provenance attached.

    Every results writer goes through here, so a new one cannot forget.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    body = dict(payload)
    body[PROVENANCE_KEY] = stamp(release, schema_version)
    path.write_text(json.dumps(body, indent=2))


# ------------------------------------------------------------------ checking
def write_refusal(path: Path, payload: dict, reasons: list[str], release: Path,
                  schema_version: str | None = None) -> None:
    """Record that this protocol refused to run, in place of its results.

    Same provenance stamp as a real result, so the refusal is attributable to
    a graph. Any 'results' key in payload is dropped: a refusal that still
    carried numbers would be read as numbers.
    """
    body = {k: v for k, v in payload.items() if k != "results"}
    body[REFUSAL_KEY] = list(reasons)
    write_results(path, body, release, schema_version)


def check_file(path: Path, release: Path) -> tuple[str, str]:
    """(status, message) for one results file against a release.

    status is 'ok', 'stale', 'unstamped' or 'refused'. 'unstamped' is reported
    rather than assumed fine: a file with no provenance is exactly the case
    that went undetected, and its numbers cannot be attributed to any graph.

    'refused' is checked AFTER staleness, so a refusal record from a graph
    that has since moved still reads stale -- it describes a refusal that may
    no longer happen.
    """
    path = Path(path)
    try:
        doc = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return "stale", f"unreadable: {exc}"
    if not isinstance(doc, dict) or PROVENANCE_KEY not in doc:
        return "unstamped", ("no provenance block; cannot tell which release "
                             "produced these numbers")
    recorded = (doc[PROVENANCE_KEY] or {}).get("release") or {}
    current = release_fingerprint(release)
    # A MISSING RELEASE IS NOT A MATCH. _sha256 returns None for a file that
    # does not exist, and two Nones compared equal -- so checking against a
    # typo'd or deleted release path reported every file as current. That is a
    # false "current" reading from the one tool whose job is to catch false
    # "current" readings, which is worse than no check at all.
    absent = [k for k in ("nodes_sha", "edges_sha") if current.get(k) is None]
    if absent:
        return "stale", (f"cannot verify: {release} has no "
                         f"{', '.join(k.replace('_sha', '.jsonl') for k in absent)}, "
                         f"so there is nothing to compare the stamp against")
    diffs = [k for k in ("nodes_sha", "edges_sha")
             if recorded.get(k) != current.get(k)]
    if diffs:
        return "stale", (
            f"computed from {recorded.get('id')} "
            f"[nodes {recorded.get('nodes_sha')}, edges {recorded.get('edges_sha')}] "
            f"but {current['id']} is now "
            f"[nodes {current['nodes_sha']}, edges {current['edges_sha']}]")
    code = (doc[PROVENANCE_KEY].get('code') or {}).get('commit', '?')[:8]
    reasons = doc.get(REFUSAL_KEY)
    if reasons:
        return "refused", (f"this protocol REFUSED to run on {current['id']} "
                           f"(code {code}); no results were produced:\n"
                           + "\n".join(f"              - {r}" for r in reasons))
    return "ok", f"matches {current['id']} (code {code})"


def check_dir(results_dir: Path, release: Path) -> list[tuple[str, str, str]]:
    """(filename, status, message) for every .json in a results directory."""
    out = []
    for p in sorted(Path(results_dir).glob("*.json")):
        status, msg = check_file(p, release)
        out.append((p.name, status, msg))
    return out


# --------------------------------------------------------- derived releases
# THE HOLE THIS CLOSES. check_file compares a RESULT against a RELEASE, which
# is the whole audit -- and temporal_split.py cuts a second release,
# v0.1-train2021, that hard_negatives.py then reads. That derived graph
# carried no stamp and nothing checked it, so after v0.1 was rebuilt the
# sequence was:
#
#   rebuild v0.1 -> hard_negatives.py loads the OLD v0.1-train2021 (it
#   exists, so it is used) -> write_results stamps the output against the NEW
#   v0.1 -> check_results reports "current".
#
# A number computed from a graph that no longer exists, certified current by
# the tool whose only job is to catch exactly that. Worse than the original
# bug, which at least left the file readable as old.
#
# The derived release now carries a results-shaped stamp of its own, written
# and checked through the same two functions as everything else rather than a
# second fingerprint implementation. docs/evaluation.md calls that the fifth
# time one calculation lived in two files; it is not becoming the sixth.
DERIVED_STAMP_NAME = "DERIVED_FROM.json"


def write_derived(derived_dir: Path, release: Path,
                  schema_version: str | None = None, **policy: Any) -> None:
    """Stamp a derived release with the release and policy that produced it.

    `policy` is recorded because the three --undated protocols write to one
    output directory, so the stamp is also the only record of WHICH split is
    currently on disk. A consumer that reports a temporal number can now say
    which protocol produced the graph it read.
    """
    write_results(Path(derived_dir) / DERIVED_STAMP_NAME, dict(policy),
                  release, schema_version)


def check_derived(derived_dir: Path, release: Path) -> tuple[str, str]:
    """(status, message) for a derived release against its parent release.

    Adds 'absent' to check_file's statuses. A derived release that was never
    built is a different situation from one cut off a graph that has since
    moved: the first needs a first run, the second needs a re-run, and a
    caller that conflates them prints the wrong instruction.
    """
    derived_dir = Path(derived_dir)
    if not derived_dir.exists():
        return "absent", f"{derived_dir} does not exist"
    stamp_path = derived_dir / DERIVED_STAMP_NAME
    if not stamp_path.exists():
        return "unstamped", (
            f"{derived_dir} carries no {DERIVED_STAMP_NAME}, so there is no "
            "way to tell which release it was cut from; rebuild it")
    return check_file(stamp_path, release)


def derived_policy(derived_dir: Path) -> dict[str, Any]:
    """The policy recorded in a derived release's stamp, or {} if unreadable."""
    try:
        doc = json.loads((Path(derived_dir) / DERIVED_STAMP_NAME).read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return {k: v for k, v in doc.items() if k != PROVENANCE_KEY}


# ---------------------------------------------------- what code wrote a layer
#
# THE FAILURE THIS IS FOR, 2026-11-02. src/kgav/gtop.py was patched to emit
# typed potency qualifiers in place of `affinity_nm`. `make rebuild` then
# failed with 4,270 UNDECLARED_QUALIFIER 'affinity_nm' -- the IDENTICAL count
# as before the patch -- and did so twice.
#
# The cause is structural, not a forgotten step. gtop is a PHASE 3 layer: it
# reads the assembled release, so it runs AFTER assemble_release and
# assemble_release runs again to pick it up. That means the phase-2 assemble
# of each run consumes the phase-3 layer output of the PREVIOUS run.
# assemble_release validates and returns 1, rebuild_release stops on a
# non-zero exit, and so the run died on last run's gtop output before this
# run's ingest_gtop.py could replace it. `rm -rf data/releases/v0.1-gtop`
# fixed it, because an absent layer is simply left out.
#
# Nothing was wrong except that no number anywhere said "this directory was
# written by code that has since changed". A layer release carries its edges
# and nothing about its origin, so a reused layer and a freshly built one are
# indistinguishable -- which pipeline.plan's docstring already calls the more
# dangerous of its two skips, for the same reason.
#
# The release-side half of this problem was solved for results files in 2026
# -09 (see the module docstring) and for derived releases by write_derived
# above. This is the same stamp for the third kind of artifact, with the same
# four statuses, because a third vocabulary for one idea is its own defect.
LAYER_STAMP_NAME = "BUILT_BY.json"


def source_closure(script: Path, root: Path) -> list[Path]:
    """The script, plus every `kgav.*` module it reaches by static import.

    WHY STATIC AND NOT sys.modules. The stamp is written by whatever ran the
    ingest, not by the ingest itself -- twelve scripts each write their own
    output and none of them would be stamping anything. Parsing the imports
    needs no cooperation from the script and gives the same answer whether
    the script ran or not, which is also what makes `check_layer` able to
    answer for a layer built months ago.

    WHY THE WHOLE TRANSITIVE SET. ingest_gtop.py imports kgav.gtop, which
    imports kgav.units; a change to unit conversion changes that layer's
    edges while leaving both the script and kgav/gtop.py untouched.

    Only `kgav.*` is followed. The standard library and third-party packages
    are pinned elsewhere and hashing site-packages would make every layer
    stale on an unrelated upgrade.
    """
    import ast

    root, script = Path(root), Path(script)
    src = root / "src"
    out: dict[str, Path] = {}
    if script.exists():
        out[str(script)] = script
    queue, seen = [script], set()
    while queue:
        path = queue.pop()
        if path in seen or not path.exists():
            continue
        seen.add(path)
        try:
            tree = ast.parse(path.read_text())
        except (OSError, SyntaxError):
            continue
        names: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level == 0 \
                    and (node.module or "").startswith("kgav"):
                names.add(node.module)
                # `from kgav.agent import tools` names a module, not an
                # attribute. Both readings are tried below; a name that is
                # not a module simply has no file.
                names.update(f"{node.module}.{a.name}" for a in node.names)
            elif isinstance(node, ast.Import):
                names.update(a.name for a in node.names
                             if a.name.startswith("kgav"))
        for name in names:
            rel = Path(*name.split("."))
            for cand in (src / rel.with_suffix(".py"), src / rel / "__init__.py"):
                if cand.exists():
                    out[str(cand)] = cand
                    queue.append(cand)
    return sorted(out.values())


def source_sha(script: Path, root: Path) -> str:
    """One hash over the layer's whole source closure.

    Keyed on the path RELATIVE to root as well as the content. The content
    hashes are concatenated in path-sorted order, so most renames move the
    hash on their own -- but `from kgav.gtop import x` resolves to either
    src/kgav/gtop.py or src/kgav/gtop/__init__.py, so a module becoming a
    package moves the file with identical content and no edit anywhere
    else. Without the path term that is one content hash in the same
    position, indistinguishable from no change. See
    tests/test_layer_stamp.py, where a mutant deleting this survived a
    test that claimed to cover it.
    """
    root = Path(root)
    h = hashlib.sha256()
    for path in source_closure(script, root):
        try:
            rel = path.relative_to(root)
        except ValueError:
            rel = path
        h.update(str(rel).encode())
        h.update(b"\0")
        h.update((_sha256(path) or "").encode())
        h.update(b"\0")
    return h.hexdigest()[:16]


def write_layer(layer_dir: Path, script: Path, root: Path) -> None:
    """Record which code wrote this layer release.

    Called by whatever drove the ingest, after it succeeded. A layer written
    by hand -- `python scripts/ingest_gtop.py` on its own -- is left
    unstamped rather than wrongly stamped, and check_layer says so.
    """
    layer_dir = Path(layer_dir)
    if not layer_dir.exists():
        return
    root = Path(root)
    closure = source_closure(Path(script), root)
    body = {
        "script": Path(script).name,
        "source_sha": source_sha(Path(script), root),
        # PER-FILE, not just the rolled-up hash. The rolled-up hash answers
        # "is this stale"; these answer "which file moved", which is the
        # question a person actually has in front of a failing rebuild.
        "file_shas": {
            (str(p.relative_to(root)) if p.is_relative_to(root) else str(p)):
            _sha256(p) for p in closure},
        "code": code_commit(root),
        "generated_utc": _dt.datetime.now(_dt.UTC).isoformat(
            timespec="seconds"),
    }
    (layer_dir / LAYER_STAMP_NAME).write_text(json.dumps(body, indent=2))


def check_layer(layer_dir: Path, script: Path, root: Path) -> tuple[str, str]:
    """(status, message) for a layer release against the code on disk.

    The four statuses are check_derived's, and they need different actions:

      absent     never built. The layer is simply left out of the release.
      unstamped  on disk with no record of what wrote it. Every layer built
                 before 2026-11-02 reads this way; it is not an error, it is
                 the absence this stamp exists to end.
      stale      written by code that has since changed. On a phase-3 layer
                 this is the state that kills a rebuild before the layer can
                 be replaced, so it names the files that moved.
      ok         written by the code now on disk.
    """
    layer_dir, root = Path(layer_dir), Path(root)
    if not layer_dir.exists():
        return "absent", f"{layer_dir.name} was never built"
    stamp_path = layer_dir / LAYER_STAMP_NAME
    if not stamp_path.exists():
        return "unstamped", (
            f"{layer_dir.name} carries no {LAYER_STAMP_NAME}, so there is no "
            f"way to tell which code wrote it; rebuild it to find out")
    try:
        doc = json.loads(stamp_path.read_text())
    except (OSError, json.JSONDecodeError):
        return "unstamped", f"{layer_dir.name}'s {LAYER_STAMP_NAME} is unreadable"
    now = source_sha(Path(script), root)
    if doc.get("source_sha") == now:
        return "ok", f"{layer_dir.name} matches the code on disk"
    # NAME THE FILES, in both directions: one whose content moved, one added
    # to the closure, and one removed from it are three different changes and
    # "the hash differs" distinguishes none of them.
    was = doc.get("file_shas") or {}
    now_shas = {
        (str(p.relative_to(root)) if p.is_relative_to(root) else str(p)):
        _sha256(p) for p in source_closure(Path(script), root)}
    moved = sorted(k for k in set(was) & set(now_shas) if was[k] != now_shas[k])
    added = sorted(set(now_shas) - set(was))
    gone = sorted(set(was) - set(now_shas))
    parts = ([f"changed: {', '.join(moved)}"] if moved else []) \
        + ([f"new in its closure: {', '.join(added)}"] if added else []) \
        + ([f"no longer imported: {', '.join(gone)}"] if gone else [])
    detail = ("; " + "; ".join(parts)) if parts else (
        "; no file in its closure differs, so the stamp itself is from "
        "another tree")
    return "stale", (
        f"{layer_dir.name} was written by different code "
        f"(stamp {doc.get('source_sha')}, now {now}){detail}; re-run "
        f"{Path(script).name} or delete the directory")
