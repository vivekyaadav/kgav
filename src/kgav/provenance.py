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


def code_commit() -> dict[str, Any]:
    """The commit that produced a result, and whether the tree was dirty.

    A dirty tree is recorded rather than refused: results are frequently
    generated mid-change, and a result that says so is more useful than one
    that cannot be produced at all.
    """
    def _git(*args: str) -> str | None:
        try:
            out = subprocess.run(("git", *args), capture_output=True, text=True,
                                 cwd=Path(__file__).resolve().parents[2],
                                 timeout=10, check=False)
            return out.stdout.strip() if out.returncode == 0 else None
        except (OSError, subprocess.SubprocessError):
            return None

    commit = _git("rev-parse", "HEAD")
    status = _git("status", "--porcelain")
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
def check_file(path: Path, release: Path) -> tuple[str, str]:
    """(status, message) for one results file against a release.

    status is 'ok', 'stale' or 'unstamped'. 'unstamped' is reported rather
    than assumed fine: a file with no provenance is exactly the case that went
    undetected, and its numbers cannot be attributed to any graph.
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
    return "ok", (f"matches {current['id']} "
                  f"(code {(doc[PROVENANCE_KEY].get('code') or {}).get('commit', '?')[:8]})")


def check_dir(results_dir: Path, release: Path) -> list[tuple[str, str, str]]:
    """(filename, status, message) for every .json in a results directory."""
    out = []
    for p in sorted(Path(results_dir).glob("*.json")):
        status, msg = check_file(p, release)
        out.append((p.name, status, msg))
    return out
