#!/usr/bin/env python3
"""Two results directories: did a rebuild change any NUMBER, or only the stamp?

`make check-results` answers "does this file describe the graph on
disk". This answers a different question, which came up on 2026-10-29
after v0.1 was reassembled and had no tool: did the figures move?

Every volatile field lives under the top-level "provenance" key --
generated_utc, code.commit, code.dirty, release.nodes_sha and
release.edges_sha. Drop that, canonicalise, compare. Anything left is
a real change.

IT REFUSES TO PASS QUIETLY, for a reason worth keeping. The first run
of this comparison happened BEFORE anything was regenerated, so it
compared the backup against its own unchanged source and reported no
differences -- a check that passed because nothing had happened, which
is indistinguishable from a check that passed because nothing changed.
Byte-identical files are therefore a PROBLEM here, not agreement, and
an empty comparison exits non-zero saying so.

--fields also prints which values moved and a histogram of the fields
across all files. That histogram is the diagnosis: one field differing
in twenty-one files is a different finding from twenty-one files each
differing in their own way. On the reassembly it showed 64 moved
values in `top.degree[].drug` and NONE in the matching `.score` --
tied names reshuffling, not numbers changing.

Usage:
    python scripts/compare_results.py <before_dir> <after_dir> [--fields]
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

SKIP_TOP = {"provenance"}
MAX_PER_FILE = 8


def canon(path: Path) -> str:
    d = json.loads(path.read_text())
    if isinstance(d, dict):
        d.pop("provenance", None)
    return json.dumps(d, sort_keys=True, indent=1)


def walk(a, b, path: str = "", out: list | None = None) -> list:
    """Leaf-level differences between two JSON documents."""
    if out is None:
        out = []
    if type(a) is not type(b):
        out.append((path, a, b))
        return out
    if isinstance(a, dict):
        for k in sorted(set(a) | set(b)):
            if not path and k in SKIP_TOP:
                continue
            p = f"{path}.{k}" if path else k
            if k not in a:
                out.append((p, "<absent>", b[k]))
            elif k not in b:
                out.append((p, a[k], "<absent>"))
            else:
                walk(a[k], b[k], p, out)
    elif isinstance(a, list):
        if len(a) != len(b):
            out.append((f"{path}[]", f"len={len(a)}", f"len={len(b)}"))
            return out
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f"{path}[{i}]", out)
    elif a != b:
        out.append((path, a, b))
    return out


def short(v, n: int = 58) -> str:
    s = v if isinstance(v, str) else json.dumps(v)
    return s if len(s) <= n else s[:n - 1] + "…"


def generic(path: str) -> str:
    """Array indices collapsed, so the histogram groups by field."""
    out, depth = [], 0
    for ch in path:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
            out.append("[]")
        elif depth == 0:
            out.append(ch)
    return "".join(out)


def compare(before: Path, after: Path, fields: bool = False) -> int:
    """Returns the number of problems; prints the report."""
    b = {p.name for p in before.glob("*.json")}
    a = {p.name for p in after.glob("*.json")}
    problems = 0
    hist: Counter = Counter()

    for name in sorted(b - a):
        print(f"MISSING AFTER    {name}")
        problems += 1
    for name in sorted(a - b):
        print(f"NEW AFTER        {name}  (not in the backup; cannot compare)")

    stamp_only = 0
    for name in sorted(b & a):
        bp, ap = before / name, after / name
        if canon(bp) == canon(ap):
            if bp.read_bytes() == ap.read_bytes():
                print(f"BYTE-IDENTICAL   {name}  "
                      f"(not regenerated -- the stamp did not move either)")
                problems += 1
            else:
                stamp_only += 1
            continue
        print(f"NUMERIC CHANGE   {name}")
        problems += 1
        if not fields:
            continue
        d = walk(json.loads(bp.read_text()), json.loads(ap.read_text()))
        for p, x, y in d[:MAX_PER_FILE]:
            print(f"    {p}\n        before  {short(x)}\n"
                  f"        after   {short(y)}")
        if len(d) > MAX_PER_FILE:
            print(f"    ... and {len(d) - MAX_PER_FILE} more")
        for p, _, _ in d:
            hist[generic(p)] += 1

    if fields and hist:
        print("\nFIELDS BY HOW MANY VALUES MOVED")
        for p, n in hist.most_common(25):
            print(f"  {n:>7,}  {p}")

    print(f"\n{stamp_only} file(s) differ in the stamp only.")
    print(f"{len(b & a)} file(s) compared, {problems} problem(s).")
    if problems == 0 and stamp_only == 0:
        print("NOTHING WAS COMPARED. Check the paths. This is not a "
              "finding that nothing changed.")
        return 1
    return problems


def main(argv: list[str]) -> int:
    fields = "--fields" in argv
    args = [a for a in argv if not a.startswith("--")]
    if len(args) != 2:
        print(__doc__)
        return 2
    before, after = Path(args[0]), Path(args[1])
    for d in (before, after):
        if not d.is_dir():
            print(f"not a directory: {d}")
            return 2
    return 1 if compare(before, after, fields) else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
