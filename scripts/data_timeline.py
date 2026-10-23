#!/usr/bin/env python3
"""What data this graph is built from, when it arrived, and what it covers.

    python scripts/data_timeline.py              # both sections
    python scripts/data_timeline.py --acquired   # raw inputs only
    python scripts/data_timeline.py --evidence   # edge dates only

TWO DIFFERENT TIMELINES, and conflating them is the mistake this separates.

  ACQUIRED  when a file landed on this machine. Answers "is this release
            built from a current ChEMBL" and "when was ORCS last pulled".

  COVERED   when the EVIDENCE was first asserted, from first_asserted_date
            on the edges themselves. Answers "what period does this graph
            describe", which is the question a reader of the paper has and
            the one the temporal split at 2021 depends on.

A file downloaded yesterday can carry evidence from 1998. Reporting only the
download date would say this is a current graph; reporting only the edge
dates would say nothing about whether a newer release exists.

MTIME IS NOT A PUBLICATION DATE. A file's mtime is when it arrived here, or
when it was last touched -- copying a directory rewrites it. Where a source
states its own version it is read instead and the mtime is shown beside it;
where it does not, the mtime is all there is and is labelled as such rather
than dressed up as a release date.

THE INPUT LIST IS NOT MAINTAINED HERE. It comes from kgav.pipeline.LAYERS,
which the rebuild driver already uses, so a layer added there appears here
without a second edit -- and cannot be quietly missing from one of them.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.pipeline import LAYERS, layer_dir

# A version a source states about itself, read from the path rather than
# guessed. ChEMBL numbers its releases and puts the number in the directory.
VERSION_IN_PATH = re.compile(r"(chembl_\d+)")


def _when(p: Path) -> str:
    return datetime.fromtimestamp(p.stat().st_mtime,
                                  UTC).strftime("%Y-%m-%d")


def _sidecar(p: Path) -> dict | None:
    """fetch_raw and filter_intact write <file>.source.json beside output."""
    for cand in (p.with_suffix(p.suffix + ".source.json"),
                 p.with_suffix(".source.json")):
        if cand.exists():
            try:
                return json.loads(cand.read_text())
            except (OSError, ValueError):
                return None
    return None


def inputs_for(layer, root: Path) -> list[Path]:
    out = [root / f for f in layer.files]
    for g in layer.globs:
        out.extend(sorted(root.glob(g)))
    return out


def acquired(root: Path) -> list[dict]:
    """One row per raw input the pipeline declares."""
    rows = []
    for layer in LAYERS:
        paths = inputs_for(layer, root)
        present = [p for p in paths if p.exists()]
        if not present:
            rows.append({"layer": layer.name, "path": "-", "state": "absent",
                         "bytes": 0, "when": "-", "version": "-",
                         "provenance": "-"})
            continue
        newest = max(present, key=lambda p: p.stat().st_mtime)
        version = "-"
        m = VERSION_IN_PATH.search(str(newest))
        if m:
            version = m.group(1)
        side = _sidecar(newest)
        if side and side.get("written"):
            version = side.get("derived_from", version)
        rows.append({
            "layer": layer.name,
            "path": (f"{len(present)} files" if len(present) > 1
                     else newest.name),
            "state": "present",
            "bytes": sum(p.stat().st_size for p in present),
            "when": _when(newest),
            "version": version,
            "provenance": "sidecar" if side else "mtime only",
        })
    return rows


def covered(release: Path) -> tuple[Counter, dict]:
    """Edge counts by year of first_asserted_date, and by source.

    Read from the assembled edges rather than from any layer's own report,
    because assembly merges an edge from two sources and keeps the EARLIEST
    date -- so a per-layer count would double-count and date it late.
    """
    years: Counter = Counter()
    by_source: dict[str, Counter] = {}
    path = release / "edges.jsonl"
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            e = json.loads(line)
            d = str(e.get("first_asserted_date") or "")
            y = d[:4] if len(d) >= 4 and d[:4].isdigit() else "undated"
            years[y] += 1
            for src in str(e.get("primary_knowledge_source") or "?").split("|"):
                by_source.setdefault(src, Counter())[y] += 1
    return years, by_source


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--release", type=Path,
                    default=ROOT / "data/releases/v0.1")
    ap.add_argument("--acquired", action="store_true")
    ap.add_argument("--evidence", action="store_true")
    args = ap.parse_args(argv)
    both = not (args.acquired or args.evidence)

    if both or args.acquired:
        print("ACQUIRED — when each raw input landed on this machine")
        print("  A date here is a FILE MTIME unless the version column names "
              "a release.\n  It is not when the source published.\n")
        print(f"  {'layer':<14}{'input':<26}{'bytes':>14}  {'arrived':<12}"
              f"{'version':<22}{'provenance'}")
        for r in acquired(args.root):
            size = f"{r['bytes']:,}" if r["bytes"] else "-"
            print(f"  {r['layer']:<14}{r['path']:<26}{size:>14}  "
                  f"{r['when']:<12}{r['version']:<22}{r['provenance']}")
        out = layer_dir(args.root, "spine").parent
        rel = out / "v0.1"
        if (rel / "MANIFEST.json").exists():
            print(f"\n  assembled release {rel.name}: "
                  f"{_when(rel / 'MANIFEST.json')}")
        print()

    if both or args.evidence:
        if not (args.release / "edges.jsonl").exists():
            print(f"no release at {args.release}; "
                  f"run make rebuild for the evidence section")
            return 0 if args.acquired else 2
        years, by_source = covered(args.release)
        total = sum(years.values())
        print("COVERED — when the evidence was first asserted")
        print("  from first_asserted_date on the assembled edges. Assembly "
              "keeps the\n  EARLIEST date when two sources carry one fact, "
              "so this is the date the\n  claim entered the literature, not "
              "the date it was downloaded.\n")
        named = sorted(y for y in years if y != "undated")
        if named:
            print(f"  range {named[0]}–{named[-1]} over {total:,} edges")
        undated = years.get("undated", 0)
        if undated:
            print(f"  {undated:,} edges ({undated / total:.1%}) carry NO "
                  f"usable date. The temporal split\n  at 2021 decides these "
                  f"by its --undated policy, not by evidence.")
        print()
        width = max((v for k, v in years.items() if k != "undated"),
                    default=1)
        for y in named:
            n = years[y]
            bar = "#" * max(1, round(40 * n / width))
            print(f"    {y}  {n:>8,}  {bar}")
        print("\n  by source, earliest dated edge first")
        order = sorted(by_source,
                       key=lambda s: min((y for y in by_source[s]
                                          if y != "undated"), default="9999"))
        print(f"    {'source':<34}{'edges':>9}  {'first':<7}{'last':<7}"
              f"{'undated':>9}")
        for src in order:
            c = by_source[src]
            dated = sorted(y for y in c if y != "undated")
            print(f"    {src:<34}{sum(c.values()):>9,}  "
                  f"{(dated[0] if dated else '-'):<7}"
                  f"{(dated[-1] if dated else '-'):<7}"
                  f"{c.get('undated', 0):>9,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
