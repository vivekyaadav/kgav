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
from kgav.temporal import PLACEHOLDER_DATE, edge_year

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


def covered(release: Path) -> tuple[Counter, dict, dict]:
    """Edge counts by year, by source, and the distinct dates per source.

    Read from the assembled edges rather than from any layer's own report,
    because assembly merges an edge from two sources and keeps the EARLIEST
    date -- so a per-layer count would double-count and date it late.

    UNDATED IS temporal.edge_year's DEFINITION, NOT A SECOND ONE. The first
    version of this read the year off the string, so the 1970-01-01
    placeholder -- which src/kgav/temporal.py names PLACEHOLDER_DATE and
    every ingest writes when it has no date -- appeared as 6,287 edges of
    1970 evidence, and `undated` read 0 for every source. The repository
    already decides this question in one place and assemble._earliest
    already special-cases the same string; a second answer here would be
    the duplication that comment was written about.
    """
    years: Counter = Counter()
    by_source: dict[str, Counter] = {}
    dates: dict[str, set] = {}
    path = release / "edges.jsonl"
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            e = json.loads(line)
            y = edge_year(e)
            key = str(y) if y else "undated"
            years[key] += 1
            raw = str(e.get("first_asserted_date") or "")
            for src in str(e.get("primary_knowledge_source") or "?").split("|"):
                by_source.setdefault(src, Counter())[key] += 1
                if y:
                    dates.setdefault(src, set()).add(raw)
    return years, by_source, dates


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
        years, by_source, dates = covered(args.release)
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
                  f"usable date -- {PLACEHOLDER_DATE}, which every ingest "
                  f"writes when it\n  has none. The temporal split at 2021 "
                  f"decides these by its --undated policy,\n  not by "
                  f"evidence.")
        print()
        # A SOURCE THAT STAMPS EVERY EDGE WITH ONE DATE IS NOT DATED.
        # ingest_host writes string_release_date to all 141,019 STRING
        # edges and ingest_vhppi defaults all of its own to 2024-01-01 --
        # the FILE's release date, not when any fact was established.
        # Drawn in the histogram those make a single year tower over every
        # real one and read as a discovery burst, so they are separated
        # out and reported as what they are.
        # TWO OR MORE DATED EDGES, or there is nothing to observe: one
        # edge trivially has one date, and calling that a stamp said
        # chembl stamps its edges on a fixture holding a single chembl row.
        def _n_dated(src):
            return sum(v for k, v in by_source[src].items() if k != "undated")

        stamped = {s_: next(iter(d)) for s_, d in dates.items()
                   if len(d) == 1 and _n_dated(s_) >= 2}
        stamp_years = {v[:4] for v in stamped.values()}
        per_edge = {y: n for y, n in years.items()
                    if y != "undated" and y not in stamp_years}
        if stamped:
            print("  NOT DRAWN BELOW -- every dated edge from these sources "
                  "carries the SAME date.\n  That is typically the file's "
                  "release date rather than when any fact was\n  "
                  "established, and in a histogram it towers over every real "
                  "year:")
            for s_, d in sorted(stamped.items(), key=lambda kv: kv[1]):
                n = sum(v for k, v in by_source[s_].items() if k != "undated")
                print(f"    {s_:<34}{n:>9,}  all dated {d}")
            print()
        if per_edge:
            print("  per-edge dates only")
            width = max(per_edge.values(), default=1)
            for y in sorted(per_edge):
                n = per_edge[y]
                bar = "#" * max(1, round(40 * n / width))
                print(f"    {y}  {n:>8,}  {bar}")
        print("\n  by source, earliest dated edge first")
        order = sorted(by_source,
                       key=lambda s: min((y for y in by_source[s]
                                          if y != "undated"), default="9999"))
        print(f"    {'source':<34}{'edges':>9}  {'first':<7}{'last':<7}"
              f"{'undated':>9}  dating")
        for src in order:
            c = by_source[src]
            dated = sorted(y for y in c if y != "undated")
            kind = ("none" if not dated
                    else "ONE DATE" if src in stamped
                    else "per-edge")
            print(f"    {src:<34}{sum(c.values()):>9,}  "
                  f"{(dated[0] if dated else '-'):<7}"
                  f"{(dated[-1] if dated else '-'):<7}"
                  f"{c.get('undated', 0):>9,}  {kind}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
