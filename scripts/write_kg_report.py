#!/usr/bin/env python3
"""Regenerate the release census inside docs/kg_audit_report.txt.

WHY GENERATED RATHER THAN WRITTEN. The first version of that report
carried no node or edge totals at all, on the reasoning that a hardcoded
count goes stale silently -- which is true, and the wrong conclusion.
An audit report IS a stamped snapshot of one release; omitting its size
makes it useless to a reader and protects nothing. The protection is to
generate the numbers from the release and stamp them with the content
hashes they were computed from, so a stale block is detectable rather
than merely regrettable.

So this script owns everything between the two markers and nothing else
in the file. It reads the assembled release, counts, and writes the
block with the release id and the nodes/edges hashes `check_results.py`
audits against. tests/test_kg_audit_report.py re-derives the same counts
and fails if the block disagrees, so a hand-edit or a rebuild without a
regenerate is caught.

    python scripts/write_kg_report.py           # rewrite the block
    python scripts/write_kg_report.py --check   # fail if it is stale
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.provenance import release_fingerprint
from kgav.schema import load_schema

REPORT = ROOT / "docs" / "kg_audit_report.txt"
BEGIN = "<<<BEGIN GENERATED: release census — python scripts/write_kg_report.py>>>"
END = "<<<END GENERATED>>>"


def _assembled(release: Path):
    """Counts straight off the two jsonl files.

    Deliberately not imported from audit_build: that module prints as it
    goes and takes an argv. The counting is six lines and a second
    copy of six lines is cheaper than a refactor that couples a report
    writer to a CLI.
    """
    from collections import Counter
    classes: Counter = Counter()
    preds: Counter = Counter()
    sources: Counter = Counter()
    n_nodes = n_edges = 0
    for line in (release / "nodes.jsonl").read_text().splitlines():
        if line.strip():
            classes[json.loads(line)["class"]] += 1
            n_nodes += 1
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        preds[e["predicate"]] += 1
        sources[e.get("primary_knowledge_source") or "(none)"] += 1
        n_edges += 1
    return n_nodes, n_edges, classes, preds, sources


def block(release: Path) -> list[str]:
    n_nodes, n_edges, classes, preds, sources = _assembled(release)
    held = set(load_schema().held_out_predicates())
    fp = release_fingerprint(release)
    out = [BEGIN, ""]
    out.append(f"    release          {release.name}")
    out.append(f"    nodes_sha        {fp['nodes_sha']}")
    out.append(f"    edges_sha        {fp['edges_sha']}")
    out.append("")
    out.append(f"    {n_nodes:,} nodes")
    out.append(f"    {n_edges:,} edges")
    out.append("")
    out.append("    NODES BY CLASS")
    for k, v in classes.most_common():
        out.append(f"      {k:<24} {v:>12,}")
    out.append("")
    out.append("    EDGES BY PREDICATE   (L = held out of every metapath)")
    for k, v in preds.most_common():
        flag = " L" if k in held else "  "
        out.append(f"      {k:<34}{flag} {v:>12,}")
    out.append("")
    out.append("    EDGES BY SOURCE")
    for k, v in sources.most_common():
        out.append(f"      {k:<34} {v:>12,}")
    man = release / "MANIFEST.json"
    if man.exists():
        m = json.loads(man.read_text())
        layers = m.get("layers") or {}
        out.append("")
        out.append(f"    {len(layers)} layer(s) assembled"
                   + (f", hub cut degree {m['hub_cut_degree']}"
                      if "hub_cut_degree" in m else ""))
        for name in sorted(layers):
            out.append(f"      {name}")
    out.append("")
    out.append(END)
    return out


def splice(text: str, lines: list[str]) -> str:
    assert BEGIN in text and END in text, \
        "the report lost its generated-block markers"
    head = text.split(BEGIN, 1)[0]
    tail = text.split(END, 1)[1]
    return head + "\n".join(lines) + tail


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release", type=Path,
                    default=ROOT / "data/releases/v0.1")
    ap.add_argument("--report", type=Path, default=REPORT)
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if the block is stale, write nothing")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"no assembled release at {args.release}", file=sys.stderr)
        return 2
    text = args.report.read_text()
    want = splice(text, block(args.release))
    if args.check:
        if want != text:
            print("docs/kg_audit_report.txt census is STALE for "
                  f"{args.release}. Run: python scripts/write_kg_report.py",
                  file=sys.stderr)
            return 1
        print(f"census is current for {args.release}")
        return 0
    args.report.write_text(want)
    n = len(block(args.release))
    print(f"wrote {n} line(s) of census into {args.report} "
          f"from {args.release}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
