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

import textwrap

from kgav.provenance import release_fingerprint
from kgav.schema import load_schema

REPORT = ROOT / "docs" / "kg_audit_report.txt"
BEGIN = "<<<BEGIN GENERATED: release census — python scripts/write_kg_report.py>>>"
END = "<<<END GENERATED>>>"


def _wrap(text: str, indent: str = "      ") -> list[str]:
    """Wrapped to the report's 80 columns.

    A summary line naming eight absent classes runs past the margin,
    and the rest of this file is hard-wrapped, so an unwrapped line is
    the one thing in a generated block a reader notices before the
    content.
    """
    return textwrap.wrap(text, width=78, initial_indent=indent,
                         subsequent_indent=indent + "  ") or [indent]


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
    # EVERY DECLARED CLASS, INCLUDING THE EMPTY ONES. Listing only what
    # is present reports presence and stays silent on absence, which is
    # the thing this project refuses everywhere else: "a source this
    # ingest cannot read is not a source with no data". v0.1 holds 7 of
    # 11 declared node classes and 14 of 24 predicates, and the first
    # version of this census showed neither fact.
    schema = load_schema()
    out.append("    NODES BY CLASS   (0 = declared in the schema, "
               "absent from this release)")
    empty_c = []
    for k in sorted(schema.node_classes, key=lambda n: (-classes[n], n)):
        out.append(f"      {k:<24} {classes[k]:>12,}")
        if not classes[k]:
            empty_c.append(k)
    total_c = sum(classes.values())
    if total_c != n_nodes:
        out.append(f"      !! class counts sum to {total_c:,}, "
                   f"not {n_nodes:,}")
    if empty_c:
        out.extend(_wrap(f"{len(empty_c)} of {len(schema.node_classes)} "
                         f"declared classes are empty: "
                         + ", ".join(empty_c)))
    out.append("")
    held = set(schema.held_out_predicates())
    labels = set(schema.evaluation_label_predicates())
    out.append("    EDGES BY PREDICATE   (L = held out of every metapath, "
               "0 = declared, absent)")
    empty_p = []
    for e in sorted(schema.edge_classes,
                    key=lambda e: (-preds[e.predicate], e.predicate)):
        k = e.predicate
        flag = " L" if k in held else "  "
        out.append(f"      {k:<34}{flag} {preds[k]:>12,}")
        if not preds[k]:
            empty_p.append(k)
    total_p = sum(preds.values())
    if total_p != n_edges:
        out.append(f"      !! predicate counts sum to {total_p:,}, "
                   f"not {n_edges:,}")
    if empty_p:
        out.append(f"      {len(empty_p)} of {len(schema.edge_classes)} "
                   f"declared predicates are empty:")
        for k in empty_p:
            why = "  <- EVALUATION LABEL" if k in labels else ""
            out.append(f"        {k}{why}")
        absent_labels = sorted(k for k in empty_p if k in labels)
        if absent_labels:
            kept = sorted(labels - set(absent_labels))
            tail = (f"so every evaluation rests on {', '.join(kept)} alone"
                    if kept else
                    "so this release carries NO evaluation label at all "
                    "and nothing in it can be scored")
            out.extend(_wrap(
                f"{len(absent_labels)} of {len(labels)} evaluation label "
                f"predicate(s) carry no edges, {tail}."))
    out.append("")
    # PER CURATOR, NOT PER STRING. primary_knowledge_source is a "|"
    # separated list once an edge has been merged, so counting the raw
    # values reports "infores:intact|infores:virhostnet" as a source of
    # its own. v0.1 showed thirty such permutations of six real
    # curators, which made the real distribution unreadable in the one
    # place it is supposed to be legible.
    from collections import Counter as _C
    per_curator: _C = _C()
    multi = 0
    repeated = 0
    for value, n in sources.items():
        parts = [q for q in str(value).split("|") if q]
        if len(parts) > 1:
            multi += n
        if len(parts) != len(set(parts)):
            repeated += n
        for q in set(parts):
            per_curator[q] += n
    out.append("    EDGES BY CURATOR   (an edge with several curators "
               "counts once for each)")
    for k, v in per_curator.most_common():
        out.append(f"      {k:<34} {v:>12,}")
    if multi:
        out.extend(_wrap(f"{multi:,} edge(s) cite more than one curator, so "
                         f"these counts sum to more than the edge total."))
    else:
        out.extend(_wrap("every edge cites exactly one curator."))
    if repeated:
        out.extend(_wrap(
            f"!! {repeated:,} edge(s) carry a REPEATED curator in one "
            f"value. merge_edges unioned source STRINGS rather than the "
            f"curators inside them, so a fact with several witnesses "
            f"accumulated duplicates. Fixed in the assembler; this "
            f"release was assembled before the fix and reassembling "
            f"rewrites the field."))
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
