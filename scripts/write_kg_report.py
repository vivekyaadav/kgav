#!/usr/bin/env python3
"""Generate docs/kg_release_census.txt from the assembled release.

WHY GENERATED RATHER THAN WRITTEN. The audit report first carried no
node or edge totals at all, on the reasoning that a hardcoded count goes
stale silently -- which is true, and the wrong conclusion. An audit
report IS a stamped snapshot of one release; omitting its size makes it
useless to a reader and protects nothing. The protection is to generate
the numbers and stamp them with the content hashes they were computed
from, so a stale figure is detectable rather than merely regrettable.

WHY ITS OWN FILE. The census lived inside docs/kg_audit_report.txt
between two markers, and the generated block and the hand-written prose
around it have different lifecycles: one changes on every rebuild, the
other on every correction. Sharing a file guaranteed they collided --
twice, a regenerated census left uncommitted made the next prose patch
fail to apply, which is friction with no upside. This file is generated
in full, so it has no markers to lose and no prose to protect.

    make kg-report          # regenerate
    make kg-report-check    # fail if it is stale, write nothing
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

CENSUS = ROOT / "docs" / "kg_release_census.txt"

#: Printed at the top of the generated file. A reader who opens it
#: without context has to learn it is not hand-maintained before they
#: think about editing it.
HEADER = (
    "=" * 78,
    "kgav — RELEASE CENSUS",
    "",
    "GENERATED FILE. Do not edit: `make kg-report` overwrites it, and",
    "`make kg-report-check` fails if it does not describe the assembled",
    "release. The prose that interprets these figures is section 4 of",
    "docs/kg_audit_report.txt.",
    "=" * 78,
    "",
)


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
    out: list[str] = []
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
    return out


def render(release: Path) -> str:
    """The whole file, header included."""
    return "\n".join(list(HEADER) + block(release)) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release", type=Path,
                    default=ROOT / "data/releases/v0.1")
    ap.add_argument("--out", type=Path, default=CENSUS)
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if the file is stale, write nothing")
    ap.add_argument("--placeholder", action="store_true",
                    help="with no assembled release, write a file that "
                         "says so instead of failing")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        if args.placeholder:
            # WRITTEN BY THE GENERATOR, not typed into the repository.
            # A fresh clone has no assembled release, and the committed
            # file has to say that rather than describe a graph nobody
            # can check. --check still fails on it, correctly: a
            # placeholder is stale until someone regenerates.
            args.out.write_text("\n".join(list(HEADER) + [
                "    NO ASSEMBLED RELEASE was on disk when this file was",
                "    written, so it describes nothing. Run `make kg-report`",
                "    against an assembled release to fill it in.",
                "",
            ]) + "\n")
            print(f"wrote a placeholder into {args.out}")
            return 0
        print(f"no assembled release at {args.release}", file=sys.stderr)
        return 2
    want = render(args.release)
    if args.check:
        # MISSING IS STALE, not absent-and-therefore-fine. A check that
        # passed when the file did not exist would certify the one state
        # it is least able to vouch for.
        have = args.out.read_text() if args.out.exists() else None
        if have != want:
            why = "is missing" if have is None else "is STALE"
            print(f"{args.out} {why} for {args.release}. "
                  f"Run: make kg-report", file=sys.stderr)
            return 1
        print(f"census is current for {args.release}")
        return 0
    args.out.write_text(want)
    print(f"wrote {want.count(chr(10))} line(s) into {args.out} "
          f"from {args.release}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
