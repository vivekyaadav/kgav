#!/usr/bin/env python3
"""Generate docs/figures/meta.dot from the assembled release.

    make metagraph-figure         # regenerate the .dot and re-render
    make metagraph-figure-check   # fail if it is stale, write nothing

WHY THIS FILE EXISTS. meta.dot was hand-written on 2026-11-03 and nothing
ever checked it. By 2026-11-20 it carried THIRTEEN stale figures and
omitted TWO predicates:

    nodes                     78,789  ->  78,649
    unique directed pairs    390,684  ->  393,927
    density                    6.3e-5 ->  6.4e-5
    dense float32              24.8GB ->  24.7GB
    CSR                        3.4MB  ->  3.5MB
    SmallMolecule             25,383  ->  25,243
    SmallMolecule mean degree      2.2 ->  2.3
    OrganismTaxon mean degree  2,033  ->  2,247.9
    SARS-CoV-2 degree         17,493  ->  17,342
    HAS_ANTIVIRAL_ACTIVITY    24,496  ->  27,928
    TARGETS                   17,982  ->  17,635
    TARGETS reach        4,130/1,946  ->  3,977/1,944
    INHIBITS                   3,943  ->  3,874
    MEASURED_INACTIVE          2,421  ->  2,420
    CHEMICALLY_SIMILAR_TO      3,384  ->  3,035

    missing entirely: IN_TRIAL_FOR (568) and TREATS (9)

THE THREE DEGREE FIGURES ABOVE ARE AGAINST THE PAIR CONVENTION. The
first version of this generator counted node degree in edges.jsonl
LINES while every arrow counted distinct (subject, object) pairs, so one
picture used two conventions and "mean deg 2.4" was not comparable to
the nnz beside it. The gap is the 2,296 lines two sources assert for one
edge. Both halves now count distinct (subject, predicate, object), which
moves SmallMolecule 2.4 -> 2.3, OrganismTaxon 2,383.7 -> 2,247.9 and
SARS-CoV-2 19,498 -> 17,342 -- the last landing near the hand-written
17,493, which suggests the original author counted pairs too.

The two missing ones matter most. They are the clinical axis, added on
2026-11-17, and they are one of the four findings v1 now publishes --
absent from the figure that is supposed to show what the graph holds.
A hand-drawn figure cannot acquire a predicate.

THE TOPOLOGY IS DERIVED, NOT JUST THE NUMBERS. Every arrow is a
(subject class, predicate, object class) block actually present in
edges.jsonl, so a predicate that appears in a later release draws itself
and one that disappears stops being drawn. That is the part a hand-
written .dot could never do, and it is why this is a generator rather
than a one-off correction.

WHAT IS NOT CHECKED. The rendered PNG/SVG/PDF cannot be byte-compared --
graphviz versions differ in their output - so --check verifies the .dot
only, and the make target re-renders. A rendered image older than the
.dot is therefore possible; `make metagraph-figure` is one command and
regenerates both.
"""
from __future__ import annotations

import argparse
import collections
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.pipeline import release_dir
from kgav.provenance import release_fingerprint
from kgav.schema import load_schema, print_capped

DOT = ROOT / "docs" / "figures" / "meta.dot"

#: Short node ids in the .dot, and the fill each class gets. A class not
#: listed here still draws -- with a generated id and the default fill --
#: so a new node class appears rather than being silently dropped.
STYLE = {
    "SmallMolecule": ("SM", "#fbe5d6"),
    "Protein": ("P", "#deebf7"),
    "Gene": ("G", "#e2f0d9"),
    "Pathway": ("PW", "#fff2cc"),
    "OrganismTaxon": ("T", "#f8cbad"),
    "TargetClass": ("TC", "#eeeeee"),
    "Disease": ("D", "#eeeeee"),
}
#: Display names where the class name alone understates what is in it.
DISPLAY = {"Pathway": "Pathway (GO+Reactome)", "Protein": "Protein"}

#: The one taxon named in the figure, because its degree is the reason
#: the OrganismTaxon mean is not a useful summary on its own.
HEADLINE_TAXON = "NCBITaxon:2697049"

#: Diff lines shown when --check finds the file stale. The remainder is
#: counted rather than dropped; print_capped is what guarantees that.
DIFF_LINES = 60

#: Drawn in the held-out colour with a label saying so. Read from the
#: schema rather than listed, so a predicate that becomes held out is
#: recoloured without anyone editing this file.
HELD_OUT_COLOUR = "#c00000"


class Census:
    """Everything the figure states, counted once from the release."""

    def __init__(self, release: Path) -> None:
        self.release = release
        self.node_class: dict[str, str] = {}
        self.is_viral: set[str] = set()
        for line in (release / "nodes.jsonl").read_text().splitlines():
            if not line.strip():
                continue
            n = json.loads(line)
            self.node_class[n["id"]] = n.get("class")
            if n.get("properties", {}).get("is_viral"):
                self.is_viral.add(n["id"])

        # nnz per block is DISTINCT (subject, object): the figure calls
        # each arrow a sparse block matrix, and a matrix entry is a pair,
        # not a row of a file. edges.jsonl can hold two lines for one
        # pair under one predicate when two sources assert it.
        self.block: dict[tuple, set] = collections.defaultdict(set)
        self.degree: collections.Counter = collections.Counter()
        self.pairs: set = set()
        # KEYED BY BLOCK, NOT BY PREDICATE. These were {predicate: set},
        # so a predicate linking two class pairs -- which the docstring
        # above promises will draw itself as two arrows -- would print
        # the same pooled counts on both, mixing classes and able to
        # exceed the class size it names. No predicate does that on this
        # release, which is exactly why it had to be fixed by reading
        # rather than by a failure.
        self.subjects: dict[tuple, set] = collections.defaultdict(set)
        self.objects: dict[tuple, set] = collections.defaultdict(set)
        self.viral_objects: dict[tuple, set] = collections.defaultdict(set)
        self.lines = 0
        self._counted: set = set()
        for line in (release / "edges.jsonl").read_text().splitlines():
            if not line.strip():
                continue
            e = json.loads(line)
            s, o, p = e["subject"], e["object"], e["predicate"]
            self.lines += 1
            key = (self.node_class.get(s), p, self.node_class.get(o))
            self.block[key].add((s, o))
            # DEGREE COUNTS THE SAME THING THE ARROWS DO: distinct
            # (subject, object) pairs under one predicate. It counted
            # edges.jsonl LINES, so the node labels and the arrow labels
            # used two conventions in one picture -- a reader comparing
            # "mean deg 2.4" against the nnz beside it was comparing
            # different quantities, and the gap is the 2,296 duplicate
            # lines two sources assert.
            if (s, p, o) not in self._counted:
                self._counted.add((s, p, o))
                self.degree[s] += 1
                self.degree[o] += 1
            self.pairs.add((s, o))
            self.subjects[key].add(s)
            self.objects[key].add(o)
            if o in self.is_viral:
                self.viral_objects[key].add(o)

        self.by_class: dict[str, list[int]] = collections.defaultdict(list)
        for nid, cls in self.node_class.items():
            self.by_class[cls].append(self.degree[nid])

    # ------------------------------------------------------------ totals
    @property
    def n_nodes(self) -> int:
        return len(self.node_class)

    @property
    def n_pairs(self) -> int:
        return len(self.pairs)

    @property
    def density(self) -> float:
        return self.n_pairs / (self.n_nodes ** 2)

    @property
    def dense_gb(self) -> float:
        """float32 adjacency, if anyone materialised it."""
        return self.n_nodes ** 2 * 4 / 1e9

    @property
    def csr_mb(self) -> float:
        """indices + values at 4 bytes each, plus the row pointer."""
        return (self.n_pairs * 8 + (self.n_nodes + 1) * 4) / 1e6

    def count(self, cls: str) -> int:
        return len(self.by_class.get(cls, []))

    def mean_degree(self, cls: str) -> float:
        v = self.by_class.get(cls) or [0]
        return sum(v) / len(v)

    def isolated(self, cls: str) -> int:
        return sum(1 for d in self.by_class.get(cls, []) if d == 0)


def _fmt_density(nnz: int, rows: int, cols: int) -> str:
    if not rows or not cols:
        return ""
    return f"{nnz / (rows * cols):.1e}"


#: What to call the members of each class in an arrow's parenthetical.
#: NAMED PER CLASS, because the first version of this keyed the phrasing
#: on the OBJECT class alone and so labelled CAUSES -- which runs
#: OrganismTaxon -> Disease -- as "16 compounds -> 10 diseases". A
#: generated figure that invents the wrong noun is worse than a stale one:
#: the number is right, so nothing looks wrong.
MEMBERS = {
    "SmallMolecule": ("compound", "compounds"),
    "Protein": ("protein", "proteins"),
    "Gene": ("gene", "genes"),
    "Pathway": ("pathway", "pathways"),
    "OrganismTaxon": ("taxon", "taxa"),
    "TargetClass": ("class", "classes"),
    "Disease": ("disease", "diseases"),
}


def _members(cls: str, n: int) -> str:
    one, many = MEMBERS.get(cls, (cls, cls))
    return f"{n:,} {one if n == 1 else many}"


def annotation(c: Census, pred: str, sc: str, oc: str) -> str:
    """The parenthetical under an arrow, computed rather than remembered.

    It answers "reaches how much of what it points at", which a block's
    nnz does not: 17,635 TARGETS edges could be 17,635 drugs or 400.

    PER BLOCK. The counts belong to this arrow's (subject class,
    predicate, object class), not to the predicate pooled across every
    class pair it might link.
    """
    key = (sc, pred, oc)
    subs, objs = len(c.subjects[key]), len(c.objects[key])
    # Into a small class, the useful figure is coverage of it.
    if c.count(oc) <= 100:
        reach = f"{objs} of {c.count(oc)} {MEMBERS.get(oc, (oc, oc))[1]}"
        return f"({_members(sc, subs)} -> {reach})"
    note = f"({_members(sc, subs)} -> {_members(oc, objs)})"
    viral = len(c.viral_objects.get(key, ()))
    if viral and oc == "Protein":
        note = note[:-1] + (", all viral)" if viral == objs
                            else f", {viral} viral)")
    return note


def render_dot(c: Census) -> str:
    """The .dot text. Deterministic: blocks sorted by descending nnz."""
    fp = release_fingerprint(c.release)
    schema = load_schema()
    held_out = set(schema.held_out_predicates())

    lines = [
        "// GENERATED FILE. Do not edit.",
        "//   make metagraph-figure         regenerate and re-render",
        "//   make metagraph-figure-check   fail if stale, write nothing",
        "//",
        "// Every figure below is counted from the release named in the",
        "// label. The predecessor of this file was hand-written and had",
        "// drifted by thirteen figures and two whole predicates.",
        f"// release {fp['id']}  nodes_sha {fp['nodes_sha']}  "
        f"edges_sha {fp['edges_sha']}",
        "digraph meta {",
    ]

    label = (
        f"kgav {fp['id']} metagraph — each arrow is one sparse block "
        f"matrix: rows×cols, nnz, density\\n"
        f"whole graph: A ∈ {{0,1}}^({c.n_nodes:,} × {c.n_nodes:,}), "
        f"{c.n_pairs:,} unique directed pairs, density {c.density:.1e} "
        f"(dense float32 = {c.dense_gb:.1f} GB; CSR ≈ {c.csr_mb:.1f} MB)\\n"
        f"nnz and every degree count DISTINCT (subject, predicate, object); "
        f"edges.jsonl holds {c.lines:,} lines\\n"
        f"edges_sha {fp['edges_sha']}"
    )
    lines += [
        f'  graph [fontname="Helvetica", label="{label}", labelloc=t, '
        f'fontsize=14, pad=0.5, nodesep=1.4, ranksep=1.3];',
        '  node [fontname="Helvetica", shape=ellipse, style=filled, '
        'fontsize=11];',
        '  edge [fontname="Helvetica", fontsize=9];',
    ]

    # ---- nodes, in descending population so the picture reads top-down
    ids: dict[str, str] = {}
    classes = sorted((cl for cl in c.by_class if cl),
                     key=lambda cl: -c.count(cl))
    for i, cls in enumerate(classes):
        short, fill = STYLE.get(cls, (f"N{i}", "#eeeeee"))
        ids[cls] = short
        bits = [DISPLAY.get(cls, cls), f"{c.count(cls):,}"]
        if c.count(cls) > 1:
            bits.append(f"mean deg {c.mean_degree(cls):,.1f}")
        if c.isolated(cls):
            bits.append(f"{c.isolated(cls):,} isolated")
        if cls == "OrganismTaxon" and HEADLINE_TAXON in c.degree:
            bits[-1] = (bits[-1] + f" (SARS-CoV-2: "
                                   f"{c.degree[HEADLINE_TAXON]:,})")
        lines.append(f'  {short:<3} [label="' + "\\n".join(bits)
                     + f'", fillcolor="{fill}"];')

    # ---- one arrow per block present in the release
    for (sc, pred, oc), pset in sorted(c.block.items(),
                                       key=lambda kv: (-len(kv[1]), kv[0][1])):
        if sc not in ids or oc not in ids:
            continue
        nnz = len(pset)
        rows, cols = c.count(sc), c.count(oc)
        # SELF-LOOPS CARRY NO SHAPE. Three Protein -> Protein blocks are
        # drawn at the same place, and graphviz stacks their labels: with
        # "20,888²" on each they overlapped and the figure was unreadable
        # where it is densest. The shape of a self block is the class
        # count, which the node already prints.
        shape = "" if sc == oc else f"{rows:,}×{cols:,} "
        label = f"{pred} {shape}· {nnz:,} · {_fmt_density(nnz, rows, cols)}"
        note = annotation(c, pred, sc, oc)
        if note:
            label += "\\n" + note
        attrs = []
        if pred in held_out:
            label = f"{pred} (HELD OUT) {shape}· {nnz:,} · " \
                    f"{_fmt_density(nnz, rows, cols)}"
            if note:
                label += "\\n" + note
            attrs += [f'color="{HELD_OUT_COLOUR}"',
                      f'fontcolor="{HELD_OUT_COLOUR}"']
            attrs.append("style=dashed" if nnz < 10_000 else "penwidth=2.5")
        elif nnz >= 100_000:
            attrs.append("penwidth=3")
        elif nnz >= 10_000:
            attrs.append("penwidth=2")
        tail = (", " + ", ".join(attrs)) if attrs else ""
        lines.append(f'  {ids[sc]:<3} -> {ids[oc]:<3} '
                     f'[label="{label}"{tail}];')

    lines.append("}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--out", type=Path, default=DOT)
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if the file on disk is not what this "
                         "would write; write nothing")
    ap.add_argument("--render", action="store_true",
                    help="also run graphviz for the PNG and SVG")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"release not found at {args.release}", file=sys.stderr)
        return 2

    text = render_dot(Census(args.release))

    # ONE STREAM. These three messages used to go to stderr while
    # print_capped writes to stdout, so the footer appeared ABOVE the
    # diff whenever the two were read together -- "regenerates it" before
    # the reader had seen what moved. The sibling generators
    # (write_scorer_table.py, write_kg_report.py) report staleness on
    # stdout and their tests read it there.
    if args.check:
        if not args.out.exists():
            print(f"{args.out} does not exist; run `make metagraph-figure`")
            return 1
        current = args.out.read_text()
        if current != text:
            import difflib
            print(f"{args.out.name} is STALE. It does not describe "
                  f"{args.release.name}:")
            diff = list(difflib.unified_diff(
                current.splitlines(), text.splitlines(),
                "on disk", "from the release", lineterm="", n=0))
            # THROUGH print_capped, not a hand-capped slice. Caught by
            # tests/test_no_duplicate_logic.py's AST sweep, which exists
            # because fourteen loops in this repository printed a
            # constant-capped head and exactly one said what it had left
            # out. A staleness report that shows ten diff lines and
            # stops reads the same whether ten moved or four hundred.
            print_capped(diff, DIFF_LINES, indent="  ", unit="more diff line(s)")
            print("\n  `make metagraph-figure` regenerates it.")
            return 1
        print(f"{args.out.name} is current for {args.release.name}")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    print(f"wrote {args.out}")

    if args.render:
        stem = args.out.parent / "kgav_metagraph"
        for fmt, dest in (("png", stem.with_suffix(".png")),
                          ("svg", stem.with_suffix(".svg"))):
            r = subprocess.run(["dot", f"-T{fmt}", str(args.out),
                                "-o", str(dest)], capture_output=True,
                               text=True)
            if r.returncode != 0:
                print(f"graphviz failed for {fmt}: {r.stderr}",
                      file=sys.stderr)
                return 1
            print(f"rendered {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
