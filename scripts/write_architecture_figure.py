#!/usr/bin/env python3
"""Generate the architecture figures from the code that defines the pipeline.

    make arch-figure         # regenerate both .dot files and re-render
    make arch-figure-check   # fail if either .dot is stale, write nothing

WHY THIS FILE EXISTS. `arch.dot` and `arch_a4.dot` were HAND-DRAWN, and
the Makefile defended that on the grounds that "the grouping, the
ordering and the choice of what to draw are editorial". The grouping is.
Nothing else in them was:

  * the layer list, the ingest script behind each layer, and which
    layers run before and after the first merge are
    `kgav.pipeline.LAYERS` -- a tuple with a `phase` field;
  * the raw inputs each layer needs are that layer's own globs;
  * the merge order is `kgav.assemble.LAYER_PRECEDENCE`;
  * every count is in `MANIFEST.json` or is a directory listing.

Drawn by hand, all of that went stale: twelve wrong items sat in those
two files for seventeen days, including three defects drawn in red that
were fixed or false and a box describing a directory that had been
deleted. tests/test_architecture_figure.py was written after that and
caught the figures up, which is a check on a hand-drawn file rather than
a reason for one.

WHAT IS STILL EDITORIAL, and is therefore declared here rather than
derived: the four STAGES, the human name each raw directory gets, and
the column each layer sits in. `_check_editorial_tables_are_complete`
fails when a layer, a raw directory or a node class has no entry, so a
thirteenth layer cannot be silently left out of the picture -- which is
the one failure mode a generated figure still has.

TWO FILES, ONE TOPOLOGY. `arch.dot` is wide (for a screen) and
`arch_a4.dot` is portrait and sized to print on one A4 sheet. They are
rendered from one description, so they cannot disagree, which the two
hand-drawn ones did: `arch.dot` carried seven red defect boxes and
`arch_a4.dot` carried eight with three marked closed.

THE SIZE FLAGS ARE LOAD-BEARING. graphviz adds about 0.2in of margin on
top of -Gsize, so asking for A4 exactly produces a page 11.89in tall
that does not fit A4. 7.9x11.3 lands at 11.68 against A4's 11.69.
300 dpi because a careless re-render once dropped the PNG from
2480x3508 to 930x1726.

WHAT IS NOT CHECKED. The rendered PNG/SVG/PDF cannot be byte-compared --
graphviz versions differ in their output -- so --check verifies the .dot
files only and `make arch-figure` re-renders both.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.agent import guards                              # noqa: E402
from kgav.assemble import LAYER_PRECEDENCE                 # noqa: E402
from kgav.pipeline import LAYERS, release_dir              # noqa: E402
from kgav.provenance import release_fingerprint            # noqa: E402
from kgav.schema import load_schema                        # noqa: E402

FIGURES = ROOT / "doc2" / "figures"
WIDE = FIGURES / "arch.dot"
A4 = FIGURES / "arch_a4.dot"

#: Human name for each raw input directory, keyed on the directory under
#: `data/raw/`. EDITORIAL -- a database's own name is not in its path --
#: and complete, which `_check_editorial_tables_are_complete` enforces:
#: a layer reading a directory absent from here fails the generator
#: rather than drawing a box labelled with a path.
RAW_SOURCE = {
    "proteomes": "UniProt\nreference proteomes",
    "host": "UniProt · STRING\nReactome · GO",
    "vhppi": "VirHostNet\nIntAct · MINT · DIP",
    "orcs": "BioGRID ORCS\nCRISPR screens",
    "curated_screens": "curated screens\n(hand-entered)",
    "chembl": "ChEMBL 37\nSQLite dump",
    "ncats": "NCATS OpenData\nSARS-CoV-2 plates",
    "folds": "UniProt · InterPro\ncatalytic chains",
    "gtop": "Guide to\nPHARMACOLOGY",
    "mondo": "MONDO\ndisease ontology",
}

#: What each layer contributes, in the words a reader needs rather than
#: the predicate names (those are in the metagraph figure). EDITORIAL and
#: complete -- a layer missing from here fails the generator.
LAYER_GIVES = {
    "spine": "viruses, viral proteins,\nthe taxonomy spine",
    "host": "human proteins, genes,\npathways, host–host PPI",
    "vhppi": "virus–host PPI",
    "orcs": "host dependency factors",
    "curated": "host factors from\nhand-entered screens",
    "chembl": "compound activity,\nviral-protein inhibition",
    "hosttargets": "compound → human target",
    "ncats": "plate activity and\ncytotoxicity",
    "folds": "protein classes,\nfold similarity",
    "selectivity": "CC50/EC50 selectivity",
    "gtop": "drug → target,\nwith direction",
    "similarity": "ECFP4 Tanimoto\nneighbours",
    "indications": "trials and approvals,\ndiseases",
}

#: Fill colours per stage. The palette is deliberately pale: these are
#: printed on A4 in black and white as often as they are read on screen.
FILL = {
    "raw": "#e8eef8",
    "layer1": "#e9f5e1",
    "layer2": "#fdebd9",
    "merge": "#fff2cc",
    "release": "#ffe599",
    "eval": "#f0eaf9",
    "results": "#e2d8f2",
    "agent": "#e3f4f4",
    "doc": "#f2f2f2",
}


def _manifest() -> dict:
    return json.loads((release_dir(ROOT) / "MANIFEST.json").read_text())


def _raw_dirs(layer) -> list[str]:
    """The `data/raw/` subdirectories one layer reads, in order.

    A layer declares its inputs three ways -- `files`, `globs` and
    `any_globs` (alternatives, any one of which satisfies it) -- and a
    figure that read only `globs` would draw `ncats` and `indications`
    as sourceless.
    """
    pats: list[str] = []
    pats += list(layer.files)
    pats += list(layer.globs)
    for alts in layer.any_globs:
        pats += list(alts)
    out: list[str] = []
    for p in pats:
        parts = Path(p).parts
        if len(parts) >= 3 and parts[0] == "data" and parts[1] == "raw":
            if parts[2] not in out:
                out.append(parts[2])
    return out


def _check_editorial_tables_are_complete() -> None:
    """The one failure mode a generated figure still has.

    A derived topology cannot forget a layer, but a hand-written LABEL
    for that layer can be absent -- and the obvious handling, falling
    back to the layer's name, draws a box that looks finished. So this
    raises instead, naming what to add.
    """
    missing = [L.name for L in LAYERS if L.name not in LAYER_GIVES]
    if missing:
        raise SystemExit(
            f"LAYER_GIVES has no entry for {missing}; add one (what the "
            f"layer contributes, in a reader's words) and re-run")
    dirs = {d for L in LAYERS for d in _raw_dirs(L)}
    unnamed = sorted(dirs - set(RAW_SOURCE))
    if unnamed:
        raise SystemExit(
            f"RAW_SOURCE has no name for data/raw/{unnamed}; add the "
            f"database's own name and re-run")
    extra = sorted(set(LAYER_GIVES) - {L.name for L in LAYERS})
    if extra:
        raise SystemExit(
            f"LAYER_GIVES names {extra}, which kgav.pipeline.LAYERS does "
            f"not; a layer was removed and the figure still draws it")


# ------------------------------------------------------------- the figure
def _esc(s: str) -> str:
    return s.replace('"', '\\"').replace("\n", "\\n")


def _box(nid: str, label: str, fill: str, **attrs) -> str:
    bits = [f'label="{_esc(label)}"', f'fillcolor="{fill}"']
    bits += [f'{k}={v}' for k, v in attrs.items()]
    return f'  {nid} [' + ", ".join(bits) + '];'


def topology() -> dict:
    """Everything the two figures draw, derived once.

    Returned as plain data so the wide and the A4 renderer cannot
    disagree about any of it -- which the two hand-drawn files did.
    """
    _check_editorial_tables_are_complete()
    m = _manifest()
    schema = load_schema()
    built = sorted(m["layers"])
    phase1 = [L for L in LAYERS if L.phase == 1]
    phase2 = [L for L in LAYERS if L.phase == 2]
    results = sorted((ROOT / "data" / "results").glob("*.json"))
    doc2 = sorted((ROOT / "doc2").glob("*.md"))

    def layer_node(L):
        return {
            "id": f"L_{L.name}",
            "name": L.name,
            "script": L.script,
            "gives": LAYER_GIVES[L.name],
            "sources": [RAW_SOURCE[d] for d in _raw_dirs(L)],
            "raw_dirs": _raw_dirs(L),
            "built": L.name in built,
            "required": L.required,
        }

    return {
        "release": m["release"],
        "schema_version": m["schema_version"],
        "n_nodes": m["nodes"],
        "n_edges": m["edges"],
        "n_node_classes": len(m["nodes_by_class"]),
        "n_predicates": len(m["edges_by_predicate"]),
        "n_predicates_declared": len(schema.edge_classes),
        "n_merged": m["merges"].get("edge_merged", 0),
        "fingerprint": release_fingerprint(release_dir(ROOT)),
        "phase1": [layer_node(L) for L in phase1],
        "phase2": [layer_node(L) for L in phase2],
        "precedence": list(LAYER_PRECEDENCE),
        "n_layers_declared": len(LAYERS),
        "n_layers_built": len(built),
        "absent": [L.name for L in LAYERS if L.name not in built],
        "n_results": len(results),
        "n_doc2": len(doc2),
        "n_intents": len(guards.ROUTED_KINDS),
        "n_caveats": len(guards.CAVEAT_LABELS),
        "n_scripts": len(sorted((ROOT / "scripts").glob("*.py"))),
    }


def _header(t: dict, wide: bool) -> str:
    title = "kgav — how the graph is built, measured and queried"
    sub = (f"release {t['release']} · schema {t['schema_version']} · "
           f"{t['n_nodes']:,} nodes · {t['n_edges']:,} edges · "
           f"{t['n_node_classes']} node classes · {t['n_predicates']} predicates")
    note = (f"GENERATED from kgav.pipeline.LAYERS, kgav.assemble."
            f"LAYER_PRECEDENCE and MANIFEST.json — nodes_sha "
            f"{t['fingerprint']['nodes_sha']}")
    size = '  size="16,9";\n' if wide else '  size="7.9,11.3"; margin=0.1;\n'
    rank = "LR" if wide else "TB"
    return (
        f'digraph kgav {{\n'
        f'  rankdir={rank};\n'
        f'{size}'
        f'  graph [label=<<font point-size="{18 if not wide else 20}">'
        f'<b>{title}</b></font><br/>'
        f'<font point-size="10">{sub}</font><br/>'
        f'<font point-size="9" color="#555555">{note}</font>>,\n'
        f'         labelloc=t, labeljust=c, fontname="Helvetica", '
        f'fontsize=11, bgcolor=white,\n'
        f'         nodesep={0.3 if not wide else 0.4}, '
        f'ranksep={0.35 if not wide else 0.6}, pad=0.3, '
        f'outputorder=edgesfirst];\n'
        f'  node [fontname="Helvetica", fontsize={8.5 if not wide else 10}, '
        f'shape=box, style="rounded,filled", fillcolor=white,\n'
        f'        penwidth=0.8, margin="0.07,0.045"];\n'
        f'  edge [fontname="Helvetica", fontsize={7.5 if not wide else 9}, '
        f'color="#666666", arrowsize=0.55, penwidth=0.9];\n'
    )


def _stage1(t: dict, wide: bool) -> list[str]:
    """Sources and the nine phase-1 ingest layers.

    One source box per raw directory, shared by the layers that read it:
    `chembl` is read by three layers and drawing it three times would
    imply three copies of a 30 GB dump.
    """
    out = ['  subgraph cluster_1 {',
           '    label=<<b>1 · INGEST</b>  <font point-size="8">'
           'one script per layer · scripts/ingest_*.py · '
           'each writes a standalone sub-release</font>>;',
           '    style="rounded"; color="#9db3d6"; fontname="Helvetica"; '
           'fontsize=10; margin=8;']
    cols = 2 if not wide else 9
    # `rank=same` CAN ONLY MERGE RANKS, NEVER SPLIT ONE. Nine layers
    # with no edge between them are one rank whatever is asked for, so
    # the first version of this drew a single row nine boxes wide and
    # `size` scaled the A4 page down to a 2in strip. Rows are forced by
    # an invisible edge from each row to the next.
    lrows = [t["phase1"][i:i + cols] for i in range(0, len(t["phase1"]), cols)]
    #: Each raw directory is drawn ONCE, in the row of the first layer
    #: that reads it -- `chembl` feeds three layers and three boxes
    #: would imply three copies of a 30 GB dump.
    placed: dict[str, int] = {}
    for r, row in enumerate(lrows):
        for L in row:
            for d in L["raw_dirs"]:
                placed.setdefault(d, r)
    srows: list[list[str]] = [[] for _ in lrows]
    for d, r in sorted(placed.items(), key=lambda kv: (kv[1], kv[0])):
        srows[r].append(f"S_{d}")
        out.append("  " + _box(f"S_{d}", RAW_SOURCE[d], FILL["raw"],
                               width=1.45, height=0.4))
    for row in lrows:
        for L in row:
            label = f"{L['name']}\n{L['gives']}"
            if not L["built"]:
                label += "\n(declared, not in this release)"
            out.append("  " + _box(
                L["id"], label,
                FILL["layer1"] if L["built"] else "#f4f4f4",
                width=1.55, height=0.55,
                style=('"rounded,filled"' if L["built"]
                       else '"rounded,filled,dashed"')))
    for r, (srow, lrow) in enumerate(zip(srows, lrows)):
        if srow:
            out.append("    {rank=same; " + " ".join(srow) + "}")
        out.append("    {rank=same; "
                   + " ".join(L["id"] for L in lrow) + "}")
        if r + 1 < len(lrows) and srows[r + 1]:
            out.append(f'    {lrow[0]["id"]} -> {srows[r + 1][0]} '
                       f'[style=invis, weight=20];')
    for L in t["phase1"]:
        for d in L["raw_dirs"]:
            out.append(f'    S_{d} -> {L["id"]} [arrowsize=0.45];')
    out.append('  }')
    return out


def _stage2(t: dict, wide: bool) -> list[str]:
    """Two merge passes with the release between them.

    The second pass is what makes this a cycle on paper and not in
    practice: the phase-2 layers READ the merged graph (similarity needs
    every compound's fingerprint, selectivity needs the activity edges)
    and are then merged back in. Drawn as two passes rather than one
    arrow, because "assemble runs twice" is the fact a reader needs.
    """
    out = ['  subgraph cluster_2 {',
           '    label=<<b>2 · ASSEMBLE</b>  <font point-size="8">'
           'scripts/assemble_release.py, run twice · '
           f'{t["n_merged"]:,} edges merged from more than one layer</font>>;',
           '    style="rounded"; color="#d6b656"; fontname="Helvetica"; '
           'fontsize=10; margin=8;']
    out.append("  " + _box("M1", "merge · pass 1\nphase-1 layers",
                           FILL["merge"], width=1.7, height=0.5))
    out.append("  " + _box(
        "REL", f"data/releases/{t['release']}\n{t['n_nodes']:,} nodes · "
        f"{t['n_edges']:,} edges\nnodes.jsonl · edges.jsonl · MANIFEST.json",
        FILL["release"], shape="cylinder", style="filled",
        width=2.4, height=0.7, penwidth=1.2))
    for L in t["phase2"]:
        out.append("  " + _box(L["id"], f"{L['name']}\n{L['gives']}",
                               FILL["layer2"], width=1.7, height=0.5))
    out.append("  " + _box("M2", "merge · pass 2\nwrites the final release",
                           FILL["merge"], width=1.7, height=0.5))
    per = len(t["phase2"]) if wide else 2
    prows = [t["phase2"][i:i + per]
             for i in range(0, len(t["phase2"]), per)]
    for r, row in enumerate(prows):
        out.append("    {rank=same; " + " ".join(L["id"] for L in row) + "}")
        if r + 1 < len(prows):
            out.append(f'    {row[0]["id"]} -> {prows[r + 1][0]["id"]} '
                       f'[style=invis, weight=20];')
    out.append('    M1 -> REL [label=" first pass"];')
    for L in t["phase2"]:
        out.append(f'    REL -> {L["id"]} [color="#c55a11", '
                   f'label=" reads the\\n merged graph"];'
                   if L is t["phase2"][0] else
                   f'    REL -> {L["id"]} [color="#c55a11"];')
        out.append(f'    {L["id"]} -> M2;')
    out.append('  }')
    out.append('  M2 -> REL [constraint=false, weight=0, '
               'color="#c55a11", style=dashed, label=" second pass"];')
    for L in t["phase1"]:
        out.append(f'  {L["id"]} -> M1 [arrowsize=0.45];')
    return out


def _stage3(t: dict, wide: bool) -> list[str]:
    out = ['  subgraph cluster_3 {',
           '    label=<<b>3 · MEASURE</b>  <font point-size="8">'
           'src/kgav/{labels,baselines,temporal,similarity,landscape}.py · '
           'every run stamped with the release it read</font>>;',
           '    style="rounded"; color="#b9a6d9"; fontname="Helvetica"; '
           'fontsize=10; margin=8;']
    for nid, label in (
            ("E_lab", "labels\nactive / measured-inactive,\nheld out of every scorer"),
            ("E_split", "splits\npublication-disjoint,\nnear-neighbour, temporal 2021"),
            ("E_score", "scorers\n8 metapath channels,\nstructure, degree baseline"),
            ("E_stat", "statistics\nAUC with interval,\npermutation null, correction")):
        out.append("  " + _box(nid, label, FILL["eval"],
                               width=1.75, height=0.6))
    out.append("  " + _box(
        "RES", f"data/results/\n{t['n_results']} stamped JSON artifacts",
        FILL["results"], shape="folder", width=1.9, height=0.5))
    out.append("  " + _box(
        "DOC", f"doc2/\n{t['n_doc2']} documents, generated\nfrom the artifacts",
        FILL["doc"], shape="note", width=1.9, height=0.5))
    out.append('    {rank=same; E_lab E_split}')
    out.append('    {rank=same; E_score E_stat}')
    out.append('    {rank=same; RES DOC}')
    out.append('    E_lab -> E_stat [label=" held-out\\n truth"];')
    out.append('    E_split -> E_score [label=" training\\n graph"];')
    out.append('    E_score -> E_stat;')
    out.append('    E_stat -> RES;')
    out.append('    RES -> DOC [style=dashed, label=" make doc2"];')
    out.append('  }')
    # THE VISIBLE EDGE AND THE RANKING EDGE ARE DIFFERENT EDGES. Drawing
    # `REL -> E_lab` with constraint=true would stretch the release box
    # across the whole measurement stage; drawing it with
    # constraint=false, which is what this did first, leaves stage 3
    # with NO ranking constraint at all and graphviz lays it BESIDE
    # stage 2 -- the A4 page came out 8.10 x 2.17in, four stages wide
    # and one deep, and `size` then scaled the whole thing to a strip.
    out.append('  REL -> E_lab [constraint=false, weight=0];')
    out.append('  REL -> E_split [constraint=false, weight=0];')
    out.append('  M2 -> E_split [style=invis, weight=20];')
    return out


def _stage4(t: dict, wide: bool) -> list[str]:
    out = ['  subgraph cluster_4 {',
           '    label=<<b>4 · ANSWER</b>  <font point-size="8">'
           'scripts/kgav_ask.py · src/kgav/agent/ · '
           f'{t["n_intents"]} question kinds · {t["n_caveats"]} caveat '
           f'labels</font>>;',
           '    style="rounded"; color="#8cc5c5"; fontname="Helvetica"; '
           'fontsize=10; margin=8;']
    for nid, label, fill in (
            ("A_q", "a question,\nin English", FILL["agent"]),
            ("A_o", "orchestrator\nintent · entities · guards", FILL["agent"]),
            ("A_t", "graph tools\nexplain · evidence · profile\n"
                    "pathways · shortlist", FILL["agent"]),
            ("A_m", "local model\nOllama, no network", "#c9eaea"),
            ("A_v", "verifier\nevery sentence against\na retrieved fact",
             FILL["agent"]),
            ("A_a", "answer\nevidence → caveats → prose", FILL["agent"])):
        out.append("  " + _box(nid, label, fill, width=1.7, height=0.55))
    out.append('    {rank=same; A_q A_o A_t}')
    out.append('    {rank=same; A_a A_v A_m}')
    out.append('    A_q -> A_o -> A_t;')
    out.append('    A_t -> A_m [label=" retrieved\\n facts"];')
    out.append('    A_m -> A_v;')
    out.append('    A_v -> A_a [label=" or WITHHELD"];')
    out.append('  }')
    out.append('  REL -> A_t [constraint=false, weight=0, '
               'color="#2a7f7f", label=" the graph"];')
    out.append('  RES -> A_o [constraint=false, weight=0, '
               'color="#2a7f7f", label=" which viruses\\n are validated"];')
    out.append('  DOC -> A_q [style=invis, weight=20];')
    return out


def _key(t: dict) -> list[str]:
    """What the reader has to know to read the picture.

    A legend rather than prose, because the figure is printed on its own
    and a caption does not travel with a PDF.
    """
    rows = [
        ("#e8eef8", "raw input", "a public database, as downloaded"),
        ("#e9f5e1", "phase-1 layer",
         "reads raw files only — can run in any order"),
        ("#fdebd9", "phase-2 layer",
         "reads the merged graph, so it runs after pass 1"),
        ("#ffe599", "the release",
         "nodes.jsonl, edges.jsonl and a MANIFEST that fingerprints both"),
        ("#f0eaf9", "measurement", "writes a stamped artifact, never a claim"),
        ("#e3f4f4", "agent", "answers from the release, refuses past it"),
    ]
    cells = []
    for fill, name, what in rows:
        cells.append(
            f'      <tr><td bgcolor="{fill}" width="14"> </td>'
            f'<td align="left"><b>{name}</b></td>'
            f'<td align="left">{what}</td></tr>')
    absent = (", ".join(t["absent"]) if t["absent"] else "none")
    note = (f'{t["n_layers_built"]} of {t["n_layers_declared"]} declared '
            f'layers are in this release (absent: {absent}). '
            f'Merge order is LAYER_PRECEDENCE: '
            f'{" → ".join(t["precedence"])}.')
    return ['  key [shape=plaintext, style="", label=<',
            '    <table border="1" color="#999999" cellborder="0" '
            'cellspacing="0" cellpadding="3" bgcolor="#fcfcfc">',
            '      <tr><td colspan="3" align="left">'
            '<font point-size="10"><b>How to read this</b></font></td></tr>',
            *cells,
            f'      <tr><td colspan="3" align="left">'
            f'<font point-size="7.5" color="#555555">{note}</font>'
            f'</td></tr>',
            '    </table>>];',
            '  A_a -> key [style=invis, weight=10];']


def build(wide: bool) -> str:
    t = topology()
    lines = ["// GENERATED FILE. Do not edit.",
             "//   make arch-figure         regenerate and re-render",
             "//   make arch-figure-check   fail if stale, write nothing",
             "//",
             "// The topology is derived: the layers, their scripts, their raw",
             "// inputs, which side of the merge each runs on, and every count.",
             f"// release {t['release']}  "
             f"nodes_sha {t['fingerprint']['nodes_sha']}  "
             f"edges_sha {t['fingerprint']['edges_sha']}",
             _header(t, wide).rstrip("\n")]
    lines += _stage1(t, wide)
    lines += _stage2(t, wide)
    lines += _stage3(t, wide)
    lines += _stage4(t, wide)
    lines += _key(t)
    lines.append("}")
    return "\n".join(lines) + "\n"


def _render() -> None:
    jobs = [
        (A4, "pdf", FIGURES / "kgav_architecture_a4.pdf", []),
        (A4, "png", FIGURES / "kgav_architecture_a4.png", ["-Gdpi=300"]),
        (WIDE, "png", FIGURES / "kgav_architecture.png", ["-Gdpi=150"]),
        (WIDE, "svg", FIGURES / "kgav_architecture.svg", []),
    ]
    for src, fmt, dst, extra in jobs:
        cmd = ["dot", f"-T{fmt}", *extra, str(src), "-o", str(dst)]
        subprocess.run(cmd, check=True)
        print(f"rendered {dst.relative_to(ROOT)}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if either .dot is stale, write nothing")
    ap.add_argument("--render", action="store_true",
                    help="also re-render the PNG, SVG and A4 PDF")
    args = ap.parse_args(argv)

    want = {WIDE: build(wide=True), A4: build(wide=False)}
    if args.check:
        bad = [p.name for p, text in want.items()
               if not p.exists() or p.read_text() != text]
        if bad:
            print(f"stale: {', '.join(sorted(bad))}\n"
                  f"run `make arch-figure` to regenerate.", file=sys.stderr)
            return 1
        print(f"architecture .dot files are current "
              f"({len(want)} files match the code and the release).")
        return 0

    FIGURES.mkdir(parents=True, exist_ok=True)
    for p, text in want.items():
        p.write_text(text)
        print(f"wrote {p.relative_to(ROOT)}  "
              f"({len(text.splitlines()):,} lines)")
    if args.render:
        _render()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
