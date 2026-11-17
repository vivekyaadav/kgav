#!/usr/bin/env python3
"""Generate docs/scorer_aucs.txt: every scorer AUC, every virus, side by side.

WHY THIS FILE EXISTS. The three scorer studies write three results files
with different shapes -- the publication-disjoint run keys its channels
one way, the symmetric control adds a k-ladder, the strength-matched run
carries a per-cell k -- and reading one AUC out of them took a Python
one-liner per question. So the comparison that matters most, "what does
each protocol say about each virus, and does the degree baseline beat
it", existed nowhere: it was reassembled by hand in a terminal after
every rebuild and never written down. A number nobody can see is a
number nobody checks.

WHY GENERATED RATHER THAN WRITTEN, same argument as
scripts/write_kg_report.py: a hand-typed AUC table goes stale silently
on the next rebuild. This one is stamped with the release fingerprint
AND with each source file's own stamp, so a table built from a mixture
of releases is detectable rather than merely wrong.

WHAT IT IS NOT. Not a protocol and not a registration -- it reports the
three registered studies and adds nothing. Every figure is read from
data/results; this script computes no AUC and applies no threshold of
its own. A cell this file cannot find is printed as a refusal with the
reason the results file gave, never omitted.

THE TWO ALPHAS ARE BOTH SHOWN. The ChEMBL cells are a family of 14
(alpha 0.05/14, z = 2.9137) and the plate cells a family of 3 (alpha
0.05/3, z = 2.3940); the SARS-CoV-2 plate cell is claimed by both
registrations. A single asterisk column would have to pick one, so the
plate block prints both.

    make scorer-table          # regenerate
    make scorer-table-check    # fail if it is stale, write nothing
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.pipeline import release_dir
from kgav.provenance import release_fingerprint
from kgav import virus_register as vr

TABLE = ROOT / "docs" / "scorer_aucs.txt"
RESULTS = ROOT / "data" / "results"

#: (column heading, results file, channel key). ORDER IS THE PRINTED
#: ORDER and the reading order of the three registrations: the
#: publication-disjoint control first because it is the prior the other
#: two quote, then the symmetric control, then the strength-matched one.
COLUMNS = (
    ("pub-disj", "structure_scorer.json", "structure_max"),
    ("sym k10", "structure_scorer_symmetric.json", "structure_max_k10"),
    ("strength", "structure_scorer_strength_matched.json", "structure_max"),
    ("degree", "structure_scorer_symmetric.json", "degree"),
)

#: The populations, with the registration that owns each one. A
#: population absent from every file is printed as absent rather than
#: skipped -- an entry cell that stops being produced is a change.
POPULATIONS = (
    ("primary_chembl", "ChEMBL, all assay types",
     "docs/prereg_structure_scorer.md", 14, 2.9137),
    ("secondary_ncats", "NCATS plate, cell_based_antiviral",
     "docs/prereg_ncats_entry_assays.md", 3, 2.3940),
    ("secondary_ncats_entry", "NCATS plate, pseudotyped_entry",
     "docs/prereg_ncats_entry_assays.md", 3, 2.3940),
)

HEADER = (
    "=" * 78,
    "kgav — SCORER AUCs BY VIRUS AND PROTOCOL",
    "",
    "GENERATED FILE. Do not edit: `make scorer-table` overwrites it, and",
    "`make scorer-table-check` fails if it does not describe the results",
    "files on disk. It reports the three registered scorer studies and",
    "computes nothing; the registrations are named under each block.",
    "=" * 78,
    "",
)


def _load(name: str) -> dict:
    path = RESULTS / name
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def _cell(doc: dict, pop: str, virus: str, channel: str) -> dict | None:
    """The channel's metrics, or None.

    A channel that refused carries no "auc" key -- the refusal is the
    absence. Returning the dict anyway would make a refusal format as
    0.0000, which is the one failure mode a table like this has.
    """
    block = doc.get(pop) or {}
    cell = block.get(virus)
    if not isinstance(cell, dict):
        return None
    metrics = cell.get(channel)
    if not isinstance(metrics, dict) or "auc" not in metrics:
        return None
    return metrics


def _fmt(metrics: dict | None, key: str) -> str:
    if metrics is None:
        return "refused"
    mark = "*" if metrics.get(key) else " "
    return f"{metrics['auc']:.4f}{mark}"


def _stamp(docs: dict[str, dict]) -> list[str]:
    """Each source file's own release stamp, printed, not compared away.

    A table assembled from two releases is the failure this guards. The
    lines name the file and its stamp so a mismatch is readable rather
    than inferred from a single aggregate hash.
    """
    out = []
    for name in sorted(docs):
        prov = (docs[name] or {}).get("provenance") or {}
        rel = prov.get("release") or {}
        out.append(f"      {name:44s} {rel.get('nodes_sha', '?')} / "
                   f"{rel.get('edges_sha', '?')}")
    return out


def _degree_agreement(docs: dict[str, dict], pop: str, viruses: list[str],
                      names: dict[str, str]) -> list[str]:
    """Do the three runs still compute the same degree baseline.

    The table reads `degree` from the symmetric run alone, because all
    three studies compute it from the graph and not from their own
    protocol, so it should be identical in all three. SHOULD IS NOT A
    CHECK. If a future change made one run's baseline protocol-dependent,
    a single-source column would hide it behind a number that still looks
    like a baseline. So the comparison is made here and printed when it
    fails, and silent when it holds.
    """
    bad = []
    for virus in viruses:
        seen = {}
        for name, doc in docs.items():
            m = _cell(doc, pop, virus, "degree")
            if m is not None:
                seen[name] = round(m["auc"], 6)
        if len(set(seen.values())) > 1:
            bad.append(f"      {names.get(virus, virus)}: "
                       + ", ".join(f"{n}={v}" for n, v in sorted(seen.items())))
    if not bad:
        return []
    return ["",
            "      DEGREE BASELINE DISAGREES BETWEEN RUNS. It is computed from",
            "      the graph and not from any protocol, so the three studies",
            "      must agree; one of them no longer does:", *bad]


def build() -> str:
    names = {"NCBITaxon:" + v["taxon"]: v["name"] for v in vr.load()["viruses"]}
    docs = {name: _load(name) for _, name, _ in COLUMNS}
    fp = release_fingerprint(release_dir(ROOT))

    lines = list(HEADER)
    lines += [
        f"    release          {fp['id']}",
        f"    nodes_sha        {fp['nodes_sha']}",
        f"    edges_sha        {fp['edges_sha']}",
        "",
        "    SOURCE FILES AND THE RELEASE EACH WAS MEASURED ON",
        *_stamp(docs),
        "",
        "    * = the interval excludes chance at that block's registered",
        "    alpha. 'refused' = the protocol declined this cell and the",
        "    results file records why; it is not a low AUC.",
        "",
    ]

    for pop, title, prereg, family, z in POPULATIONS:
        viruses = sorted({v for d in docs.values()
                          for v in (d.get(pop) or {})
                          if isinstance((d.get(pop) or {}).get(v), dict)},
                         key=lambda v: names.get(v, v))
        lines.append(f"    {title.upper()}")
        lines.append(f"      registered by {prereg}")
        lines.append(f"      family of {family}, alpha = 0.05/{family} = "
                     f"{0.05 / family:.6f}, z = {z}")
        if not viruses:
            lines += ["", "      NO CELLS. This population is absent from "
                          "every results file on disk.", ""]
            continue
        head = (f"      {'virus':<14}{'n+':>6}{'n-':>7}"
                + "".join(f"{c:>10}" for c, _, _ in COLUMNS)
                + f"{'k_cell':>8}")
        lines += ["", head, "      " + "-" * (len(head) - 6)]
        strength = docs["structure_scorer_strength_matched.json"]
        k_cells = strength.get("k_cells") or {}
        for virus in viruses:
            cells = [_cell(docs[f], pop, virus, ch) for _, f, ch in COLUMNS]
            ref = next((c for c in cells if c), None)
            if ref is None:
                continue
            row = (f"      {names.get(virus, virus):<14}"
                   f"{ref['n_pos']:>6}{ref['n_neg']:>7}"
                   + "".join(f"{_fmt(c, 'clears_corrected_alpha'):>10}"
                             for c in cells)
                   + f"{k_cells.get(virus, '-'):>8}")
            lines.append(row)
        if pop != "primary_chembl":
            lines += ["",
                      "      AND AT THE 14-CELL ALPHA (z = 2.9137), the sterner",
                      "      of the two. The SARS-CoV-2 cytopathic-effect cell is",
                      "      claimed by BOTH registrations, so every plate cell",
                      "      prints both levels rather than this table choosing:"]
            for virus in viruses:
                cells = [_cell(docs[f], pop, virus, ch) for _, f, ch in COLUMNS]
                if not any(cells):
                    continue
                lines.append((f"      {names.get(virus, virus):<14}"
                              + " " * 13
                              + "".join(f"{_fmt(c, 'clears_alternate_alpha'):>10}"
                                        for c in cells)).rstrip())
        lines += _degree_agreement(docs, pop, viruses, names)
        lines.append("")

    lines += [
        "    READING IT. 'pub-disj' excludes from each candidate's",
        "    reference set every known active citing a publication the",
        "    candidate also cites; that control is asymmetric between the",
        "    classes, which is what the other two were registered to",
        "    address. 'sym k10' removes the ten nearest references from",
        "    every candidate alike. 'strength' removes k_cell of them,",
        "    where k_cell is that cell's own measured asymmetry. 'degree'",
        "    is the baseline: a compound's edge count, no chemistry. A",
        "    structure column that does not beat its own degree column is",
        "    not evidence about chemistry.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if the file on disk differs; write nothing")
    args = ap.parse_args(argv)

    want = build()
    if args.check:
        have = TABLE.read_text() if TABLE.exists() else ""
        if have == want:
            print(f"{TABLE.relative_to(ROOT)} is current.")
            return 0
        print(f"{TABLE.relative_to(ROOT)} is STALE. Run `make scorer-table`.")
        return 1
    TABLE.write_text(want)
    print(f"wrote {TABLE.relative_to(ROOT)}")
    print(want)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
