#!/usr/bin/env python3
"""Where the AUC actually comes from, as a share per approach.

    make auc-contribution

THE QUESTION THIS ANSWERS. "Percentage-wise, what did each approach
contribute to the AUC?" -- asked of a project that has four scorer
families, eight metapath channels, three leakage controls and a degree
baseline.

WHY A SCRIPT AND NOT A PARAGRAPH. Every figure below is read from
data/results. A decomposition typed into a document is a decomposition
that was true on the day it was typed; this one is re-derived from the
files on disk, and prints the release fingerprint each file was measured
on so a mixture is visible rather than merely wrong.

THE ARITHMETIC, STATED BECAUSE IT IS A CHOICE AND NOT A FACT.

  excess(x) = AUC(x) - 0.5

An AUC of 0.5 is chance, so the distance from 0.5 is the whole of what a
scorer contributes; a share of the AUC itself would credit chance with
50 percentage points and flatter every method equally.

AUC IS NOT ADDITIVE AND THESE APPROACHES ARE NOT STACKED LAYERS. There
is no sequence of steps whose excesses sum to a total, so a single pie
chart would be a fabrication. What CAN be measured, and is measured
here, is three different things, reported separately:

  1. SUBSTITUTION. Each approach's own excess on the same population, so
     they can be compared. The share column is that approach's excess
     over the sum of the positive excesses -- a statement about relative
     strength, not about addition.

  2. DECOMPOSITION OF THE STRONGEST CLAIM. For every cell that carries a
     structure AUC, how much of its excess the degree baseline alone
     reproduces. Degree is edge count with no chemistry in it, so that
     share is the part of the result that is NOT about chemistry.

  3. WHAT THE CONTROLS COST. The leakage controls do not raise AUC, they
     lower it, and reporting only the uncontrolled figure is how a
     project talks itself into a result. Those rows are negative on
     purpose.

NOTHING HERE IS A NEW MEASUREMENT. No AUC is computed, no interval is
widened, no threshold is applied. Every number is read from a results
file and every verdict ('clears') is the one that file already recorded.
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

RESULTS = ROOT / "data" / "results"
CHANCE = 0.5

#: The structure studies, sternest control last. Each is (label, file, key).
STRUCTURE = (
    ("publication-disjoint", "structure_scorer.json", "structure_max"),
    ("symmetric k=10", "structure_scorer_symmetric.json", "structure_max_k10"),
    ("strength-matched", "structure_scorer_strength_matched.json",
     "structure_max"),
)
#: The baseline, which lives in the symmetric run's file.
DEGREE = ("structure_scorer_symmetric.json", "degree")

RULE = "=" * 78
THIN = "-" * 78


def load(name: str) -> dict:
    path = RESULTS / name
    if not path.exists():
        raise SystemExit(f"{path} is missing; `make check-results` says which "
                         f"results files are current")
    return json.loads(path.read_text())


def stamp(doc: dict) -> str:
    """The release a results file says it was measured on."""
    rel = (doc.get("provenance") or {}).get("release") or {}
    return f"{rel.get('nodes_sha', '?')} / {rel.get('edges_sha', '?')}"


def cells(doc: dict) -> dict[str, dict]:
    """Every virus cell in a structure-scorer file, ChEMBL and plate alike.

    The plate blocks are keyed by assay type under `secondary_ncats*`, so
    a cell is identified by (section, taxon) rather than by taxon alone --
    SARS-CoV-2 appears in three of them and they are different
    populations with different answers.
    """
    out: dict[str, dict] = {}
    for section in ("primary_chembl", "secondary_ncats",
                    "secondary_ncats_entry"):
        block = doc.get(section) or {}
        for taxon, cell in block.items():
            out[f"{section}|{taxon}"] = cell
    return out


def auc(cell: dict, key: str) -> float | None:
    """The AUC of one channel of one cell, or None if it was refused.

    `measurable` false means the protocol declined the cell and recorded
    why. That is not a low AUC and must never be read as one -- a refused
    cell whose arithmetic AUC is 0.5 would otherwise enter the average as
    'exactly chance', which is a measurement nobody made.
    """
    ch = (cell or {}).get(key)
    if not isinstance(ch, dict):
        return None
    if ch.get("measurable") is False or ch.get("refusal"):
        return None
    val = ch.get("auc")
    return float(val) if isinstance(val, (int, float)) else None


def clears(cell: dict, key: str) -> bool:
    ch = (cell or {}).get(key) or {}
    return bool(ch.get("clears_corrected_alpha")) and ch.get(
        "direction") != "below"


#: taxon CURIE -> the display name, from the register. Spelled the way
#: scripts/write_scorer_table.py spells it: a virus name comes from
#: config/viruses.yaml and never from prose or a graph node property.
VIRUS_NAMES = {"NCBITaxon:" + v["taxon"]: v["name"] for v in vr.load()["viruses"]}


def name_of(label: str) -> str:
    section, taxon = label.split("|", 1)
    nice = VIRUS_NAMES.get(taxon, taxon)
    tag = {"primary_chembl": "ChEMBL",
           "secondary_ncats": "plate/CPE",
           "secondary_ncats_entry": "plate/entry"}.get(section, section)
    return f"{nice} ({tag})"


def pct(x: float) -> str:
    return f"{100 * x:6.1f}%"


# ======================================================================
def section_1_substitution(out) -> None:
    """Each approach's own excess, on the populations it was measured on.

    SCORERS AND BASELINES ARE SEPARATED, and only scorers are in the
    denominator. The first version of this table put the degree baseline
    in it, which handed 36% of "the contribution" to a control -- a row
    whose entire purpose is to be the thing the approaches are measured
    AGAINST. A baseline that scores well is evidence against the approach
    above it, not a contribution of its own, so baselines print below the
    rule with their share column blank.
    """
    out(RULE)
    out("1. WHAT EACH APPROACH CONTRIBUTES ON ITS OWN")
    out(RULE)
    out("")
    out("   Each row is one approach's best measured cell. 'share' is that")
    out("   approach's excess over the summed excess of the APPROACHES --")
    out("   a comparison of strength, not a decomposition of a total,")
    out("   because these are separate scorers on overlapping populations")
    out("   and they do not stack.")
    out("")
    out("   The baselines below the rule are controls, not approaches.")
    out("   They are excluded from the denominator: a baseline that scores")
    out("   well is evidence against the approach it is compared with.")
    out("")

    # (label, where, auc, note, is_baseline)
    rows: list[tuple[str, str, float, str, bool]] = []

    # ---- the graph: metapath channels, on the assay-matched population
    hn = load("hard_negatives.json")
    for virus, chans in (hn.get("cross_sectional") or {}).items():
        where = f"{virus}, plate, publication-disjoint"
        chan_cells = {k: c for k, c in chans.items()
                      if k not in ("degree", "COMBINED")
                      and isinstance(c, dict)}
        scored = {k: c["auc"] for k, c in chan_cells.items()
                  if c.get("evaluable") and isinstance(c.get("auc"),
                                                       (int, float))}
        if scored:
            best_key = max(scored, key=lambda k: scored[k])
            rows.append((f"KG metapath, best single channel ({best_key})",
                         where,
                         scored[best_key],
                         f"{len(scored)} of {len(chan_cells)} channels scored",
                         False))
        comb = chans.get("COMBINED") or {}
        if isinstance(comb.get("auc"), (int, float)):
            rows.append(("KG metapaths, all channels COMBINED", where,
                         comb["auc"], "", False))
        deg = chans.get("degree") or {}
        if isinstance(deg.get("auc"), (int, float)):
            rows.append(("degree baseline — edge count, no chemistry", where,
                         deg["auc"], "the control", True))

    # ---- chemistry: the structure scorer, under its sternest control
    sm = load("structure_scorer_strength_matched.json")
    sm_cells = cells(sm)
    scored = {lbl: auc(c, "structure_max") for lbl, c in sm_cells.items()
              if auc(c, "structure_max") is not None}
    best = max(scored, key=lambda k: scored[k])
    n_clear = sum(1 for c in sm_cells.values() if clears(c, "structure_max"))
    rows.append(("structure ECFP4 max-Tanimoto (strength-matched)",
                 f"best of {len(scored)} cells: {name_of(best)}",
                 scored[best],
                 f"{n_clear} of {len(sm_cells)} cells clear corrected alpha",
                 False))

    # ---- chemistry, forward: the temporal scorer.
    # THE AUC IS UNDER `primary`, AND `clears` IS A LIST OF TAXA, NOT A
    # COUNT. Read as cell["auc"] and len()-less, this row vanished from
    # the table in silence -- the only out-of-time result in the project,
    # absent from a summary of what the project's approaches contribute,
    # with nothing printed to say so.
    ts = load("temporal_structure_scorer.json")
    tcells = ts.get("cells") or {}
    tvals = {k: c["primary"]["auc"] for k, c in tcells.items()
             if isinstance((c.get("primary") or {}).get("auc"), (int, float))}
    if tvals:
        tbest = max(tvals, key=lambda k: tvals[k])
        n_clear = len(ts.get("clears") or [])
        n_comp = len(ts.get("comparator_usable") or [])
        rows.append((f"structure, PROSPECTIVE (train<={ts['cutoff']}, "
                     f"test>{ts['cutoff']})",
                     f"best of {len(tvals)} cells: "
                     f"{tcells[tbest].get('name', tbest)}",
                     tvals[tbest],
                     f"{n_clear} of {ts.get('n_evaluable', len(tvals))} "
                     f"clear; {n_comp} usable degree comparators",
                     False))

    # ---- the clinical axis
    cl = load("clinical_labels.json")
    disease = cl.get("disease", "?")
    for key, label, is_base in (
            ("structure_max_k10", "structure, clinical labels", False),
            ("degree", "degree baseline, clinical labels", True)):
        ch = cl.get(key) or {}
        if isinstance(ch.get("auc"), (int, float)):
            rows.append((f"{label} ({disease})",
                         f"{cl.get('virus', '?')}, TREATS/IN_TRIAL_FOR edges",
                         ch["auc"], "the control" if is_base else "",
                         is_base))

    approaches = [r for r in rows if not r[4]]
    baselines = [r for r in rows if r[4]]
    total = sum(max(0.0, r[2] - CHANCE) for r in approaches)

    out(f"   {'approach':<48}{'AUC':>8}{'excess':>9}{'share':>9}")
    out("   " + THIN[:74])
    for label, where, a, note, _ in sorted(approaches, key=lambda r: -r[2]):
        ex = a - CHANCE
        share = pct(max(0.0, ex) / total) if total else "     --"
        out(f"   {label:<48}{a:8.4f}{ex:+9.4f}{share:>9}")
        out(f"     {where}" + (f"  [{note}]" if note else ""))
    out("   " + THIN[:74])
    out(f"   {'sum of approach excesses':<48}{'':>8}{total:+9.4f}   100.0%")
    out("")
    out("   BASELINES — controls, not in the denominator above")
    out("   " + THIN[:74])
    for label, where, a, note, _ in sorted(baselines, key=lambda r: -r[2]):
        out(f"   {label:<48}{a:8.4f}{a - CHANCE:+9.4f}{'':>9}")
        out(f"     {where}" + (f"  [{note}]" if note else ""))
    out("")


def section_2_how_much_is_chemistry(out) -> None:
    """Per cell: the share of the structure excess that degree reproduces."""
    out(RULE)
    out("2. HOW MUCH OF THE STRUCTURE RESULT IS NOT CHEMISTRY")
    out(RULE)
    out("")
    out("   The degree baseline ranks a compound by its edge count. It")
    out("   contains no chemistry at all, so whatever it reproduces of a")
    out("   structure cell's excess is the part of that cell NOT earned by")
    out("   chemical similarity. Per cell, under the sternest control:")
    out("")
    out("   'deg/struct' is capped for display at >100%. Uncapped, a cell")
    out("   whose structure excess is a few thousandths prints a ratio in")
    out("   the thousands of percent -- a number that says only that the")
    out("   denominator is near zero, while reading as a measurement.")
    out("")

    sm = cells(load("structure_scorer_strength_matched.json"))
    deg = cells(load(DEGREE[0]))

    out(f"   {'cell':<30}{'struct':>8}{'degree':>8}"
        f"{'struct ex':>11}{'deg/struct':>12}")
    out("   " + THIN[:71])
    shares, counted, skipped, degree_wins = [], 0, [], []
    for label in sorted(sm, key=lambda k: name_of(k)):
        s_auc = auc(sm[label], "structure_max")
        d_auc = auc(deg.get(label, {}), DEGREE[1])
        if s_auc is None or d_auc is None:
            skipped.append(name_of(label))
            continue
        s_ex, d_ex = s_auc - CHANCE, d_auc - CHANCE
        counted += 1
        if s_ex <= 0:
            shown = "  below .5"
        else:
            frac = max(0.0, d_ex) / s_ex
            shares.append(min(frac, 1.0))
            if frac >= 1.0:
                shown = "     >100%"
                degree_wins.append(name_of(label))
            else:
                shown = pct(frac)
        out(f"   {name_of(label):<30}{s_auc:8.4f}{d_auc:8.4f}{s_ex:+11.4f}"
            f"{shown:>12}")
    out("   " + THIN[:71])
    if shares:
        shares.sort()
        n = len(shares)
        mid = shares[n // 2] if n % 2 else (shares[n // 2 - 1]
                                            + shares[n // 2]) / 2
        out(f"   {counted} cells compared; {n} have a positive structure "
            f"excess.")
        out("   Shares are capped at 100% before averaging, so a cell where")
        out("   degree wins counts as 'all of it', not as 34 times all of it.")
        out(f"   median share reproduced by degree alone: {pct(mid)}")
        out(f"   mean:                                    "
            f"{pct(sum(shares) / n)}")
        out(f"   cells where degree matches or beats structure: "
            f"{len(degree_wins)} of {n}")
        if degree_wins:
            out("     " + ", ".join(degree_wins))
    if skipped:
        out(f"   not compared ({len(skipped)}): " + ", ".join(skipped))
        out("     (a cell the protocol refused, with the reason in the "
            "results file; not a low AUC)")
    out("")


def section_3_what_the_controls_cost(out) -> None:
    out(RULE)
    out("3. WHAT THE LEAKAGE CONTROLS COST")
    out(RULE)
    out("")
    out("   These rows are NEGATIVE contributions and they are the point.")
    out("   Each control removes a way the scorer could have been reading")
    out("   its own answer sheet. The uncontrolled column is what the same")
    out("   method would have reported without them.")
    out("")

    docs = {f: load(f) for _, f, _ in STRUCTURE}
    per = {lbl: cells(docs[f]) for lbl, f, _ in STRUCTURE}
    keys = {lbl: k for lbl, _, k in STRUCTURE}

    out(f"   {'cell':<30}" + "".join(f"{lbl.split()[0][:9]:>11}"
                                     for lbl, _, _ in STRUCTURE))
    out("   " + THIN[:63])
    common = set.intersection(*(set(c) for c in per.values()))
    deltas: dict[str, list[float]] = {lbl: [] for lbl, _, _ in STRUCTURE}
    base_lbl = STRUCTURE[1][0]          # symmetric k=10, the loosest numeric
    for label in sorted(common, key=name_of):
        vals = {lbl: auc(per[lbl][label], keys[lbl]) for lbl in per}
        if vals[base_lbl] is None:
            continue
        line = f"   {name_of(label):<30}"
        for lbl, _, _ in STRUCTURE:
            v = vals[lbl]
            line += f"{(f'{v:.4f}' if v is not None else 'refused'):>11}"
            if v is not None and vals[base_lbl] is not None:
                deltas[lbl].append(v - vals[base_lbl])
        out(line)
    out("   " + THIN[:63])
    out("")
    out(f"   mean AUC change against '{base_lbl}', the loosest numeric")
    out("   control, over the cells all three could measure:")
    for lbl, _, _ in STRUCTURE:
        d = deltas[lbl]
        if d:
            out(f"     {lbl:<26}{sum(d) / len(d):+8.4f}   ({len(d)} cells)")
    out("")


def section_4_the_graph(out) -> None:
    out(RULE)
    out("4. THE KNOWLEDGE GRAPH'S OWN CONTRIBUTION")
    out(RULE)
    out("")

    mp = load("metapath_comparison.json")
    per_virus = mp.get("per_virus") or {}
    total = readable = beats = both = 0
    for chans in per_virus.values():
        for key, ch in chans.items():
            if key in ("degree", "COMBINED") or not isinstance(ch, dict):
                continue
            total += 1
            if ch.get("vs_degree_readable"):
                readable += 1
            if ch.get("beats_degree"):
                beats += 1
            if ch.get("beats_degree") and ch.get("vs_degree_readable"):
                both += 1
    out(f"   metapath channel-cells measured:                    {total:4d}")
    out(f"     beat the degree baseline on lift:                 {beats:4d}")
    out(f"     readable against degree at all:                   {readable:4d}")
    out(f"     BOTH beat degree and are readable against it:     {both:4d}")
    out("")

    hn = load("hard_negatives.json")
    for virus, chans in (hn.get("cross_sectional") or {}).items():
        comb = (chans.get("COMBINED") or {}).get("auc")
        deg = (chans.get("degree") or {}).get("auc")
        if isinstance(comb, (int, float)) and isinstance(deg, (int, float)):
            out(f"   {virus}, assay-matched plate, publication-disjoint:")
            out(f"     all eight channels COMBINED      AUC {comb:.4f}"
                f"  ({comb - CHANCE:+.4f})")
            out(f"     degree baseline alone            AUC {deg:.4f}"
                f"  ({deg - CHANCE:+.4f})")
            out(f"     the graph's reasoning adds       {comb - deg:+.4f}")
            out("")

    # THE PATHWAY CHANNEL. Two channels in one file: `taxa` is the
    # dependency channel and `interactome` the PPI one. Counted together
    # and reported separately, because the 62->63 move of 2026-11-18 was
    # misattributed to the interactome for exactly the want of that split.
    land = load("pathway_landscape.json")
    for channel, key in (("dependency", "taxa"), ("interactome",
                                                  "interactome")):
        rows = land.get(key) or {}
        if isinstance(rows, dict):
            rows = list(rows.values())
        rows = [r for r in rows if isinstance(r, dict)]
        if not rows:
            continue
        cleared = sum(1 for r in rows if r.get("clears_corrected_alpha")
                      or r.get("clears"))
        out(f"   pathway overlap, {channel} channel:")
        out(f"     rows measured:                                    "
            f"{len(rows):4d}")
        out(f"     clearing the corrected alpha:                     "
            f"{cleared:4d}")
        out("")


def section_5_prospective(out) -> None:
    """The only out-of-time measurement in the project, cell by cell.

    Everything in sections 1 to 4 is cross-sectional: the labels were
    split by publication or by assay, never by time, so a scorer could
    be reading chemistry that was already known when the test compound
    was published. This section trains on what existed at the cutoff and
    tests on what came after, which is the only arrangement that answers
    "would this have worked".
    """
    ts = load("temporal_structure_scorer.json")
    cellmap = ts.get("cells") or {}
    clears_set = set(ts.get("clears") or [])
    usable = set(ts.get("comparator_usable") or [])

    out(RULE)
    out("5. THE ONLY PROSPECTIVE MEASUREMENT")
    out(RULE)
    out("")
    out(f"   split {ts.get('split', '?')}, cutoff {ts.get('cutoff')}, "
        f"undated {ts.get('undated_policy')}, "
        f"overlap {ts.get('overlap_policy')}")
    out(f"   family of {ts.get('family_size')}, alpha "
        f"{ts.get('alpha'):.6f}, z {ts.get('z_corrected')}")
    out("")
    out("   Trained on chemistry published at or before the cutoff, tested")
    out("   on activity published after it. A cell 'clears' when its")
    out("   corrected interval excludes chance ABOVE it.")
    out("")
    out(f"   {'cell':<16}{'n+':>6}{'n-':>7}{'AUC':>9}{'excess':>9}"
        f"{'degree':>9}  verdict")
    out("   " + THIN[:71])

    rows = []
    for taxon, cell in cellmap.items():
        a = (cell.get("primary") or {}).get("auc")
        if not isinstance(a, (int, float)):
            continue
        rows.append((taxon, cell, a))
    for taxon, cell, a in sorted(rows, key=lambda r: -r[2]):
        d = (cell.get("degree") or {}).get("auc")
        d_txt = f"{d:.4f}" if isinstance(d, (int, float)) else "--"
        if taxon in clears_set:
            verdict = "clears ABOVE chance"
        elif a < CHANCE and (cell.get("primary") or {}).get(
                "clears_corrected_alpha"):
            verdict = "clears BELOW chance"
        else:
            verdict = "does not clear"
        out(f"   {cell.get('name', taxon):<16}"
            f"{cell['n_test_positives']:>6}{cell['n_test_negatives']:>7}"
            f"{a:9.4f}{a - CHANCE:+9.4f}{d_txt:>9}  {verdict}")
    out("   " + THIN[:71])
    out(f"   {len(rows)} cells evaluable of {ts.get('n_cells_total')}; "
        f"{len(clears_set)} clear above chance.")
    out(f"   usable degree comparators: {len(usable)} of {len(rows)}.")
    out("   A cell with no usable comparator has no baseline to beat: the")
    out("   graph reaches too few of the compounds being ranked to score")
    out("   them, which is a fact about the graph and not about the scorer.")
    out("")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path,
                    help="also write the report to this file")
    args = ap.parse_args()

    lines: list[str] = []

    def out(s: str = "") -> None:
        lines.append(s)

    fp = release_fingerprint(release_dir(ROOT))
    out(RULE)
    out("kgav — AUC CONTRIBUTION BY APPROACH")
    out(RULE)
    out("")
    out(f"   release    {fp['id']}   nodes {fp['nodes_sha']} / "
        f"edges {fp['edges_sha']}")
    out("")
    out("   SOURCE FILES AND THE RELEASE EACH RECORDS")
    for name in sorted({f for _, f, _ in STRUCTURE}
                       | {"hard_negatives.json", "metapath_comparison.json",
                          "temporal_structure_scorer.json",
                          "clinical_labels.json",
                          "pathway_landscape.json"}):
        try:
            out(f"     {name:<44}{stamp(load(name))}")
        except SystemExit:
            out(f"     {name:<44}MISSING")
    out("")
    out("   excess(x) = AUC(x) - 0.5. Chance contributes nothing, so the")
    out("   distance from 0.5 is the whole of what a scorer supplies.")
    out("")

    section_1_substitution(out)
    section_2_how_much_is_chemistry(out)
    section_3_what_the_controls_cost(out)
    section_4_the_graph(out)
    section_5_prospective(out)

    out(RULE)
    out("READING IT")
    out(RULE)
    out("")
    out("   A share in section 1 is a share of the summed positive excess,")
    out("   which is a comparison of strength and NOT a decomposition of a")
    out("   total: these scorers were measured on overlapping populations")
    out("   and do not stack. Section 2 is the one genuine decomposition")
    out("   here, because degree and structure rank the same candidates in")
    out("   the same cell.")
    out("")

    text = "\n".join(lines)
    print(text)
    if args.out:
        args.out.write_text(text + "\n")
        print(f"\nwrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
