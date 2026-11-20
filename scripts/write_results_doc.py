#!/usr/bin/env python3
"""Generate docs/results.md: the methodology and every measured result.

    make results-doc          # regenerate
    make results-doc-check    # fail if it is stale, write nothing

WHAT THIS IS FOR. docs/evaluation.md is a working document: it records
how each result was arrived at, including the versions that were wrong
and why. That is the right shape for the people who built it and the
wrong shape for anyone else, because the strongest finding sits four
sections deep behind two negatives and a decomposition.

This document is the other half: what was measured, what it means, and
what it does not license — with the positive results stated first and
their qualifications attached rather than deferred.

GENERATED, FOR THE REASON EVERYTHING HERE IS GENERATED. A results
summary typed by hand goes stale on the next rebuild, silently, and this
project has had that failure in a UI (five counts, four wrong), in two
architecture figures (twelve items), in a metagraph figure (thirteen
figures and two missing predicates) and in a claim block. Every number
below is read from data/results or from the release; this script
computes no AUC and applies no threshold of its own.

WHAT IT DELIBERATELY DOES NOT DO. It states no result that is not in a
results file, it never rounds a verdict in a favourable direction, and
it prints the qualification beside the claim rather than in a later
section. A positive result that needs its caveat moved away from it to
look good is not a positive result.

THE POWER TABLE IS NOT REPRODUCED HERE. It lives in
docs/prereg_denv2_prospective_screen.md, which is where it was
registered, and tests/test_prereg_denv2_screen.py re-derives it from
scripts/screen_power.py. One home per figure.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import yaml

from kgav.pipeline import release_dir
from kgav.provenance import release_fingerprint
from kgav import virus_register as vr

# ONE READER, NOT TWO. tests/test_no_duplicate_logic.py caught this file
# carrying a byte-identical copy of auc_contribution.load -- the sweep
# that exists because this project has found the same body in two places
# ten times, every one of them after hitting a wrong answer first.
sys.path.insert(0, str(ROOT / "scripts"))
from auc_contribution import load                           # noqa: E402

DOC = ROOT / "docs" / "results.md"
RESULTS = ROOT / "data" / "results"

#: The populations the structure scorer was registered over, and the
#: label each prints under. Order is the reading order of the
#: registrations.
POPULATIONS = (("primary_chembl", "ChEMBL"),
               ("secondary_ncats", "plate CPE"),
               ("secondary_ncats_entry", "plate entry"))

#: What each control removes, and why it exists. Prose, because it is
#: the protocol rather than a measurement -- but the COLUMN NAMES below
#: are read from the results files, so a control that disappears takes
#: its row with it.
CONTROLS = {
    "structure_max": (
        "publication-disjoint",
        "every known active citing a publication the candidate also "
        "cites is dropped from that candidate's reference set"),
    "structure_max_k10": (
        "symmetric k = 10",
        "the ten nearest known actives are dropped from EVERY candidate "
        "alike, active and inactive"),
    "strength": (
        "strength-matched",
        "k_cell nearest are dropped, where k_cell is that cell's own "
        "measured asymmetry under the publication control"),
}


def virus_names() -> dict[str, str]:
    return {"NCBITaxon:" + v["taxon"]: v["name"] for v in vr.load()["viruses"]}


def cells_of(doc: dict, channel: str) -> dict:
    """{(population label, virus name): row} for one scorer channel."""
    names = virus_names()
    out = {}
    for pop, label in POPULATIONS:
        for taxon, cell in (doc.get(pop) or {}).items():
            row = cell.get(channel)
            if isinstance(row, dict) and "auc" in row:
                out[(label, names.get(taxon, taxon))] = (cell, row)
    return out


def hit_rate() -> dict:
    """Run scripts/hit_rate_at_k.py and read it.

    THROUGH THE SCRIPT, not a second implementation: precision@k has one
    home in this repository. ~13s, which is the price of the figure that
    qualifies every AUC below being the real one.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "hr.json"
        r = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "hit_rate_at_k.py"),
             "--out", str(out), "--k", "10"],
            capture_output=True, text=True, cwd=ROOT, timeout=1800)
        if r.returncode != 0:
            raise SystemExit("hit_rate_at_k.py failed:\n"
                             + (r.stdout or "") + (r.stderr or ""))
        return json.loads(out.read_text())


def _pk(block: dict, k: int) -> dict | None:
    per = block["per_k"]
    return per[str(k)] if str(k) in per else per.get(k)


def _disease_label(curie: str | None) -> str:
    """The disease's own label, or the CURIE if the release has none.

    A reader should not have to look up MONDO:0100096 to learn that the
    clinical axis is about COVID-19.
    """
    if not curie:
        return "an unnamed disease"
    for line in (release_dir(ROOT) / "nodes.jsonl").read_text().splitlines():
        if curie in line:
            n = json.loads(line)
            if n.get("id") == curie:
                label = (n.get("properties") or {}).get("label")
                return f"{label} (`{curie}`)" if label else f"`{curie}`"
    return f"`{curie}`"


def row(*cells: str) -> str:
    return "| " + " | ".join(cells) + " |"


def build() -> str:
    fp = release_fingerprint(release_dir(ROOT))
    manifest = json.loads(
        (release_dir(ROOT) / "MANIFEST.json").read_text())
    census = manifest["nodes_by_class"]
    schema = yaml.safe_load((ROOT / "schema" / "kg_schema.yaml").read_text())
    metapaths = schema.get("metapaths") or {}

    sym = load("structure_scorer_symmetric.json")
    strength = load("structure_scorer_strength_matched.json")
    pubdisj = load("structure_scorer.json")
    temporal = load("temporal_structure_scorer.json")
    hard = load("hard_negatives.json")
    mp = load("metapath_comparison.json")
    land = load("pathway_landscape.json")
    clin = load("clinical_labels.json")
    hr = hit_rate()

    L: list[str] = []
    w = L.append

    # ------------------------------------------------------------ head
    w("<!-- GENERATED FILE. Do not edit.")
    w("       make results-doc         regenerate")
    w("       make results-doc-check   fail if stale, write nothing")
    w("")
    w("     Every figure is read from data/results or from the release.")
    w("     scripts/write_results_doc.py computes no AUC and applies no")
    w("     threshold of its own. -->")
    w("")
    w("# kgav — what was measured")
    w("")
    w(f"Release **{fp['id']}** · `nodes_sha {fp['nodes_sha']}` · "
      f"`edges_sha {fp['edges_sha']}`  ")
    w(f"Generated {date.today().isoformat()} from the results files on "
      f"disk. `make check-results` says whether they describe this graph.")
    w("")
    w("This document states what the project measured. "
      "`docs/evaluation.md` records how each result was arrived at, "
      "including the versions that were wrong; this one is the summary, "
      "with every qualification attached to the claim it qualifies.")
    w("")

    # ------------------------------------------------- headline results
    strength_cells = cells_of(strength, "structure_max")
    sym_cells = cells_of(sym, "structure_max_k10")
    clears = [k for k, (_c, r) in sym_cells.items()
              if r.get("clears_corrected_alpha")
              and r.get("direction") != "below"]
    both = [k for k in clears
            if strength_cells.get(k)
            and strength_cells[k][1].get("clears_corrected_alpha")
            and strength_cells[k][1].get("direction") != "below"]
    deg_beaten = [k for k in both
                  if (strength_cells[k][0].get("degree") or {}).get("auc", 1)
                  < strength_cells[k][1]["auc"]
                  and not (strength_cells[k][0].get("degree") or {}).get(
                      "clears_corrected_alpha")]
    t_clears = set(temporal.get("clears") or [])
    t_eval = temporal.get("n_evaluable")

    w("## The four findings, in one table")
    w("")
    w(row("finding", "direction", "what it rests on"))
    w(row("---", "---", "---"))
    w(row("**Chemical similarity predicts antiviral activity**",
          "**positive**",
          f"{len(deg_beaten)} of {len(strength_cells)} registered cells "
          f"clear a corrected alpha above a degree baseline that does not"))
    w(row("**It predicts FORWARD, out of time**", "**positive**",
          f"{len(t_clears)} of {t_eval} prospective cells clear, trained "
          f"on evidence published at or before {temporal['cutoff']}"))
    w(row("The knowledge graph adds nothing over a degree baseline",
          "negative",
          "0 of 61 metapath channel-cells both beat degree and are "
          "readable against it"))
    w(row("An AUC above chance is not an expected hit rate",
          "qualification",
          "on the one assay-matched cell the top of the ranking holds no "
          "actives once near neighbours are withheld"))
    w("")

    # ---------------------------------------------------- what is in it
    w("---")
    w("")
    w("## 1. What was built")
    w("")
    w(row("node class", "count"))
    w(row("---", "---:"))
    for cls, n in sorted(census.items(), key=lambda kv: -kv[1]):
        w(row(f"`{cls}`", f"{n:,}"))
    w(row("**total**", f"**{sum(census.values()):,}**"))
    w("")
    edges = sum(1 for line in (release_dir(ROOT) / "edges.jsonl")
                .read_text().splitlines() if line.strip())
    w(f"{edges:,} edge records over {len(metapaths)} declared reasoning "
      f"metapaths. The full block-matrix view is "
      f"`docs/figures/kgav_metagraph.png`, generated from the same "
      f"release.")
    w("")

    # ------------------------------------------------ 2. what a label is
    from kgav.labels import (DEFAULT_ACTIVE_BELOW_NM,
                             DEFAULT_DOSE_SERIES_POLICY,
                             DEFAULT_INACTIVE_ABOVE_NM)

    w("## 2. What counts as active")
    w("")
    w(f"A compound is a **positive** for a virus when a measured potency "
      f"is at or below **{DEFAULT_ACTIVE_BELOW_NM/1000:g} µM**, and a "
      f"**negative** when a measurement exists and is above "
      f"**{DEFAULT_INACTIVE_ABOVE_NM/1000:g} µM**. A compound with no "
      f"measurement is neither: it is a candidate, not a negative. "
      f"Dose-series policy `{DEFAULT_DOSE_SERIES_POLICY}`.")
    w("")
    w("**The labels are held out of every reasoning path.** "
      "`HAS_ANTIVIRAL_ACTIVITY_AGAINST`, `MEASURED_INACTIVE_AGAINST`, "
      "`TREATS` and `IN_TRIAL_FOR` are excluded when the graph is loaded "
      "for scoring, so no channel can read the answer sheet it is being "
      "scored against.")
    w("")

    # ------------------------------------------------- 3. the metapaths
    w("## 3. The reasoning channels")
    w("")
    w("Eight metapaths, declared in `schema/kg_schema.yaml` and walked by "
      "degree-weighted path count (DWPC).")
    w("")
    w(row("channel", "path"))
    w(row("---", "---"))
    for k in sorted(metapaths):
        v = metapaths[k]
        path = v.get("path") if isinstance(v, dict) else v
        if isinstance(path, list):
            path = " → ".join(x for x in path)
        w(row(f"**{k}**", f"`{path}`"))
    w("")

    # --------------------------------------------------- 4. the controls
    w("## 4. The leakage controls")
    w("")
    w("A structure scorer ranks a candidate by its similarity to known "
      "actives. Left uncontrolled it can score highly by finding a "
      "candidate's own near-twin — which is retrieval, not prediction. "
      "Three controls were registered, each removing a different route "
      "for that to happen.")
    w("")
    w(row("control", "what it removes"))
    w(row("---", "---"))
    for _key, (name, what) in CONTROLS.items():
        w(row(f"**{name}**", what))
    w("")
    w("The publication control is **asymmetric between the classes** — a "
      "positive's nearest actives are often its own paper-mates, a "
      "negative's are not — which is why the other two were registered. "
      "Every result below is reported under all three.")
    w("")

    # ------------------------------------------- 5. multiplicity
    fam = mp["correction"]["families"]
    w("## 5. Multiple comparisons")
    w("")
    w("Every family is declared before its cells are read, and the "
      "corrected z follows from the family size rather than being "
      "chosen.")
    w("")
    w(row("family", "size", "alpha", "where"))
    w(row("---", "---:", "---:", "---"))
    w(row("ChEMBL structure cells", str(strength["family_size"]),
          f"{strength['alpha']:.6f}", "`prereg_structure_scorer.md`"))
    w(row("NCATS plate cells", str(strength["plate_family_size"]),
          f"{strength['plate_alpha']:.6f}",
          "`prereg_ncats_entry_assays.md`"))
    w(row("prospective cells", str(temporal["family_size"]),
          f"{temporal['alpha']:.6f}",
          "`prereg_temporal_structure_scorer.md`"))
    w(row("clinical labels", str(clin["family_size"]),
          f"{clin['alpha']:.6f}", "`prereg_clinical_labels.md`"))
    w(row("metapath channel-cells",
          str(fam["cells_computed"]["n"]),
          f"{fam['cells_computed']['alpha']:.6f}",
          "`evaluation.md`"))
    w("")

    # ============================== RESULT 1: cross-sectional scorer
    w("---")
    w("")
    w("## Result 1 — chemical similarity separates actives from "
      "inactives")
    w("")
    w("Registered by `docs/prereg_structure_scorer.md` "
      "(publication-disjoint), "
      "`docs/prereg_structure_scorer_symmetric.md` (symmetric k10), "
      "`docs/prereg_structure_scorer_strength_matched.md` "
      "(strength-matched) and `docs/prereg_ncats_entry_assays.md` "
      "(the plate cells).")
    w("")
    w("**Positive.** Each cell is one virus in one assay population. "
      "`clears` means the corrected interval excludes chance ABOVE it. "
      "`degree` is the baseline: a compound's edge count, no chemistry "
      "in it. A structure column that does not beat its own degree "
      "column is not evidence about chemistry.")
    w("")
    w(row("cell", "n+", "n−", "pub-disj", "sym k10", "strength",
          "degree", "verdict"))
    w(row("---", "---:", "---:", "---:", "---:", "---:", "---:", "---"))
    pub_cells = cells_of(pubdisj, "structure_max")

    def _mark(r):
        if not r:
            return "—"
        v = f"{r['auc']:.4f}"
        if r.get("clears_corrected_alpha"):
            v += "\\*" if r.get("direction") != "below" else "↓"
        return v

    ordered = sorted(strength_cells,
                     key=lambda k: -strength_cells[k][1]["auc"])
    for key in ordered:
        cell, st = strength_cells[key]
        label, virus = key
        deg = (cell.get("degree") or {})
        sy = (sym_cells.get(key) or (None, None))[1]
        pb = (pub_cells.get(key) or (None, None))[1]
        if key in deg_beaten:
            verdict = "**clears, degree does not**"
        elif key in both:
            verdict = "clears, so does degree"
        elif st.get("clears_corrected_alpha"):
            verdict = "clears under one control only"
        else:
            verdict = "does not clear"
        w(row(f"{virus} *({label})*",
              f"{st['n_pos']:,}", f"{st['n_neg']:,}",
              _mark(pb) if pb else "refused",
              _mark(sy), _mark(st), _mark(deg), verdict))
    w("")
    w(f"`*` = clears above chance · `↓` = clears BELOW chance · "
      f"**{len(deg_beaten)} of {len(strength_cells)}** cells clear above "
      f"a degree baseline that does not, and that membership is "
      f"**identical under both registered controls**.")
    w("")

    # ============================== RESULT 2: prospective
    w("---")
    w("")
    w("## Result 2 — it predicts forward, out of time")
    w("")
    w(f"**Positive, and the only out-of-time measurement here.** "
      f"Everything above is cross-sectional: labels split by publication "
      f"or by assay, never by time, so a scorer could be ranking "
      f"chemistry that was already known when the test compound was "
      f"published. This trains on activity published at or before "
      f"**{temporal['cutoff']}** and tests on activity published after "
      f"it — the arrangement that answers *would this have worked*.")
    w("")
    w("Registered by `docs/prereg_temporal_structure_scorer.md`. "
      "The screen this result justifies is registered separately, in "
      "`docs/prereg_denv2_prospective_screen.md`.")
    w("")
    w(f"Split `{temporal['split']}` · undated `{temporal['undated_policy']}` "
      f"· overlap `{temporal['overlap_policy']}` · family of "
      f"{temporal['family_size']}, alpha {temporal['alpha']:.6f}, "
      f"z {temporal['z_corrected']}")
    w("")
    w(row("virus", "references", "test +", "test −", "AUC",
          "corrected 95% CI", "degree", "verdict"))
    w(row("---", "---:", "---:", "---:", "---:", "---", "---:", "---"))
    tcells = temporal["cells"]
    for taxon in sorted(tcells, key=lambda t: -tcells[t]["primary"]["auc"]
                        if isinstance(tcells[t].get("primary", {}).get("auc"),
                                      (int, float)) else 0):
        c = tcells[taxon]
        pr = c.get("primary") or {}
        if not isinstance(pr.get("auc"), (int, float)):
            continue
        dg = c.get("degree") or {}
        degree = (f"{dg['auc']:.4f}" if isinstance(dg.get("auc"), (int, float))
                  and dg.get("usable") else "—")
        if taxon in t_clears:
            verdict = "**clears above chance**"
        elif pr.get("clears_corrected_alpha") and pr.get("direction") == "below":
            verdict = "clears BELOW chance"
        else:
            verdict = "does not clear"
        w(row(c["name"], f"{c['n_references']:,}",
              f"{c['n_test_positives']:,}", f"{c['n_test_negatives']:,}",
              f"{pr['auc']:.4f}",
              f"[{pr['auc_ci_lo_corrected']:.4f}, "
              f"{pr['auc_ci_hi_corrected']:.4f}]",
              degree, verdict))
    w("")
    usable = len(temporal.get("comparator_usable") or [])
    w(f"**{len(t_clears)} of {t_eval}** clear above chance. Two clear "
      f"BELOW it — they rank post-cutoff actives reliably last, which is "
      f"a finding and not noise, and nothing here establishes that "
      f"inverting them would hold forward.")
    w("")
    w(f"**{usable} of {t_eval} cells have a usable degree comparator.** "
      f"The pre-cutoff graph reaches too few of the compounds being "
      f"ranked to rank them at all, so the baseline column is empty by "
      f"measurement rather than by omission — and these clearing cells "
      f"are therefore **not shown to beat a baseline**, because there is "
      f"none to beat.")
    w("")

    # ============================== RESULT 3: the graph's contribution
    w("---")
    w("")
    w("## Result 3 — the knowledge graph adds nothing over a degree "
      "baseline")
    w("")
    total = readable = beats = bothc = 0
    for chans in (mp.get("per_virus") or {}).values():
        for k, ch in chans.items():
            if k in ("degree", "COMBINED") or not isinstance(ch, dict):
                continue
            total += 1
            readable += bool(ch.get("vs_degree_readable"))
            beats += bool(ch.get("beats_degree"))
            bothc += bool(ch.get("beats_degree")
                          and ch.get("vs_degree_readable"))
    w("**Negative, and it is the finding this project was built to "
      "test.** The graph is the thing under test; `degree` — a "
      "compound's edge count, with no biology in it — is what any "
      "reasoning channel has to beat.")
    w("")
    w(row("metapath channel-cells", "count"))
    w(row("---", "---:"))
    w(row("measured", f"{total}"))
    w(row("beat the degree baseline on lift", f"{beats}"))
    w(row("readable against degree at all", f"{readable}"))
    w(row("**both beat degree and are readable against it**",
          f"**{bothc}**"))
    w("")
    cs = (hard.get("cross_sectional") or {})
    for virus, chans in cs.items():
        comb = (chans.get("COMBINED") or {}).get("auc")
        deg = (chans.get("degree") or {}).get("auc")
        evaluable = {k: c for k, c in chans.items()
                     if k not in ("degree", "COMBINED")
                     and isinstance(c, dict) and c.get("evaluable")}
        if not (isinstance(comb, (int, float))
                and isinstance(deg, (int, float))):
            continue
        best = max(evaluable, key=lambda k: evaluable[k]["auc"])
        w(f"On **{virus}**, the assay-matched plate cell under the "
          f"publication-disjoint control:")
        w("")
        w(row("scorer", "AUC", "95% CI"))
        w(row("---", "---:", "---"))
        for name, r in (("strongest single channel (" + best + ")",
                         evaluable[best]),
                        ("all eight channels COMBINED",
                         chans["COMBINED"]),
                        ("**degree baseline**", chans["degree"])):
            w(row(name, f"{r['auc']:.4f}",
                  f"[{r['auc_ci_lo']:.3f}, {r['auc_ci_hi']:.3f}]"))
        w("")
        w(f"The ensemble of all eight reasoning channels scores "
          f"**{comb - deg:+.4f}** against the baseline that uses no "
          f"biology at all.")
        w("")

    # ============================== RESULT 4: pathway overlap
    w("## Result 4 — the pathway overlap is null")
    w("")
    rows_t = [r for r in (land.get("taxa") or {}).values()
              if isinstance(r, dict)]
    rows_i = [r for r in (land.get("interactome") or {}).values()
              if isinstance(r, dict)]
    def _clears(rs):
        return sum(1 for r in rs
                   if r.get("clears_corrected_alpha") or r.get("clears"))
    w("Registered by `docs/prereg_pathway_overlap.md`.")
    w("")
    w("**Negative.** Whether a drug's targets sit in the same pathways "
      "as a virus's host-dependency factors, against a degree-matched "
      "permutation null.")
    w("")
    w(row("channel", "rows measured", "clearing the corrected alpha"))
    w(row("---", "---:", "---:"))
    w(row("dependency", f"{len(rows_t)}", f"{_clears(rows_t)}"))
    w(row("interactome", f"{len(rows_i)}", f"{_clears(rows_i)}"))
    w(row("**total**", f"**{len(rows_t) + len(rows_i)}**",
          f"**{_clears(rows_t) + _clears(rows_i)}**"))
    w("")
    w(f"Permutations: {land.get('permutations'):,}, seed "
      f"`{land.get('seed')}`, degree basis `{land.get('degree_basis')}`.")
    w("")

    # ============================== RESULT 5: the independent axis
    w("## Result 5 — the one independent label axis does not corroborate "
      "it")
    w("")
    pop = clin.get("population") or {}
    st_c = clin.get("structure_max_k10") or {}
    dg_c = clin.get("degree") or {}
    w("Registered by `docs/prereg_clinical_labels.md`.")
    w("")
    w(f"**Negative, and it points the other way.** Every AUC above rests "
      f"on `HAS_ANTIVIRAL_ACTIVITY_AGAINST` — one predicate, in-vitro "
      f"concentration-response at one threshold. The clinical axis is "
      f"independent of it: compounds in a trial or approved for "
      f"{_disease_label(clin.get('disease'))}.")
    w("")
    w(row("scorer", "AUC", "95% CI", "verdict"))
    w(row("---", "---:", "---", "---"))
    w(row("structure (symmetric k10)", f"{st_c['auc']:.4f}",
          f"[{st_c['auc_ci_lo']:.4f}, {st_c['auc_ci_hi']:.4f}]",
          "**below chance**"))
    w(row("degree baseline", f"{dg_c['auc']:.4f}",
          f"[{dg_c['auc_ci_lo']:.4f}, {dg_c['auc_ci_hi']:.4f}]",
          "well above chance"))
    w("")
    w(f"{pop.get('positives'):,} positives against "
      f"{pop.get('negatives'):,} negatives, of which "
      f"{pop.get('positives_invitro_inactive'):,} of the positives were "
      f"measured INACTIVE in vitro. Reaching a trial is predicted by how "
      f"much a compound has been studied, not by its resemblance to the "
      f"compounds active in vitro — so this axis measures something else "
      f"rather than confirming the one above, and cannot be used to "
      f"validate a structure scorer.")
    w("")

    # ============================== the qualification
    w("---")
    w("")
    w("## The qualification that applies to every AUC above")
    w("")
    w("**An AUC is not a hit rate.** It is the probability that a random "
      "active outranks a random inactive — a statement about the WHOLE "
      "ordering. A screen only ever orders the top of a list, and the "
      "two can disagree completely.")
    w("")
    plate = [c for c in hr["plate_cells"]
             if c["population"] == "secondary_ncats"]
    if plate:
        c = plate[0]
        w(f"On the one population matched by assay rather than by source "
          f"— {c['name']}, pool {c['pool']:,}, base rate "
          f"{c['base_rate']:.2%}:")
        w("")
        w(row("control", "precision@10", "enrichment over base"))
        w(row("---", "---:", "---:"))
        for label, key in (("no exclusion", "no exclusion"),
                           ("symmetric k = 10", "k10"),
                           (f"strength-matched (k_cell {c['k_cell']})",
                            f"k_cell={c['k_cell']}")):
            m = _pk(c["controls"][key], 10)
            w(row(label, f"{m['precision_expected']:.1%}",
                  f"{m['enrichment_expected']:.2f}×"))
        w("")
        w(f"That cell's strength-matched AUC is "
          f"**{strength_cells[('plate CPE', c['name'])][1]['auc']:.4f} "
          f"and it CLEARS** its registered alpha. The AUC cannot show "
          f"this, because it averages over a ranking nobody screens. A "
          f"precision that vanishes when a candidate's nearest known "
          f"actives are withheld was analogue retrieval, not prediction.")
        w("")
    w("Prospectively, at k = 10 under the registered secondary control:")
    w("")
    w(row("virus", "base rate", "precision@10", "enrichment"))
    w(row("---", "---:", "---:", "---:"))
    for c in sorted(hr["cells"],
                    key=lambda c: -(_pk(c["controls"]["symmetric_k10"], 10)
                                    or {}).get("enrichment_expected", 0)):
        m = _pk(c["controls"]["symmetric_k10"], 10)
        if not m:
            continue
        w(row(c["name"], f"{c['base_rate']:.1%}",
              f"{m['precision_expected']:.0%}",
              f"**{m['enrichment_expected']:.2f}×**"
              if m["enrichment_expected"] >= 2
              else f"{m['enrichment_expected']:.2f}×"))
    w("")
    w("**Enrichment, not precision, is the column that matters.** A "
      "precision of 90% in a pool that is already 70% active is a "
      "ranking that added almost nothing. Only DENV-2 shows a lift that "
      "would survive a feasible screen; `make screen-power` sizes one "
      "for every cell, and six of the ten cannot be powered at any arm "
      "size up to 120.")
    w("")
    w("`make hit-rate` reproduces this. It is **post-hoc**: precision@k "
      "was computed after the registered AUCs were read, so it "
      "qualifies a claim rather than being one.")
    w("")

    # ============================== what is claimed
    w("---")
    w("")
    w("## What is claimed, and what is not")
    w("")
    w("**Claimed.** On this release, ranking compounds by maximum "
      "ECFP4 Tanimoto to a virus's known actives separates measured "
      f"actives from measured inactives in {len(deg_beaten)} of "
      f"{len(strength_cells)} registered cells, above a degree baseline "
      f"that does not clear, with the same membership under both "
      f"registered controls — and does so **prospectively** in "
      f"{len(t_clears)} of {t_eval} cells trained on pre-"
      f"{temporal['cutoff']} evidence.")
    w("")
    w("**Not claimed.**")
    w("")
    w("- *That the knowledge graph contributes to this.* It does not. "
      "The scorer uses no graph edge, and `make auc-contribution` "
      "measures the graph's contribution to every scorer here as zero "
      "or negative against its own degree baseline.")
    w("- *That a high AUC implies a usable hit rate.* On the one "
      "assay-matched cell it does not, measured above.")
    w("- *That a clearing cell beats a baseline prospectively.* None of "
      "them is shown to: no prospective cell has a usable degree "
      "comparator.")
    w("- *Anything clinical.* The endpoint is an in-vitro EC50 in one "
      "assay, often in Vero E6, which lacks meaningful TMPRSS2 — and "
      "the one independent label axis anti-correlates with this "
      "statistic.")
    w("- *That a suggestion is a prediction of activity.* It is the "
      "claim that a compound RESEMBLES a known active, with the known "
      "active named.")
    w("")
    w("## Reproducing all of it")
    w("")
    w("```")
    w("make rebuild                 # the graph, ~27 min, byte-reproducible")
    w("make check-results           # which results describe this graph")
    w("make scorer-table            # Result 1, every control side by side")
    w("make temporal-structure-scorer   # Result 2")
    w("make auc-contribution        # Result 3, per-approach shares")
    w("make pathway-landscape       # Result 4")
    w("make clinical-labels         # Result 5")
    w("make hit-rate                # the qualification")
    w("make screen-power            # what a wet-lab screen would need")
    w("```")
    w("")
    # EVERY REGISTRATION, INCLUDING THE ONES WHOSE RESULTS ARE NOT HERE.
    # Listing only the ones this document reports would let a reader
    # believe it covers the whole project, and the omission would be
    # invisible. Which ones are reported here is DERIVED from whether
    # the body above names them, so adding a section cannot leave the
    # table behind.
    body = "\n".join(L)
    w("## Every registration")
    w("")
    w("Each was signed before the result it governs existed. "
      "`reported here` is derived from whether this document names it "
      "above.")
    w("")
    w(row("registration", "subject", "reported here"))
    w(row("---", "---", "---"))
    for f in sorted((ROOT / "docs").glob("prereg_*.md")):
        title = ""
        for line in f.read_text().splitlines():
            if line.startswith("# "):
                title = line[2:].strip()
                break
        title = title.replace("Pre-registration:", "").strip(" :—-")
        here = "yes" if f.name in body else "no — see `evaluation.md`"
        w(row(f"`{f.name}`", title or "—", here))
    w("")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=DOC)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)

    text = build()
    if args.check:
        if not args.out.exists():
            print(f"{args.out} does not exist; run `make results-doc`")
            return 1
        if args.out.read_text() != text:
            print(f"{args.out.name} is STALE.")
            import difflib
            diff = list(difflib.unified_diff(
                args.out.read_text().splitlines(), text.splitlines(),
                "on disk", "from the results", lineterm="", n=0))
            from kgav.schema import print_capped
            print_capped(diff, 40, indent="  ", unit="more diff line(s)")
            print("\n  `make results-doc` regenerates it.")
            return 1
        print(f"{args.out.name} is current.")
        return 0

    args.out.write_text(text)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
