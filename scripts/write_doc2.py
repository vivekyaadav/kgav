#!/usr/bin/env python3
"""Generate doc2/: the whole project described from the code and the artifacts.

    python scripts/write_doc2.py            # write doc2/*.md
    python scripts/write_doc2.py --check    # fail if any file is out of date

WHY GENERATED AND NOT WRITTEN. doc2 is a from-scratch description of this
repository, and the requirement on it was "nothing should be wrong". A
hand-written document cannot meet that for longer than one rebuild: this
project has already shipped a hand-typed 0.7841 for a cell reading 0.7822, a
"none below chance" for a family with two below, thirteen stale figures in a
hand-drawn metagraph, and a pool size of 21,081 that a filter moved to 21,039
the same day. Every one of those was true when typed.

So the prose lives in doc2/_templates/*.md.in with `{{placeholders}}`, and
every figure is computed here from:

    schema/kg_schema.yaml          the declared graph model
    config/viruses.yaml            the virus register
    data/releases/v0.1/MANIFEST.json   what was actually built
    data/results                   every registered run, read never written
    src/kgav/*.py                  the constants, imported not copied
    tests/, scripts/, Makefile     counted, not estimated

IT WRITES NOTHING UNDER data/results, and the docstring above names that
directory without a file glob on purpose. tests/test_provenance.py flags any
script that writes a file AND carries a string that looks like a path to a
stamped artifact, because it cannot tell a path a script reads from one it
writes -- and the match is over the WHOLE string constant, so naming the
suffix anywhere in this docstring would trip it from thirty lines away. That
conservatism is the right direction: the rule exists because one script
hand-rolled its own stamp while write_results' docstring claimed "every
results writer goes through here".

THREE RULES, each of which this repository learned the hard way:

  A PLACEHOLDER WITH NO VALUE IS AN ERROR, not an empty string. A figure
  that renders blank reads as prose that was always meant to read that way.

  A VALUE NO TEMPLATE USES IS REPORTED. A derivation nothing renders is
  dead code that looks like coverage -- the same defect as a check that
  cannot fail.

  --check DIFFS, and names the first line that moved. `make doc2-check` is
  the half that belongs in a gate.

WHAT THIS DOES NOT DO. It does not reproduce docs/prereg_*.md. Those are
signed pre-registrations, amended by appending and never rewritten; doc2
describes them, quotes their verdicts from the results files they registered,
and links to them. Rewriting a signed registration from scratch would destroy
the one property that makes it one.
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav import virus_register as vr                      # noqa: E402
from kgav.agent import guards                              # noqa: E402
from kgav.agent import tools as agent_tools                # noqa: E402
from kgav.assemble import LAYER_PRECEDENCE                 # noqa: E402
from kgav.calibration import (                             # noqa: E402
    CHANCE, MIN_REACHED_POSITIVES)
from kgav.labels import (                                  # noqa: E402
    DEFAULT_ACTIVE_BELOW_NM, DEFAULT_INACTIVE_ABOVE_NM,
    DEFAULT_DOSE_SERIES_POLICY, DOSE_SERIES_POLICIES)
from kgav.landscape import ALPHA as LANDSCAPE_ALPHA        # noqa: E402
from kgav.pipeline import LAYERS, release_dir              # noqa: E402
from kgav.schema import load_schema                        # noqa: E402
from kgav.similarity import FP_BITS, FP_RADIUS             # noqa: E402

OUT = ROOT / "doc2"
TEMPLATES = OUT / "_templates"
# THE RELEASE PATH HAS ONE OWNER, and this file spelling the literal is
# what tests/test_release_path_has_one_owner.py exists to catch -- it
# found this line within a minute of the file being written, having been
# built after the same literal was found at 39 sites.
RELEASE = release_dir(ROOT)
RESULTS = ROOT / "data" / "results"
PLACEHOLDER = re.compile(r"\{\{([a-z0-9_]+)\}\}")


# --------------------------------------------------------------- helpers
def _j(path: Path) -> dict:
    return json.loads(path.read_text())


def table(header: list[str], rows: list[list[str]],
          align: str | None = None) -> str:
    """A markdown table. `align` is one character per column: l, r or c.

    RENDERED HERE RATHER THAN IN THE TEMPLATE because a row count is a
    figure like any other: docs/ has a test of its own
    (test_docs_tables_render.py) for the case where a cell past the
    header's column count is silently not displayed at all.
    """
    align = align or "l" * len(header)
    assert len(align) == len(header), (header, align)
    bar = {"l": ":---", "r": "---:", "c": ":---:"}
    lines = ["| " + " | ".join(header) + " |",
             "|" + "|".join(bar[a] for a in align) + "|"]
    for r in rows:
        assert len(r) == len(header), f"row {r} has {len(r)} of {len(header)}"
        lines.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(lines)


def n(x) -> str:
    return f"{x:,}" if isinstance(x, int) else str(x)


def pct(x: float, places: int = 1) -> str:
    return f"{100 * x:.{places}f}%"


def auc(x) -> str:
    return "—" if x is None else f"{x:.4f}"


# ------------------------------------------------------------- the values
def derive() -> dict[str, str]:
    v: dict[str, str] = {}
    schema = load_schema()
    reg = vr.load()
    man = _j(RELEASE / "MANIFEST.json")
    names = {f"NCBITaxon:{e['taxon']}": e["name"] for e in reg["viruses"]}

    # ---------------------------------------------------- release census
    v["schema_version"] = schema.version
    v["release_id"] = man["release"]
    v["n_nodes"] = n(man["nodes"])
    v["n_edges"] = n(man["edges"])
    v["nodes_sha"] = man["layers"]["spine"]["nodes_sha"][:0] or _release_sha("nodes")
    v["edges_sha"] = _release_sha("edges")
    v["build_commit"] = man["code"]["commit"][:12]
    v["build_dirty"] = str(man["code"]["dirty"])
    v["hub_cut_degree"] = n(man["hub_cut_degree"])
    v["degree_median"] = n(man["degree"]["median"])
    v["degree_p99"] = n(man["degree"]["p99"])
    v["degree_p99_9"] = n(man["degree"]["p99_9"])
    v["degree_max"] = n(man["degree"]["max"])

    v["nodes_by_class_table"] = table(
        ["node class", "nodes", "id prefixes"],
        [[c, n(k), ", ".join(schema.node_classes[c]["id_prefixes"])]
         for c, k in sorted(man["nodes_by_class"].items(), key=lambda x: -x[1])],
        "lrl")

    held = schema.held_out_predicates()
    labelp = schema.evaluation_label_predicates()
    sym = schema.symmetry()
    rows = []
    for p, k in sorted(man["edges_by_predicate"].items(), key=lambda x: -x[1]):
        ec = next((e for e in schema.edge_classes if e.predicate == p), None)
        flags = []
        if p in held:
            flags.append("held out")
        if p in labelp:
            flags.append("label")
        if sym.get(p):
            flags.append("symmetric")
        rows.append([f"`{p}`", n(k),
                     f"{ec.subject} → {ec.object}" if ec else "?",
                     ", ".join(flags) or "—"])
    v["edges_by_predicate_table"] = table(
        ["predicate", "edges", "subject → object", "flags"], rows, "lrll")

    v["n_predicates_built"] = str(len(man["edges_by_predicate"]))
    v["n_predicates_declared"] = str(len({e.predicate for e in schema.edge_classes}))
    v["n_node_classes_declared"] = str(len(schema.node_classes))
    v["n_node_classes_built"] = str(len(man["nodes_by_class"]))
    v["held_out_predicates"] = ", ".join(f"`{p}`" for p in sorted(held))
    v["label_predicates"] = ", ".join(f"`{p}`" for p in sorted(labelp))
    v["n_prefixes"] = str(len(schema.prefixes))
    v["n_enums"] = str(len(schema.enums))

    # metapaths
    rows = []
    for code, mp in sorted(schema.metapaths.items()):
        path = mp.get("path") or []
        pretty = " → ".join(path)
        cons = mp.get("constraints") or {}
        rows.append([f"**{code}**", mp.get("doc", ""), f"`{pretty}`",
                     "; ".join(f"`{k}` = `{x}`" for k, x in cons.items()) or "—",
                     "yes" if mp.get("exclude_hubs") else "no"])
    v["metapath_table"] = table(
        ["code", "what it is", "path", "constraints", "hub exclusion"], rows)
    v["n_metapaths"] = str(len(schema.metapaths))
    rl = schema.retrieval_limits
    v["max_pathway_members"] = n(rl["max_pathway_members"])
    v["max_neighbors_per_node"] = n(rl["max_neighbors_per_node"])
    v["max_subgraph_edges"] = n(rl["max_subgraph_edges"])
    v["hub_exclude_classes"] = ", ".join(rl["hub_exclude_classes"])

    # ------------------------------------------------------- the viruses
    rows = []
    for e in reg["viruses"]:
        curie = f"NCBITaxon:{e['taxon']}"
        cs = e.get("chain_source") or {}
        rows.append([
            e["name"], curie, vr.family_of(e, reg) or "?",
            cs.get("accession") or "—",
            str(cs.get("chains") if cs.get("chains") is not None else "—"),
            "yes" if str(cs.get("taxon")) != str(e["taxon"]) else "no",
            e.get("mondo") or "—"])
    v["virus_table"] = table(
        ["name", "taxon", "family", "polyprotein", "chains",
         "chains from a strain", "disease"], rows, "lllrrll")
    v["n_viruses_built"] = str(len(reg["viruses"]))
    v["n_viruses_surveyed"] = str(len(reg.get("surveyed") or []))
    v["n_register_violations"] = str(len(vr.validate(reg)))
    fams = sorted(set(vr.families(reg).values()))
    v["virus_families"] = (" and ".join(fams) if len(fams) < 3
                           else ", ".join(fams[:-1]) + " and " + fams[-1])

    # --------------------------------------------------------- the layers
    rows = []
    for layer in LAYERS:
        src = (man["layers"].get(layer.name) or {}).get("sources")
        inputs = list(layer.files) + list(layer.globs) + [
            " or ".join(a) for a in layer.any_globs]
        rows.append([
            layer.name, f"`scripts/{layer.script}`", str(layer.phase),
            "required" if layer.required else "optional",
            "peers" if layer.reads_peers else "raw",
            ", ".join(f"`{x}`" for x in inputs) or "—",
            _sources_text(src)])
    v["layer_table"] = table(
        ["layer", "script", "phase", "status", "reads", "raw inputs",
         "source version"], rows)
    v["n_layers"] = str(len(LAYERS))
    # DECLARED IS NOT BUILT, and on this release they differ by one. The
    # first draft of doc2 said "assembled from 13 layers" because the
    # pipeline declares 13; the MANIFEST holds 12. `curated` is skipped
    # -- correctly, and by name -- because its raw input directory is
    # empty, and the flavivirus host factors it exists to supply are
    # therefore absent from this graph. A document that counted the
    # declaration would have overstated the build.
    from kgav.pipeline import availability
    built = set(man["layers"])
    absent = [lyr for lyr in LAYERS if lyr.name not in built]
    v["n_layers_built"] = str(len(built))
    v["layers_absent_why"] = "; ".join(
        f"`{lyr.name}` — {availability(lyr, ROOT)[1] or 'no reason recorded'}"
        for lyr in absent) or "none absent"
    v["layer_precedence"] = " > ".join(LAYER_PRECEDENCE)
    v["layers_no_version"] = ", ".join(
        man["sources"]["layers_declaring_no_source_version"])

    merges = man["merges"]
    v["n_edges_merged"] = n(merges["edge_merged"])
    v["merge_table"] = table(
        ["merge statistic", "count"],
        [[f"`{k}`", n(x)] for k, x in sorted(merges.items(),
                                             key=lambda kv: -kv[1])], "lr")
    v["n_precedence_decided"] = n(sum(
        x for k, x in merges.items()
        if k.startswith("qualifier_precedence_decided")))
    v["n_irreconcilable"] = n(sum(
        x for k, x in merges.items()
        if k.startswith("qualifier_irreconcilable")))

    # --------------------------------------------------------- the labels
    v["inactive_above_nm"] = n(int(DEFAULT_INACTIVE_ABOVE_NM))
    v["active_below_nm"] = n(int(DEFAULT_ACTIVE_BELOW_NM))
    v["dose_series_default"] = DEFAULT_DOSE_SERIES_POLICY
    v["dose_series_policies"] = ", ".join(f"`{p}`" for p in DOSE_SERIES_POLICIES)
    v["fp_bits"] = n(FP_BITS)
    v["fp_radius"] = str(FP_RADIUS)
    v["chance"] = str(CHANCE)
    v["min_reached_positives"] = str(MIN_REACHED_POSITIVES)
    v["landscape_alpha"] = str(LANDSCAPE_ALPHA)

    # ------------------------------------------------- the scorer runs
    for tag, fname in (("pub", "structure_scorer.json"),
                       ("sym", "structure_scorer_symmetric.json"),
                       ("str", "structure_scorer_strength_matched.json")):
        d = _j(RESULTS / fname)
        v[f"{tag}_control"] = d["control"]
        v[f"{tag}_family"] = str(d["family_size"])
        v[f"{tag}_alpha"] = f"{d['alpha']:.6f}"
        v[f"{tag}_z"] = str(d["z_corrected"])
        v[f"{tag}_n_labelled"] = n(d["n_labelled"])
        rows, clears, below = [], [], []
        for curie, cell in d["primary_chembl"].items():
            m = cell.get("structure_max") or {}
            if m.get("auc") is None:
                continue
            nm = names.get(curie, curie)
            ok = bool(m.get("clears_corrected_alpha")) and m.get("direction") != "below"
            if ok:
                clears.append(nm)
            if m.get("direction") == "below" and m.get("clears_corrected_alpha"):
                below.append(nm)
            # THREE STATES, NOT TWO, and the first version of this
            # printed two. A cell whose AUC is under 0.5 and whose
            # interval still SPANS 0.5 has no verdict; a cell whose
            # interval EXCLUDES 0.5 on the wrong side is a measured
            # anti-ranking. Rendering both as "below chance" made this
            # document report four anti-rankings where the run records
            # two, in the same table as the count that says two.
            verdict = ("**clears**" if ok else
                       "**clears BELOW chance**"
                       if (m.get("direction") == "below"
                           and m.get("clears_corrected_alpha"))
                       else "no verdict")
            rows.append([nm, n(m["n_pos"]), n(m["n_neg"]), f"{m['auc']:.4f}",
                         f"[{m['auc_ci_lo_corrected']:.4f}, "
                         f"{m['auc_ci_hi_corrected']:.4f}]", verdict])
        rows.sort(key=lambda r: -float(r[3]))
        v[f"{tag}_table"] = table(
            ["virus", "active", "inactive", "AUC",
             f"corrected {100 * (1 - d['alpha']):.3f}% interval", "verdict"],
            rows, "lrrllc")
        v[f"{tag}_n_clears"] = str(len(clears))
        v[f"{tag}_clears"] = ", ".join(sorted(clears))
        v[f"{tag}_n_below"] = str(len(below))
        v[f"{tag}_below"] = ", ".join(sorted(below)) or "none"
        v[f"{tag}_n_cells"] = str(len(rows))
    # TWO DEFINITIONS OF "BELOW", AND THE READ-OUT GATES ON THE OTHER ONE.
    # `{tag}_below` above is the RUN's verdict: the corrected interval
    # excludes chance on the wrong side. `repurposing_candidates.tier`
    # is stricter -- `auc < 0.5 OR direction == below` -- because a
    # ranking of unlabelled compounds from a cell whose point estimate is
    # under chance is not something to print under a caveat. DENV-4 is
    # the cell where the two differ, and conflating them is how a
    # document comes to disagree with the tool it describes.
    sm = _j(RESULTS / "structure_scorer_strength_matched.json")
    from repurposing_candidates import tier as _tier_fn                # noqa
    tiers = {}
    for curie, cell in sm["primary_chembl"].items():
        row = cell.get("structure_max") or {}
        if row.get("auc") is None:
            continue
        tiers[names.get(curie, curie)] = _tier_fn({
            "auc": row["auc"], "direction": row.get("direction"),
            "clears_corrected_alpha": row.get("clears_corrected_alpha")})
    v["str_tier_table"] = table(
        ["virus", "AUC", "interval verdict", "read-out tier"],
        [[nm, f"{(sm['primary_chembl'][c].get('structure_max') or {})['auc']:.4f}",
          ("clears" if (sm["primary_chembl"][c]["structure_max"].get(
              "clears_corrected_alpha")
              and sm["primary_chembl"][c]["structure_max"].get("direction")
              != "below")
           else "clears BELOW chance"
           if sm["primary_chembl"][c]["structure_max"].get(
               "clears_corrected_alpha") else "no verdict"),
          f"`{tiers[nm]}`"]
         for c, nm in sorted(((c, names.get(c, c)) for c in sm["primary_chembl"]
                              if (sm["primary_chembl"][c].get("structure_max")
                                  or {}).get("auc") is not None),
                             key=lambda cn: -sm["primary_chembl"][cn[0]][
                                 "structure_max"]["auc"])], "llll")
    v["str_n_tier_below"] = str(sum(1 for t in tiers.values() if t == "below"))
    v["str_tier_below"] = ", ".join(
        sorted(k for k, t in tiers.items() if t == "below")) or "none"
    v["str_n_tier_clears"] = str(sum(1 for t in tiers.values() if t == "clears"))
    v["str_n_tier_weak"] = str(sum(1 for t in tiers.values() if t == "weak"))
    v["k_cell_table"] = table(
        ["virus", "k_cell"],
        [[names.get(k, k), str(x)] for k, x in sorted(
            sm["k_cells"].items(), key=lambda kv: -kv[1])], "lr")
    v["sym_k"] = str(_j(RESULTS / "structure_scorer_symmetric.json")["k_primary"])

    # plate cells (the assay-matched secondary populations)
    rows = []
    for tag, fname in (("publication-disjoint", "structure_scorer.json"),
                       ("symmetric k=10", "structure_scorer_symmetric.json"),
                       ("strength-matched", "structure_scorer_strength_matched.json")):
        d = _j(RESULTS / fname)
        for pop, label in (("secondary_ncats", "NCATS CPE"),
                           ("secondary_ncats_entry", "NCATS entry")):
            for curie, cell in (d.get(pop) or {}).items():
                m = cell.get("structure_max") or {}
                if m.get("auc") is None:
                    continue
                rows.append([tag, label, names.get(curie, curie),
                             n(m["n_pos"]), n(m["n_neg"]), f"{m['auc']:.4f}",
                             "**clears**"
                             if m.get("clears_corrected_alpha")
                             and m.get("direction") != "below"
                             else ("**clears BELOW chance**"
                                   if m.get("clears_corrected_alpha")
                                   else "no verdict")])
    v["plate_table"] = table(
        ["control", "population", "virus", "active", "inactive", "AUC",
         "verdict"], rows, "lllrrlc")
    v["plate_family"] = str(_j(RESULTS / "structure_scorer.json")["plate_family_size"])
    v["plate_z"] = str(_j(RESULTS / "structure_scorer.json")["plate_z_corrected"])

    # ------------------------------------------------------- prospective
    t = _j(RESULTS / "temporal_structure_scorer.json")
    v["temporal_split"] = t["split"]
    v["temporal_cutoff"] = str(t["cutoff"])
    v["temporal_undated"] = t["undated_policy"]
    v["temporal_overlap"] = t["overlap_policy"]
    v["temporal_family"] = str(t["family_size"])
    v["temporal_alpha"] = str(t["alpha"])
    v["temporal_z"] = str(t["z_corrected"])
    v["temporal_n_evaluable"] = str(t["n_evaluable"])
    v["temporal_n_candidates"] = str(t["n_cells_total"])
    v["temporal_min_comparator_reach"] = pct(t["min_comparator_reach"], 0)
    rows = []
    for curie, c in sorted(t["cells"].items(), key=lambda kv: kv[1]["name"]):
        p = c.get("primary") or {}
        s = c.get("symmetric_k10") or {}
        dg = c.get("degree") or {}
        state = ("**clears**" if curie in t["clears"] else
                 "**clears BELOW chance**"
                 if p.get("clears_corrected_alpha")
                 and p.get("direction") == "below" else "no verdict")
        rows.append([c["name"], n(c["n_references"]), n(c["n_test_positives"]),
                     n(c["n_test_negatives"]), f"{p['auc']:.4f}",
                     f"[{p['auc_ci_lo_corrected']:.4f}, "
                     f"{p['auc_ci_hi_corrected']:.4f}]",
                     f"{s['auc']:.4f}" if s.get("auc") is not None else "—",
                     auc(dg.get("auc")), state])
    v["temporal_table"] = table(
        ["virus", "references", "test +", "test −", "AUC",
         "corrected interval", "k=10", "degree", "verdict"], rows,
        "lrrrllllc")
    v["temporal_clears"] = ", ".join(sorted(names.get(c, c) for c in t["clears"]))
    v["temporal_n_clears"] = str(len(t["clears"]))
    below = [c["name"] for k, c in t["cells"].items()
             if (c.get("primary") or {}).get("direction") == "below"
             and (c.get("primary") or {}).get("clears_corrected_alpha")]
    v["temporal_below"] = ", ".join(sorted(below))
    v["temporal_n_below"] = str(len(below))
    v["temporal_comparator_usable"] = str(len(t["comparator_usable"]))
    # THE ONE CELL THAT DISAGREES WITH ITSELF ACROSS PROTOCOLS, derived
    # from both files rather than typed into the prose that explains it.
    d4 = "NCBITaxon:11070"
    v["denv4_name"] = names.get(d4, d4)
    v["denv4_cross_auc"] = f"{sm_rows['primary_chembl'][d4]['structure_max']['auc']:.4f}" \
        if (sm_rows := _j(RESULTS / "structure_scorer_strength_matched.json")) \
        else "—"
    v["denv4_prospective_auc"] = f"{t['cells'][d4]['primary']['auc']:.4f}"
    v["denv4_k_cell"] = str(sm_rows["k_cells"][d4])

    # --------------------------------------------------- clinical labels
    c = _j(RESULTS / "clinical_labels.json")
    v["clinical_disease"] = c["disease"]
    v["clinical_virus"] = names.get(c["virus"], c["virus"])
    v["clinical_family"] = str(c["family_size"])
    v["clinical_alpha"] = str(c["alpha"])
    v["clinical_treats"] = n(c["clinical_edges"]["TREATS"])
    v["clinical_trials"] = n(c["clinical_edges"]["IN_TRIAL_FOR"])
    p = c["population"]
    v["clinical_pop_table"] = table(
        ["population", "compounds"],
        [[k.replace("_", " "), n(x)] for k, x in p.items()], "lr")
    v["clinical_structure_auc"] = f"{c['structure_max_k10']['auc']:.4f}"
    v["clinical_structure_ci"] = (
        f"[{c['structure_max_k10']['auc_ci_lo_corrected']:.4f}, "
        f"{c['structure_max_k10']['auc_ci_hi_corrected']:.4f}]")
    v["clinical_structure_dir"] = c["structure_max_k10"]["direction"]
    v["clinical_degree_auc"] = f"{c['degree']['auc']:.4f}"
    v["clinical_degree_ci"] = (f"[{c['degree']['auc_ci_lo_corrected']:.4f}, "
                               f"{c['degree']['auc_ci_hi_corrected']:.4f}]")
    v["clinical_degree_wins"] = str(c["p2_degree_at_or_above_structure"])

    # ---------------------------------------------------- hard negatives
    h = _j(RESULTS / "hard_negatives.json")
    v["hn_label_source"] = ", ".join(h["label_source"])
    v["hn_assay_type"] = ", ".join(h["label_assay_type"])
    v["hn_publication_disjoint"] = str(h["publication_disjoint"])
    v["hn_selectivity"] = f"{h['selectivity_filter']} / {h['selectivity_applies']}"
    v["hn_inactive_above"] = n(int(h["inactive_above_nm"]))
    rows = []
    for virus, chans in h["cross_sectional"].items():
        for name_, row in chans.items():
            if not isinstance(row, dict) or row.get("auc") is None:
                continue
            rows.append([virus, name_, n(row["n_pos"]), n(row["n_neg"]),
                         f"{row['auc']:.4f}",
                         f"[{row['auc_ci_lo']:.4f}, {row['auc_ci_hi']:.4f}]",
                         pct(row["coverage_pos"]),
                         f"{row['auc_floor']:.4f}",
                         "—" if row.get("auc_skill") is None
                         else f"{row['auc_skill']:.3f}"])
    v["hn_table"] = table(
        ["virus", "channel", "+", "−", "AUC", "95% interval",
         "positives reached", "floor", "skill"], rows, "llrrllrrr")
    v["hn_temporal"] = ("skipped — the section is empty in the artifact"
                        if not h["temporal"] else "present")
    # WHICH SCORERS CLEAR, counted rather than read off the table by a
    # reader. The first version of 01-overview said "none has an interval
    # excluding chance" and TWO do: the ensemble and the degree baseline.
    # The true statement is about the eight metapath CHANNELS, and the
    # distinction is the finding.
    metapaths = [c for c in ("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")]
    # THE CEILING AND FLOOR EXAMPLES, DERIVED. kgav.labels' docstrings
    # quote these from the release they were written on, and doc2 quoted
    # the docstring -- a figure one release away from being wrong with
    # full confidence, which is the defect this folder exists to avoid.
    worst = min(
        ((c, (h["cross_sectional"]["SARS-CoV-2"].get(c) or {}))
         for c in metapaths),
        key=lambda cx: cx[1].get("coverage_pos", 1.0))
    v["hn_thinnest_channel"] = worst[0]
    v["hn_thinnest_auc"] = f"{worst[1]['auc']:.4f}"
    v["hn_thinnest_ceiling"] = f"{worst[1]['auc_ceiling']:.4f}"
    v["hn_thinnest_cov"] = pct(worst[1]["coverage_pos"])
    v["hn_thinnest_reached"] = str(round(
        worst[1]["coverage_pos"] * worst[1]["n_pos"]))
    v["hn_n_pos"] = n(worst[1]["n_pos"])
    on_ceiling = [c for c in metapaths
                  if abs((h["cross_sectional"]["SARS-CoV-2"][c]["auc"]
                          - h["cross_sectional"]["SARS-CoV-2"][c]
                          ["auc_ceiling"])) < 5e-4]
    v["hn_on_ceiling"] = ", ".join(on_ceiling) or "none"
    v["hn_n_on_ceiling"] = str(len(on_ceiling))
    best = max(metapaths,
               key=lambda c: h["cross_sectional"]["SARS-CoV-2"][c]["auc"])
    v["hn_best_channel"] = best
    v["hn_best_channel_auc"] = \
        f"{h['cross_sectional']['SARS-CoV-2'][best]['auc']:.4f}"
    below_ceiling = [c for c in metapaths
                     if h["cross_sectional"]["SARS-CoV-2"][c]["auc_ceiling"]
                     < CHANCE]
    v["hn_ceiling_below_chance"] = ", ".join(below_ceiling) or "none"
    v["hn_combined_margin_share"] = pct(
        (cells_c := h["cross_sectional"]["SARS-CoV-2"]["COMBINED"])
        and (cells_c["auc_floor"] - CHANCE) / (cells_c["auc"] - CHANCE), 0)
    cells = h["cross_sectional"]["SARS-CoV-2"]
    clearing = [c for c in metapaths
                if (cells.get(c) or {}).get("auc_ci_lo", 0) > CHANCE]
    v["hn_n_channels"] = str(len(metapaths))
    v["hn_channels_clearing"] = ", ".join(clearing) or "none"
    v["hn_n_channels_clearing"] = str(len(clearing))
    for name_, key in (("degree", "degree"), ("COMBINED", "combined")):
        row = cells[name_]
        v[f"hn_{key}_auc"] = f"{row['auc']:.4f}"
        v[f"hn_{key}_ci"] = f"[{row['auc_ci_lo']:.4f}, {row['auc_ci_hi']:.4f}]"
        v[f"hn_{key}_floor"] = f"{row['auc_floor']:.4f}"
        v[f"hn_{key}_skill"] = pct(row["auc_skill"])
        v[f"hn_{key}_cov"] = pct(row["coverage_pos"])

    # ---------------------------------------------- metapath comparison
    mc = _j(RESULTS / "metapath_comparison.json")
    corr = mc["correction"]
    rows = [[k, n(f["n"]), f"{f['alpha']:.7f}",
             ", ".join(f["survivors"]) or "none"]
            for k, f in corr["families"].items()]
    v["correction_table"] = table(
        ["family", "tests", "corrected alpha", "survivors"], rows, "lrll")
    v["correction_note"] = corr["note"]
    v["mc_n_viruses"] = str(len(mc["per_virus"]))
    v["mc_channels"] = ", ".join(list(next(iter(mc["per_virus"].values()))))

    # ------------------------------------------------- pathway landscape
    pl = _j(RESULTS / "pathway_landscape.json")
    v["pl_permutations"] = n(pl["permutations"])
    v["pl_n_tests"] = str(pl["n_tests"])
    v["pl_seed_scheme"] = pl["seed_scheme"]
    v["pl_degree_basis"] = pl["degree_basis"]
    v["pl_drugged_proteome"] = n(pl["operands"]["drugged_proteome"]["n"])
    rows = []
    for curie, cell in sorted(pl["taxa"].items(),
                              key=lambda kv: names.get(kv[0], kv[0])):
        rows.append([names.get(curie, curie), n(cell["n_dependency"]),
                     n(cell["n_observed"]), f"{cell['expected']:.2f}",
                     f"{cell['p_naive']:.4f}", f"{cell['p_matched']:.4f}",
                     f"{cell['p_matched_less']:.4f}"])
    v["pl_table"] = table(
        ["virus", "dependency factors", "observed overlap", "expected",
         "naive p", "degree-matched p", "p (depletion)"], rows, "lrrrrrr")
    v["pl_uncovered"] = ", ".join(sorted(names.get(x, x)
                                         for x in pl["uncovered"]))
    # WHY A CELL IS UNCOVERED, FROM THE GRAPH RATHER THAN FROM A GUESS.
    # The first draft attributed all of the landscape's uncovered cells
    # to the absent `curated` layer. Three of them are nothing to do with
    # it, and the graph answers the question directly: a virus with no
    # HOST_FACTOR_FOR edge has no dependency-factor set to overlap.
    import collections
    from kgav.labels import ACTIVITY_PREDICATE, classify
    hf = collections.Counter()
    # ONE PASS, TWO CENSUSES. The schema and the labels module both quote
    # the share of the label predicate that is a measurement of
    # INACTIVITY, and both quote it from a comment. It is a property of
    # the release, so it is counted here.
    verdicts = collections.Counter()
    bounds = collections.Counter()
    for line in (RELEASE / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        if '"HOST_FACTOR_FOR"' in line:
            e = json.loads(line)
            if e["predicate"] == "HOST_FACTOR_FOR":
                hf[e["object"]] += 1
                continue
        if ACTIVITY_PREDICATE not in line:
            continue
        e = json.loads(line)
        if e["predicate"] != ACTIVITY_PREDICATE:
            continue
        verdicts[classify(e) or "undecidable"] += 1
        q = e.get("qualifiers") or {}
        if (q.get("relation") or "=") in (">", ">="):
            for f in ("ec50_nm", "ic50_nm", "ki_nm", "kd_nm"):
                if q.get(f) is not None:
                    bounds[float(q[f])] += 1
                    break
    total = sum(verdicts.values())
    v["n_label_edges"] = n(total)
    v["label_verdict_table"] = table(
        ["verdict", "edges", "share"],
        [[k, n(x), pct(x / total)]
         for k, x in verdicts.most_common()], "lrr")
    v["label_inactive_share"] = pct(verdicts["inactive"] / total)
    n_cens = sum(bounds.values())
    at_or_above = sum(x for b, x in bounds.items()
                      if b >= DEFAULT_INACTIVE_ABOVE_NM)
    v["n_censored_rows"] = n(n_cens)
    v["censored_at_or_above_share"] = pct(at_or_above / n_cens)
    v["censored_below_share"] = pct((n_cens - at_or_above) / n_cens)
    v["n_distinct_bounds"] = n(len(bounds))
    top = bounds.most_common(4)
    v["censored_bound_modes"] = ", ".join(
        f"{b / 1000:g} µM ({n(x)} rows, {pct(x / n_cens)})" for b, x in top)
    v["host_factor_table"] = table(
        ["virus", "HOST_FACTOR_FOR edges"],
        [[names.get(k, k), n(x)] for k, x in hf.most_common()], "lr")
    without = sorted(names[k] for k in names if k not in hf)
    v["n_taxa_without_host_factors"] = str(len(without))
    v["taxa_without_host_factors"] = ", ".join(without)
    v["pl_n_uncovered"] = str(len(pl["uncovered"]))

    # ----------------------------------------------------- label overlap
    lo = _j(RESULTS / "label_overlap.json")
    v["lo_both"] = n(lo["both"])
    v["lo_share_a_pub"] = n(lo["share_a_pub"])
    v["lo_share_pct"] = pct(lo["share_a_pub"] / lo["both_sides_have_pubs"])
    v["lo_unwithholdable"] = n(lo["unwithholdable"])
    v["lo_viral_proteins"] = n(lo["viral_proteins_reached"])

    # ------------------------------------------- degree weighting arms
    dw = _j(RESULTS / "degree_weighting_sensitivity.json")
    v["dw_cells"] = n(dw["cells"])
    v["dw_moved"] = str(dw["cells_moved"])
    v["dw_crossing"] = str(dw["cells_crossing_alpha"])
    v["dw_claim_cell"] = (f"{names.get(dw['claim_cell']['virus'], '?')} "
                          f"{dw['claim_cell']['channel']}")

    # ---------------------------------------------- repurposing read-out
    rc = _j(RESULTS / "repurposing_candidates.json")
    cf = rc["chemistry_filter"]
    v["rc_top"] = str(rc["top"])
    v["rc_n_molecules"] = n(cf["n_molecules"])
    v["rc_n_annotated"] = n(cf["n_annotated"])
    v["rc_n_by_inchikey"] = n(cf["n_annotated_by_inchikey"])
    v["rc_n_unchecked"] = n(cf["n_unchecked"])
    v["rc_n_diagnostic"] = n(cf["n_diagnostic_identified"])
    v["rc_n_diagnostic_keyed"] = n(cf["n_diagnostic_identified_by_inchikey"])
    v["rc_n_cytotoxic"] = n(cf["n_cytotoxic"])
    v["rc_n_withdrawn"] = n(cf["n_withdrawn"])
    v["rc_atc"] = ", ".join(f"`{k}` ({x})" for k, x
                            in sorted(cf["atc_non_therapeutic"].items()))
    rows, shown = [], 0
    for _curie, cell in sorted(rc["viruses"].items(),
                              key=lambda kv: kv[1]["name"]):
        fl = cell["chemistry_filter"]["resemblance_floor"]
        shown += len(cell["suggestions"])
        rows.append([cell["name"], cell["scorer_verdict"]["tier"],
                     n(cell["n_known_actives"]), n(cell["n_candidates"]),
                     auc(fl["floor"]), n(fl["n_inactives_scored"]),
                     str(fl["n_rows_displaced"]), str(len(cell["suggestions"]))])
    v["rc_table"] = table(
        ["virus", "verdict", "known actives", "candidates",
         "resemblance floor", "failures scored", "rows the floor cost",
         "rows shown"], rows, "llrrlrrr")
    v["rc_rows_shown"] = n(shown)

    # ------------------------------------------------------- the capability
    sys.path.insert(0, str(ROOT / "scripts"))
    from repurposing_candidates import tier as _tier       # noqa: E402
    from kgav.calibration import capabilities              # noqa: E402
    caps = capabilities(RESULTS, names=dict(names), cross_tier=_tier)
    rows = []
    for _curie, cap in sorted(caps.items(), key=lambda kv: kv[1].name):
        rows.append([cap.name, cap.cross_sectional, auc(cap.cross_auc),
                     cap.prospective, auc(cap.prospective_auc),
                     "**yes**" if cap.shortlist_licensed else "no",
                     "**yes**" if cap.contradicted else "no"])
    v["capability_table"] = table(
        ["virus", "cross-sectional", "AUC", "prospective", "AUC",
         "shortlist licensed", "contradicted"], rows, "lllllcc")
    v["n_licensed"] = str(sum(1 for c in caps.values() if c.shortlist_licensed))
    v["licensed_viruses"] = ", ".join(sorted(
        c.name for c in caps.values() if c.shortlist_licensed))
    v["contradicted_viruses"] = ", ".join(sorted(
        c.name for c in caps.values() if c.contradicted)) or "none"

    # ------------------------------------------------------------- agent
    v["agent_intents"] = ", ".join(f"`{k.value}`" for k in guards.QueryKind)
    v["agent_caveats_table"] = table(
        ["key", "label a synthesised answer must carry"],
        [[f"`{k}`", f"**{x}**"] for k, x in guards.CAVEAT_LABELS.items()])
    v["n_caveats"] = str(len(guards.CAVEAT_LABELS))
    v["agent_routed_kinds"] = ", ".join(
        f"`{k.value}`" for k in sorted(guards.ROUTED_KINDS, key=lambda x: x.value))
    v["agent_max_neighbours"] = str(agent_tools.MAX_NEIGHBOURS)
    v["agent_max_paths"] = str(agent_tools.MAX_PATHS)
    v["agent_max_subgraph"] = str(agent_tools.MAX_SUBGRAPH_EDGES)
    v["agent_thin_route"] = str(guards.THIN_ROUTE)
    v["agent_si_threshold"] = str(guards.SI_THRESHOLD)
    v["agent_biology_sources"] = table(
        ["claim", "reference"],
        [[f"`{k}`", s[1]] for k, s in guards.BIOLOGY_SOURCES.items()])
    v["n_anaphora"] = str(len(guards.ANAPHORA))
    v["agent_tools_table"] = table(
        ["tool", "what it answers"],
        [[f"`{name}`", (getattr(agent_tools.GraphTools, name).__doc__ or ""
                        ).strip().splitlines()[0]]
         for name in sorted(
             x for x in dir(agent_tools.GraphTools)
             if not x.startswith("_")
             and callable(getattr(agent_tools.GraphTools, x)))])

    # ------------------------------------------------------------- tests
    counts = {}
    for p in sorted((ROOT / "tests").glob("test_*.py")):
        src = p.read_text()
        counts[p.name] = src.count("\ndef test_") + src.count("\n    def test_")
    v["n_test_files"] = str(len(counts))
    v["n_test_functions"] = n(sum(counts.values()))
    v["n_slow_tests"] = n(sum(
        p.read_text().count("pytest.mark.slow")
        for p in (ROOT / "tests").glob("test_*.py")))
    v["biggest_test_files"] = table(
        ["test file", "test functions", "what it holds"],
        [[f"`tests/{k}`", str(x), _first_docline(ROOT / "tests" / k)]
         for k, x in sorted(counts.items(), key=lambda kv: -kv[1])[:20]], "lrl")

    # ----------------------------------------------------------- scripts
    rows = []
    for p in sorted((ROOT / "scripts").glob("*.py")):
        rows.append([f"`{p.name}`", str(len(p.read_text().splitlines())),
                     _first_docline(p)])
    v["script_table"] = table(["script", "lines", "what it does"], rows, "lrl")
    v["n_scripts"] = str(len(rows))
    v["n_src_modules"] = str(len(list((ROOT / "src" / "kgav").rglob("*.py"))))
    v["n_src_lines"] = n(sum(len(p.read_text().splitlines())
                             for p in (ROOT / "src" / "kgav").rglob("*.py")))
    v["n_script_lines"] = n(sum(len(p.read_text().splitlines())
                                for p in (ROOT / "scripts").glob("*.py")))
    v["n_test_lines"] = n(sum(len(p.read_text().splitlines())
                              for p in (ROOT / "tests").glob("*.py")))

    # ------------------------------------------------------------ make
    v["make_table"] = table(
        ["target", "the command it runs"],
        [[f"`make {t}`", d] for t, d in _make_targets()], "ll")
    v["n_make_targets"] = str(len(_make_targets()))

    # --------------------------------------------------- results census
    rows = []
    for p in sorted(RESULTS.glob("*.json")):
        d = _j(p)
        prov = d.get("provenance") or {}
        rel = prov.get("release") or {}
        rows.append([
            f"`{p.name}`", n(p.stat().st_size),
            "REFUSED" if "refused_by" in d else "result",
            (rel.get("edges_sha") or "—")[:8],
            (prov.get("generated_utc") or "—")[:10]])
    v["results_table"] = table(
        ["file", "bytes", "kind", "edges sha", "generated"], rows, "lrlll")
    v["n_results_files"] = str(len(rows))
    v["n_refused_results"] = str(sum(
        1 for p in RESULTS.glob("*.json") if "refused_by" in _j(p)))
    v["refusal_reasons"] = "\n".join(
        f"- `{p.name}`: " + "; ".join(_j(p)["refused_by"])
        for p in sorted(RESULTS.glob("*.json")) if "refused_by" in _j(p))

    return v


def _release_sha(which: str) -> str:
    from kgav.provenance import release_fingerprint
    return release_fingerprint(RELEASE)[f"{which}_sha"]


def _sources_text(src) -> str:
    """What a layer declared about the release of the source it read.

    A SOURCE THAT DECLARED NO RELEASE STRING SAYS SO. The first version
    fell back to the first 40 characters of the JSON, which rendered the
    indications layer's table names as a truncated object -- a cell that
    looks like a version and is not one, in the column whose whole job is
    to say which version was read.
    """
    if not src:
        return "—"
    out = []
    for k, x in sorted(src.items()):
        if isinstance(x, dict):
            ver = x.get("release") or x.get("published")
            if ver:
                out.append(f"{k} `{ver}`")
            else:
                out.append(f"{k}: no release string "
                           f"(declares {', '.join(sorted(x))})")
        else:
            out.append(f"{k} `{x}`")
    return ", ".join(out)


def _first_docline(path: Path) -> str:
    try:
        doc = ast.get_docstring(ast.parse(path.read_text())) or ""
    except SyntaxError:
        return "—"
    first = (doc.strip().splitlines() or [""])[0]
    return first.replace("|", "\\|")


def _make_targets() -> list[tuple[str, str]]:
    """Every target in the Makefile, with the COMMAND it runs.

    NOT WITH THE COMMENT ABOVE IT, which is what the first version did
    and which attributed the wrong description to four targets. A
    Makefile comment block is not delimited per target: the block above
    `lint-fast` discusses `make test`'s worker distribution first and
    `lint-fast`'s own rule-subset last, so taking the first line of the
    block put a sentence about one target beside another, and left the
    target it belonged to with no description at all.

    A recipe cannot be mis-attributed. It is also the thing a reader
    needs: what the target actually runs.
    """
    out: list[tuple[str, str]] = []
    lines = (ROOT / "Makefile").read_text().splitlines()
    for i, line in enumerate(lines):
        m = re.match(r"^([a-zA-Z][a-zA-Z0-9_-]*):", line)
        if not m or m.group(1) == ".PHONY":
            continue
        recipe: list[str] = []
        for nxt in lines[i + 1:]:
            if not nxt.startswith(("\t", " " * 4)):
                break
            cmd = nxt.strip()
            if cmd.startswith("@test -n") or not cmd:
                continue
            recipe.append(cmd.lstrip("@").replace("|", "\\|"))
        joined = "<br>".join(f"`{c}`" for c in recipe) or "—"
        out.append((m.group(1), joined))
    return out


# ------------------------------------------------------------- rendering
def render(values: dict[str, str]) -> dict[str, str]:
    """{filename: rendered text} for every template.

    An unknown placeholder RAISES. A figure that renders as an empty
    string reads like prose that was always meant to read that way, and
    this project has shipped that defect under four other names.
    """
    out, used = {}, set()
    for tpl in sorted(TEMPLATES.glob("*.md.in")):
        text = tpl.read_text()
        missing = sorted({m.group(1) for m in PLACEHOLDER.finditer(text)}
                         - set(values))
        if missing:
            raise SystemExit(
                f"{tpl.name} uses {missing}, which write_doc2.derive() does "
                f"not produce")
        used |= {m.group(1) for m in PLACEHOLDER.finditer(text)}
        rendered = PLACEHOLDER.sub(lambda m: values[m.group(1)], text)
        out[tpl.name[:-len(".in")]] = rendered.rstrip("\n") + "\n"
    unused = sorted(set(values) - used)
    if unused:
        print(f"note: {len(unused)} derived value(s) no template uses: "
              + ", ".join(unused), file=sys.stderr)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if any generated file differs, write nothing")
    args = ap.parse_args(argv)

    if not TEMPLATES.exists():
        print(f"no templates in {TEMPLATES}", file=sys.stderr)
        return 2
    files = render(derive())

    if args.check:
        bad = []
        for name, text in sorted(files.items()):
            path = OUT / name
            if not path.exists():
                bad.append(f"{name}: missing")
                continue
            have = path.read_text()
            if have != text:
                bad.append(f"{name}: {_first_diff(have, text)}")
        if bad:
            print("doc2 is out of date:")
            for b in bad:
                print(f"  {b}")
            print("\nrun `make doc2` to regenerate it.")
            return 1
        print(f"doc2 is current: {len(files)} file(s) match the code and the "
              f"artifacts.")
        return 0

    OUT.mkdir(parents=True, exist_ok=True)
    for name, text in sorted(files.items()):
        (OUT / name).write_text(text)
        print(f"wrote doc2/{name}  ({len(text.splitlines()):,} lines)")
    return 0


def _first_diff(have: str, want: str) -> str:
    a, b = have.splitlines(), want.splitlines()
    for i, (x, y) in enumerate(zip(a, b), 1):
        if x != y:
            return f"line {i} is {x[:60]!r}, should be {y[:60]!r}"
    return f"{len(a)} lines on disk, {len(b)} generated"


if __name__ == "__main__":
    raise SystemExit(main())
