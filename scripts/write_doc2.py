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

WHAT THIS DOES NOT DO. It does not reproduce doc2/registrations/prereg_*.md. Those are
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
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav import virus_register as vr                      # noqa: E402
from kgav.agent import guards                              # noqa: E402
from kgav.agent import llm                              # noqa: E402
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


def _analogy() -> dict:
    """The near-neighbour exclusion ladder, by running the diagnostic.

    THROUGH THE SCRIPT, one implementation, and it costs about ten
    seconds of `make doc2`. That is the difference between "6 of 10
    cells clear" being a result and being a summary that hides two
    different kinds of result inside one count: three of the six hold
    with forty of their nearest known actives withheld, and two stop
    clearing before ten.

    DELIBERATELY NOT IN data/results. It is a post-hoc diagnostic on a
    registered run, so it is labelled as one rather than stamped as a
    registered artifact -- the same decision `make hit-rate` records.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "analogy.json"
        r = subprocess.run(
            [sys.executable,
             str(ROOT / "scripts" / "diagnose_prospective_analogy.py"),
             "--out", str(out)],
            capture_output=True, text=True, cwd=ROOT, timeout=1800)
        if r.returncode != 0:
            raise SystemExit("diagnose_prospective_analogy.py failed:\n"
                             + (r.stdout or "") + (r.stderr or ""))
        return json.loads(out.read_text())


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
    # DECLARED AND EMPTY. A schema that declares more than the graph
    # holds is not a defect, but a reader has to be told which questions
    # therefore cannot be asked -- and naming them is the only way the
    # list stays true when one of them acquires edges, as TREATS and
    # IN_TRIAL_FOR did on 2026-11-17 while the scope statement went on
    # saying they were empty.
    unbuilt_c = sorted(set(schema.node_classes) - set(man["nodes_by_class"]))
    unbuilt_p = sorted({e.predicate for e in schema.edge_classes}
                       - set(man["edges_by_predicate"]))
    v["unbuilt_classes"] = ", ".join(f"`{x}`" for x in unbuilt_c) or "none"
    v["n_unbuilt_classes"] = str(len(unbuilt_c))
    v["unbuilt_predicates"] = ", ".join(f"`{x}`" for x in unbuilt_p) or "none"
    v["n_unbuilt_predicates"] = str(len(unbuilt_p))
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
    # WITH UNITS, because a bare "10,000" in a sentence about potency is
    # a number a reader has to guess the units of -- and in µM too,
    # because that is the unit the literature states the cut in.
    _ab = (f"{int(DEFAULT_ACTIVE_BELOW_NM):,} nM "
           f"({DEFAULT_ACTIVE_BELOW_NM / 1000:g} µM)")
    _ia = (f"{int(DEFAULT_INACTIVE_ABOVE_NM):,} nM "
           f"({DEFAULT_INACTIVE_ABOVE_NM / 1000:g} µM)")
    # AND WHETHER THEY ARE ONE CUT OR TWO. They are configurable
    # separately and are equal by default, so a document describing a
    # band between them would be describing a configuration nobody runs
    # -- which an earlier draft of 06-methodology.md did.
    v["threshold_shape"] = (
        f"a single cut at {_ab}: a measured potency at or "
        f"below it is **active** and anything above it is "
        f"**measured-inactive**. The two bounds are separate parameters "
        f"and are equal by default, so there is no band between them "
        f"on this release"
        if DEFAULT_ACTIVE_BELOW_NM == DEFAULT_INACTIVE_ABOVE_NM else
        f"two cuts — **active** at or below {_ab}, "
        f"**measured-inactive** at or above {_ia} — with "
        f"everything between them assigned to neither class")
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

    # EVERY REGISTERED CELL BESIDE ITS OWN DEGREE COMPARATOR, in four
    # groups. The headline claim is not "the scorer clears" but "the
    # scorer clears WHERE DEGREE DOES NOT", and those are different
    # counts: on this release 17 cells are readable, the scorer clears
    # in more of them than that claim uses, and the group where degree
    # clears too is not evidence about chemistry. A table that printed
    # only the first group would be the failure this project is built
    # against -- and so would one that is almost complete, which is why
    # the four group sizes are derived and summed below.
    GROUPS = ("clears, degree does not", "clears, degree clears too",
              "does not clear")
    cells_sd, group_of, degree_above = [], {}, 0
    for sect, pop in (("primary_chembl", "ChEMBL"),
                      ("secondary_ncats", "plate CPE"),
                      ("secondary_ncats_entry", "plate entry")):
        for curie, cell in sm.get(sect, {}).items():
            if not cell.get("readable"):
                continue
            s_ = cell.get("structure_max") or {}
            g = cell.get("degree") or {}
            if s_.get("auc") is None:
                continue
            s_clears = bool(s_.get("clears_corrected_alpha")) \
                and s_.get("direction") != "below"
            g_clears = bool(g.get("clears_corrected_alpha"))
            grp = (GROUPS[0] if s_clears and not g_clears
                   else GROUPS[1] if s_clears else GROUPS[2])
            if not s_clears and g_clears and (g.get("auc") or 0) > s_["auc"]:
                degree_above += 1
            group_of.setdefault(grp, []).append(1)
            cells_sd.append([
                names.get(curie, curie), pop, f"{s_['auc']:.4f}",
                "**yes**" if s_clears else "no",
                "—" if g.get("auc") is None else f"{g['auc']:.4f}",
                "yes" if g_clears else "no", grp])
    cells_sd.sort(key=lambda r: (GROUPS.index(r[6]), -float(r[2])))
    v["scorer_vs_degree_table"] = table(
        ["virus", "population", "structure AUC", "clears",
         "degree AUC", "clears", "group"], cells_sd, "lllcllc")
    v["n_scorer_cells"] = str(len(cells_sd))
    for i, grp in enumerate(GROUPS):
        v[f"n_cells_group{i + 1}"] = str(len(group_of.get(grp, ())))
    v["n_cells_degree_above"] = str(degree_above)

    # DO THE TWO CONTROLS AGREE, AND WHERE EXACTLY. The claim's main
    # defence is that the clearing set is not an artifact of one
    # exclusion rule -- and stated over the whole 17-cell family that is
    # FALSE on this release, which tests/test_doc2_claims.py found the
    # day the sentence was written. They agree on every ChEMBL cell and
    # on the cytopathic plate cell; they disagree on one pseudotyped-
    # entry cell, where degree clears and sits ABOVE the structure
    # statistic under both readings, so that cell is not evidence about
    # chemistry either way and is outside the claim in both.
    sym_doc = _j(RESULTS / "structure_scorer_symmetric.json")

    def _clearing(doc, field):
        out = set()
        for sect, pop in (("primary_chembl", "ChEMBL"),
                          ("secondary_ncats", "plate CPE"),
                          ("secondary_ncats_entry", "plate entry")):
            for curie, cell in (doc.get(sect) or {}).items():
                if not cell.get("readable"):
                    continue
                r = cell.get(field) or {}
                if (r.get("clears_corrected_alpha")
                        and r.get("direction") != "below"):
                    out.add((names.get(curie, curie), pop))
        return out

    sm_clear = _clearing(sm, "structure_max")
    sy_clear = _clearing(sym_doc, "structure_max_k10")
    disagree = sorted(sm_clear ^ sy_clear)
    v["n_controls_disagree"] = str(len(disagree))
    v["controls_disagree"] = ", ".join(
        f"{nm} ({pop})" for nm, pop in disagree) or "none"
    if disagree:
        v["controls_agreement"] = (
            f"The two registered near-neighbour controls — symmetric "
            f"k = {{sym_k}} and strength-matched — agree on every cell the "
            f"claim rests on, and disagree on "
            f"{len(disagree)}: {v['controls_disagree']}. In that cell the "
            f"degree baseline clears and sits above the structure "
            f"statistic under both readings, so it is not evidence about "
            f"chemistry either way and is outside the claim in both. The "
            f"membership of the {v['n_cells_group1']} cells that clear "
            f"above a non-clearing baseline is **identical under both "
            f"controls**.")
    else:
        v["controls_agreement"] = (
            f"The membership of the clearing set is **identical under both "
            f"registered controls** — symmetric k = {{sym_k}} and "
            f"strength-matched — in all {v['n_scorer_cells']} cells, which "
            f"is the result's main defence: it is not an artifact of one "
            f"exclusion rule.")
    v["controls_agreement"] = v["controls_agreement"].replace(
        "{sym_k}", v["sym_k"])

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

    # THE METAPATH CHANNELS UNDER THE SAME SPLIT, which is a different
    # study from the structure scorer above and must not be read as the
    # same one. This table had NO test until 2026-11-08 and every cell
    # in it was stale -- figures from a release two rebuilds earlier,
    # quoted at the top of the document as the summary of the whole
    # temporal protocol.
    ti = _j(RESULTS / "temporal_2021_include.json")
    rows = []
    for scorer, m in ti["results"]["SARS-CoV-2"].items():
        if not isinstance(m, dict) or "auc" not in m:
            continue
        marks = []
        if not m.get("hits_determined"):
            marks.append("~ties")
        if m.get("trustworthy") is False:
            marks.append("(!)")
        rows.append([
            scorer, n(m["pool"]), n(m["pos"]), f"{m['lift']:.2f}",
            "—" if m.get("p_value") is None else f"{m['p_value']:.4f}",
            f"{m['auc']:.3f} [{m['auc_ci_lo']:.3f}, {m['auc_ci_hi']:.3f}]",
            " ".join(marks)])
    v["temporal_metapath_table"] = table(
        ["scorer", "pool", "pos", "lift", "p", "AUC [95% CI]", ""],
        rows, "lrrrrll")
    # AND THE DECOMPOSITION OF IT, which is where this protocol stops
    # being a measurement of ranking at all: a cell whose AUC equals its
    # own imputation floor to four decimals is reporting its coverage.
    # THE FLOOR IS COMPUTED, NOT READ. temporal_split.py records the
    # coverages (auc_coverage_pos / auc_coverage_neg) and not the floor,
    # which is why the decomposition cost nothing new to measure -- and
    # why a first version of this, looking for an `auc_floor` field,
    # found none and would have published "no cell equals its own
    # floor": true of the field, false of the table.
    from kgav.labels import auc_ceiling as _ceiling     # noqa
    from kgav.labels import auc_floor as _floor         # noqa
    from kgav.labels import auc_skill as _skill         # noqa
    on_floor = []
    for k, m in ti["results"]["SARS-CoV-2"].items():
        if not isinstance(m, dict) or m.get("auc") is None:
            continue
        np_, nn_ = m.get("auc_n_pos"), m.get("auc_n_neg")
        cp, cn = m.get("auc_coverage_pos"), m.get("auc_coverage_neg")
        if None in (np_, nn_, cp, cn):
            continue
        fl = _floor(round(cp * np_), round(cn * nn_), np_, nn_)
        if round(m["auc"], 4) == round(fl, 4):
            on_floor.append(k)
    v["temporal_n_on_floor"] = str(len(on_floor))
    v["temporal_on_floor"] = ", ".join(sorted(on_floor)) or "none"
    clearing = [k for k, m in ti["results"]["SARS-CoV-2"].items()
                if isinstance(m, dict)
                and (m.get("auc_ci_lo") or 0) > CHANCE]
    v["temporal_metapath_clearing"] = ", ".join(sorted(clearing)) or "none"
    v["temporal_n_metapath_clearing"] = str(len(clearing))
    # ACROSS EVERY TEMPORAL CELL, not just SARS-CoV-2 under `include`.
    # The sharpest statement this protocol supports is a census, and a
    # census of one virus under one policy is not it.
    n_cells = n_floor = 0
    worst = span = 0.0
    for f in sorted(RESULTS.glob("temporal_2021_*.json")):
        for _virus, per in (_j(f).get("results") or {}).items():
            for _k, m in per.items():
                if not isinstance(m, dict) or m.get("auc") is None:
                    continue
                np_, nn_ = m.get("auc_n_pos"), m.get("auc_n_neg")
                cp, cn = m.get("auc_coverage_pos"), m.get("auc_coverage_neg")
                if None in (np_, nn_, cp, cn):
                    continue
                n_cells += 1
                fl = _floor(round(cp * np_), round(cn * nn_), np_, nn_)
                worst = max(worst, abs(m["auc"] - fl))
                span = max(span, _ceiling(round(cp * np_), round(cn * nn_),
                                          np_, nn_) - fl)
                if round(m["auc"], 4) == round(fl, 4):
                    n_floor += 1
    # THE BAR THE INTERVAL HAS TO BEAT AT THIS SAMPLE SIZE. The AUC is
    # computed over the WHOLE labelled test set, not over what a scorer
    # reaches -- the difference between the ensemble clearing the bar
    # and missing it.
    sizes = {(m["auc_n_pos"], m["auc_n_neg"])
             for m in ti["results"]["SARS-CoV-2"].values()
             if isinstance(m, dict) and "auc_n_pos" in m}
    if len(sizes) != 1:
        raise SystemExit(
            f"temporal scorers disagree on the AUC sample size: {sizes}")
    tn_pos, tn_neg = sizes.pop()
    from kgav.labels import auc_interval as _ival
    bar = CHANCE
    while _ival(bar, tn_pos, tn_neg)[0] <= CHANCE:
        bar += 0.0005
    # THE SAME DECOMPOSITION ON THE TEMPORAL CELLS, which on this
    # release leaves nothing behind. EMPTY IS A STATE THE DOCUMENT MUST
    # ASSERT, not one it may fall into by leaving a table out: five
    # cells were here on 2026-11-08 and all five were the similarity
    # leak, so a reader who finds no table has to be told which it is.
    dec_rows = []
    for pol in ("include", "computed_only", "exclude"):
        f = RESULTS / f"temporal_2021_{pol}.json"
        if not f.exists():
            continue
        for virus, chans in (_j(f).get("results") or {}).items():
            for scorer, m in chans.items():
                if not isinstance(m, dict) or m.get("auc_ci_lo") is None:
                    continue
                if m["auc_ci_lo"] <= CHANCE:
                    continue
                np_, nn_ = m["auc_n_pos"], m["auc_n_neg"]
                cp, cn = m["auc_coverage_pos"], m["auc_coverage_neg"]
                rp, rn = round(cp * np_), round(cn * nn_)
                fl = _floor(rp, rn, np_, nn_)
                ce = _ceiling(rp, rn, np_, nn_)
                sk = _skill(m["auc"], fl, ce)
                dec_rows.append([pol, virus, scorer, pct(cp), pct(cn),
                                 f"{fl:.4f}", f"{m['auc']:.4f}",
                                 f"{ce:.4f}",
                                 "n/a" if sk is None else pct(sk)])
    v["temporal_decomposition"] = (
        table(["policy", "virus", "scorer", "cov_pos", "cov_neg", "floor",
               "AUC", "ceiling", "skill"], dec_rows, "lllrrrrrr")
        if dec_rows else
        "**No temporal cell has an AUC interval excluding 0.5 on this "
        "release**, so this table is empty and the paragraph below is what "
        "the decomposition says instead. Five cells were here on 2026-11-08 "
        "and all five were the similarity leak.")

    v["temporal_auc_n"] = f"{tn_pos:,} and {tn_neg:,}"
    v["temporal_auc_bar"] = f"{bar:.3f}"
    v["temporal_n_auc_cells"] = str(n_cells)
    v["temporal_n_cells_on_floor"] = str(n_floor)
    v["temporal_max_floor_gap"] = f"{worst:.4f}"
    # THE ROOM A RANKING HAD TO SHOW ITSELF IN, which is the
    # other half of "the temporal AUC is a coverage statistic".
    v["temporal_max_span"] = f"{span:.4f}"

    # THE EXCLUSION LADDER: the six clearing cells are not one result.
    # Recomputing each cell's AUC with its k nearest pre-cutoff actives
    # withheld separates a chemotype signal from retrieval of a known
    # series, and the two are indistinguishable in the headline number.
    an = _analogy()
    ks = list(an["ks"])
    lad_rows = []
    for cell in sorted(an["cells"],
                       key=lambda c: -c["ladder"][str(ks[0])]["auc"]):
        row = [cell["name"]]
        last = None
        for k in ks:
            e = cell["ladder"][str(k)]
            row.append(f"{e['auc']:.4f}" + ("\\*" if e["clears"] else ""))
            if e["clears"]:
                last = k
        row.append(f"k≤{last}" if last is not None else "none")
        lad_rows.append(row)
    v["analogy_table"] = table(
        ["virus", *[f"k={k}" for k in ks], "clears up to"], lad_rows,
        "l" + "r" * len(ks) + "l")
    v["analogy_ks"] = ", ".join(str(k) for k in ks)
    v["analogy_k_max"] = str(ks[-1])
    survive = sorted(c["name"] for c in an["cells"]
                     if c["ladder"][str(ks[-1])]["clears"])
    v["analogy_survivors"] = ", ".join(survive) or "none"
    v["analogy_n_survivors"] = str(len(survive))
    fragile = sorted(c["name"] for c in an["cells"]
                     if not c["ladder"]["10"]["clears"])
    v["analogy_fragile"] = ", ".join(fragile) or "none"
    v["analogy_n_fragile"] = str(len(fragile))
    v["analogy_n_cells"] = str(len(an["cells"]))
    _W = {0: "None", 1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
          6: "Six", 7: "Seven", 8: "Eight", 9: "Nine", 10: "Ten"}
    v["analogy_heading"] = (
        f"{_W[len(survive)]} of those {_W[len(an['cells'])].lower()} are "
        f"analogy, and {_W[len(an['cells']) - len(survive)].lower()} are not")
    # THE RANK REVERSAL, derived rather than asserted: the cell with the
    # highest published AUC is the one that loses the most.
    by_primary = sorted(an["cells"],
                        key=lambda c: -c["ladder"][str(ks[0])]["auc"])
    by_robust = sorted(an["cells"],
                       key=lambda c: -c["ladder"][str(ks[-1])]["auc"])
    v["analogy_by_primary"] = " > ".join(c["name"] for c in by_primary)
    v["analogy_by_robust"] = " > ".join(c["name"] for c in by_robust)
    top = by_primary[0]
    v["analogy_top_name"] = top["name"]
    v["analogy_top_drop"] = (
        f"{top['ladder'][str(ks[0])]['auc']:.4f} to "
        f"{top['ladder'][str(ks[-1])]['auc']:.4f}")

    # EVERY REGISTRATION, FROM THE DIRECTORY. The predecessor of this
    # index was a hand-written table in doc2/08-results.md; a registration
    # added without an edit to it was simply absent, which is the one
    # thing an index of registrations must not be.
    regs = sorted((ROOT / "doc2" / "registrations").glob("prereg_*.md"))
    v["n_registrations"] = str(len(regs))
    v["registration_table"] = table(
        ["registration", "subject"],
        [[f"[`{p.name}`](registrations/{p.name})",
          _first_heading(p)] for p in regs])
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
        # THE TWO COMPONENTS OF THE MARGIN, so the share above it can be
        # checked by subtraction rather than believed.
        v[f"hn_{key}_coverage_component"] = f"{row['auc_floor'] - CHANCE:.4f}"
        v[f"hn_{key}_ranking_component"] = f"{row['auc'] - row['auc_floor']:.4f}"

    # THE DECOMPOSITION, WHICH IS WHAT MAKES TWO AUCs COMPARABLE. floor
    # is what a scorer's coverage asymmetry scores on its own, ceiling
    # is the best that coverage could reach, and skill is the fraction
    # of the span between them that the ranking used. All three are
    # recorded by the protocol; none is recomputed here, because a
    # document that re-derives a statistic can disagree with the
    # artifact it is describing and look right doing it.
    dec = []
    for name_ in [*metapaths, "COMBINED", "degree"]:
        row = cells.get(name_)
        if not isinstance(row, dict) or row.get("auc") is None:
            continue
        bold = "**" if name_ in ("COMBINED", "degree") else ""
        dec.append([f"{bold}{name_}{bold}",
                    pct(row["coverage_pos"]), pct(row["coverage_neg"]),
                    f"{bold}{row['auc_floor']:.4f}{bold}",
                    f"{row['auc']:.4f}",
                    f"{row['auc_ceiling']:.4f}",
                    "n/a" if row.get("auc_skill") is None
                    else f"{bold}{pct(row['auc_skill'])}{bold}"])
    v["auc_decomposition_table"] = table(
        ["scorer", "reaches + ", "reaches − ", "floor", "AUC", "ceiling",
         "skill"], dec, "lrrrrrr")

    # THE SEQUENCE: one control added per row. Every row writes its OWN
    # file, because all five once wrote one hard_negatives.json and only
    # the last survived -- which is how a five-row table in docs/ came
    # to describe four protocols the graph had stopped producing while
    # check-results certified the one file it could see.
    seq, seq_rows = [], [
        ("hard_negatives_row1.json", "`selective-only`, positives only"),
        ("hard_negatives_row2.json", "`all` — no selectivity filter"),
        ("hard_negatives_row3.json", "`verified-only`, **both** sides"),
        ("hard_negatives_row4.json", "+ `--publication-disjoint`"),
        ("hard_negatives.json", "+ `--label-source` (one screen)"),
    ]
    for i, (fname, what) in enumerate(seq_rows, 1):
        d = _j(RESULTS / fname)
        row = d["cross_sectional"]["SARS-CoV-2"]["M1"]
        reached = round(row["coverage_pos"] * row["n_pos"])
        cell = (f"{row['auc']:.4f} [{row['auc_ci_lo']:.4f}, "
                f"{row['auc_ci_hi']:.4f}]"
                if reached >= MIN_REACHED_POSITIVES
                else "**not evaluable**")
        seq.append([str(i), what, n(row["n_pos"]), n(row["n_neg"]),
                    pct(row["coverage_pos"]), pct(row["coverage_neg"]),
                    str(reached), cell])
    v["sequence_table"] = table(
        ["#", "protocol", "+", "−", "reaches + ", "reaches − ",
         "actives reached", "M1"], seq, "llrrrrrl")

    # ---------------------------------------------- metapath comparison
    mc = _j(RESULTS / "metapath_comparison.json")
    # THE UNCONTROLLED VIEW, WHICH IS WHERE THE DISCOVERIES LIVE. This
    # protocol applies none of the controls: it runs on unfiltered
    # build_labels output, so a p-value here is the probability of that
    # many hits under a RANDOM RANKING and cannot see a confound that
    # inflates the scorer and the labels together. Printing it is the
    # point -- it is the table a reader would have produced without the
    # rest of this repository, and the reason to distrust it is visible
    # in its own columns.
    from kgav.baselines import confound_warning as _confound
    CHANS = ("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
    rows, unreadable = [], []
    for virus, per in sorted(mc["per_virus"].items()):
        readable = [c for c in CHANS
                    if (per.get(c) or {}).get("trustworthy")
                    and (per.get(c) or {}).get("hits_determined")]
        best = max(readable, key=lambda c: per[c]["lift"], default=None)
        sig = [c for c in readable if per[c]["p_value"] < 0.05]
        lowest = min(sig, key=lambda c: per[c]["p_value"], default=None)
        deg = per.get("degree") or {}
        if best is None:
            unreadable.append(virus)
        rows.append([
            virus,
            best or "—",
            f"{per[best]['lift']:.2f}" if best else "—",
            f"{per[best]['p_value']:.4f}" if best else "—",
            f"{lowest} {per[lowest]['p_value']:.4f}" if lowest else "—",
            f"{deg.get('lift', 0):.2f}",
            # <.0001 RATHER THAN 0.0000. A p rendered as zero reads as a
            # certainty, and the two degree rows that produce it are the
            # two this table exists to warn about.
            ("<.0001" if deg.get("p_value", 1) < 1e-4
             else f"{deg.get('p_value', 1):.4f}"),
            # THE VERDICT, NOT THE PARAGRAPH. confound_warning returns
            # the whole explanation; the first clause is the verdict and
            # the rest belongs under the table, not in a cell.
            ((_confound(per) or "—").split(":")[0]).strip()])
    v["mc_table"] = table(
        ["virus", "best readable channel", "lift", "p", "significant",
         "degree", "p", "baseline verdict"], rows, "llrrlrrl")
    v["mc_unreadable"] = ", ".join(sorted(unreadable))
    # A ROW IS UNREADABLE FOR ONE OF TWO RECORDED REASONS, and both are
    # severe here. `hits_determined: False` means the published hit
    # count is one draw from a tie block straddling the cut, so the lift
    # computed from it is one draw too; `trustworthy: False` means the
    # expectation at the cut is a fraction of a hit, so a single hit
    # produces an enormous lift. Every eye-catching double-digit lift in
    # this project is of the second shape.
    ties, thin = [], []
    for virus, per in sorted(mc["per_virus"].items()):
        for c in CHANS:
            cell = per.get(c) or {}
            if cell.get("pos", 0) == 0:
                continue
            if not cell.get("hits_determined"):
                ties.append((virus, c, cell))
            elif not cell.get("trustworthy"):
                thin.append((virus, c, cell))
    ties.sort(key=lambda x: x[2].get("p_value", 1))
    thin.sort(key=lambda x: -x[2].get("lift", 0))
    v["mc_tie_table"] = table(
        ["virus", "channel", "nominal lift", "p", "hits range", "positives"],
        [[vv, c, f"{x['lift']:.2f}", f"{x['p_value']:.6f}",
          f"[{x['hits_lo']}, {x['hits_hi']}]", n(x["pos"])]
         for vv, c, x in ties[:6]], "llrrlr")
    v["mc_thin_table"] = table(
        ["virus", "channel", "lift", "hits expected at the cut"],
        [[vv, c, f"{x['lift']:.2f}", f"{x['expected']:.2f}"]
         for vv, c, x in thin[:6]], "llrr")
    v["mc_n_ties"] = str(len(ties))
    v["mc_n_thin"] = str(len(thin))
    no_reach = sorted(vv for vv, per in mc["per_virus"].items()
                      if (per.get("COMBINED") or {}).get("pos", 0) == 0)
    v["mc_no_reachable_positive"] = ", ".join(no_reach)
    empty_pool = sorted(vv for vv, per in mc["per_virus"].items()
                        if (per.get("COMBINED") or {}).get("pool", 0) == 0)
    v["mc_n_empty_pool"] = str(len(empty_pool))
    # "A and B" RATHER THAN "A, B": the two shapes of unreachability are
    # not alike and the paragraph has to distinguish them in prose, not
    # in a comma-separated list a reader skims past.
    v["mc_empty_pool"] = (" and ".join(empty_pool) if len(empty_pool) < 3
                          else ", ".join(empty_pool[:-1])
                          + " and " + empty_pool[-1]) or "none"
    # SPELLED OUT, because the prose does and a reader counting rows in
    # a fifteen-row table should not have to.
    WORD = {0: "No", 1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
            6: "Six", 7: "Seven", 8: "Eight", 9: "Nine", 10: "Ten",
            11: "Eleven", 12: "Twelve", 13: "Thirteen", 14: "Fourteen",
            15: "Fifteen", 16: "Sixteen"}
    nv = len(mc["per_virus"])
    v["mc_unreadable_sentence"] = (
        f"{WORD[len(unreadable)]} of {WORD[nv].lower()} viruses have no "
        f"readable channel at all")
    v["mc_no_reach_sentence"] = (
        f"{WORD[len(no_reach)]} of {WORD[nv].lower()} viruses have no "
        f"positive reachable by any metapath")
    readable = [(vv, c, m) for vv, per in mc["per_virus"].items()
                for c, m in per.items()
                if c not in ("degree", "COMBINED")
                and m.get("hits_determined") and m.get("trustworthy")]
    sig = [(vv, c, m) for vv, c, m in readable
           if m.get("p_value") is not None and m["p_value"] < 0.05]
    v["mc_readable_sentence"] = (
        f"{WORD[len(readable)].lower()} cells are both `trustworthy` and "
        f"`hits_determined`")
    v["mc_significant_sentence"] = (
        f"{WORD[len(sig)]} of them {'is' if len(sig) == 1 else 'are'} "
        f"nominally significant")
    v["mc_significant_cells"] = ", ".join(
        f"{vv} {c}, lift {m['lift']:.2f}, p={m['p_value']:.6f}"
        for vv, c, m in sorted(sig, key=lambda x: x[2]["p_value"]))
    v["mc_smallest_p"] = (f"{min(m['p_value'] for _, _, m in sig):.6f}"
                          if sig else "—")

    # EVERY CELL IN THE PROJECT WHOSE INTERVAL EXCLUDES CHANCE, with the
    # floor beside each. THE CENSUS THAT EXISTED IN NO TEST: the claim
    # block enumerated these in prose -- "two scorers do, both only
    # cross-sectionally and both on SARS-CoV-2" -- while five of the
    # seven were temporal and two were a different virus. An interval
    # excluding chance is not a result when the imputation floor sits on
    # top of it, so a census that prints the interval and not the floor
    # is the sentence this project kept having to retract, one column
    # short.
    CLEAR_RUNS = {
        "cross-sectional": ("hard_negatives.json", "cross_sectional"),
        "temporal include": ("temporal_2021_include.json", "results"),
        "temporal computed_only": ("temporal_2021_computed_only.json",
                                   "results"),
        "temporal exclude": ("temporal_2021_exclude.json", "results"),
    }
    census = []
    for label, (fname, key) in CLEAR_RUNS.items():
        f = RESULTS / fname
        if not f.exists():
            continue
        for virus, chans in (_j(f).get(key) or {}).items():
            if not isinstance(chans, dict):
                continue
            for scorer, m in chans.items():
                if not isinstance(m, dict):
                    continue
                lo, hi = m.get("auc_ci_lo"), m.get("auc_ci_hi")
                if lo is None or hi is None or (lo <= CHANCE <= hi):
                    continue
                if "coverage_pos" in m:
                    np_, nn_ = m["n_pos"], m["n_neg"]
                    cp_, cn_ = m["coverage_pos"], m["coverage_neg"]
                else:
                    np_, nn_ = m["auc_n_pos"], m["auc_n_neg"]
                    cp_, cn_ = m["auc_coverage_pos"], m["auc_coverage_neg"]
                fl = _floor(round(cp_ * np_), round(cn_ * nn_), np_, nn_)
                share = ((fl - CHANCE) / (m["auc"] - CHANCE)
                         if m["auc"] != CHANCE else 0.0)
                census.append([
                    label, virus, scorer,
                    f"{scorer} {m['auc']:.3f} [{lo:.3f}, {hi:.3f}]",
                    f"{fl:.4f}", pct(share, 0)])
    census.sort()
    # REBUILT SO THE RENDERED FORM THE TEST READS IS LITERALLY PRESENT:
    # `{scorer} {auc} [lo, hi]` is one string, and splitting it across
    # table cells would make the census unreadable to anything that
    # checks it.
    v["chance_census"] = (
        "\n".join(f"    {r[0]:<24} {r[1]:<12} {r[3]}  floor {r[4]}  "
                   f"{r[5]} of the margin is reach" for r in census)
        if census else "    (no cell on this release)")
    v["n_chance_census"] = str(len(census))
    v["chance_census_viruses"] = ", ".join(
        sorted({r[1] for r in census})) or "none"

    # THE BASELINE VERDICT CENSUS, over every (run, virus) pair in the
    # project. These counts went stale three times in six days while the
    # prose around them was rewritten by hand, and only one of the three
    # classes had a test -- in the direction that assumed it held
    # somewhere. A verdict that holds NOWHERE is as assertable as one
    # that names a holder, so the block below says which.
    VERDICT_RUNS = {"cross-sectional": ("metapath_comparison.json",
                                        "per_virus")}
    for policy in ("computed_only", "include", "exclude"):
        VERDICT_RUNS[f"temporal {policy}"] = (
            f"temporal_2021_{policy}.json", "results")
    vc = {"CONFOUNDED": [], "BASELINE AT CHANCE": [],
          "BASELINE NOT DETERMINED": []}
    for label, (fname, key) in VERDICT_RUNS.items():
        f = RESULTS / fname
        if not f.exists():
            continue
        for virus, row in (_j(f).get(key) or {}).items():
            w = _confound(row)
            if w:
                vc[w.split(":")[0]].append(f"{label} / {virus}")
    WORDS = {0: "No", 1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
             6: "Six", 7: "Seven", 8: "Eight", 9: "Nine", 10: "Ten",
             11: "Eleven", 12: "Twelve", 13: "Thirteen", 14: "Fourteen",
             15: "Fifteen", 16: "Sixteen", 17: "Seventeen",
             18: "Eighteen", 19: "Nineteen", 20: "Twenty"}
    MEANS = {
        "CONFOUNDED": "the degree baseline is itself significant, so a "
                      "channel beating random says nothing about mechanism",
        "BASELINE AT CHANCE": "degree is indistinguishable from random, "
                              "which is the only condition under which a "
                              "channel's p means what it appears to mean",
        "BASELINE NOT DETERMINED": "the degree row's own hit count is one "
                                   "draw from a tie block, so no comparison "
                                   "against it can be read",
    }
    lines = []
    for k in ("CONFOUNDED", "BASELINE AT CHANCE", "BASELINE NOT DETERMINED"):
        pairs = sorted(vc[k])
        count = (f"{WORDS.get(len(pairs), len(pairs))} (run, virus) "
                 f"{'pair' if len(pairs) == 1 else 'pairs'}")
        where = (", ".join(pairs) if pairs
                 else "NO (run, virus) pair on this release")
        lines.append(f"    {k}\n"
                     f"        {count} — {MEANS[k]}.\n"
                     f"        {where}.")
    lines.append("    (nothing)\n"
                 "        every other pair: the comparison is readable and "
                 "no verdict applies.")
    v["verdict_block"] = "\n".join(lines)
    # THE TWO PROSE FIGURES BESIDE THE TABLE, and the ratio only while
    # it exists: when degree's own Hits@k is an interval the artifact
    # records `lift_vs_degree: null`, and a document that recomputed one
    # from the point draw would be quoting a number nobody measured.
    m4 = mc["per_virus"]["SARS-CoV-2"]["M4"]
    deg = mc["per_virus"]["SARS-CoV-2"]["degree"]
    v["mc_degree_lift"] = f"**{deg['lift']:.2f}**"
    # THE AT-CHANCE CELLS, SPELLED OUT WITH BOTH SIDES' FIGURES. The
    # verdict is a statement about the baseline, so it is only readable
    # beside the channel it licenses -- and the warning text itself
    # quotes one, so the document has to agree with it.
    at_chance = []
    for virus, per in sorted(mc["per_virus"].items()):
        w = _confound(per) or ""
        if not w.startswith("BASELINE AT CHANCE"):
            continue
        d = per["degree"]
        carried = w.split("so ", 1)[1].split("'s", 1)[0]
        cm = per[carried]
        at_chance.append(
            f"On {virus}, degree lift {d['lift']:.2f} "
            f"(p={d['p_value']:.4f}) is indistinguishable from random, so "
            f"{carried}'s {cm['lift']:.2f} (p={cm['p_value']:.4f}) is not "
            f"explained by study volume")
    # SARS-CoV M1: A RANGE AND NO POINT. The hit count is one draw from
    # a tie block, and writing the release in a canonical order moved
    # the published lift from 2.08 to 0.83 without anything about the
    # graph changing -- because rank() sorts on score with no tiebreak
    # and Python's sort is stable, so which of the tied compounds fall
    # inside Hits@100 follows the order edges.jsonl was written in.
    und = mc["per_virus"]["SARS-CoV"]["M1"]
    v["mc_undetermined_cell"] = (
        f"SARS-CoV M1's Hits@100 is [{und['hits_lo']}, {und['hits_hi']}] of "
        f"{und['pos']} positives and its lift is "
        f"[{und['lift_lo']:.2f}, {und['lift_hi']:.2f}], which "
        f"{'contains' if und['lift_lo'] <= 1.0 <= und['lift_hi'] else 'excludes'}"
        f" 1.00 — no enrichment at all is inside the range"
        if und.get("hits_determined") is False else
        f"SARS-CoV M1 is determined on this release, at lift "
        f"{und['lift']:.2f}")
    v["mc_at_chance"] = (". ".join(at_chance) + "."
                         if at_chance else
                         "No (run, virus) pair on this release reads "
                         "BASELINE AT CHANCE.")
    v["mc_m4_vs_degree"] = (
        f"{m4['lift_vs_degree']:.2f}x a degree baseline that is itself at "
        f"chance here"
        if m4.get("lift_vs_degree") is not None else
        "no ratio against degree at all, because degree's own Hits@k on "
        "this cell is an interval rather than a count")
    corr = mc["correction"]
    #: What each family IS, in a reader's words, and how its size is
    #: arrived at. EDITORIAL and complete: a family the artifact adds
    #: without an entry here fails the generator rather than being
    #: printed under its own key, because "cells_computed" names the
    #: field and not the thing, and the difference between two of these
    #: families is a five-fold change in the threshold the one
    #: significant cell is judged against.
    FAMILY = {
        "readable_cells": ("readable cells only",
                           "cells that are both `trustworthy` and "
                           "`hits_determined`"),
        "cells_computed": ("every cell computed",
                           "every (virus, channel) the run produced"),
        "cells_planned": ("every cell planned",
                          f"({len(mc['per_virus'])} viruses x "
                          f"{len(schema.metapaths)} metapaths)"),
        "planned_plus_baselines": (
            "planned, plus COMBINED and degree",
            f"({len(mc['per_virus'])} viruses x "
            f"{len(schema.metapaths) + 2} scorers)"),
    }
    unnamed = sorted(set(corr["families"]) - set(FAMILY))
    if unnamed:
        raise SystemExit(
            f"metapath_comparison.json declares correction families "
            f"{unnamed} that write_doc2.FAMILY does not name")
    # THE VERDICT IS COMPUTED, and it is the sentence a reader takes
    # away: "clears" beside an alpha smaller than p is the whole point
    # of the table. The factor is computed too -- one row read "by a
    # factor of 10" against a threshold twelve times the p-value.
    best_p = min((m["p_value"] for per in mc["per_virus"].values()
                  for c, m in per.items()
                  if c not in ("degree", "COMBINED")
                  and m.get("hits_determined") and m.get("trustworthy")
                  and m.get("p_value") is not None), default=None)
    rows = []
    for k, f in corr["families"].items():
        a = f["alpha"]
        if best_p is None:
            verdict = "no readable cell is significant"
        elif best_p < a:
            verdict = f"clears, by a factor of {round(a / best_p)}"
        else:
            # NO FACTOR ON A ROW THAT DOES NOT CLEAR. Stating one here
            # means stating alpha/p, which is below 1 and rounds to
            # nothing readable; stating p/alpha instead puts a second
            # convention in one column.
            verdict = "does not clear"
        rows.append([FAMILY[k][0], n(f["n"]), f"{a:.5f}", verdict,
                     FAMILY[k][1], ", ".join(f["survivors"]) or "none"])
    # THE JUDGED CELL'S p IS IN THE HEADER. A verdict column with no
    # p beside it is an opinion: the reader has to be able to do the
    # comparison, and a test has to be able to find what was compared.
    v["correction_table"] = table(
        ["family", "tests", "corrected alpha",
         f"the smallest readable cell, at p={best_p:.6f}"
         if best_p is not None else "the smallest readable cell",
         "how the size is arrived at", "survivors"], rows, "lrllll")
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
    # BOTH CHANNELS, OR THE TABLE CONTRADICTS THE SENTENCE BESIDE IT.
    # This printed the `taxa` channel alone -- 8 rows -- under a
    # paragraph saying 44 tests. 44 is 22 rows tested in two directions,
    # and the 14 missing rows are the whole interactome channel: the
    # softer, larger, more interesting half, including the smallest
    # p-value in the analysis. A table that is almost complete is the
    # failure this project is built against.
    pl_alpha = 0.05 / pl["n_tests"]
    rows = []
    for channel in ("taxa", "interactome"):
        for curie, cell in sorted(pl[channel].items(),
                                  key=lambda kv: names.get(kv[0], kv[0])):
            hi, lo = cell["p_matched"], cell["p_matched_less"]
            rows.append([channel, names.get(curie, curie),
                         n(cell["n_dependency"]),
                         n(cell["n_observed"]), f"{cell['expected']:.2f}",
                         f"{cell['p_naive']:.4f}", f"{hi:.4f}", f"{lo:.4f}",
                         # THE VERDICT IS COMPUTED FROM THE FILE'S OWN
                         # ALPHA, not restated beside it. A table whose
                         # verdict column is typed can disagree with the
                         # threshold printed above it.
                         "ENRICHED" if hi < pl_alpha else
                         "DEPLETED" if lo < pl_alpha else "—"])
    v["pl_table"] = table(
        ["channel", "virus", "dependency factors", "observed overlap",
         "expected", "naive p", "degree-matched p", "p (depletion)",
         "verdict"],
        rows, "llrrrrrrl")
    v["pl_n_rows"] = str(len(rows))
    v["pl_alpha_exact"] = f"{pl_alpha:.5f}"
    pl_clearing = [r for r in rows if r[8] != "—"]
    v["pl_n_clearing"] = str(len(pl_clearing))
    # THE UNCORRECTED CONTRAST, FOR THE TWO ROWS IT IS SHARPEST ON. The
    # gap between the naive and the degree-matched p IS the finding of
    # this section, and quoting it for a row where the two agree would
    # hide that.
    # EVERY ROW THE UNCORRECTED TEST WOULD HAVE CALLED A DISCOVERY, in
    # both channels, rather than two chosen ones: a selection would be
    # the thing this section is about.
    naive_gap = sorted(
        ((names.get(c, c), ch, r["p_naive"], r["p_matched"])
         for ch in ("taxa", "interactome")
         for c, r in pl[ch].items()
         if "p_naive" in r and r["p_naive"] < 0.05 and r["p_matched"] >= 0.05),
        key=lambda t: t[2])
    v["pl_n_naive_discoveries"] = str(len(naive_gap))
    v["pl_naive_contrast"] = "; ".join(
        f"{nm}'s {ch} overlap is {pn:.4f} naive and {pm:.4f} matched"
        for nm, ch, pn, pm in naive_gap)
    v["pl_n_taxa_rows"] = str(len(pl["taxa"]))
    v["pl_n_interactome_rows"] = str(len(pl["interactome"]))
    v["pl_uncovered"] = ", ".join(sorted(names.get(x, x)
                                         for x in pl["uncovered"]))

    # WHAT THE NULL COULD HAVE DETECTED. A null result is only a finding
    # if the reader knows what would have survived: failing to detect a
    # 2x excess is not evidence against a 2x excess when the test needed
    # 7x. This table did not exist anywhere in the project until
    # 2026-11-02, and "nothing survived correction" was being reported
    # without it.
    from kgav import landscape as _L
    land = _L.load(RESULTS / "pathway_landscape.json", release=man["release"])
    dead = {(r.channel, r.virus) for r in land.powerless_depletion()}
    rows = []
    for r in sorted(land.family, key=lambda r: r.expected):
        k, ratio, _ = r.minimum_detectable
        rows.append([r.channel, names.get(r.virus, r.virus),
                     f"{r.expected:.2f}", str(r.n_observed),
                     f"{k} hits = **{ratio:.2f}x**",
                     "**no rejection region**"
                     if (r.channel, r.virus) in dead else "testable"])
    v["pl_power_table"] = table(
        ["channel", "virus", "expected", "observed", "min. detectable",
         "downward test"], rows, "llrrll")
    _PW = {0: "No", 1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
           6: "Six", 7: "Seven", 8: "Eight", 9: "Nine", 10: "Ten",
           20: "Twenty", 21: "Twenty-one", 22: "Twenty-two",
           23: "Twenty-three", 24: "Twenty-four"}
    v["pl_powerless_sentence"] = (
        f"{_PW.get(len(dead), len(dead))} of the "
        f"{_PW.get(len(land.family), len(land.family)).lower()} downward "
        f"tests have no rejection region at all")
    # AND HOW COARSE THE REST OF THE FAMILY IS. "No enrichment was
    # found" over a family where most rows could only have detected a
    # 2.4-fold excess is a different sentence from the same words over
    # a family that could see 1.2.
    coarse = sum(1 for r in land.family if r.minimum_detectable[1] >= 2.4)
    v["pl_coarse_sentence"] = (
        f"{_PW.get(coarse, coarse)} of the "
        f"{_PW.get(len(land.family), len(land.family)).lower()} need a 2.4x "
        f"excess or more before any observation could clear the bar")
    v["pl_n_powerless"] = str(len(dead))
    v["pl_n_family_rows"] = str(len(land.family))
    # THE ROWS THAT INVITE A READING AS DEPLETION AND MUST NOT GET ONE.
    # A ratio computed on an expectation of three is one Poisson draw.
    thin = sorted((r for r in land.family if r.expected < 10),
                  key=lambda r: r.expected)
    v["pl_thin_rows"] = ", ".join(
        f"{names.get(r.virus, r.virus)} ({r.channel}: {r.n_observed} observed, "
        f"{r.expected:.2f} expected)" for r in thin) or "none"
    v["pl_n_thin_rows"] = str(len(thin))

    # TEN DISEASE NODES FOR SIXTEEN VIRUSES, which was for ten days the
    # stated reason the clinical labels were empty, and was false when
    # it was stated. Derived both ways: the register's distinct MONDO
    # ids and the release's Disease nodes have to agree.
    mondos = {e.get("mondo") for e in reg["viruses"] if e.get("mondo")}
    no_mondo = [e["name"] for e in reg["viruses"] if not e.get("mondo")]
    n_disease = int(man["nodes_by_class"].get("Disease", 0))
    if len(mondos) != n_disease or no_mondo:
        raise SystemExit(
            f"{len(mondos)} distinct MONDO ids in the register, "
            f"{n_disease} Disease nodes in the release, {len(no_mondo)} "
            f"entries with none: 09-limits.md states one number for both")
    v["n_disease_nodes"] = str(n_disease)
    # UNCOVERED IN WHICH CHANNEL. `uncovered` is the taxa channel's, and
    # reading it as "could not be tested at all" is wrong for six of the
    # eight: they have an interactome row. The pair missing from BOTH is
    # the one the agent refuses to rank, arrived at independently.
    both = sorted(set(pl["uncovered"]) & set(pl["interactome_uncovered"]))
    v["pl_n_uncovered_interactome"] = str(len(pl["interactome_uncovered"]))
    v["pl_uncovered_interactome"] = ", ".join(
        sorted(names.get(x, x) for x in pl["interactome_uncovered"])) or "none"
    v["pl_n_uncovered_both"] = str(len(both))
    v["pl_uncovered_both"] = ", ".join(
        sorted(names.get(x, x) for x in both)) or "none"
    # WHY A CELL IS UNCOVERED, FROM THE GRAPH RATHER THAN FROM A GUESS.
    # The first draft attributed all of the landscape's uncovered cells
    # to the absent `curated` layer. Three of them are nothing to do with
    # it, and the graph answers the question directly: a virus with no
    # HOST_FACTOR_FOR edge has no dependency-factor set to overlap.
    import collections
    from kgav.labels import ACTIVITY_PREDICATE, classify
    hf = collections.Counter()
    # THE SCOPE CENSUS, in the same pass. `HOST_FACTOR_FOR` splits by its
    # `direction` qualifier into dependency and restriction factors, and
    # only the dependency side enters the overlap analysis -- so a virus
    # with restriction factors and no dependency factors is in the graph
    # and absent from that table, which is a distinction the earlier
    # scope statement made and doc2 had dropped.
    dep = collections.Counter()
    res = collections.Counter()
    av = collections.Counter()
    # WHAT AN EVIDENCE DATE MEANS DIFFERS BY SOURCE, and the temporal
    # protocol rests on it. `1970-01-01` is the placeholder every ingest
    # writes when it has no date; counting it by source is the only way
    # a reader can tell "undated" from "dated 1970".
    undated = collections.Counter()
    dated = 0
    n_targets = n_undirected = 0
    # ONE PASS, TWO CENSUSES. The schema and the labels module both quote
    # the share of the label predicate that is a measurement of
    # INACTIVITY, and both quote it from a comment. It is a property of
    # the release, so it is counted here.
    verdicts = collections.Counter()
    bounds = collections.Counter()
    for line in (RELEASE / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        dated += 1
        if '"first_asserted_date": "1970-01-01"' in line:
            undated[json.loads(line).get("primary_knowledge_source")
                    or "unrecorded"] += 1
        if '"TARGETS"' in line:
            e = json.loads(line)
            if e["predicate"] == "TARGETS":
                n_targets += 1
                # "unknown" IS NO DIRECTION. The qualifier is present on
                # every edge and carries the string `unknown` where the
                # source states none, so a test looking for a missing
                # key counts zero and a document typing the figure by
                # hand is the only thing that ever had it right.
                if (e.get("qualifiers") or {}).get("direction") in (
                        None, "", "unknown"):
                    n_undirected += 1
        if '"HOST_FACTOR_FOR"' in line:
            e = json.loads(line)
            if e["predicate"] == "HOST_FACTOR_FOR":
                hf[e["object"]] += 1
                d = (e.get("qualifiers") or {}).get("direction")
                (dep if d == "dependency" else res)[e["object"]] += 1
                continue
        if ACTIVITY_PREDICATE not in line:
            continue
        e = json.loads(line)
        if e["predicate"] != ACTIVITY_PREDICATE:
            continue
        if str(e.get("object", "")).startswith("NCBITaxon:"):
            av[e["object"]] += 1
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
    # THE ANALYSABLE SET, and what the untouched remainder is made of.
    # The host-directed channels can only reach a compound carrying a
    # human TARGETS edge; the rest of the library is out of their reach
    # by construction rather than by scoring badly. The phase breakdown
    # of the remainder is what says how much of it is even worth
    # reaching: a preclinical compound is not a repurposing candidate.
    phase = {}
    for line in (RELEASE / "nodes.jsonl").read_text().splitlines():
        if '"SmallMolecule"' not in line:
            continue
        nd = json.loads(line)
        if nd["class"] == "SmallMolecule":
            phase[nd["id"]] = (nd.get("properties") or {}).get("max_phase")
    targeted = set()
    for line in (RELEASE / "edges.jsonl").read_text().splitlines():
        if '"TARGETS"' not in line:
            continue
        e = json.loads(line)
        if e["predicate"] == "TARGETS":
            targeted.add(e["subject"])

    def _bucket(ph):
        return ("unknown" if ph is None else "approved" if ph == 4
                else "clinical" if ph in (2, 3) else "phase 1" if ph == 1
                else "preclinical")

    without = {k: 0 for k in ("approved", "clinical", "phase 1",
                              "preclinical", "unknown")}
    for cid, ph in phase.items():
        if cid not in targeted:
            without[_bucket(ph)] += 1
    routed = len(set(phase) & targeted)
    v["n_small_molecules"] = n(len(phase))
    v["n_routed_small_molecules"] = n(routed)
    v["routed_share"] = f"{round(100 * routed / len(phase), 1)}"
    for k in without:
        v[f"n_unrouted_{k.replace(' ', '')}"] = n(without[k])
    v["n_unrouted_ceiling"] = n(without["approved"] + without["clinical"])
    v["n_unrouted"] = n(len(phase) - routed)

    v["n_targets"] = n(n_targets)
    v["n_targets_undirected"] = n(n_undirected)
    v["targets_undirected_share"] = f"{100 * n_undirected / n_targets:.1f}%"
    n_undated = sum(undated.values())
    # WHAT THE DOSE-SERIES POLICY COSTS, both ways. Four seconds of
    # `make doc2`, and the alternative is a paragraph stating what an
    # exclusion would cost without anyone having run the exclusion.
    from kgav.labels import build_labels as _build
    _inc = _build(RELEASE)
    _exc = _build(RELEASE, dose_series_policy="exclude")

    def _tot(L, which):
        return sum(len(x) for x in getattr(L, which).values())

    v["ds_flagged"] = n(_inc.stats["outside_dose_series_kept"])
    v["ds_flagged_positive"] = n(_tot(_inc, "positives")
                                 - _tot(_exc, "positives"))
    v["ds_pos_include"] = n(_tot(_inc, "positives"))
    v["ds_pos_exclude"] = n(_tot(_exc, "positives"))
    v["ds_neg_include"] = n(_tot(_inc, "negatives"))
    v["ds_neg_exclude"] = n(_tot(_exc, "negatives"))
    # THE OVERLAP ACROSS THE CUTOFF, counted. A compound measured on
    # both sides of the split is the thing that makes a temporal
    # protocol not temporal, and the census that decides the policy is
    # cheap: two more build_labels calls, about four seconds.
    from kgav.temporal import drop_already_known as _dak
    _cut = int(_j(RESULTS / "temporal_structure_scorer.json")["cutoff"])
    _pre = _build(RELEASE, year_to=_cut)
    _post = _build(RELEASE, year_from=_cut + 1)
    _cen = {}
    for _virus in sorted(set(_post.positives) | set(_post.negatives)):
        bef = {"active": _pre.positives.get(_virus, set()),
               "inactive": _pre.negatives.get(_virus, set())}
        aft = {"active": _post.positives.get(_virus, set()),
               "inactive": _post.negatives.get(_virus, set())}
        seen = bef["active"] | bef["inactive"]
        for _c in sorted(aft["active"] | aft["inactive"]):
            if _c not in seen:
                continue
            # A SIDE'S POLARITY IS A SET, NOT A FLAG: a compound can be
            # in both the positive and the negative set of one side --
            # different cell lines, strains, labs -- and collapsing that
            # here would invent the label the loader refuses to invent.
            _cen[(_virus, _c)] = (
                frozenset(k for k, vv in bef.items() if _c in vv),
                frozenset(k for k, vv in aft.items() if _c in vv))
    _agree = sum(1 for b, a in _cen.values() if b == a)
    _stats, _ = _dak({k: set(x) for k, x in _post.positives.items()},
                     {k: set(x) for k, x in _post.negatives.items()},
                     _pre.positives, _pre.negatives, _cut)
    v["ov_pairs"] = str(len(_cen))
    v["ov_agree"] = str(_agree)
    v["ov_flip"] = str(len(_cen) - _agree)
    v["ov_dropped_pos"] = str(
        _stats.get("test_positive_already_known_dropped", 0))
    v["ov_dropped_neg"] = str(
        _stats.get("test_negative_already_known_dropped", 0))
    # M8'S EVIDENCE, gated here exactly as the driver enumerates it.
    # The channel's name says "fold similarity" and what it can actually
    # walk on this release is a different and weaker thing, so the four
    # counts are derived rather than described.
    import importlib
    _ing = importlib.import_module("ingest_folds")
    from kgav.folds import gate as _gate
    _raw = _j(ROOT / "data" / "raw" / "folds" / "phase2_chains.json")
    _chains = _ing.load_cache(_raw)
    _elig = [c for c in _chains
             if c.ok and c.spec.role in _ing.SCORE_ELIGIBLE_ROLES]
    _allowed = []
    for _i, _a in enumerate(_elig):
        for _b in _elig[_i + 1:]:
            if _a.spec.family != _b.spec.family:
                g = _gate(_a, _b)
                if g.allowed:
                    _allowed.append(g)
    _strong = [g for g in _allowed
               if (g.mode, g.strength) == ("nucleophile", "strong")]
    _weak = [g for g in _allowed if (g.mode, g.strength) == ("motif", "weak")]
    if len(_allowed) != len(_strong) + len(_weak):
        raise SystemExit(
            "an allowed fold pair is neither strong-nucleophile nor "
            "weak-motif; 08-results.md states only those two")
    _walk = []
    for line in (RELEASE / "edges.jsonl").read_text().splitlines():
        if '"FOLD_SIMILAR_TO"' not in line:
            continue
        e = json.loads(line)
        q = e.get("qualifiers") or {}
        if e["predicate"] == "FOLD_SIMILAR_TO" and q.get(
                "catalytic_type_match"):
            _walk.append(q)
    v["m8_allowed"] = str(len(_allowed))
    v["m8_strong"] = str(len(_strong))
    v["m8_weak"] = str(len(_weak))
    v["m8_strong_in_release"] = str(
        sum(1 for q in _walk if q.get("fold_strength") == "strong"))
    v["m8_walkable"] = str(len(_walk))
    v["n_undated_edges"] = n(n_undated)
    v["undated_share"] = pct(n_undated / dated)
    v["undated_table"] = table(
        ["source", "undated edges"],
        [[f"`{k}`", n(x)] for k, x in undated.most_common()], "lr")
    v["host_factor_table"] = table(
        ["virus", "HOST_FACTOR_FOR edges"],
        [[names.get(k, k), n(x)] for k, x in hf.most_common()], "lr")
    without = sorted(names[k] for k in names if k not in hf)
    v["n_taxa_without_host_factors"] = str(len(without))
    v["taxa_without_host_factors"] = ", ".join(without)
    v["pl_n_uncovered"] = str(len(pl["uncovered"]))

    # THE SCOPE TABLE: what the graph holds for each virus it declares.
    # This is the part a reader uses to decide what a result applies to,
    # and every figure in it moves whenever the graph is rebuilt. Its
    # predecessor's header said "generated" from the day it was written
    # and nothing generated it; the rows were maintained by hand.
    rows = []
    for e in reg["viruses"]:
        curie = f"NCBITaxon:{e['taxon']}"
        cs = e.get("chain_source") or {}
        rows.append([
            e["name"], e["taxon"], vr.family_of(e, reg) or "?",
            n(dep[curie]), n(res[curie]), n(av[curie]),
            ("species" if str(cs.get("taxon")) in (None, "None",
                                                   str(e["taxon"]))
             else f"strain {cs['taxon']}"),
            "yes" if dep[curie] else "**no**"])
    rows.sort(key=lambda r: -int(r[3].replace(",", "")))
    v["scope_table"] = table(
        ["virus", "taxon", "family", "dependency", "restriction",
         "antiviral", "chain source", "in the overlap table"],
        rows, "lllrrrll")
    v["n_dependency_factors"] = n(sum(dep.values()))
    v["n_restriction_factors"] = n(sum(res.values()))
    v["n_antiviral_edges"] = n(sum(av.values()))
    v["n_in_overlap"] = str(sum(1 for r in rows if r[7] == "yes"))
    # THE EXTREME CASE, NAMED RATHER THAN LEFT IN A NODE COUNT. One
    # virus is in the register and answers no question at all.
    thinnest = min(rows, key=lambda r: (int(r[3].replace(",", ""))
                                        + int(r[4].replace(",", ""))
                                        + int(r[5].replace(",", ""))))
    v["thinnest_virus"] = thinnest[0]
    v["thinnest_virus_edges"] = str(
        sum(int(thinnest[i].replace(",", "")) for i in (3, 4, 5)))

    # ----------------------------------------------------- label overlap
    lo = _j(RESULTS / "label_overlap.json")
    v["lo_both"] = n(lo["both"])
    v["lo_share_a_pub"] = n(lo["share_a_pub"])
    v["lo_share_pct"] = pct(lo["share_a_pub"] / lo["both_sides_have_pubs"])
    v["lo_unwithholdable"] = n(lo["unwithholdable"])
    v["lo_viral_proteins"] = n(lo["viral_proteins_reached"])
    # THE DENOMINATOR OF THE SHARE IS NOT THE POPULATION. `both` is
    # every compound carrying an INHIBITS edge AND a label; of those,
    # `both_sides_have_pubs` cite a publication on both sides, and the
    # share is over THAT -- which is the qualifier the figure was
    # quoted without.
    v["lo_both_pubs"] = n(lo["both_sides_have_pubs"])
    v["lo_summary_row"] = (
        f"{v['lo_share_pct']} shared PMIDs "
        f"({lo['share_a_pub']:,}/{lo['both_sides_have_pubs']:,})")
    v["lo_unwithholdable_pct"] = (
        f"{100 * lo['unwithholdable'] / lo['both']:.1f}%")

    # ------------------------------------------- degree weighting arms
    dw = _j(RESULTS / "degree_weighting_sensitivity.json")
    v["dw_cells"] = n(dw["cells"])
    v["dw_moved"] = str(dw["cells_moved"])
    v["dw_crossing"] = str(dw["cells_crossing_alpha"])
    v["dw_claim_cell"] = (f"{names.get(dw['claim_cell']['virus'], '?')} "
                          f"{dw['claim_cell']['channel']}")
    # THE SENSITIVITY, IN THE ORDER THAT KEEPS IT ONE. data/results holds
    # arm A, so arm A is the result and arm B is the sensitivity;
    # quoting B first would restate a figure no results file carries.
    cc = dw["claim_cell"]
    v["dw_p_a"] = f"{cc['p_a']:.6f}"
    v["dw_p_b"] = f"{cc['p_b']:.6f}"
    v["dw_lift_a"] = f"{cc['lift_a']:.2f}"
    v["dw_lift_b"] = f"{cc['lift_b']:.2f}"
    v["dw_hub_set"] = str(dw["hub_set_size"]["A"])
    readable_n = len([1 for per in mc["per_virus"].values()
                      for c, m in per.items()
                      if c not in ("degree", "COMBINED")
                      and m.get("hits_determined") and m.get("trustworthy")])
    v["dw_alpha_factor"] = str(round((0.05 / readable_n) / cc["p_b"]))
    v["dw_arms"] = table(
        ["arm", "what it changes"],
        [[f"`{k}`", x] for k, x in dw["arms"].items()])

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
    # THE SPREAD OF THE FLOOR IS THE ARGUMENT AGAINST A FLAT THRESHOLD,
    # so the two ends of it are derived rather than described. They are
    # what makes "shortlist anything above 0.55" wrong in both
    # directions at once.
    floors = {c["name"]: c["chemistry_filter"]["resemblance_floor"]["floor"]
              for c in rc["viruses"].values()
              if c["chemistry_filter"]["resemblance_floor"].get("floor")
              is not None}
    lo_v = min(floors, key=floors.get)
    hi_v = max(floors, key=floors.get)
    v["rc_floor_min"] = f"{floors[lo_v]:.4f}"
    v["rc_floor_min_virus"] = lo_v
    v["rc_floor_max"] = f"{floors[hi_v]:.4f}"
    v["rc_floor_max_virus"] = hi_v
    v["rc_floor_spread"] = f"{floors[hi_v] / floors[lo_v]:.1f}"
    # THE ROW THE REVIEW ASKED TO HAVE REMOVED, kept with its position.
    # The floor was not bent until the disliked row fell out, and the
    # only way to show that is to publish the row and the cell's own
    # standard beside each other.
    named = None
    for cell in rc["viruses"].values():
        for row in cell["suggestions"]:
            if (row.get("name") or "").upper() == "TEMAZEPAM":
                named = (cell, row)
                break
    if named:
        cell, row = named
        nm = (row.get("name") or "").upper()
        fl = cell["chemistry_filter"]["resemblance_floor"]["floor"]
        v["rc_named_note"] = (
            f"**The threshold was not bent until the disliked row fell "
            f"out.** The review of 2026-11-19 named {nm} on "
            f"{cell['name']} as tail noise. That cell's own failures have a "
            f"median of {fl:.4f}, so the cell's own standard keeps the row, "
            f"and it is still on the list with "
            f"{row['beats_failures']:.0%} of this virus's failures beside "
            f"it.")
    else:
        v["rc_named_note"] = (
            "The row the 2026-11-19 review named is no longer on any list; "
            "the floor that kept it has not changed.")

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
    contra = sorted(c.name for c in caps.values() if c.contradicted)
    v["contradicted_note"] = (
        f"**{len(contra)} "
        f"{'virus clears' if len(contra) == 1 else 'viruses clear'} "
        f"cross-sectionally and anti-rank "
        f"prospectively — {', '.join(contra)} — and are refused by name.** "
        f"A cross-sectional verdict alone does not license a shortlist, "
        f"because a shortlist is a prospective act."
        if contra else
        "No virus currently clears cross-sectionally while anti-ranking "
        "prospectively.")

    # ----------------------------------------------- how to name a virus
    # THE ALIASES ARE THE REGISTER'S, not a list retyped into a guide.
    # The predecessor of 11-user-guide.md listed them by hand, and a
    # reader who typed a spelling the register accepts but the page
    # omitted had no way to know the page was the thing that was wrong.
    reg = vr.load()
    names = vr.search_names(reg)
    dropped = vr.dropped_short_names(reg)
    v["n_virus_aliases"] = str(len(names))
    by_taxon: dict[str, list[str]] = {}
    for alias, taxon in names.items():
        by_taxon.setdefault(taxon, []).append(alias)
    rows = []
    for entry in reg["viruses"]:
        taxon = f"NCBITaxon:{entry['taxon']}"
        # A VIRUS WITH NO ACCEPTED SPELLING IS A BUG, NOT AN EMPTY CELL.
        # JEV is the near miss: its own name is the one string the
        # resolver refuses, so an earlier version of this keyed on the
        # name and dropped the row -- a guide silently listing fifteen
        # of sixteen viruses.
        aliases = sorted(by_taxon.get(taxon, ()), key=lambda s: (len(s), s))
        if not aliases:
            raise SystemExit(
                f"{entry['name']} ({taxon}) has no accepted spelling at "
                f"all; the resolver cannot be asked about it")
        no = dropped.get(taxon) or []
        rows.append([entry["name"],
                     ", ".join(f"`{a}`" for a in aliases),
                     ", ".join(f"`{a}`" for a in sorted(no)) or "—"])
    v["virus_naming_table"] = table(
        ["virus", "write it as (any of these)", "deliberately not accepted"],
        rows)
    n_dropped = sum(len(x) for x in dropped.values())
    v["virus_naming_note"] = (
        f"**{n_dropped} short "
        f"{'name is' if n_dropped == 1 else 'names are'} refused on "
        f"purpose**, in the last "
        f"column. A name shorter than {vr.MIN_SEARCH_NAME} characters, or "
        f"one that is a substring of ordinary English, matches text that "
        f"is not a virus — and a wrong match is worse than no match, "
        f"because it answers confidently about the wrong organism."
        if n_dropped else
        "Every name in the register is accepted; none is short enough to "
        "need refusing.")

    # ------------------------------------------------------------- agent
    v["agent_model"] = f"`{llm.DEFAULT_MODEL}`"
    # THE VIRUS AN UNRECOGNISED NAME FALLS BACK TO, by its human name.
    # Printing the taxon id would be accurate and useless to the reader
    # this page is for, and hard-coding "SARS-CoV-2" would survive the
    # fallback being changed.
    fallback = [e["name"] for e in reg["viruses"]
                if f"NCBITaxon:{e['taxon']}" == guards.SUPPORTED_VIRUS]
    if len(fallback) != 1:
        raise SystemExit(
            f"guards.SUPPORTED_VIRUS is {guards.SUPPORTED_VIRUS}, which "
            f"matches {len(fallback)} register entries; 11-user-guide.md "
            f"has to name it")
    v["default_virus"] = fallback[0]
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


def _first_heading(path: Path) -> str:
    """A registration's subject, out of its own title line.

    `# Pre-registration: <subject>` is the form every one of them uses.
    Falling back to the file name would make a registration that broke
    the convention look described when it was not, so this raises.
    """
    first = path.read_text().splitlines()[0].strip()
    m = re.match(r"#\s*Pre-registration:\s*(.+)$", first)
    if not m:
        raise SystemExit(
            f"{path.name} does not open with "
            f"'# Pre-registration: <subject>'; its first line is {first!r}")
    return m.group(1).strip()


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
    # A COPY TO CHECK AGAINST, so a test of the staleness check does not
    # have to mutate a tracked file that every other test is reading.
    # The predecessor of this flag did not exist and the mutation test
    # edited doc2/01-overview.md in place: under `pytest -n auto` another
    # worker reads that file mid-edit, which is the same race the
    # rewrite-tables gate had and the same repair.
    ap.add_argument("--out", type=Path, default=OUT,
                    help="the directory to compare against or write into "
                         "(default: doc2/)")
    args = ap.parse_args(argv)
    out_dir = args.out

    if not TEMPLATES.exists():
        print(f"no templates in {TEMPLATES}", file=sys.stderr)
        return 2
    files = render(derive())

    if args.check:
        bad = []
        for name, text in sorted(files.items()):
            path = out_dir / name
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

    out_dir.mkdir(parents=True, exist_ok=True)
    for name, text in sorted(files.items()):
        (out_dir / name).write_text(text)
        print(f"wrote {out_dir.name}/{name}  "
              f"({len(text.splitlines()):,} lines)")
    return 0


def _first_diff(have: str, want: str) -> str:
    a, b = have.splitlines(), want.splitlines()
    for i, (x, y) in enumerate(zip(a, b), 1):
        if x != y:
            return f"line {i} is {x[:60]!r}, should be {y[:60]!r}"
    return f"{len(a)} lines on disk, {len(b)} generated"


if __name__ == "__main__":
    raise SystemExit(main())
