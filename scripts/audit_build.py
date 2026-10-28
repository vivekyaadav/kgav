#!/usr/bin/env python3
"""What has been ingested, what has been built, and what is missing.

    python scripts/audit_build.py
    python scripts/audit_build.py --virus NCBITaxon:11053

WHY THIS EXISTS. The facts needed to decide what to work on next are
spread across five places: pipeline.LAYERS says what a layer needs,
data/raw says what arrived, data/releases/v0.1-* says what was built,
the assembled release says what joined, and the register says which
viruses were supposed to be covered. Answering "can the agent handle
DENV-1 yet" meant reading all five and holding the join in your head,
which is how a virus with a taxon in the map and no file on disk went
unnoticed for weeks.

EVERY NUMBER HERE IS READ, NONE IS REMEMBERED. The per-virus table
counts nodes and edges in the release on disk at the moment it runs. A
figure quoted from a previous run is a figure about a previous graph,
and this project has already shipped two reports that did exactly that.

WHAT IT IS NOT. It does not decide whether a gap is worth closing, and
it does not say a virus is "ready": a channel can be structurally
reachable and still rest on four edges. The counts are there so that
judgement has something to stand on.

Read-only. Writes nothing. Degrades per section: a missing release is
reported and the rest still runs.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.pipeline import LAYERS, availability, layer_dir
from kgav.provenance import check_file
from kgav.schema import load_schema, print_capped
from kgav.virus_register import declared_taxa, family_of

VIRAL_PPI = "PHYSICALLY_INTERACTS_WITH"
HOST_FACTOR = "HOST_FACTOR_FOR"
INHIBITS = "INHIBITS"
# THE OTHER CHEMISTRY EDGE. chembl.py writes four predicates and this
# counted one. A cell-based assay measuring "compound inhibits yellow
# fever virus" has no protein to point at, so it lands on the
# OrganismTaxon instead -- and Yellow fever has 1,070 ChEMBL compounds
# while reading zero here. "No direct-acting compound" was true and
# "no compound" was what it looked like.
ANTIVIRAL = "HAS_ANTIVIRAL_ACTIVITY_AGAINST"
INACTIVE = "MEASURED_INACTIVE_AGAINST"
CAUSES = "CAUSES"
TARGETS = "TARGETS"


def _count_lines(path: Path) -> int:
    if not path.exists():
        return 0
    with open(path, "rb") as fh:
        return sum(1 for ln in fh if ln.strip())


def raw_inputs(root: Path) -> list[tuple[str, int, bool, str]]:
    """(layer, phase, runnable, why not) for every declared layer."""
    out = []
    for layer in LAYERS:
        ok, why = availability(layer, root)
        out.append((layer.name, layer.phase, ok, why))
    return out


def layer_releases(root: Path) -> list[tuple[str, int, int, bool]]:
    """(layer, nodes, edges, built) per declared layer.

    A layer release with zero edges is reported rather than skipped: it
    is the shape a layer takes when its ingest ran and joined nothing,
    which looks identical to "not built" in a directory listing and
    needs completely different work.
    """
    out = []
    for layer in LAYERS:
        d = layer_dir(root, layer.name)
        built = (d / "edges.jsonl").exists() or (d / "nodes.jsonl").exists()
        out.append((layer.name,
                    _count_lines(d / "nodes.jsonl"),
                    _count_lines(d / "edges.jsonl"),
                    built))
    return out


def assembled(release: Path) -> dict:
    """Node classes and edge predicates in the assembled release."""
    classes: Counter = Counter()
    preds: Counter = Counter()
    sources: Counter = Counter()
    nodes: dict[str, dict] = {}
    nf, ef = release / "nodes.jsonl", release / "edges.jsonl"
    if nf.exists():
        for line in nf.read_text().splitlines():
            if not line.strip():
                continue
            n = json.loads(line)
            nodes[n["id"]] = n
            classes[n["class"]] += 1
    edges = []
    if ef.exists():
        for line in ef.read_text().splitlines():
            if not line.strip():
                continue
            e = json.loads(line)
            edges.append(e)
            preds[e["predicate"]] += 1
            sources[e.get("primary_knowledge_source") or "?"] += 1
    return {"nodes": nodes, "edges": edges, "classes": classes,
            "predicates": preds, "sources": sources}


def per_virus(a: dict, taxa: list[str]) -> dict[str, Counter]:
    """What the graph actually holds for each declared virus.

    THE JOIN THAT WAS BEING DONE BY HAND. A virus needs proteins to be
    a target, chains to be a SPECIFIC target, host PPI to bridge to the
    interactome, and dependency factors to reach the CRISPR channels.
    Having three of the four is not visible from any single layer's
    output, and it decides which questions the agent can answer.
    """
    nodes = a["nodes"]
    prot_of: dict[str, set[str]] = defaultdict(set)
    for nid, n in nodes.items():
        if n.get("class") != "Protein":
            continue
        p = n.get("properties") or {}
        if p.get("is_viral") and p.get("taxon_id"):
            prot_of[str(p["taxon_id"])].add(nid)

    out: dict[str, Counter] = {t: Counter() for t in taxa}
    for t in taxa:
        out[t]["proteins"] = len(prot_of.get(t, ()))
        out[t]["chains"] = sum(1 for p in prot_of.get(t, ())
                               if "PRO_" in p)

    def viral_taxon(nid: str) -> str | None:
        n = nodes.get(nid) or {}
        p = n.get("properties") or {}
        return str(p["taxon_id"]) if p.get("is_viral") and p.get("taxon_id") \
            else None

    partners: dict[str, set[str]] = defaultdict(set)
    for e in a["edges"]:
        pred, s, o = e["predicate"], e["subject"], e["object"]
        if pred == VIRAL_PPI:
            for viral, host in ((s, o), (o, s)):
                t = viral_taxon(viral)
                if t and t in out and not viral_taxon(host):
                    out[t]["vh_ppi_edges"] += 1
                    partners[t].add(host)
        elif pred == HOST_FACTOR:
            if o in out:
                d = (e.get("qualifiers") or {}).get("direction") or "unstated"
                out[o][f"hf_{d}"] += 1
        elif pred == INHIBITS:
            t = viral_taxon(o)
            if t and t in out:
                out[t]["inhibits"] += 1
        elif pred == INACTIVE:
            t = viral_taxon(o)
            if t and t in out:
                out[t]["inactive"] += 1
        elif pred == ANTIVIRAL and o in out:
            # Object is the TAXON, not a protein: that is the whole
            # point of the predicate.
            out[o]["antiviral"] += 1
        elif pred == CAUSES and s in out:
            out[s]["causes_disease"] += 1
    for t, ps in partners.items():
        out[t]["host_partners"] = len(ps)
    return out


def measured_reach(results: Path, release: Path) -> dict[str, int] | None:
    """Held-out positives each virus's metapaths actually reach.

    `{virus label: COMBINED.pos}`, or None when that cannot be known.

    WHY THIS IS NOT DERIVED FROM THE GRAPH. Whether a metapath reaches a
    labelled compound is the result of walking it against a label set,
    which is what compare_metapaths.py does. Edge counts cannot answer
    it: DENV-1 has 349 labels and 53 virus-host interactions and reaches
    none of the 349.

    ONLY IF THE FILE IS CURRENT. provenance.check_file is the same audit
    `make check-results` runs, so a comparison measured on a graph that
    has since moved is reported as unmeasured rather than quoted here.
    Returning None for "cannot be known" rather than 0 for "reaches
    nothing" is the whole point: those are opposite claims, and this
    script has twice shipped a figure from a previous graph.
    """
    path = Path(results) / "metapath_comparison.json"
    if not path.exists():
        return None
    status, _msg = check_file(path, release)
    if status != "ok":
        return None
    try:
        per = json.loads(path.read_text())["per_virus"]
    except (OSError, ValueError, KeyError):
        return None
    out = {}
    for virus, channels in per.items():
        pos = (channels.get("COMBINED") or {}).get("pos")
        if isinstance(pos, int):
            out[virus] = pos
    return out or None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", type=Path, default=ROOT / "config/viruses.yaml")
    ap.add_argument("--release", type=Path,
                    default=ROOT / "data/releases/v0.1")
    ap.add_argument("--root", type=Path, default=ROOT)
    # EXPLICIT, not derived from --release. release.parent.parent happens
    # to give data/results for the default layout and walks out of the
    # tree for any other, which is how a test would silently exercise a
    # path no run uses.
    ap.add_argument("--results", type=Path, default=ROOT / "data/results")
    args = ap.parse_args(argv)

    doc = yaml.safe_load(args.config.read_text())
    register = {str(v["taxon"]): v for v in doc["viruses"]}
    taxa = [f"NCBITaxon:{t}" for t in register]
    label = {f"NCBITaxon:{t}": v["name"] for t, v in register.items()}
    # DEFAULTS APPLIED. Reading the per-entry key alone printed "?" for
    # every coronavirus, because they inherit Coronaviridae from the
    # register's defaults block -- and the audit's own family column was
    # the fourth copy of that bug.
    family = {f"NCBITaxon:{t}": (family_of(v, doc) or "?")
              for t, v in register.items()}

    print("=" * 74)
    print("BUILD AUDIT")
    print("=" * 74)
    print(f"register: {len(register)} viruses, "
          f"{len(declared_taxa(doc))} declared taxids")

    # ------------------------------------------------- 1. raw inputs
    print(f"\n{'-' * 74}\n1. RAW INPUTS — can each layer run from what is on "
          f"disk?\n{'-' * 74}")
    rows = raw_inputs(args.root)
    ready = [r for r in rows if r[2]]
    print(f"  {len(ready)} of {len(rows)} layers can run\n")
    print(f"    {'layer':<14} {'phase':>5}  {'status'}")
    for name, phase, ok, why in rows:
        print(f"    {name:<14} {phase:>5}  "
              + ("ready" if ok else f"BLOCKED: {why}"))

    # --------------------------------------------- 2. layer releases
    print(f"\n{'-' * 74}\n2. LAYERS BUILT — what each ingest actually "
          f"wrote\n{'-' * 74}")
    print(f"    {'layer':<14} {'nodes':>10} {'edges':>10}  note")
    unbuilt = []
    for name, n, e, built in layer_releases(args.root):
        if not built:
            unbuilt.append(name)
            note = "NOT BUILT"
        elif not e:
            note = "built, no edges — ran and joined nothing"
        else:
            note = ""
        print(f"    {name:<14} {n:>10,} {e:>10,}  {note}")

    # ------------------------------------------ 3. assembled release
    print(f"\n{'-' * 74}\n3. ASSEMBLED RELEASE — {args.release}\n{'-' * 74}")
    if not (args.release / "nodes.jsonl").exists():
        print("  NOT ASSEMBLED. Nothing below can be computed.")
        return 1
    a = assembled(args.release)
    print(f"  {len(a['nodes']):,} nodes, {len(a['edges']):,} edges\n")
    print("  nodes by class")
    for k, v in a["classes"].most_common():
        print(f"    {k:<22} {v:>10,}")
    print("\n  edges by predicate")
    for k, v in a["predicates"].most_common():
        print(f"    {k:<32} {v:>10,}")
    print("\n  edges by source")
    print_capped(a["sources"].most_common(), 8,
                 lambda kv: f"{kv[0]:<32} {kv[1]:>10,}",
                 unit="more source(s)")
    man = args.release / "MANIFEST.json"
    if man.exists():
        m = json.loads(man.read_text())
        layers = m.get("layers") or {}
        print(f"\n  manifest: {len(layers)} layer(s) assembled"
              + (f", hub cut degree {m['hub_cut_degree']}"
                 if "hub_cut_degree" in m else ""))

    # --------------------------------------------- 4. per-virus matrix
    print(f"\n{'-' * 74}\n4. PER-VIRUS COVERAGE — what the graph holds for "
          f"each declared virus\n{'-' * 74}")
    pv = per_virus(a, taxa)
    hdr = (f"    {'virus':<14} {'fam':<5} {'prot':>5} {'chain':>6} "
           f"{'vhPPI':>7} {'partn':>6} {'dep':>6} {'restr':>6} {'inhib':>6} "
           f"{'antiv':>6} {'inact':>6} {'dis':>4}")
    print(hdr)
    for t in sorted(taxa, key=lambda x: (family[x], label[x])):
        c = pv[t]
        print(f"    {label[t]:<14} {family[t][:4]:<5} {c['proteins']:>5,} "
              f"{c['chains']:>6,} {c['vh_ppi_edges']:>7,} "
              f"{c['host_partners']:>6,} {c['hf_dependency']:>6,} "
              f"{c['hf_restriction']:>6,} {c['inhibits']:>6,} "
              f"{c['antiviral']:>6,} {c['inactive']:>6,} "
              f"{c['causes_disease']:>4,}")
    print("\n    prot/chain  viral proteins and mature peptides in the spine")
    print("    vhPPI       virus-host edges; partn, distinct host proteins")
    print("    dep/restr   CRISPR host factors, by direction")
    print("    inhib       compounds with an INHIBITS edge to a viral protein")
    print("    antiv       organism-level antivirals: a cell-based assay has")
    print("                no protein to point at, so the edge lands on the")
    print("                TAXON. A virus can have none of the first and")
    print("                hundreds of these.")
    print("    inact       measured and found INACTIVE; evidence, not absence")

    print("    dis         CAUSES edges to a Disease node")

    # ---------------------------------------- 4b. labels versus features
    # THE DISTINCTION THE TABLE ABOVE DOES NOT MAKE. antiv is
    # held_out: true and evaluation_label: true in the schema, and
    # appears in NO metapath -- it is the ground truth being predicted,
    # not evidence to predict from. Yellow fever has 1,030 of them and
    # no PPI, no CRISPR and no INHIBITS, so nothing can reach a
    # compound for it: every prediction is a miss by construction, and
    # the row above reads like its best-covered column.
    reach = measured_reach(args.results, args.release)
    held = load_schema().held_out_predicates()
    print(f"\n    LABELS vs FEATURES\n    held out of every metapath: "
          f"{', '.join(sorted(held))}")
    print("    `labels` counts POSITIVES only "
          "(HAS_ANTIVIRAL_ACTIVITY_AGAINST).\n    "
          "MEASURED_INACTIVE_AGAINST is the negative label and is held out "
          "too,\n    so the `inact` column above is also answer sheet, not "
          "evidence.")
    print(f"\n    {'virus':<14} {'labels':>7} {'feature routes':<34} "
          f"verdict")
    verdicts: dict[str, list[str]] = defaultdict(list)
    for t in sorted(taxa, key=lambda x: (family[x], label[x])):
        c = pv[t]
        # WITH THEIR SIZE. "evaluable" on seven PPI edges is true and
        # reads like "evaluable" on eight thousand. No threshold is
        # invented here: the number is shown and the reader judges.
        routes = [f"{n}({v:,})" for n, v in (("inhibits", c["inhibits"]),
                                             ("crispr", c["hf_dependency"]),
                                             ("ppi", c["vh_ppi_edges"])) if v]
        # MEASURED REACH, WHERE IT IS AVAILABLE AND CURRENT. A route
        # that exists and reaches no labelled compound is invisible in
        # edge counts: DENV-1 carries 349 labels and 53 PPI edges and
        # this column read "evaluable", while the scorer reached 0 of
        # those 349. Seven of fifteen viruses are in that state.
        hit = reach.get(label[t]) if reach is not None else None
        if c["antiviral"] and routes and hit == 0:
            verdict = "ROUTES REACH NO LABEL"
        elif c["antiviral"] and routes and reach is None:
            verdict = "evaluable (reach unmeasured)"
        elif c["antiviral"] and routes:
            verdict = "evaluable"
        elif c["antiviral"]:
            verdict = "LABELS, NO FEATURE ROUTE"
        elif routes:
            verdict = "features, nothing to score"
        else:
            verdict = "neither"
        verdicts[verdict].append(label[t])
        print(f"    {label[t]:<14} {c['antiviral']:>7,} "
              f"{'+'.join(routes) or '-':<34} {verdict}")
    for v in ("LABELS, NO FEATURE ROUTE", "features, nothing to score",
              "ROUTES REACH NO LABEL"):
        if verdicts[v]:
            print(f"\n    {v}: {', '.join(sorted(verdicts[v]))}")
    if verdicts["LABELS, NO FEATURE ROUTE"]:
        print("      No metapath terminates at these, so every prediction "
              "is a miss by\n      construction. Scoring them measures the "
              "graph's reach, not the model.")
    if verdicts["ROUTES REACH NO LABEL"]:
        print("      A route exists and walks it, and reaches none of the "
              "held-out\n      positives. A ranking for these is complete, "
              "and nothing in this\n      graph can mark it right or wrong. "
              "Edge counts cannot see this,\n      which is why the column "
              "beside them read `evaluable`.")
    if reach is None:
        print("\n      reach unmeasured: metapath_comparison.json is "
              "absent, stale or\n      unstamped for this release, so "
              "whether a route reaches a label is\n      unknown rather "
              "than assumed. `make check-results` says which.")

    # --------------------------------------------------- 5. what is left
    print(f"\n{'-' * 74}\n5. WHAT IS MISSING\n{'-' * 74}")
    blocked = [(n, w) for n, _p, ok, w in rows if not ok]
    if blocked:
        print(f"  {len(blocked)} layer(s) cannot run:")
        for n, w in blocked:
            print(f"    {n:<14} {w}")
    if unbuilt:
        print(f"  {len(unbuilt)} layer(s) never built: {', '.join(unbuilt)}")

    gaps = defaultdict(list)
    for t in taxa:
        c = pv[t]
        if not c["proteins"]:
            gaps["no viral protein in the spine"].append(label[t])
        if not c["vh_ppi_edges"]:
            gaps["no virus-host PPI"].append(label[t])
        if not c["hf_dependency"]:
            gaps["no CRISPR dependency factor"].append(label[t])
        if not c["inhibits"] and not c["antiviral"]:
            gaps["no compound of any kind"].append(label[t])
        elif not c["inhibits"]:
            gaps["no protein-level inhibitor (organism-level only)"].append(
                label[t])
        if not c["causes_disease"]:
            gaps["no disease link"].append(label[t])
    for why, names in sorted(gaps.items(), key=lambda kv: -len(kv[1])):
        print(f"  {len(names):>2} virus(es) with {why}:")
        print(f"       {', '.join(sorted(names))}")
    if not gaps:
        print("  No per-virus gap in any measured axis.")
    print("\n  A gap here is a fact about the SOURCES, not a defect. Whether "
          "it is\n  worth closing depends on what the release is for, and "
          "this script\n  does not decide that.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
