#!/usr/bin/env python3
"""What the schema declares against what the release actually holds.

    make data-audit
    python scripts/data_audit.py --release data/releases/v0.1

WHY A SCRIPT RATHER THAN A ONE-OFF READ. The 2026-10-24 audit found four
declared node classes with no nodes and ten predicates with no edges,
including BOTH predicates the schema marks evaluation_label -- so the
project's only label was HAS_ANTIVIRAL_ACTIVITY_AGAINST and nothing said
so. None of that is visible in assemble_release's output, which lists what
IS there: a class with no nodes prints no line, and a reader counting
lines sees a complete graph.

Everything here is derived from schema/kg_schema.yaml and the release on
disk. Nothing is hardcoded, so a predicate added to the schema and never
ingested shows up here on the next run rather than when someone notices.

THE ONE EXCEPTION IS `NOTES`, which records what to DO about a gap. That
is judgement, not derivation, so each entry carries the date it was formed
and the audit prints it as such. It is here rather than in a commit
message because a remedy worked out once and lost is a remedy worked out
twice.

WHAT THIS DOES NOT DO. It does not validate the graph -- scripts/
validate_release.py does that and is a different question ("is what is
here well-formed") from this one ("is what is here everything"). It does
not judge whether a count is big enough, only whether it is zero, except
for enum qualifiers where the share of `unknown` is reported because a
required qualifier can be present on every edge and still say nothing.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

# RECORDED JUDGEMENTS, not derived. Each says what would fill a gap and
# when that was worked out. Keyed by the class or predicate it concerns.
NOTES: dict[str, tuple[str, str]] = {
    "CellLine": ("2026-10-24", (
        "No nodes, and none expected: cell lines are stored as a `cell_line` "
        "qualifier holding a CVCL CURIE. The class is vestigial unless "
        "something needs to hang data off a cell line itself.")),
    "Publication": ("2026-10-24", (
        "No nodes, and none expected: every edge carries `publications` as "
        "a list of PMID CURIEs. Vestigial in the same way as CellLine.")),
    "VirusStrain": ("2026-10-24", (
        "NOT vestigial the way CellLine is, despite the qualifier above. "
        "HAS_ANTIVIRAL_ACTIVITY_AGAINST records a strain as free text -- "
        "'USA-WA1/2020' -- which names it and cannot be joined on, so "
        "CONFERS_RESISTANCE_TO has nothing to hang off and the resistance "
        "axis is absent. The qualifier is where a future ingest would get "
        "the strain names; it is not a substitute for the nodes. "
        "Stanford's CoV-RDB covers SARS-CoV-2 resistance; flaviviruses "
        "have far less.")),
    "Variant": ("2026-10-24", (
        "Genuinely absent, and with it HAS_VARIANT and OCCURS_IN. Needed "
        "before any 'this mutation confers resistance' claim.")),
    "TREATS": ("2026-10-24", (
        "EMPTY, and the schema marks it evaluation_label. ChEMBL's "
        "drug_indication table fills it (max_phase_for_ind 4) and is "
        "already on disk in chembl_37.db -- the layer queries activities "
        "and assays and never touches it. ChEMBL indexes by MeSH/EFO while "
        "Disease nodes are MONDO, and label matching is banned in "
        "normalize_bio, so it needs EFO->MONDO xrefs. Blocked behind the "
        "Disease axis: there is nothing for an indication to attach to.")),
    "IN_TRIAL_FOR": ("2026-10-24", (
        "EMPTY, and the schema marks it evaluation_label. Same source: "
        "drug_indication at max_phase_for_ind 1-3, with indication_refs "
        "supplying the nct_id this predicate requires.")),
    "CAUSES": ("2026-10-24", (
        "Three edges for sixteen viruses. 13 of 16 register entries carry "
        "mondo: null, so there is no Disease node to point at. This is the "
        "prerequisite for TREATS and IN_TRIAL_FOR and the cheapest of the "
        "three -- curation, no download.")),
    "SUBCLASS_OF": ("2026-10-24", (
        "Empty. NCBI Taxonomy gives parent-child cheaply, and it is what "
        "would let a question be asked at the Flaviviridae level instead "
        "of once per species.")),
    "SUBCLASS_OF_DISEASE": ("2026-10-24", (
        "Empty, and follows the Disease axis.")),
    "COMBINED_WITH": ("2026-10-24", (
        "Empty. Drug synergy, with no obvious redistributable source; "
        "NCATS matrix screens are the nearest.")),
    "PREDICTED_INDICATION": ("2026-10-24", (
        "Empty BY DESIGN -- is_model_output. It is written by a scorer, "
        "and no scorer is adopted.")),
}


def snake(name: str) -> str:
    return re.sub(r"(?<!^)(?=[A-Z])", "_", name).lower()


def load_release(release: Path) -> tuple[Counter, Counter, dict, dict, set]:
    """One pass: node classes, predicates, qualifier values, taxon reach."""
    classes: Counter = Counter()
    preds: Counter = Counter()
    quals: dict[str, dict[str, Counter]] = defaultdict(
        lambda: defaultdict(Counter))
    # EVERY TAXON NODE, not only those an edge touches. Collected here so
    # a virus with no edge at all appears in the coverage table with a row
    # of zeros -- which is the single most useful line in it, and the one
    # a table built from edges alone cannot contain.
    taxa: dict[str, Counter] = defaultdict(Counter)
    edge_fields: set[str] = set()

    with open(release / "nodes.jsonl") as fh:
        for line in fh:
            if not line.strip():
                continue
            n = json.loads(line)
            classes[n.get("class", "?")] += 1
            if n.get("class") == "OrganismTaxon":
                taxa.setdefault(n["id"], Counter())
    with open(release / "edges.jsonl") as fh:
        for line in fh:
            if not line.strip():
                continue
            e = json.loads(line)
            p = e.get("predicate", "?")
            preds[p] += 1
            edge_fields |= set(e)
            for k, v in (e.get("qualifiers") or {}).items():
                if isinstance(v, str):
                    quals[p][k][v] += 1
            for end in (e.get("subject"), e.get("object")):
                if isinstance(end, str) and end.startswith("NCBITaxon:"):
                    taxa[end][p] += 1
    return classes, preds, quals, taxa, edge_fields


def note_for(key: str) -> str:
    if key not in NOTES:
        return ""
    when, text = NOTES[key]
    return f"        [{when}] {text}"


def wrap(text: str, width: int = 74, indent: str = "        ") -> str:
    out, line = [], indent
    for word in text.split():
        if len(line) + len(word) + 1 > width and line.strip():
            out.append(line)
            line = indent + word
        else:
            line = f"{line} {word}" if line.strip() else indent + word
    if line.strip():
        out.append(line)
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path,
                    default=ROOT / "data/releases/v0.1")
    ap.add_argument("--schema", type=Path,
                    default=ROOT / "schema/kg_schema.yaml")
    ap.add_argument("--unknown-threshold", type=float, default=0.5,
                    help="flag an enum qualifier whose `unknown` share "
                         "exceeds this")
    args = ap.parse_args(argv)

    if not (args.release / "edges.jsonl").exists():
        print(f"no release at {args.release}")
        return 2
    schema = yaml.safe_load(args.schema.read_text())
    classes, preds, quals, taxa, edge_fields = load_release(args.release)
    gaps = 0

    print(f"DATA AUDIT — {args.release.name}, schema "
          f"{schema.get('schema_version')}")
    print(f"{sum(classes.values()):,} nodes, {sum(preds.values()):,} edges\n")

    # ------------------------------------------------------ node classes
    declared_c = list(schema["node_classes"])
    print(f"NODE CLASSES — {len(declared_c)} declared, "
          f"{sum(1 for c in declared_c if classes.get(c))} populated")
    for c in declared_c:
        n = classes.get(c, 0)
        if n:
            print(f"  {c:<16}{n:>10,}")
    for c in declared_c:
        if classes.get(c):
            continue
        gaps += 1
        # Is its data kept somewhere other than a node? A class whose
        # snake_case name is an edge qualifier or an edge field is stored
        # as a property, which is a design choice rather than a gap.
        as_qual = sorted({p for p, kk in quals.items() if snake(c) in kk})
        as_field = [f for f in edge_fields
                    if f in (snake(c), snake(c) + "s")]
        print(f"  {c:<16}{'EMPTY':>10}")
        if as_qual:
            print(f"        stored as a `{snake(c)}` qualifier on "
                  f"{len(as_qual)} predicate(s): {', '.join(as_qual[:4])}")
        if as_field:
            print(f"        stored as the edge field `{as_field[0]}`")
        if note_for(c):
            print(wrap(NOTES[c][1], indent="        ")
                  .replace("        ", f"        [{NOTES[c][0]}] ", 1))
    for c in sorted(set(classes) - set(declared_c)):
        gaps += 1
        print(f"  {c:<16}{classes[c]:>10,}  <- IN THE RELEASE, NOT IN THE "
              f"SCHEMA")

    # -------------------------------------------------------- predicates
    edges = {e["predicate"]: e for e in schema["edge_classes"]}
    print(f"\nPREDICATES — {len(edges)} declared, "
          f"{sum(1 for p in edges if preds.get(p))} populated")
    for p in sorted(edges, key=lambda x: -preds.get(x, 0)):
        if preds.get(p):
            print(f"  {p:<32}{preds[p]:>10,}")
    for p, spec in edges.items():
        if preds.get(p):
            continue
        why = []
        if spec.get("is_model_output"):
            why.append("model output")
        if spec.get("evaluation_label"):
            why.append("EVALUATION LABEL")
        for end in ("subject", "object"):
            cls = spec.get(end)
            if cls in declared_c and not classes.get(cls):
                why.append(f"{end} class {cls} is empty")
        if not spec.get("is_model_output"):
            gaps += 1
        tail = f"  ({'; '.join(why)})" if why else ""
        print(f"  {p:<32}{'EMPTY':>10}{tail}")
        if p in NOTES:
            print(wrap(NOTES[p][1], indent="        ")
                  .replace("        ", f"        [{NOTES[p][0]}] ", 1))
    for p in sorted(set(preds) - set(edges)):
        gaps += 1
        print(f"  {p:<32}{preds[p]:>10,}  <- IN THE RELEASE, NOT IN THE "
              f"SCHEMA")

    # -------------------------------------------------- evaluation labels
    labels = [p for p, e in edges.items() if e.get("evaluation_label")]
    live = [p for p in labels if preds.get(p)]
    print(f"\nEVALUATION LABELS — {len(live)} of {len(labels)} have data")
    for p in labels:
        n = preds.get(p, 0)
        print(f"  {p:<32}{(f'{n:,}' if n else 'EMPTY'):>10}")
    if not live:
        print(wrap(
            "NO label carries data. Nothing in this release can be "
            "evaluated against anything.", indent="  "))
    elif len(live) == 1:
        print(wrap(
            f"Only {live[0]} carries data, so every evaluation in this "
            "project rests on it. A second, independent label is what "
            "would show whether a result is about the biology or about "
            "that one source.", indent="  "))

    # ------------------------------------------------------- the channels
    print("\nMETAPATHS")
    for name, mp in schema["metapaths"].items():
        path = [str(x) for x in mp["path"]]
        missing = [x.lstrip("<>") for i, x in enumerate(path)
                   if i % 2 == 1 and not preds.get(x.lstrip("<>"))]
        if missing:
            gaps += 1
            print(f"  {name:<5} DEAD  needs {', '.join(sorted(set(missing)))}")
        else:
            print(f"  {name:<5} walks")

    # ----------------------------------------- a required enum that says nothing
    print(f"\nENUM QUALIFIERS over {args.unknown_threshold:.0%} unknown")
    flagged = 0
    for p, spec in edges.items():
        for q, qspec in (spec.get("qualifiers") or {}).items():
            if not isinstance(qspec, dict) or "enum" not in qspec:
                continue
            counts = quals.get(p, {}).get(q)
            if not counts:
                continue
            total = sum(counts.values())
            unknown = counts.get("unknown", 0)
            if total and unknown / total > args.unknown_threshold:
                flagged += 1
                gaps += 1
                req = ("required" if q in (spec.get("required_qualifiers")
                                           or []) else "optional")
                print(f"  {p}.{q:<22}{unknown:>9,}/{total:<9,} "
                      f"{unknown / total:>6.1%} unknown  ({req})")
    if not flagged:
        print("  none")
    else:
        print(wrap(
            "A required qualifier present on every edge and set to "
            "`unknown` passes validation and carries no information. This "
            "is the gap a curated source closes.", indent="  "))

    # ------------------------------------------------- per-taxon coverage
    print("\nPER-TAXON COVERAGE — predicates reaching each built virus")
    want = [p for p, e in edges.items()
            if "OrganismTaxon" in (e.get("subject"), e.get("object"))
            and not e.get("is_model_output")]
    print(f"  {'taxon':<20}" + "".join(f"{p[:11]:>13}" for p in want))
    for t in sorted(taxa, key=lambda x: -sum(taxa[x].values())):
        row = "".join(f"{(taxa[t].get(p) or 0):>13,}" for p in want)
        empty = sum(1 for p in want if not taxa[t].get(p))
        flag = "   <- no edge on " + str(empty) + " of these" if empty else ""
        print(f"  {t:<20}{row}{flag}")

    print(f"\n{gaps} gap(s). A gap is a declared thing with nothing in it, "
          f"or a\nrequired enum that mostly says `unknown` -- not a judgement "
          f"that a\ncount is too small.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
