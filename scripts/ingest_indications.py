#!/usr/bin/env python3
"""Driver: drug -> Disease TREATS / IN_TRIAL_FOR from ChEMBL indications.

    python scripts/ingest_indications.py --inspect   # counts, no write
    python scripts/ingest_indications.py             # build the layer

REGISTERED IN ADVANCE by doc2/registrations/prereg_clinical_labels.md, signed
2026-11-17 before this script existed. The registration counted the
operands -- 578 placeable (compound, disease) pairs over 294 compounds
already in the graph -- and fixed the analysis that will read them: one
evaluable cell, COVID-19 against SARS-CoV-2, family of one, alpha 0.05.
This script ingests; it computes no AUC and ranks nothing.

--inspect first, always. It performs every join and every refusal and
writes nothing, so the number that decides whether the layer is worth
assembling -- how many indications survive the node filter -- is
readable before anything lands on disk.

WHY THIS LAYER ADDS NO NODES. The registration refuses it: a new
SmallMolecule would move every other population in the repository and
make `make auc-delta` uninterpretable for the release that introduces
the independent label. 111 of the 405 indicated compounds that have a
structure in ChEMBL are not in this graph and are dropped by this rule,
counted as `compound_not_in_graph`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.emit import Emit  # noqa: E402
from kgav.indications import (  # noqa: E402
    EVIDENCE_TIER,
    INDICATION_SOURCE,
    IndicationError,
    disease_xrefs,
    ingest_indications,
)
from kgav.pipeline import CHEMBL_DB, layer_dir, peer_nodes  # noqa: E402
from kgav.provenance import SOURCE_STAMP_NAME, write_source  # noqa: E402
from kgav.schema import load_schema, report_violations  # noqa: E402
from kgav import virus_register as vr  # noqa: E402

SELF_LAYER = "indications"


def graph_ids(all_nodes) -> tuple[set[str], set[str]]:
    """(SmallMolecule ids, Disease ids) from the layers built before this.

    Both endpoints have to already exist: the compounds because the
    registration forbids adding any, and the diseases because
    ingest_spine builds them from the register and this layer must not
    invent a Disease node for a MONDO id the spine declined to create.
    """
    compounds: set[str] = set()
    diseases: set[str] = set()
    for n in all_nodes:
        if n.get("class") == "SmallMolecule":
            compounds.add(n["id"])
        elif n.get("class") == "Disease":
            diseases.add(n["id"])
    return compounds, diseases


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--db", type=Path, default=ROOT / CHEMBL_DB,
                    help="ChEMBL SQLite (default: the one the pipeline names)")
    ap.add_argument("--config", type=Path,
                    default=ROOT / "config" / "viruses.yaml")
    ap.add_argument("--mondo", type=Path,
                    default=ROOT / "data" / "raw" / "mondo" / "mondo.json")
    # THE LAYER NAME IS WRITTEN OUT, not passed as SELF_LAYER, and the
    # duplication is deliberate: assemble_release reads
    # `v0.1-indications`, and tests/test_assemble.py::
    # test_every_ingest_script_writes_where_assemble_reads checks that
    # this default NAMES that directory. A constant here reads as
    # tidier and makes the guard unable to see the name, which is how a
    # driver comes to write where nothing looks.
    ap.add_argument("--out", type=Path,
                    default=layer_dir(ROOT, "indications"),
                    help="layer directory (default: the pipeline's)")
    ap.add_argument("--inspect", action="store_true",
                    help="report every join and refusal; write nothing")
    args = ap.parse_args(argv)

    if not args.db.exists():
        print(f"no ChEMBL database at {args.db}")
        return 1
    if not args.mondo.exists():
        print(f"no MONDO release at {args.mondo}")
        return 1

    viruses = vr.load(args.config)["viruses"]
    try:
        xrefs = disease_xrefs(viruses, args.mondo)
    except IndicationError as exc:
        print(f"REFUSED: {exc}")
        return 1

    compounds, diseases = graph_ids(peer_nodes(ROOT, SELF_LAYER))
    print(f"peers hold {len(compounds):,} SmallMolecule and "
          f"{len(diseases)} Disease node(s)")
    print(f"\n{len(xrefs)} registered disease(s), xrefs derived from "
          f"{args.mondo.name}:")
    for mondo in sorted(xrefs):
        mesh, efo = xrefs[mondo]
        mark = "" if mondo in diseases else "   NOT IN GRAPH"
        print(f"  {mondo:16s} mesh={sorted(mesh) or '-'} "
              f"efo={sorted(efo)}{mark}")

    em = Emit()
    s = ingest_indications(args.db, xrefs, compounds, diseases, em)

    print(f"\n{s['rows']:,} indication row(s) over the registered diseases")
    print(f"  {s['treats']:7,} TREATS        (max_phase_for_ind 4)")
    print(f"  {s['in_trial']:7,} IN_TRIAL_FOR  (0.5 through 3)")
    print(f"  {s['edges']:7,} edges")
    print("\n  per disease:")
    for mondo in sorted(xrefs):
        t, i = s.get(f"treats:{mondo}", 0), s.get(f"in_trial:{mondo}", 0)
        if t or i:
            print(f"    {mondo:16s} TREATS {t:4,}   IN_TRIAL_FOR {i:4,}")
    print("\n  dropped:")
    for key, label in (
            ("compound_not_in_graph",
             "indicated compound is not a node here -- the registration "
             "forbids adding one"),
            ("in_trial_without_nct",
             "no ClinicalTrials reference, and nct_id is a REQUIRED "
             "qualifier"),
            ("row_without_structure", "no structure in ChEMBL"),
            ("row_without_phase", "no max_phase_for_ind"),
            ("disease_not_in_graph", "no Disease node for the MONDO id")):
        if s[key]:
            print(f"    {key:<24} {s[key]:7,}  {label}")
    for key, label in (("pair_with_several_trials",
                        "pair(s) cite more than one trial; support_count "
                        "carries the number"),
                       ("pair_with_several_phases",
                        "pair(s) are filed at several phases; the furthest "
                        "is the edge")):
        if s[key]:
            print(f"\n  {s[key]:7,} {label}")

    if args.inspect:
        print("\n--inspect: nothing written.")
        return 0
    if not s["edges"]:
        print("\nNo indication survived the joins. Nothing written: an empty")
        print("layer directory would assemble silently and look like a label.")
        return 1

    # VALIDATED HERE, NOT ONLY AT ASSEMBLY. The assembler validates the
    # whole graph, so a bad edge from this layer surfaces as a schema
    # failure over 396,000 edges with no indication of which layer wrote
    # it. IN_TRIAL_FOR requires nct_id AND status, and this layer's
    # entire reason for dropping a row is that requirement, so the check
    # that the requirement is met belongs next to the drop rule.
    #
    # The node scope is the PEER nodes, because this layer emits none
    # and both of its endpoints live elsewhere -- without them every
    # edge would read as DANGLING, which would be a statement about the
    # validation scope and not about the edges.
    peers = peer_nodes(ROOT, SELF_LAYER)
    violations = load_schema().validate_batch(peers, em.edges)
    if report_violations(violations):
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("")
    (args.out / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in em.edges))
    write_source(args.out, chembl_indications={
        "table": "drug_indication",
        "refs_table": "indication_refs",
        "ref_type": "ClinicalTrials",
        "db": str(args.db.relative_to(ROOT)),
        "evidence_tier": EVIDENCE_TIER,
    })
    print(f"\nwrote {args.out}  (source {INDICATION_SOURCE})")
    print(f"recorded the table and the ref type in {SOURCE_STAMP_NAME}")
    print("NO NODES. nodes.jsonl is written EMPTY rather than omitted, so a")
    print("layer that adds no node is distinguishable from one that failed")
    print("before writing.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
