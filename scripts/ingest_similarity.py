#!/usr/bin/env python3
"""Build CHEMICALLY_SIMILAR_TO edges so M6 (cold start) has something to walk.

Reports the similarity distribution before committing to a threshold, because
0.7 is conventional rather than justified.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.emit import Emit
from kgav.schema import load_schema
from kgav.similarity import (
    build,
    fingerprints,
    load_compounds,
    similarity_distribution,
)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--out", type=Path, default=root / "data/releases/v0.1-similarity")
    ap.add_argument("--threshold", type=float, default=0.7)
    ap.add_argument("--top-k", type=int, default=10)
    ap.add_argument("--survey-only", action="store_true")
    args = ap.parse_args()

    # The schema pins tanimoto_ecfp4 at min 0.7. Passing a lower threshold here
    # would emit edges the validator rejects, so fail loudly instead of
    # producing an invalid release -- and changing the cut means changing the
    # schema, which forces a version bump and a record of the decision.
    schema = load_schema()
    for ec in schema.edge_classes:
        if ec.predicate == "CHEMICALLY_SIMILAR_TO":
            floor = (ec.qualifiers.get("tanimoto_ecfp4") or {}).get("min")
            if floor is not None and args.threshold < floor:
                raise SystemExit(
                    f"--threshold {args.threshold} is below the schema minimum "
                    f"{floor} for tanimoto_ecfp4. Change the schema (with a "
                    f"version bump) if the cut should move.")

    smiles = load_compounds(args.release)
    print(f"{len(smiles):,} compounds in the release")

    # Reference set = compounds that already have target annotations, since M6
    # exists to borrow those. Query set = compounds that have none.
    #
    # MEASURED_INACTIVE_AGAINST COUNTS AS AN ANNOTATION. A compound assayed
    # against a viral protein and found not to inhibit it is not a cold-start
    # compound -- it has been measured, and the measurement was negative.
    #
    # It was excluded here until schema 0.11.0 made no difference, because
    # those rows were emitted as INHIBITS and landed in `annotated` anyway.
    # ce133f1 corrected the predicate and, without anyone tracing it, moved
    # 1,625 pairs into the QUERY set. M7 walks
    # similar -> INHIBITS -> protein -> gene -> virus, so those compounds
    # acquired similarity edges that let M7 propose activity against the exact
    # protein the graph records them as measured inactive against -- a
    # contradiction sourced from a measurement, which nothing checks.
    #
    # THIS IS THE COARSE FIX, deliberately. It removes a compound from every
    # cold-start route because it was measured inactive against ONE protein,
    # so a compound inactive on nsp5 also loses its host-directed M6 route.
    # The precise fix is a per-(compound, target) exclusion at traversal, which
    # the metapath language has no way to express -- there are no negative
    # constraints. Erring toward not contradicting a measurement.
    ANNOTATING = ("TARGETS", "INHIBITS", "MEASURED_INACTIVE_AGAINST")
    annotated: set[str] = set()
    measured_inactive: set[str] = set()
    for line in (args.release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] in ANNOTATING:
            annotated.add(e["subject"])
            if e["predicate"] == "MEASURED_INACTIVE_AGAINST":
                measured_inactive.add(e["subject"])
    query = set(smiles) - annotated
    print(f"  {len(annotated):,} have target annotations (reference set)")
    print(f"  {len(query):,} have none (query set -- unreachable without M6)")
    if measured_inactive:
        print(f"  {len(measured_inactive):,} of the reference set are there "
              f"ONLY via a measured-inactive\n  result; they are excluded from "
              f"the query set so M7 cannot propose activity\n  against a "
              f"protein they were measured not to inhibit")

    q_fps, q_stats = fingerprints({i: smiles.get(i, "") for i in query})
    r_fps, r_stats = fingerprints({i: smiles.get(i, "") for i in annotated})
    print(f"\nfingerprints: {q_stats['fingerprinted']:,} query, "
          f"{r_stats['fingerprinted']:,} reference")
    bad = q_stats["unparseable"] + r_stats["unparseable"]
    missing = q_stats["no_smiles"] + r_stats["no_smiles"]
    if bad or missing:
        print(f"  {bad:,} unparseable, {missing:,} without SMILES "
              f"(unreachable by any method)")

    dist = similarity_distribution(q_fps, r_fps)
    print(f"\nquery compounds with at least one neighbour at each cut "
          f"(of {len(q_fps):,}):")
    for cut in sorted(dist, reverse=True):
        print(f"  >= {cut:.2f}   {dist[cut]:>6,}  ({100*dist[cut]/max(len(q_fps),1):>5.1f}%)")
    if args.survey_only:
        return 0

    em = Emit()
    s = build(em, smiles, query, annotated, "infores:kgav-computed",
              args.threshold, args.top_k)
    print(f"\n{s['edges']:,} CHEMICALLY_SIMILAR_TO edges at tanimoto >= "
          f"{args.threshold} (top {args.top_k} per query)")
    print(f"  {s['queries_connected']:,} previously unreachable compounds now "
          f"have a route")

    # Endpoints all exist in the release, so validate against its nodes.
    nodes = [json.loads(x) for x in (args.release / "nodes.jsonl").read_text().splitlines()
             if x.strip()]
    violations = load_schema().validate_batch(nodes, em.edges)
    if violations:
        print(f"\nSCHEMA FAIL — {len(violations)} violations")
        for v in violations[:8]:
            print(f"  {v}")
        return 1
    print("\nschema: PASS")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("")
    (args.out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in em.edges))
    print(f"{len(em.edges):,} edges -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
