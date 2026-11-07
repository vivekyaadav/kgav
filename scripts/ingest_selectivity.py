#!/usr/bin/env python3
"""Build the selectivity layer from same-document EC50/CC50 pairs."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.emit import Emit
from kgav.schema import load_schema, report_violations
from kgav.selectivity import SI_THRESHOLD, ingest_selectivity, query_pairs
from kgav.virus_register import declared_taxa


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path,
                    default=root / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--config", type=Path, default=root / "config" / "viruses.yaml")
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--out", type=Path, default=root / "data/releases/v0.1-selectivity")
    ap.add_argument("--threshold", type=float, default=SI_THRESHOLD)
    args = ap.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    # ONE MAP, FROM THE REGISTER. virus_register.declared_taxa exists for
    # exactly this and its docstring says so -- "the same shape the ingests
    # build ad hoc from isolate_taxa, derived here so chain_source cannot
    # drift out of the map" -- and nothing called it: four scripts built the
    # map by hand and the canonical function was reached only by its own test.
    # Five implementations of one mapping, which is how AUC 0.826, the
    # coverage floor and lift's trustworthy rule each ended up with two.
    #
    # declared_taxa also folds in chain_source.taxon, which the hand-built
    # copies omitted. Identical on the current register -- verified before
    # this change -- because validate()'s UNDECLARED_CHAIN_TAXON already
    # forces every differing chain taxon into isolate_taxa. So this is
    # behaviour-neutral today and stops being luck tomorrow.
    taxon_map = declared_taxa(cfg)

    print("pairing antiviral EC50/IC50 with CC50/TC50 within the same document")
    # EVERY GATE REPORTED. Four of them lived in SQL WHERE clauses and
    # discarded input with nothing recorded, which is fix-list item E7:
    # only the ChEMBL layer asserts `emitted + dropped == rows`.
    q: Counter = Counter()
    rows = query_pairs(args.db, sorted(taxon_map), stats=q)
    print(f"  {len(rows):,} matched (compound, document, virus) rows")
    print(f"\n  antiviral rows read {q['av_rows']:,} -> "
          f"{q['av_groups']:,} (compound, doc, virus) group(s)")
    for k, label in (("av_dropped_censored_or_null_relation",
                      "dropped: relation is not '=' (censored or absent)"),
                     ("av_dropped_unconvertible_units",
                      "dropped: units not convertible to nM"),
                     ("av_converted", "converted from a non-nM unit")):
        if q[k]:
            print(f"    {q[k]:>8,}  {label}")
    print(f"  toxicity rows read {q['tox_rows']:,} -> "
          f"{q['tox_groups']:,} group(s)")
    for k, label in (("tox_dropped_relation",
                      "dropped: relation is neither '=' nor '>'"),
                     ("tox_dropped_unconvertible_units",
                      "dropped: units not convertible to nM"),
                     ("tox_converted", "converted from a non-nM unit"),
                     ("tox_mixed_group_exact_preferred",
                      "group holds both exact and censored CC50s; the "
                      "exact minimum is used")):
        if q[k]:
            print(f"    {q[k]:>8,}  {label}")
    if q["no_matched_counterscreen"]:
        print(f"    {q['no_matched_counterscreen']:>8,}  antiviral group(s) "
              f"with no counterscreen in the same document")
    # THE ACCOUNTING CLOSES, or the run says so rather than being believed.
    for side in ("av", "tox"):
        dropped = sum(v for k, v in q.items()
                      if k.startswith(f"{side}_dropped_"))
        accounted = dropped + q[f"{side}_kept"]
        if accounted != q[f"{side}_rows"]:
            print(f"\n  ! ACCOUNTING DOES NOT CLOSE on the {side} side: "
                  f"{q[f'{side}_rows']:,} rows read, {accounted:,} accounted "
                  f"for ({dropped:,} dropped + {q[f'{side}_kept']:,} kept). "
                  f"A row is being discarded by a gate with no counter.")
            return 1
    print(f"\n  accounting closes: "
          f"{q['av_rows']:,} antiviral and {q['tox_rows']:,} toxicity rows "
          f"read, every one either dropped by a named gate or kept")
    # The units actually seen, so a new one shows up as a line rather than
    # as a silently smaller number.
    for k in sorted(k for k in q if k.startswith(("av_unit:", "tox_unit:"))):
        print(f"    {q[k]:>8,}  {k}")

    em = Emit()
    s = ingest_selectivity(em, rows, taxon_map, "infores:chembl", args.threshold)

    total = s["selective"] + s["cytotoxic"]
    print(f"\n{s['edges']:,} compounds with a verified selectivity index")
    if total:
        print(f"  {s['selective']:,} selective (SI >= {args.threshold:g})  "
              f"{100*s['selective']/total:.0f}%")
        print(f"  {s['cytotoxic']:,} cytotoxic (SI < {args.threshold:g})  "
              f"{100*s['cytotoxic']/total:.0f}%")
    for k in ("no_inchikey", "unmapped_taxon", "si_uncomputable"):
        if s[k]:
            print(f"  dropped {k}: {s[k]:,}")

    per_virus: Counter = Counter()
    for e in em.edges:
        per_virus[e["object"]] += 1
    if per_virus:
        print("\nper virus:")
        for tid, n in per_virus.most_common():
            print(f"  {tid:<22} {n:>6,}")

    nodes = [json.loads(x) for x in (args.release / "nodes.jsonl").read_text().splitlines()
             if x.strip()]
    known = {n["id"] for n in nodes}
    dangling = [e for e in em.edges if e["subject"] not in known or e["object"] not in known]
    if dangling:
        # A compound paired here but absent from the release means the
        # chemistry layer filtered it out; emitting the edge would dangle.
        print(f"\n  {len(dangling):,} edges reference compounds absent from the "
              f"release and are dropped")
        em.edges = [e for e in em.edges if e not in dangling]

    violations = load_schema().validate_batch(nodes, em.edges)
    if report_violations(violations):
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("")
    (args.out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in em.edges))
    print(f"{len(em.edges):,} edges -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
