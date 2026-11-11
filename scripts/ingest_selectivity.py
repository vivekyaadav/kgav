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
from kgav.provenance import chembl_release, write_source
from kgav.schema import load_schema, report_violations
from kgav.selectivity import SI_THRESHOLD, ingest_selectivity, query_pairs
from kgav.virus_register import declared_taxa
from kgav.pipeline import layer_dir, peer_nodes

#: This layer's own name, so peer_nodes can exclude it.
SELF_LAYER = "selectivity"


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path,
                    default=root / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--config", type=Path, default=root / "config" / "viruses.yaml")
    # --peers, NOT --release. This layer no longer reads the assembled
    # graph; it reads the layer directories built before it (L20 (b)).
    # Same flag name and same default as ingest_ncats, which has read
    # its peers this way since 2026-10-09.
    ap.add_argument("--peers", type=Path, default=root / "data/releases",
                    help="directory holding the per-layer releases")
    ap.add_argument("--out", type=Path, default=layer_dir(root, "selectivity"))
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
                     # THE GATE WAS IN THE ACCOUNTING AND NOT IN THE REPORT.
                     # The loop below sums every av_dropped_* key, so the
                     # total closed from the day this gate was added while
                     # the breakdown a reader sees did not mention it --
                     # which is E7's defect with the counter present.
                     ("av_dropped_not_cell_based",
                      "dropped: not a cell-based antiviral readout "
                      "(biochemical target_type)"),
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
    # HOW THE RATIO WAS PAIRED, per pair, because "same document" is not
    # "same cell system" and the difference is the whole meaning of the
    # number. selectivity.query_pairs says what each one is.
    for k, label in (("ec50_paired_on_the_cc50s_cell_line",
                      "EC50 taken from the CC50's own cell line"),
                     ("ec50_not_pairable_on_a_cell_line",
                      "no EC50 in the CC50's cell line (or either unnamed): "
                      "the conservative row is used"),
                     ("ec50_conservative_over_several_rows",
                      "several candidate EC50s; the largest is used")):
        if q[k]:
            print(f"    {q[k]:>8,}  {label}")
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

    # THE COMPOUNDS THE RELEASE DOES NOT HOLD ARE DROPPED BEFORE THE EMIT,
    # NOT AFTER IT, so every counter below describes the layer that ships.
    # Until 2026-11-11 the dangling edges were removed at the END and the
    # breakdown printed above them, so the run reported "2,758 emitted edges:
    # 1,694 agree, 253 ..., 811 ..." for a layer holding 2,638 edges with
    # 1,609 verified -- the three classes summed to the PRE-drop set and were
    # read as describing the file. The emitted edge set is unchanged by moving
    # the filter: `dangling` tested subject-or-object, the objects are taxa
    # and always present, so compound presence is the whole of it. Verified
    # byte-identical against the post-filter run before this landed.
    #
    # AND THEY ARE NOT AN L27 SYMPTOM, which was the first thing to rule out:
    # this layer keys compounds on ChEMBL's `standard_inchi_key` exactly as
    # chembl.py did before the identity repair. Measured on this release --
    # of the 121 paired compound ids absent from the release, the identity
    # policy maps ZERO of them onto a node that is present. They are
    # compounds the chemistry layer filtered out, which is what the comment
    # below always said.
    nodes = peer_nodes(root, SELF_LAYER, releases=args.peers)
    known = {n["id"] for n in nodes}
    before = len(rows)
    rows = [r for r in rows
            if not r.get("inchikey") or f"INCHIKEY:{r['inchikey']}" in known]
    if before != len(rows):
        print(f"\n  {before - len(rows):,} paired row(s) name a compound the "
              f"release does not hold and are dropped before the emit")

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
    # WHAT selectivity_verified MEANS ON THIS RUN, which until now was
    # computed, written onto every edge and never reported.
    for k, label in (("context_matched", "EC50 and CC50 name the SAME cell line"),
                     ("context_cell_line_mismatch",
                      "they name DIFFERENT cell lines (not verifiable)"),
                     ("context_cell_line_unnamed",
                      "at least one side names no cell line (not verifiable)")):
        if s[k]:
            print(f"  {s[k]:>7,}  {label}")

    per_virus: Counter = Counter()
    for e in em.edges:
        per_virus[e["object"]] += 1
    if per_virus:
        print("\nper virus:")
        for tid, n in per_virus.most_common():
            print(f"  {tid:<22} {n:>6,}")

    # READS THE LAYERS BEFORE IT, not the assembled release (fix-list
    # L20, option (b)). pipeline.peer_nodes dedupes by id in
    # LAYER_PRECEDENCE order, which is assemble_release's own rule, and
    # on v0.1 reproduces the assembled node set exactly -- 78,789 both
    # ways, zero difference. Nothing needs an assembly to exist before
    # this layer runs, so the pipeline assembles once, at the end.
    dangling = [e for e in em.edges
                if e["subject"] not in known or e["object"] not in known]
    # KEPT AS A GUARD RATHER THAN AS THE FILTER. The rows were filtered
    # above, so this must now find nothing; if it finds something, a
    # dangling endpoint is arriving by some route the row filter does not
    # see and the run REFUSES rather than quietly trimming.
    if dangling:
        print(f"\n  ! {len(dangling):,} emitted edge(s) still dangle after "
              f"the row filter, e.g. {dangling[0]['subject']} -> "
              f"{dangling[0]['object']}. An endpoint is reaching the emit by "
              f"a route the filter does not cover.")
        return 1

    violations = load_schema().validate_batch(nodes, em.edges)
    if report_violations(violations):
        return 1

    # DATE PRECISION, REPORTED (fix-list C10). A year-only date is written
    # as a year; the placeholder is counted apart, because "undated" spread
    # across seven sources is seven gaps with seven remedies.
    if s["dated_year_only"] or s["undated_placeholder"]:
        tot = s["dated_year_only"] + s["undated_placeholder"]
        print(f"\n  dates: {s['dated_year_only']:,} at YEAR precision, "
              f"{s['undated_placeholder']:,} undated "
              f"({s['undated_placeholder'] / tot:.1%})")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("")
    (args.out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in em.edges))
    # WHICH RELEASE OF THE SOURCE (fix-list V2/V3/V5). Read from the
    # file on disk, never from a register kept by hand: a release named in
    # a path or a filename is evidence, one typed into a config is a claim.
    write_source(args.out, chembl=chembl_release(args.db))
    print(f"{len(em.edges):,} edges -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
