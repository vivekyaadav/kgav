#!/usr/bin/env python3
"""Ingest the NCATS OpenData SARS-CoV-2 CPE screen and its counterscreen.

    mkdir -p data/raw/ncats
    base=https://opendata.ncats.nih.gov/public/odp/assay
    curl -L -o data/raw/ncats/cpe_assay.csv \\
      "$base/SARS-CoV-2_Cytopathic_Effect_(CPE).csv"
    curl -L -o data/raw/ncats/cpe_tox.csv \\
      "$base/SARS-CoV-2_Cytopathic_Effect_(Host_Tox_Counterscreen).csv"

    python scripts/ingest_ncats.py --inspect      # mapping and call breakdown
    python scripts/ingest_ncats.py                # build the layer

NOT the portal's root cpe.tsv: it pools three assays run on different
protocols and carries no cytotoxicity column.

--inspect first, always. The mapping and the activity call are the two things
a schema check cannot verify: a mis-mapped potency column, or a direction read
the wrong way in a gain-of-signal assay, produces a release that validates,
counts plausibly and encodes the opposite conclusion.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.assemble import LAYER_PRECEDENCE
from kgav.emit import Emit
from kgav.ncats import (
    CANONICAL_ASSAY,
    COLUMN_CANDIDATES,
    COUNTERSCREEN_ASSAY,
    CPE_PUBLICATION,
    DEFAULT_DATE,
    MIN_EFFICACY_PCT,
    SARS_COV_2,
    SOURCE,
    ColumnError,
    _key,
    classify_row,
    index_counterscreen,
    ingest_cpe,
    join_report,
    read_table,
    top_concentration_nm,
    unit_for,
)
from kgav.schema import load_schema

# Below this, assume the normalization is wrong rather than the library being
# unusual. These are approved and annotated compound collections; ChEMBL's
# coronavirus set and DrugCentral overlap them heavily, so a few per cent join
# rate means the keys are being built differently, not that the overlap is
# genuinely tiny.
MIN_PLAUSIBLE_JOIN_RATE = 0.05

BASE = "https://opendata.ncats.nih.gov/public/odp/assay"

# This layer's own name in LAYER_PRECEDENCE, excluded from the peer
# set so ownership never depends on whether this script already ran.
SELF_LAYER = "ncats"


def _mapping_report(name: str, path: Path, cols: dict, rows: list[dict]) -> None:
    print(f"\n{name}: {path}  ({len(rows):,} rows)")
    for field in COLUMN_CANDIDATES:
        got = cols.get(field)
        print(f"  {' ' if got else '!'} {field:<13} <- {got or '(not found)'}")
    print(f"    {'conc_cols':<13} <- {len(cols.get('conc_cols', []))} columns")
    if "assay" in cols:
        c = Counter(r.get(cols["assay"], "") for r in rows)
        print(f"  assays present ({len(c)}):")
        for k, n in c.most_common(6):
            print(f"      {k!r:<52} {n:>6,}")
    if "library" in cols:
        c = Counter(r.get(cols["library"], "") for r in rows)
        print(f"  libraries present ({len(c)}), top 5:")
        for k, n in c.most_common(5):
            print(f"      {k!r:<52} {n:>6,}")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpe", type=Path, default=root / "data/raw/ncats/cpe_assay.csv")
    ap.add_argument("--tox", type=Path, default=root / "data/raw/ncats/cpe_tox.csv",
                    help="host-tox counterscreen; without it no selectivity "
                         "index is computed and actives carry EC50 only")
    ap.add_argument("--peers", type=Path, default=root / "data/releases",
                    help="directory holding the per-layer releases. Node "
                         "ownership is decided against the OTHER layers, "
                         "never against the assembled graph: this layer "
                         "contributes to that graph, so reading it back makes "
                         "the decision depend on whether this script has "
                         "already run.")
    ap.add_argument("--out", type=Path, default=root / "data/releases/v0.1-ncats")
    ap.add_argument("--date", default=DEFAULT_DATE)
    ap.add_argument("--ac50-units", default=None, choices=["M", "uM", "nM", "logM"])
    ap.add_argument("--min-efficacy-pct", type=float, default=MIN_EFFICACY_PCT,
                    help=f"minimum POSITIVE efficacy for an active (default "
                         f"{MIN_EFFICACY_PCT}; Chen et al. selected hits at 55)")
    ap.add_argument("--assay", default=CANONICAL_ASSAY,
                    help="assay_name to keep; 'all' disables the filter")
    ap.add_argument("--inspect", action="store_true")
    args = ap.parse_args()

    if not args.cpe.exists():
        print(f"not found: {args.cpe}\n\nDownload both files with:\n"
              f"  mkdir -p {args.cpe.parent}\n"
              f'  base={BASE}\n'
              f'  curl -L -o {args.cpe} \\\n'
              f'    "$base/SARS-CoV-2_Cytopathic_Effect_(CPE).csv"\n'
              f'  curl -L -o {args.tox} \\\n'
              f'    "$base/SARS-CoV-2_Cytopathic_Effect_(Host_Tox_Counterscreen).csv"')
        return 1

    try:
        rows, cols = read_table(args.cpe)
        tox_rows, tox_cols = ([], None)
        if args.tox.exists():
            tox_rows, tox_cols = read_table(args.tox)
    except ColumnError as e:
        print(f"HEADER NOT UNDERSTOOD\n\n{e}")
        return 1

    ac50_units = args.ac50_units or (unit_for(cols["ac50"]) if "ac50" in cols else "uM")
    tox_units = unit_for(tox_cols["ac50"]) if tox_cols and "ac50" in tox_cols else "uM"
    assay_filter = None if args.assay == "all" else args.assay

    _mapping_report("CPE", args.cpe, cols, rows)
    if tox_cols:
        _mapping_report("counterscreen", args.tox, tox_cols, tox_rows)
    else:
        print(f"\n  ! {args.tox} not found -- no selectivity index will be "
              f"computed.\n    Actives will carry EC50 with no CC50 beside it.")

    tox_index = index_counterscreen(tox_rows, tox_cols, COUNTERSCREEN_ASSAY) \
        if tox_cols else {}
    print(f"\nAC50 read as {ac50_units}, counterscreen AC50 as {tox_units}")
    print(f"assay filter: {assay_filter or '(none -- all assays)'}")
    print(f"active requires POSITIVE efficacy >= {args.min_efficacy_pct}%")
    print(f"counterscreen rows indexed by sample: {len(tox_index):,}")
    print(f"every edge cites {CPE_PUBLICATION}, dated {args.date}")

    # top concentration actually present, so the censoring point is reported
    # rather than assumed
    tops = [t for r in rows if (t := top_concentration_nm(r, cols)) is not None]
    if tops:
        print(f"top concentration across rows: min {min(tops):,.0f} nM, "
              f"max {max(tops):,.0f} nM")
    else:
        print("  ! NO CONCENTRATION COLUMNS PARSED — every inactive will be "
              "refused,\n    because its censoring point is unknown.")

    # ---------------------------------------------------- dry-run the calls
    counts: Counter = Counter()
    for r in rows:
        if assay_filter and "assay" in cols and \
                _key(r.get(cols["assay"], "")) != _key(assay_filter):
            counts["other_assay"] += 1
            continue
        sid = (r.get(cols["sample"]) or "").strip() if "sample" in cols else ""
        v, _, why = classify_row(r, cols, ac50_units=ac50_units,
                                 tox_row=tox_index.get(sid), tox_cols=tox_cols,
                                 tox_ac50_units=tox_units,
                                 min_efficacy_pct=args.min_efficacy_pct)
        counts[f"{v or 'none'}: {why}"] += 1
    print("\ncall breakdown (before structure normalization)")
    for k, n in counts.most_common():
        print(f"  {k:<46} {n:>7,}")

    if args.inspect:
        print("\n--inspect: nothing written. Re-run without it to build.")
        return 0

    try:
        from kgav.normalize_chem import normalize_chemical
    except ImportError as e:
        print(f"\nrdkit is required to normalize structures: {e}")
        return 1

    # NODE OWNERSHIP IS DECIDED AGAINST THE PEER LAYERS, NOT THE ASSEMBLED
    # GRAPH. Reading the assembled graph was circular and self-destructive:
    # this layer contributes its new compounds to it, so a second run saw
    # those compounds "already present", declined to emit them, and the only
    # nodes its 5,254 edges had stopped existing. The assembler reported 5,254
    # DANGLING edges and the graph fell back to 66,007 nodes. Running an
    # ingest twice must be a no-op, and against peers it is.
    #
    # Full node dicts, not just ids, because this layer deliberately does not
    # re-emit a node it does not own (see ingest_cpe on is_approved) and so
    # cannot validate alone -- the OrganismTaxon belongs to the spine.
    # ingest_chembl.py validates against index.nodes + em.nodes for the same
    # reason.
    peer_dirs = [args.peers / f"v0.1-{name}" for name in LAYER_PRECEDENCE
                 if name != SELF_LAYER]
    peer_dirs = [d for d in peer_dirs if (d / "nodes.jsonl").exists()]
    existing: set[str] = set()
    # DEDUPED BY ID. Peer layers overlap on purpose -- a compound appears in
    # chembl, selectivity and hosttargets alike, and the assembler merges them
    # into one node. Concatenating the layer files hands validate_batch the
    # same id several times and it correctly reports DUPLICATE_NODE: 67,360
    # peer rows against 66,007 assembled nodes, so 1,353 violations that say
    # nothing about this layer. First writer wins, exactly as the assembler's
    # own precedence does.
    by_id: dict[str, dict] = {}
    for d in peer_dirs:
        for line in (d / "nodes.jsonl").read_text().splitlines():
            if line.strip():
                n = json.loads(line)
                by_id.setdefault(n["id"], n)
                if n["class"] == "SmallMolecule":
                    existing.add(n["id"])
    release_nodes = list(by_id.values())
    if peer_dirs:
        print(f"\n{len(existing):,} compounds across {len(peer_dirs)} peer "
              f"layers ({len(release_nodes):,} distinct nodes): "
              f"{', '.join(d.name.replace('v0.1-', '') for d in peer_dirs)}")
    else:
        print(f"\n  ! no peer layers found under {args.peers} -- every compound "
              f"treated as\n    new, so this layer will assert is_approved on "
              f"compounds another layer\n    knows better. Build the other "
              f"layers first.")

    em = Emit()
    print("normalizing and ingesting (rdkit, this takes a minute)")
    stats = ingest_cpe(em, rows, cols, normalize_chemical,
                       tox_by_sample=tox_index, tox_cols=tox_cols,
                       existing_compounds=existing,
                       taxon=SARS_COV_2, source=SOURCE, date=args.date,
                       ac50_units=ac50_units, tox_ac50_units=tox_units,
                       assay_filter=assay_filter,
                       min_efficacy_pct=args.min_efficacy_pct)

    print(f"\n{stats['active']:,} active | {stats['inactive']:,} measured INACTIVE")
    print(f"  {stats['with_same_plate_si']:,} actives carry a same-plate "
          f"selectivity index")
    print(f"  {stats['no_cytotoxicity_detected']:,} showed no cytotoxicity up to "
          f"the top concentration")
    print("\n  per-row outcome:")
    for k in sorted(k for k in stats if k.startswith("reason_")):
        print(f"    {k.replace('reason_', ''):<38} {stats[k]:>7,}")
    for k in ("other_assay", "no_structure", "unnormalizable",
              "duplicate_structure", "upgraded_inactive_to_active",
              "paired_with_counterscreen", "joined_existing_compound",
              "new_compound_node"):
        if stats[k]:
            print(f"    {k:<38} {stats[k]:>7,}")

    labelled = {e["subject"] for e in em.edges}
    jr = join_report(labelled, existing)
    print(f"\njoin against {len(peer_dirs)} peer layers ({len(existing):,} compounds)")
    print(f"  compounds labelled {jr['emitted']:>7,}")
    print(f"  already in graph   {jr['already_in_graph']:>7,}   {jr['join_rate']:.1%}")
    print(f"  new compounds      {jr['new_compounds']:>7,}")
    print("\n  A compound that joins carries the graph's existing TARGETS edges, "
          "so it can\n  be reached by M2-M6. One that does not is a label with "
          "no path to it.")

    if existing and jr["join_rate"] < MIN_PLAUSIBLE_JOIN_RATE:
        print(f"\nREFUSING TO WRITE — join rate {jr['join_rate']:.1%} is below "
              f"{MIN_PLAUSIBLE_JOIN_RATE:.0%}.\n"
              f"  These are approved and annotated collections, which ChEMBL and "
              f"DrugCentral\n  cover heavily. A rate this low means the keys are "
              f"being built differently\n  from the rest of the graph, not that "
              f"the overlap is small -- which is the\n  silent failure this layer "
              f"exists to avoid. Compare a few of these\n  InChIKeys against the "
              f"release before overriding.")
        return 1

    violations = load_schema().validate_batch(
        release_nodes + list(em.nodes.values()), em.edges)
    if violations:
        print(f"\nSCHEMA FAIL — {len(violations)} violations")
        for v in violations[:8]:
            print(f"  {v}")
        if not release_nodes:
            print(f"\n  No peer layer was read, so the node set is this "
                  f"layer's own nodes\n  only -- which omits the OrganismTaxon "
                  f"the spine owns. DANGLING here\n  means the validation "
                  f"scope is wrong, not the edges.")
        return 1
    print("\nschema: PASS")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in em.nodes.values()))
    (args.out / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in em.edges))
    print(f"{len(em.nodes):,} new compounds, {len(em.edges):,} edges -> {args.out}")
    print("\nNext: add this layer to the assembler's LAYER_PRECEDENCE, "
          "reassemble, then\n  python scripts/hard_negatives.py --selectivity "
          "verified-only \\\n      --selectivity-applies both "
          "--publication-disjoint")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
