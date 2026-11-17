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
    ASSAYS,
    NCATS_UNITS,
    COLUMN_CANDIDATES,
    MIN_EFFICACY_PCT,
    SOURCE,
    PolarityError,
    check_polarity,
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
from kgav.schema import load_schema, report_violations
from kgav.pipeline import layer_dir

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
    ap.add_argument("--out", type=Path, default=layer_dir(root, "ncats"))
    # DEFAULT None, NOT DEFAULT_DATE. Each assay declares its own
    # first_asserted_date from its own cited paper, and a flag defaulting to
    # the CPE date would have silently stamped a 2020-08-18 preprint date on
    # a 2022 SLAS assay and on two MERS/SARS-CoV screens -- overriding the
    # per-assay value the registry exists to hold. An explicit --date still
    # wins, for all assays, which is what an override is for.
    ap.add_argument("--date", default=None,
                    help="override first_asserted_date for EVERY selected "
                         "assay; by default each uses its own cited paper's "
                         "date")
    ap.add_argument("--ac50-units", default=None, choices=["M", "uM", "nM", "logM"])
    ap.add_argument("--min-efficacy-pct", type=float, default=MIN_EFFICACY_PCT,
                    help=f"minimum POSITIVE efficacy for an active (default "
                         f"{MIN_EFFICACY_PCT}; Chen et al. selected hits at 55)")
    # THE CONCENTRATION UNIT, DECLARED WITH ITS EVIDENCE (fix-list P7).
    # cpe.tsv names it -- CONC0(M) -- and the two files this layer reads,
    # cpe_assay.csv and cpe_tox.csv, use a bare `conc0`. Measured across
    # them: cpe_assay.csv's conc3 runs 4e-08 to 2e-04 with median 2e-05 and
    # cpe.tsv's CONC3(M) runs 4e-08 to 6e-05 with the same median, so the
    # unlabelled columns are molar by evidence. A header that states its
    # own unit overrides this.
    ap.add_argument("--conc-units", default="M", choices=sorted(NCATS_UNITS),
                    help="unit of the conc* columns WHEN THE HEADER DOES NOT "
                         "STATE ONE. cpe.tsv says CONC0(M); cpe_assay.csv "
                         "and cpe_tox.csv do not, and their value "
                         "distributions match cpe.tsv's, which is the "
                         "evidence for the default.")
    ap.add_argument("--assay", default=None,
                    help="assay_name to keep WITHIN the selected assays; "
                         "'all' disables the filter. Normally unnecessary -- "
                         "each assay in --assays declares its own name.")
    # WHICH ASSAYS, FROM THE REGISTRY RATHER THAN FROM PATHS. Before
    # 2026-11-17 this script took one file pair and one hardcoded taxon, so
    # adding an assay meant another run with another --cpe and remembering
    # to change the taxon with it -- three chances to file a MERS screen
    # against SARS-CoV-2. Now the assay is named and everything about it
    # (virus, polarity, cell line, counterscreen, citation, filenames)
    # comes from kgav.ncats.ASSAYS.
    ap.add_argument("--assays", default="all",
                    help=f"comma-separated assay keys, or 'all' (default). "
                         f"Available: {', '.join(ASSAYS)}")
    ap.add_argument("--raw", type=Path, default=root / "data/raw/ncats",
                    help="directory holding the assay exports")
    ap.add_argument("--inspect", action="store_true")
    args = ap.parse_args()

    keys = list(ASSAYS) if args.assays == "all" else \
        [k.strip() for k in args.assays.split(",") if k.strip()]
    unknown = [k for k in keys if k not in ASSAYS]
    if unknown:
        print(f"unknown assay(s): {', '.join(unknown)}\n"
              f"available: {', '.join(ASSAYS)}")
        return 2
    specs = [ASSAYS[k] for k in keys]

    # --cpe/--tox stay as explicit overrides for the CPE pair, so every
    # existing invocation keeps working; the other assays resolve from the
    # registry under --raw.
    paths = {}
    for s in specs:
        if s.key == "cpe":
            paths[s.key] = (args.cpe, args.tox)
        else:
            paths[s.key] = (args.raw / s.file, args.raw / s.counterscreen_file)

    missing = [(s, a) for s in specs for a in (paths[s.key][0],)
               if not a.exists()]
    if missing:
        print("not found:\n" + "\n".join(f"  {a}" for _s, a in missing))
        print(f"\nDownload them with:\n  mkdir -p {args.raw}\n  base={BASE}")
        for s in specs:
            print(f'  curl -L -o {paths[s.key][0]} \\\n'
                  f'    "$base/{s.remote}"')
            print(f'  curl -L -o {paths[s.key][1]} \\\n'
                  f'    "$base/{s.remote_counterscreen}"')
        return 1

    # ---------------------------------------------------- read every assay
    # ONE PASS PER ASSAY, and the polarity of each is CHECKED against its
    # own rows before anything is counted. A file that contradicts its
    # declaration stops the run: applying the CPE sign rule to a
    # loss-of-signal export yields zero actives, which looks exactly like a
    # screen that found nothing.
    loaded = []
    for s in specs:
        apath, tpath = paths[s.key]
        try:
            rows, cols = read_table(apath)
            tox_rows, tox_cols = ([], None)
            if tpath.exists():
                tox_rows, tox_cols = read_table(tpath)
        except ColumnError as e:
            print(f"HEADER NOT UNDERSTOOD in {apath.name}\n\n{e}")
            return 1

        in_assay = [r for r in rows
                    if "assay" not in cols
                    or _key(r.get(cols["assay"], "")) == _key(s.assay_name)]
        try:
            neg, pos = check_polarity(in_assay, cols, s, args.min_efficacy_pct)
        except PolarityError as e:
            print(f"POLARITY REFUSED in {apath.name}\n\n{e}")
            return 1

        ac50_units = args.ac50_units or (
            unit_for(cols["ac50"]) if "ac50" in cols else "uM")
        tox_units = unit_for(tox_cols["ac50"]) \
            if tox_cols and "ac50" in tox_cols else "uM"
        tox_index = index_counterscreen(tox_rows, tox_cols,
                                        s.counterscreen_assay) if tox_cols else {}

        _mapping_report(s.key, apath, cols, rows)
        if tox_cols:
            _mapping_report(f"{s.key} counterscreen", tpath, tox_cols, tox_rows)
        else:
            print(f"\n  ! {tpath} not found -- no selectivity index, and for "
                  f"a loss-of-signal\n    assay no cytotoxicity refusal "
                  f"either, which is the confound itself.")

        print(f"\n  {s.key}: {s.taxon}, "
              f"{'LOSS' if s.polarity < 0 else 'GAIN'}-of-signal, "
              f"{s.assay_type}")
        print(f"    polarity verified: {neg:,} responders at or beyond "
              f"-{args.min_efficacy_pct:g}% vs {pos:,} at or beyond "
              f"+{args.min_efficacy_pct:g}%")
        print(f"    {len(in_assay):,} of {len(rows):,} rows in this assay; "
              f"AC50 as {ac50_units}, counterscreen as {tox_units}; "
              f"{len(tox_index):,} counterscreen rows indexed")
        print(f"    cites {s.publication}, dated {args.date or s.date}")

        tops = [x for r in in_assay
                if (x := top_concentration_nm(r, cols, args.conc_units))
                is not None]
        if tops:
            print(f"    top concentration: {min(tops):,.0f} .. "
                  f"{max(tops):,.0f} nM")
        else:
            print("    ! NO CONCENTRATION COLUMNS PARSED -- every inactive "
                  "will be refused")

        counts: Counter = Counter()
        for r in in_assay:
            sid = (r.get(cols["sample"]) or "").strip() if "sample" in cols else ""
            v, _q, why = classify_row(
                r, cols, ac50_units=ac50_units, conc_units=args.conc_units,
                tox_row=tox_index.get(sid), tox_cols=tox_cols,
                tox_ac50_units=tox_units,
                min_efficacy_pct=args.min_efficacy_pct, spec=s)
            counts[f"{v or 'none'}: {why}"] += 1
        print("    call breakdown (before structure normalization)")
        for k, n in counts.most_common():
            print(f"      {k:<44} {n:>7,}")

        loaded.append((s, rows, cols, tox_index, tox_cols,
                       ac50_units, tox_units))

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
    stats: Counter = Counter()
    # ONE Emit FOR ALL ASSAYS, because a compound measured in two of them is
    # ONE node carrying two label edges -- the assays answer different
    # questions about it and `assay_type` keeps them apart. Per-assay totals
    # are printed too, since a single sum would hide an assay that
    # contributed nothing.
    for s, rows, cols, tox_index, tox_cols, ac50_units, tox_units in loaded:
        print(f"\nnormalizing and ingesting {s.key} (rdkit)")
        got = ingest_cpe(em, rows, cols, normalize_chemical,
                         tox_by_sample=tox_index, tox_cols=tox_cols,
                         existing_compounds=existing,
                         source=SOURCE, date=args.date,
                         ac50_units=ac50_units, tox_ac50_units=tox_units,
                         conc_units=args.conc_units,
                         assay_filter=(None if args.assay == "all"
                                       else args.assay or s.assay_name),
                         min_efficacy_pct=args.min_efficacy_pct,
                         spec=s)
        print(f"  {s.key:<14} {got['active']:>6,} active | "
              f"{got['inactive']:>7,} inactive | "
              f"{got['reason_weak_hit_above_active_threshold']:>5,} weak hits "
              f"refused")
        for k, v in got.items():
            stats[k] += v

    print(f"\n{stats['active']:,} active | {stats['inactive']:,} measured INACTIVE")
    print(f"  {stats['with_same_plate_si']:,} actives carry a same-plate "
          f"selectivity index")
    print(f"  {stats['no_cytotoxicity_detected']:,} showed no cytotoxicity up to "
          f"the top concentration")
    print("\n  per-row outcome:")
    for k in sorted(k for k in stats if k.startswith("reason_")):
        print(f"    {k.replace('reason_', ''):<38} {stats[k]:>7,}")
    ORDERED = ("other_assay", "no_structure", "unnormalizable",
               "duplicate_structure", "upgraded_inactive_to_active",
               "paired_with_counterscreen", "joined_existing_compound",
               "new_compound_node")
    for k in ORDERED:
        if stats[k]:
            print(f"    {k:<38} {stats[k]:>7,}")

    # A HARDCODED LIST OF KEYS MAKES A NEW COUNTER INVISIBLE, which is how
    # fix-list L19's `tautomer_conflated_nodes` was added to ingest_screens,
    # populated, and printed nowhere -- the first run after wiring it looked
    # exactly like a run that had found none. The curated list above stays,
    # because its ORDER is the narrative of the filter; anything else the
    # ingest counted is printed after it rather than dropped.
    shown = set(ORDERED) | {"rows", "in_assay", "active", "inactive", "none"}
    extra = sorted(k for k, v in stats.items()
                   if v and k not in shown and not k.startswith("reason_"))
    for k in extra:
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
    if report_violations(violations, note=(
            "\n  No peer layer was read, so the node set is this layer's own "
            "nodes\n  only -- which omits the OrganismTaxon the spine owns. "
            "DANGLING here\n  means the validation scope is wrong, not the "
            "edges." if not release_nodes else "")):
        return 1

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
