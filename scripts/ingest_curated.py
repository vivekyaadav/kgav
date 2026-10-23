#!/usr/bin/env python3
"""Driver: host factors from published screens ORCS does not carry.

    mkdir -p data/raw/curated_screens
    # one TSV per paper, columns:
    #   ncbi_gene_id  taxon  direction  pmid  [gene_symbol cell_line
    #                                          screen_type effect_size date]
    python scripts/ingest_curated_screens.py --inspect   # validate, emit nothing
    python scripts/ingest_curated_screens.py             # build the layer

--inspect first, always. Every row that cannot become an edge is named with
its row number and its reason, and nothing is written. A table that half
parses is the failure this reports rather than ingests.

DIRECTION IS THE CURATOR'S, NOT THIS SCRIPT'S. src/kgav/curated.py refuses a
row that does not state one. The flavivirus CRISPR literature runs both ways
within one paper -- a resistance screen's enriched genes are dependency
factors, a susceptibility screen's are restriction factors -- and a default
here would silently file half of them backwards.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.curated import (  # noqa: E402
    CURATED_SOURCE,
    CuratedRowError,
    ingest_curated,
    read_table,
    resolve_symbols,
    validate_row,
)
from kgav.emit import Emit  # noqa: E402
from kgav.orcs import build_gene_to_proteins  # noqa: E402


def built_taxa(spine: Path) -> set[str]:
    """Taxon ids the spine actually built, bare."""
    out = set()
    p = spine / "nodes.jsonl"
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        if not line.strip():
            continue
        n = json.loads(line)
        if n.get("class") == "OrganismTaxon":
            out.add(str(n["id"]).removeprefix("NCBITaxon:"))
    return out


def gene_symbols(host: Path, spine: Path) -> dict[str, str]:
    """{NCBIGene CURIE: symbol} from the built Gene nodes."""
    out = {}
    for rel in (host, spine):
        p = rel / "nodes.jsonl"
        if not p.exists():
            continue
        for line in p.read_text().splitlines():
            if not line.strip():
                continue
            n = json.loads(line)
            if n.get("class") == "Gene":
                sym = (n.get("properties") or {}).get("symbol")
                if sym:
                    out[n["id"]] = sym
    return out


def resolve_mode(args) -> int:
    """Turn a list of gene symbols into a TSV skeleton with real ids.

    Prints the skeleton to stdout so it can be redirected into
    data/raw/curated_screens/<paper>.tsv and filled in.
    """
    symbols = [l for l in args.resolve.read_text().splitlines() if l.strip()]
    genes = gene_symbols(args.host, args.spine)
    if not genes:
        print(f"no Gene nodes under {args.host} or {args.spine}; build the "
              f"host layer first")
        return 2
    resolved, unresolved, ambiguous = resolve_symbols(symbols, genes)
    print(f"# {len(resolved)}/{len(symbols)} symbols resolved against "
          f"{len(genes):,} Gene nodes", file=sys.stderr)
    if ambiguous:
        print(f"# {len(ambiguous)} AMBIGUOUS, left out:", file=sys.stderr)
        for s, ids in sorted(ambiguous.items()):
            print(f"#   {s}: {ids}", file=sys.stderr)
    if unresolved:
        print(f"# {len(unresolved)} UNRESOLVED, left out. A hit whose gene is "
              f"not in the host layer can never become an edge; aliases are "
              f"not handled, so check these by hand:", file=sys.stderr)
        for s in unresolved:
            print(f"#   {s}", file=sys.stderr)
    print("ncbi_gene_id\ttaxon\tdirection\tpmid\tgene_symbol\tcell_line\t"
          "screen_type\teffect_size")
    for sym, gid in sorted(resolved.items(), key=lambda kv: kv[0]):
        print(f"{gid.removeprefix('NCBIGene:')}\t\t\t\t{sym}\t\t\t")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tables", type=Path,
                    default=ROOT / "data" / "raw" / "curated_screens")
    ap.add_argument("--host", type=Path,
                    default=ROOT / "data" / "releases" / "v0.1-host")
    ap.add_argument("--spine", type=Path,
                    default=ROOT / "data" / "releases" / "v0.1-spine")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "releases" / "v0.1-curated")
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--resolve", type=Path, default=None,
                    help="a file of gene SYMBOLS, one per line; prints a TSV "
                         "skeleton with their NCBIGene ids filled in")
    args = ap.parse_args(argv)

    if args.resolve:
        return resolve_mode(args)

    tables = sorted(args.tables.glob("*.tsv")) if args.tables.exists() else []
    if not tables:
        print(f"no .tsv tables in {args.tables}")
        print("Nothing to ingest. This is not an error: the layer is empty")
        print("until a screen is curated into it.")
        return 1

    print(f"{len(tables)} table(s) in {args.tables}\n")
    rows, bad = [], 0
    for path in tables:
        try:
            raw = read_table(path)
        except CuratedRowError as e:
            print(f"  {path.name}: REFUSED — {e}")
            bad += 1
            continue
        ok = 0
        for i, r in enumerate(raw, 2):          # 1 is the header
            try:
                rows.append(validate_row(r, i))
                ok += 1
            except CuratedRowError as e:
                print(f"  {path.name}: {e}")
                bad += 1
        print(f"  {path.name}: {ok}/{len(raw)} rows usable")

    print(f"\n{len(rows):,} validated rows, {bad} refused")
    by_taxon: dict[str, int] = {}
    by_pmid: dict[str, int] = {}
    for r in rows:
        by_taxon[r["taxon"]] = by_taxon.get(r["taxon"], 0) + 1
        by_pmid[r["pmid"]] = by_pmid.get(r["pmid"], 0) + 1
    print(f"  {len(by_taxon)} taxa, {len(by_pmid)} publications")
    for t, n in sorted(by_taxon.items(), key=lambda kv: -kv[1]):
        print(f"    NCBITaxon:{t:<10} {n:6,} rows")

    if bad:
        print("\nREFUSED ROWS ABOVE. Fix them and re-run; a partly parsed")
        print("table is the thing this reports rather than ingests.")
    if args.inspect:
        print("\n--inspect: nothing written.")
        return 0 if not bad else 1
    if bad:
        return 1

    taxa = built_taxa(args.spine)
    gene_to_protein = build_gene_to_proteins(
        [args.host / "edges.jsonl", args.spine / "edges.jsonl"])
    print(f"\n{len(gene_to_protein):,} genes map to host proteins; "
          f"{len(taxa)} taxa built")
    em = Emit()
    s = ingest_curated(em, rows, gene_to_protein, known_taxa=taxa)
    print(f"\n{s['edges']:,} edges  "
          f"({s['dependency']:,} dependency, {s['restriction']:,} restriction)")
    for k in ("taxon_not_built", "gene_not_in_graph", "duplicate"):
        if s[k]:
            print(f"  skipped {k:<20} {s[k]:6,}")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in em.nodes.values()))
    (args.out / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in em.edges))
    print(f"\nwrote {args.out}  (source {CURATED_SOURCE}, evidence tier 2)")
    print("`curated` is already in LAYER_PRECEDENCE, directly below orcs;")
    print("re-assemble to bring this layer into the release.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
