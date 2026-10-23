#!/usr/bin/env python3
"""Driver: drug -> human protein TARGETS from the Guide to PHARMACOLOGY.

    mkdir -p data/raw/gtop
    # from https://www.guidetopharmacology.org/download.jsp
    #   interactions.csv   and   ligands.csv  (or ligand_id_mapping.csv)
    python scripts/ingest_gtop.py --inspect   # columns, coverage, no write
    python scripts/ingest_gtop.py             # build the layer

--inspect first, always. This reader was built against the published
column list and not against a downloaded file, so the first thing worth
knowing is whether the release you have matches it. --inspect reports the
join rate too: how many interactions survive the human filter, the
UniProt filter and the InChIKey join, which is the number that says
whether the layer is worth assembling.

WHAT IT ADDS, in order of value: a DIRECTION on edges that have none
(12,155 of 14,381 in v0.1 state nothing), then reach for approved and
clinical ligands the chemistry layer missed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.emit import Emit  # noqa: E402
from kgav.gtop import (  # noqa: E402
    GTOP_SOURCE,
    INTERACTION_COLUMNS,
    GtopFormatError,
    _rows,
    ingest_interactions,
    ligand_inchikeys,
)


def graph_ids(release: Path) -> tuple[set[str], set[str]]:
    """(human protein ids, small molecule ids) already in the graph."""
    proteins, compounds = set(), set()
    p = release / "nodes.jsonl"
    if not p.exists():
        return proteins, compounds
    for line in p.read_text().splitlines():
        if not line.strip():
            continue
        n = json.loads(line)
        if n.get("class") == "Protein" and not (
                n.get("properties") or {}).get("is_viral"):
            proteins.add(n["id"])
        elif n.get("class") == "SmallMolecule":
            compounds.add(n["id"])
    return proteins, compounds


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--raw", type=Path, default=ROOT / "data" / "raw" / "gtop")
    ap.add_argument("--release", type=Path,
                    default=ROOT / "data" / "releases" / "v0.1")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "releases" / "v0.1-gtop")
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--mint", action="store_true",
                    help="also emit edges for ligands not already in the "
                         "graph; off by default, see src/kgav/gtop.py")
    args = ap.parse_args(argv)

    inter = args.raw / "interactions.csv"
    ligands = next((args.raw / n for n in
                    ("ligands.csv", "ligand_id_mapping.csv")
                    if (args.raw / n).exists()), None)
    if not inter.exists() or ligands is None:
        print(f"need interactions.csv and ligands.csv (or "
              f"ligand_id_mapping.csv) in {args.raw}")
        print("Download them from "
              "https://www.guidetopharmacology.org/download.jsp")
        return 1

    try:
        rows = _rows(inter, INTERACTION_COLUMNS)
        keys = ligand_inchikeys(ligands)
    except GtopFormatError as e:
        print(f"REFUSED: {e}")
        return 2

    print(f"{len(rows):,} interactions, {len(keys):,} ligands with an "
          f"InChIKey\n")
    proteins, compounds = graph_ids(args.release)
    if not proteins:
        print(f"no nodes under {args.release}; build the release first")
        return 2
    print(f"graph: {len(proteins):,} human proteins, "
          f"{len(compounds):,} small molecules")

    em = Emit()
    s = ingest_interactions(
        em, rows, keys, known_proteins=proteins,
        known_compounds=None if args.mint else compounds)

    print(f"\n{s['edges']:,} TARGETS edges")
    for k in ("direction_inhibitor", "direction_agonist",
              "direction_modulator", "direction_unknown"):
        if s[k]:
            print(f"  {k.removeprefix('direction_'):<12} {s[k]:7,}")
    known = s["edges"] - s["direction_unknown"]
    if s["edges"]:
        print(f"\n  {known:,} of {s['edges']:,} state a direction "
              f"({known / s['edges']:.1%}). v0.1's existing TARGETS layer "
              f"states one on 2,226 of 14,381 (15.5%).")
    print("\n  dropped:")
    for k in ("not_human", "no_uniprot", "multi_accession",
              "protein_not_in_graph", "ligand_without_inchikey",
              "compound_not_in_graph", "duplicate", "unparseable_affinity"):
        if s[k]:
            print(f"    {k:<26} {s[k]:7,}")

    if args.inspect:
        print("\n--inspect: nothing written.")
        return 0
    if not s["edges"]:
        print("\nNo edge survived the joins. Nothing written: an empty layer")
        print("directory would assemble silently and look like coverage.")
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in em.nodes.values()))
    (args.out / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in em.edges))
    print(f"\nwrote {args.out}  (source {GTOP_SOURCE})")
    print("`gtop` sits above `hosttargets` in LAYER_PRECEDENCE, so its")
    print("curated action wins a direction conflict with ChEMBL's bulk rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
