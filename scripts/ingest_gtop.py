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
(12,155 of 14,381 stated nothing before this layer; 10,395 of 17,982 do
on the assembled v0.1, which contains it), then reach for approved and
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
    CompoundIndex,
    GtopFormatError,
    _rows,
    file_version,
    ingest_interactions,
    ligand_chembl_ids,
    ligand_inchikeys,
    parse_version,
)
from kgav.provenance import SOURCE_STAMP_NAME, write_source  # noqa: E402


def graph_ids(release: Path) -> tuple[set[str], CompoundIndex]:
    """(human protein ids, the compound index) for a built release."""
    proteins: set[str] = set()
    nodes: list[dict] = []
    p = release / "nodes.jsonl"
    if not p.exists():
        return proteins, CompoundIndex(frozenset())
    for line in p.read_text().splitlines():
        if not line.strip():
            continue
        n = json.loads(line)
        if n.get("class") == "Protein" and not (
                n.get("properties") or {}).get("is_viral"):
            proteins.add(n["id"])
        elif n.get("class") == "SmallMolecule":
            nodes.append(n)
    return proteins, CompoundIndex.from_nodes(nodes)


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

    # BOTH FILES, AND THEY MUST AGREE. The layer joins interactions to
    # ligands, so two downloads from different releases would silently
    # produce a graph that is neither -- ligand ids are reused across
    # releases, so the join succeeds and the result is wrong. Only
    # interactions.csv was ever read for the version, which made the
    # mismatch undetectable as well as unreported.
    vers = {f.name: parse_version(file_version(f)) for f in (inter, ligands)}
    declared = {v["raw"] for v in vers.values()}
    if len(declared) > 1:
        print("REFUSED: the two downloads declare different GtoPdb releases, "
              "and ligand ids are reused across releases, so the join would "
              "succeed and mean nothing:")
        for name, v in vers.items():
            print(f"  {name:<20} {v['raw'] or 'no version declared'}")
        print("Re-download both from the same release.")
        return 2
    version = vers[inter.name]
    print(version["raw"] or "no GtoPdb version declared in either file")

    try:
        rows = _rows(inter, INTERACTION_COLUMNS)
        keys = ligand_inchikeys(ligands)
        chembl = ligand_chembl_ids(ligands)
    except GtopFormatError as e:
        print(f"REFUSED: {e}")
        return 2

    print(f"{len(rows):,} interactions, {len(keys):,} ligands with an "
          f"InChIKey, {len(chembl):,} with a ChEMBL id\n")
    proteins, compounds = graph_ids(args.release)
    if not proteins:
        print(f"no nodes under {args.release}; build the release first")
        return 2
    print(f"graph: {len(proteins):,} human proteins, "
          f"{len(compounds.exact):,} small molecules "
          f"({len(compounds.by_skeleton):,} unambiguous skeletons, "
          f"{len(compounds.by_chembl):,} ChEMBL ids)")

    em = Emit()
    s = ingest_interactions(
        em, rows, keys, known_proteins=proteins,
        compounds=None if args.mint else compounds,
        ligand_chembl=chembl)

    print(f"\n{s['edges']:,} TARGETS edges")
    for k in ("direction_inhibitor", "direction_agonist",
              "direction_modulator", "direction_unknown"):
        if s[k]:
            print(f"  {k.removeprefix('direction_'):<12} {s[k]:7,}")
    known = s["edges"] - s["direction_unknown"]
    if s["edges"]:
        # THE BASELINE IS PRE-GTOP, and saying "v0.1's existing TARGETS
        # layer" was wrong twice over: v0.1 now CONTAINS this layer, so
        # the comparison was against a release that already includes the
        # thing being measured, and the figure it quoted was from the
        # release before. Labelled by what it is now.
        print(f"\n  {known:,} of {s['edges']:,} state a direction "
              f"({known / s['edges']:.1%}). Before this layer existed, "
              f"the ChEMBL-only TARGETS layer stated one on 2,226 of "
              f"14,381 (15.5%); on the assembled v0.1 it is 7,587 of "
              f"17,982 (42.2%).")
    for k in ("compound_by_skeleton", "compound_by_chembl_id"):
        if s[k]:
            print(f"  {k.removeprefix('compound_by_'):<12} {s[k]:7,} edges "
                  f"joined by a looser key than the exact InChIKey")
    for k, label in (("duplicate_changed_direction",
                      "pairs whose direction the second row changed"),
                     ("qualifier_irreconcilable:direction",
                      "pairs where GtoP's own rows contradict each other")):
        if s[k]:
            print(f"  {s[k]:7,} {label}")
    print(f"\n  {s['no_affinity_reported']:7,} row(s) report no affinity at "
          f"all, which is not the same as\n          an affinity this reader "
          f"could not express")
    if s["affinity_unit_unmapped"]:
        print(f"  {s['affinity_unit_unmapped']:7,} potency value(s) dropped: "
              f"their affinity unit is not in AFFINITY_UNIT_TO_QUALIFIER.")
        print("          The edge still ships, with a direction and no "
              "number. By unit:")
        for k in sorted(k for k in s if k.startswith("affinity_unit:")):
            print(f"            {s[k]:>5,}  {k.removeprefix('affinity_unit:')}")
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
    # THE RELEASE, WRITTEN DOWN. Parsed since the quoted-comment fix and
    # printed to a terminal that nobody kept, so the v0.1 MANIFEST names
    # no GtoPdb release and no figure derived from this layer can say
    # which one it came from. `published` is a date the source states
    # about itself -- the only one in this project that is not a download
    # mtime. Absent fields are omitted, never written as "unrecorded".
    write_source(args.out, guidetopharmacology={
        k: v for k, v in version.items() if v is not None} or None)
    print(f"\nwrote {args.out}  (source {GTOP_SOURCE})")
    print(f"recorded GtoPdb release "
          f"{version['release'] or 'unknown'} in {SOURCE_STAMP_NAME}")
    print("`gtop` sits above `hosttargets` in LAYER_PRECEDENCE, so its")
    print("curated action wins a direction conflict with ChEMBL's bulk rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
