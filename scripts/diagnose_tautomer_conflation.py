#!/usr/bin/env python3
"""Does ONE layer file two tautomers under one InChIKey? (fix-list L19)

    python scripts/diagnose_tautomer_conflation.py --layer ncats
    python scripts/diagnose_tautomer_conflation.py --layer chembl --db <path>

WHY THIS QUERY EXISTS AND WHY IT HAD NOT BEEN RUN. L19 found 32 nodes
carrying the drug-target edges of two tautomers -- `ncats` and
`hosttargets` disagreeing about which annular nitrogen holds the proton,
both spellings resolving to the node's key, because standard InChI's
mobile-H layer treats annular tautomers as one compound. It was a
CROSS-LAYER diagnostic, and the only build with two layers describing the
same compound was B1, which was declined. So on the shipped release that
diagnostic reports nothing, and the row says in as many words that whether
a SINGLE layer also does it "is a different query and has not been run".

The two answers are opposite claims about v0.1:

    nothing found   the conflation is LATENT -- a property of merging two
                    sources, which this release does not do. L19 is a
                    limit to state in scope_and_limits.md.
    something found the conflation is IN v0.1 TODAY. One node is
                    accumulating two molecules' edges with one source, and
                    L19 is a defect to fix rather than a limit to declare.

HOW IT ASKS. A layer writes one node per id -- measured, every per-layer
nodes.jsonl has zero ids carrying two SMILES -- so the collapse cannot be
seen in the output. It happens when the layer DERIVES one id from two raw
rows. So this reads the raw input, computes `normalize_chemical(smiles)`
exactly as the ingest does, groups the raw spellings by the resulting
InChIKey, and inside each group asks `same_structure`, which is
deliberately weaker than the id: it normalises aromaticity, atom order and
stereo NOTATION while keeping tautomers apart.

A group whose members share an InChIKey and which `same_structure` splits
is the thing L19 describes. Spelling differences are not reported, because
`same_structure` calls them one molecule, which is what they are.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.normalize_chem import normalize_chemical, same_structure  # noqa: E402
from kgav.schema import print_capped  # noqa: E402

#: Layers that supply a SMILES, and therefore can conflate on their own.
#: A layer that supplies no structure cannot derive an id from one.
SMILES_LAYERS = ("ncats", "chembl", "hosttargets")

#: Where ncats keeps its structures. Two spellings of one header, because
#: the CSV and the TSV disagree about case, which is `ncats.COLUMNS`' job
#: everywhere else.
NCATS_FILES = ("cpe_assay.csv", "cpe_tox.csv", "cpe.tsv")


def ncats_rows(root: Path) -> list[tuple[str, str]]:
    """(source id, raw SMILES) for every ncats row that carries one."""
    out: list[tuple[str, str]] = []
    for name in NCATS_FILES:
        p = root / "data/raw/ncats" / name
        if not p.exists():
            continue
        delim = "\t" if p.suffix == ".tsv" else ","
        with p.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh, delimiter=delim):
                keys = {k.lower().replace("_", ""): k for k in row if k}
                sk = keys.get("smiles") or keys.get("structure")
                ik = keys.get("sampleid") or keys.get("sampleid")
                if not sk:
                    continue
                smi = (row.get(sk) or "").strip()
                if smi:
                    out.append((f"{name}:{(row.get(ik) or '?').strip()}", smi))
    return out


def chembl_rows(db: Path, layer: str) -> list[tuple[str, str]]:
    """(chembl_id, canonical_smiles) as each layer's own query selects it."""
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        cur = con.execute(
            "SELECT md.chembl_id, cs.canonical_smiles "
            "FROM molecule_dictionary md "
            "JOIN compound_structures cs ON cs.molregno = md.molregno "
            "WHERE cs.canonical_smiles IS NOT NULL")
        return [(a, b) for a, b in cur if b]
    finally:
        con.close()


#: What distinguishes the spellings inside one conflated key. The first
#: measurement of this on ncats came back dominated by the second class,
#: which is NOT what L19 is about, and reporting one number would have
#: answered a question nobody asked.
SALT_OR_CHARGE = "salt_or_charge"
TAUTOMER = "tautomer"


def classify(spellings: list[str], ids: dict[str, object]) -> str:
    """Why `same_structure` split these, given their ChemIds.

    `normalize_chemical` strips counterions (largest fragment) and
    neutralises (uncharge) BEFORE computing the InChIKey, and
    `same_structure` deliberately does neither. So two spellings can share
    a key and be called different molecules for two unrelated reasons:

      SALT_OR_CHARGE  one is a salt or a charged form of the other --
                      `Cl.CN(C)CC...` against the free base. This is
                      fix-list L18's class, not L19's: the id is right, the
                      stored SMILES describes a different substance, and
                      nothing about the edges is conflated, because both
                      rows ARE the same parent compound.

      TAUTOMER        both are single neutral fragments and the proton has
                      moved. This is L19: standard InChI's mobile-H layer
                      treats them as one compound, so one node carries two
                      molecules' edges.

    The flags come from ChemId, which already records exactly this:
    `had_multiple_fragments` is set when the input had a counterion to
    strip, `was_charged` when it had to be neutralised.
    """
    salty = any(getattr(ids[s], "had_multiple_fragments", False)
                or getattr(ids[s], "was_charged", False) for s in spellings)
    return SALT_OR_CHARGE if salty else TAUTOMER


def conflations(rows: list[tuple[str, str]]) -> tuple[dict, dict]:
    """Group raw spellings by derived InChIKey; split each by same_structure.

    Returns (conflated, stats). `conflated` maps InChIKey -> a dict with the
    representative spellings, their source ids, and the `kind` that says
    whether this is L19's tautomer conflation or L18's salt class.
    """
    by_key: dict[str, dict[str, list[str]]] = collections.defaultdict(
        lambda: collections.defaultdict(list))
    chem: dict[str, object] = {}
    failed = 0
    for src, smi in rows:
        if smi not in chem:
            try:
                chem[smi] = normalize_chemical(smi)
            except Exception:
                chem[smi] = None
        cid = chem[smi]
        if cid is None:
            failed += 1
            continue
        by_key[cid.inchikey][smi].append(src)

    conflated = {}
    for key, spellings in by_key.items():
        if len(spellings) < 2:
            continue
        groups: list[list[str]] = []
        for s in spellings:
            for g in groups:
                if same_structure(g[0], s):
                    g.append(s)
                    break
            else:
                groups.append([s])
        if len(groups) > 1:
            reps = [g[0] for g in groups]
            conflated[key] = {
                "kind": classify(reps, chem),
                "spellings": {r: spellings[r] for r in reps},
            }
    kinds = collections.Counter(v["kind"] for v in conflated.values())
    return conflated, {
        "rows": len(rows), "keys": len(by_key), "unparseable": failed,
        "keys_with_several_spellings": sum(1 for v in by_key.values() if len(v) > 1),
        "tautomer": kinds[TAUTOMER], "salt_or_charge": kinds[SALT_OR_CHARGE],
    }


def _render_key(item) -> str:
    """One conflated key and its spellings, for print_capped.

    print_capped rather than a `[:40]` slice and a print: this repository
    found 14 capped printing loops of which one said what it had left out,
    and `tests/test_no_duplicate_logic.py` now rejects the other shape. It
    rejected this script on its first run, which is the test working.
    """
    key, v = item
    lines = [f"\n  {key}"]
    for smi, srcs in v["spellings"].items():
        shown = ", ".join(srcs[:3]) + (" +more" if len(srcs) > 3 else "")
        lines.append(f"    {smi}   <- {shown}")
    return "\n".join(lines)


def build_parser(root: Path | None = None) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--layer", choices=SMILES_LAYERS, default="ncats")
    ap.add_argument("--db", type=Path, default=None,
                    help="chembl sqlite, required for the chembl layers")
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after N raw rows; 0 reads all of them")
    ap.add_argument("--out", type=Path, default=None,
                    help="write the two classes as JSON")
    return ap


def main(argv=None) -> int:
    args = build_parser(ROOT).parse_args(argv)
    if args.layer == "ncats":
        rows = ncats_rows(ROOT)
    else:
        if args.db is None or not args.db.exists():
            print(f"--db is required for layer {args.layer} and must exist",
                  file=sys.stderr)
            return 1
        rows = chembl_rows(args.db, args.layer)
    if args.limit:
        rows = rows[:args.limit]
    if not rows:
        print(f"no raw SMILES found for layer {args.layer}", file=sys.stderr)
        return 1

    conflated, stats = conflations(rows)
    taut = {k: v for k, v in conflated.items() if v["kind"] == TAUTOMER}
    salt = {k: v for k, v in conflated.items() if v["kind"] == SALT_OR_CHARGE}

    print(f"layer {args.layer}")
    print(f"  {stats['rows']:,} raw row(s) with a SMILES")
    print(f"  {stats['keys']:,} distinct InChIKey(s) derived")
    print(f"  {stats['unparseable']:,} unparseable")
    print(f"  {stats['keys_with_several_spellings']:,} key(s) written more "
          f"than one way at all (spelling INCLUDED)")
    print(f"  {len(conflated):,} key(s) where same_structure calls the "
          f"spellings DIFFERENT MOLECULES, of which")
    print(f"      {len(taut):,}  TAUTOMER       <- this is L19")
    print(f"      {len(salt):,}  SALT_OR_CHARGE <- this is L18's class, "
          f"and the id is right")
    print()

    if args.out:
        args.out.write_text(json.dumps({
            "layer": args.layer, "stats": stats,
            "tautomer": {k: v["spellings"] for k, v in taut.items()},
            "salt_or_charge": {k: v["spellings"] for k, v in salt.items()},
        }, indent=1, sort_keys=True) + "\n")
        print(f"wrote {args.out}\n")

    if not taut:
        print("NO TAUTOMER CONFLATION FROM A SINGLE LAYER. On this input "
              "every id this layer\nderives comes from one molecule, so "
              "L19's conflation is a property of MERGING\ntwo sources and "
              "is LATENT in v0.1 rather than present in it. L19 is then a "
              "limit\nto state, not a defect to fix.")
    else:
        print(f"L19 IS PRESENT IN v0.1, NOT LATENT. {len(taut)} node(s) "
              f"below take their edges from\ntwo tautomers supplied by THIS "
              f"LAYER ALONE, so no merge is needed to produce the\n"
              f"conflation and declining B1 did not avoid it.")
        print_capped(sorted(taut.items()), 40, _render_key, indent="  ",
                     unit="more node(s) with a tautomer conflation")

    if salt:
        print(f"\nAND SEPARATELY, {len(salt)} key(s) are L18 AT SCALE. "
              f"That row is filed at 4 nodes,\nmeasured on the declined B1 "
              f"build with a cross-layer diagnostic. Within ncats\nalone "
              f"the layer supplies both a salt and its free base for "
              f"{len(salt)} compounds. The\nnode id is correct in every "
              f"case -- normalize_chemical strips the counterion --\nso "
              f"nothing is conflated; what is wrong is the SMILES stored "
              f"beside it, for\nwhichever row happened to be written last.")
        print_capped(sorted(salt.items()), 5, _render_key, indent="  ",
                     unit="more key(s) in the salt class")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
