#!/usr/bin/env python3
"""Which layers disagree about a node property, and whether it is a real
disagreement or two spellings of the same thing.

    python scripts/diagnose_property_conflicts.py
    python scripts/diagnose_property_conflicts.py --property smiles

WHY. assemble_release.py reports `property_conflict:smiles 7` and prints the
pairs truncated to 34 characters:

    INCHIKEY:HJORMJIFDVBMOB-LBPRGKRZSA-N.smiles:
        'COc1ccc([C@@H]2CNC(=O)C2)cc1OC1CCC' (chembl)
     vs 'COC1=C(OC2CCCC2)C=C(C=C1)[C@@H]3CN' (ncats)

That is unreadable as chemistry. One writes aromatic rings lowercase and the
other Kekule, the atom order differs, and SMILES stereo descriptors (@ / @@)
are ORDER-DEPENDENT -- the same centre is written either way depending on
which neighbour comes first. So two SMILES that look different by eye may be
the same molecule, and two that look similar may not be. Counting them as
seven conflicts says nothing about whether anything is wrong.

WHAT THIS DECIDES. For a SmallMolecule the node id IS the structure: an
InChIKey is derived from the molecule. So each candidate SMILES can be
checked against the id it was filed under, which turns an unreadable string
comparison into a verdict:

    agree              same molecule, two spellings -- nothing is wrong, and
                       the merge is choosing between equals
    stereo only        same skeleton, different stereochemistry. For a drug
                       that is a different compound, and the first layer in
                       LAYER_PRECEDENCE silently wins today
    different molecule one of them is filed under an InChIKey that is not
                       its own structure, which is a broken identifier, not
                       a formatting difference

The last case cannot be read off the truncated print at all, which is the
reason this exists.

NEEDS rdkit (the `chem` extra). Without it the layer disagreements are still
listed in full, with no verdict -- a partial answer, said to be partial,
rather than none.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.assemble import LAYER_PRECEDENCE

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.inchi import MolToInchiKey
    RDLogger.DisableLog("rdApp.*")      # a bad SMILES is this script's output
    HAVE_RDKIT = True
except ImportError:                      # pragma: no cover - env dependent
    HAVE_RDKIT = False


def _canonical(smiles: str) -> tuple[str | None, str | None, str | None]:
    """(canonical SMILES, canonical SMILES without stereo, InChIKey)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    flat = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    Chem.RemoveStereochemistry(flat)
    return Chem.MolToSmiles(mol), Chem.MolToSmiles(flat), MolToInchiKey(mol)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--releases", type=Path, default=root / "data" / "releases")
    ap.add_argument("--prefix", default="v0.1-")
    ap.add_argument("--property", default=None,
                    help="only this property (default: every one that clashes)")
    args = ap.parse_args()

    # Read the LAYER releases, not the assembled one: by the time layers are
    # merged the losing value is gone. _merge_properties fills gaps and keeps
    # the FIRST non-None value, so the assembled release holds whichever
    # layer came first in LAYER_PRECEDENCE and nothing records the other.
    per_node: dict[str, dict[str, dict]] = defaultdict(dict)
    found = []
    for layer in LAYER_PRECEDENCE:
        nf = args.releases / f"{args.prefix}{layer}" / "nodes.jsonl"
        if not nf.exists():
            continue
        found.append(layer)
        for line in nf.read_text().splitlines():
            if not line.strip():
                continue
            n = json.loads(line)
            per_node[n["id"]][layer] = n.get("properties") or {}
    if not found:
        print(f"no {args.prefix}* layer releases in {args.releases}",
              file=sys.stderr)
        return 2
    print(f"layers read, in precedence order: {', '.join(found)}")
    print(f"{len(per_node):,} nodes carry properties in at least one layer")
    if not HAVE_RDKIT:
        print("\nrdkit NOT INSTALLED -- disagreements are listed, but no "
              "verdict is given.\n  pip install 'rdkit>=2024.3'")

    clashes: dict[str, list[tuple[str, dict[str, str]]]] = defaultdict(list)
    for nid, by_layer in per_node.items():
        vals: dict[str, dict[str, str]] = defaultdict(dict)
        for layer, props in by_layer.items():
            for k, v in props.items():
                if v is not None:
                    vals[k][layer] = v
        for k, per_layer in vals.items():
            if args.property and k != args.property:
                continue
            if len({json.dumps(v, sort_keys=True) for v in per_layer.values()}) > 1:
                clashes[k].append((nid, per_layer))

    if not clashes:
        print("\nNo property is given two different values by two layers.")
        return 0

    verdicts: dict[str, int] = defaultdict(int)
    for prop in sorted(clashes, key=lambda k: -len(clashes[k])):
        rows = clashes[prop]
        print(f"\n{'=' * 74}\n{prop}: {len(rows)} node(s)\n{'=' * 74}")
        for nid, per_layer in sorted(rows):
            print(f"\n  {nid}")
            for layer, v in per_layer.items():
                print(f"    {layer:<14} {v}")
            if not (HAVE_RDKIT and prop == "smiles"):
                continue

            # THE NODE ID IS THE REFERENCE. An InChIKey is derived from the
            # structure, so "does this SMILES belong under this id" is a
            # question with an answer, unlike "do these two strings differ".
            claimed = nid.split(":", 1)[1] if ":" in nid else nid
            seen: dict[str, tuple] = {}
            for layer, v in per_layer.items():
                canon, flat, key = _canonical(v)
                if canon is None:
                    print(f"    -> {layer}: rdkit cannot parse this SMILES")
                    verdicts["unparseable SMILES"] += 1
                    continue
                seen[layer] = (canon, flat, key)
                mark = "matches the node id" if key == claimed else \
                       f"NOT this node's key -- resolves to {key}"
                print(f"    -> {layer}: {mark}")
            if len(seen) < 2:
                continue
            canons = {c for c, _, _ in seen.values()}
            flats = {f for _, f, _ in seen.values()}
            if len(canons) == 1:
                v = "agree (two spellings of one molecule)"
            elif len(flats) == 1:
                v = "STEREOCHEMISTRY DIFFERS (same skeleton)"
            else:
                v = "DIFFERENT MOLECULES under one InChIKey"
            verdicts[v] += 1
            print(f"    => {v}")

    if verdicts:
        print(f"\n{'=' * 74}\nVERDICTS\n{'=' * 74}")
        for v, n in sorted(verdicts.items(), key=lambda kv: -kv[1]):
            print(f"  {v:<42} {n:>4}")
        print("""
  agree                 the merge is choosing between equals. Canonicalising
                        at ingest would remove the conflict entirely -- it
                        is a formatting difference being counted as a data
                        disagreement.
  STEREOCHEMISTRY       a different compound for any pharmacological
                        purpose. _merge_properties keeps the FIRST layer in
                        LAYER_PRECEDENCE and reports the clash; which layer
                        is right is a sourcing question, not a merge rule.
  DIFFERENT MOLECULES   an identifier bug. The InChIKey is derived from the
                        structure, so two structures cannot share one. The
                        '-> matches the node id' lines above say which
                        layer filed it wrongly.""")

    print("""
WHAT THIS DOES NOT DECIDE
  Nothing here changes a value. Which source to believe is a curation
  question; this reports what the layers say and what the identifier says
  about them.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
