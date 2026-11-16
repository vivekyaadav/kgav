#!/usr/bin/env python3
"""Is every node id what the normaliser derives from that node's own SMILES?

FIX-LIST L27, AND THE ONLY CHECK THAT SHOWS THE REPAIR HELD. The defect was
that node ids came from ChEMBL's `standard_inchi_key` verbatim while every
emit site declared `salt_collapsed: True`. Measured before the repair: 549 of
25,383 SmallMolecule nodes carried a multi-fragment `smiles`, all 549
declared the flag, and NOT ONE was a fixed point of `normalize_chemical`.
287 of them normalised onto a node already in the release, so 287 compounds
were filed twice, 980 edges sat on the duplicate copies, and 31 compounds
were a positive AND a negative for the same virus under two ids -- which
`build_labels` could not see, because the two sides were different subjects.

WHY THIS IS A SCRIPT AND NOT ONLY A TEST. The full pass normalises every
stored structure: an RDKit cleanup, largest-fragment, uncharge and
canonical-tautomer pass per molecule, which is **minutes, not seconds**, on
25,094 compounds. `make test-all` is `pytest -n auto --dist loadfile`, so
every test in one file lands on one worker and the longest file sets the
wall clock -- putting the whole pass in the suite would take the full gate
from four minutes to a quarter of an hour. That is a bad trade for a check
whose job is to confirm a one-time repair.

So the division is explicit: this script does the whole release on demand
(`make identity-fixed-points`), and `tests/test_identity_policy.py` calls the
same function on a deterministic sample plus every node where the source and
the parent disagreed -- the 213 most likely to regress. ONE implementation,
because two would drift, which is the defect this repository keeps finding.

    make identity-fixed-points              # the whole release
    make identity-fixed-points ARGS="--limit 2000"
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.pipeline import release_dir                    # noqa: E402
from kgav.schema import print_capped                    # noqa: E402


def small_molecules(release: Path) -> list[tuple[str, str, bool]]:
    """(node id, stored smiles, source and parent disagreed) per compound."""
    out = []
    for line in (release / "nodes.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        node = json.loads(line)
        if node.get("class") != "SmallMolecule":
            continue
        props = node.get("properties") or {}
        smiles = props.get("smiles") or ""
        if not smiles:
            continue
        out.append((node["id"], smiles,
                    bool(props.get("smiles_as_supplied"))))
    return out


def _keys(batch: list[tuple[str, str]]) -> list[tuple[str, str | None]]:
    """(id, normalised key or None) for a batch, in a worker process.

    Strings in, strings out: an RDKit molecule is not picklable and does not
    need to be, since every id this answers about is derived from a SMILES
    and compared as a string.
    """
    from kgav.normalize_chem import normalize_chemical as N

    out = []
    for nid, smiles in batch:
        try:
            out.append((nid, f"INCHIKEY:{N(smiles).inchikey}"))
        except Exception:
            # BROADLY, AND COUNTED BY THE CALLER. normalize_chemical raises
            # NormalizationError, but RDKit's sanitisation and tautomer
            # passes raise their own types out of C++, and a structure this
            # pass cannot read is not evidence about the ids of the ones it
            # can.
            out.append((nid, None))
    return out


def fixed_point_census(release: Path, limit: int | None = None,
                       seed: int = 0, jobs: int = 0) -> dict:
    """How many node ids the normaliser does NOT re-derive, and which ids
    collide on one normalised key.

    `limit` samples, and the sample is NOT just the first N: it is every
    compound whose source and parent disagreed -- the ones that exercise the
    salt and charge paths at all -- plus a seeded random draw from the rest.
    Taking the first N would scan whichever ids sort lowest, which is a
    property of the InChIKey alphabet and not of the chemistry.

    `collisions` is only meaningful over the WHOLE release: two ids sharing
    a key cannot be seen if the sample holds one of them. It is reported
    with `complete` so a caller cannot read a sampled zero as a clean bill.
    """
    mols = small_molecules(release)
    complete = limit is None or limit >= len(mols)
    if not complete:
        disagreed = [m for m in mols if m[2]]
        rest = [m for m in mols if not m[2]]
        take = max(0, limit - len(disagreed))
        rng = random.Random(seed)
        mols = disagreed + rng.sample(rest, min(take, len(rest)))

    # ACROSS THE CORES, because the whole release took 21m58s on one and
    # this is the check that closes L27 -- a 22-minute confirmation is one
    # nobody runs. The work is per molecule and shares nothing, so it
    # chunks cleanly. `jobs=1` keeps the serial path for debugging, and the
    # result is order-independent: `not_fixed` is sorted before returning,
    # so two runs at different widths produce the same output.
    pairs = [(nid, smiles) for nid, smiles, _ in mols]
    n_jobs = jobs or min(os.cpu_count() or 1, 32)
    if n_jobs > 1 and len(pairs) > 256:
        chunk = max(16, len(pairs) // (n_jobs * 8))
        batches = [pairs[i:i + chunk] for i in range(0, len(pairs), chunk)]
        results: list[tuple[str, str | None]] = []
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            for got in pool.map(_keys, batches):
                results.extend(got)
    else:
        results = _keys(pairs)

    by_key: dict[str, set[str]] = defaultdict(set)
    not_fixed, unparseable = [], 0
    for nid, key in results:
        if key is None:
            unparseable += 1
            continue
        by_key[key].add(nid)
        if key != nid:
            not_fixed.append((nid, key))
    not_fixed.sort()
    return {"checked": len(mols), "unparseable": unparseable,
            "not_fixed": not_fixed,
            "collisions": {k: sorted(v) for k, v in by_key.items()
                           if len(v) > 1},
            "complete": complete, "total_small_molecules": len(
                mols if complete else small_molecules(release))}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=None)
    ap.add_argument("--limit", type=int, default=None,
                    help="sample this many compounds instead of all of them")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=0,
                    help="worker processes; 0 = all cores, 1 = serial")
    args = ap.parse_args(argv)

    release = args.release or release_dir(ROOT)
    if not (release / "nodes.jsonl").exists():
        print(f"no release at {release}")
        return 2

    c = fixed_point_census(release, args.limit, args.seed, args.jobs)
    scope = "the whole release" if c["complete"] else \
            f"a seeded sample of {c['checked']:,}"
    print(f"identity fixed points on {release.name}, {scope}")
    print(f"  SmallMolecule with a stored structure  "
          f"{c['total_small_molecules']:>8,}")
    print(f"  normalised                             {c['checked']:>8,}")
    print(f"  unparseable                            "
          f"{c['unparseable']:>8,}")
    print(f"  NOT a fixed point of the normaliser    "
          f"{len(c['not_fixed']):>8,}   (549 before the L27 repair)")
    # print_capped, not a [:10] slice: ten offenders and four hundred must
    # not read identically, which tests/test_no_duplicate_logic.py sweeps
    # the whole of scripts/ for -- and caught here on the first run.
    print_capped(c["not_fixed"], 10,
                 render=lambda x: f"{x[0]}  ->  {x[1]}", indent="      ")
    if c["complete"]:
        print(f"  normalised keys carried by >1 node     "
              f"{len(c['collisions']):>8,}   (287 before the repair)")
        print_capped(sorted(c["collisions"].items()), 10,
                     render=lambda kv: f"{kv[0]}  <-  {', '.join(kv[1])}",
                     indent="      ")
    else:
        print("  normalised keys carried by >1 node          n/a   "
              "(needs the whole release: a sample can hold one of a pair)")
    ok = not c["not_fixed"] and (not c["complete"] or not c["collisions"])
    print(f"\n  {'PASS' if ok else 'FAIL'}: every id checked is "
          f"{'' if ok else 'NOT '}what the normaliser derives from its own "
          f"stored structure")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
