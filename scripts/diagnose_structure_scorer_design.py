#!/usr/bin/env python3
"""The design facts a structure-derived scorer has to be registered against.

WHY THIS EXISTS, AND WHAT IT DELIBERATELY DOES NOT DO.

Every evaluated channel in this project has almost no room to rank. Across
56 cells the largest margin above the imputed floor is +0.0050 and the
median is +0.0000, because 74.8% of labelled compounds carry no traversable
edge, so an AUC is mostly `auc_floor` -- the coverage gap -- and not the
ranking. A structure-derived scorer is the one candidate that reaches every
labelled compound, because SMILES coverage is complete. So it has a span,
and a span is exactly what makes it worth pre-registering rather than
trying.

**THIS SCRIPT COMPUTES NO SCORE AND NO AUC, AND THAT IS THE POINT.** A
pre-registration is written before the measurement; a diagnostic that
printed the thing being registered would make the registration a
description. What it measures is the DESIGN: how large each labelled class
is, whether the compounds can be fingerprinted at all, how much of each
label set is near-duplicate chemistry, how many labels carry a publication,
and what AUC each cell is powered to detect. Those are properties of the
labels and the release. None of them is an outcome.

THE ONE LINE THAT MATTERS MOST. The obvious hazard measurement -- the
distribution of maximum Tanimoto from each compound to the known actives --
is NOT taken here, although the code to take it is two imports away. For a
positive, that quantity IS its leave-one-out score; for a negative it IS
its score. Printing both distributions is printing the ranking, and the
separation between them is the AUC to four decimal places. So the hazard is
quantified by a stricter and cruder proxy that cannot be turned into a
score: EXACT SKELETON IDENTITY, `inchikey_skel`, which asks how many
actives are stereoisomers or tautomers of another active rather than how
similar they are. It understates analog redundancy, and understating it
before signing is the correct direction to be wrong in.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.calibration import MIN_REACHED_POSITIVES     # noqa: E402
from kgav.labels import (                              # noqa: E402
    ACTIVITY_PREDICATE,
    auc_interval,
    build_labels,
)
from kgav.pipeline import release_dir                   # noqa: E402
from kgav.schema import print_capped                    # noqa: E402
from kgav.virus_register import load as load_register   # noqa: E402

# Imported rather than restated: the fingerprint any scorer would use is the
# one the similarity layer already builds, ECFP4 at radius 2 over 2048 bits.
from kgav.similarity import FP_BITS, FP_RADIUS          # noqa: E402


def virus_names() -> dict[str, str]:
    """NCBITaxon id -> the register's name for it.

    NOT a convenience. Writing this registration by hand put the wrong name
    against four taxa -- 11082 as a generic flavivirus when it is West Nile,
    11072 as an enterovirus when it is JEV, 11084 as a morbillivirus when it
    is TBEV -- and named HCoV-NL63 as the virus behind a PMID figure that
    belongs to West Nile. A table of bare taxon ids invites exactly that, so
    the name comes off the register.
    """
    doc = load_register()
    out = {}
    for entry in doc.get("viruses") or []:
        taxon = entry.get("taxon")
        if taxon and entry.get("name"):
            out[f"NCBITaxon:{taxon}"] = entry["name"]
    return out


def small_molecules(release: Path) -> dict[str, dict]:
    """id -> properties, for every SmallMolecule in the release."""
    out = {}
    for line in (release / "nodes.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        n = json.loads(line)
        if n["class"] == "SmallMolecule":
            out[n["id"]] = n.get("properties") or {}
    return out


def label_publications(release: Path) -> dict[str, set[str]]:
    """compound -> PMIDs cited by its activity-label edges."""
    idx: dict[str, set[str]] = defaultdict(set)
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] == ACTIVITY_PREDICATE:
            pubs = {p for p in (e.get("publications") or []) if p}
            if pubs:
                idx[e["subject"]] |= pubs
    return dict(idx)


def skeleton_redundancy(ids: set[str], props: dict[str, dict]) -> dict:
    """How much of a label set is near-duplicate chemistry.

    EXACT skeleton identity, not similarity -- see the module docstring for
    why the similarity version is withheld. `inchikey_skel` is the
    stereo-agnostic grouping key the identity policy already defines, so
    two ids sharing one are the same constitution differing in stereo or
    tautomer.
    """
    skels = Counter(props.get(i, {}).get("inchikey_skel", "") for i in ids)
    unknown = skels.pop("", 0)
    members_in_groups = sum(c for c in skels.values() if c > 1)
    return {
        "n": len(ids),
        "no_skeleton": unknown,
        "distinct_skeletons": len(skels),
        "groups": sum(1 for c in skels.values() if c > 1),
        "members_in_groups": members_in_groups,
        "largest_group": max(skels.values(), default=0),
    }


def detectable_auc(n_pos: int, n_neg: int, alpha_z: float = 1.96) -> float | None:
    """The smallest AUC whose interval excludes 0.5 for this cell.

    Walks up from 0.5 using the project's own Hanley & McNeil interval, so
    the figure is the one the evaluation would actually print. Returns None
    if no AUC below 1.0 can clear chance at this class size, which is the
    "the cell cannot reject" verdict fix-list B3 made load-bearing.
    """
    if n_pos < 1 or n_neg < 1:
        return None
    auc = 0.5
    while auc < 1.0:
        auc = round(auc + 0.0001, 4)
        lo, _hi = auc_interval(auc, n_pos, n_neg, z=alpha_z)
        if lo > 0.5:
            return auc
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--sources", default=None,
                    help="restrict labels to one screen, e.g. infores:ncats")
    args = ap.parse_args()

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    print("STRUCTURE-SCORER DESIGN FACTS. No score is computed here and no")
    print("AUC; see the module docstring for the line this script does not")
    print("cross and why.\n")
    print(f"release {args.release.name}")
    print(f"fingerprint any scorer would reuse: ECFP4, radius {FP_RADIUS}, "
          f"{FP_BITS} bits (kgav.similarity)\n")

    props = small_molecules(args.release)
    sources = {s.strip() for s in args.sources.split(",")} if args.sources else None
    labels = build_labels(args.release, sources=sources)
    pubs = label_publications(args.release)

    viruses = sorted(set(labels.positives) | set(labels.negatives))
    names = virus_names()

    def label(v: str) -> str:
        n = names.get(v)
        return f"{n} ({v.split(':')[1]})" if n else v

    # ---------------------------------------------------- 1. can it even score
    print("1. CAN A STRUCTURE SCORER REACH THE LABELS AT ALL")
    print("   The premise of the whole exercise. A compound with no SMILES")
    print("   cannot be fingerprinted, and is unreachable for this scorer in")
    print("   exactly the way 74.8% of labels are unreachable for the graph.")
    labelled = set()
    for v in viruses:
        labelled |= labels.positives.get(v, set())
        labelled |= labels.negatives.get(v, set())
    with_smiles = {i for i in labelled if (props.get(i, {}).get("smiles") or "")}
    missing = sorted(labelled - with_smiles)
    print(f"\n   {len(labelled):,} distinct labelled compounds")
    print(f"   {len(with_smiles):,} carry a SMILES  "
          f"({len(with_smiles) / len(labelled):.1%})")
    print(f"   {len(missing):,} do not")
    if missing:
        print_capped(missing, 5, indent="     ",
                     unit="more with no SMILES")
    print("\n   A SMILES is not a fingerprint. RDKit can reject a SMILES it")
    print("   cannot parse, and this script does not call it -- the parse")
    print("   rate belongs to the run and is registered as a figure to")
    print("   report, not one to know first.")

    # -------------------------------------------- 2. the power of each cell
    print("\n\n2. WHAT EACH CELL IS POWERED TO DETECT")
    print("   Uncorrected, z=1.96, by the project's own Hanley & McNeil")
    print(f"   interval. A cell needs >= {MIN_REACHED_POSITIVES} of each")
    print("   class to be read at all; a structure scorer reaches every")
    print("   compound it can fingerprint, so 'reached' is the class size.")
    print(f"\n   {'virus':<26} {'pos':>6} {'neg':>7} {'min AUC':>9}  "
          f"{'readable':>8}")
    rows = []
    for v in viruses:
        pos = labels.positives.get(v, set()) & with_smiles
        neg = labels.negatives.get(v, set()) & with_smiles
        det = detectable_auc(len(pos), len(neg))
        readable = (len(pos) >= MIN_REACHED_POSITIVES
                    and len(neg) >= MIN_REACHED_POSITIVES)
        rows.append((v, len(pos), len(neg), det, readable))
        shown = f"{det:.4f}" if det else "none"
        print(f"   {label(v):<26} {len(pos):>6,} {len(neg):>7,} {shown:>9}  "
              f"{'yes' if readable else 'no':>8}")
    n_readable = sum(1 for r in rows if r[4])
    print(f"\n   {n_readable} of {len(rows)} cells are readable, and the "
          f"correction family")
    print("   has to be fixed over those before any of them is computed.")

    # ------------------------------------------------- 3. the analog hazard
    print("\n\n3. HOW MUCH OF EACH LABEL SET IS THE SAME CHEMISTRY")
    print("   THE CENTRAL HAZARD. A max-similarity scorer ranks a compound")
    print("   by its nearest known active, so a label set built from a")
    print("   congeneric series scores itself. Measured by EXACT skeleton")
    print("   identity, which understates it -- a close analog has a")
    print("   different skeleton and still carries the leak.")
    print(f"\n   {'virus':<26} {'class':<4} {'n':>6} {'skels':>7} "
          f"{'in groups':>10} {'largest':>8}")
    for v, _npos, _nneg, _det, readable in rows:
        if not readable:
            continue
        for name, ids in (("pos", labels.positives.get(v, set()) & with_smiles),
                          ("neg", labels.negatives.get(v, set()) & with_smiles)):
            r = skeleton_redundancy(ids, props)
            share = (r["members_in_groups"] / r["n"]) if r["n"] else 0
            print(f"   {label(v):<26} {name:<4} {r['n']:>6,} "
                  f"{r['distinct_skeletons']:>7,} "
                  f"{r['members_in_groups']:>6,} ({share:>3.0%}) "
                  f"{r['largest_group']:>8,}")

    # -------------------------------------------- 4. can we split on papers
    print("\n\n4. IS A PUBLICATION-DISJOINT SPLIT AVAILABLE")
    print("   The only split that separates a congeneric series from its")
    print("   own paper. baselines.label_publication_index already feeds")
    print("   Graph.load(drop_pubs_shared_with=...), so the mechanism")
    print("   exists; what is unknown is whether the labels carry enough")
    print("   PMIDs to use it.")
    print(f"\n   {'virus':<26} {'pos w/ PMID':>12} {'neg w/ PMID':>12} "
          f"{'papers':>8}")
    for v, _pos_n, _neg_n, _det, readable in rows:
        if not readable:
            continue
        pos = labels.positives.get(v, set()) & with_smiles
        neg = labels.negatives.get(v, set()) & with_smiles
        pp = sum(1 for i in pos if pubs.get(i))
        pn = sum(1 for i in neg if pubs.get(i))
        papers = set()
        for i in pos | neg:
            papers |= pubs.get(i, set())
        print(f"   {label(v):<26} {pp:>5,}/{len(pos):<6,} {pn:>5,}/{len(neg):<6,} "
              f"{len(papers):>8,}")

    print("\n\nWHAT THIS DOES NOT DECIDE")
    print("  Nothing here says whether a structure scorer works. It says")
    print("  which cells could show that it does, what the labels would")
    print("  leak into it, and which splits are buildable. The scorer, the")
    print("  reference set, the split and the correction family are")
    print("  decisions, and they belong in docs/prereg_structure_scorer.md")
    print("  BEFORE the first score is computed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
