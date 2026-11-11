#!/usr/bin/env python3
"""What the graph suggests, by name, for each virus.

THIS IS THE READ-OUT, NOT A HYPOTHESIS TEST. Everything upstream of it has
been about whether the structure-derived scorer measures anything: three
pre-registrations, their verdicts, and the controls in
docs/prereg_structure_scorer_strength_matched.md. This script does not add a
protocol. It points the scorer that survived those tests at the question the
graph exists for and prints compound names.

WHAT IT RANKS. For each virus with enough known actives, every small
molecule in the release that is NOT already a measured active or inactive
for that virus, scored by maximum Tanimoto to that virus's known actives
under the same exclusions the registered protocol uses. The nearest known
active is printed beside each suggestion, because a similarity score with
no neighbour named is unauditable: the whole claim is "this looks like that",
and a reader has to see the "that".

WHAT IT IS WORTH, IN THE SAME BREATH AS THE LIST. On the labelled
populations this scorer reads AUC 0.78 on SARS-CoV-2 and clears a corrected
threshold in ten of fourteen cells, with zero below chance -- but the
plate-matched cell is marginal and the collection confound between
source-matched and plate-matched populations measures about 0.19 to 0.22
AUC. So a high rank here is a chemical-analogy suggestion with a measured
hit rate, not a prediction of activity and not a mechanism. The header
prints that every time; fix-list A4 exists because a ranking that does not
say so is making a claim the evidence refuses.

MEASURED ACTIVES AND INACTIVES ARE BOTH EXCLUDED FROM THE SUGGESTIONS. An
active is not a repurposing candidate, it is a known answer; an inactive was
tested and failed, so suggesting it would be proposing an experiment that
has already been done. Both counts are printed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import Graph                        # noqa: E402
from kgav.labels import build_labels                    # noqa: E402
from kgav.pipeline import release_dir                   # noqa: E402
from kgav.provenance import write_results               # noqa: E402
from kgav.schema import load_schema                     # noqa: E402
from kgav.similarity import FP_BITS, FP_RADIUS, fingerprints  # noqa: E402
from kgav.virus_register import load as load_register   # noqa: E402

CHEMBL_DB = ROOT / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db"

MIN_ACTIVES = 10          # calibration.MIN_REACHED_POSITIVES
DEFAULT_TOP = 15


def virus_names() -> dict[str, str]:
    doc = load_register()
    return {f"NCBITaxon:{e['taxon']}": e["name"]
            for e in (doc.get("viruses") or [])
            if e.get("taxon") and e.get("name")}


def preferred_names(chembl_ids: set[str]) -> dict[str, str]:
    """ChEMBL id -> pref_name, for the compounds being suggested.

    THE GRAPH HOLDS NO NAME AND THAT IS WHY THIS IS HERE. A SmallMolecule
    node carries `chembl_id`, `smiles`, `inchikey_skel` and the identity
    flags -- nothing a person can read. A candidate list of accession
    numbers cannot be judged by anyone, which is the difference between a
    result and a read-out. The names are looked up from the raw ChEMBL
    database rather than added to the release, because a display label is
    not a fact the graph asserts.

    Absent database is not fatal: the rows print their ChEMBL ids.
    """
    if not CHEMBL_DB.exists() or not chembl_ids:
        return {}
    import sqlite3
    con = sqlite3.connect(f"file:{CHEMBL_DB}?mode=ro", uri=True)
    out = {}
    ids = sorted(chembl_ids)
    for i in range(0, len(ids), 500):
        chunk = ids[i:i + 500]
        q = ("SELECT chembl_id, pref_name FROM molecule_dictionary "
             f"WHERE chembl_id IN ({','.join('?' * len(chunk))})")
        for cid, name in con.execute(q, chunk):
            if name:
                out[cid] = name
    con.close()
    return out


def parent_of(smiles: str) -> str | None:
    """The node id this structure WOULD have under the identity policy.

    WHY A SUGGESTION NEEDS THIS. The first run of this script ranked
    ERLOTINIB HYDROCHLORIDE first for DENV-1 at Tanimoto 0.979, resembling
    ERLOTINIB -- the same drug, filed as a second node because its SMILES
    kept the `.Cl` and so hashed to its own InChIKey. 549 nodes in the
    release are multi-fragment like that, 287 of them duplicate a node
    already present, and none is a fixed point of `normalize_chemical`.
    Suggesting a salt of a known active is not a repurposing candidate, so
    those rows are dropped here and counted.

    The defect itself is fix-list L27 and is not repaired by this script.
    """
    from kgav.normalize_chem import normalize_chemical
    cid = normalize_chemical(smiles)
    return f"INCHIKEY:{cid.inchikey}" if cid else None


def molecules(release: Path) -> dict[str, dict]:
    out = {}
    for line in (release / "nodes.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        n = json.loads(line)
        if n["class"] == "SmallMolecule":
            out[n["id"]] = n.get("properties") or {}
    return out


def suggest(candidates: list[str], actives: list[str], fps: dict,
            props: dict, top: int, degree: dict | None = None) -> list[dict]:
    """Top `top` candidates by nearest known active.

    The two exclusions that survive from the registered protocol and mean
    something here: a candidate is never scored against itself (it is not in
    `actives` by construction), and a reference sharing its `inchikey_skel`
    is dropped -- a stereoisomer or tautomer of a known active is the same
    compound in this graph's identity policy, so scoring against it would
    report a Tanimoto of 1.0 and suggest a compound already known.
    """
    from rdkit import DataStructs

    ref_ids = [a for a in actives if a in fps]
    ref_fps = [fps[a] for a in ref_ids]
    ref_skel = [props.get(a, {}).get("inchikey_skel", "") for a in ref_ids]
    rows = []
    for c in candidates:
        fp = fps.get(c)
        if fp is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
        c_skel = props.get(c, {}).get("inchikey_skel", "")
        best, best_i = -1.0, None
        for i, s in enumerate(sims):
            if c_skel and ref_skel[i] == c_skel:
                continue
            if s > best:
                best, best_i = s, i
        if best_i is None:
            continue
        rows.append({"candidate": c,
                     "chembl_id": props.get(c, {}).get("chembl_id"),
                     "is_approved": bool(props.get(c, {}).get("is_approved")),
                     # WHETHER THE GRAPH HAS ANYTHING ELSE TO SAY ABOUT IT.
                     # The structure score reaches every compound; the graph
                     # channels reach a quarter of them. A suggestion with
                     # degree 0 rests on chemistry alone, and one with a
                     # target carries a second, independent kind of
                     # evidence. Label edges are excluded from this count
                     # (degree_unlabelled), so it cannot be read off the
                     # answer sheet -- fix-list L7.
                     "graph_degree": int((degree or {}).get(c, 0)),
                     "tanimoto": round(best, 4),
                     "nearest_active": ref_ids[best_i],
                     "nearest_active_chembl":
                         props.get(ref_ids[best_i], {}).get("chembl_id")})
    rows.sort(key=lambda r: -r["tanimoto"])
    return rows[:top]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--top", type=int, default=DEFAULT_TOP)
    ap.add_argument("--approved-only", action="store_true",
                    help="restrict suggestions to compounds the release "
                         "flags is_approved -- the repurposing question in "
                         "its narrowest form")
    ap.add_argument("--virus", default=None)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data/results/repurposing_candidates.json")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    print("REPURPOSING CANDIDATES — the structure scorer, pointed at the")
    print("question rather than at its own validation.\n")
    print("HOW MUCH THIS IS WORTH, BEFORE THE NAMES:")
    print("  On labelled compounds this scorer reads AUC 0.7841 on")
    print("  SARS-CoV-2 and clears a corrected threshold above chance in 10")
    print("  of 14 viruses, none below -- see")
    print("  docs/prereg_structure_scorer_strength_matched.md. But the one")
    print("  population matched by assay rather than by source is MARGINAL")
    print("  (0.5898, interval containing chance) and the gap between the")
    print("  two is 0.19-0.22 AUC. So every row below is a chemical-analogy")
    print("  suggestion with a measured hit rate. It is NOT a prediction of")
    print("  activity and NOT a mechanism: the evidence is 'resembles a")
    print("  known active', and the known active is named so you can judge.\n")

    props = molecules(args.release)
    schema = load_schema()
    g = Graph.load(args.release, symmetry=schema.symmetry())
    labels = build_labels(args.release)          # every source: more actives
    names = virus_names()

    fps, stats = fingerprints(
        {i: (p.get("smiles") or "") for i, p in sorted(props.items())})
    print(f"fingerprint ECFP4 radius {FP_RADIUS}, {FP_BITS} bits: "
          f"{stats['fingerprinted']:,} of {len(props):,} molecules "
          f"({stats['unparseable']:,} unparseable, "
          f"{stats['no_smiles']:,} no SMILES)")

    pool = sorted(i for i in props if i in fps)
    if args.approved_only:
        pool = [i for i in pool if props[i].get("is_approved")]
        print(f"restricted to {len(pool):,} approved compounds")
    print(f"candidate pool: {len(pool):,}\n")

    payload = {"measured_by": "scripts/repurposing_candidates.py",
               "scorer": "max Tanimoto to known actives (ECFP4 2048/r2)",
               "registered_by":
                   "docs/prereg_structure_scorer_strength_matched.md",
               "approved_only": bool(args.approved_only),
               "top": args.top, "viruses": {}}

    for virus in sorted(set(labels.positives) | set(labels.negatives)):
        if args.virus and args.virus not in (virus, names.get(virus)):
            continue
        actives = sorted(labels.positives.get(virus, set()))
        inactives = labels.negatives.get(virus, set())
        if len(actives) < MIN_ACTIVES:
            continue
        known = set(actives) | inactives
        # AND NOT A SALT OF SOMETHING KNOWN, IN EITHER DIRECTION.
        #
        # The first version only normalised CANDIDATES, and the mirror case
        # went straight through: SIMEPREVIR ranked second for SARS-CoV-2
        # resembling SIMEPREVIR SODIUM at Tanimoto 0.832 -- the salt is the
        # known active and the parent is the "suggestion", so the row says
        # nothing but "you already know this drug works". Caught by reading
        # the names, which is the whole argument for printing them.
        #
        # So both sides are mapped: a candidate is dropped if its own
        # parent is known, and if it IS the parent of something known.
        known_parents = set()
        for k in known:
            s = props.get(k, {}).get("smiles") or ""
            if "." in s:
                par = parent_of(s)
                if par:
                    known_parents.add(par)
        candidates, salts = [], 0
        for c in pool:
            if c in known:
                continue
            if c in known_parents:
                salts += 1
                continue
            s = props[c].get("smiles") or ""
            if "." in s:
                par = parent_of(s)
                if par in known or (par in pool and par != c):
                    salts += 1
                    continue
            candidates.append(c)
        rows = suggest(candidates, actives, fps, props, args.top,
                       degree=g.degree_unlabelled)
        payload["viruses"][virus] = {
            "name": names.get(virus),
            "n_known_actives": len(actives),
            "n_known_inactives": len(inactives),
            "n_candidates": len(candidates),
            "n_salts_of_known_compounds_dropped": salts,
            "suggestions": rows,
        }
        label = names.get(virus, virus)
        print(f"{'=' * 72}")
        print(f"{label}  ({virus})")
        print(f"  {len(actives):,} known actives, {len(inactives):,} known "
              f"inactives excluded, {len(candidates):,} candidates")
        print(f"  {salts:,} dropped as the same compound as something "
              f"already known, in either direction (fix-list L27)")
        print(f"{'=' * 72}")
        wanted = {r["chembl_id"] for r in rows if r["chembl_id"]}
        wanted |= {r["nearest_active_chembl"] for r in rows
                   if r["nearest_active_chembl"]}
        pref = preferred_names(wanted)
        for r in rows:
            r["name"] = pref.get(r["chembl_id"] or "")
            r["nearest_active_name"] = pref.get(
                r["nearest_active_chembl"] or "")
        print(f"  {'#':>2}  {'compound':<30} {'appr':<5} {'Tc':>6} "
              f"{'deg':>5}  resembles")
        for n, r in enumerate(rows, 1):
            cid = r["chembl_id"] or r["candidate"].split(":", 1)[1][:14]
            shown = r["name"] or cid
            near = (r["nearest_active_name"] or r["nearest_active_chembl"]
                    or r["nearest_active"].split(":", 1)[1][:14])
            print(f"  {n:>2}  {shown[:30]:<30} "
                  f"{'YES' if r['is_approved'] else '-':<5}"
                  f" {r['tanimoto']:>6.3f} {r['graph_degree']:>5,}  "
                  f"{near[:28]}")
        if rows:
            n_appr = sum(1 for r in rows if r["is_approved"])
            n_reach = sum(1 for r in rows if r["graph_degree"])
            print(f"\n  {n_appr} of the top {len(rows)} are approved drugs; "
                  f"{n_reach} carry any non-label graph edge at all.")
            print(f"  Highest Tc {rows[0]['tanimoto']:.3f}, lowest shown "
                  f"{rows[-1]['tanimoto']:.3f}.")

    write_results(args.out, payload, args.release, schema.version)
    print(f"\n{'=' * 72}")
    print("A Tanimoto of 1.0 would mean the same ECFP4 bits, not the same")
    print("compound: the identity policy already collapses salts, charges")
    print("and tautomers, and `inchikey_skel` matches are excluded above, so")
    print("a near-1.0 row is a close analog rather than a duplicate -- which")
    print("is also the row least likely to be novel.")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
