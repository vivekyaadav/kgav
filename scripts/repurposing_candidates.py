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

#: The registered runs this read-out quotes its own trustworthiness from.
#: The strength-matched control is the PRIMARY because it is the sternest
#: one that ran -- it removes each candidate's k_cell nearest references --
#: so quoting it is the conservative choice. The other two are carried
#: beside it, because a reader deciding whether to chase a row deserves the
#: range rather than the number I picked.
SCORER_RUNS = (
    ("strength_matched", "data/results/structure_scorer_strength_matched.json"),
    ("symmetric_k10", "data/results/structure_scorer_symmetric.json"),
    ("publication", "data/results/structure_scorer.json"),
)
PRIMARY_RUN = SCORER_RUNS[0][0]


def verdicts(root: Path = ROOT) -> dict[str, dict]:
    """virus curie -> what the registered runs say about THAT cell.

    WHY THIS EXISTS. Until 2026-11-12 this read-out emitted an
    identically-shaped ranked list for every virus with ten actives, and
    the file carried `n_known_actives` and `n_candidates` and NO measure of
    whether the scorer works on that virus. So Yellow fever at AUC 0.8117,
    which clears its corrected threshold, and DENV-4 at **0.4792**, where
    the scorer ANTI-RANKS, printed the same way. A list from a cell below
    chance is worse than random and there was nothing in the output to say
    so -- the header's "clears in ten of fourteen" is a statement about the
    FAMILY, and a reader wants the statement about the row in front of them.

    Each cell gets its primary AUC, its corrected interval, whether it
    clears, the direction, its k_cell, and the other two controls' AUCs.
    """
    out: dict[str, dict] = {}
    for tag, rel in SCORER_RUNS:
        path = root / rel
        if not path.exists():
            continue
        doc = json.loads(path.read_text())
        for virus, cell in doc.get("primary_chembl", {}).items():
            m = cell.get("structure_max") or {}
            if m.get("auc") is None:
                continue
            v = out.setdefault(virus, {})
            v[f"auc_{tag}"] = round(m["auc"], 4)
            if tag != PRIMARY_RUN:
                continue
            v.update({
                "auc": round(m["auc"], 4),
                "auc_ci_lo": round(m["auc_ci_lo_corrected"], 4),
                "auc_ci_hi": round(m["auc_ci_hi_corrected"], 4),
                "clears_corrected_alpha": bool(
                    m.get("clears_corrected_alpha")),
                "direction": m.get("direction"),
                "k_cell": cell.get("k_cell"),
                "control": PRIMARY_RUN,
            })
    return out


def tier(v: dict | None) -> str:
    """What a cell's verdict licenses: 'clears', 'weak', 'below', 'unknown'.

    THREE TIERS AND NOT TWO, because "does not clear" covers two situations
    a reader must not confuse. A cell above chance that misses the
    corrected threshold is a weak suggestion; a cell BELOW chance means the
    scorer ranks known actives *under* known inactives there, so its
    ranking of unlabelled compounds is actively misleading and no list
    should be printed from it at all.
    """
    if not v or v.get("auc") is None:
        return "unknown"
    if v["auc"] < 0.5 or v.get("direction") == "below":
        return "below"
    return "clears" if v.get("clears_corrected_alpha") else "weak"


def headline(vs: dict[str, dict], root: Path = ROOT) -> list[str]:
    """The 'what this is worth' banner, DERIVED.

    Every figure in it was typed by hand until 2026-11-12 -- 0.7841 for a
    cell that now reads 0.7822, "none below" for a family that now has two
    below chance -- in the one paragraph whose job is to stop a reader
    over-reading the list.
    """
    lines = []
    sars = vs.get("NCBITaxon:2697049", {})
    clears = sum(1 for v in vs.values() if tier(v) == "clears")
    below = [k for k, v in vs.items() if tier(v) == "below"]
    weak = [k for k, v in vs.items() if tier(v) == "weak"]
    if sars.get("auc") is not None:
        lines.append(f"  On labelled compounds this scorer reads AUC "
                     f"{sars['auc']:.4f} on SARS-CoV-2")
        lines.append(f"  and clears a corrected threshold above chance in "
                     f"{clears} of {len(vs)} viruses,")
        lines.append(f"  with {len(weak)} above chance but not clearing and "
                     f"{len(below)} BELOW chance -- see")
        lines.append("  docs/prereg_structure_scorer_strength_matched.md.")
    path = root / "data/results/structure_scorer_strength_matched.json"
    if path.exists():
        doc = json.loads(path.read_text())
        mm = (doc.get("secondary_ncats", {})
              .get("NCBITaxon:2697049", {}).get("structure_max") or {})
        knn = root / "data/results/structure_scorer_symmetric.json"
        gap10 = None
        if knn.exists():
            kd = json.loads(knn.read_text())
            a = (kd["primary_chembl"]["NCBITaxon:2697049"]
                 ["structure_max"]["auc"])
            b = (kd["secondary_ncats"]["NCBITaxon:2697049"]
                 ["structure_max"]["auc"])
            gap10 = a - b
        if mm.get("auc") is not None and sars.get("auc") is not None:
            gap = sars["auc"] - mm["auc"]
            lo, hi = mm["auc_ci_lo_corrected"], mm["auc_ci_hi_corrected"]
            contains = "containing chance" if lo < 0.5 < hi else \
                       "excluding chance"
            lines.append("  But the one population matched by assay "
                         "rather than by source")
            lines.append(f"  is MARGINAL ({mm['auc']:.4f}, interval "
                         f"{contains}) and the gap")
            lines.append(f"  between the two is "
                         f"{min(gap, gap10 or gap):.2f}-"
                         f"{max(gap, gap10 or gap):.2f} AUC.")
    lines.append("  So every row below is a chemical-analogy suggestion "
                 "with a")
    lines.append("  measured hit rate. It is NOT a prediction of activity "
                 "and NOT a")
    lines.append("  mechanism: the evidence is 'resembles a known active', "
                 "and the")
    lines.append("  known active is named so you can judge.")
    return lines


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


def _heavy_delta(a: str | None, b: str | None) -> int | None:
    """Heavy atoms in `a` minus heavy atoms in `b`, or None if either is
    unparseable. A ring homolog reads 1; a residue permutation reads 0, and
    that zero is informative rather than missing -- it says the two differ
    in a way neither the fingerprint NOR the atom count can see."""
    from rdkit import Chem

    try:
        ma = Chem.MolFromSmiles(a) if a else None
        mb = Chem.MolFromSmiles(b) if b else None
    except Exception:
        return None
    if ma is None or mb is None:
        return None
    return ma.GetNumHeavyAtoms() - mb.GetNumHeavyAtoms()


def suggest(candidates: list[str], actives: list[str], fps: dict,
            props: dict, top: int, degree: dict | None = None) -> list[dict]:
    """Top `top` candidates by nearest known active.

    The two exclusions that survive from the registered protocol and mean
    something here: a candidate is never scored against itself (it is not in
    `actives` by construction), and a reference sharing its `inchikey_skel`
    is dropped -- a stereoisomer or tautomer of a known active is the same
    compound in this graph's identity policy, so scoring against it would
    report a Tanimoto of 1.0 and suggest a compound already known.

    AND THE SKELETON EXCLUSION DOES NOT CATCH EVERY TANIMOTO OF 1.0, which
    reading the output showed on 2026-11-12: 14 of 210 suggestions sat at
    **1.0** with a DIFFERENT skeleton from the active they resemble -- so not
    salts, which L27 collapsed, and not stereoisomers, which the line above
    removes. They are distinct molecules with an identical ECFP4 bit vector,
    because radius 2 folded into 2048 bits cannot see:

        a saturated RING SIZE   cycloheptyl vs cyclohexyl (CHEMBL4639263 /
                                CHEMBL4526128) -- every atom's radius-2
                                environment is CH2 flanked by CH2 either way
        a peptide's LENGTH      one extra leucine (CHEMBL4857345 /
                                CHEMBL1821856)
        a residue ORDER         Arg-Arg-Lys vs Arg-Lys-Lys (CHEMBL468049 /
                                CHEMBL522355)

    A suggestion at 1.0 therefore reads "this is the same molecule as a known
    active" when it may be one methylene away, which is a different and much
    weaker statement. `fingerprint_identical` and `heavy_atom_delta` now say
    which it is, so the top of the list cannot be read as an identity it is
    not. **The registered AUC is not exposed to this:** both symmetric
    controls drop each candidate's k NEAREST references, and a near-identity
    is by definition among them -- measured at zero candidates above 0.90
    under either, against 0.0033 of AUC under the publication control, by
    `scripts/diagnose_fingerprint_blindspots.py`. This read-out applies no
    k-nearest exclusion, which is why it is the thing that shows the defect.
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
                     # WHAT A 1.0 ACTUALLY MEANS HERE. Identical bits with a
                     # different skeleton is a fingerprint limit, not an
                     # identity, and the heavy-atom difference is the
                     # cheapest honest measure of how far apart they are:
                     # 1 for a ring homolog, 0 for a residue permutation.
                     "fingerprint_identical": best >= 0.9999,
                     "heavy_atom_delta": _heavy_delta(
                         props.get(c, {}).get("smiles"),
                         props.get(ref_ids[best_i], {}).get("smiles")),
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

    # THE VERDICTS FIRST, AND NO LIST WITHOUT THEM. A ranked list whose
    # trustworthiness cannot be stated is the thing this script was
    # criticised for printing, so an absent registered run is a REFUSAL
    # rather than a list with a missing column.
    vs = verdicts(ROOT)
    if not vs:
        print("REFUSED: no registered scorer run on disk, so no cell's "
              "AUC can be\nquoted beside its suggestions -- and a ranked "
              "list with no measure of\nwhether the scorer works on that "
              "virus is what this script exists not\nto print. Run "
              "`make structure-scorer-strength` first.")
        return 2

    print("REPURPOSING CANDIDATES — the structure scorer, pointed at the")
    print("question rather than at its own validation.\n")
    print("HOW MUCH THIS IS WORTH, BEFORE THE NAMES:")
    for line in headline(vs, ROOT):
        print(line)
    print()

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
        # THIS CELL'S OWN VERDICT DECIDES WHETHER THERE IS A LIST. Below
        # chance the scorer ranks known actives UNDER known inactives, so
        # its ranking of unlabelled compounds is actively misleading and
        # the suggestions are refused rather than printed with a caveat --
        # a caveat beside fifteen named drugs loses to the drugs.
        v = vs.get(virus)
        t_ = tier(v)
        # THE PAYLOAD AND THE PRINTOUT MUST AGREE. The first version set
        # `rows = []` only for "below" and then skipped printing for
        # "below" OR "unknown" -- so an unknown cell wrote fifteen
        # suggestions into the JSON under a banner saying NO LIST IS
        # PRINTED. A file disagreeing with the run that wrote it is worse
        # than either being wrong.
        if t_ in ("below", "unknown"):
            rows = []
        else:
            rows = suggest(candidates, actives, fps, props, args.top,
                           degree=g.degree_unlabelled)
        payload["viruses"][virus] = {
            "name": names.get(virus),
            "n_known_actives": len(actives),
            "n_known_inactives": len(inactives),
            "n_candidates": len(candidates),
            "n_salts_of_known_compounds_dropped": salts,
            # WHAT THE REGISTERED RUNS SAY ABOUT THIS CELL, carried in the
            # file beside the list so the two cannot be read apart.
            "scorer_verdict": dict(v or {}, tier=t_),
            "suggestions": rows,
        }
        label = names.get(virus, virus)
        print(f"{'=' * 72}")
        print(f"{label}  ({virus})")
        print(f"  {len(actives):,} known actives, {len(inactives):,} known "
              f"inactives excluded, {len(candidates):,} candidates")
        print(f"  {salts:,} dropped as the same compound as something "
              f"already known, in either direction (fix-list L27)")
        if v and v.get("auc") is not None:
            others = "  ".join(
                f"{tag.replace('_', ' ')} {v['auc_' + tag]:.4f}"
                for tag, _ in SCORER_RUNS
                if 'auc_' + tag in v and tag != PRIMARY_RUN)
            k = v.get("k_cell")
            kbit = f" (k_cell={k})" if k else ""
            print(f"  THIS CELL: AUC {v['auc']:.4f} "
                  f"[{v['auc_ci_lo']:.4f}, {v['auc_ci_hi']:.4f}] under the "
                  f"strength-matched control{kbit}")
            if others:
                print(f"             other controls: {others}")
            # WHEN THE CONTROLS DISAGREE ABOUT THE DIRECTION, SAY SO. DENV-4
            # reads 0.4792 under the strength-matched control and 0.5868
            # under the symmetric one, so "below chance" is a statement
            # about the sternest control and not about the chemistry. The
            # refusal still stands -- the conservative control decides --
            # but a reader who is told only the refusal would not know the
            # question is open.
            sides = {tag for tag, _ in SCORER_RUNS
                     if 'auc_' + tag in v
                     and (v['auc_' + tag] >= 0.5) != (v['auc'] >= 0.5)}
            if sides:
                print(f"             ! THE CONTROLS DISAGREE about which "
                      f"side of chance this cell is on\n"
                      f"               ({', '.join(sorted(sides))} put it "
                      f"the other way). The sternest one decides here.")
        verdict_line = {
            "clears": "  VERDICT: clears the corrected alpha above chance. "
                      "The list below is usable.",
            "weak": "  VERDICT: ABOVE chance but does NOT clear the "
                    "corrected alpha.\n           Treat the list as weak: "
                    "the cell is readable and not significant.",
            "below": "  VERDICT: BELOW CHANCE. The scorer ranks this "
                     "virus's known actives\n           UNDER its known "
                     "inactives, so any ranking of unlabelled\n"
                     "           compounds from it is worse than random. "
                     "NO LIST IS PRINTED.",
            "unknown": "  VERDICT: no registered run covers this cell, so "
                       "nothing is known\n           about whether the "
                       "scorer works here. NO LIST IS PRINTED.",
        }[t_]
        print(verdict_line)
        print(f"{'=' * 72}")
        if t_ in ("below", "unknown"):
            continue
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
            # A 1.0 IS MARKED, WITH WHAT ACTUALLY DIFFERS. Unmarked it
            # reads as an identity, and 14 of these are a ring homolog or a
            # peptide variant that ECFP4 at radius 2 cannot resolve.
            flag = ""
            if r.get("fingerprint_identical"):
                d = r.get("heavy_atom_delta")
                flag = ("  <- SAME BITS, "
                        + (f"{abs(d)} heavy atom(s) apart"
                           if d else "same atom count: a permutation")
                        + "; NOT the same compound")
            print(f"  {n:>2}  {shown[:30]:<30} "
                  f"{'YES' if r['is_approved'] else '-':<5}"
                  f" {r['tanimoto']:>6.3f} {r['graph_degree']:>5,}  "
                  f"{near[:28]}{flag}")
        if rows:
            n_appr = sum(1 for r in rows if r["is_approved"])
            n_reach = sum(1 for r in rows if r["graph_degree"])
            print(f"\n  {n_appr} of the top {len(rows)} are approved drugs; "
                  f"{n_reach} carry any non-label graph edge at all.")
            n_same = sum(1 for r in rows if r.get("fingerprint_identical"))
            print(f"  Highest Tc {rows[0]['tanimoto']:.3f}, lowest shown "
                  f"{rows[-1]['tanimoto']:.3f}.")
            if n_same:
                print(f"  {n_same} of them share the active's ECFP4 bits "
                      f"EXACTLY while being a different structure -- read "
                      f"those as analogs, not identities.")

    # THE CENSUS, so the run's own shape is readable without scrolling.
    seen = [(r["name"] or k, r["scorer_verdict"]["tier"], len(r["suggestions"]))
            for k, r in payload["viruses"].items()]
    counts = {name: [n for n, tr, _ in seen if tr == name]
              for name in ("clears", "weak", "below", "unknown")}
    print(f"{'=' * 72}")
    print("WHICH VIRUSES THIS RUN CAN AND CANNOT ANSWER FOR")
    print(f"  usable list (clears the corrected alpha): "
          f"{len(counts['clears'])}")
    for name in sorted(counts["clears"]):
        print(f"      {name}")
    for tier_name, why in (
            ("weak", "above chance, does NOT clear -- weak"),
            ("below", "BELOW chance -- refused, no list"),
            ("unknown", "no registered run -- refused, no list")):
        if counts[tier_name]:
            print(f"  {why}: {len(counts[tier_name])}")
            for name in sorted(counts[tier_name]):
                print(f"      {name}")
    # AND THE VIRUSES THAT NEVER REACHED THE LOOP, which is the half of
    # the question a census of the output cannot answer: a virus with fewer
    # than MIN_ACTIVES known actives is skipped before any scoring, so it
    # appears nowhere above and its absence looks like an oversight.
    # HCoV-NL63 has 5 actives and 7 inactives; HCoV-HKU1 has none of
    # either in any source. Named with their counts, because "needs five
    # more actives" and "has no antiviral data at all" are different
    # problems with different remedies.
    absent = sorted(set(names) - set(payload["viruses"]))
    if absent and not args.virus:
        print(f"  not scored at all (fewer than {MIN_ACTIVES} known "
              f"actives): {len(absent)}")
        for v in absent:
            print(f"      {names.get(v, v):<14} "
                  f"{len(labels.positives.get(v, set())):>3} actives, "
                  f"{len(labels.negatives.get(v, set())):>3} measured "
                  f"inactives")
    print(f"{'=' * 72}")

    write_results(args.out, payload, args.release, schema.version)
    print(f"\n{'=' * 72}")
    print("A Tanimoto of 1.0 means the same ECFP4 bits, NOT the same")
    print("compound, and on 2026-11-12 that stopped being a caveat and")
    print("became a measurement: the identity policy collapses salts,")
    print("charges and tautomers, `inchikey_skel` matches are excluded")
    print("above, and 14 rows STILL read 1.0 against a different skeleton.")
    print("Radius 2 folded into 2048 bits cannot see a saturated ring's")
    print("size (cycloheptyl vs cyclohexyl), a peptide's length (one extra")
    print("leucine) or a residue's position (Arg-Arg-Lys vs Arg-Lys-Lys).")
    print("Those rows are marked SAME BITS with the heavy-atom difference,")
    print("so a 1.0 cannot be read as an identity it is not.")
    print("The registered AUC is NOT exposed to this: both symmetric")
    print("controls drop each candidate's k nearest references and a")
    print("near-identity is among them -- zero candidates above 0.90 under")
    print("either, against 0.0033 of AUC under the publication control.")
    print("Measured by scripts/diagnose_fingerprint_blindspots.py.")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
