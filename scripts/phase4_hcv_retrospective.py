#!/usr/bin/env python3
"""Phase 4: the broad set, measured exactly as docs/prereg_phase4_hcv.md
registered it.

THE PREREG IS THE SPECIFICATION AND THIS SCRIPT HAS NO FREEDOM. Every
definition below is quoted from it, and where the registered wording did not
anticipate a case, the wording decides the PRIMARY result and the
alternative is reported beside it as a declared sensitivity. A protocol
improved after seeing the data is not the protocol that was registered.

WHAT IS BEING ASKED. The paired arms are refused: arm A is twelve compounds
over four independent observations, arm B is four. What remains is the
retrospective as the 2020 literature actually ran it -- of the compounds
known to inhibit HCV NS3, how many stopped a coronavirus, by any mechanism?
Seventeen compounds over seventeen independent observations, every pair
cross-paper, because the HCV enzyme literature and the coronavirus
literature do not overlap.

WHAT IT CANNOT SETTLE, in the prereg's own words: "This test can refute the
gate. It cannot confirm it." The gate predicts a null. A significant
enrichment says the catalytic-type gate is deleting real signal; a
non-significant result says nothing was detected at this n, which is not the
same as nothing being there. The registered power table says the test sees
roughly a three-fold enrichment and nothing finer.

    python scripts/phase4_hcv_retrospective.py

Read-only. Writes nothing: the verdict belongs beneath the prereg, in the
shape docs/prereg_temporal_auc.md uses, not in a results file stamped
against a release this measurement does not read.

EXIT CODES
  0  the test ran; read the verdict
  1  the set is below the registered floor, so it did not run
  2  insufficient inputs
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
from collections import Counter
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.calibration import MIN_REACHED_POSITIVES  # noqa: E402
from kgav.chembl import QUANT_TYPES, to_nm  # noqa: E402
from kgav.labels import classify  # noqa: E402

# ONE READER. The organism patterns, the role patterns, the polyprotein
# description rule, the potency rule and the identity rule all live in the
# survey that produced the registered counts. A second copy here would be a
# second definition of the set the prereg is about.
from phase4_survey_chemistry import (  # noqa: E402
    IDENTITIES,
    broad_coronavirus_compounds,
    broad_publications,
    compound_id,
    independent_units,
    load_activities,
    load_targets,
    missing_columns,
    role_publications,
    usable_potency,
)

ALPHA = 0.05


def coronavirus_edge(a: dict) -> dict | None:
    """One activity row in the shape labels.classify reads.

    THE PREREG SAYS THE IMPLEMENTATION CALLS THAT FUNCTION AND DOES NOT
    RESTATE THE RULE. So the 10,000 nM threshold, and what a ">" or a "<"
    proves, appear nowhere in this file: classify owns both, with its
    declared defaults, and this function only translates a ChEMBL row into
    the dict it expects.
    """
    field = QUANT_TYPES.get(a["standard_type"])
    nm = to_nm(a["standard_value"], a["standard_units"])
    if field is None or nm is None:
        return None
    return {"qualifiers": {field: nm,
                           "relation": a["standard_relation"] or "="}}


def coronavirus_verdicts(acts: list[dict], targets: dict[int, dict],
                         *, skeleton: bool = False
                         ) -> tuple[dict[str, set[str]], Counter]:
    """compound -> the set of verdicts its coronavirus rows produced."""
    out: dict[str, set[str]] = {}
    stats: Counter = Counter()
    for a in acts:
        t = targets.get(a["tid"])
        if t is None or t["family"] != "Coronaviridae":
            continue
        cid = compound_id(a, skeleton=skeleton)
        if cid is None:
            stats["no_inchikey"] += 1
            continue
        edge = coronavirus_edge(a)
        if edge is None:
            stats["unconvertible"] += 1
            continue
        verdict = classify(edge)
        stats[f"row_{verdict}"] += 1
        out.setdefault(cid, set()).add(verdict)
    return out, stats


def fisher_one_sided(a: int, b: int, c: int, d: int, *,
                     greater: bool = True) -> float:
    """Exact hypergeometric tail for the 2x2 [[a, b], [c, d]].

    BY HAND, NOT scipy. pyproject declares pyyaml, pandas and pyarrow;
    labels.auc_interval implements Hanley & McNeil the same way and for the
    same reason. A statistical dependency taken for one p-value is a
    dependency the whole project then carries.

    `greater` is the registered alternative -- enrichment of the first row.
    The other tail is computed only to satisfy P2, which requires a
    depletion to be reported and claimed as evidence for nothing.
    """
    n1, n2, k, n = a + b, c + d, a + c, a + b + c + d
    if n1 == 0 or n2 == 0 or k == 0 or k == n:
        return 1.0
    lo, hi = max(0, k - n2), min(n1, k)
    rng = range(a, hi + 1) if greater else range(lo, a + 1)
    return sum(comb(n1, i) * comb(n2, k - i) for i in rng) / comb(n, k)


def two_by_two(hcv: set[str], verdicts: dict[str, set[str]],
               *, undecidable_is_inactive: bool) -> tuple[int, int, int, int, int]:
    """(a, b, c, d, excluded) for the registered 2x2.

    THE PREREG DEFINES "ACTIVE" AND NOT "INACTIVE". Its words are: "A
    compound is active if classify returns 'active' for at least one
    coronavirus row." Read literally, everything else is not-active,
    including a compound whose every row was undecidable -- and an
    undecidable row is not evidence of inactivity. The literal reading is
    the PRIMARY because it is what was registered; the three-way split is
    reported beside it, and the gap between them is the size of a case the
    registered wording did not consider.
    """
    a = b = c = d = excluded = 0
    for cid, vs in verdicts.items():
        active = "active" in vs
        decided = active or ("inactive" in vs)
        if not decided and not undecidable_is_inactive:
            excluded += 1
            continue
        if cid in hcv:
            a, b = (a + 1, b) if active else (a, b + 1)
        else:
            c, d = (c + 1, d) if active else (c, d + 1)
    return a, b, c, d, excluded


def report(label: str, hcv: set[str], verdicts: dict[str, set[str]],
           *, undecidable_is_inactive: bool) -> tuple[float, float, int]:
    a, b, c, d, excluded = two_by_two(
        hcv, verdicts, undecidable_is_inactive=undecidable_is_inactive)
    n1, n2 = a + b, c + d
    p_up = fisher_one_sided(a, b, c, d, greater=True)
    p_dn = fisher_one_sided(a, b, c, d, greater=False)
    print(f"\n  {label}")
    print(f"    {'':<22} {'active':>8} {'inactive':>9} {'n':>7} {'rate':>8}")
    print(f"    {'HCV NS3 inhibitors':<22} {a:>8,} {b:>9,} {n1:>7,} "
          f"{(a / n1 if n1 else 0):>8.1%}")
    print(f"    {'background':<22} {c:>8,} {d:>9,} {n2:>7,} "
          f"{(c / n2 if n2 else 0):>8.1%}")
    if excluded:
        print(f"    {excluded:,} compound(s) excluded: every coronavirus "
              f"row undecidable")
    print(f"    one-sided Fisher, enrichment (REGISTERED): p = {p_up:.4f}")
    print(f"    one-sided Fisher, depletion  (P2 only):    p = {p_dn:.4f}")
    return p_up, p_dn, n1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    args = ap.parse_args()

    if not args.db.exists():
        print(f"no ChEMBL database at {args.db}\nExiting 2.")
        return 2
    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    gone = missing_columns(con)
    if gone:
        print(f"missing columns: {', '.join(gone[:6])}\nExiting 2.")
        return 2
    targets = load_targets(con)
    acts = load_activities(con, [t["tid"] for t in targets])
    con.close()
    by_tid = {t["tid"]: t for t in targets}
    print(f"{args.db}\n")
    print("Protocol: docs/prereg_phase4_hcv.md, registered 2026-10-16 before "
          "this script\nexisted. Nothing below departs from it.\n")

    verdict_rows: list[tuple] = []
    for ident, skel in IDENTITIES:
        pot, _, _ = usable_potency(acts, by_tid, skeleton=skel)
        broad = broad_coronavirus_compounds(acts, by_tid, skeleton=skel)
        hcv = set(pot.get("hcv_ns3pro", {})) & broad
        pubs = {"hcv_ns3pro": role_publications(acts, by_tid, pot,
                                                skeleton=skel)
                .get("hcv_ns3pro", {}),
                "_coronavirus": broad_publications(acts, by_tid,
                                                   skeleton=skel)}
        units = independent_units(hcv, "hcv_ns3pro", "_coronavirus", pubs)
        verdicts, stats = coronavirus_verdicts(acts, by_tid, skeleton=skel)

        print(f"=== {ident}")
        print(f"  the set: {len(hcv)} compound(s), {units} independent "
              f"observation(s)")
        print("  coronavirus rows classified: "
              + ", ".join(f"{k.replace('row_', '')} {v:,}"
                          for k, v in sorted(stats.items())))
        if units < MIN_REACHED_POSITIVES:
            print(f"  BELOW THE REGISTERED FLOOR of "
                  f"{MIN_REACHED_POSITIVES}. Not run at this identity.")
            continue
        p_up, p_dn, n1 = report(
            "PRIMARY — the prereg's literal wording (undecidable counts as "
            "not-active)", hcv, verdicts, undecidable_is_inactive=True)
        report("SENSITIVITY — undecidable excluded, a case the registered "
               "wording did not consider", hcv, verdicts,
               undecidable_is_inactive=False)
        verdict_rows.append((ident, p_up, p_dn, n1))

    # ----------------------------------------------------- the decision rule
    print("\n" + "=" * 70)
    if not verdict_rows:
        print("NOT RUN. No identity reached the registered floor.\nExiting 1.")
        return 1
    ident, p_up, p_dn, n1 = verdict_rows[0]
    print(f"DECISION, per the rule registered in the prereg. Primary "
          f"identity:\n{ident}, n = {n1}.\n")
    if p_up < ALPHA:
        print(f"P1 REJECTED (p = {p_up:.4f} < {ALPHA}). HCV NS3 inhibitors "
              f"ARE enriched for\ncoronavirus activity. The conservation "
              f"premise's nucleophile rule is refuted\nat the level this "
              f"evidence reaches: the catalytic-type gate is removing "
              f"signal,\nand phases 5-10 must not be built on it as it "
              f"stands. This is the outcome\ningest_folds names.")
    elif p_dn < ALPHA:
        print(f"NO RESULT, and a significant DEPLETION (p = {p_dn:.4f}).\n"
              f"Per P2 this is reported and claimed as evidence for nothing: "
              f"HCV NS3\ninhibitors are a narrow chemotype of large "
              f"macrocyclic peptidomimetics, and\nunder-performing a "
              f"repurposing library is what that chemotype predicts\n"
              f"regardless of any fold or nucleophile.")
    else:
        print(f"NO RESULT (p = {p_up:.4f}). Nothing is confirmed. The "
              f"catalytic-type gate\nremains an UNTESTED safety measure, and "
              f"that wording belongs wherever the\nproject describes it.")
    print("\nIn every branch the paired arms stay refused. The registered "
          "power table\nsays this test sees roughly a three-fold enrichment "
          "and nothing finer, so a\nnull here is not a tight bound.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
