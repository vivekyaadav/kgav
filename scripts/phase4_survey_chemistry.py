#!/usr/bin/env python3
"""Phase 4, step 0 — can the HCV retrospective be measured at all?

WHAT PHASE 4 IS FOR. The folds layer ships a catalytic-type gate whose whole
job is one refusal. Coronavirus 3CLpro and HCV NS3 share the trypsin-like
protease fold -- phase 0 measured it, and the shipped layer records 54
refused-despite-shared-fold pairs -- and differ in the catalytic residue,
cysteine against serine. kg_schema.yaml states the consequence of removing
the gate: M8 "traverses from SARS-CoV-2 Mpro to HCV NS3 -- same superfamily,
opposite nucleophile -- and recommends nirmatrelvir for hepatitis C".

NOTHING IN THIS REPOSITORY HAS TESTED WHETHER THAT REFUSAL IS RIGHT.
phase0_verify_folds.py and ingest_folds.py both say phases 5-10 must not be
built if the conservation premise is refuted, and both measure ANNOTATION --
which InterPro entries two chains share. Annotation cannot say whether a
compound that inhibits one protease inhibits the other. Only measured
chemistry can, and phase 4 is where it is asked.

THE TEST PHASE 4 WILL RUN, in three claims:

  A  ALLOWED arm.   CoV 3CLpro vs picornavirus 3C -- shared fold, cysteine
     both, gate PERMITS. The layer emits 20 such transfers and calls them
     strong. If compounds measured on both ends are unrelated in potency,
     the conservation premise fails on its own best case.
  B  REFUSED arm.   CoV 3CLpro vs HCV NS3 -- shared fold, cysteine against
     serine, gate BLOCKS. If compounds measured on both ends ARE related,
     the gate is deleting real signal.
  C  The contrast.  A must exceed B. If they are indistinguishable, the
     catalytic type predicts nothing and the gate is decoration.

THIS SCRIPT DOES NOT RUN THAT TEST. It asks the prior question -- does
ChEMBL hold enough compounds measured on BOTH ends of either arm for the
answer to mean anything -- and exits non-zero if it does not. A protocol
written after the counts are known, and committed before the code that
measures it, is what docs/prereg_m1_potency.md means by "a prediction made
afterwards is a description". The counts come first; the prediction second;
the measurement third.

WHY THIS READS ChEMBL AND NOT THE ASSEMBLED GRAPH. Two reasons, each
sufficient on its own:

  1  The graph does not contain the data. INHIBITS edges are built by
     ingest_chembl from the virus register, and chembl.query_activities
     filters on `t.tax_id IN (declared taxa)`. config/viruses.yaml declares
     8 viruses and 8 surveyed flaviviruses; HCV is in neither list, and
     neither is any picornavirus. Both arms of the test are invisible to
     the graph by construction -- which is also why 14 of the 24 resolved
     fold chains are comparators the spine does not build.
  2  The gate is the thing under test. Measuring a filter with the output
     of that filter answers nothing.

HOW TARGETS ARE FOUND, AND WHAT THAT MISSES. Every ChEMBL target is read and
classified in Python -- no target id is written down here, because an id
copied from memory is an assumption. Family comes from the organism STRING,
so a virus whose ChEMBL organism name contains none of the words in FAMILIES
is not seen at all; the distinct organism strings actually matched are
printed, with counts, so what was pulled in is visible rather than asserted.
Role comes from the target's pref_name, MOST SPECIFIC FIRST, and every
target in a matched organism is printed with the role it got or a dash.

TWO REJECTIONS ARE NOT INCIDENTAL AND ARE COUNTED SEPARATELY.

  POLYPROTEIN TARGETS. ChEMBL files a great deal of viral protease data
  against "Genome polyprotein". Counting those as protease measurements is
  the polyprotein problem this project has now hit at least three times --
  v1 chapter 9, the first phase 0 script comparing whole replicases, and
  the flavivirus NS3 rows that fell through to UniProtKB:<polyprotein>.
  They are rejected and reported, because how much data is lost that way is
  a number phase 4 needs.

  ORGANISM-LEVEL TARGETS. A target whose target_type is ORGANISM carries
  antiviral readouts -- replicon, CPE, cell culture -- not enzyme
  measurements. For HCV that is most of the chemistry in ChEMBL. A replicon
  EC50 says a compound stops the virus; it does not say it inhibits NS3.
  Rejected for the paired arms, and counted into the broad set below.

ON ASSAY READOUT. Arms A and B are claims about ENZYMES, so they take only
rows chembl.assay_readout calls "protein". A cellular antiviral number is
not a measurement of the protease it is filed against -- the same reasoning
that moved chloroquine's 160 nM off spike.

    python scripts/phase4_survey_chemistry.py
    python scripts/phase4_survey_chemistry.py --db path/to/chembl_37.db

Read-only. Writes nothing. Recommends nothing.

EXIT CODES, following phase0_verify_folds:
  0  measurable -- at least one arm clears the floor; go write the prereg
  1  not measurable -- the arms are too thin; phase 4 needs another design
  2  insufficient inputs -- no database, or a family with no targets at all
"""
from __future__ import annotations

import argparse
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.calibration import MIN_REACHED_POSITIVES
from kgav.chembl import QUANT_TYPES, assay_readout, to_nm
from kgav.schema import print_capped

# THE FLOOR, BORROWED RATHER THAN CHOSEN. calibration.MIN_REACHED_POSITIVES
# is 10 and its comment says why: below ten measured compounds an AUC
# describes the imputation rather than the channel, and 10 is already
# hard_negatives' minimum class size. An arm of phase 4 is a class in
# exactly that sense -- the compounds measured on both ends -- so it takes
# the same floor instead of inventing a third 10.
FLOOR = MIN_REACHED_POSITIVES

# Family is read from target_dictionary.organism. These are the words that
# appear in ChEMBL organism strings for the three families the folds layer
# compares; anything else is out of scope and is never fetched.
FAMILIES: dict[str, tuple[str, ...]] = {
    "Coronaviridae": (r"coronavirus", r"\bSARS\b", r"\bMERS\b", r"HCoV"),
    "Picornaviridae": (r"rhinovirus", r"enterovirus", r"poliovirus",
                       r"coxsackie", r"echovirus", r"hepatitis A virus",
                       r"picornavir", r"aphthovirus", r"cardiovirus"),
    "Flaviviridae": (r"hepatitis C virus", r"hepacivirus", r"dengue",
                     r"zika", r"west nile", r"yellow fever",
                     r"japanese encephalitis", r"tick-?borne encephalitis",
                     r"flavivir"),
}

# HCV is pulled out of Flaviviridae because arm B is specifically about it:
# it is the one serine protease in the comparison with a large, approved
# inhibitor class, and the one the schema names.
HCV = (r"hepatitis C", r"hepacivirus")

# Role from pref_name, MOST SPECIFIC FIRST. "3C" is a substring of
# "3C-like"; the project's own note calls the substring hazard its sixth
# instance. Each role is (name, require, reject, despite): the target takes
# the role when a require pattern matches and no reject pattern does -- and
# a reject is cancelled when a `despite` pattern is also present.
#
# Only 3clpro, 3cpro and ns3pro feed the arms. plpro is classified because
# it is one of the ten chains the folds layer ships and the one cysteine
# protease M8 cannot reach, so its chemistry is worth seeing beside them.
ROLES: dict[str, tuple[tuple[str, tuple[str, ...], tuple[str, ...],
                             tuple[str, ...]], ...]] = {
    "Coronaviridae": (
        ("3clpro", (r"3C-?like", r"\b3CL\w*", r"main protease", r"\bMpro\b",
                    r"\bnsp5\b"), (), ()),
        ("plpro", (r"papain-?like", r"\bPLpro\b", r"\bnsp3\b"), (), ()),
    ),
    "Picornaviridae": (
        ("3cpro", (r"picornain\s*3C\b", r"\b3C\s*protease", r"protease\s*3C\b",
                   r"\b3C\b"), (r"\b3CD\b",), ()),
    ),
    "Flaviviridae": (
        # NS3 is a protease AND a helicase in one chain -- the split
        # ingest_folds makes by domain coordinates. "NS3 helicase" is the
        # wrong enzyme. "NS3 protease/helicase" is the bifunctional
        # construct the protease data is actually filed against, so naming
        # the protease cancels the helicase rejection. Writing that as a
        # lookahead instead got it backwards: in "NS3 protease/helicase"
        # the word protease comes FIRST, so `helicase(?!.*protease)` fired
        # and dropped the main HCV target.
        ("ns3pro", (r"NS2B[-/ ]?NS3", r"NS3[-/ ]?4A", r"\bNS3\b"),
         (r"helicase", r"NTPase", r"polymerase", r"methyltransferase"),
         (r"protease", r"proteinase")),
    ),
}

# Rejected before role matching, and counted. See the module docstring.
POLYPROTEIN = re.compile(r"polyprotein", re.I)

# target_type values that denote a protein measurement. Anything else that
# matched a role is printed as excluded, with its type, rather than dropped
# silently -- for HCV the excluded set is most of the chemistry.
PROTEIN_TARGET_TYPES = frozenset({"SINGLE PROTEIN", "PROTEIN COMPLEX",
                                  "CHIMERIC PROTEIN"})

# Columns this script reads. Checked against the live schema before any
# query runs, so a ChEMBL release that renamed one fails with the column
# list in hand rather than as a SQL error three joins deep.
REQUIRED_COLUMNS = {
    "target_dictionary": ("tid", "chembl_id", "pref_name", "organism",
                          "tax_id", "target_type"),
    "activities": ("assay_id", "molregno", "standard_type", "standard_value",
                   "standard_units", "standard_relation"),
    "assays": ("assay_id", "tid", "description", "assay_type", "doc_id"),
    "molecule_dictionary": ("molregno", "chembl_id", "pref_name", "max_phase"),
}


def _any(patterns: tuple[str, ...], text: str) -> bool:
    return any(re.search(p, text or "", re.I) for p in patterns)


def classify_family(organism: str | None) -> str | None:
    """Which of the three compared families this organism belongs to.

    Read from the string, which is the limitation stated in the docstring:
    a coronavirus whose ChEMBL organism name says neither "coronavirus" nor
    SARS/MERS/HCoV is invisible here.
    """
    for family, patterns in FAMILIES.items():
        if _any(patterns, organism or ""):
            return family
    return None


def is_hcv(organism: str | None) -> bool:
    return _any(HCV, organism or "")


def classify_role(family: str, pref_name: str | None) -> str | None:
    """Protease role from the target's name, or None.

    A polyprotein target has no role: it names thirteen proteins at once,
    and a measurement filed against it does not say which one was inhibited.
    """
    name = pref_name or ""
    if POLYPROTEIN.search(name):
        return None
    for role, require, reject, despite in ROLES.get(family, ()):
        if not _any(require, name):
            continue
        if _any(reject, name) and not _any(despite, name):
            continue
        return role
    return None


def missing_columns(con: sqlite3.Connection) -> list[str]:
    """Every required column this database does not have, as table.column."""
    gone = []
    for table, cols in REQUIRED_COLUMNS.items():
        have = {r[1] for r in con.execute(f"PRAGMA table_info({table})")}
        if not have:
            gone.append(f"{table} (no such table)")
            continue
        gone += [f"{table}.{c}" for c in cols if c not in have]
    return gone


def load_targets(con: sqlite3.Connection) -> list[dict]:
    """Every target in a matched organism, with family, role and acceptance.

    Reads the whole table -- about fifteen thousand rows -- rather than
    filtering in SQL, so the family patterns live in one testable place
    instead of being half in Python and half in a LIKE clause.
    """
    out = []
    for r in con.execute("""SELECT tid, chembl_id, pref_name, organism,
                                   tax_id, target_type
                            FROM target_dictionary"""):
        tid, chembl_id, pref_name, organism, tax_id, ttype = r
        family = classify_family(organism)
        if family is None:
            continue
        role = classify_role(family, pref_name)
        out.append({"tid": tid, "chembl_id": chembl_id,
                    "pref_name": pref_name or "", "organism": organism or "",
                    "tax_id": tax_id, "target_type": ttype or "",
                    "family": family, "hcv": is_hcv(organism), "role": role,
                    "polyprotein": bool(POLYPROTEIN.search(pref_name or "")),
                    "protein_level": (ttype or "") in PROTEIN_TARGET_TYPES})
    return out


ACTIVITY_FIELDS = ("tid", "compound", "drug_name", "max_phase",
                   "standard_type", "standard_value", "standard_units",
                   "standard_relation", "description", "assay_type")


def load_activities(con: sqlite3.Connection, tids: list[int]) -> list[dict]:
    """Quantified potency rows for these targets, with the assay text.

    Fetched for EVERY matched target, not only the usable ones: the arms
    need protein-level enzyme rows, and the broad set needs the cellular
    coronavirus rows the arms reject. Filtering happens in Python, where
    both readings come from the same fetch and cannot disagree.

    AN EMPTY tids LIST IS REFUSED. SQLite accepts `IN ()` as always-false
    and returns zero rows, which is indistinguishable from these targets
    having no chemistry -- the failure chembl.query_activities documents
    costing a DENV-2 measurement of 8,748 rows a reading of zero.
    """
    if not tids:
        raise ValueError(
            "load_activities called with no targets: `IN ()` returns zero "
            "rows in SQLite rather than failing, so a caller bug would read "
            "as these proteases having no measured chemistry at all.")
    qmarks = ",".join("?" * len(tids))
    types = ",".join("?" * len(QUANT_TYPES))
    rows = con.execute(f"""
        SELECT a.tid, md.chembl_id, md.pref_name, md.max_phase,
               act.standard_type, act.standard_value, act.standard_units,
               act.standard_relation, a.description, a.assay_type
        FROM activities act
        JOIN assays a ON act.assay_id = a.assay_id
        JOIN molecule_dictionary md ON act.molregno = md.molregno
        WHERE a.tid IN ({qmarks})
          AND act.standard_value IS NOT NULL
          AND act.standard_type IN ({types})
    """, list(tids) + list(QUANT_TYPES)).fetchall()
    return [dict(zip(ACTIVITY_FIELDS, r)) for r in rows]


def compounds_by_role(acts: list[dict], targets: dict[int, dict],
                      *, protein_only: bool = True) -> dict[str, set[str]]:
    """role key -> the compounds measured on it, one entry per compound.

    The role key is the role name for every family but Flaviviridae, where
    HCV is split out as `hcv_ns3pro` because arm B is about HCV specifically
    and the other six NS3 proteases are a separate comparison.
    """
    out: dict[str, set[str]] = defaultdict(set)
    for a in acts:
        t = targets.get(a["tid"])
        if t is None or t["role"] is None or not t["protein_level"]:
            continue
        if to_nm(a["standard_value"], a["standard_units"]) is None:
            continue
        if protein_only:
            where, _ = assay_readout(a["description"] or "", a["assay_type"])
            if where != "protein":
                continue
        key = role_key(t)
        out[key].add(a["compound"])
    return out


def role_key(t: dict) -> str:
    if t["family"] == "Flaviviridae" and t["role"] == "ns3pro":
        return "hcv_ns3pro" if t["hcv"] else "flavi_ns3pro"
    return t["role"]


def broad_coronavirus_compounds(acts: list[dict], targets: dict[int, dict]) -> set[str]:
    """Compounds with ANY coronavirus measurement, enzyme or antiviral.

    Arms A and B need both ends measured on the enzyme. This set is the
    fallback the survey reports when they are thin: a compound with an HCV
    NS3 enzyme number and a coronavirus antiviral number is still the
    retrospective's real-world subject -- "the HCV protease inhibitors that
    were tested against SARS-CoV-2" -- it just answers a weaker question.
    """
    out: set[str] = set()
    for a in acts:
        t = targets.get(a["tid"])
        if t is not None and t["family"] == "Coronaviridae":
            out.add(a["compound"])
    return out


ARMS = (
    ("A  ALLOWED  CoV 3CLpro + picornavirus 3C", "3clpro", "3cpro",
     "gate PERMITS: shared fold, cysteine both"),
    ("B  REFUSED  CoV 3CLpro + HCV NS3", "3clpro", "hcv_ns3pro",
     "gate BLOCKS: shared fold, cysteine vs serine"),
    ("B' REFUSED  CoV 3CLpro + other flavivirus NS3", "3clpro", "flavi_ns3pro",
     "gate BLOCKS: shared fold, cysteine vs serine"),
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    args = ap.parse_args()

    if not args.db.exists():
        print(f"no ChEMBL database at {args.db}\n"
              f"Phase 4 reads the same file ingest_chembl reads; pass --db if "
              f"it lives elsewhere.\nExiting 2: nothing was measured.")
        return 2

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    gone = missing_columns(con)
    if gone:
        print("this database does not have the columns this script reads:")
        print_capped(gone, 20, indent="    ", unit="more column(s)")
        print("\nExiting 2: nothing was measured.")
        return 2

    targets = load_targets(con)
    by_tid = {t["tid"]: t for t in targets}
    print(f"{args.db}\n")

    # ---------------------------------------------------- what was matched
    print("organisms matched, by family — family is read from this string, "
          "so\na virus absent here was never fetched:")
    for family in FAMILIES:
        fam = [t for t in targets if t["family"] == family]
        orgs = Counter(t["organism"] for t in fam)
        print(f"  {family}  —  {len(orgs)} organism(s), {len(fam)} target(s)")
        print_capped(sorted(orgs.items(), key=lambda kv: (-kv[1], kv[0])),
                     8, lambda kv: f"{kv[0][:58]:<60} {kv[1]:>4}",
                     indent="      ", unit="more organism(s)")
    missing_family = [f for f in FAMILIES
                      if not any(t["family"] == f for t in targets)]

    # ------------------------------------------------- targets and roles
    print("\nprotease targets found, by family. Every target in a matched "
          "organism is\nlisted; a dash means no role pattern claimed it:")
    for family in FAMILIES:
        fam = sorted((t for t in targets if t["family"] == family),
                     key=lambda t: (t["role"] or "~", t["pref_name"]))
        roled = [t for t in fam if t["role"]]
        print(f"  {family}  —  {len(roled)} with a protease role, "
              f"{len(fam) - len(roled)} without")
        print_capped(roled, 14,
                     lambda t: (f"{t['chembl_id']:<14} {role_key(t):<13} "
                                f"{'' if t['protein_level'] else 'EXCLUDED '}"
                                f"{t['target_type']:<16} {t['pref_name'][:44]}"),
                     indent="      ", unit="more with a role")

    # ------------------------------------------------------- the chemistry
    usable = [t["tid"] for t in targets if t["role"] and t["protein_level"]]
    if not usable:
        print("\nno protein-level protease target in any of the three "
              "families.\nExiting 2: nothing was measured.")
        return 2
    acts = load_activities(con, [t["tid"] for t in targets])
    con.close()
    rows_per_tid = Counter(a["tid"] for a in acts)

    # WHAT WAS REJECTED, IN ROWS AS WELL AS TARGETS. An earlier draft of this
    # block counted only targets that had already taken a role, and printed
    # "rejected as NOT PROTEIN-LEVEL: 0" on a database holding an
    # ORGANISM-level HCV target -- true of the set it counted, and read as
    # nothing having been dropped. The two classes below are disjoint and
    # cover every matched target that cannot reach an arm.
    def _rows(ts):
        return sum(rows_per_tid[t["tid"]] for t in ts)

    poly = [t for t in targets if t["polyprotein"]]
    print(f"\n  rejected as POLYPROTEIN targets: {len(poly)} target(s), "
          f"{_rows(poly):,} quantified row(s)")
    print("  A measurement filed against a polyprotein does not say which of "
          "its\n  proteins was inhibited. This is the polyprotein problem, "
          "counted rather\n  than absorbed.")
    print_capped(sorted(poly, key=lambda t: -rows_per_tid[t["tid"]]), 8,
                 lambda t: f"{t['chembl_id']:<14} {t['organism'][:40]:<42} "
                           f"{rows_per_tid[t['tid']]:>7,} rows",
                 indent="      ", unit="more polyprotein target(s)")

    nonprotein = [t for t in targets
                  if not t["protein_level"] and not t["polyprotein"]]
    print(f"\n  rejected as NOT PROTEIN-LEVEL: {len(nonprotein)} target(s), "
          f"{_rows(nonprotein):,} quantified row(s)")
    print("  An ORGANISM target carries replicon and cell-culture antiviral "
          "numbers.\n  Those are evidence about the virus, not about the "
          "enzyme they sit beside,\n  so they cannot stand at either end of "
          "arms A or B. Counted whether or not\n  a role pattern would have "
          "claimed them, because the rejection is the\n  target's type and "
          "not its name.")
    print_capped(sorted(nonprotein, key=lambda t: -rows_per_tid[t["tid"]]), 8,
                 lambda t: f"{t['chembl_id']:<14} {t['target_type']:<16} "
                           f"{t['organism'][:28]:<30} "
                           f"{rows_per_tid[t['tid']]:>7,} rows",
                 indent="      ", unit="more excluded target(s)")

    protein = compounds_by_role(acts, by_tid, protein_only=True)
    anywhere = compounds_by_role(acts, by_tid, protein_only=False)
    print(f"\n{len(acts):,} quantified rows on those targets "
          f"({', '.join(QUANT_TYPES)})\n")
    print("compounds per role — protein-level readout is what the arms use; "
          "`any`\nincludes cellular rows filed against the same enzyme:")
    print(f"  {'role':<16} {'protein':>9} {'any readout':>13}")
    for key in sorted(set(protein) | set(anywhere)):
        print(f"  {key:<16} {len(protein.get(key, ())):>9,} "
              f"{len(anywhere.get(key, ())):>13,}")

    # ------------------------------------------------------------ the arms
    print("\nTHE ARMS — compounds measured at BOTH ends, on the enzyme:")
    print(f"  {'arm':<44} {'n':>6}   what the gate says")
    results = {}
    for label, left, right, says in ARMS:
        both = protein.get(left, set()) & protein.get(right, set())
        results[label] = both
        print(f"  {label:<44} {len(both):>6}   {says}")

    broad = broad_coronavirus_compounds(acts, by_tid)
    hcv_broad = protein.get("hcv_ns3pro", set()) & broad
    print(f"\n  the broad set: HCV NS3 enzyme inhibitors with ANY coronavirus "
          f"measurement\n  (enzyme or antiviral, any target_type): "
          f"{len(hcv_broad)}")
    print("  This is the retrospective as the literature ran it in 2020. It "
          "answers a\n  weaker question than arm B -- whether the compounds "
          "stopped the virus, not\n  whether they inhibited Mpro -- and is "
          "reported so a thin arm B has a\n  fallback rather than a dead end.")

    # --------------------------------------------------------- the verdict
    print()
    if missing_family:
        print(f"NOT MEASURABLE. No target at all for: "
              f"{', '.join(missing_family)}.\nA family with no targets is a "
              f"pattern problem or a ChEMBL gap, not a\nfinding about "
              f"chemistry.\nExiting 2.")
        return 2

    clears = {k: v for k, v in results.items() if len(v) >= FLOOR}
    if not clears:
        print(f"NOT MEASURABLE. No arm reaches {FLOOR} compounds measured at "
              f"both ends.")
        print(f"{FLOOR} is calibration.MIN_REACHED_POSITIVES, the floor this "
              f"project already\nuses for a class worth evaluating; below it "
              f"a figure describes its own\nsmall sample. Phase 4 needs "
              f"another design -- the broad set above is the\ncandidate, and "
              f"it tests a weaker claim.")
        print("\nExiting 1 so this is not mistaken for a refutation of the "
              "premise.\nAn unmeasurable test is not a failed test.")
        return 1

    print(f"MEASURABLE. {len(clears)} arm(s) reach the {FLOOR}-compound floor:")
    for label in clears:
        print(f"  {label}  n={len(clears[label])}")
    print("\nNext: docs/prereg_phase4_hcv.md, written against these counts "
          "and committed\nBEFORE the code that measures the arms. The "
          "prediction the premise makes is\nthat arm A relates and arm B "
          "does not; writing it afterwards would make it a\ndescription.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
