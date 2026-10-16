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

POLYPROTEIN TARGETS ARE READ, NOT REJECTED -- AND THE FIRST VERSION OF THIS
SCRIPT GOT THAT WRONG. ChEMBL files most viral protease data against the
polyprotein target and names the enzyme in the ASSAY DESCRIPTION. Rejecting
those targets by name, which is what the first run did, threw out 31
targets and 17,102 rows including CHEMBL4523582's 6,690, reported
"Coronaviridae: 0 with a protease role", and returned NOT MEASURABLE. That
verdict was an artifact of the filter, not a fact about chemistry.

This project already knows how these rows are filed: ingest_activities
calls chembl.resolve_nsp on the description and stamps the resulting edge
target_resolution="chain_from_assay_description". The same tables resolve
the role here, per ROW. What they do not resolve is printed verbatim, in
full, and that listing is what the tables get written from: it is how
chembl.PICORNA_PATTERNS came to exist, out of all 19 descriptions on the
six picornavirus polyprotein targets, rather than out of what the proteins
are called in a textbook. The flavivirus table was written the same way.

ORGANISM-LEVEL TARGETS ARE REJECTED. A target whose target_type is ORGANISM
carries antiviral readouts -- replicon, CPE, cell culture -- not enzyme
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
import math
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from statistics import median
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.calibration import MIN_REACHED_POSITIVES
from kgav.chembl import (PATTERNS_BY_FAMILY, QUANT_TYPES, assay_readout,
                         resolve_nsp, to_nm)
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

# A POLYPROTEIN TARGET IS WHERE THE DATA IS, NOT WHERE IT IS LOST.
#
# The first version of this script rejected every polyprotein-named target
# outright, on the reasoning that a measurement filed against a polyprotein
# does not say which of its proteins was inhibited. Run on chembl_37 that
# threw out 31 targets and 17,102 rows -- among them CHEMBL4523582 with
# 6,690 -- and reported "Coronaviridae: 0 with a protease role". The verdict
# was NOT MEASURABLE, and it was an artifact of my own filter.
#
# ChEMBL files coronavirus and flavivirus enzyme data against the
# polyprotein target and says which enzyme in the ASSAY DESCRIPTION. This
# project already knows that: chembl.ingest_activities calls resolve_nsp on
# the description and stamps the edge target_resolution=
# "chain_from_assay_description". So the role is resolved per ROW here too,
# from the same tables, and a polyprotein target is a target whose rows need
# reading rather than one whose rows are gone.
POLYPROTEIN = re.compile(r"polyprotein", re.I)

# (family, nsp, domain) from resolve_nsp -> the arm role it supplies.
POLYPROTEIN_ROLE = {
    ("Coronaviridae", "nsp5", None): "3clpro",
    ("Coronaviridae", "nsp3", "PLpro"): "plpro",
    ("Flaviviridae", "NS3", "protease"): "ns3pro",
    ("Picornaviridae", "3C", None): "3cpro",
}

# WHICH FAMILIES HAVE A MEASURED DESCRIPTION TABLE.
#
# DERIVED FROM chembl.PATTERNS_BY_FAMILY, not listed again here. The earlier
# version wrote the set out by hand, which was correct on the day and is a
# second copy of a fact the library already holds: adding Picornaviridae to
# the library while this set still said {Coronaviridae, Flaviviridae} would
# have left the new table unreachable, with the survey still reporting "NO
# DESCRIPTION TABLE EXISTS" beside a table that does.
#
# Coronaviridae is unioned in because it has no key of its own:
# resolve_nsp's fallback IS the coronavirus table. Every other family falls
# through to that fallback too, and letting a picornavirus row read with
# CORONAVIRUS patterns is the cross-family naming bug
# diagnose_flavivirus_targets.py was written to find -- so a family absent
# from this set is refused BY NAME rather than by trusting the coronavirus
# patterns not to match it.
DESCRIPTION_TABLE = frozenset(PATTERNS_BY_FAMILY) | {"Coronaviridae"}

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
    "compound_structures": ("molregno", "standard_inchi_key"),
}

# TWO COMPOUND IDENTITIES, BOTH THE PROJECT'S OWN.
#
# The first version counted compounds by molecule_dictionary.chembl_id and
# found 4 compounds across arm B's 3,049 x 2,125 ends. A ChEMBL id is a
# REGISTRATION, not a molecule: a free base and its hydrochloride, or the
# same compound entered twice, carry different ids. This project never uses
# that id as identity -- ingest_activities keys every SmallMolecule node
# `INCHIKEY:<standard_inchi_key>` and drops a row with no key.
#
# The skeleton is the first 14 characters, the InChIKey's connectivity
# block. normalize_chem.same_as_rule already declares two compounds sharing
# it to be one node, with merge_rule
# "shared_inchikey_skeleton:salt_or_stereo_or_tautomer". So the arms are
# reported at both: the exact key is the graph's node, the skeleton is the
# graph's own declared merge, and the gap between them is how much of an
# arm is salt or stereo bookkeeping rather than chemistry.
IDENTITIES = (("exact InChIKey  (the graph's node id)", False),
              ("InChIKey skeleton  (its declared SAME_AS merge)", True))


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


ACTIVITY_FIELDS = ("tid", "compound", "inchikey", "drug_name", "max_phase",
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
        SELECT a.tid, md.chembl_id, cs.standard_inchi_key, md.pref_name,
               md.max_phase, act.standard_type, act.standard_value,
               act.standard_units, act.standard_relation, a.description,
               a.assay_type
        FROM activities act
        JOIN assays a ON act.assay_id = a.assay_id
        JOIN molecule_dictionary md ON act.molregno = md.molregno
        LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
        WHERE a.tid IN ({qmarks})
          AND act.standard_value IS NOT NULL
          AND act.standard_type IN ({types})
    """, list(tids) + list(QUANT_TYPES)).fetchall()
    return [dict(zip(ACTIVITY_FIELDS, r)) for r in rows]


def row_role(t: dict, description: str | None) -> str | None:
    """The enzyme one activity row measures, or None.

    A target named for its enzyme supplies the role directly. A polyprotein
    target supplies it only through the assay description, read with the
    family's own table -- the mechanism ingest_activities already uses.
    """
    if not t["protein_level"]:
        return None
    if t["role"]:
        return t["role"]
    if not t["polyprotein"] or t["family"] not in DESCRIPTION_TABLE:
        return None
    nsp, domain = resolve_nsp(description or "", t["family"])
    return POLYPROTEIN_ROLE.get((t["family"], nsp, domain))


def role_key(t: dict, role: str | None = None) -> str | None:
    """The role name the arms use, with HCV split out of the flaviviruses.

    Arm B is about HCV specifically: it is the one serine protease in the
    comparison with a large approved inhibitor class, and the one the schema
    names. The other six NS3 proteases are a separate comparison.
    """
    role = role if role is not None else t["role"]
    if role is None:
        return None
    if t["family"] == "Flaviviridae" and role == "ns3pro":
        return "hcv_ns3pro" if t["hcv"] else "flavi_ns3pro"
    return role


def compound_id(a: dict, *, skeleton: bool = False) -> str | None:
    """The compound this row is about, as the GRAPH identifies compounds.

    None when ChEMBL holds no structure for it -- ingest_activities drops
    those rows too, under the counter `no_inchikey`, so a survey that
    counted them would be counting compounds the graph will never build.
    """
    key = a.get("inchikey")
    if not key:
        return None
    return key[:14] if skeleton else f"INCHIKEY:{key}"


def compounds_by_role(acts: list[dict], targets: dict[int, dict],
                      *, protein_only: bool = True,
                      skeleton: bool = False) -> dict[str, set[str]]:
    """role key -> the compounds measured on it, one entry per compound."""
    out: dict[str, set[str]] = defaultdict(set)
    for a in acts:
        t = targets.get(a["tid"])
        if t is None:
            continue
        key = role_key(t, row_role(t, a["description"]))
        if key is None:
            continue
        cid = compound_id(a, skeleton=skeleton)
        if cid is None:
            continue
        if to_nm(a["standard_value"], a["standard_units"]) is None:
            continue
        if protein_only:
            where, _ = assay_readout(a["description"] or "", a["assay_type"])
            if where != "protein":
                continue
        out[key].add(cid)
    return out


def polyprotein_rows(acts: list[dict], targets: dict[int, dict]
                     ) -> tuple[Counter, Counter, dict[str, Counter]]:
    """Polyprotein rows in THREE buckets, not two.

    (arm role -> rows, other enzyme -> rows, family -> unknown text -> rows)

    THE EARLIER TWO-BUCKET SPLIT WAS WRONG, in the way this session keeps
    finding. "rows it did NOT resolve: 4,533" put together two different
    things, and 1,438 of the coronavirus share were a single description
    reading "...targets SARS-COV-2 nsp3 macrodomain...". resolve_nsp reads
    that perfectly -- NSP_PATTERNS' first entry is the macrodomain rule --
    and returns ("nsp3", "macrodomain"), which is simply not a protease and
    so not an arm role. Filing it beside the picornavirus rows, which the
    tables genuinely cannot read, made the coronavirus table look an order
    of magnitude more broken than it is.

    A row is UNKNOWN only when the family's table returned nothing, or when
    the family has no table at all. That bucket, and only that bucket, is
    evidence that a pattern is missing.
    """
    arm: Counter = Counter()
    other: Counter = Counter()
    unknown: dict[str, Counter] = defaultdict(Counter)
    for a in acts:
        t = targets.get(a["tid"])
        if t is None or not t["polyprotein"] or not t["protein_level"]:
            continue
        key = role_key(t, row_role(t, a["description"]))
        if key:
            arm[f"{t['family']}  {key}"] += 1
            continue
        nsp, domain = ((None, None) if t["family"] not in DESCRIPTION_TABLE
                       else resolve_nsp(a["description"] or "", t["family"]))
        if nsp:
            label = f"{nsp} {domain}" if domain else nsp
            other[f"{t['family']}  {label}"] += 1
        else:
            unknown[t["family"]][a["description"] or "(no description)"] += 1
    return arm, other, unknown


def usable_potency(acts: list[dict], targets: dict[int, dict],
                   *, skeleton: bool = False
                   ) -> tuple[dict[str, dict[str, float]], Counter, dict]:
    """(role -> compound -> median log10 nM, drops, assay-type mix).

    A COMPOUND WITH A PROTEIN-LEVEL ROW IS NOT A COMPOUND WITH A USABLE
    NUMBER, and the arms were sized on the first. A row still has to carry
    an exact relation and a convertible unit before its potency can enter a
    correlation, so the arm a protocol can actually run on is this one, not
    the membership count printed above it.
    """
    per: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list))
    drops: Counter = Counter()
    for a in acts:
        t = targets.get(a["tid"])
        if t is None:
            continue
        key = role_key(t, row_role(t, a["description"]))
        if key is None:
            drops["no_role"] += 1
            continue
        cid = compound_id(a, skeleton=skeleton)
        if cid is None:
            drops["no_inchikey"] += 1
            continue
        where, _ = assay_readout(a["description"] or "", a["assay_type"])
        if where != "protein":
            drops["not_a_protein_readout"] += 1
            continue
        if (a["standard_relation"] or "") != EXACT_RELATION:
            drops["censored_relation"] += 1
            continue
        nm = to_nm(a["standard_value"], a["standard_units"])
        if nm is None or nm <= 0:
            drops["unconvertible_or_nonpositive"] += 1
            continue
        per[(key, cid)][a["standard_type"]].append(math.log10(nm))

    out: dict[str, dict[str, float]] = defaultdict(dict)
    mix: Counter = Counter()
    for (key, cid), by_type in per.items():
        for stype in POTENCY_PREFERENCE:
            if by_type.get(stype):
                out[key][cid] = median(by_type[stype])
                mix[f"{key}  {stype}"] += 1
                break
    return out, drops, mix


def role_organisms(acts: list[dict], targets: dict[int, dict],
                   usable: dict[str, dict[str, float]],
                   *, skeleton: bool = False) -> dict[str, Counter]:
    """role -> organism -> how many of the USABLE compounds it supplied.

    Arm A pools four picornaviruses, and 59 of the 199 rows are one
    Coxsackievirus B3 paper. Whether to pool is a prereg decision and it
    cannot be taken without seeing the split.

    Restricted to `usable` on purpose. Counted over every role-resolved row
    instead, this table said 3clpro drew on 31 compounds while the arm
    beside it used 26 -- a split describing a set the protocol does not
    read.
    """
    seen: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for a in acts:
        t = targets.get(a["tid"])
        if t is None:
            continue
        key = role_key(t, row_role(t, a["description"]))
        cid = compound_id(a, skeleton=skeleton)
        if key is None or cid is None or cid not in usable.get(key, {}):
            continue
        seen[key][t["organism"]].add(cid)
    return {k: Counter({org: len(c) for org, c in v.items()})
            for k, v in seen.items()}


def broad_coronavirus_compounds(acts: list[dict], targets: dict[int, dict],
                                *, skeleton: bool = False) -> set[str]:
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
        cid = compound_id(a, skeleton=skeleton)
        if cid is not None and t is not None and t["family"] == "Coronaviridae":
            out.add(cid)
    return out


# THE DESIGN IS A 2x2, AND THE FIRST VERSION HAD TWO OF ITS FOUR CELLS.
#
# The gate makes its decision on two facts: do the chains share a fold, and
# is the catalytic residue the same. A and B hold the fold constant and vary
# the nucleophile, which tests the nucleophile rule and says nothing about
# whether the FOLD is doing any work. If compounds overlap just as much
# between two cysteine proteases with NO shared fold, then what A measures
# is compound promiscuity and screening-panel habit, not conservation.
#
# So both of the missing cells are here, as controls. They are not claims
# the premise makes; they are what makes A's number readable.
#
#   shared fold   same nucleophile   pair
#   yes           yes                A, A2   gate ALLOWS
#   yes           no                 B, B'   gate BLOCKS on catalytic type
#   no            yes                C1      gate blocks for want of a fold
#   no            no                 C2      nothing links them at all
#
# A2 applies the gate's own rule WITHIN Flaviviridae rather than across
# families: HCV NS3 and the other flavivirus NS3 proteases are serine both
# and share the fold, so gate() allows them. It is not one of the 36
# cross-family transfers the layer emits, and it is here because arm A is
# twelve compounds and a second allowed pair with more of them is worth
# having before anything is concluded from twelve.
ARMS = (
    ("A   ALLOWED  CoV 3CLpro + picornavirus 3C", "3clpro", "3cpro", "allowed",
     "shared fold, cysteine both — gate PERMITS"),
    ("A2  ALLOWED  HCV NS3 + other flavivirus NS3", "hcv_ns3pro",
     "flavi_ns3pro", "allowed",
     "shared fold, serine both — gate PERMITS (within family)"),
    ("B   REFUSED  CoV 3CLpro + HCV NS3", "3clpro", "hcv_ns3pro", "refused",
     "shared fold, cysteine vs serine — gate BLOCKS"),
    ("B'  REFUSED  CoV 3CLpro + other flavivirus NS3", "3clpro",
     "flavi_ns3pro", "refused",
     "shared fold, cysteine vs serine — gate BLOCKS"),
    ("C1  CONTROL  CoV 3CLpro + CoV PLpro", "3clpro", "plpro", "control",
     "NO shared fold, cysteine both — isolates the fold"),
    ("C2  CONTROL  CoV PLpro + HCV NS3", "plpro", "hcv_ns3pro", "control",
     "NO shared fold, cysteine vs serine — the floor"),
)

# POTENCY, AND WHY THE RULE IS A MEDIAN RATHER THAN THE BEST VALUE.
#
# A compound usually carries several rows on one enzyme. Taking the most
# potent is what this project already knows not to do: chembl.py's own note
# records chloroquine's "IC50 160 nM on Spike" winning a best-potency
# selection by 44x over its real 7,000 nM measurement, because the 160 came
# from a pseudovirus entry assay. A median over the protein-level rows
# cannot be captured by one outlier that way.
#
# ONE ASSAY TYPE AT A TIME. An IC50 and a Ki are not the same quantity and
# averaging them produces a number that is neither. The first type present
# in this order is used and the rest are ignored, with the mix reported so
# the prereg can argue with the order rather than inherit it silently.
POTENCY_PREFERENCE = ("IC50", "Ki", "Kd", "EC50")

# Censored rows -- ">", "<" -- are a bound, not a measurement, and a
# correlation cannot use a bound. chembl.VALID_RELATIONS keeps them with
# their relation and counts them; so does this.
EXACT_RELATION = "="


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--samples", type=int, default=8,
                    help="distinct descriptions to print per family that HAS "
                         "a table; a family with none prints all of them")
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
    print("\nprotease targets NAMED for their enzyme, by family. A "
          "polyprotein target\nnames no enzyme and is read from its assay "
          "descriptions instead, below:")
    for family in FAMILIES:
        fam = sorted((t for t in targets if t["family"] == family),
                     key=lambda t: (t["role"] or "~", t["pref_name"]))
        roled = [t for t in fam if t["role"]]
        poly_here = [t for t in fam if t["polyprotein"]]
        neither = [t for t in fam if not t["role"] and not t["polyprotein"]]
        print(f"  {family}  —  {len(roled)} named for a protease, "
              f"{len(poly_here)} polyprotein, {len(neither)} neither")
        print_capped(roled, 14,
                     lambda t: (f"{t['chembl_id']:<14} {role_key(t):<13} "
                                f"{'' if t['protein_level'] else 'EXCLUDED '}"
                                f"{t['target_type']:<16} {t['pref_name'][:44]}"),
                     indent="      ", unit="more named for a protease")
        # A FAMILY WITH NO ROLED TARGET PRINTED NOTHING AT ALL, which is how
        # "Picornaviridae -- 0 named for a protease" arrived with 74 targets
        # standing behind it and no way to see what they were called. If the
        # patterns claimed none of them, the names are the evidence for
        # whether that is right.
        if not roled and neither:
            print(f"      none claimed. What the {len(neither)} unclaimed "
                  f"targets are called:")
            print_capped(sorted(neither, key=lambda t: t["pref_name"]), 10,
                         lambda t: f"  {t['target_type']:<16} "
                                   f"{t['pref_name'][:56]}",
                         indent="      ", unit="more unclaimed target(s)")

    # ------------------------------------------------------- the chemistry
    acts = load_activities(con, [t["tid"] for t in targets])
    con.close()
    if not acts:
        print("\nnot one quantified row on any matched target.\n"
              "Exiting 2: nothing was measured.")
        return 2
    rows_per_tid = Counter(a["tid"] for a in acts)

    def _rows(ts):
        return sum(rows_per_tid[t["tid"]] for t in ts)

    # POLYPROTEIN TARGETS, READ RATHER THAN REJECTED. See the note beside
    # POLYPROTEIN above: rejecting these by name threw out 6,690 coronavirus
    # protease rows and produced a NOT MEASURABLE verdict that was an
    # artifact of the filter.
    poly = [t for t in targets if t["polyprotein"] and t["protein_level"]]
    arm, other, unknown = polyprotein_rows(acts, by_tid)
    print(f"\n  POLYPROTEIN targets: {len(poly)} target(s), "
          f"{_rows(poly):,} quantified row(s)")
    print("  The target name says polyprotein; the ASSAY DESCRIPTION says "
          "which enzyme,\n  and that is how ingest_activities already "
          "resolves these rows\n  (target_resolution=chain_from_assay_"
          "description). Read the same way here:")
    print_capped(sorted(poly, key=lambda t: -rows_per_tid[t["tid"]]), 8,
                 lambda t: f"{t['chembl_id']:<14} {t['pref_name'][:26]:<28} "
                           f"{t['organism'][:30]:<32} "
                           f"{rows_per_tid[t['tid']]:>7,} row(s)",
                 indent="      ", unit="more polyprotein target(s)")
    print(f"\n    rows the description resolved to an ARM ROLE: "
          f"{sum(arm.values()):,}")
    for k, v in arm.most_common():
        print(f"      {k:<34} {v:>7,}")
    print(f"\n    rows it resolved to SOME OTHER ENZYME — read correctly, "
          f"just not a\n    protease, so no arm: {sum(other.values()):,}")
    for k, v in other.most_common():
        print(f"      {k:<34} {v:>7,}")
    unk_total = sum(sum(c.values()) for c in unknown.values())
    print(f"\n    rows the table could NOT READ AT ALL: {unk_total:,}")
    print("    Only this bucket is evidence that a pattern is missing.")
    for family in FAMILIES:
        c = unknown.get(family)
        if not c:
            continue
        has_table = family in DESCRIPTION_TABLE
        table = ("the family's own table"
                 if has_table else "NO DESCRIPTION TABLE EXISTS for this family")
        print(f"      {family}  {sum(c.values()):,} row(s) over "
              f"{len(c):,} distinct description(s)  —  {table}")
        # VERBATIM AND WHOLE. Truncated to 96 characters, the commonest
        # coronavirus description read "Functional biochemical assay to
        # identify treatments for viral infectious disease that targets SA"
        # -- cut exactly where it names the enzyme. A gap in a pattern table
        # cannot be judged from a description clipped before the word that
        # matters.
        #
        # AND EVERY ONE OF THEM, for a family with no table. That family is
        # precisely the one whose table has to be written, the text is the
        # deliverable, and a capped sample is a table written from a guess
        # about the lines not shown. The flavivirus table says it was closed
        # by reading all 19 unresolved rows; this is the same situation.
        print_capped(c.most_common(), args.samples if has_table else len(c),
                     lambda kv: f"{kv[1]:>6,}  {kv[0]}",
                     indent="        ", unit="more distinct description(s)")

    # Not "and not polyprotein": a polyprotein target that is ORGANISM-typed
    # belongs here, and `poly` above already requires protein_level. The two
    # lists stay disjoint and together cover every matched target.
    nonprotein = [t for t in targets if not t["protein_level"]]
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
                           f"{rows_per_tid[t['tid']]:>7,} row(s)",
                 indent="      ", unit="more excluded target(s)")

    no_key = sum(1 for a in acts if not a.get("inchikey"))
    print(f"\n{len(acts):,} quantified rows on those targets "
          f"({', '.join(QUANT_TYPES)}); {no_key:,} carry no InChIKey and are "
          f"dropped,\nas ingest_activities drops them under `no_inchikey`.")

    # ------------------------------------------------------------ the arms
    # AT BOTH IDENTITIES. A ChEMBL id is a registration, not a molecule; the
    # graph's node id is INCHIKEY:<key> and its declared SAME_AS merge is the
    # 14-character skeleton. Counting by chembl_id, as the first version did,
    # reported 4 compounds across a 3,049 x 2,125 pair of ends.
    results: dict[str, set[str]] = {}
    broadest: set[str] = set()
    for ident, skel in IDENTITIES:
        protein = compounds_by_role(acts, by_tid, protein_only=True,
                                    skeleton=skel)
        anywhere = compounds_by_role(acts, by_tid, protein_only=False,
                                     skeleton=skel)
        print(f"\n=== {ident}")
        print("  compounds per role — protein-level readout is what the arms "
              "use; `any`\n  includes cellular rows filed against the same "
              "enzyme:")
        print(f"    {'role':<16} {'protein':>9} {'any readout':>13}")
        for key in sorted(set(protein) | set(anywhere)):
            print(f"    {key:<16} {len(protein.get(key, ())):>9,} "
                  f"{len(anywhere.get(key, ())):>13,}")

        pot, drops, mix = usable_potency(acts, by_tid, skeleton=skel)
        print("\n  THE ARMS — `both` is compounds measured at both ends; "
              "`usable` is\n  those with a potency a correlation can read, "
              "after the rules below:")
        print(f"    {'arm':<46} {'both':>5} {'usable':>7}   what the gate says")
        for label, left, right, kind, says in ARMS:
            both = protein.get(left, set()) & protein.get(right, set())
            pair = set(pot.get(left, {})) & set(pot.get(right, {}))
            # The widest identity wins: a compound the graph would merge is
            # one compound, and refusing on the narrower count would refuse
            # a measurable arm over salt bookkeeping. USABLE is what the
            # verdict reads, because membership without a readable number
            # sizes a protocol that cannot run.
            if kind != "control":
                results[label] = max(results.get(label, set()), pair, key=len)
            print(f"    {label:<46} {len(both):>5} {len(pair):>7}   {says}")
        if not skel:
            print(f"\n    rows dropped before a potency could be read: "
                  f"{sum(drops.values()):,}")
            for k, v in drops.most_common():
                print(f"      {k:<34} {v:>8,}")
            print("\n    assay type each compound's number came from "
                  f"(preference {' > '.join(POTENCY_PREFERENCE)}):")
            print_capped(mix.most_common(), 10,
                         lambda kv: f"{kv[0]:<30} {kv[1]:>7,}",
                         indent="      ", unit="more role/type pair(s)")
            orgs = role_organisms(acts, by_tid, pot, skeleton=skel)
            print("\n    which organisms supplied each role's compounds — "
                  "arm A pools four\n    picornaviruses, and whether to pool "
                  "is a decision the prereg has to take:")
            for key in sorted(orgs):
                print(f"      {key}")
                print_capped(orgs[key].most_common(), 5,
                             lambda kv: f"  {kv[0][:44]:<46} {kv[1]:>6,}",
                             indent="      ", unit="more organism(s)")

        broad = broad_coronavirus_compounds(acts, by_tid, skeleton=skel)
        hcv_broad = protein.get("hcv_ns3pro", set()) & broad
        broadest = max(broadest, hcv_broad, key=len)
        print(f"\n    the broad set: HCV NS3 enzyme inhibitors with ANY "
              f"coronavirus\n    measurement (enzyme or antiviral, any "
              f"target_type): {len(hcv_broad)}")
    print("\n  The broad set is the retrospective as the literature ran it in "
          "2020. It\n  answers a weaker question than arm B -- whether the "
          "compounds stopped the\n  virus, not whether they inhibited Mpro -- "
          "and is reported so a thin arm B\n  has a fallback rather than a "
          "dead end.")

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
        print(f"NO PAIRED ARM IS MEASURABLE. None reaches {FLOOR} compounds "
              f"measured at both\nends. {FLOOR} is "
              f"calibration.MIN_REACHED_POSITIVES, the floor this project "
              f"already\nuses for a class worth evaluating; below it a figure "
              f"describes its own\nsmall sample.")
        # THE SCRIPT CALLED THE BROAD SET THE FALLBACK AND THEN DID NOT LOOK
        # AT IT. The first version checked only ARMS against the floor and
        # printed "Phase 4 needs another design -- the broad set above is
        # the candidate" while that candidate stood at 17, above the floor,
        # one screen higher. A verdict that names an alternative owes it the
        # same test.
        if len(broadest) >= FLOOR:
            print(f"\nTHE BROAD SET DOES: {len(broadest)} compounds, at or "
                  f"above the same floor.\nIt is the weaker claim -- whether "
                  f"HCV NS3 inhibitors stopped a coronavirus,\nnot whether "
                  f"they inhibited its protease -- and it is the design phase "
                  f"4\nhas left. The prereg has to say that plainly, because "
                  f"a reader will take\n'the HCV retrospective' to mean the "
                  f"enzyme test unless told otherwise.")
        else:
            print(f"\nThe broad set does not reach it either: "
                  f"{len(broadest)} compounds.")
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
