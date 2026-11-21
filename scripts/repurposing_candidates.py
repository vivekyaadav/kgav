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

AND SO ARE DIAGNOSTIC AGENTS, AND SO IS RESEMBLANCE WEAKER THAN RESEMBLANCE
THAT ALREADY FAILED. A similarity score knows what a compound LOOKS like and
nothing about what it is for, so without the two filters in the "chemistry"
section below this list offers an ophthalmic dye as an antiviral, and offers
compounds less like the known actives than the median compound that was
screened against this very virus and came back negative. Both filters name
and count what they removed, and the floor is derived per cell rather than
set by a constant: the resemblance a failure achieves runs from 0.21 to 0.73
across these fourteen viruses, so one typed threshold is wrong in both
directions at once. Every row that survives carries its own position in
that distribution, because a bare Tanimoto cannot be read across cells.
"""
from __future__ import annotations

import argparse
import bisect
import json
import statistics
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
    # "WITH A MEASURED HIT RATE" WAS FALSE AND SAT HERE FOR WEEKS.
    # This line read "every row below is a chemical-analogy suggestion
    # with a measured hit rate", and no hit rate had ever been measured
    # -- checked across structure_scorer.json, its two controls,
    # temporal_structure_scorer.json and this script's own artifact. What
    # it quoted was AUC, which is the probability a random active
    # outranks a random inactive: a statement about the WHOLE ordering,
    # and the one translation that does not yield a hit rate. A screen
    # only ever orders the top of the list.
    #
    # Measured on 2026-11-20 by scripts/hit_rate_at_k.py, and the reason
    # this text now points at that report instead of claiming a number:
    # on the one population matched by assay rather than by source, the
    # precision of the top of this ranking falls to zero once each
    # candidate's nearest known actives are withheld. The uncontrolled
    # figure is large. Both registered controls take it to nothing. No
    # figure is typed here, because a hit rate copied into a banner is a
    # hit rate that goes stale silently -- which is how the false claim
    # survived.
    lines.append("  So every row below is a chemical-analogy suggestion. "
                 "It is NOT a")
    lines.append("  prediction of activity, NOT a mechanism and NOT a hit "
                 "rate: the")
    lines.append("  evidence is 'resembles a known active', and the known "
                 "active is")
    lines.append("  named so you can judge.")
    lines.append("")
    lines.append("  WHAT THAT RESEMBLANCE IS WORTH: `make hit-rate`. It "
                 "measures how")
    lines.append("  many of the top-ranked compounds actually came back "
                 "active, against")
    lines.append("  the base rate of the population they were ranked in, "
                 "at three")
    lines.append("  exclusion strengths. Read it before ordering anything. "
                 "On the one")
    lines.append("  assay-matched cell the precision of this ranking's "
                 "head does not")
    lines.append("  survive the near-neighbour control, and the AUC above "
                 "cannot show")
    lines.append("  that, because it averages over a ranking a screen "
                 "never sees.")
    return lines


def _is(n: int) -> str:
    """"1 are antineoplastics" announces that a template wrote the line,
    and this project fixed the same defect in the agent's node summaries
    two days ago. Verb agreement is cheap."""
    return "is" if n == 1 else "are"


def _does(n: int) -> str:
    return "does" if n == 1 else "do"


def _displaced_names(rows: list[dict], ann: dict, limit: int = 6) -> str:
    """The rows the floor removed, by name where ChEMBL has one.

    NAMED, NOT COUNTED. "the floor removed 7 rows" cannot be checked
    against the previous list; "TEMAZEPAM 0.517, FLUORESCEIN 0.522" can.
    """
    shown = [
        ((ann.get(r["candidate"]) or {}).get("name")
         or r["chembl_id"] or r["candidate"].split(":", 1)[-1][:14])
        + f" {r['tanimoto']:.3f}" for r in rows[:limit]]
    more = f", and {len(rows) - limit} more" if len(rows) > limit else ""
    return ", ".join(shown) + more


def virus_names() -> dict[str, str]:
    doc = load_register()
    return {f"NCBITaxon:{e['taxon']}": e["name"]
            for e in (doc.get("viruses") or [])
            if e.get("taxon") and e.get("name")}


# --------------------------------------------------------------- chemistry
#
# WHAT A MAX-TANIMOTO RANKING CANNOT KNOW. The score says "this looks like
# that" and nothing else, so the list it produces contains, by construction,
# every compound that looks like a known active for ANY reason. Reading the
# committed artifact on 2026-11-21 found three kinds of row a referee
# reaches for first:
#
#   A DYE MATCHED TO A DYE. FLUORESCEIN, Tanimoto 0.5217, printed under
#   "APPROVED DRUGS -- orderable now" for West Nile and Yellow fever,
#   resembling ERYTHROSINE. Fluorescein is an ophthalmic dye; erythrosine
#   is a food colourant. Nobody doses either for a flavivirus, and the row
#   says more about ECFP4 and xanthene rings than about either virus.
#
#   CYTOTOXIC ANTIMETABOLITES, UNMARKED. GEMCITABINE, CYTARABINE,
#   CLOFARABINE, NELARABINE, AZACITIDINE, FLOXURIDINE, CLADRIBINE. They are
#   nucleoside analogues, so they resemble nucleoside antivirals by
#   construction, and in a cell-based antiviral assay a compound that kills
#   the host cell lowers the viral readout without touching the virus. That
#   is a confound to declare, NOT a reason to drop the row: nucleoside
#   analogues are an antiviral family, and oncology-to-antiviral is a real
#   repurposing direction.
#
#   A TANIMOTO NOBODY CAN READ. TEMAZEPAM, a benzodiazepine hypnotic, at
#   0.5167 for Yellow fever, called tail noise on the strength of the
#   number alone -- and a Tanimoto has no absolute meaning. 0.52 sits near
#   the top of what MERS-CoV's screened-and-failed compounds reach and
#   under the median of HCoV-OC43's. Measured, Yellow fever's own failures
#   have a median of 0.3008, so its own standard KEEPS that row: it
#   resembles the actives better than 77% of the compounds screened
#   against Yellow fever that came back negative. The fix is therefore not
#   a threshold that removes it but a column that reads it -- see
#   `failed_resemblance` and `>fail`.
#
# THREE MECHANISMS, DELIBERATELY DIFFERENT: one drops the row and names
# it, one labels it, one puts a number in units a reader can use. Which
# rule got which was decided by measuring whether the field it rests on
# means what its name says -- see NON_THERAPEUTIC_ATC and the note on
# `withdrawn_flag` below.

#: ATC classes whose compounds exist to be MEASURED rather than to act:
#: contrast media, radiotracers, provocation agents, ophthalmic dyes. A
#: compound all of whose ATC classes fall in here is not a repurposing
#: candidate in any sense, because it is not a therapy at all. Those are
#: dropped from the pool and named in the output.
#:
#: EVERY CLASS, NOT ANY, and that is the whole rule. TOLBUTAMIDE is V04CA01,
#: a provocation test for insulin secretion, AND A10BB03, a sulfonylurea
#: people take daily. So are FOLIC ACID (V04CX02 + B03BB01), DEXTROSE,
#: SORBITOL and GONADORELIN. "Any class is diagnostic" dropped all five;
#: "every class is diagnostic" drops 25 compounds of the 19,884 annotated,
#: and they are contrast media, F-18 tracers and fluorescein.
NON_THERAPEUTIC_ATC = {
    "V04": "diagnostic agents",
    "V07": "all other non-therapeutic products",
    "V08": "contrast media",
    "V09": "diagnostic radiopharmaceuticals",
    "S01J": "ophthalmological diagnostic agents",
}

#: ATC L01, antineoplastic agents. FLAGGED AND KEPT. The confound is real
#: -- a cytotoxic nucleoside lowers a cell-based viral readout by killing
#: the cell -- and so is the class, so the row carries the label and the
#: reader decides.
CYTOTOXIC_ATC = "L01"

ROLE_DIAGNOSTIC = "diagnostic_agent"
ROLE_CYTOTOXIC = "cytotoxic"
ROLE_WITHDRAWN = "withdrawn"
ROLE_BLACK_BOX = "black_box_warning"
ROLE_NO_THERAPEUTIC_USE = "no_therapeutic_use"

#: How a role prints beside a row. ROLE_DIAGNOSTIC has no entry because a
#: diagnostic agent never reaches a row: it leaves the pool before anything
#: is scored, and is named where it would have appeared.
_MARK = {
    ROLE_CYTOTOXIC: f"cytotoxic (ATC {CYTOTOXIC_ATC} antineoplastic)",
    ROLE_WITHDRAWN: "withdrawn somewhere",
    ROLE_BLACK_BOX: "black-box warning",
    ROLE_NO_THERAPEUTIC_USE: "ChEMBL records no therapeutic use",
}

#: Fewer failures than this is not a distribution, so the cell gets no
#: resemblance floor and the output says so rather than quietly applying
#: a median of three numbers. Every cell in v0.1 has at least 27.
MIN_FLOOR_INACTIVES = 10


def annotations(chembl_ids: set[str]) -> dict[str, dict]:
    """ChEMBL id -> name, ATC classes, withdrawal and black-box status.

    THE GRAPH HOLDS NEITHER A NAME NOR A ROLE, AND THAT IS WHY THIS IS
    HERE. A SmallMolecule node carries `chembl_id`, `smiles`,
    `inchikey_skel` and the identity flags -- nothing a person can read,
    and nothing about what the compound is for. A list of accession
    numbers cannot be judged by anyone, which is the difference between a
    result and a read-out; and a list of NAMES with no roles is judged
    WRONGLY, because FLUORESCEIN reads like a drug. Both are looked up
    from the raw ChEMBL database rather than added to the release,
    because a display label and a clinical class are not facts this graph
    asserts.

    ABSENT DATABASE IS NOT FATAL AND MUST NOT BE SILENT: the rows print
    their ChEMBL ids and `main` states that no compound's role could be
    checked. A filter that passes everything when its input is missing,
    without saying so, is the shape of defect this project keeps finding.
    """
    if not CHEMBL_DB.exists() or not chembl_ids:
        return {}
    import sqlite3
    con = sqlite3.connect(f"file:{CHEMBL_DB}?mode=ro", uri=True)
    out: dict[str, dict] = {}
    ids = sorted(chembl_ids)
    for i in range(0, len(ids), 500):
        chunk = ids[i:i + 500]
        q = ("SELECT md.chembl_id, md.pref_name, md.withdrawn_flag, "
             "md.black_box_warning, md.therapeutic_flag, "
             "group_concat(mac.level5) "
             "FROM molecule_dictionary md "
             "LEFT JOIN molecule_atc_classification mac "
             "  ON mac.molregno = md.molregno "
             f"WHERE md.chembl_id IN ({','.join('?' * len(chunk))}) "
             "GROUP BY md.molregno")
        for cid, name, withdrawn, bbw, thera, atc in con.execute(q, chunk):
            out[cid] = {"name": name or None,
                        "atc": sorted((atc or "").split(",")) if atc else [],
                        "withdrawn": bool(withdrawn),
                        "black_box": bool(bbw),
                        "therapeutic": bool(thera)}
    con.close()
    return out


def preferred_names(chembl_ids: set[str]) -> dict[str, str]:
    """Just the names, for callers that want nothing else."""
    return {c: a["name"] for c, a in annotations(chembl_ids).items()
            if a["name"]}


def annotations_by_inchikey(keys: set[str]) -> dict[str, dict]:
    """The same annotation, for compounds ChEMBL knows only by structure.

    ChEMBL is keyed on the standard InChIKey, so where the release's
    normalisation left the key intact the record comes back; where it did
    not, the compound stays unannotated and is counted as such.
    """
    if not CHEMBL_DB.exists() or not keys:
        return {}
    import sqlite3
    con = sqlite3.connect(f"file:{CHEMBL_DB}?mode=ro", uri=True)
    out: dict[str, dict] = {}
    ks = sorted(keys)
    for i in range(0, len(ks), 500):
        chunk = ks[i:i + 500]
        q = ("SELECT cs.standard_inchi_key, md.chembl_id, md.pref_name, "
             "md.withdrawn_flag, md.black_box_warning, md.therapeutic_flag, "
             "group_concat(mac.level5) "
             "FROM compound_structures cs "
             "JOIN molecule_dictionary md ON md.molregno = cs.molregno "
             "LEFT JOIN molecule_atc_classification mac "
             "  ON mac.molregno = md.molregno "
             f"WHERE cs.standard_inchi_key IN ({','.join('?' * len(chunk))}) "
             "GROUP BY md.molregno")
        for key, cid, name, withdrawn, bbw, thera, atc in con.execute(
                q, chunk):
            out[key] = {"name": name or None, "chembl_id": cid,
                        "atc": sorted((atc or "").split(",")) if atc else [],
                        "withdrawn": bool(withdrawn),
                        "black_box": bool(bbw),
                        "therapeutic": bool(thera)}
    con.close()
    return out


def pool_annotations(props: dict) -> dict[str, dict]:
    """node id -> what ChEMBL knows about that compound.

    KEYED ON THE NODE AND NOT ON THE CHEMBL ID, and that is the whole
    reason this function exists rather than a call to `annotations`.
    5,359 of the 25,243 small molecules in v0.1 carry NO chembl_id --
    they came from sources that identify a compound by structure -- so a
    role filter reading the chembl_id column alone would have inspected
    79% of the pool and reported itself as applied.

    MEASURED, BECAUSE THE GAP IS NOT UNIFORM. Of those 5,359, ChEMBL
    recognises 3,121 by standard InChIKey, and among THEM are 17 further
    diagnostic agents -- PHENOLSULFONPHTHALEIN (phenol red), INDOCYANINE
    GREEN, SULFOBROMOPHTHALEIN and most of the iodinated contrast series
    -- and 20 further antineoplastics, including MECHLORETHAMINE and
    THIOTEPA. The dyes this filter exists to remove were concentrated in
    the part of the pool it could not see, which is what a filter keyed
    on the convenient column usually looks like.

    Compounds neither route resolves stay unannotated, and `main` prints
    how many those are: unchecked is not the same as clean.
    """
    by_id = annotations({p.get("chembl_id") for p in props.values()
                         if p.get("chembl_id")})
    out: dict[str, dict] = {}
    unresolved = []
    for node, prop in props.items():
        hit = by_id.get(prop.get("chembl_id") or "")
        if hit:
            out[node] = dict(hit, chembl_id=prop.get("chembl_id"),
                             resolved_by="chembl_id")
        elif node.startswith("INCHIKEY:"):
            unresolved.append(node)
    by_key = annotations_by_inchikey(
        {n.split(":", 1)[1] for n in unresolved})
    for node in unresolved:
        hit = by_key.get(node.split(":", 1)[1])
        if hit:
            out[node] = dict(hit, resolved_by="inchikey")
    return out


def compound_flags(ann: dict | None,
                   approved: bool = False) -> tuple[str, ...]:
    """What ChEMBL says about this compound that changes how a row reads.

    ROLE_DIAGNOSTIC is the only one this script filters on, because it is
    the only one whose field means what its name says.

    `withdrawn_flag` DOES NOT, and the check is cheap enough that it was
    worth running before writing a filter on it: 195 of the 1,676 approved
    compounds in this pool carry it, and the set contains CODEINE and
    THALIDOMIDE -- withdrawn somewhere, some time, for some formulation,
    both still prescribed today. So withdrawal is REPORTED and never
    filtered on. "A text marker is not an enforcement mechanism" is a
    defect this project has already named twice; a filter resting on a
    field that does not mean what its name says is the same defect wearing
    the other hat.
    """
    if not ann:
        return ()
    atc = ann.get("atc") or []
    out = []
    if atc and all(any(c.startswith(p) for p in NON_THERAPEUTIC_ATC)
                   for c in atc):
        out.append(ROLE_DIAGNOSTIC)
    if any(c.startswith(CYTOTOXIC_ATC) for c in atc):
        out.append(ROLE_CYTOTOXIC)
    if ann.get("withdrawn"):
        out.append(ROLE_WITHDRAWN)
    if ann.get("black_box"):
        out.append(ROLE_BLACK_BOX)
    # ONLY ON AN APPROVED COMPOUND, AND THAT RESTRICTION IS THE WHOLE
    # VALUE OF THE FIELD. `therapeutic_flag` is 0 for 18,251 of the
    # 19,884 annotated compounds, including BALAPIRAVIR and OBELDESIVIR
    # -- antivirals in development -- so as a general rule it means
    # little more than "not a marketed dosed product". On a compound the
    # release flags APPROVED it says something precise and worth saying:
    # the product is approved and ChEMBL records no therapeutic use for
    # it, which is what PHENOLPHTHALEIN, ERYTHROSINE and BETA CAROTENE
    # have in common and what makes "APPROVED DRUGS -- orderable now" the
    # wrong banner to meet them under.
    if approved and ann.get("therapeutic") is False:
        out.append(ROLE_NO_THERAPEUTIC_USE)
    return tuple(out)


def _max_tc(fp, ref_fps: list, ref_skel: list[str],
            skel: str) -> tuple[float, int] | None:
    """Best Tanimoto from one fingerprint to the references, and which one.

    SHARED BY THE SUGGESTIONS AND BY THE FLOOR ON PURPOSE. The floor is
    the resemblance a known INACTIVE achieves, and it is only comparable
    to a candidate's resemblance if both are computed the same way --
    same fingerprint, same skeleton exclusion. Two implementations of
    "max Tanimoto to the actives" would drift, and the drift would move
    the threshold rather than announce itself.
    """
    from rdkit import DataStructs

    if not ref_fps:
        return None
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    best, best_i = -1.0, None
    for i, s in enumerate(sims):
        if skel and ref_skel[i] == skel:
            continue
        if s > best:
            best, best_i = s, i
    return None if best_i is None else (best, best_i)


def failed_resemblance(actives: list[str], inactives: list[str], fps: dict,
                       props: dict) -> list[float]:
    """Sorted max-Tanimoto from each known INACTIVE to the known actives.

    THE REFERENCE DISTRIBUTION THIS READ-OUT HAS AND NEVER USED. Every
    compound in it was screened against this virus and came back
    negative, so it is the resemblance that has already been tried and
    failed -- the only yardstick on the page that is in the same units as
    a candidate's score. `resemblance_floor` reads a threshold off it and
    `suggest` reads each row's position in it; both have to be the same
    distribution, computed with the same skeleton exclusion, or the
    threshold and the percentile would describe different things.
    """
    ref_ids = [a for a in actives if a in fps]
    ref_fps = [fps[a] for a in ref_ids]
    ref_skel = [props.get(a, {}).get("inchikey_skel", "") for a in ref_ids]
    vals = []
    for c in inactives:
        if c not in fps:
            continue
        hit = _max_tc(fps[c], ref_fps, ref_skel,
                      props.get(c, {}).get("inchikey_skel", ""))
        if hit:
            vals.append(hit[0])
    vals.sort()
    return vals


def resemblance_floor(values: list[float], n_inactives: int) -> dict:
    """How much this virus's actives are resembled BY COMPOUNDS THAT FAILED.

    A FLAT 0.55 WAS ASKED FOR AND THE DATA REFUSES A CONSTANT. Measured on
    v0.1: the median max-Tanimoto from a screened-and-inactive compound to
    the same virus's known actives runs from 0.212 on MERS-CoV to 0.726 on
    HCoV-OC43, a factor of three. One typed threshold is therefore wrong
    in both directions at once -- far above MERS-CoV's entire inactive
    distribution, where it discards rows that do carry information, and
    far below HCoV-OC43's median, where it keeps rows that carry none.

    So the floor is the cell's own: the MEDIAN max-Tanimoto its known
    INACTIVES reach against its known actives, under the same skeleton
    exclusion the suggestions use. A candidate below it resembles the
    actives less than half the compounds that were actually screened
    against this virus and came back negative.

    WHAT THAT IS AND IS NOT. It is not a prediction that the compound is
    inactive, and it is not a second protocol -- no score changes and the
    cell's registered verdict still decides whether any list is printed.
    It is the statement that the row's only evidence, resemblance, is
    weaker here than resemblance that has already been tried and failed.

    The quartile and decile above the median ride along so the choice of
    statistic is visible rather than buried: a reader who wants a sterner
    floor can see where it would sit and how many rows it would cost.

    AND THE FLOOR IS NOT THE WHOLE ANSWER, which the measurement showed
    the moment it ran. TEMAZEPAM -- the row the review called tail noise
    -- sits at 0.5167 for Yellow fever, whose failures have a median of
    0.3008, so the cell's own standard KEEPS it. Rather than bend the
    threshold until the disliked row falls out, every row now carries
    where it sits in this distribution: `beats_failures`, printed as the
    `>fail` column. The floor is simply the 50% line on that column.
    """
    out = {"n_inactives_scored": len(values),
           "n_inactives": n_inactives,
           "floor": None, "statistic": "median", "p75": None, "p90": None,
           "min_inactives_required": MIN_FLOOR_INACTIVES}
    if len(values) < MIN_FLOOR_INACTIVES:
        return out
    out["floor"] = round(statistics.median(values), 4)
    out["p75"] = round(statistics.quantiles(values, n=4)[2], 4)
    out["p90"] = round(statistics.quantiles(values, n=10)[8], 4)
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


#: Clinical maturity, most actionable first. A lab asking "what can we
#: test next week" and a chemist asking "what should we make" are
#: different questions about the SAME ranking, and one list answers both
#: only if it says which rows are which.
TIERS = ("approved", "clinical", "research")

#: Printed above each tier's block, so a reader knows what changed
#: between them is clinical maturity and NOT the score.
_TIER_BANNER = {
    "approved": "APPROVED DRUGS -- orderable now",
    "clinical": "IN CLINICAL DEVELOPMENT (ChEMBL max_phase 1-3)",
    "research": "RESEARCH COMPOUNDS -- mostly analogs of the known actives",
}


def _phase_tier(prop: dict) -> str:
    """Which tier a compound belongs to.

    `is_approved` and `max_phase == 4` are not quite the same set (1,676
    against 1,680 on v0.1), so the approved flag wins where they differ:
    it is the field the rest of this script and its tests already read.
    `max_phase` is absent on 20,626 of 25,243 compounds, and absent means
    research -- not unknown-and-therefore-promoted.
    """
    if prop.get("is_approved"):
        return "approved"
    try:
        phase = float(prop.get("max_phase"))
    except (TypeError, ValueError):
        return "research"
    return "clinical" if phase >= 1 else "research"


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
            props: dict, top: int, degree: dict | None = None,
            flat: bool = False, floor: float | None = None,
            flags: dict[str, tuple[str, ...]] | None = None,
            failures: list[float] | None = None,
            ) -> tuple[list[dict], dict]:
    """Top `top` candidates by nearest known active, and what the
    resemblance floor cost.

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

    `floor` drops candidates whose best resemblance is under the value
    `resemblance_floor` measured for this cell. The second return value
    says what that cost: the rows the floor DISPLACED -- the ones the same
    selection would have printed without it -- plus how much of the pool
    sits below the floor. Returned rather than printed from inside, so the
    printout and the written file cannot disagree about it.

    `flags` is the compound-role map from `compound_flags`: it labels rows
    and never filters them, because the one role worth filtering on is
    removed from the pool before anything is scored.
    """
    ref_ids = [a for a in actives if a in fps]
    ref_fps = [fps[a] for a in ref_ids]
    ref_skel = [props.get(a, {}).get("inchikey_skel", "") for a in ref_ids]
    rows: list[dict] = []
    for c in candidates:
        fp = fps.get(c)
        if fp is None:
            continue
        hit = _max_tc(fp, ref_fps, ref_skel,
                      props.get(c, {}).get("inchikey_skel", ""))
        if hit is None:
            continue
        best, best_i = hit
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
        rows[-1]["max_phase"] = props.get(c, {}).get("max_phase")
        rows[-1]["phase_tier"] = _phase_tier(props.get(c, {}))
        # WHAT THE COMPOUND IS FOR, beside what it looks like. A
        # cytotoxic antimetabolite and an antiviral nucleoside are the
        # same row to a fingerprint; this is the column that separates
        # them, and the reference's flags are carried too because a
        # suggestion is only as readable as the thing it resembles.
        rows[-1]["flags"] = list((flags or {}).get(c, ()))
        rows[-1]["nearest_active_flags"] = list(
            (flags or {}).get(ref_ids[best_i], ()))
        rows[-1]["below_resemblance_floor"] = bool(
            floor is not None and best < floor)
        # WHERE THIS ROW SITS AMONG THE COMPOUNDS THAT FAILED. A bare
        # Tanimoto is not readable across cells -- 0.52 is near the top
        # of MERS-CoV's failure distribution and under the median of
        # HCoV-OC43's -- and this is the same number in units a reader
        # can use. The floor is the 50% line on it.
        rows[-1]["beats_failures"] = (
            None if not failures else
            round(bisect.bisect_left(failures, best) / len(failures), 4))
    rows.sort(key=lambda r: -r["tanimoto"])

    # TOP `top` WITHIN EACH TIER, NOT top overall. THE POOL IS NOT
    # RESTRICTED -- every small molecule is still scored and still
    # eligible -- but a single max-Tanimoto ranking never surfaces a drug
    # anyone can order. Measured on SARS-CoV-2 on 2026-11-19: of the top
    # 60 by Tanimoto, 35 are >= 0.90, the lowest is 0.874, and ZERO are
    # approved or even named. ~25,000 research analogues outrank every
    # approved compound, so a reader of the top of this list sees the
    # same scaffold sixty times and concludes the method found nothing.
    #
    # Tiering is not a second protocol and changes no score: the ranking
    # inside each tier is the identical statistic, and the cell's verdict
    # still gates whether any list is printed at all.
    out: list[dict] = []
    # FLAT IS ONE RANKING OF `top`, NOT `top` PER TIER.
    #
    # The tiers exist so a browsing reader sees approved drugs without
    # the pool being restricted to them. A SCREEN is a different use:
    # docs/prereg_denv2_prospective_screen.md section 3 fixes arm A as
    # "the top 20 by max-Tanimoto", one ranking -- and the command the
    # registration records, `--top 20`, returns SIXTY under tiering,
    # twenty from each of three tiers. Measured on DENV-2 before this
    # flag existed.
    def pick(rs: list[dict]) -> list[dict]:
        if flat:
            return sorted(rs, key=lambda r: -r["tanimoto"])[:top]
        picked: list[dict] = []
        for tier in TIERS:
            picked.extend([r for r in rs if r["phase_tier"] == tier][:top])
        picked.sort(key=lambda r: (TIERS.index(r["phase_tier"]),
                                   -r["tanimoto"]))
        return picked

    # THE FLOOR IS APPLIED BEFORE EACH TIER'S CAP, and that ordering is
    # the point: filtering after the cap would let a row with no evidence
    # behind it eat one of the fifteen slots a reader actually sees, and
    # the list would come back short with nothing to say why. Filtering
    # first refills the slot from the next candidate that clears.
    out = pick([r for r in rows if not r["below_resemblance_floor"]])
    # AND WHAT IT COST, NAMED. The number of candidates under the floor is
    # most of the pool and says little; the rows that WOULD have been
    # printed and are not is the figure a reader of the old list wants,
    # and it is the one that can be checked against the old list.
    displaced = [r for r in pick(rows) if r["below_resemblance_floor"]]
    return out, {"displaced": displaced,
                 "n_candidates_below_floor": sum(
                     1 for r in rows if r["below_resemblance_floor"]),
                 "n_candidates_scored": len(rows)}


#: Predicates that mean somebody has screened this compound against this
#: virus. Both directions: an inactive result is a screening result.
SCREENED_PREDICATES = ("HAS_ANTIVIRAL_ACTIVITY_AGAINST",
                       "MEASURED_INACTIVE_AGAINST")


def screened_against(release: Path, virus: str) -> set[str]:
    """Every compound with ANY activity edge to `virus`, either way.

    STRICTER THAN THE LABEL SETS, and that is the point. labels.positives
    and labels.negatives are what the potency thresholds could classify;
    a row whose value the thresholds could not read is in neither, and
    its compound has still been screened. A prospective screen that drew
    such a compound would be re-running somebody's experiment and
    calling the result a prediction.
    """
    out: set[str] = set()
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] not in SCREENED_PREDICATES:
            continue
        if e["object"] == virus:
            out.add(e["subject"])
        elif e["subject"] == virus:
            out.add(e["object"])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--top", type=int, default=DEFAULT_TOP)
    ap.add_argument("--approved-only", action="store_true",
                    help="restrict suggestions to compounds the release "
                         "flags is_approved -- the repurposing question in "
                         "its narrowest form")
    ap.add_argument("--virus", default=None)
    ap.add_argument("--exclude-labelled", action="store_true",
                    help="drop every candidate carrying ANY activity edge "
                         "to this virus, in either direction -- not only "
                         "those the potency thresholds classified. "
                         "docs/prereg_denv2_prospective_screen.md section 4 "
                         "requires 'never screened against DENV-2 in either "
                         "direction', which is stricter than the label sets: "
                         "a measurement the thresholds could not classify "
                         "still means somebody has screened it.")
    ap.add_argument("--flat", action="store_true",
                    help="one ranking of --top compounds by max-Tanimoto, "
                         "instead of --top from each of the three phase "
                         "tiers. What a two-arm screen needs: the "
                         "registration fixes arm A as 'the top 20', and "
                         "tiered that command returns sixty.")
    ap.add_argument("--random-arm", type=int, default=0, metavar="N",
                    help="also draw N compounds uniformly at random from "
                         "the same pool, as the comparator arm of a two-arm "
                         "screen. Requires --seed.")
    ap.add_argument("--seed", type=int, default=None,
                    help="the seed for --random-arm, recorded in the "
                         "artifact. Required with --random-arm so the draw "
                         "is reproducible and cannot be re-rolled.")
    # THE OFF-SWITCHES EXIST SO THE FILTERS CAN BE MEASURED. A filter
    # whose effect cannot be computed is a claim, and both of these change
    # which compounds a reader is asked to order. The tests run the script
    # both ways and assert on the difference.
    ap.add_argument("--keep-diagnostic-agents", action="store_true",
                    help="do NOT drop compounds all of whose ATC classes "
                         "are diagnostic, contrast or radiotracer classes. "
                         "Only useful for measuring what that filter "
                         "removes.")
    ap.add_argument("--no-resemblance-floor", action="store_true",
                    help="do NOT drop candidates whose best Tanimoto is "
                         "below the median a known INACTIVE of the same "
                         "virus reaches. Only useful for measuring what "
                         "that floor removes.")
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
    # skip_predicates IS LOAD-BEARING HERE, and its absence was a defect
    # until 2026-11-19. This script reports `graph_degree` per suggestion
    # from g.degree_unlabelled, and on a graph loaded without this
    # argument that counter is identical to `degree` -- so the column
    # said how connected a compound is while counting the held-out
    # labels the shortlist exists to predict. Measured on the committed
    # artifact: 163 of 195 reported degrees were inflated, and 78 of 195
    # compounds the graph does not reach at all were printed as reached.
    # The ranking is max-Tanimoto and never moved; what moved is what a
    # reader of the shortlist believes about each row.
    g = Graph.load(args.release,
                   skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
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

    # WHAT EACH COMPOUND IS FOR, ONCE, FOR THE WHOLE POOL. One ChEMBL
    # pass over the 19,884 annotated ids; the roles are needed BEFORE
    # scoring because the diagnostic filter has to remove a compound from
    # the pool rather than from the printed window -- see `suggest`.
    ann = pool_annotations(props)
    flags = {i: compound_flags(ann.get(i),
                               approved=bool(props[i].get("is_approved")))
             for i in props}
    # IDENTIFIED AND APPLIED ARE TWO DIFFERENT SETS. --keep-diagnostic-agents
    # must change what is DROPPED without changing what was FOUND, or the
    # census would report the filter's effect as the filter's input and the
    # off switch could not be used to measure anything.
    identified = {i for i, f in flags.items() if ROLE_DIAGNOSTIC in f}
    diagnostic = set() if args.keep_diagnostic_agents else identified
    if not ann:
        # NOT SILENT. A filter that passes everything because its input is
        # missing, and says nothing, is worse than no filter: the reader
        # believes the list was checked.
        print("NO CHEMISTRY ANNOTATION: the raw ChEMBL database is not on "
              "this checkout, so\n  NO compound's clinical role could be "
              "checked. Every row below is unfiltered\n  for role, and "
              f"may be a dye, a contrast agent or a radiotracer.\n"
              f"  Expected at {CHEMBL_DB}.\n")
    else:
        n_cyt = sum(1 for f in flags.values() if ROLE_CYTOTOXIC in f)
        n_wdn = sum(1 for f in flags.values() if ROLE_WITHDRAWN in f)
        n_key = sum(1 for a in ann.values()
                    if a.get("resolved_by") == "inchikey")
        print(f"clinical roles from ChEMBL: {len(ann):,} of "
              f"{len(props):,} molecules annotated "
              f"({n_key:,} of them by InChIKey, carrying no chembl_id in "
              f"the release),\n  so {len(props) - len(ann):,} could NOT "
              f"be role-checked at all -- unchecked, not clean.")
        print(f"  {len(identified):,} are diagnostic, contrast or "
              f"radiotracer agents (every ATC class in "
              f"{'/'.join(sorted(NON_THERAPEUTIC_ATC))})")
        print("  -- not therapies, so DROPPED from every list below and "
              "named where they would have appeared."
              if not args.keep_diagnostic_agents else
              "  -- KEPT: --keep-diagnostic-agents was passed.")
        print(f"  {n_cyt:,} are antineoplastics (ATC {CYTOTOXIC_ATC}) and "
              f"{n_wdn:,} are withdrawn somewhere: both are LABELLED, "
              f"never dropped.\n")
    payload = {"measured_by": "scripts/repurposing_candidates.py",
               "scorer": "max Tanimoto to known actives (ECFP4 2048/r2)",
               "registered_by":
                   "docs/prereg_structure_scorer_strength_matched.md",
               "approved_only": bool(args.approved_only),
               "top": args.top, "viruses": {}}
    # THE CENSUS IN THE FILE, so a document quoting these numbers can be
    # re-derived from the artifact instead of from a terminal somebody
    # scrolled past. docs/scope_and_limits.md quotes five of them and
    # tests/test_chemistry_filter.py checks each against this block.
    payload["chemistry_filter"] = {
        "atc_non_therapeutic": dict(NON_THERAPEUTIC_ATC),
        "atc_cytotoxic": CYTOTOXIC_ATC,
        "n_molecules": len(props),
        "n_annotated": len(ann),
        "n_annotated_by_inchikey": sum(
            1 for a in ann.values() if a.get("resolved_by") == "inchikey"),
        "n_unchecked": len(props) - len(ann),
        "n_diagnostic_identified": len(identified),
        "n_diagnostic_identified_by_inchikey": sum(
            1 for i in identified
            if (ann.get(i) or {}).get("resolved_by") == "inchikey"),
        "n_cytotoxic": sum(1 for f in flags.values()
                           if ROLE_CYTOTOXIC in f),
        "n_withdrawn": sum(1 for f in flags.values()
                           if ROLE_WITHDRAWN in f),
        "diagnostic_agents": sorted(
            (ann.get(i) or {}).get("name") or i for i in identified),
        "diagnostic_filter_applied": not args.keep_diagnostic_agents,
        "resemblance_floor_applied": not args.no_resemblance_floor,
    }

    # A DRAW THAT CANNOT BE REPRODUCED IS NOT A DRAW. The registration
    # records the seed in the artifact and in the commit message, so a
    # random arm without one is refused rather than given an arbitrary
    # default -- which would let the draw be re-rolled until it looked
    # right, with nothing on the record to show it had been.
    if args.random_arm and args.seed is None:
        print("--random-arm requires --seed: the draw has to be "
              "reproducible from the artifact", file=sys.stderr)
        return 2
    if args.seed is not None and not args.random_arm:
        print("--seed without --random-arm draws nothing", file=sys.stderr)
        return 2

    screened: set[str] = set()
    if args.exclude_labelled:
        if not args.virus:
            print("--exclude-labelled needs --virus: 'screened against "
                  "what' has no answer otherwise", file=sys.stderr)
            return 2

    for virus in sorted(set(labels.positives) | set(labels.negatives)):
        if args.virus and args.virus not in (virus, names.get(virus)):
            continue
        if args.exclude_labelled:
            screened = screened_against(args.release, virus)
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
        candidates, salts, screened_out = [], 0, 0
        dropped_diagnostic: list[str] = []
        for c in pool:
            if c in known:
                continue
            if c in screened:
                # Has an activity edge to this virus that the thresholds
                # did not classify. Known is a subset of this.
                screened_out += 1
                continue
            # NOT A THERAPY AT ALL. A contrast agent or an ophthalmic dye
            # cannot be repurposed as an antiviral, whatever it resembles,
            # so it leaves the pool here rather than being printed with a
            # label -- and it is named below, because a drop nobody can
            # see is indistinguishable from a compound that never scored.
            if c in diagnostic:
                dropped_diagnostic.append(
                    (ann.get(c) or {}).get("name")
                    or props[c].get("chembl_id") or c)
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
        # THE FLOOR IS THIS CELL'S OWN, measured from its own failures.
        # Computed even when no list will be printed, because "how much do
        # this virus's inactives resemble its actives" is a property of
        # the cell and belongs in the artifact either way.
        failures = failed_resemblance(actives, sorted(inactives), fps,
                                      props)
        fl = resemblance_floor(failures, len(inactives))
        floor = None if args.no_resemblance_floor else fl["floor"]
        if t_ in ("below", "unknown"):
            rows, cost = [], {"displaced": [], "n_candidates_below_floor": 0,
                              "n_candidates_scored": 0}
        else:
            rows, cost = suggest(
                candidates, actives, fps, props, args.top,
                flat=args.flat, floor=floor, flags=flags,
                failures=failures, degree=g.degree_unlabelled)
        # ----------------------------------------------- the comparator arm
        #
        # DRAWN FROM THE POOL ARM A WAS RANKED FROM, AND DISJOINT FROM IT.
        #
        # docs/prereg_denv2_prospective_screen.md section 3 says "k drawn
        # at random from the same pool". Read strictly that admits a
        # compound into BOTH arms, which makes the 2x2 incoherent -- the
        # same well would be a hit in both columns. Disjoint is the only
        # coherent reading and is recorded here as a choice rather than
        # assumed, because it is a detail the registration did not fix.
        #
        # The seed is recorded in the artifact. The draw is over a SORTED
        # pool so it depends on the seed and the release, not on dict
        # ordering.
        arm_b: list[str] = []
        if args.random_arm:
            import random as _random

            # "candidate" is the key suggest() emits the compound id
            # under. Guessing "inchikey" or "id" here would silently make
            # the arms non-disjoint: the set would be empty, every
            # filter would pass, and arm B could redraw arm A.
            chosen = {r["candidate"] for r in rows}
            eligible = sorted(c for c in candidates if c not in chosen)
            rng = _random.Random(args.seed)
            arm_b = (rng.sample(eligible, args.random_arm)
                     if len(eligible) >= args.random_arm else [])
            if not arm_b:
                print(f"  REFUSED the random arm: {len(eligible):,} "
                      f"eligible, {args.random_arm} needed")

        payload["viruses"][virus] = {
            "name": names.get(virus),
            "n_known_actives": len(actives),
            "n_known_inactives": len(inactives),
            "n_candidates": len(candidates),
            "n_salts_of_known_compounds_dropped": salts,
            "n_screened_against_this_virus_dropped": screened_out,
            "exclude_labelled": bool(args.exclude_labelled),
            # WHAT THE CHEMISTRY FILTERS REMOVED FROM THIS CELL, in the
            # file and not only on the terminal. The diagnostic agents are
            # named rather than counted: 25 compounds is a list a reader
            # can check, and a count is a claim they cannot.
            "chemistry_filter": {
                "diagnostic_agents_dropped": sorted(dropped_diagnostic),
                "n_diagnostic_agents_dropped": len(dropped_diagnostic),
                "diagnostic_filter_applied":
                    not args.keep_diagnostic_agents,
                "resemblance_floor": dict(
                    fl, applied=floor is not None,
                    n_candidates_below_floor=cost[
                        "n_candidates_below_floor"],
                    n_candidates_scored=cost["n_candidates_scored"],
                    n_rows_displaced=len(cost["displaced"]),
                    rows_displaced=[
                        {"chembl_id": r["chembl_id"],
                         "candidate": r["candidate"],
                         "tanimoto": r["tanimoto"],
                         "phase_tier": r["phase_tier"]}
                        for r in cost["displaced"]]),
                "annotated": bool(ann),
            },
            "random_arm": ({
                "seed": args.seed,
                "n_requested": args.random_arm,
                "n_drawn": len(arm_b),
                "disjoint_from_arm_a": True,
                "drawn_from": len(candidates),
                "compounds": arm_b,
            } if args.random_arm else None),
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
        if dropped_diagnostic:
            head = ", ".join(sorted(dropped_diagnostic)[:6])
            more = (f", and {len(dropped_diagnostic) - 6} more"
                    if len(dropped_diagnostic) > 6 else "")
            print(f"  {len(dropped_diagnostic):,} dropped as diagnostic, "
                  f"contrast or radiotracer agents, which are not "
                  f"therapies:\n      {head}{more}")
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
        if fl["floor"] is None:
            print(f"  NO RESEMBLANCE FLOOR: {fl['n_inactives']:,} measured "
                  f"inactives is under the {MIN_FLOOR_INACTIVES} this "
                  f"needs,\n           so this cell has no "
                  f"failed-resemblance distribution to compare a row to.")
        elif floor is None:
            print(f"  RESEMBLANCE FLOOR NOT APPLIED "
                  f"(--no-resemblance-floor). It would be "
                  f"{fl['floor']:.4f}.")
        else:
            print(f"  RESEMBLANCE FLOOR {fl['floor']:.4f}: the MEDIAN "
                  f"max-Tanimoto this virus's "
                  f"{fl['n_inactives_scored']:,}\n           "
                  f"screened-and-INACTIVE compounds reach against its "
                  f"known actives. A candidate\n           below it "
                  f"resembles the actives less than half the compounds "
                  f"that were\n           screened against this very "
                  f"virus and came back negative.\n"
                  + (f"\n           "
                     f"{cost['n_candidates_below_floor']:,} of "
                     f"{cost['n_candidates_scored']:,} candidates sit "
                     f"below it. It cost the list "
                     f"{len(cost['displaced'])} row(s): "
                     f"{_displaced_names(cost['displaced'], ann) or 'none'}"
                     if cost["n_candidates_scored"] else "")
                  + f"\n           A sterner floor sits at "
                    f"{fl['p75']:.4f} (upper quartile) or "
                    f"{fl['p90']:.4f} (upper decile).")
        print(f"{'=' * 72}")
        if t_ in ("below", "unknown"):
            continue
        # NAMED FROM THE SAME MAP THE ROLES CAME FROM. 15 rows in the
        # committed artifact read "resembles NMDYYWFGPIMTKO", because that
        # reference node carried no chembl_id -- and a reference nobody
        # can look up defeats the one audit this read-out offers. The
        # InChIKey route that fixed the role gap fixes the naming gap with
        # it, on both sides of the row.
        for r in rows:
            r["name"] = (ann.get(r["candidate"]) or {}).get("name")
            r["nearest_active_name"] = (
                (ann.get(r["nearest_active"]) or {}).get("name"))
        print(f"  {'#':>2}  {'compound':<30} {'tier':<9} {'Tc':>6} "
              f"{'>fail':>6} {'deg':>5}  resembles")
        print("      >fail = the share of this virus's screened-and-INACTIVE "
              "compounds whose own\n      best resemblance to the known "
              "actives is WEAKER than this row's. The floor is the 50% "
              "line.")
        last_tier = None
        for n, r in enumerate(rows, 1):
            if r["phase_tier"] != last_tier:
                last_tier = r["phase_tier"]
                print(f"  -- {_TIER_BANNER[last_tier]}")
            cid = r["chembl_id"] or r["candidate"].split(":", 1)[1][:14]
            shown = r["name"] or cid
            near = (r["nearest_active_name"] or r["nearest_active_chembl"]
                    or r["nearest_active"].split(":", 1)[1][:14])
            # A 1.0 IS MARKED, WITH WHAT ACTUALLY DIFFERS. Unmarked it
            # reads as an identity, and 14 of these are a ring homolog or a
            # peptide variant that ECFP4 at radius 2 cannot resolve.
            notes = []
            if r.get("fingerprint_identical"):
                d = r.get("heavy_atom_delta")
                notes.append(
                    "SAME BITS, "
                    + (f"{abs(d)} heavy atom(s) apart"
                       if d else "same atom count: a permutation")
                    + "; NOT the same compound")
            # WHAT IT IS FOR, BESIDE WHAT IT LOOKS LIKE. An unmarked
            # GEMCITABINE under "orderable now" is the row a referee
            # quotes: it resembles a nucleoside antiviral because it IS a
            # nucleoside, and it lowers a cell-based viral readout by
            # killing the cell.
            notes += [_MARK[m] for m in r.get("flags") or () if m in _MARK]
            # ONLY THE CONFOUND TRAVELS FROM THE REFERENCE. A black-box
            # warning on the compound a suggestion RESEMBLES says nothing
            # about the suggestion, and printing it on every
            # prochlorperazine analogue buried the one marker that does
            # matter: that the thing it resembles is a cytotoxic.
            nearnote = "".join(
                f" [{_MARK[m]}]" for m in r.get("nearest_active_flags") or ()
                if m == ROLE_CYTOTOXIC)
            bf = r.get("beats_failures")
            print(f"  {n:>2}  {shown[:30]:<30} "
                  f"{r['phase_tier']:<9}"
                  f" {r['tanimoto']:>6.3f} "
                  f"{'—' if bf is None else format(bf, '.0%'):>6} "
                  f"{r['graph_degree']:>5,}  "
                  f"{near[:28]}{nearnote}"
                  + (f"  <- {'; '.join(notes)}" if notes else ""))
        if rows:
            n_appr = sum(1 for r in rows if r["is_approved"])
            n_reach = sum(1 for r in rows if r["graph_degree"])
            by = {t: sum(1 for r in rows if r["phase_tier"] == t)
                  for t in TIERS}
            print(f"\n  {len(rows)} shown: "
                  + ", ".join(f"{by[t]} {t}" for t in TIERS) + ".")
            # WHAT THE GRAPH ADDS, COUNTED RATHER THAN IMPLIED. Until
            # 2026-11-19 this line read the held-out labels and reported
            # every candidate as reached; the true figure on the committed
            # artifact was 117 of 195. A suggestion with degree 0 rests on
            # chemistry alone and the graph has nothing to add about it.
            print(f"  {n_appr} are approved drugs; {n_reach} of {len(rows)} "
                  f"carry any non-label graph edge at all.")
            n_cy = sum(1 for r in rows if ROLE_CYTOTOXIC in (r["flags"] or []))
            n_wd = sum(1 for r in rows if ROLE_WITHDRAWN in (r["flags"] or []))
            if n_cy or n_wd:
                print(f"  {n_cy} {_is(n_cy)} antineoplastic(s) (ATC "
                      f"{CYTOTOXIC_ATC}) and {n_wd} {_is(n_wd)} withdrawn "
                      f"in some jurisdiction.\n  Both are shown on "
                      f"purpose: the first is an antiviral class as well "
                      f"as a cytotoxic\n  one, and the second flag reads "
                      f"1 for CODEINE and THALIDOMIDE, so it cannot be a "
                      f"filter.")
            n_same = sum(1 for r in rows if r.get("fingerprint_identical"))
            # PER TIER. rows is now sorted by (tier, -Tc), so rows[0]
            # and rows[-1] are the top of the first tier and the bottom
            # of the last -- which printed "Highest Tc 0.723, lowest
            # shown 0.925" the first time this ran. A summary line that
            # contradicts the table above it is worse than no summary.
            for t in TIERS:
                tc = [r["tanimoto"] for r in rows if r["phase_tier"] == t]
                if tc:
                    print(f"  {t:<9} Tc {max(tc):.3f} down to {min(tc):.3f}")
            if n_same:
                print(f"  {n_same} of them {_does(n_same)} share the "
                      f"active's ECFP4 bits EXACTLY while being a "
                      f"different structure -- read those as analogs, not "
                      f"identities.")

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
