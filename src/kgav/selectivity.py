"""Selectivity index from same-document EC50/CC50 pairs.

CORRECTING AN EARLIER ERROR. This project previously recorded that ChEMBL
contains almost no cytotoxicity data for coronaviruses -- 2 CC50 records across
21,900 activities. That was wrong, and the error was in the query rather than
the database: cytotoxicity assays are filed against the CELL LINE they were run
in, not against the virus, so a search restricted to coronavirus organism
targets finds nothing. ChEMBL 37 holds 74,450 CC50 records under cell-line
targets and a further 10,780 under unchecked targets.

WHY SAME-DOCUMENT PAIRING. Selectivity index is CC50 / EC50, and the ratio is
only meaningful if both numbers come from the same experimental context. A CC50
measured in HeLa by one group and an EC50 measured in Calu-3 by another
produces a number with the right units and no interpretation. Pairing on
(compound, document) gives both measurements from one paper, which in practice
means one cell system, one protocol, one set of controls.

Cross-document pairing would yield far more pairs and would be much weaker. It
is not done here; compounds without a same-document pair are marked
`selectivity_verified: false` rather than being given a fabricated ratio.

WHAT THE PAIRING FOUND. 1,634 compounds with matched data across 1,819
(compound, document) pairs. Of those pairs, 691 -- 38% -- have a selectivity
index below 10, the conventional threshold below which a cell-based antiviral
hit is not pursued. Roughly two in five compounds that look active in a
coronavirus cell assay are cytotoxic enough that the activity is suspect.

CONSEQUENCE FOR THE EVALUATION. The positive set used in every earlier
evaluation treats those compounds as actives. Re-running the hard-negative
comparison with them excluded is a real test of whether the direct-acting
signal survives, and it can fail: if the signal partly reflects the graph
recognising well-characterised compounds that happen to be cytotoxic, the AUC
will fall.
"""
from __future__ import annotations

import sqlite3
from collections import Counter
from pathlib import Path

# ONE UNIT IMPLEMENTATION AND ONE ALLOW-SET. units.py:17 names this file as
# a place the conversion was missing entirely -- the SQL demanded literal
# 'nM' -- and warns that a second copy is the worst-placed instance of a
# pattern this project has found ten times.
from kgav.chembl import CHEMBL_UNITS
from kgav.temporal import PLACEHOLDER_DATE
from kgav.units import to_nm

# Conventional cut below which a cell-based antiviral hit is not pursued.
SI_THRESHOLD = 10.0

# Cytotoxicity readouts. GI50 (growth inhibition) is deliberately EXCLUDED:
# it measures growth arrest rather than death, is dominated by oncology
# panels, and mixing it with CC50 would conflate two different endpoints.
TOX_TYPES = ("CC50", "TC50")
ACTIVITY_TYPES = ("EC50", "IC50")


#: Units this layer will convert, the same frozen set chembl and
#: host_targets use. units.py says why pinning it matters: adding a unit to
#: CONCENTRATION_FACTORS must not silently change what an existing ingest
#: accepts.
SELECTIVITY_UNITS = CHEMBL_UNITS


def query_pairs(db: Path, taxa: list[str],
                stats: Counter | None = None) -> list[dict]:
    """Matched (compound, document) antiviral/cytotoxicity measurements.

    EVERY GATE IS COUNTED, which is why the grouping is here rather than in
    SQL. The previous version did it in four CTE WHERE clauses, so a row
    outside IC50/EC50, a censored EC50, and -- the one that mattered -- any
    row whose units were not LITERALLY 'nM' vanished with nothing recorded.
    Fix-list item E7 is that this layer, hosttargets and ncats have no
    `emitted + dropped == rows` accounting, which only ChEMBL has; a gate
    inside a WHERE clause cannot be given one without moving it out.

    UNITS GO THROUGH kgav.units.to_nm. They were not converted at all: the
    SQL demanded 'nM'. On the antiviral side that costs almost nothing --
    ChEMBL normalises potencies, so 37,776 of 38,023 rows are already nM
    and 240 of the remaining 247 are ug.mL-1, which needs a molecular
    weight and is genuinely unconvertible. On the TOXICITY side it dropped
    4,327 uM and 116 mM rows, worth 72 extra (compound, document) pairs.
    A CASE expression in SQL would have been a fourth implementation of the
    conversion, and units.py calls that "the worst placed instance of a
    pattern this project has now found ten times".

    THE RETAINED CC50 AND ITS exact FLAG DESCRIBE THE SAME ROW. The flag
    was MIN(CASE WHEN relation = '=' THEN 0 ELSE 1 END) -- a minimum, so it
    read "at least one row is exact" while being named `all_exact` -- and
    cc50_nm was MIN over every row including censored ones. 51 of the 3,388
    SI-contributing groups therefore reported an exact CC50 taken from a
    censored row, and in all 51 the real exact value is HIGHER, so the
    index was understated and the compound looked more cytotoxic than its
    own exact measurement says. When any exact row exists the CC50 is now
    the minimum over the exact rows only.
    """
    counts: Counter = stats if stats is not None else Counter()
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    placeholders = ",".join("?" * len(taxa))
    act_types = ",".join(f"'{t}'" for t in ACTIVITY_TYPES)
    tox_types = ",".join(f"'{t}'" for t in TOX_TYPES)

    # ---------------------------------------------------- antiviral side
    av_rows = con.execute(f"""
        SELECT act.molregno, a.doc_id, t.tax_id,
               act.standard_value AS value, act.standard_units AS units,
               act.standard_relation AS relation, a.description AS descr
        FROM activities act
        JOIN assays a ON act.assay_id = a.assay_id
        JOIN target_dictionary t ON a.tid = t.tid
        WHERE t.tax_id IN ({placeholders})
          AND act.standard_type IN ({act_types})
    """, taxa).fetchall()

    antiviral: dict[tuple, dict] = {}
    for r in av_rows:
        counts["av_rows"] += 1
        if r["relation"] != "=":
            counts["av_dropped_censored_or_null_relation"] += 1
            continue
        nm = to_nm(r["value"], r["units"], allow=SELECTIVITY_UNITS)
        if nm is None or nm <= 0:
            counts["av_dropped_unconvertible_units"] += 1
            counts[f"av_unit:{r['units']}"] += 1
            continue
        if r["units"] != "nM":
            counts["av_converted"] += 1
        # ROWS KEPT, separately from the groups they collapse into. The
        # accounting identity is dropped + kept == read; groups are fewer
        # than kept rows by construction, and using them closed nothing.
        counts["av_kept"] += 1
        key = (r["molregno"], r["doc_id"], r["tax_id"])
        cur = antiviral.get(key)
        # Tie-break on the description so two runs over one database cannot
        # disagree about which row was retained.
        cand = (nm, r["descr"] or "")
        if cur is None:
            antiviral[key] = {"ec50_nm": nm, "descr": r["descr"],
                              "n_activity": 1, "_k": cand}
        else:
            cur["n_activity"] += 1
            if cand < cur["_k"]:
                cur.update(ec50_nm=nm, descr=r["descr"], _k=cand)
    counts["av_groups"] = len(antiviral)

    # ---------------------------------------------------- toxicity side
    tox_rows = con.execute(f"""
        SELECT act.molregno, a.doc_id,
               act.standard_value AS value, act.standard_units AS units,
               act.standard_relation AS relation, a.description AS descr
        FROM activities act
        JOIN assays a ON act.assay_id = a.assay_id
        WHERE act.standard_type IN ({tox_types})
    """).fetchall()

    toxicity: dict[tuple, dict] = {}
    for r in tox_rows:
        counts["tox_rows"] += 1
        if r["relation"] not in ("=", ">"):
            counts["tox_dropped_relation"] += 1
            continue
        nm = to_nm(r["value"], r["units"], allow=SELECTIVITY_UNITS)
        if nm is None or nm <= 0:
            counts["tox_dropped_unconvertible_units"] += 1
            counts[f"tox_unit:{r['units']}"] += 1
            continue
        if r["units"] != "nM":
            counts["tox_converted"] += 1
        counts["tox_kept"] += 1
        key = (r["molregno"], r["doc_id"])
        exact = r["relation"] == "="
        rec = toxicity.setdefault(key, {"n_tox": 0, "exact": [], "bound": []})
        rec["n_tox"] += 1
        rec["exact" if exact else "bound"].append((nm, r["descr"] or ""))
    counts["tox_groups"] = len(toxicity)

    # ---------------------------------------------------- the join
    meta = {}
    for row in con.execute("""
        SELECT md.molregno, md.chembl_id AS compound,
               cs.standard_inchi_key AS inchikey
        FROM molecule_dictionary md
        LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
    """):
        meta[row["molregno"]] = (row["compound"], row["inchikey"])
    docs = {row["doc_id"]: (row["year"], row["pubmed_id"])
            for row in con.execute("SELECT doc_id, year, pubmed_id FROM docs")}
    con.close()

    out: list[dict] = []
    for (mol, doc, tax), av in antiviral.items():
        tx = toxicity.get((mol, doc))
        if tx is None:
            counts["no_matched_counterscreen"] += 1
            continue
        # THE CC50 AND ITS FLAG FROM ONE ROW. Exact rows win outright: a
        # censored bound and an exact measurement are different claims and
        # the minimum across both belongs to neither.
        pool = tx["exact"] or tx["bound"]
        cc50_nm, tox_descr = min(pool)
        cc50_exact = bool(tx["exact"])
        if tx["exact"] and tx["bound"]:
            counts["tox_mixed_group_exact_preferred"] += 1
        compound, inchikey = meta.get(mol, (None, None))
        year, pubmed_id = docs.get(doc, (None, None))
        out.append({"molregno": mol, "doc_id": doc, "tax_id": tax,
                    "ec50_nm": av["ec50_nm"], "n_activity": av["n_activity"],
                    "cc50_nm": cc50_nm, "n_tox": tx["n_tox"],
                    "cc50_exact": cc50_exact, "tox_descr": tox_descr,
                    "compound": compound, "inchikey": inchikey,
                    "year": year, "pubmed_id": pubmed_id,
                    "assay_text": av["descr"]})
    counts["pairs"] = len(out)
    return out


# Cell lines appearing in coronavirus assays, longest patterns first so
# "Vero E6" is not swallowed by "Vero".
CELL_LINE_PATTERNS = [
    ("Vero E6", ("vero e6", "veroe6")),
    ("Vero", ("vero",)),
    ("Calu-3", ("calu-3", "calu3")),
    ("Caco-2", ("caco-2", "caco2")),
    ("Huh-7.5", ("huh-7.5", "huh7.5")),
    ("Huh-7", ("huh-7", "huh7")),
    ("A549", ("a549",)),
    ("HEK293T", ("hek293t", "hek 293t")),
    ("HEK293", ("hek293", "hek 293")),
    ("BHK-21", ("bhk-21", "bhk21")),
    ("MRC-5", ("mrc5", "mrc-5")),
    ("HEL 299", ("hel 299",)),
    ("HeLa", ("hela",)),
    ("primary airway", ("air-liquid", "primary human airway", "hae ")),
]


def extract_cell_line(description: str | None) -> str | None:
    """Cell line named in an assay description.

    The same compound's EC50 varies thirtyfold across cell lines in this data
    -- chloroquine ranges 330 nM in BHK-21 to 10,900 nM in a Vero line -- so a
    potency or selectivity figure without its cell line is not interpretable.
    """
    d = (description or "").lower()
    for label, patterns in CELL_LINE_PATTERNS:
        if any(pat in d for pat in patterns):
            return label
    return None


def selectivity_index(ec50_nm: float, cc50_nm: float) -> float | None:
    if not ec50_nm or ec50_nm <= 0 or not cc50_nm or cc50_nm <= 0:
        return None
    return cc50_nm / ec50_nm


def ingest_selectivity(em, rows: list[dict], taxon_map: dict[str, str],
                       source: str, threshold: float = SI_THRESHOLD) -> Counter:
    """Emit activity edges carrying CC50 and a computed selectivity index.

    These share (subject, predicate, object) with the chemistry layer's
    activity edges, so assembly merges them and the selectivity qualifiers fill
    gaps the chemistry layer left empty.
    """
    stats: Counter = Counter()
    best: dict[tuple[str, str], dict] = {}

    for r in rows:
        stats["rows"] += 1
        inchikey = r.get("inchikey")
        if not inchikey:
            stats["no_inchikey"] += 1
            continue
        tax = taxon_map.get(str(r["tax_id"]))
        if tax is None:
            stats["unmapped_taxon"] += 1
            continue
        si = selectivity_index(r["ec50_nm"], r["cc50_nm"])
        if si is None:
            stats["si_uncomputable"] += 1
            continue

        key = (f"INCHIKEY:{inchikey}", f"NCBITaxon:{tax}")
        # A compound may be paired in several papers. Keep the WORST
        # selectivity, not the best: a compound shown cytotoxic anywhere is
        # cytotoxic, and choosing the flattering measurement is how a
        # selectivity filter stops filtering.
        # Keep the worst for filtering, but retain the spread: reporting
        # nirmatrelvir as cytotoxic at SI 1.8 from one unflattering paper,
        # with no sign that other studies disagree, states something false
        # about a licensed drug.
        # Keep the BEST-EVIDENCED pair, not the worst. The worst-case rule was
        # right for a filter and is wrong here: with a censored CC50 a low
        # ratio is an uninformative lower bound, not a cytotoxicity finding,
        # so minimising over them systematically defames selective compounds.
        # Exact measurements outrank censored ones; among equals the higher
        # ratio wins.
        # A censored CC50 is a LOWER BOUND: it can establish that a compound
        # is selective, never that it is cytotoxic. So censored pairs are kept
        # only when the bound already clears the threshold; below it they are
        # uninformative and must not vote either way.
        #
        # Taking the worst ratio defamed nirmatrelvir (SI 1.8) and remdesivir
        # (0.9). Taking the best exonerated chloroquine (17.7) and
        # hydroxychloroquine (21.2). Both were the same error -- selecting a
        # study by its answer -- so exact measurements decide.
        #
        # THIS COMMENT USED TO END "and among several exact measurements the
        # WORST is kept, since a compound shown cytotoxic in a real
        # measurement is cytotoxic." That was superseded in the same function
        # and never corrected here: the emit loop below takes the MEDIAN of
        # the exact values, for the reason stated there -- the worst of 34
        # remdesivir measurements is itself an outlier. The line immediately
        # under this one has said "collect all; median chosen below" the
        # whole time, and tests/test_selectivity.py pins the median.
        #
        # It is corrected rather than deleted because the stale clause was
        # CITED as authority elsewhere -- baselines.potency_tiebreak took
        # "keep the worst" from it -- and a comment that outlived its code
        # does its damage by being quoted.
        # cc50_exact must be present. Defaulting a missing value to "censored"
        # silently discards the row; defaulting to "exact" silently trusts a
        # bound as a measurement. Both are the same class of error that
        # produced SI 1.8 for nirmatrelvir, so an absent field is counted and
        # dropped rather than interpreted.
        if r.get("cc50_exact") is None:
            stats["missing_relation_flag"] += 1
            continue
        exact = bool(r["cc50_exact"])
        if not exact and si < threshold:
            stats["censored_uninformative"] += 1
            continue
        prior = best.get(key)
        if prior is None:
            better = True
        elif exact and not prior["exact"]:
            better = True                       # exact beats censored
        elif exact and prior["exact"]:
            better = False                      # collect all; median chosen below
        elif not exact and not prior["exact"]:
            better = si > prior["si"]           # strongest bound wins
        else:
            better = False                      # never demote exact to censored
        if exact:
            if key in best:
                best[key].setdefault("exact_all", []).append(si)
            else:
                best[key] = {"si": si, "row": r, "exact": True,
                             "all": [si], "exact_all": [si]}
                continue
        if better:
            best[key] = {"si": si, "row": r, "exact": exact,
                         "all": (prior["all"] if prior else []) + [si],
                         "exact_all": (prior.get("exact_all", []) if prior else [])}
        else:
            best[key]["all"].append(si)

    for (drug, virus), payload in best.items():
        r = payload["row"]
        # MEDIAN across exact measurements, not the minimum. Taking the worst
        # of 34 remdesivir measurements selected an outlier -- one paper
        # reporting CC50 168 nM against EC50 55.6 nM -- and reported a
        # licensed antiviral as cytotoxic at SI 0.9. With one measurement the
        # minimum IS the measurement; with 34 it is selection on the outcome.
        exact_vals = sorted(payload.get("exact_all") or [])
        if exact_vals:
            mid = len(exact_vals) // 2
            si = (exact_vals[mid] if len(exact_vals) % 2
                  else (exact_vals[mid - 1] + exact_vals[mid]) / 2)
        else:
            si = payload["si"]
        selective = si >= threshold
        stats["selective" if selective else "cytotoxic"] += 1
        # THE PLACEHOLDER, COUNTED (fix-list C10).
        stats["dated_year_only" if r.get("year")
              else "undated_placeholder"] += 1

        # ONE EXPERIMENTAL CONTEXT, WHICH A SHARED DOCUMENT IS NOT.
        #
        # selectivity_verified was written True on every edge, and the
        # schema declares it means "CC50 and EC50 come from the SAME
        # document, so the selectivity index reflects one experimental
        # context". The join does guarantee the document. A paper routinely
        # measures antiviral activity in one cell line and cytotoxicity in
        # another, and then the ratio is not a selectivity index -- it is a
        # potency in Calu-3 divided by a toxicity in Vero.
        #
        # Measured on chembl_37, over the 4,167 (compound, document, taxon)
        # pairs: 2,244 agree, 337 name DIFFERENT cell lines, and 1,586 have
        # at least one side unnamed. Only the first group can be verified.
        # The other two are not errors in the source; they are the limit of
        # what a document-level join can establish, and the flag now says so
        # rather than asserting the strongest reading of all three.
        ec_cl = extract_cell_line(r.get("assay_text"))
        tox_cl = extract_cell_line(r.get("tox_descr"))
        if ec_cl and tox_cl:
            verified = ec_cl == tox_cl
            stats["context_matched" if verified
                  else "context_cell_line_mismatch"] += 1
        else:
            verified = False
            stats["context_cell_line_unnamed"] += 1
        em.edge(drug, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", virus,
                source=source,
                # A YEAR IS WRITTEN AS A YEAR (fix-list C10).
                date=(str(int(r["year"])) if r.get("year")
                      else PLACEHOLDER_DATE),
                tier=1,
                quals={"document_id": str(r["doc_id"]),
                       "ec50_nm": float(r["ec50_nm"]),
                       "cc50_nm": float(r["cc50_nm"]),
                       "selectivity_index": round(si, 3),
                       "selectivity_verified": verified,
                       "selectivity_index_max": round(max(payload["all"]), 3),
                       "selectivity_n_studies": len(payload["all"]),
                       "relation": "=",
                       "assay_type": "cell_based_antiviral",
                       # The EC50's OWN cell line. Before the query used a
                       # window function this came from whichever row sorted
                       # first alphabetically, which named a different cell
                       # line than the potency on 794 groups.
                       **({"cell_line": ec_cl} if ec_cl else {}),
                       # `unquantified` is NOT set here either. Writing false
                       # to mean "a selectivity index exists" was a third
                       # definition of the field, and it lost the merge to
                       # chembl's true on the 1,264 pairs that had both. The
                       # cc50_nm and ec50_nm above are the evidence; the flag
                       # is derived from them at assembly.
                       },
                pmids=[f"PMID:{r['pubmed_id']}"] if r.get("pubmed_id") else None)
        stats["edges"] += 1
    return stats


def selectivity_index_by_compound(release: Path) -> dict[tuple[str, str], float]:
    """(compound, virus) -> selectivity index, for filtering an evaluation."""
    import json
    out: dict[tuple[str, str], float] = {}
    for line in (Path(release) / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] != "HAS_ANTIVIRAL_ACTIVITY_AGAINST":
            continue
        si = (e.get("qualifiers") or {}).get("selectivity_index")
        if si is not None:
            out[(e["subject"], e["object"])] = float(si)
    return out
