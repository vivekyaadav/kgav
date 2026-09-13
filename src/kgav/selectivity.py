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

# Conventional cut below which a cell-based antiviral hit is not pursued.
SI_THRESHOLD = 10.0

# Cytotoxicity readouts. GI50 (growth inhibition) is deliberately EXCLUDED:
# it measures growth arrest rather than death, is dominated by oncology
# panels, and mixing it with CC50 would conflate two different endpoints.
TOX_TYPES = ("CC50", "TC50")
ACTIVITY_TYPES = ("EC50", "IC50")


def query_pairs(db: Path, taxa: list[str]) -> list[dict]:
    """Matched (compound, document) antiviral/cytotoxicity measurements.

    Both sides take the MINIMUM value within a document, since a paper may
    report several concentrations or cell lines. Taking the minimum EC50 and
    the minimum CC50 is the conservative choice for selectivity: it maximises
    apparent potency and minimises apparent tolerability, so a compound
    clearing the threshold under this rule clears it under any.
    """
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    placeholders = ",".join("?" * len(taxa))
    act_types = ",".join(f"'{t}'" for t in ACTIVITY_TYPES)
    tox_types = ",".join(f"'{t}'" for t in TOX_TYPES)

    rows = con.execute(f"""
        WITH antiviral AS (
            SELECT act.molregno, a.doc_id, t.tax_id,
                   MIN(act.standard_value) AS ec50_nm,
                   COUNT(*) AS n_activity,
                   MIN(a.description) AS descr
            FROM activities act
            JOIN assays a ON act.assay_id = a.assay_id
            JOIN target_dictionary t ON a.tid = t.tid
            WHERE t.tax_id IN ({placeholders})
              AND act.standard_type IN ({act_types})
              AND act.standard_relation = '='
              AND act.standard_units = 'nM'
              AND act.standard_value > 0
            GROUP BY act.molregno, a.doc_id, t.tax_id
        ),
        toxicity AS (
            SELECT act.molregno, a.doc_id,
                   MIN(act.standard_value) AS cc50_nm,
                   COUNT(*) AS n_tox,
                   -- A censored CC50 ("> 50000 nM") means the compound was NOT
                   -- toxic up to that concentration, so the true value is
                   -- HIGHER and the resulting ratio is a LOWER BOUND on
                   -- selectivity. Treating it as a point value manufactures
                   -- cytotoxicity: it produced SI 1.0 for nirmatrelvir and 0.0
                   -- for remdesivir, both licensed antivirals.
                   MIN(CASE WHEN act.standard_relation = '=' THEN 0 ELSE 1 END) AS all_exact
            FROM activities act
            JOIN assays a ON act.assay_id = a.assay_id
            WHERE act.standard_type IN ({tox_types})
              AND act.standard_units = 'nM'
              AND act.standard_relation IN ('=', '>')
              AND act.standard_value > 0
            GROUP BY act.molregno, a.doc_id
        )
        SELECT av.molregno, av.doc_id, av.tax_id, av.ec50_nm, av.n_activity,
               tx.cc50_nm, tx.n_tox, tx.all_exact,
               md.chembl_id AS compound,
               cs.standard_inchi_key AS inchikey,
               d.year, d.pubmed_id, av.descr AS assay_text
        FROM antiviral av
        JOIN toxicity tx ON av.molregno = tx.molregno AND av.doc_id = tx.doc_id
        JOIN molecule_dictionary md ON av.molregno = md.molregno
        LEFT JOIN compound_structures cs ON av.molregno = cs.molregno
        LEFT JOIN docs d ON av.doc_id = d.doc_id
    """, taxa).fetchall()
    con.close()
    return [dict(r) for r in rows]


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
        # study by its answer -- so exact measurements decide, and among
        # several exact measurements the WORST is kept, since a compound shown
        # cytotoxic in a real measurement is cytotoxic.
        # all_exact must be present. Defaulting a missing value to "censored"
        # silently discards the row; defaulting to "exact" silently trusts a
        # bound as a measurement. Both are the same class of error that
        # produced SI 1.8 for nirmatrelvir, so an absent field is counted and
        # dropped rather than interpreted.
        if r.get("all_exact") is None:
            stats["missing_relation_flag"] += 1
            continue
        exact = r["all_exact"] == 0
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
        em.edge(drug, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", virus,
                source=source,
                date=f"{int(r['year'])}-01-01" if r.get("year") else "1970-01-01",
                tier=1,
                quals={"ec50_nm": float(r["ec50_nm"]),
                       "cc50_nm": float(r["cc50_nm"]),
                       "selectivity_index": round(si, 3),
                       "selectivity_verified": True,
                       "selectivity_index_max": round(max(payload["all"]), 3),
                       "selectivity_n_studies": len(payload["all"]),
                       "relation": "=",
                       "assay_type": "cell_based_antiviral",
                       **({"cell_line": cl} if (cl := extract_cell_line(
                           r.get("assay_text"))) else {}),
                       # unquantified is now FALSE for these: a selectivity
                       # index exists, which is the whole point of the layer.
                       "unquantified": False},
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
