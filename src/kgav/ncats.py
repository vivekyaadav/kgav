"""NCATS OpenData SARS-CoV-2 CPE screen: labels in the population that matters.

    SmallMolecule -HAS_ANTIVIRAL_ACTIVITY_AGAINST-> OrganismTaxon

WHY THIS LAYER EXISTS. Four evaluation protocols have now scored at chance or
below once their confounds were controlled (docs/evaluation.md). Two of the six
measured causes are properties of the LABEL SET, not of the graph:

  * Negatives cannot reach the host-directed metapaths. 1.5% of
    MEASURED_INACTIVE compounds reach a dependency host factor against 12.1%
    expected -- 7.9x under-represented. ChEMBL negatives are antiviral research
    compounds, which carry viral-target annotations and almost no host ones.
  * 94.9% of compounds carrying both an INHIBITS edge and an antiviral label
    cite the same PMID on both, because one paper reports the target IC50 and
    the cell-based EC50 together.

This panel fixes both. It reports what did NOT work -- 9,848 of 10,633 rows are
curve class 4, tested and inactive -- drawn from approved and annotated
libraries, which is the population holding the TARGETS edges. It is cell-based
and generated independently of ChEMBL target annotation, so no binding edge
cites its paper and --publication-disjoint withholds nothing when scoring
against it.

TWO FILES, NOT ONE. The portal's root cpe.tsv pools three different assays
(SARS-CoV-2 cytopathic effect 5,862 rows; SARS-CoV-2_CPE_SRI 4,722;
DPI_SARS-CoV-2_General 49) and carries no cytotoxicity column. The canonical
export is the pair under /public/odp/assay/:

    SARS-CoV-2_Cytopathic_Effect_(CPE).csv                  activity
    SARS-CoV-2_Cytopathic_Effect_(Host_Tox_Counterscreen).csv  cytotoxicity

Same Vero E6 protocol without virus, same sample_ids, same concentrations. The
pair is what makes the selectivity index same-plate rather than same-document.

THE POLARITY IS THE DANGEROUS PART, and it is the BioGRID ORCS problem in a
new costume. CPE is a GAIN-of-signal assay: a compound that protects cells
raises the readout. A compound that kills them lowers it. Both produce a fitted
AC50, and nothing in the AC50 says which happened -- the sign is carried by
`efficacy`. The first row of the real CPE export is ac50 7.08 uM with efficacy
-46.9: a cytotoxic compound that an AC50-only rule would file as an antiviral
with EC50 7 uM. So activity requires a POSITIVE efficacy of at least
MIN_EFFICACY_PCT, never a fitted value alone.

A FITTED AC50 IS ALSO NOT A CALL. 785 of 10,633 rows carry an AC50 but only
~190 sit in a real curve class; the other 586 are class 3, single-point
activity, which NCATS treats as inconclusive. Class is read first, value
second.

IDENTITY IS THE OTHER RISK. The graph's compound keys came out of
normalize_chemical(); a key taken from this file instead would differ on salt,
stereo or tautomer and the new nodes would sit beside the existing ones without
joining to a single TARGETS edge. Nothing would error and no count would look
wrong -- the README calls a discarded identifier the dominant bug class, found
four times, silent every time.
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

ANTIVIRAL_PREDICATE = "HAS_ANTIVIRAL_ACTIVITY_AGAINST"

# Chen CZ et al., Front Pharmacol 2021;11:592737. Vero E6, MOI 0.002, 72h.
CPE_PUBLICATION = "PMID:33708112"
SOURCE = "infores:ncats-opendata"
# Vero E6, the line this screen used. Cellosaurus: Vero E6 is CVCL_0574, which
# is NOT Vero (CVCL_0059) -- a distinct subclone selected for ACE2 expression.
VERO_E6 = "CVCL:0574"
SARS_COV_2 = "NCBITaxon:2697049"

# First public assertion: bioRxiv 2020.08.18.255877. Not the portal's file
# mtime (2023-05-22), which is a re-export date and would put a 2020 screen
# after both PPI layers. Overridable; never inferred.
DEFAULT_DATE = "2020-08-18"

# Only the canonical assay. cpe.tsv pools three, and SRI ran a different
# protocol: merging them would average two concentration ranges into one
# censoring point.
CANONICAL_ASSAY = "sars-cov-2 cytopathic effect (cpe)"
COUNTERSCREEN_ASSAY = "sars-cov-2 cytopathic effect (host tox counterscreen)"

# NCATS qHTS curve classes (Inglese 2006). 1.x complete curve, 2.x incomplete
# with one asymptote -- both are real concentration-response fits. 3 is
# single-point activity at the top concentration only, 5 is a poor fit: both
# inconclusive, and 586 + n rows here are class 3, so admitting them would
# quadruple the active set with unreplicated single points. 4 is inactive, and
# it is the class this layer is here for.
CURVE_CLASS_FITTED = frozenset({"1.1", "1.2", "1.3", "1.4",
                                "2.1", "2.2", "2.3", "2.4"})
CURVE_CLASS_INACTIVE = frozenset({"4"})

# Chen et al. selected hits at >55% efficacy. 50 is the round floor below that;
# raise it with --min-efficacy-pct to reproduce the paper's hit list exactly.
MIN_EFFICACY_PCT = 50.0

# A compound is only callable INACTIVE if the assay reached a concentration at
# or above the threshold the rest of the pipeline uses for inactivity. Tested
# to 5 uM and found inactive says nothing about 10 uM, and labels.classify()
# would silently return None for it; counting it as a negative anyway would
# invent evidence.
DEFAULT_INACTIVE_AT_NM = 10_000.0


class ColumnError(ValueError):
    """Raised when the header cannot be mapped. Carries the whole header.

    Deliberately fatal. Guessing a column here produces a release that
    validates, counts plausibly, and encodes the wrong measurement.
    """


# field -> header spellings seen across NCATS ODP exports. The assay/ CSVs are
# lowercase and the root TSV is uppercase, so _key() flattens both. Order
# matters: sample_id is the stable NCGC identifier, sample_name is a label and
# is sometimes the identifier repeated.
COLUMN_CANDIDATES: dict[str, tuple[str, ...]] = {
    "smiles": ("smiles", "structure", "canonicalsmiles"),
    "sample": ("sampleid", "ncgcsid", "sid", "ncgcid", "samplename", "sample"),
    "assay": ("assayname", "assay"),
    "library": ("library",),
    "ac50": ("ac50", "ac50um", "ac50m", "ac50nm"),
    "lac50": ("logac50", "lac50"),
    "efficacy": ("efficacy", "efficacypct"),
    "curve": ("curveclass2", "curveclass", "cclass"),
    "max_response": ("maxresponse",),
    "r2": ("r2",),
    "moa": ("primarymoa", "moa"),
}


def _key(s: str) -> str:
    return "".join(c for c in s.strip().lower() if c.isalnum())


def resolve_columns(header: list[str],
                    required: tuple[str, ...] = ("smiles", "curve", "efficacy"),
                    ) -> dict:
    """Map our field names onto this file's actual header.

    `smiles`, `curve` and `efficacy` are all required, not just the structure:
    without the curve class there is no activity call, and without efficacy
    there is no direction -- and a direction guessed in a gain-of-signal assay
    files cytotoxic compounds as antivirals.

    Also collects conc_cols, in file order, so the top concentration is read
    from the row rather than assumed.
    """
    index = {_key(h): h for h in header}
    found: dict = {}
    for field, cands in COLUMN_CANDIDATES.items():
        for c in cands:
            if _key(c) in index:
                found[field] = index[_key(c)]
                break
    found["conc_cols"] = [h for h in header if _key(h).startswith("conc")]
    missing = [f for f in required if f not in found]
    if missing:
        raise ColumnError(
            f"could not find {missing} in header.\n"
            f"  header seen ({len(header)} columns): {header}\n"
            f"  add the real spelling to COLUMN_CANDIDATES in src/kgav/ncats.py "
            f"rather than renaming the file's columns, so the next export of "
            f"this dataset keeps working."
        )
    return found


def read_table(path: Path) -> tuple[list[dict], dict]:
    """Read a .tsv or .csv export and resolve its header."""
    p = Path(path)
    delim = "\t" if p.suffix.lower() in (".tsv", ".tab") else ","
    with p.open(newline="", encoding="utf-8-sig") as fh:
        rdr = csv.DictReader(fh, delimiter=delim)
        if not rdr.fieldnames:
            raise ColumnError(f"{p} has no header row")
        cols = resolve_columns(list(rdr.fieldnames))
        return list(rdr), cols


def _f(value) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def to_nm(value, unit: str) -> float | None:
    """Concentration to nanomolar. `unit` is one of M, uM, nM, logM.

    logM is NCATS's log_ac50: log10 of molar. -5.0 is 10 uM, i.e. 10,000 nM.
    """
    v = _f(value)
    if v is None:
        return None
    if unit == "nM":
        return v
    if unit == "uM":
        return v * 1_000.0
    if unit == "M":
        return v * 1e9
    if unit == "logM":
        return (10.0 ** v) * 1e9
    raise ValueError(f"unknown concentration unit {unit!r}")


def unit_for(column_name: str, default: str = "uM") -> str:
    """Infer the unit from the column's own name, since NCATS varies it."""
    k = _key(column_name)
    if k.startswith(("lac50", "logac50")):
        return "logM"
    if k.endswith("nm"):
        return "nM"
    if k.endswith("um"):
        return "uM"
    if k.endswith("m"):
        return "M"
    return default


def top_concentration_nm(row: dict, cols: dict) -> float | None:
    """The highest concentration this row was actually tested at.

    Read from the row's own conc columns, which are in molar. The alternative
    -- a single --max-conc flag -- writes one censoring point across libraries
    that were screened over different ranges, and every negative then carries
    a concentration it was never tested at. The real export tops out at 2.0E-5
    (20 uM) on the 4-point primary and lower on some plates.
    """
    best: float | None = None
    for c in cols.get("conc_cols", ()):
        v = _f(row.get(c))
        if v is None or v <= 0:
            continue
        nm = v * 1e9
        if best is None or nm > best:
            best = nm
    return best


def curve_class(row: dict, cols: dict) -> tuple[str, int]:
    """(magnitude, sign) of CURVE_CLASS2, e.g. '-2.4' -> ('2.4', -1).

    The sign is returned but NOT used to decide direction. NCATS's convention
    for it is not stated in the export, and efficacy carries the same
    information as a measured quantity -- so the sign is reported for
    provenance and efficacy decides. One inferred polarity convention per
    project is already one too many.
    """
    raw = (row.get(cols["curve"]) or "").strip()
    if not raw:
        return "", 0
    sign = -1 if raw.startswith("-") else 1
    return raw.lstrip("+-"), sign


def classify_row(row: dict, cols: dict, *, ac50_units: str,
                 tox_row: dict | None = None, tox_cols: dict | None = None,
                 tox_ac50_units: str = "uM",
                 min_efficacy_pct: float = MIN_EFFICACY_PCT,
                 inactive_at_nm: float = DEFAULT_INACTIVE_AT_NM,
                 ) -> tuple[str | None, dict, str]:
    """One screen row -> (verdict, qualifiers, reason).

    verdict is "active", "inactive" or None. The qualifiers are written so that
    labels.classify() reaches the SAME verdict from the emitted edge alone: an
    active carries relation "=" with its fitted EC50, an inactive carries
    relation ">" with the top concentration tested. That symmetry is the
    contract.
    """
    quals: dict = {"assay_type": "cell_based_antiviral", "cell_line": VERO_E6,
                   "virus_strain": "USA-WA1/2020"}
    cls, sign = curve_class(row, cols)
    eff = _f(row.get(cols["efficacy"]))
    ac50 = to_nm(row.get(cols["ac50"]), ac50_units) if "ac50" in cols else None
    top = top_concentration_nm(row, cols)
    if cls:
        quals["source_curve_class"] = f"{'-' if sign < 0 else ''}{cls}"
    if eff is not None:
        quals["source_efficacy_pct"] = eff

    # ---------------------------------------------------------------- inactive
    if cls in CURVE_CLASS_INACTIVE:
        if top is None:
            return None, quals, "inactive_without_a_concentration"
        if top < inactive_at_nm:
            return None, quals, "below_threshold"
        quals["ec50_nm"] = top
        quals["relation"] = ">"
        return "inactive", quals, "curve_class_4"

    # ------------------------------------------------------------ not a fit
    if cls not in CURVE_CLASS_FITTED:
        # Class 3 is single-point activity, class 5 a poor fit, blank is
        # unscored. 586 class-3 rows carry an AC50; admitting them on the
        # strength of that value alone would multiply the active set with
        # unreplicated single points.
        return None, quals, f"curve_class_{cls or 'blank'}_not_a_fit"

    # ------------------------------------------- a real fit: which direction?
    if eff is None:
        return None, quals, "fitted_without_efficacy"
    if eff < 0:
        # THE TRAP. Gain-of-signal assay: a negative efficacy means the
        # compound REDUCED viability. It has a clean AC50 and it is cytotoxic,
        # not antiviral. Counted, never emitted.
        return None, quals, "negative_efficacy_cytotoxic"
    if eff < min_efficacy_pct:
        return None, quals, "below_min_efficacy"
    if ac50 is None or ac50 <= 0:
        # Fit and direction agree but no usable potency. `unquantified` is
        # DERIVED at assembly (see chembl.py), so emit without a value rather
        # than invent one.
        quals["relation"] = "="
        return "active", quals, "fitted_unquantified"

    quals["ec50_nm"] = ac50
    quals["relation"] = "="

    # ------------------------------------------------- same-plate selectivity
    if tox_row is not None and tox_cols is not None:
        cc50 = to_nm(tox_row.get(tox_cols["ac50"]), tox_ac50_units) \
            if "ac50" in tox_cols else None
        tcls, _ = curve_class(tox_row, tox_cols)
        if cc50 and cc50 > 0 and tcls in CURVE_CLASS_FITTED:
            quals["cc50_nm"] = cc50
            quals["selectivity_index"] = cc50 / ac50
            quals["selectivity_verified"] = True
            quals["selectivity_n_studies"] = 1
        elif tcls in CURVE_CLASS_INACTIVE:
            # Not cytotoxic up to the top concentration. A real and useful
            # result, but CC50 is censored, so no finite index is computed --
            # writing top_conc/ec50 would understate a compound that is
            # cleaner than the assay can measure.
            quals["cytotoxicity_not_detected"] = True
    return "active", quals, "fitted_active"


def ingest_cpe(em, rows: list[dict], cols: dict, normalizer, *,
               tox_by_sample: dict[str, dict] | None = None,
               tox_cols: dict | None = None,
               existing_compounds: frozenset[str] | set[str] = frozenset(),
               taxon: str = SARS_COV_2, source: str = SOURCE,
               date: str = DEFAULT_DATE, ac50_units: str = "uM",
               tox_ac50_units: str = "uM",
               assay_filter: str | None = CANONICAL_ASSAY,
               min_efficacy_pct: float = MIN_EFFICACY_PCT,
               inactive_at_nm: float = DEFAULT_INACTIVE_AT_NM) -> Counter:
    """Emit one antiviral-activity edge per usable row.

    `normalizer` is injected so tests need no rdkit and so the function cannot
    quietly fall back to the file's own identifiers: it takes a structure
    string and returns an object with .curie and .inchikey, or raises.

    `existing_compounds` is why this layer does not become the authority on
    `is_approved`. That property is required on SmallMolecule and boolean, with
    no unknown, and this file does not say whether a compound is approved. A
    layer that asserts it anyway is the H2 bug again -- a true statement about
    one layer that becomes false after the merge, resolved by layer precedence
    rather than by evidence. So a compound the graph already holds gets NO node
    from here, only its label edge. New compounds get is_approved False, a
    floor rather than a claim: no metapath or filter reads the field.
    """
    stats: Counter = Counter()
    seen: dict[str, str] = {}
    tox_by_sample = tox_by_sample or {}

    for row in rows:
        stats["rows"] += 1
        if assay_filter is not None and "assay" in cols:
            if _key(row.get(cols["assay"], "")) != _key(assay_filter):
                stats["other_assay"] += 1
                continue
        stats["in_assay"] += 1

        struct = (row.get(cols["smiles"]) or "").strip()
        if not struct:
            stats["no_structure"] += 1
            continue
        try:
            cid = normalizer(struct)
        except Exception:
            stats["unnormalizable"] += 1
            continue

        sample = (row.get(cols["sample"]) or "").strip() if "sample" in cols else ""
        tox_row = tox_by_sample.get(sample) if sample else None
        if tox_row is not None:
            stats["paired_with_counterscreen"] += 1

        verdict, quals, reason = classify_row(
            row, cols, ac50_units=ac50_units, tox_row=tox_row,
            tox_cols=tox_cols, tox_ac50_units=tox_ac50_units,
            min_efficacy_pct=min_efficacy_pct, inactive_at_nm=inactive_at_nm)
        stats[f"reason_{reason}"] += 1
        if verdict is None:
            continue

        drug = cid.curie
        if drug in seen:
            # One structure, two plate records. Keeping both would double its
            # DWPC weight and let one compound vote twice in the AUC. An
            # active beats an inactive: a compound that worked on any plate is
            # not evidence of inactivity.
            if seen[drug] == "inactive" and verdict == "active":
                for i, e in enumerate(em.edges):
                    if e["subject"] == drug and e["predicate"] == ANTIVIRAL_PREDICATE:
                        em.edges[i] = dict(e, qualifiers=quals)
                        break
                seen[drug] = "active"
                stats["upgraded_inactive_to_active"] += 1
            else:
                stats["duplicate_structure"] += 1
            continue
        seen[drug] = verdict

        if drug in existing_compounds:
            stats["joined_existing_compound"] += 1
        else:
            em.node(drug, "SmallMolecule",
                    smiles=struct,
                    inchikey_skel=cid.inchikey[:14],
                    # See the docstring: a floor, not a claim.
                    is_approved=False,
                    salt_collapsed=False, stereo_collapsed=False)
            stats["new_compound_node"] += 1
        em.edge(drug, ANTIVIRAL_PREDICATE, taxon,
                source=source, date=date, tier=1, quals=quals,
                pmids=[CPE_PUBLICATION])
        stats[verdict] += 1
        if quals.get("selectivity_verified"):
            stats["with_same_plate_si"] += 1
        if quals.get("cytotoxicity_not_detected"):
            stats["no_cytotoxicity_detected"] += 1

    return stats


def index_counterscreen(rows: list[dict], cols: dict,
                        assay_filter: str | None = COUNTERSCREEN_ASSAY) -> dict:
    """sample_id -> counterscreen row, for the same-plate CC50."""
    out: dict[str, dict] = {}
    for r in rows:
        if assay_filter is not None and "assay" in cols:
            if _key(r.get(cols["assay"], "")) != _key(assay_filter):
                continue
        sid = (r.get(cols["sample"]) or "").strip() if "sample" in cols else ""
        if sid:
            out.setdefault(sid, r)
    return out


def join_report(labelled: set[str], existing: set[str]) -> dict:
    """How much of this layer lands on compounds the graph already holds.

    `labelled` is every compound this layer wrote an edge for -- taken from the
    edges, not from em.nodes, because a compound that joined deliberately gets
    no node from this layer (see ingest_cpe).

    The number that matters is not how many nodes were created but how many
    joined, because a node that joins nothing contributes no path. Reported so
    a bad normalization shows up as a count rather than as a null result three
    scripts later -- cross-entity count comparison is this project's primary
    detector and it only works if the count is printed.
    """
    joined = labelled & existing
    return {"emitted": len(labelled), "already_in_graph": len(joined),
            "new_compounds": len(labelled - existing),
            "join_rate": (len(joined) / len(labelled)) if labelled else 0.0}
