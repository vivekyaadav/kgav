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
import re
from collections import Counter
from pathlib import Path

from kgav.normalize_chem import tautomer_conflations
from kgav.units import to_nm_strict

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


#: NCGC sample ids are `NCGC<substance>-<lot>`, and the lot is a BATCH of
#: the same substance, not a different compound. Measured on cpe_assay.csv:
#: 8,596 distinct substance ids, of which 815 appear under more than one lot.
LOT_SUFFIX = re.compile(r"^(NCGC\d+)-\d+$")


def substance_id(sample: str) -> str:
    """The substance a sample id names, with its lot suffix removed.

    WHY fix-list L19 NEEDS THIS. A conflation is two distinct SUBSTANCES
    pooled onto one node. Two lots of one substance are one compound, so
    counting `NCGC00263533-01` and `NCGC00263533-12` as two members
    overstates it -- and that is not hypothetical: penciclovir appears under
    both, with inconsistent SMILES between them, which is a spelling defect
    in one row rather than two compounds sharing a node.

    A sample id that does not match the pattern is returned unchanged rather
    than guessed at; the count then treats it as its own substance, which is
    the conservative direction for a detector.
    """
    if not sample:
        return "?"
    m = LOT_SUFFIX.match(sample.strip())
    return m.group(1) if m else sample.strip()


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


#: Exactly what unit_for can produce. Pinned so that a unit added to
#: CONCENTRATION_FACTORS cannot quietly become acceptable here: outside
#: this set the inference is broken, which is the error this raises.
NCATS_UNITS = frozenset({"M", "uM", "nM", "logM"})


def to_nm(value, unit: str) -> float | None:
    """Concentration to nanomolar. `unit` is one of M, uM, nM, logM.

    STRICT, and deliberately unlike chembl.to_nm: `unit` here is derived by
    unit_for from a column name this codebase recognises, so an unexpected
    value means the inference broke. Returning None would drop the row and
    the evidence of the bug together. See kgav.units for both contracts.
    """
    return to_nm_strict(value, unit, allow=NCATS_UNITS)


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


#: Bounds on the TOP of a dose series, which is what this function returns.
#:
#: MEASURED, not guessed. Over the 10,707 rows of cpe_assay.csv that carry
#: any concentration, the per-row maximum runs 50 nM to 200,000 nM with a
#: median of 20,000 nM (20 uM, the 4-point primary); 8,424 rows sit exactly
#: at 20,000 and only 48 fall below 1,000, the lowest being 50. So 10 nM
#: to 1 M clears every real row by at least 5x.
#:
#: It still catches both unit errors a bare `conc0` column invites, because
#: a factor of a million does not land near a bound: micromolar read as
#: molar sends 20 to 2e10 nM, above the ceiling, and molar read as
#: micromolar sends 2e-5 to 0.02 nM, below the floor.
#:
#: TWO EARLIER VERSIONS OF THIS BOUND WERE WRONG, both rejected by the real
#: panel on its first row, and both for the same reason this fix exists --
#: a plausible-sounding range written without looking at the file:
#:
#:   every conc column, 1 nM - 1 M   conc0 is 4e-10 M = 0.4 nM, the BOTTOM
#:                                   of a series running to 40 uM. A
#:                                   dose-response series is meant to start
#:                                   below the active range; only the top
#:                                   has a plausible range, and only the
#:                                   top leaves this function.
#:   the top, 100 nM - 1 M           48 rows top out below 1,000 nM and one
#:                                   at 50 nM.
PLAUSIBLE_TOP_NM = (10.0, 1e9)


def top_concentration_nm(row: dict, cols: dict,
                         conc_units: str | None = None) -> float | None:
    """The highest concentration this row was actually tested at.

    Read from the row's own conc columns. The alternative -- a single
    --max-conc flag -- writes one censoring point across libraries that
    were screened over different ranges, and every negative then carries a
    concentration it was never tested at. The real export tops out at
    2.0E-5 (20 uM) on the 4-point primary and lower on some plates.

    THE UNIT IS READ, NOT ASSUMED (fix-list P7). This was `nm = v * 1e9`,
    hardcoding molar on the strength of a docstring.

    Three of the portal's files, three answers. cpe.tsv names the unit --
    CONC0(M) through CONC14(M) -- and unit_for() in this module already
    parses exactly that. The two files this layer ACTUALLY reads,
    cpe_assay.csv and cpe_tox.csv, use a bare `conc0`, so the header
    cannot answer it there and the assumption really was load bearing.

    WHAT THE UNIT IS, AND HOW THAT IS KNOWN. Measured across the files:
    cpe_assay.csv's conc3 runs 4e-08 to 2e-04 with median 2e-05, and
    cpe.tsv's CONC3(M) runs 4e-08 to 6e-05 with the same median. Identical
    magnitudes, and one of the two states molar -- so the unlabelled
    columns are molar by evidence rather than by assertion. The driver
    passes that as `--conc-units M` so the claim is declared at the call
    site next to the measurement supporting it, instead of being a
    multiplier nobody can question.

    IT WAS ALSO INCONSISTENT WITH unit_for. That function defaults to "uM"
    for a column whose name carries no unit, while this multiplier assumed
    M: the same unlabelled column would have been read two ways by two
    functions in one file, differing by a factor of a million. A conc
    column with no stated unit is now REFUSED rather than defaulted,
    because this value is the censoring point for every inactive label the
    layer emits -- `inactive_at_nm` compares against it -- so a silent
    factor of 1e6 would move every negative in the panel.
    """
    best: float | None = None
    best_col = best_raw = best_unit = None
    for c in cols.get("conc_cols", ()):
        v = _f(row.get(c))
        if v is None or v <= 0:
            continue
        # THE HEADER FIRST, THE CALLER SECOND, NEVER A SILENT DEFAULT.
        unit = unit_for(c, default="") or (conc_units or "")
        if not unit:
            raise ValueError(
                f"concentration column {c!r} states no unit in its name and "
                f"no conc_units was given. unit_for would default it to uM "
                f"and the previous hardcoded multiplier assumed M -- a "
                f"factor of a million on the censoring point of every "
                f"inactive label. Either name the unit in the header, as "
                f"cpe.tsv's CONC0(M) does, or pass conc_units with evidence "
                f"for it.")
        nm = to_nm(v, unit)
        if nm is None:
            continue
        if best is None or nm > best:
            best = nm
            best_col, best_raw, best_unit = c, v, unit
    # THE TOP, ONCE, not every point of the series.
    if best is not None:
        lo, hi = PLAUSIBLE_TOP_NM
        if not lo <= best <= hi:
            raise ValueError(
                f"top concentration {best_col} = {best_raw} {best_unit} is "
                f"{best:.3g} nM, outside the plausible range "
                f"{lo:.0f}-{hi:.0g} nM for the top of a dose series. A unit "
                f"error lands here rather than slightly off; check the "
                f"column header and --conc-units before relaxing this.")
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
                 conc_units: str | None = None,
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
    top = top_concentration_nm(row, cols, conc_units)
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
               conc_units: str | None = None,
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
    # fix-list L19. node id -> raw SMILES -> the sample ids that supplied it.
    # The id is derived with the canonical tautomer, so two library members
    # that differ only in where a proton sits land on ONE node and pool their
    # screening results. Nothing recorded that, because a layer writes one
    # node per id and the collapse is only visible from the raw rows.
    derivations: dict[str, dict[str, set[str]]] = {}

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
        sample_col = cols.get("sample")
        sample = (row.get(sample_col) or "").strip() if sample_col else ""
        derivations.setdefault(f"INCHIKEY:{cid.inchikey}", {}) \
            .setdefault(struct, set()).add(substance_id(sample))

        sample = (row.get(cols["sample"]) or "").strip() if "sample" in cols else ""
        tox_row = tox_by_sample.get(sample) if sample else None
        if tox_row is not None:
            stats["paired_with_counterscreen"] += 1

        verdict, quals, reason = classify_row(
            row, cols, ac50_units=ac50_units, conc_units=conc_units,
            tox_row=tox_row,
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
                # THE DISAGREEMENT SURVIVES THE OVERWRITE (fix-list L14).
                # This replaced the inactive edge's qualifiers in place and
                # the only trace was a counter, so an edge that the panel
                # itself contradicted was indistinguishable downstream from
                # one measured active on every plate. labels.py's stated
                # policy is that an ambiguous compound is EXCLUDED rather
                # than assigned, and build_labels implements it -- but it
                # cannot see a disagreement the ingest already settled,
                # and identity_includes_source cannot separate two records
                # from one source.
                #
                # Marking it is what makes that decision reachable. Whether
                # build_labels should then drop these compounds is a
                # protocol question, so it is NOT changed here: the label
                # set is untouched and the fact is now on the edge.
                for i, e in enumerate(em.edges):
                    if e["subject"] == drug and e["predicate"] == ANTIVIRAL_PREDICATE:
                        em.edges[i] = dict(
                            e, qualifiers={**quals,
                                           "replicate_disagreement": True})
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
            # fix-list L18. `smiles` DESCRIBES THE NODE, not the row.
            # ncats supplies salt and deprotonated forms for compounds whose
            # id is the free base's InChIKey -- `...fumarate`, `Cl.X`, a
            # carboxylate `C([O-])=O` -- so storing the raw string put a
            # structure that resolves to a DIFFERENT key on the node. 108 of
            # them within this layer alone, measured 2026-11-11; the row had
            # it filed at 4, from a cross-layer diagnostic on a build that
            # was declined. cid.smiles is the de-salted, uncharged parent
            # that the id was derived from, so the two now agree.
            #
            # THE SOURCE'S OWN STRING IS KEPT, not overwritten -- it is
            # ncats's statement about its sample and losing it would be the
            # opposite error -- but only where normalisation changed the
            # SUBSTANCE, not merely its spelling.
            #
            # `struct != cid.smiles` WAS THE WRONG TEST AND THE REBUILD SAID
            # SO: cid.smiles is RDKit's canonical SMILES, which differs from
            # the source's spelling for almost every compound, so the
            # property landed on 5,134 of 5,234 ncats nodes. A field present
            # on 98% of a layer carries almost no information, and it buried
            # the ~108 cases this row is about -- the salt and deprotonated
            # forms -- among 5,026 cosmetic reorderings. The flags below are
            # the ones that mean "a different form of the compound was
            # supplied": a counterion was stripped or a charge neutralised.
            #
            # tautomer_changed IS NOT ONE OF THEM, and the first attempt
            # included it. Measured on the rebuild: of 1,352 nodes it
            # selected, 1,130 (84%) were tautomer canonicalisation ALONE
            # against 124 salts and 98 charged forms -- so it buried the
            # cases this row is about for the second time. Under InChI's
            # mobile-H layer two tautomers are the SAME COMPOUND, which is
            # why these nodes exist at all (fix-list L19) and why the graph
            # now declares `tautomer_collapsed: true`. A canonicalised
            # tautomer is a different spelling of what the node holds, not
            # a different substance.
            differs = cid.had_multiple_fragments or cid.was_charged
            em.node(drug, "SmallMolecule",
                    smiles=cid.smiles,
                    smiles_as_supplied=(struct if differs else None),
                    inchikey_skel=cid.inchikey[:14],
                    # See the docstring: a floor, not a claim.
                    is_approved=False,
                                        # THE IDENTITY POLICY, DECLARED (fix-list L24). Not a
                    # per-compound measurement: normalize_chem's docstring
                    # says these record WHICH KEY IS THE NODE ID, so a
                    # constant is correct and both layers must state the
                    # same one. salt_collapsed was False on all 25,383
                    # nodes and that is wrong -- largest_fragment.choose
                    # runs before the InChIKey, so `Cl.X` and `X` really
                    # do share a node. stereo_collapsed stays False: the
                    # id is the FULL key and L- and D-alanine differ in it.
                    # tests/test_identity_policy.py derives both from what
                    # normalize_chemical actually does.
                    salt_collapsed=True, stereo_collapsed=False,
                    tautomer_collapsed=True)
            stats["new_compound_node"] += 1
        em.edge(drug, ANTIVIRAL_PREDICATE, taxon,
                source=source, date=date, tier=1, quals=quals,
                pmids=[CPE_PUBLICATION])
        stats[verdict] += 1
        if quals.get("selectivity_verified"):
            stats["with_same_plate_si"] += 1
        if quals.get("cytotoxicity_not_detected"):
            stats["no_cytotoxicity_detected"] += 1

    for members in tautomer_conflations(derivations).values():
        stats["tautomer_conflated_nodes"] += 1
        stats["tautomer_conflated_members"] += members

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
