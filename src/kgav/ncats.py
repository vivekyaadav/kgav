"""NCATS OpenData SARS-CoV-2 CPE screen: labels in the population that matters.

    SmallMolecule -HAS_ANTIVIRAL_ACTIVITY_AGAINST-> OrganismTaxon

WHY THIS LAYER EXISTS. Four evaluation protocols have now scored at chance or
below once their confounds were controlled (doc2/08-results.md). Two of the six
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
from dataclasses import dataclass
from pathlib import Path

# THE THRESHOLD HAS ONE OWNER. Importing it rather than restating
# it is the point: a second 10_000 here would be the "check written
# with the literals it is checking" shape, and the two would drift.
from kgav.labels import DEFAULT_ACTIVE_BELOW_NM
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


# --------------------------------------------------------------- assay specs
#: GAIN vs LOSS of signal, and why this is a declared field rather than an
#: inference. The CPE assay reads cell viability: a compound that protects
#: cells RAISES the signal, so antiviral activity is POSITIVE efficacy. The
#: pseudotyped-particle entry assays read a reporter carried INTO the cell:
#: a compound that blocks entry LOWERS the signal, so the same activity is
#: NEGATIVE efficacy and a negative curve class.
#:
#: Measured on the ODP exports before this was written: applying the CPE
#: rule to the three entry files yields ZERO actives, because all 1,004
#: fitted inhibitors sit at negative efficacy. A layer with one global
#: convention emits three all-negative cells that validate, count plausibly
#: and encode the opposite conclusion -- the failure this module's own
#: ingest script warns about in those words.
#:
#: So polarity is DECLARED per assay and CHECKED against the file
#: (`observed_polarity`), never inferred at runtime: a truncated or empty
#: download would otherwise redefine what activity means.
#: Distinguishes "caller said nothing" from "caller said None",
#: which means ALL assays. A plain None default would have made
#: the two indistinguishable the moment assay_filter started
#: defaulting from the spec.
_UNSET = "\x00unset"

GAIN_OF_SIGNAL = 1
LOSS_OF_SIGNAL = -1

#: A pseudotyped particle carries the named virus's spike on another viral
#: core, so inhibiting its entry is evidence about ENTRY and not about
#: replication of the authentic virus. It is weaker evidence than a
#: cytopathic-effect assay and carries its own assay_type so a consumer can
#: separate the two. Registered in doc2/registrations/prereg_ncats_entry_assays.md.
CELL_BASED_ANTIVIRAL = "cell_based_antiviral"
CELL_BASED_PSEUDOTYPE_ENTRY = "cell_based_pseudotype_entry"

#: THE ENTRY ASSAYS ARE NOT THE CPE PAPER, and the first wiring of them
#: cited PMID:33708112 for all three -- including a MERS screen. Verified
#: against PubMed rather than inferred:
#:
#:   PMID:35086793  Xu M et al., SLAS Discov 2022 Mar. "A high throughput
#:                  screening assay for inhibitors of SARS-CoV-2
#:                  pseudotyped particle entry."
#:   PMID:33330839  Chen CZ et al., ACS Pharmacol Transl Sci 2020-12-11.
#:                  "Identifying SARS-CoV-2 Entry Inhibitors through Drug
#:                  Repurposing Screens of SARS-S and MERS-S Pseudotyped
#:                  Particles." ONE paper covers both of those assays.
PP_SARS_COV_2_PUBLICATION = "PMID:35086793"
PP_SARS_S_MERS_S_PUBLICATION = "PMID:33330839"

#: `first_asserted_date` is the cited paper's publication date, which is
#: VERIFIABLE. The CPE layer uses a bioRxiv preprint date because one is
#: known for it; no preprint date is known for these, and the portal's file
#: mtime is a re-export date -- the thing DEFAULT_DATE's comment refuses.
#: So the journal date is used and the choice is stated rather than the
#: earlier-but-unknown date being guessed. It feeds the temporal split, so
#: it is deliberately the conservative (later) end.
PP_SARS_COV_2_DATE = "2022-03-01"
PP_SARS_S_MERS_S_DATE = "2020-12-11"

SARS_COV = "NCBITaxon:694009"
MERS_COV = "NCBITaxon:1335626"
#: Cellosaurus: HEK293T (CVCL_0063) and Huh-7 (CVCL_0336).
HEK293T = "CVCL:0063"
HUH7 = "CVCL:0336"


@dataclass(frozen=True)
class AssaySpec:
    """One NCATS assay, and everything about it that is NOT in its rows.

    The file states efficacy, curve class and AC50. It does NOT state which
    virus the assay is about, which direction means activity, which cell
    line, or which counterscreen pairs with it -- and every one of those
    changes the meaning of the numbers. They are declared here, beside each
    other, so a new assay cannot be added by passing a path.
    """

    key: str
    assay_name: str                  # as it appears in `assay_name`, lowercased
    taxon: str
    polarity: int                    # GAIN_OF_SIGNAL or LOSS_OF_SIGNAL
    assay_type: str
    cell_line: str
    counterscreen_assay: str
    #: Local basenames under data/raw/ncats, and the portal filenames they
    #: came from. THE REMOTE NAME LIVES HERE because cpe_assay.csv's
    #: provenance record says `"url": null` with a note that it "is not
    #: recoverable" -- while this repository's own ingest script documents
    #: the URL in its header. A download route that is known but not
    #: machine-readable is the same shape as one that is lost, so every
    #: assay added from 2026-11-17 carries its own.
    file: str = "cpe_assay.csv"
    counterscreen_file: str = "cpe_tox.csv"
    remote: str = "SARS-CoV-2_Cytopathic_Effect_(CPE).csv"
    remote_counterscreen: str = (
        "SARS-CoV-2_Cytopathic_Effect_(Host_Tox_Counterscreen).csv")
    publication: str = CPE_PUBLICATION
    date: str = DEFAULT_DATE
    virus_strain: str | None = None


CPE = AssaySpec(
    key="cpe",
    assay_name=CANONICAL_ASSAY,
    taxon=SARS_COV_2,
    polarity=GAIN_OF_SIGNAL,
    assay_type=CELL_BASED_ANTIVIRAL,
    cell_line=VERO_E6,
    counterscreen_assay=COUNTERSCREEN_ASSAY,
    virus_strain="USA-WA1/2020",
)

#: The three entry assays, registered 2026-11-17 before this code existed.
#: The enzymatic and binding assays on the same portal (3CL, RdRp, TMPRSS2,
#: ACE2, Spike-ACE2) are REFUSED: they measure target engagement, and a
#: compound that inhibits a protease in a tube has not been shown to inhibit
#: a virus. They are absent from this table on purpose.
PP_SARS_COV_2 = AssaySpec(
    key="pp_sars_cov_2",
    publication=PP_SARS_COV_2_PUBLICATION,
    date=PP_SARS_COV_2_DATE,
    file="pp_sars_cov_2.csv",
    counterscreen_file="pp_sars_cov_2_tox.csv",
    remote="SARS-CoV-2_Pseudotyped_Particle_Entry.csv",
    remote_counterscreen="SARS-CoV-2_Pseudotyped_Particle_Entry_(Tox_Counterscreen).csv",
    assay_name="sars-cov-2 pseudotyped particle entry",
    taxon=SARS_COV_2,
    polarity=LOSS_OF_SIGNAL,
    assay_type=CELL_BASED_PSEUDOTYPE_ENTRY,
    cell_line=HEK293T,
    counterscreen_assay="sars-cov-2 pseudotyped particle entry (tox counterscreen)",
)
PP_SARS_COV = AssaySpec(
    key="pp_sars_cov",
    publication=PP_SARS_S_MERS_S_PUBLICATION,
    date=PP_SARS_S_MERS_S_DATE,
    file="pp_sars_cov.csv",
    counterscreen_file="pp_sars_cov_tox.csv",
    remote="SARS-CoV_Pseudotyped_Particle_Entry.csv",
    remote_counterscreen="SARS-CoV_Pseudotyped_Particle_Entry_(Vero_E6_Tox_Counterscreen).csv",
    assay_name="sars-cov pseudotyped particle entry",
    taxon=SARS_COV,
    polarity=LOSS_OF_SIGNAL,
    assay_type=CELL_BASED_PSEUDOTYPE_ENTRY,
    cell_line=VERO_E6,
    counterscreen_assay="sars-cov pseudotyped particle entry (veroe6 tox counterscreen)",
)
PP_MERS = AssaySpec(
    key="pp_mers",
    publication=PP_SARS_S_MERS_S_PUBLICATION,
    date=PP_SARS_S_MERS_S_DATE,
    file="pp_mers.csv",
    counterscreen_file="pp_mers_tox.csv",
    remote="MERS_Pseudotyped_Particle_Entry.csv",
    remote_counterscreen="MERS_Pseudotyped_Particle_Entry_(Huh7_Tox_Counterscreen).csv",
    assay_name="mers pseudotyped particle entry",
    taxon=MERS_COV,
    polarity=LOSS_OF_SIGNAL,
    assay_type=CELL_BASED_PSEUDOTYPE_ENTRY,
    cell_line=HUH7,
    counterscreen_assay="mers pseudotyped particle entry (huh7 tox counterscreen)",
)

ASSAYS: dict[str, AssaySpec] = {
    s.key: s for s in (CPE, PP_SARS_COV_2, PP_SARS_COV, PP_MERS)}

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


class PolarityError(ValueError):
    """A file's own rows disagree with the polarity declared for its assay.

    Deliberately fatal, like ColumnError and for the same reason: a wrong
    direction produces a release that validates, counts plausibly and
    encodes the opposite conclusion. The remedy is to look at the file, not
    to let the data pick a convention at runtime -- a truncated download
    would then silently redefine activity.
    """


def observed_polarity(rows: list[dict], cols: dict,
                      min_efficacy_pct: float = MIN_EFFICACY_PCT,
                      ) -> tuple[int, int, int]:
    """(sign, n_negative, n_positive) over this file's FITTED rows.

    The sign is whichever direction the strong responders actually go, read
    only from rows with a real concentration-response fit and an efficacy
    at or beyond the threshold -- class 4 and the inconclusive classes say
    nothing about direction and are excluded.

    This does NOT decide the convention. It checks the declared one, so a
    file that arrives re-exported with its sign flipped stops the build
    rather than inverting every label in the layer.
    """
    neg = pos = 0
    for row in rows:
        cls, _sign = curve_class(row, cols)
        if cls not in CURVE_CLASS_FITTED:
            continue
        eff = _f(row.get(cols["efficacy"]))
        if eff is None or abs(eff) < min_efficacy_pct:
            continue
        if eff < 0:
            neg += 1
        else:
            pos += 1
    if neg == pos:
        return 0, neg, pos
    return (LOSS_OF_SIGNAL if neg > pos else GAIN_OF_SIGNAL), neg, pos


def check_polarity(rows: list[dict], cols: dict, spec: AssaySpec,
                   min_efficacy_pct: float = MIN_EFFICACY_PCT) -> tuple[int, int]:
    """Refuse a file whose strong responders contradict its declared sign."""
    got, neg, pos = observed_polarity(rows, cols, min_efficacy_pct)
    if got == 0:
        raise PolarityError(
            f"{spec.key}: {neg} negative and {pos} positive strong responders, "
            f"so the file states no direction. Declared "
            f"{'loss' if spec.polarity < 0 else 'gain'}-of-signal; look at the "
            f"file rather than letting it choose")
    if got != spec.polarity:
        raise PolarityError(
            f"{spec.key} is declared "
            f"{'loss' if spec.polarity < 0 else 'gain'}-of-signal and its own "
            f"rows say the opposite: {neg:,} fitted responders at or beyond "
            f"-{min_efficacy_pct:g}% efficacy against {pos:,} at or beyond "
            f"+{min_efficacy_pct:g}%. One of the two is wrong and guessing "
            f"here inverts every label in this assay")
    return neg, pos


def classify_row(row: dict, cols: dict, *, ac50_units: str,
                 conc_units: str | None = None,
                 tox_row: dict | None = None, tox_cols: dict | None = None,
                 tox_ac50_units: str = "uM",
                 min_efficacy_pct: float = MIN_EFFICACY_PCT,
                 inactive_at_nm: float = DEFAULT_INACTIVE_AT_NM,
                 active_below_nm: float = DEFAULT_ACTIVE_BELOW_NM,
                 spec: AssaySpec = CPE,
                 ) -> tuple[str | None, dict, str]:
    """One screen row -> (verdict, qualifiers, reason).

    verdict is "active", "inactive" or None. The qualifiers are written so that
    labels.classify() reaches the SAME verdict from the emitted edge alone: an
    active carries relation "=" with its fitted EC50, an inactive carries
    relation ">" with the top concentration tested. That symmetry is the
    contract.
    """
    quals: dict = {"assay_type": spec.assay_type, "cell_line": spec.cell_line}
    if spec.virus_strain:
        quals["virus_strain"] = spec.virus_strain
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
    # THE DIRECTION IS THE ASSAY'S AND NOT THIS MODULE'S. `eff` is multiplied
    # by the declared polarity, so "active" always means beyond the threshold
    # in the direction that assay calls inhibition: RAISING viability for
    # CPE, LOWERING the entry reporter for a pseudotype. A single global sign
    # test read all 1,004 fitted entry inhibitors as no effect, which is how
    # this field came to exist.
    signal = eff * spec.polarity
    if signal < 0:
        # THE TRAP, and it stays one under either polarity. For CPE a
        # negative efficacy means the compound REDUCED viability -- a clean
        # AC50, cytotoxic rather than antiviral. For an entry assay the
        # mirror case is a compound that INCREASED entry. Counted, never
        # emitted, and the reason keeps its original name so a run of either
        # assay stays comparable with every run before this change.
        return None, quals, "negative_efficacy_cytotoxic"
    if signal < min_efficacy_pct:
        return None, quals, "below_min_efficacy"
    if ac50 is None or ac50 <= 0:
        # Fit and direction agree but no usable potency. `unquantified` is
        # DERIVED at assembly (see chembl.py), so emit without a value rather
        # than invent one.
        quals["relation"] = "="
        return "active", quals, "fitted_unquantified"

    # ------------------------------- a fitted hit the threshold calls neither
    # A COMPOUND THE ASSAY FITTED AS AN INHIBITOR IS NOT EVIDENCE OF
    # INACTIVITY, and until 2026-11-17 it was counted as exactly that. This
    # function returned "active" for any fitted hit, and `labels.classify`
    # then read its own threshold off the emitted edge and returned
    # "inactive" -- so 17 CPE compounds with EC50 between 11,220 and 25,119
    # nM were negatives in every cell this layer feeds, which is a KNOWN
    # weak inhibitor sitting in the class of non-inhibitors.
    #
    # It also broke this function's stated contract, two paragraphs up: the
    # qualifiers are written so that labels.classify reaches the SAME
    # verdict from the edge alone. It did not, and nothing checked it;
    # tests/test_ncats.py now does, over every row of every assay.
    #
    # The refusal is NEITHER class, counted and named -- not "active"
    # (above the threshold the project uses everywhere) and not "inactive"
    # (the assay fitted a curve). Decided 2026-11-17 and recorded in
    # doc2/registrations/prereg_ncats_entry_assays.md as the amendment that followed P1
    # failing.
    #
    # SCOPE, STATED: this is the NCATS layer only. The same question exists
    # for ChEMBL -- a published "= 25 uM" becomes an inactive label there
    # too -- over a much larger population, and it is NOT changed here and
    # NOT measured. Filed rather than folded in.
    if ac50 > active_below_nm:
        quals["ec50_nm"] = ac50
        quals["relation"] = "="
        return None, quals, "weak_hit_above_active_threshold"

    quals["ec50_nm"] = ac50
    quals["relation"] = "="

    # ------------------------------------------------- same-plate selectivity
    if tox_row is not None and tox_cols is not None:
        cc50 = to_nm(tox_row.get(tox_cols["ac50"]), tox_ac50_units) \
            if "ac50" in tox_cols else None
        tcls, _ = curve_class(tox_row, tox_cols)
        if cc50 and cc50 > 0 and tcls in CURVE_CLASS_FITTED:
            # CYTOTOXICITY IS A CONFOUND IN ONE DIRECTION ONLY, and that is
            # why this refusal is tied to polarity rather than switched on
            # per assay. In a LOSS-of-signal assay the readout is a reporter
            # carried into the cell, so a compound that kills cells lowers
            # it by exactly the route a compound that blocks entry does:
            # cytotoxicity MIMICS the activity and is indistinguishable from
            # it in the fitted curve. In the GAIN-of-signal CPE assay the
            # two have opposite signs, so a cytotoxic compound is already
            # refused above as `negative_efficacy_cytotoxic` and a second
            # test here would be redundant -- which is why the CPE cell
            # registered in prereg_structure_scorer.md does not move.
            #
            # Registered in doc2/registrations/prereg_ncats_entry_assays.md before this
            # code: "a compound whose counterscreen AC50 is at or below its
            # antiviral AC50 is killing cells at the concentration where it
            # works, and is REFUSED rather than called active."
            if spec.polarity == LOSS_OF_SIGNAL and cc50 <= ac50:
                quals["cc50_nm"] = cc50
                return None, quals, "cytotoxic_at_or_below_ec50"
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
               taxon: str | None = None, source: str = SOURCE,
               date: str | None = None, ac50_units: str = "uM",
               tox_ac50_units: str = "uM",
               conc_units: str | None = None,
               assay_filter: str | None = _UNSET,
               min_efficacy_pct: float = MIN_EFFICACY_PCT,
               inactive_at_nm: float = DEFAULT_INACTIVE_AT_NM,
               active_below_nm: float = DEFAULT_ACTIVE_BELOW_NM,
               spec: AssaySpec = CPE) -> Counter:
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
    # EVERY ASSAY-DEPENDENT VALUE COMES FROM ONE PLACE. taxon, the assay
    # name to keep and the citation used to be three arguments a caller had
    # to keep consistent with each other; three assays later that is three
    # chances to file a MERS screen against SARS-CoV-2. They now default
    # from `spec`, and an explicit argument still wins so existing callers
    # and tests are unaffected.
    taxon = spec.taxon if taxon is None else taxon
    date = spec.date if date is None else date
    if assay_filter is _UNSET:
        assay_filter = spec.assay_name
    publication = spec.publication
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
            min_efficacy_pct=min_efficacy_pct, inactive_at_nm=inactive_at_nm,
            active_below_nm=active_below_nm, spec=spec)
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
                pmids=[publication])
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
