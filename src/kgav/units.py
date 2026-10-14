"""Concentration conversion, in one place.

WHY THIS FILE EXISTS. `to_nm` was defined three times -- chembl.py,
host_targets.py and ncats.py -- and the three did not agree:

    unit      chembl        host_targets   ncats
    nM        yes           yes            yes
    uM        yes           yes            yes
    pM        yes           yes            RAISES
    mM        yes           yes            RAISES
    M         None          None           yes
    logM      None          None           yes
    unknown   None          None           RAISES

Three tables, three coverages, two opposite failure modes. This is the
calculation every potency threshold in the project runs through -- the 10 uM
inactive cut, active_below_nm, the censored-relation rule, every AUC computed
against the labels those produce -- so three versions of it is the worst
placed instance of a pattern this project has now found ten times.

THE TWO THAT LOOKED IDENTICAL WERE NOT. chembl divided by 1000 for pM;
host_targets multiplied by 0.001. Those differ in the last bit for 13.75% of
values (measured over 200,000 random magnitudes), e.g. 511274.72136909724 ->
511.2747213690972 against 511.2747213690973. Nothing downstream can feel one
ULP at pM scale, and that is not the point: two copies believed to be the
same were never checked, and the units they support had already diverged.

TWO CONTRACTS, BOTH DELIBERATE. They are not an accident of history and
collapsing them would be wrong:

  to_nm         the unit comes from a source's own column -- ChEMBL
                standard_units is an open set including things no conversion
                can handle (ug.mL-1 needs a molecular weight, and guessing
                one produces a number that looks measured). Unrecognised ->
                None, and the CALLER counts it.

  to_nm_strict  the unit was derived by this codebase, from a column name
                ncats.unit_for recognises. An unexpected value there means
                the inference broke, which is a bug and must not be silently
                dropped. Unrecognised -> ValueError.
"""
from __future__ import annotations

#: Multiply by this to reach nanomolar. ug.mL-1 is deliberately absent: it
#: needs the compound's molecular weight, so there is no factor to put here.
CONCENTRATION_FACTORS: dict[str, float] = {
    "pM": 0.001,
    "nM": 1.0,
    "uM": 1_000.0,
    "mM": 1_000_000.0,
    "M": 1_000_000_000.0,
}


def _f(value) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def to_nm(value, unit: str | None, *,
          allow: frozenset[str] | None = None) -> float | None:
    """Concentration to nanomolar, or None when it cannot be converted.

    `allow` restricts the accepted units. It exists so that adding a unit to
    CONCENTRATION_FACTORS cannot silently change what an existing ingest
    accepts: a source that has never been checked for 'M' rows keeps
    returning None for them until someone measures whether any exist.
    to_nm_strict takes the same argument for the same reason -- pinning the
    set is what makes the shared table safe to extend.
    """
    v = _f(value)
    if v is None:
        return None
    u = (unit or "").strip()
    if allow is not None and u not in allow:
        return None
    f = CONCENTRATION_FACTORS.get(u)
    return v * f if f is not None else None


def to_nm_strict(value, unit: str, *,
                 allow: frozenset[str] | None = None) -> float | None:
    """As to_nm, but an unrecognised unit is an error rather than a miss.

    For a unit this codebase derived rather than read from a source: there
    the unexpected value is a broken inference, and returning None would
    drop the row and the evidence of the bug together.

    logM is NCATS's log_ac50: log10 of molar. -5.0 is 10 uM, i.e. 10,000 nM.
    """
    v = _f(value)
    if v is None:
        return None
    if allow is not None and unit not in allow:
        raise ValueError(f"unknown concentration unit {unit!r}")
    if unit == "logM":
        return (10.0 ** v) * CONCENTRATION_FACTORS["M"]
    f = CONCENTRATION_FACTORS.get(unit)
    if f is None:
        raise ValueError(f"unknown concentration unit {unit!r}")
    return v * f
