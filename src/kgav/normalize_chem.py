"""Chemical identifier normalisation.

Turns any structure representation into a stable node key: an InChIKey, plus a
14-character skeleton that groups salts, stereoisomers and tautomers together.

Two keys, deliberately:
  inchikey       full key — distinguishes stereoisomers
  inchikey_skel  first block — the grouping key

Which one is the node id is a project decision recorded on every node as
salt_collapsed / stereo_collapsed. The graph uses the FULL key as id and carries
the skeleton as a property, so a bad merge is always reversible.
"""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

_UNCHARGER = rdMolStandardize.Uncharger()
_LFC = rdMolStandardize.LargestFragmentChooser()
_TE = rdMolStandardize.TautomerEnumerator()


@dataclass(frozen=True)
class ChemId:
    inchikey: str
    inchikey_skel: str
    smiles: str
    had_multiple_fragments: bool
    was_charged: bool
    tautomer_changed: bool

    @property
    def curie(self) -> str:
        return f"INCHIKEY:{self.inchikey}"


class NormalizationError(ValueError):
    pass


def _parse(structure: str) -> Chem.Mol:
    """Accept SMILES, InChI, or molblock. Never guess silently."""
    s = (structure or "").strip()
    if not s:
        raise NormalizationError("empty structure")

    if s.startswith("InChI="):
        mol = inchi.MolFromInchi(s, sanitize=True)
    elif "\n" in s and ("V2000" in s or "V3000" in s):
        mol = Chem.MolFromMolBlock(s, sanitize=True)
    else:
        mol = Chem.MolFromSmiles(s, sanitize=True)

    if mol is None:
        raise NormalizationError(f"could not parse structure: {s[:60]!r}")
    return mol


def normalize_chemical(structure: str, *, canonical_tautomer: bool = True) -> ChemId:
    """Standardise a structure and return its identifiers.

    Pipeline, in order:
      1. parse and sanitize
      2. cleanup (normalize functional groups, disconnect metals, reionize)
      3. largest fragment  — strips counterions
      4. uncharge          — neutralises the parent
      5. canonical tautomer (optional; slow on large molecules)
      6. InChIKey

    Order matters. Uncharging before fragment selection can neutralise a
    counterion into something that then wins the largest-fragment contest.
    """
    mol = _parse(structure)

    n_frags_before = len(Chem.GetMolFrags(mol))
    charge_before = Chem.GetFormalCharge(mol)

    mol = rdMolStandardize.Cleanup(mol)
    mol = _LFC.choose(mol)
    mol = _UNCHARGER.uncharge(mol)

    # InChIKey is computed from the standardised-but-NOT-tautomer-canonicalised
    # molecule. RDKit's TautomerEnumerator can strip stereocentres adjacent to a
    # tautomerisable group, which would silently merge enantiomers into one node.
    # For drugs that difference is frequently the whole pharmacology.
    key = inchi.MolToInchiKey(mol)
    if not key:
        raise NormalizationError("InChIKey generation returned empty")

    # The canonical tautomer is used only as a GROUPING signal, never as the id.
    tautomer_changed = False
    if canonical_tautomer:
        before = Chem.MolToSmiles(mol)
        try:
            taut = _TE.Canonicalize(mol)
        except Exception as exc:  # RDKit raises on pathological inputs
            raise NormalizationError(f"tautomer canonicalisation failed: {exc}") from exc
        tautomer_changed = Chem.MolToSmiles(taut) != before

    return ChemId(
        inchikey=key,
        inchikey_skel=key[:14],
        smiles=Chem.MolToSmiles(mol),
        had_multiple_fragments=n_frags_before > 1,
        was_charged=charge_before != 0,
        tautomer_changed=tautomer_changed,
    )


def same_as_rule(a: ChemId, b: ChemId) -> str | None:
    """Why two ChemIds would be merged, or None if they would not be.

    Returned string goes in the SAME_AS edge's merge_rule qualifier, so every
    collapse in the graph carries its own justification.
    """
    if a.inchikey == b.inchikey:
        return "identical_inchikey"
    if a.inchikey_skel == b.inchikey_skel:
        return "shared_inchikey_skeleton:salt_or_stereo_or_tautomer"
    return None
