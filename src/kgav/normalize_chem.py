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
from typing import TYPE_CHECKING

if TYPE_CHECKING:                     # annotations only; never imported at runtime
    from rdkit import Chem

# rdkit is an OPTIONAL extra ([chem]), so it is imported inside the functions
# that need it rather than at module scope. Imported at module scope, a missing
# rdkit makes this module unimportable -- which takes down every test that
# imports it, and a collection error aborts the WHOLE suite, not just the tests
# that need chemistry. Same pattern as similarity.py.
_STANDARDIZERS: tuple | None = None


def _standardizers() -> tuple:
    """RDKit's standardiser objects, built once on first use.

    Cached rather than constructed per call: a TautomerEnumerator is expensive
    to build and normalize_chemical runs once per compound.
    """
    global _STANDARDIZERS
    if _STANDARDIZERS is None:
        from rdkit import RDLogger
        from rdkit.Chem.MolStandardize import rdMolStandardize

        RDLogger.DisableLog("rdApp.*")
        _STANDARDIZERS = (rdMolStandardize.Uncharger(),
                          rdMolStandardize.LargestFragmentChooser(),
                          rdMolStandardize.TautomerEnumerator())
    return _STANDARDIZERS


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
    from rdkit import Chem, RDLogger
    from rdkit.Chem import inchi

    RDLogger.DisableLog("rdApp.*")

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
    from rdkit import Chem
    from rdkit.Chem import inchi
    from rdkit.Chem.MolStandardize import rdMolStandardize

    uncharger, largest_fragment, tautomers = _standardizers()

    mol = _parse(structure)

    n_frags_before = len(Chem.GetMolFrags(mol))
    charge_before = Chem.GetFormalCharge(mol)

    mol = rdMolStandardize.Cleanup(mol)
    mol = largest_fragment.choose(mol)
    mol = uncharger.uncharge(mol)

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
            taut = tautomers.Canonicalize(mol)
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


#: Length of an InChIKey's first block. The grouping key: salts,
#: stereoisomers and tautomers of one compound share it, which is what
#: same_as_rule below says in words. src/kgav/gtop.py has its own
#: SKELETON_LEN for the same 14; that is a second constant for one fact.
SKELETON_LEN = 14


def skeleton_of_id(node_id: str) -> str | None:
    """The InChIKey skeleton an `INCHIKEY:` node id carries, or None.

    WHY THIS IS NEEDED AWAY FROM A ChemId. similarity.py compares node ids
    and has no ChemId to hand, and it treated id inequality as proof that
    two fingerprints are two compounds. They are not: the graph uses the
    FULL InChIKey as the id (see the module docstring) precisely so a salt
    and its parent, or two stereoisomers, stay separate nodes -- and 434
    of 3,384 CHEMICALLY_SIMILAR_TO edges in v0.1 joined two such keys,
    reporting a compound as its own neighbour at tanimoto 1.0.

    gtop.py uses the same 14 characters as PROOF two keys are one compound
    while similarity.py used id inequality as proof they are two. The
    disagreement is what this closes; one definition is how it stays
    closed.

    None for an id with no InChIKey in it -- a locally minted `KGAV:` node,
    for instance -- because nothing can be concluded about its skeleton,
    and a caller must not read that as "different skeleton".
    """
    prefix = "INCHIKEY:"
    if not node_id.startswith(prefix):
        return None
    key = node_id[len(prefix):]
    return key[:SKELETON_LEN] if len(key) >= SKELETON_LEN else None


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


def same_structure(a: str, b: str) -> bool:
    """Do two SMILES denote one molecule, written differently?

    WHAT THIS IS FOR. assemble reported `property_conflict:smiles 7` on
    v0.1, and all seven were spelling: one layer writes aromatic rings
    lowercase and the other Kekule, atom order differs, and a SMILES stereo
    descriptor is ORDER-DEPENDENT -- the same centre reads @ or @@ according
    to which neighbour comes first. Every pair canonicalises to one molecule
    and both members match the node's own InChIKey.

    DELIBERATELY NOT normalize_chemical. That collapses salts, tautomers and
    (via the skeleton key) stereoisomers, which is the right grouping for a
    node id and far too strong here: two tautomers really are two values and
    should be reported, not merged away. Canonical SMILES normalises
    aromaticity, atom order and stereo NOTATION while keeping the
    stereochemistry itself, which is exactly the question "is this the same
    molecule spelled differently".

    AND CANONICAL SMILES WAS NOT ENOUGH (fix-list B1 found it). It does not
    normalise CHARGE-SEPARATED hypervalent sulfur: a sulfone written
    `S(=O)(=O)` and the same sulfone written `[S+]([O-])(=O)` canonicalise
    differently while InChI treats them as one molecule. On the B1 rebuild
    the hosttargets layer began supplying SMILES for 1,988 compounds ncats
    had already filed, and 73 pairs were reported as conflicts -- of which
    **72 were this**, every one with BOTH spellings matching the node's own
    InChIKey. `diagnose_property_conflicts.py` labelled them "DIFFERENT
    MOLECULES under one InChIKey ... an identifier bug" while printing the
    evidence against that verdict two lines above.

    TWO SIMPLER FIXES WERE MEASURED AND REJECTED. rdkit's
    rdMolStandardize.Normalize fixes a charge-separated SULFOXIDE and not a
    sulfone, so it does not cover the case found. InChIKey equality covers
    it but collapses amide/imidol tautomers (`CC(=O)N` and `CC(O)=N` share
    one key), which is exactly the merging this function must not do. So the
    fallback is a targeted recombination of S+/O- pairs and nothing else:
    it cannot reach a tautomer, a stereocentre, a salt or a double bond,
    and it deliberately leaves N-oxides alone, since neutralising those
    would invent a pentavalent nitrogen.

    Returns False for anything rdkit cannot parse. A value that is not a
    structure is a difference like any other and belongs in the report.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    ma, mb = Chem.MolFromSmiles(a), Chem.MolFromSmiles(b)
    if ma is None or mb is None:
        return False
    if Chem.MolToSmiles(ma) == Chem.MolToSmiles(mb):
        return True
    na, nb = _recombine_sulfur_oxygen(ma), _recombine_sulfur_oxygen(mb)
    if na is None or nb is None:
        return False
    return Chem.MolToSmiles(na) == Chem.MolToSmiles(nb)


#: `[S+][O-]` written back as `S=O`. Only sulfur, and only this bond: see
#: same_structure for the two broader fixes that were measured and rejected.
_SO_SMARTS = "[S+][O-]"
_SO_RXN = "[S+:1][O-:2]>>[S+0:1]=[O+0:2]"


def _recombine_sulfur_oxygen(mol):
    """Neutralise charge-separated S-O pairs, or None if that breaks the mol.

    Loops because a sulfone has two such pairs and a reaction applies one
    match at a time. Bounded rather than `while`: a transform that failed to
    reduce the match count would otherwise spin forever on some input nobody
    has seen yet.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    patt, rxn = Chem.MolFromSmarts(_SO_SMARTS), AllChem.ReactionFromSmarts(_SO_RXN)
    for _ in range(8):
        if not mol.GetSubstructMatches(patt):
            return mol
        out = rxn.RunReactants((mol,))
        if not out:
            return mol
        mol = out[0][0]
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return None
    return mol
