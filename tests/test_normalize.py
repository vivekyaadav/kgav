"""Day 2 acceptance gate: round-trip on hand-checked identifiers.

Reference files live in data/raw/refs and are NOT committed. Tests that need
them skip when absent, so CI stays green on a clean checkout while still
running for real on the workstation.

The property under test is not "does it return something" — it is "does it
either return the right thing, or say why it could not". A resolver that
silently passes bad input through is the failure this whole module exists to
prevent.
"""
from pathlib import Path

import pytest

from kgav.normalize_bio import MondoResolver, TaxonResolver, UniProtResolver
from kgav.normalize_chem import (
    NormalizationError,
    normalize_chemical,
    same_as_rule,
    same_structure,
)

REFS = Path(__file__).resolve().parents[1] / "data" / "raw" / "refs"
needs_refs = pytest.mark.skipif(not REFS.exists(), reason="reference files not downloaded")


# ----------------------------------------------------------------- chemicals
SALT_PAIRS = [
    ("nicotine",     "CN1CCC[C@H]1c1cccnc1",                  "CN1CCC[C@H]1c1cccnc1.Cl.Cl"),
    ("acetate",      "CC(=O)O",                               "CC(=O)[O-].[Na+]"),
    ("benzoate",     "OC(=O)c1ccccc1",                        "[O-]C(=O)c1ccccc1.[K+]"),
    ("lidocaine",    "CCN(CC)CC(=O)Nc1c(C)cccc1C",            "CCN(CC)CC(=O)Nc1c(C)cccc1C.Cl"),
    ("metformin",    "CN(C)C(=N)N=C(N)N",                     "CN(C)C(=N)N=C(N)N.Cl"),
]


@pytest.mark.parametrize("name,free,salt", SALT_PAIRS, ids=[p[0] for p in SALT_PAIRS])
def test_salt_forms_collapse_to_parent(name, free, salt):
    a, b = normalize_chemical(free), normalize_chemical(salt)
    assert a.inchikey == b.inchikey, f"{name}: salt did not collapse to parent"
    assert b.had_multiple_fragments, f"{name}: multi-fragment input not flagged"
    assert same_as_rule(a, b) == "identical_inchikey"


def test_charged_form_is_neutralised():
    neutral = normalize_chemical("CC(=O)O")
    anion = normalize_chemical("CC(=O)[O-]")
    assert neutral.inchikey == anion.inchikey
    assert anion.was_charged


def test_stereoisomers_are_distinct_but_share_skeleton():
    r = normalize_chemical("C[C@H](N)C(=O)O")
    s = normalize_chemical("C[C@@H](N)C(=O)O")
    assert r.inchikey != s.inchikey, "stereoisomers must not collapse"
    assert r.inchikey_skel == s.inchikey_skel
    assert same_as_rule(r, s) == "shared_inchikey_skeleton:salt_or_stereo_or_tautomer"


def test_inchi_and_smiles_agree():
    from rdkit import Chem
    from rdkit.Chem import inchi as rdinchi
    smiles = "CC(=O)Oc1ccccc1C(=O)O"          # aspirin
    as_inchi = rdinchi.MolToInchi(Chem.MolFromSmiles(smiles))
    assert normalize_chemical(smiles).inchikey == normalize_chemical(as_inchi).inchikey


def test_idempotent():
    """Normalising an already-normalised structure must be a no-op."""
    once = normalize_chemical("CN1CCC[C@H]1c1cccnc1.Cl.Cl")
    twice = normalize_chemical(once.smiles)
    assert once.inchikey == twice.inchikey


@pytest.mark.parametrize("bad", ["", "   ", "not_a_smiles!!!", "C1CC"])
def test_unparseable_structures_raise(bad):
    with pytest.raises(NormalizationError):
        normalize_chemical(bad)


# ------------------------------------------------------------------ proteins
@pytest.fixture(scope="module")
def uniprot():
    if not REFS.exists():
        pytest.skip("reference files not downloaded")
    return UniProtResolver(REFS / "sec_ac.txt", REFS / "delac_sp.txt")


@needs_refs
@pytest.mark.parametrize("acc", [
    "P0DTD1",   # SARS-CoV-2 replicase polyprotein 1ab
    "P0DTC2",   # SARS-CoV-2 spike
    "Q9BYF1",   # human ACE2
    "O15393",   # human TMPRSS2
    "P08183",   # human ABCB1
    "P05067",   # human APP
])
def test_known_primary_accessions_resolve(uniprot, acc):
    r = uniprot.resolve(acc)
    assert r.status == "primary", f"{acc}: {r.status} {r.note}"
    assert r.curie == f"UniProtKB:{acc}"


@needs_refs
def test_isoform_suffix_is_stripped(uniprot):
    assert uniprot.resolve("O15393-2").curie == "UniProtKB:O15393"


@needs_refs
def test_curie_input_is_accepted(uniprot):
    assert uniprot.resolve("UniProtKB:P0DTD1").curie == "UniProtKB:P0DTD1"


@needs_refs
@pytest.mark.parametrize("bad", ["NOTANACC", "12345", "", "P0DTD1XXXX"])
def test_malformed_accessions_are_rejected(uniprot, bad):
    r = uniprot.resolve(bad)
    assert not r.ok and r.status == "malformed"


@needs_refs
def test_secondary_accessions_resolve_to_a_primary(uniprot):
    """Sample the real secondary map — every entry must land on a primary."""
    sample = [a for a, p in list(uniprot.secondary.items())[:200] if len(p) == 1]
    assert sample, "no single-target secondary accessions found"
    for acc in sample:
        r = uniprot.resolve(acc)
        assert r.status == "secondary" and r.ok, f"{acc}: {r.status}"


@needs_refs
def test_deleted_accessions_are_not_silently_passed(uniprot):
    for acc in list(uniprot.deleted)[:50]:
        r = uniprot.resolve(acc)
        assert not r.ok and r.status == "deleted"


# --------------------------------------------------------------------- taxa
@pytest.fixture(scope="module")
def taxon():
    if not REFS.exists():
        pytest.skip("reference files not downloaded")
    return TaxonResolver(REFS / "taxdump.tar.gz")


@needs_refs
@pytest.mark.parametrize("taxid,expect", [
    ("2697049", "Severe acute respiratory syndrome coronavirus 2"),
    ("694009",  "Severe acute respiratory syndrome-related coronavirus"),
    ("1335626", "Middle East respiratory syndrome-related coronavirus"),
    ("11118",   "Coronaviridae"),
    ("9606",    "Homo sapiens"),
])
def test_coronavirus_taxa_resolve(taxon, taxid, expect):
    r = taxon.resolve(taxid)
    assert r.status == "primary", f"{taxid}: {r.status} {r.note}"
    assert expect.lower() in r.note.lower(), f"{taxid}: got {r.note!r}"


@needs_refs
def test_nonexistent_taxon_is_unmapped(taxon):
    """The bug that shipped in the first draft: bogus ids reported primary."""
    r = taxon.resolve("999999999")
    assert not r.ok and r.status == "unmapped"


@needs_refs
def test_merged_taxa_follow_to_current_id(taxon):
    sample = list(taxon.merged.items())[:200]
    assert sample, "merged.dmp empty"
    for old, new in sample:
        r = taxon.resolve(old)
        assert r.status == "secondary" and r.curie == f"NCBITaxon:{new}"


@needs_refs
@pytest.mark.parametrize("bad", ["abc", "", "NCBITaxon:xyz"])
def test_malformed_taxa_are_rejected(taxon, bad):
    assert taxon.resolve(bad).status == "malformed"


# ------------------------------------------------------------------ diseases
@pytest.fixture(scope="module")
def mondo():
    if not REFS.exists():
        pytest.skip("reference files not downloaded")
    return MondoResolver(REFS / "mondo.json")


@needs_refs
@pytest.mark.parametrize("xref", ["MONDO:0100096", "DOID:0080600", "UMLS:C5203670", "ICD10CM:U07.1"])
def test_covid_xrefs_all_land_on_one_term(mondo, xref):
    r = mondo.resolve(xref)
    assert r.curie == "MONDO:0100096", f"{xref}: {r.status} {r.note}"


@needs_refs
def test_unknown_xref_is_unmapped_not_guessed(mondo):
    r = mondo.resolve("NOPE:123")
    assert not r.ok and r.status == "unmapped"


@needs_refs
def test_labels_are_never_used_for_matching(mondo):
    """Passing a disease NAME must fail. Label matching is banned by design."""
    r = mondo.resolve("COVID-19")
    assert not r.ok, "resolver matched on a label string"


@needs_refs
def test_obsolete_terms_are_flagged(mondo):
    if not mondo.obsolete:
        pytest.skip("no obsolete terms parsed")
    r = mondo.resolve(next(iter(mondo.obsolete)))
    assert not r.ok and r.status == "deleted"


# ------------------------------------------------- spelling vs disagreement
def test_same_structure_sees_through_spelling():
    """The real v0.1 case: chembl writes aromatic rings lowercase, ncats
    writes them Kekule, the atom order differs, and the stereo descriptor
    flips with neighbour order. One molecule (rolipram), two spellings."""
    assert same_structure(
        "COc1ccc([C@@H]2CNC(=O)C2)cc1OC1CCCC1",
        "COC1=C(OC2CCCC2)C=C(C=C1)[C@@H]3CNC(=O)C3")


def test_same_structure_keeps_stereochemistry():
    """DELIBERATELY NOT normalize_chemical, whose skeleton key groups
    stereoisomers. Two enantiomers are two compounds and must not be
    called one value spelled two ways."""
    assert not same_structure("C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O")


def test_same_structure_is_false_for_different_molecules_and_for_junk():
    """A value rdkit cannot parse is a difference like any other: it belongs
    in the conflict report, not silently merged."""
    assert not same_structure("CCO", "c1ccccc1")
    assert not same_structure("not a smiles !!!", "CCO")
    assert not same_structure("CCO", "")


# ------------------- charge-separated sulfur, found by fix-list B1's rebuild
def test_same_structure_sees_through_a_charge_separated_sulfonyl():
    """CANONICAL SMILES WAS NOT ENOUGH. A sulfone written S(=O)(=O) and the
    same sulfone written [S+]([O-])(=O) canonicalise differently while InChI
    treats them as one molecule. When B1 let hosttargets supply SMILES for
    1,988 compounds ncats had already filed, 73 pairs were reported as
    conflicts and 36 of them were this -- every one with BOTH spellings
    matching the node's own InChIKey, which is the evidence that they are
    one molecule."""
    assert same_structure(
        "Cc1cc(C)c(S(=O)(=O)n2cc(C)cc2C)c(C)c1",
        "CC1=C[N](C(=C1)C)[S+]([O-])(=O)C2=C(C)C=C(C)C=C2C")


def test_same_structure_handles_BOTH_sulfonyls_in_one_molecule():
    """The transform rewrites one match at a time, so a molecule with two
    charge-separated sulfonyls needs more than one pass. A single-shot fix
    would quietly leave these reported as conflicts."""
    assert same_structure(
        "CS(=O)(=O)c1cc(-c2ccccc2)cc(S(=O)(=O)c2ccc(CN)s2)c1",
        "C[S+]([O-])(=O)C1=CC(=CC(=C1)C2=CC=CC=C2)[S+]([O-])(=O)C3=CC=C(CN)S3")


def test_same_structure_sees_through_a_charge_separated_sulfoxide():
    assert same_structure("CS(=O)c1ccccc1", "C[S+]([O-])c1ccccc1")


def test_the_sulfur_fix_does_not_collapse_a_tautomer():
    """THE PROPERTY THAT MAKES THIS FIX NARROW ENOUGH TO BE SAFE. Two
    broader fixes were measured and rejected for failing exactly here:
    rdMolStandardize.Normalize does not cover a sulfone at all, and
    InChIKey equality covers it but gives amide and imidol one key. The
    S+/O- transform cannot reach a tautomer."""
    assert not same_structure("CC(=O)N", "CC(O)=N")
    assert not same_structure("CC(=O)CC(C)=O", "CC(O)=CC(C)=O")


def test_the_sulfur_fix_does_not_collapse_an_annular_tautomer():
    """THE 32 THAT REMAIN REPORTED, AND SHOULD BE. ncats and hosttargets
    disagree about which imidazole nitrogen carries the proton. Both
    spellings resolve to the node's InChIKey, because InChI's mobile-H layer
    treats them as one compound -- so the node cannot separate them and the
    SMILES can. Reporting that is the comparator working, not failing."""
    assert not same_structure(
        "Cc1cnc(-c2cnc(NCCNc3ccc(C#N)cn3)nc2-c2ccc(Cl)cc2Cl)[nH]1",
        "CC1=C[NH]C(=N1)C2=C(N=C(NCCNC3=CC=C(C=N3)C#N)N=C2)C4=CC=C(Cl)C=C4Cl")


def test_the_sulfur_fix_does_not_collapse_salts_stereo_or_geometry():
    """The other three things a looser comparator would have merged away."""
    assert not same_structure("CCN", "CCN.Cl")
    assert not same_structure("C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O")
    assert not same_structure("C/C=C/C", r"C/C=C\C")


def test_an_N_oxide_is_left_alone_by_the_sulfur_fix():
    """Neutralising [N+][O-] would invent a pentavalent nitrogen, so the
    transform names sulfur explicitly rather than 'any charge-separated
    pair'. An N-oxide compared with itself must still be equal, and the
    molecule must survive the pass."""
    no = "c1ccnc(c1)[N+]([O-])=O"
    assert same_structure(no, no)
    assert not same_structure(no, "c1ccnc(c1)N=O")
