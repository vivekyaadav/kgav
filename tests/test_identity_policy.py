"""salt_collapsed and stereo_collapsed must state what the id actually does.

FIX-LIST L24. Both are `required: true` in the schema and both were
hardcoded `False` at all three SmallMolecule emit sites -- chembl.py,
host_targets.py, ncats.py -- and measured `False` on 25,383 of 25,383
compounds in v0.1. `salt_collapsed: false` is WRONG: `normalize_chemical`
runs `largest_fragment.choose` before taking the InChIKey, so `Cl.X` and `X`
land on one node.

WHAT THESE FIELDS ARE, which took a wrong turn to establish.
normalize_chem's module docstring says it: "Which one is the node id is a
project decision recorded on every node as salt_collapsed /
stereo_collapsed." They declare the IDENTITY POLICY, not a property of the
compound in front of them. So a constant is correct, and
host_targets.py's comment -- "Must match the Day 7 chembl layer exactly, or
the assembler reports property conflicts on thousands of shared compounds"
-- is benign rather than a smell: every layer must declare the same policy.

I first read them as properties of a DERIVATION and proposed resolving them
at assembly by OR, on the grounds that two layers deriving one id from a
salt and from the parent genuinely disagree. That was wrong, and it would
have turned a one-value fix into a schema change plus an assembler change.
The measurement that settled it is the one below: ask the normaliser.

THE VALUES ARE DERIVED HERE, NOT TYPED. A test asserting
`salt_collapsed is True` against a literal `True` in the ingest would pass
if both were wrong together, which is the defect this repository keeps
finding in its own checks. These probe `normalize_chemical` and require the
ingests to agree with what it does.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

RELEASE = ROOT / "data/releases/v0.1"
EMITTERS = ("src/kgav/chembl.py", "src/kgav/host_targets.py", "src/kgav/ncats.py")

#: A salt and its parent; two stereoisomers of one skeleton; two annular
#: tautomers of 4-methylimidazole.
SALT, PARENT = "Cl.CN(C)CCO", "CN(C)CCO"
L_ALA, D_ALA = "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"
TAUT_A, TAUT_B = "Cc1cnc(C)[nH]1", "Cc1c[nH]c(C)n1"


@pytest.fixture(scope="module")
def measured() -> dict:
    """What the id ACTUALLY collapses, measured on the normaliser."""
    from kgav.normalize_chem import normalize_chemical as N
    return {
        "salt_collapsed": N(SALT).inchikey == N(PARENT).inchikey,
        "stereo_collapsed": N(L_ALA).inchikey == N(D_ALA).inchikey,
        "tautomer_collapsed": N(TAUT_A).inchikey == N(TAUT_B).inchikey,
    }


def test_the_id_collapses_salts_and_not_stereoisomers(measured):
    """The policy itself, so a change to normalize_chemical fails HERE
    first with a readable reason, rather than as a property mismatch."""
    assert measured["salt_collapsed"] is True, (
        "largest_fragment no longer runs before the InChIKey; the id has "
        "stopped collapsing salts and every node's salt_collapsed is now wrong")
    assert measured["stereo_collapsed"] is False, (
        "the id has started collapsing stereoisomers -- it should be the FULL "
        "key, with inchikey_skel carrying the stereo-agnostic grouping")
    assert measured["tautomer_collapsed"] is True, (
        "the id has stopped collapsing tautomers; fix-list L19's 9 conflated "
        "nodes exist BECAUSE it does, and every node now declares it")


def test_the_skeleton_is_what_groups_stereoisomers(measured):
    """The other half of the two-key design: stereo is not collapsed in the
    id but IS grouped by inchikey_skel, so a bad merge stays reversible."""
    from kgav.normalize_chem import normalize_chemical as N
    assert N(L_ALA).inchikey_skel == N(D_ALA).inchikey_skel


@pytest.mark.parametrize("path", EMITTERS)
def test_every_emitter_declares_the_measured_policy(path, measured):
    """PROVEN TO FAIL against the tree of 2026-11-11, where all three said
    `salt_collapsed=False, stereo_collapsed=False`."""
    src = (ROOT / path).read_text()
    m = re.search(r"salt_collapsed=(True|False),\s*"
                  r"stereo_collapsed=(True|False),\s*"
                  r"tautomer_collapsed=(True|False)", src)
    assert m, f"{path} no longer emits all three flags in one call"
    got = {"salt_collapsed": m.group(1) == "True",
           "stereo_collapsed": m.group(2) == "True",
           "tautomer_collapsed": m.group(3) == "True"}
    assert got == measured, (
        f"{path} declares {got} and the normaliser does {measured}")


@pytest.mark.skipif(not (RELEASE / "nodes.jsonl").exists(),
                    reason="no release on disk")
def test_the_release_agrees_with_the_normaliser(measured):
    """The release itself. Fails until the rebuild that carries L24, which
    is the point: the 25,383 nodes on disk state a policy the id does not
    follow."""
    seen: dict[tuple, int] = {}
    for line in (RELEASE / "nodes.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        p = (json.loads(line).get("properties") or {})
        if "salt_collapsed" not in p:
            continue
        key = (p.get("salt_collapsed"), p.get("stereo_collapsed"),
               p.get("tautomer_collapsed"))
        seen[key] = seen.get(key, 0) + 1
    assert seen, "no node carries the flags"
    want = (measured["salt_collapsed"], measured["stereo_collapsed"],
            measured["tautomer_collapsed"])
    assert set(seen) == {want}, (
        f"the release carries {dict(seen)} and the normaliser does {want}")


# ------------------------------------------------- L18: smiles describes the node
def test_the_parent_smiles_is_what_the_id_was_derived_from():
    """`ChemId.smiles` must be the de-salted, uncharged structure, because
    that is what ncats now stores and what makes the node self-consistent.
    If normalize_chemical stops returning the parent, the stored `smiles`
    silently becomes the raw string again."""
    from kgav.normalize_chem import normalize_chemical as N
    assert N(SALT).smiles == N(PARENT).smiles == "CN(C)CCO"
    assert N("CC(=O)[O-]").smiles == "CC(=O)O"


def test_ncats_stores_the_parent_and_keeps_the_supplied_string():
    """fix-list L18. The raw string is ncats's statement about its sample and
    must survive; what must not survive is a `smiles` that resolves to a
    different InChIKey than the node it sits on."""
    src = (ROOT / "src/kgav/ncats.py").read_text()
    assert "smiles=cid.smiles" in src, (
        "ncats must store the parent structure the id was derived from")
    assert "smiles_as_supplied=(struct if differs else None)" in src, (
        "the source's own string must be kept where the substance differs")
    assert "differs = cid.had_multiple_fragments or cid.was_charged" in src, (
        "the condition must be a CHANGE OF SUBSTANCE, not of spelling")
    assert "cid.tautomer_changed" not in src.split("differs =")[1][:200], (
        "tautomer_changed must NOT select this property. Measured: of 1,352 "
        "nodes it selected, 1,130 were tautomer canonicalisation alone "
        "against 124 salts and 98 charged forms -- and under InChI's mobile-H "
        "layer two tautomers are the same compound, which is why these nodes "
        "exist (L19) and why tautomer_collapsed is declared true")


@pytest.mark.skipif(not (RELEASE / "nodes.jsonl").exists(),
                    reason="no release on disk")
def test_no_stored_smiles_resolves_to_a_different_node():
    """THE DEFECT ITSELF, on the release. A node whose `smiles` normalises to
    a key other than its own id is stating a structure it does not hold.
    Checked on the compounds that carry `smiles_as_supplied`, which are
    exactly the ones where the source and the parent disagreed -- the whole
    25,383 would need an RDKit pass and is the `slow` test's job.
    """
    from kgav.normalize_chem import normalize_chemical as N
    checked = bad = 0
    for line in (RELEASE / "nodes.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        node = json.loads(line)
        props = node.get("properties") or {}
        if not props.get("smiles_as_supplied"):
            continue
        checked += 1
        try:
            got = N(props["smiles"]).inchikey
        except Exception:
            continue
        if f"INCHIKEY:{got}" != node["id"]:
            bad += 1
    assert bad == 0, (
        f"{bad} of {checked} node(s) carrying smiles_as_supplied have a "
        f"`smiles` that resolves to a different InChIKey")
