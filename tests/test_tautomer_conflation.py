"""One node, two compounds: L19 made countable at the point of collapse.

FIX-LIST L19. `normalize_chemical` canonicalises the tautomer before taking
the InChIKey, and standard InChI's mobile-H layer treats annular tautomers
as one compound regardless, so two library members differing only in where a
proton sits land on ONE node and pool their screening results.

WHY IT STAYED OPEN. A layer writes one node per id -- measured, every
per-layer nodes.jsonl has zero ids carrying two SMILES -- so the collapse is
invisible in the output and the validator cannot see it either. It is only
visible from the raw rows. The row said so and said the query had not been
run; run on ncats it returns 8.

WHAT IS NOT COUNTED, because the first version of the diagnostic counted it
and answered a broader question than the row asked:

  a salt or charge variant  `normalize_chemical` strips counterions and
                            neutralises before keying, so `Cl.X` and `X`
                            really are one parent and pooling is CORRECT.
                            That is L18 -- a wrong stored `smiles`, 108 of
                            them in ncats -- and not a conflation.
  one source, two spellings 3 of the 11 ncats tautomer pairs are a single
                            NCGC sample written inconsistently across files.
                            The node holds one compound; only the spelling
                            is wrong.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.ncats import substance_id  # noqa: E402
from kgav.normalize_chem import tautomer_conflations  # noqa: E402

#: 4-methylimidazole's two annular tautomers: one InChIKey, two SMILES,
#: and `same_structure` keeps them apart because collapsing tautomers is
#: the one thing its docstring says it must not do.
TAUT_A = "Cc1cnc(C)[nH]1"
TAUT_B = "Cc1c[nH]c(C)n1"


def test_two_source_compounds_on_one_node_is_a_conflation():
    got = tautomer_conflations({"INCHIKEY:A": {TAUT_A: {"S1"}, TAUT_B: {"S2"}}})
    assert got == {"INCHIKEY:A": 2}


def test_one_source_written_two_ways_is_not():
    """The node holds ONE compound. Reporting it would inflate L19 from 8
    to 11 and the extra 3 would not be conflations."""
    assert tautomer_conflations(
        {"INCHIKEY:A": {TAUT_A: {"S1"}, TAUT_B: {"S1"}}}) == {}


def test_a_salt_against_its_parent_is_not():
    """L18's class. normalize_chemical strips the counterion, so both rows
    ARE the same parent compound and pooling their results is correct."""
    assert tautomer_conflations(
        {"INCHIKEY:A": {"CN(C)CCO": {"S1"}, "Cl.CN(C)CCO": {"S2"}}}) == {}


def test_a_charged_form_against_its_parent_is_not():
    """Same reason, via `was_charged` rather than `had_multiple_fragments`:
    uncharge runs before the key is taken."""
    assert tautomer_conflations(
        {"INCHIKEY:A": {"CC(=O)O": {"S1"}, "CC(=O)[O-]": {"S2"}}}) == {}


def test_one_spelling_is_never_a_conflation():
    assert tautomer_conflations({"INCHIKEY:A": {TAUT_A: {"S1", "S2"}}}) == {}


def test_a_spelling_difference_that_is_one_molecule_is_not():
    """Aromatic vs Kekule, and atom order. same_structure calls these one
    molecule, which is what they are -- the seven `property_conflict:smiles`
    cases on v0.1 were all of this kind."""
    assert tautomer_conflations(
        {"INCHIKEY:A": {"c1ccccc1O": {"S1"}, "OC1=CC=CC=C1": {"S2"}}}) == {}


def test_a_salt_spelling_alongside_a_tautomer_pair_does_not_suppress_it():
    """THE BUG THIS FUNCTION SHIPPED WITH, for one hour on 2026-11-11.

    `salty` was computed over every representative and the whole key was
    then skipped, so {tautomer A, tautomer B, salt of A} reported nothing --
    even though A against B is a real conflation and one salt spelling does
    not make it go away.

    THE TELL WAS A COUNT THAT FELL WHEN ROWS WERE ADDED. Over the three raw
    ncats files it returned 8; over the ONE file the layer actually ingests
    it returned 11. The larger population carried more salt spellings, so
    more keys were suppressed, and a subset read higher than its superset.
    """
    taut_and_salt = {"K": {TAUT_A: {"S1"}, TAUT_B: {"S2"}, "Cl." + TAUT_A: {"S3"}}}
    assert tautomer_conflations(taut_and_salt) == {"K": 2}


def test_three_distinct_members_are_counted_as_three():
    got = tautomer_conflations({"INCHIKEY:A": {TAUT_A: {"S1", "S3"}, TAUT_B: {"S2"}}})
    assert got == {"INCHIKEY:A": 3}


# ----------------------------------------------------- the ingest wiring
def test_the_ingest_counts_them_rather_than_merging_silently():
    """ingest_screens must REPORT the collapse it performs. Before this the
    stats counter had no key for it and the eight nodes were indistinguishable
    from the other 5,226."""
    src = (ROOT / "src/kgav/ncats.py").read_text()
    assert "tautomer_conflations(derivations)" in src, (
        "the ncats ingest no longer calls the counter")
    assert 'stats["tautomer_conflated_nodes"]' in src
    assert 'stats["tautomer_conflated_members"]' in src


def test_the_count_has_one_owner():
    """THE DIAGNOSTIC MUST NOT RE-IMPLEMENT IT. Two implementations of one
    question agreeing is not corroboration when they are the same
    implementation twice -- and worse when they are not, because then the
    build and the diagnostic disagree about how many there are."""
    src = (ROOT / "scripts/diagnose_tautomer_conflation.py").read_text()
    assert "tautomer_conflations" in src, (
        "the diagnostic must defer to kgav.normalize_chem.tautomer_conflations")


# ------------------------------------------------------- the real input
@pytest.mark.slow
@pytest.mark.skipif(not (ROOT / "data/raw/ncats/cpe_assay.csv").exists(),
                    reason="data/raw is gitignored and hand-populated")
def test_the_measured_count_on_the_real_ncats_input():
    """9 nodes, 18 members, on the raw files as of 2026-11-11. This is the
    figure the row quotes; if the input changes it should fail here rather
    than age in the document.

    THE SOURCE IS A SUBSTANCE, VIA ncats.substance_id. NCGC sample ids are
    <substance>-<lot> and 815 of 8,596 base ids appear under more than one
    lot, so counting lots gives 12 nodes rather than 9 and calls two batches
    of one compound a conflation.

    MARKED slow BECAUSE IT IS, and the marker's docstring in pyproject.toml
    is about exactly this: it normalises 31,305 raw SMILES through RDKit's
    tautomer canonicaliser, which takes about nine minutes. One such test
    held half the wall time once already. `make test` excludes it and
    `make test-all` runs it.
    """
    import csv
    derivations: dict[str, dict[str, set[str]]] = {}
    from kgav.normalize_chem import normalize_chemical
    cache: dict[str, object] = {}
    for name in ("cpe_assay.csv", "cpe_tox.csv", "cpe.tsv"):
        p = ROOT / "data/raw/ncats" / name
        if not p.exists():
            continue
        delim = "\t" if p.suffix == ".tsv" else ","
        with p.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh, delimiter=delim):
                keys = {k.lower().replace("_", ""): k for k in row if k}
                sk, ik = keys.get("smiles"), keys.get("sampleid")
                if not sk:
                    continue
                smi = (row.get(sk) or "").strip()
                if not smi:
                    continue
                if smi not in cache:
                    try:
                        cache[smi] = normalize_chemical(smi)
                    except Exception:
                        cache[smi] = None
                cid = cache[smi]
                if cid is None:
                    continue
                derivations.setdefault(f"INCHIKEY:{cid.inchikey}", {}) \
                    .setdefault(smi, set()).add(
                        substance_id((row.get(ik) or "").strip()))
    got = tautomer_conflations(derivations)
    assert len(got) == 9, f"{len(got)} conflated nodes, expected 9: {sorted(got)}"
    assert sum(got.values()) == 18, sum(got.values())
