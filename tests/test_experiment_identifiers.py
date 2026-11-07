"""Every measurement edge says which experiment it came from.

Fix-list V6: no edge carried an assay, document or screen identifier,
though every layer has one at emit time. ChEMBL's assay CHEMBL id was
already in this ingest's own dedup key -- three separate key tuples are
built from it -- and BioGRID's SCREEN_ID was the loop variable. They
simply never reached an edge.

Three things need them:

  traceability   a reader cannot get from an edge back to the experiment
  fix-list E6    assemble settles 17,254 measurement conflicts by layer
                 precedence and could not say whether two disagreeing
                 potencies came from two assays or one
  fix-list L2    "keep each measurement as its own edge, keyed on the
                 assay" is not expressible until the assay is on the edge

NOT identity qualifiers. Making them so IS L2's decision and would split
every merged edge; this carries the fact so the decision can be made with
the data in front of it.

Coverage on the real inputs: chembl 31,569/31,569 both fields, hosttargets
18,954/18,954 assay_id, selectivity 3,226/3,226 document_id, orcs
11,112/11,112 screen_id.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from kgav.schema import load_schema

ROOT = Path(__file__).resolve().parents[1]
IDENTIFIER_PREDICATES = {
    "HAS_ANTIVIRAL_ACTIVITY_AGAINST": {"assay_id", "document_id"},
    "INHIBITS": {"assay_id", "document_id"},
    "MEASURED_INACTIVE_AGAINST": {"assay_id", "document_id"},
    "TARGETS": {"assay_id"},
    "HOST_FACTOR_FOR": {"screen_id"},
}


@pytest.mark.parametrize("predicate,names",
                         sorted(IDENTIFIER_PREDICATES.items()))
def test_the_identifier_is_declared(predicate, names):
    ec = next(e for e in load_schema().edge_classes
              if e.predicate == predicate)
    missing = names - set(ec.qualifiers)
    assert not missing, (
        f"{predicate} emits {sorted(missing)} and does not declare them; "
        f"the release would fail with UNDECLARED_QUALIFIER")


@pytest.mark.parametrize("predicate,names",
                         sorted(IDENTIFIER_PREDICATES.items()))
def test_an_identifier_is_not_an_identity_qualifier(predicate, names):
    """L2's decision, not this commit's. An identity qualifier splits every
    merged edge, which would change the graph rather than describe it."""
    ec = next(e for e in load_schema().edge_classes
              if e.predicate == predicate)
    ident = set(getattr(ec, "identity_qualifiers", None) or ())
    assert not (ident & names), (
        f"{predicate} keys its identity on {sorted(ident & names)}; that is "
        f"fix-list L2 and moves the graph")


def test_chembl_sets_them_before_the_branch():
    """Setting them inside the protein branch reached only 39% of the
    HAS_ANTIVIRAL edges: the direct-organism rows emit from an earlier
    branch and never saw the assignment. Measured, not assumed -- the
    first version of this fix shipped at 9,465 of 24,073."""
    src = (ROOT / "src" / "kgav" / "chembl.py").read_text()
    shared = src.index('quals = {qual_field: nm, "relation": relation}')
    branch = src.index('if r["target_type"] == "ORGANISM":', shared)
    between = src[shared:branch]
    assert 'quals["assay_id"]' in between
    assert 'quals["document_id"]' in between


def test_the_organism_allowlist_admits_them():
    """ORGANISM_QUALIFIERS is a deliberate allow-list that excludes
    protein-RESOLUTION fields, which say nothing once the edge is on the
    organism. An experiment identifier is not one of those -- and the
    organism edges are the held-out labels, so which assay produced a
    label is the question V6 exists to answer."""
    from kgav.chembl import ORGANISM_QUALIFIERS
    assert "assay_id" in ORGANISM_QUALIFIERS
    assert "document_id" in ORGANISM_QUALIFIERS
    # The exclusions it exists for are still excluded.
    assert "resolved_by" not in ORGANISM_QUALIFIERS


def test_the_identifier_is_read_not_invented():
    """Each comes from the source's own field. A minted or derived id would
    be traceable to nothing, which is the opposite of the point."""
    checks = {
        "src/kgav/chembl.py": 'str(r["assay"])',
        "src/kgav/host_targets.py": 'str(r["assay"])',
        "src/kgav/selectivity.py": 'str(r["doc_id"])',
        "src/kgav/orcs.py": '"screen_id": sid',
    }
    for path, needle in checks.items():
        assert needle in (ROOT / path).read_text(), f"{path}: {needle}"
