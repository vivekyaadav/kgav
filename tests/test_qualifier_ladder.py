"""Two witnesses to one qualifier, when precedence cannot decide it.

WHAT THIS IS FOR. `merge_edges` took the higher layer's value for every
qualifier: `merged_quals.update(a)`. For `direction` that is wrong in both
directions, and the Guide to PHARMACOLOGY download of 2026-11-01 measured
how wrong.

  49 edges would have been made WORSE. ChEMBL states `inhibitor` or
  `agonist`; GtoP states `modulator` for the same pair, because the
  compound binds an allosteric site. Both are right, `inhibitor` is the
  more useful of the two, and gtop sits above hosttargets in
  LAYER_PRECEDENCE -- so the vaguer value won, placed there by a comment
  claiming it knew more.

  AND THE GAP-FILLING NEVER WORKED BY PRECEDENCE. host_targets writes
  `direction: "unknown"` explicitly, because the schema requires the
  qualifier. "Fill the gaps, do not overwrite" is therefore not something
  a layer position can express: whichever layer is higher wins, including
  when it is the one with no claim.

So `direction` has a LADDER -- no claim, then graded, then specific -- and
the rules are a claim beats no claim, specific beats graded, and two
opposite specifics reconcile to neither.
"""
from __future__ import annotations

import itertools
import sys
from collections import Counter
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.assemble import (
    LAYER_PRECEDENCE,
    QUALIFIER_LADDER,
    merge_edges,
    resolve_qualifier,
)

DIRECTIONS = ("inhibitor", "agonist", "modulator", "unknown")


def _edge(direction, source, tier=1):
    return {"subject": "INCHIKEY:D", "predicate": "TARGETS",
            "object": "UniProtKB:P1", "qualifiers": {"direction": direction},
            "primary_knowledge_source": source, "evidence_tier": tier,
            "first_asserted_date": "2020-01-01"}


# ------------------------------------------------------------- the three rules
@pytest.mark.parametrize("other", ["inhibitor", "agonist", "modulator"])
def test_no_claim_never_wins(other):
    """Either way round. This is the rule that makes "gtop fills the gaps"
    work at all, and no layer ordering can express it."""
    assert resolve_qualifier("direction", "unknown", other) == other
    assert resolve_qualifier("direction", other, "unknown") == other


@pytest.mark.parametrize("specific", ["inhibitor", "agonist"])
def test_a_specific_direction_is_not_replaced_by_a_graded_one(specific):
    """THE 49. A negative allosteric modulator of ABL1 is both a NAM and an
    inhibitor; `modulator` is the more precise mechanism and `inhibitor` the
    more useful direction, and overwriting one with the other loses what
    M2, M4 and M5 read."""
    assert resolve_qualifier("direction", "modulator", specific) == specific
    assert resolve_qualifier("direction", specific, "modulator") == specific


def test_two_opposite_specifics_reconcile_to_neither():
    """THE 40. Two curated sources contradicting each other is not a fact,
    and keeping whichever sits higher in LAYER_PRECEDENCE would assert a
    mechanism on the strength of a list's ordering."""
    assert resolve_qualifier("direction", "inhibitor", "agonist") == "unknown"
    assert resolve_qualifier("direction", "agonist", "inhibitor") == "unknown"


@pytest.mark.parametrize("value", DIRECTIONS)
def test_two_witnesses_that_agree_pass_through(value):
    assert resolve_qualifier("direction", value, value) == value


# --------------------------------------------------------------- the invariants
@pytest.mark.parametrize("hi,lo", list(itertools.product(DIRECTIONS, repeat=2)))
def test_the_outcome_does_not_depend_on_which_witness_arrived_first(hi, lo):
    """ORDER INDEPENDENCE IS THE POINT, not tidiness. `primary_knowledge_source`
    once grew quadratically because merging was neither associative nor
    idempotent and the value depended on arrival order. Here the same flaw
    would be worse: one release built twice could state opposite directions
    on one edge, and the fingerprint would catch it without saying why."""
    assert (resolve_qualifier("direction", hi, lo)
            == resolve_qualifier("direction", lo, hi))


@pytest.mark.parametrize("hi,lo", list(itertools.product(DIRECTIONS, repeat=2)))
def test_every_outcome_is_a_value_the_schema_allows(hi, lo):
    enum = yaml.safe_load(
        (ROOT / "schema" / "kg_schema.yaml").read_text()
    )["enums"]["target_direction"]
    assert resolve_qualifier("direction", hi, lo) in enum


@pytest.mark.parametrize("hi,lo", list(itertools.product(DIRECTIONS, repeat=2)))
def test_reconciling_is_idempotent(hi, lo):
    """A third witness agreeing with the outcome must not move it."""
    once = resolve_qualifier("direction", hi, lo)
    assert resolve_qualifier("direction", once, once) == once


def test_a_qualifier_with_no_ladder_keeps_the_precedence_behaviour():
    """The ladder is opt-in. Every other qualifier must behave exactly as
    before, or this change reaches far beyond the one it was made for."""
    assert "assay_type" not in QUALIFIER_LADDER
    assert resolve_qualifier("assay_type", "binding", "biochemical") == "binding"
    assert resolve_qualifier("affinity_nm", 40.0, 90.0) == 40.0


def test_a_single_witness_is_unchanged():
    assert resolve_qualifier("direction", "inhibitor", None) == "inhibitor"
    assert resolve_qualifier("direction", None, "inhibitor") == "inhibitor"


def test_the_ladder_is_ordered_least_specific_first():
    """The rule reads a ladder, so a ladder written the other way round
    would invert every decision silently."""
    assert QUALIFIER_LADDER["direction"] == ("unknown", "modulator")


def test_every_ladder_value_is_in_the_schema_enum():
    enums = yaml.safe_load(
        (ROOT / "schema" / "kg_schema.yaml").read_text())["enums"]
    for name, ladder in QUALIFIER_LADDER.items():
        enum = enums.get(f"target_{name}") or enums.get(name)
        assert enum, f"no enum found for the {name!r} ladder"
        assert set(ladder) <= set(enum), (name, ladder, enum)


# ------------------------------------------------------- through merge_edges
def test_a_merged_edge_keeps_the_more_specific_direction():
    """End to end, in the real layer order: gtop sits ABOVE hosttargets, so
    its qualifiers are the `a` side."""
    assert LAYER_PRECEDENCE.index("gtop") < LAYER_PRECEDENCE.index("hosttargets")
    out = merge_edges(_edge("modulator", "infores:guidetopharmacology"),
                      _edge("inhibitor", "infores:chembl"))
    assert out["qualifiers"]["direction"] == "inhibitor"
    assert out["primary_knowledge_source"] == (
        "infores:chembl|infores:guidetopharmacology")


def test_a_merged_edge_fills_an_unknown_direction():
    out = merge_edges(_edge("inhibitor", "infores:guidetopharmacology"),
                      _edge("unknown", "infores:chembl"))
    assert out["qualifiers"]["direction"] == "inhibitor"


def test_a_merged_edge_refuses_an_outright_contradiction():
    out = merge_edges(_edge("agonist", "infores:guidetopharmacology"),
                      _edge("inhibitor", "infores:chembl"))
    assert out["qualifiers"]["direction"] == "unknown"


def test_both_reconciliations_are_counted_not_silent():
    """These reach the MANIFEST through Assembly.stats. A merge that
    quietly changed 89 directions would be indistinguishable from one that
    changed none."""
    st: Counter = Counter()
    merge_edges(_edge("modulator", "s1"), _edge("inhibitor", "s2"), stats=st)
    assert st["qualifier_kept_more_specific:direction"] == 1
    merge_edges(_edge("agonist", "s1"), _edge("inhibitor", "s2"), stats=st)
    assert st["qualifier_irreconcilable:direction"] == 1
    # Agreement is not an event worth counting; a counter that ticks on
    # every merge is one nobody reads.
    merge_edges(_edge("inhibitor", "s1"), _edge("inhibitor", "s2"), stats=st)
    assert sum(st.values()) == 2


def test_merging_without_a_counter_still_works():
    """Assembly passes one; merge_edges is used directly elsewhere."""
    out = merge_edges(_edge("modulator", "s1"), _edge("inhibitor", "s2"))
    assert out["qualifiers"]["direction"] == "inhibitor"


def test_a_qualifier_only_one_side_states_survives_the_merge():
    """The loop runs over `a`'s keys, so a key only `b` has must not be
    dropped -- which is how affinity_nm would vanish."""
    a = _edge("inhibitor", "s1")
    b = dict(_edge("unknown", "s2"))
    b["qualifiers"] = {"direction": "unknown", "affinity_nm": 40.0}
    out = merge_edges(a, b)
    assert out["qualifiers"] == {"direction": "inhibitor", "affinity_nm": 40.0}
