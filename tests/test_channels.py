"""Which channels a missing input blocks, derived rather than listed.

The point of this module is that "what is still needed" was prose in a comment
that said nothing about WHAT each item blocked -- so ORCS_KEY and the
VirHostNet downloads read as blockers on the whole flavivirus effort, when
neither touches M1, the direct-acting channel and the only one this project
has a defensible result on.
"""
from __future__ import annotations

import json

import pytest
import yaml

from kgav.channels import (blocked_by, channel_layers, layer_predicates,
                           predicates_in)

# The layer -> predicate map as assemble_release reported it on the
# workstation, 2026-10-11 (the connectivity-by-predicate block).
REAL_LAYERS = {
    "spine": {"ENCODED_BY", "BELONGS_TO", "CAUSES"},
    "host": {"PARTICIPATES_IN", "PHYSICALLY_INTERACTS_WITH", "ENCODED_BY"},
    "vhppi": {"PHYSICALLY_INTERACTS_WITH"},
    "orcs": {"HOST_FACTOR_FOR"},
    "chembl": {"INHIBITS", "MEASURED_INACTIVE_AGAINST",
               "HAS_ANTIVIRAL_ACTIVITY_AGAINST"},
    "selectivity": {"TARGETS"},
    "hosttargets": {"TARGETS"},
    "similarity": {"CHEMICALLY_SIMILAR_TO"},
    "ncats": {"HAS_ANTIVIRAL_ACTIVITY_AGAINST"},
}


@pytest.fixture(scope="module")
def metapaths():
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    return yaml.safe_load((root / "schema" / "kg_schema.yaml").read_text())[
        "metapaths"]


def test_the_reverse_marker_is_stripped_by_the_schemas_own_parser(metapaths):
    """BUG THIS CAUGHT IN ITS OWN FIRST VERSION. M4's path contains
    "<PARTICIPATES_IN" -- walk it backwards, Pathway to Protein. A local
    startswith("<") was omitted, so the predicate was compared verbatim
    against the layers and reported as one no layer emits, while
    PARTICIPATES_IN sits in the host layer with 157,211 edges.

    schema.parse_hop does the stripping now, because the reverse marker is the
    schema's convention and a second copy of it is a second thing to keep in
    step.
    """
    assert predicates_in(metapaths["M4"]["path"]) == [
        "TARGETS", "PARTICIPATES_IN", "PARTICIPATES_IN", "HOST_FACTOR_FOR"]
    assert not any(p.startswith("<")
                   for spec in metapaths.values()
                   for p in predicates_in(spec["path"]))


def test_m1_needs_neither_orcs_nor_virhostnet(metapaths):
    """THE FINDING. M1 is
    SmallMolecule -INHIBITS-> Protein -ENCODED_BY-> Gene -BELONGS_TO->
    OrganismTaxon, so it needs the chembl and spine layers and nothing from
    the host-mediated sources.

    That is why a flavivirus M1 evaluation is blocked only by the proteome
    file: ORCS_KEY and the VirHostNet downloads gate M2 through M6.
    """
    needs = channel_layers(metapaths, REAL_LAYERS)["M1"]
    assert needs["unsupplied"] == set()
    assert "orcs" not in needs["layers"]
    assert "vhppi" not in needs["layers"]
    assert {"chembl", "spine"} <= needs["layers"]


def test_only_the_host_mediated_channels_are_blocked_without_orcs(metapaths):
    """Stated as the whole picture rather than one channel, so a schema change
    that moves a predicate shows up here."""
    blocked = blocked_by(metapaths, REAL_LAYERS, {"orcs", "vhppi"})
    assert blocked["M1"] == set(), "direct-acting needs neither"
    assert blocked["M7"] == set(), "similarity-to-direct needs neither"
    for code in ("M2", "M4", "M6"):
        assert "orcs" in blocked[code], code
    assert "vhppi" in blocked["M3"]
    assert {"orcs", "vhppi"} <= blocked["M5"]


def test_an_unsupplied_predicate_is_a_blocker_not_a_clear_result(metapaths):
    """THE FALSE GREEN THIS MODULE PRODUCED ABOUT ITSELF.

    blocked_by first returned only `layers & unavailable`. M8's
    MEMBER_OF_CLASS and FOLD_SIMILAR_TO are emitted by NO layer -- ingest_folds
    writes a layer LAYER_PRECEDENCE does not list -- so M8 came back as an
    empty set and read as "evaluable today" while being unwalkable. A channel
    that cannot be traversed at all, reported as ready.

    The two causes stay distinguishable in the output: an unavailable layer is
    something to fetch, an unsupplied predicate is something to wire up.
    """
    blocked = blocked_by(metapaths, REAL_LAYERS, set())
    assert blocked["M1"] == set()
    assert blocked["M8"], "M8 is not evaluable and must not report as clear"
    assert any("FOLD_SIMILAR_TO" in b for b in blocked["M8"])
    assert any("no layer emits" in b for b in blocked["M8"])


def test_every_declared_metapath_is_classified(metapaths):
    """A metapath missing from the report would read as no constraint at all.
    EVERY DECLARED METAPATH IS REPORTED is already the rule in
    baselines.py:287; this is the same rule for the requirements view."""
    got = set(channel_layers(metapaths, REAL_LAYERS))
    assert got == set(metapaths), got ^ set(metapaths)


def test_layer_predicates_reads_what_a_layer_actually_emitted(tmp_path):
    """Measured from edges.jsonl rather than declared, so a layer that ran and
    produced nothing is distinguishable from one that was never run."""
    (tmp_path / "v0.1-spine").mkdir()
    (tmp_path / "v0.1-spine" / "edges.jsonl").write_text(
        json.dumps({"predicate": "ENCODED_BY"}) + "\n"
        + json.dumps({"predicate": "BELONGS_TO"}) + "\n")
    (tmp_path / "v0.1-orcs").mkdir()          # ran, wrote nothing
    (tmp_path / "v0.1-orcs" / "edges.jsonl").write_text("")

    got = layer_predicates(tmp_path, "v0.1")
    assert got["spine"] == {"ENCODED_BY", "BELONGS_TO"}
    assert got["orcs"] == set(), "an empty layer is not a missing layer"
    assert "vhppi" not in got, "a layer never run is absent, not empty"


def test_a_layer_that_was_never_run_blocks_its_channels(tmp_path, metapaths):
    """The case that matters for a promotion: orcs absent entirely, because
    ORCS_KEY was unset and ingest_orcs exited without writing."""
    (tmp_path / "v0.1-spine").mkdir()
    (tmp_path / "v0.1-spine" / "edges.jsonl").write_text(
        json.dumps({"predicate": "ENCODED_BY"}) + "\n"
        + json.dumps({"predicate": "BELONGS_TO"}) + "\n")
    (tmp_path / "v0.1-chembl").mkdir()
    (tmp_path / "v0.1-chembl" / "edges.jsonl").write_text(
        json.dumps({"predicate": "INHIBITS"}) + "\n")

    preds = layer_predicates(tmp_path, "v0.1")
    info = channel_layers(metapaths, preds)
    assert info["M1"]["unsupplied"] == set(), "M1 is walkable on these two"
    assert "HOST_FACTOR_FOR" in info["M2"]["unsupplied"], \
        "with no orcs layer, M2's predicate has no source at all"
