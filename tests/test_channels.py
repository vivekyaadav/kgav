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

    got, skipped = layer_predicates(tmp_path, "v0.1")
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

    preds, _ = layer_predicates(tmp_path, "v0.1")
    info = channel_layers(metapaths, preds)
    assert info["M1"]["unsupplied"] == set(), "M1 is walkable on these two"
    assert "HOST_FACTOR_FOR" in info["M2"]["unsupplied"], \
        "with no orcs layer, M2's predicate has no source at all"


# --------------------------------------------------------------------------
# "is the predicate wired" vs "does THIS virus have data there"
# --------------------------------------------------------------------------
def test_m3_cannot_reach_a_newly_promoted_virus_without_virhostnet(metapaths):
    """THE FALSE GREEN THIS FIXES.

    The dry run reported "M3 evaluable once the proteome lands" for Zika,
    because PHYSICALLY_INTERACTS_WITH is emitted by the host layer -- 150,432
    edges of it. But M3's path is

        SmallMolecule -TARGETS-> Protein -PHYSICALLY_INTERACTS_WITH->
        Protein -ENCODED_BY-> Gene -BELONGS_TO-> OrganismTaxon

    and the second Protein is the VIRAL one, because the Gene under it belongs
    to the virus. host's PPIs are human-human and cannot reach a viral taxon;
    only vhppi bridges host to virus, and it has no flavivirus files. So M3
    cannot reach Zika, and the report said it could.
    """
    from kgav.channels import reachable_for_virus

    got = reachable_for_virus(metapaths, REAL_LAYERS, {"spine", "chembl"})
    assert got["M3"]["entry"] == "PHYSICALLY_INTERACTS_WITH"
    assert got["M3"]["sources"] == {"host", "vhppi"}
    assert got["M3"]["reachable"] is False, \
        "host supplies the predicate but cannot reach a viral taxon"


def test_only_m1_and_m7_reach_a_virus_covered_by_spine_and_chembl(metapaths):
    """Zika's real situation after fetch_raw: proteome cached, 355 ChEMBL
    protein rows, no ORCS screens ingestable, no VirHostNet files."""
    from kgav.channels import reachable_for_virus

    got = reachable_for_virus(metapaths, REAL_LAYERS, {"spine", "chembl"})
    reaches = {c for c, i in got.items() if i["reachable"]}
    assert reaches == {"M1", "M7"}, reaches


def test_the_entry_predicate_is_the_first_non_spine_hop_from_the_taxon(metapaths):
    """Derived by walking the path backwards from OrganismTaxon: every hop the
    spine supplies stays inside the virus's own subgraph, so the first one it
    does not supply is where the virus must be entered from outside."""
    from kgav.channels import virus_entry_predicate

    spine = REAL_LAYERS["spine"]
    assert virus_entry_predicate(metapaths["M1"]["path"], spine) == "INHIBITS"
    assert virus_entry_predicate(metapaths["M7"]["path"], spine) == "INHIBITS"
    assert virus_entry_predicate(metapaths["M2"]["path"], spine) \
        == "HOST_FACTOR_FOR"
    # SCHEMA 0.13.0: M8 walks Protein -FOLD_SIMILAR_TO-> Protein instead of
    # out through TargetClass and back. The entry predicate moved with it,
    # and the new answer is the stricter one -- MEMBER_OF_CLASS is emitted
    # for every chain the folds layer resolves, so it would read "evaluable"
    # for a virus whose proteins merely have class memberships. Reaching a
    # virus through M8 requires a fold edge that actually lands on one of
    # ITS proteins, which is what this module exists to distinguish.
    assert virus_entry_predicate(metapaths["M8"]["path"], spine) \
        == "FOLD_SIMILAR_TO"


def test_adding_virhostnet_coverage_makes_m3_reachable(metapaths):
    """The counterfactual, so the test says what WOULD unblock it rather than
    only that it is blocked."""
    from kgav.channels import reachable_for_virus

    with_vhn = reachable_for_virus(metapaths, REAL_LAYERS,
                                   {"spine", "chembl", "vhppi"})
    assert with_vhn["M3"]["reachable"] is True
    assert with_vhn["M2"]["reachable"] is False, "M2 still needs orcs"


def test_every_metapath_has_an_entry_predicate(metapaths):
    """A path entirely spine-supplied would never need an outside source,
    which no declared metapath is -- each has to enter the virus somewhere. If
    one returned None the reachability verdict would be vacuous."""
    from kgav.channels import virus_entry_predicate

    for code, spec in metapaths.items():
        assert virus_entry_predicate(spec["path"], REAL_LAYERS["spine"]), code


def test_training_splits_are_not_counted_as_layers(tmp_path):
    """THE LATENT FAULT. v0.1-train2021 and friends are temporal training
    splits -- derived subsets of the assembled release -- so they carry
    INHIBITS and HOST_FACTOR_FOR exactly as the release does. Unfiltered, the
    dry run printed

        M1  enters via INHIBITS (chembl, train2021, train2021-corrected, ...)

    naming four things that supply nothing. The verdict held only because
    coverage is measured separately and no split name was ever in `covering`,
    so this was one step from being a wrong answer rather than a cosmetic one.
    """
    for name, pred in (("spine", "ENCODED_BY"), ("chembl", "INHIBITS"),
                       ("train2021", "INHIBITS"),
                       ("train2021-corrected-strict", "INHIBITS")):
        (tmp_path / f"v0.1-{name}").mkdir()
        (tmp_path / f"v0.1-{name}" / "edges.jsonl").write_text(
            json.dumps({"predicate": pred}) + "\n")

    got, skipped = layer_predicates(tmp_path, "v0.1",
                                    include={"spine", "chembl"})
    assert set(got) == {"spine", "chembl"}
    assert sorted(skipped) == ["train2021", "train2021-corrected-strict"]


def test_a_skipped_directory_is_reported_not_dropped(tmp_path):
    """A layer silently excluded is this project's own worst bug on record:
    ingest_folds writes a layer LAYER_PRECEDENCE does not list, and
    assemble_release filtered it out while every count still looked right. So
    the filter returns what it removed."""
    (tmp_path / "v0.1-folds").mkdir()
    (tmp_path / "v0.1-folds" / "edges.jsonl").write_text(
        json.dumps({"predicate": "FOLD_SIMILAR_TO"}) + "\n")
    got, skipped = layer_predicates(tmp_path, "v0.1", include={"spine"})
    assert got == {}
    assert skipped == ["folds"], \
        "an unlisted layer must be named, since that is exactly the M8 gap"
