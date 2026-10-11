"""Which metapath channels a missing input actually blocks.

WHY THIS IS NOT A LIST SOMEONE MAINTAINS

Before this, "what is still needed to build a flavivirus" was four or five
items in a comment, and nothing said which of them blocked WHAT. That reads as
one undifferentiated wall, and it hid a fact worth knowing: ORCS_KEY and the
VirHostNet downloads block the host-mediated channels, and M1 -- the
direct-acting channel, the only one this project has a defensible result on --
needs neither. A flavivirus M1 evaluation needs the proteome and nothing else.

So the mapping is DERIVED, twice over:

  metapath -> predicates   from schema/kg_schema.yaml, the same declaration
                           walk_metapath traverses.
  predicate -> layer       by reading the layer releases on disk and seeing
                           which predicates each one actually emitted.

Neither half is a hand-kept list, because this repository has watched four
hand-built taxon maps drift out of step with the register and replaced them
with declared_taxa for exactly this reason. A list of requirements that goes
stale is worse than no list: it is the stale text in dryrun_promote that still
said "there is no fetcher in this repo" after one was added.
"""
from __future__ import annotations

import json
from pathlib import Path

from kgav.schema import parse_hop


def layer_predicates(release_root: Path, prefix: str,
                     include: set[str] | None = None
                     ) -> tuple[dict[str, set[str]], list[str]]:
    """({layer name: predicates it emitted}, [directories skipped])

    Measured rather than declared: a layer's edges.jsonl is the only thing
    that knows what it produced. A layer that exists and is empty maps to an
    empty set, which is the honest answer and is what makes a missing input
    visible downstream.

    `include` RESTRICTS WHAT COUNTS AS A LAYER, and the first version had no
    such filter. It globbed <prefix>-* and so counted v0.1-train2021,
    v0.1-train2021-corrected and the two strict variants as layers -- those
    are TEMPORAL TRAINING SPLITS, derived subsets of the assembled release,
    so they contain INHIBITS and HOST_FACTOR_FOR like the release does. The
    dry run then printed

        M1  enters via INHIBITS (chembl, train2021, train2021-corrected, ...)

    naming four things that supply nothing. The verdict survived only because
    coverage is measured separately, so a split name was never in `covering` --
    which means this was a latent fault rather than a wrong answer, and it
    would have become a wrong answer the moment a split name got in.

    Skipped directories are RETURNED rather than dropped, because a layer
    silently excluded is the failure this project already has on record:
    ingest_folds writes a layer LAYER_PRECEDENCE does not list, and
    assemble_release filtered it out while every count still looked right.
    """
    out: dict[str, set[str]] = {}
    skipped: list[str] = []
    for d in sorted(release_root.glob(f"{prefix}-*")):
        if not d.is_dir():
            continue
        name = d.name[len(prefix) + 1:]
        if include is not None and name not in include:
            skipped.append(name)
            continue
        preds: set[str] = set()
        edges = d / "edges.jsonl"
        if edges.exists():
            for line in edges.read_text().splitlines():
                if line.strip():
                    pred = json.loads(line).get("predicate")
                    if pred:
                        preds.add(pred)
        out[name] = preds
    return out, skipped


def predicates_in(path: list[str]) -> list[str]:
    """The predicates in a metapath's declared path, direction stripped.

    A path alternates node class and predicate: [SmallMolecule, INHIBITS,
    Protein, ENCODED_BY, Gene, ...], so the predicates are the odd positions.

    schema.parse_hop does the stripping rather than a local startswith("<"),
    because the reverse marker is the schema's convention and a second copy of
    it is a second thing to keep in step. The first version of this module did
    reimplement it, forgot, and reported M4's "<PARTICIPATES_IN" as a
    predicate no layer emits -- while PARTICIPATES_IN sits in the host layer
    with 157,211 edges.
    """
    return [parse_hop(p)[0] for i, p in enumerate(path) if i % 2 == 1]


def channel_layers(metapaths: dict, layer_preds: dict[str, set[str]]
                   ) -> dict[str, dict[str, set[str]]]:
    """{metapath code: {"layers": needed, "unsupplied": predicates no layer has}}

    `unsupplied` is the part that matters and the reason this returns both. A
    predicate no layer emits means the channel cannot be walked at all, and
    saying only "needs layers X, Y" would hide that: M8 needed `folds`, a
    layer LAYER_PRECEDENCE does not list, and the channel read as merely empty
    for months rather than as unreachable.
    """
    out: dict[str, dict[str, set[str]]] = {}
    for code, spec in (metapaths or {}).items():
        needed: set[str] = set()
        unsupplied: set[str] = set()
        for pred in predicates_in(spec.get("path") or []):
            owners = {lyr for lyr, preds in layer_preds.items()
                      if pred in preds}
            if owners:
                needed |= owners
            else:
                unsupplied.add(pred)
        out[code] = {"layers": needed, "unsupplied": unsupplied}
    return out


def blocked_by(metapaths: dict, layer_preds: dict[str, set[str]],
               unavailable: set[str]) -> dict[str, set[str]]:
    """{metapath code: what blocks it}, empty set if it can be evaluated today.

    Two kinds of blocker, and the first version reported only one. It returned
    `info["layers"] & unavailable`, so M8 -- whose MEMBER_OF_CLASS and
    FOLD_SIMILAR_TO are emitted by NO layer, because ingest_folds writes a
    layer LAYER_PRECEDENCE does not list -- came back as an empty set and read
    as "evaluable today". A channel that cannot be walked at all reported as
    ready is the exact false green this module exists to prevent, produced by
    the module itself.

    An unsupplied predicate is returned as "no layer emits <PREDICATE>" so the
    two causes stay distinguishable: an unavailable layer is something to go
    and fetch, an unsupplied predicate is something to go and wire up.
    """
    out: dict[str, set[str]] = {}
    for code, info in channel_layers(metapaths, layer_preds).items():
        blockers = set(info["layers"] & unavailable)
        blockers |= {f"no layer emits {p}" for p in info["unsupplied"]}
        out[code] = blockers
    return out


def virus_entry_predicate(path: list[str],
                          spine_predicates: set[str]) -> str | None:
    """The predicate that must attach something to THIS virus's own proteins.

    WHY "IS THE PREDICATE WIRED" IS THE WRONG QUESTION FOR A NEW VIRUS.

    channel_layers answers whether some layer anywhere emits a predicate. For
    a virus being promoted that is not enough, and the first version of this
    module reported M3 as "evaluable once the proteome lands" because
    PHYSICALLY_INTERACTS_WITH is emitted by the host layer -- 150,432 edges of
    it. But M3 is

        SmallMolecule -TARGETS-> Protein -PHYSICALLY_INTERACTS_WITH->
        Protein -ENCODED_BY-> Gene -BELONGS_TO-> OrganismTaxon

    and that second Protein is the VIRAL one, because the Gene below it
    belongs to the virus. host's PPIs are human-human and can never reach a
    viral taxon; only the vhppi layer bridges host to virus, and it has no
    flavivirus files. So M3 cannot reach Zika, and reporting it as evaluable
    was a false green about the very thing this module exists to report.

    Derived by walking the path BACKWARDS from the OrganismTaxon endpoint.
    Every step the SPINE supplies (ENCODED_BY, BELONGS_TO) stays inside the
    virus's own subgraph, so the node before it is still viral. The first
    predicate the spine does not supply is the bridge: the one hop that has to
    come from a source carrying data about this particular virus.

    Returns None when the whole path is spine-supplied, which no declared
    metapath is -- every one of them has to enter the virus from somewhere.
    """
    preds = predicates_in(path)
    for pred in reversed(preds):
        if pred not in spine_predicates:
            return pred
    return None


def reachable_for_virus(metapaths: dict, layer_preds: dict[str, set[str]],
                        covering: set[str]) -> dict[str, dict[str, object]]:
    """{code: {"entry": predicate, "sources": layers, "reachable": bool}}

    `covering` is the set of layers that will hold data about the virus being
    promoted -- spine once its proteome is cached, chembl when the ChEMBL
    database has rows for its taxa, vhppi when its PPI files are downloaded,
    orcs when its screens are ingestable. A channel is reachable for that
    virus when at least one layer supplying its entry predicate is in there.

    This is a weaker claim than "the channel will produce paths": it says the
    hop into the virus has a source. Whether that source has enough rows is
    what diagnose_flavivirus_targets measures.
    """
    spine_preds = layer_preds.get("spine", set())
    out: dict[str, dict[str, object]] = {}
    for code, spec in (metapaths or {}).items():
        entry = virus_entry_predicate(spec.get("path") or [], spine_preds)
        sources = {lyr for lyr, preds in layer_preds.items()
                   if entry and entry in preds}
        out[code] = {"entry": entry, "sources": sources,
                     "reachable": bool(sources & covering)}
    return out
