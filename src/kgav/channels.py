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


def layer_predicates(release_root: Path, prefix: str) -> dict[str, set[str]]:
    """{layer name: predicates it emitted}, read from the layer releases.

    Measured rather than declared: a layer's edges.jsonl is the only thing
    that knows what it produced. A layer directory that does not exist, or
    exists and is empty, maps to an empty set -- which is the honest answer
    and is what makes a missing input visible downstream.
    """
    out: dict[str, set[str]] = {}
    for d in sorted(release_root.glob(f"{prefix}-*")):
        if not d.is_dir():
            continue
        name = d.name[len(prefix) + 1:]
        preds: set[str] = set()
        edges = d / "edges.jsonl"
        if edges.exists():
            for line in edges.read_text().splitlines():
                if line.strip():
                    p = json.loads(line).get("predicate")
                    if p:
                        preds.add(p)
        out[name] = preds
    return out


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
