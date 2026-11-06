"""MONDO disease labels, read from the ontology release on disk.

WHY THIS EXISTS. Every Disease node in v0.1 was labelled `f"{virus} disease"`
-- a string built from the virus being ingested -- and two of the ten MONDO
ids are shared by four viruses each, so the surviving label was whichever
virus the config listed last:

    MONDO:0005709   graph said "HCoV-HKU1 disease"   MONDO says "common cold"
    MONDO:0005502   graph said "DENV-4 disease"      MONDO says "dengue disease"
    MONDO:0100096   graph said "SARS-CoV-2 disease"  MONDO says "COVID-19"

Six of sixteen viruses pointed at a disease node named after a different
virus, and the other four labels were invented where a real one existed.

THE DATA WAS ALREADY HERE. data/raw/mondo/mondo.json has been on disk since
2026-10-28 and the 2026-11-07 audit found that no ingest reads it: 107 MB
of ontology carried along and never opened, while the labels it holds were
being made up one directory away. It parses in about a second.

ONLY THE LABELS, and only for ids that were asked for. The full release is
36,131 terms; this project pins ten. Returning the whole map would invite
callers to resolve ids the config never declared, which is how an
unvalidated ontology lookup becomes a silent source of nodes.
"""
from __future__ import annotations

import json
from pathlib import Path

#: Where fetch/record puts it. data/raw/refs/mondo.json is a second,
#: older copy; fix list P14 is about deleting one of them.
DEFAULT_PATH = Path("data/raw/mondo/mondo.json")

_IRI_MARK = "MONDO_"


def _curie(iri: str) -> str | None:
    """"http://purl.obolibrary.org/obo/MONDO_0005709" -> "MONDO:0005709"."""
    if _IRI_MARK not in iri:
        return None
    return "MONDO:" + iri.rsplit(_IRI_MARK, 1)[1]


def load_labels(path: Path, wanted: set[str]) -> dict[str, str]:
    """{MONDO curie: label} for the ids in `wanted` that the release defines.

    An id the release does not define is simply absent from the result.
    The caller decides what that means -- ingest_spine refuses rather than
    falling back to an invented label, because falling back is the bug
    this function exists to fix.
    """
    doc = json.loads(Path(path).read_text())
    out: dict[str, str] = {}
    for graph in doc.get("graphs", ()):
        for node in graph.get("nodes", ()):
            curie = _curie(node.get("id") or "")
            if curie in wanted and node.get("lbl"):
                out[curie] = node["lbl"]
    return out
