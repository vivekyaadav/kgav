"""Shared node/edge accumulator used by every ingest script."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Any


@dataclass
class Emit:
    """Accumulates nodes and edges, deduping node ids."""
    nodes: dict[str, dict] = field(default_factory=dict)
    edges: list[dict] = field(default_factory=list)
    notes: Counter = field(default_factory=Counter)

    def node(self, nid: str, cls: str, **props: Any) -> str:
        props = {k: v for k, v in props.items() if v is not None}
        if nid in self.nodes:
            self.nodes[nid]["properties"].update(props)
            self.notes[f"node_merged:{cls}"] += 1
        else:
            self.nodes[nid] = {"id": nid, "class": cls, "properties": props}
        return nid

    def edge(self, subj: str, pred: str, obj: str, *, source: str, date: str,
             tier: int = 1, quals: dict | None = None, pmids: list | None = None,
             score: float | None = None, support: int | None = None) -> None:
        e = {
            "subject": subj, "predicate": pred, "object": obj,
            "qualifiers": quals or {},
            "primary_knowledge_source": source,
            "evidence_tier": tier,
            "first_asserted_date": date,
        }
        if pmids:
            e["publications"] = pmids
        if support is not None:
            # HOW MANY WITNESSES BEHIND ONE FACT, for a source that cites
            # several and an identity that merges them. The indications
            # layer cites up to two ClinicalTrials registrations for one
            # (compound, disease) pair and nct_id is not an identity
            # qualifier, so without this the second registration would
            # vanish without trace rather than being counted.
            e["support_count"] = support
        if score is not None:
            # Raw score as the source reports it. Deliberately NOT called
            # confidence: a VirHostNet miscore and a STRING combined score are
            # not commensurable and neither is a probability.
            e["source_score"] = score
        self.edges.append(e)
