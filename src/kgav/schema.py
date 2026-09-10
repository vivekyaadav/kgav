"""Schema loading and validation for the antiviral repurposing KG.

The single job of this module: make it impossible for a triple that is not
declared in kg_schema.yaml to enter the graph, and impossible for an edge to
enter without provenance and a date.

Every validate_* function returns a list of Violation. Empty list == valid.
"""
from __future__ import annotations

import datetime as _dt
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

CURIE_RE = re.compile(r"^([A-Za-z][A-Za-z0-9_.]*):(\S+)$")
ANY = "any"


@dataclass(frozen=True)
class Violation:
    code: str
    message: str
    ref: str = ""

    def __str__(self) -> str:
        return f"[{self.code}] {self.ref}: {self.message}" if self.ref else f"[{self.code}] {self.message}"


@dataclass
class EdgeClass:
    predicate: str
    subject: str
    object: str
    subject_constraints: dict = field(default_factory=dict)
    object_constraints: dict = field(default_factory=dict)
    qualifiers: dict = field(default_factory=dict)
    required_qualifiers: list = field(default_factory=list)
    held_out: bool = False
    is_model_output: bool = False

    @property
    def key(self) -> tuple:
        return (self.subject, self.predicate, self.object)


def _is_iso_date(v: Any) -> bool:
    if not isinstance(v, str):
        return False
    try:
        _dt.date.fromisoformat(v[:10])
        return True
    except ValueError:
        return False


_TYPE_CHECKS = {
    "string": lambda v: isinstance(v, str),
    "int": lambda v: isinstance(v, int) and not isinstance(v, bool),
    "float": lambda v: isinstance(v, (int, float)) and not isinstance(v, bool),
    "bool": lambda v: isinstance(v, bool),
    "date": lambda v: isinstance(v, (_dt.date, _dt.datetime)) or _is_iso_date(v),
    "list[string]": lambda v: isinstance(v, (list, tuple)) and all(isinstance(x, str) for x in v),
}


class Schema:
    def __init__(self, doc: Mapping[str, Any]):
        self.raw = doc
        self.version: str = doc["schema_version"]
        self.prefixes: dict = doc["prefixes"]
        self.enums: dict = doc["enums"]
        self.node_classes: dict = doc["node_classes"]
        self.metapaths: dict = doc.get("metapaths", {})
        self.retrieval_limits: dict = doc.get("retrieval_limits", {})

        prov = doc["edge_provenance"]
        self.required_provenance: list = list(prov["required"])

        self.edge_classes: list[EdgeClass] = [
            EdgeClass(
                predicate=e["predicate"],
                subject=e["subject"],
                object=e["object"],
                subject_constraints=e.get("subject_constraints") or {},
                object_constraints=e.get("object_constraints") or {},
                qualifiers=e.get("qualifiers") or {},
                required_qualifiers=list(e.get("required_qualifiers") or []),
                held_out=bool(e.get("held_out", False)),
                is_model_output=bool(e.get("is_model_output", False)),
            )
            for e in doc["edge_classes"]
        ]
        self._by_key = {ec.key: ec for ec in self.edge_classes}
        self._by_predicate: dict[str, list[EdgeClass]] = {}
        for ec in self.edge_classes:
            self._by_predicate.setdefault(ec.predicate, []).append(ec)

    @classmethod
    def load(cls, path: str | Path) -> Schema:
        with open(path, encoding="utf-8") as fh:
            return cls(yaml.safe_load(fh))

    def held_out_predicates(self) -> set[str]:
        """Predicates that must be stripped from the TRAINING graph."""
        return {ec.predicate for ec in self.edge_classes if ec.held_out}

    def model_output_predicates(self) -> set[str]:
        return {ec.predicate for ec in self.edge_classes if ec.is_model_output}

    def find_edge_class(self, subj_class: str, predicate: str, obj_class: str) -> EdgeClass | None:
        for key in ((subj_class, predicate, obj_class),
                    (ANY, predicate, ANY),
                    (subj_class, predicate, ANY),
                    (ANY, predicate, obj_class)):
            if key in self._by_key:
                return self._by_key[key]
        return None

    def _check_value(self, name: str, spec: Mapping, value: Any, ref: str) -> list[Violation]:
        out: list[Violation] = []
        typ = spec.get("type", "string")
        check = _TYPE_CHECKS.get(typ)
        if check is None:
            return [Violation("SCHEMA_BAD_TYPE", f"unknown declared type {typ!r} for {name}", ref)]
        if not check(value):
            return [Violation("TYPE", f"{name}={value!r} is not {typ}", ref)]
        if "enum" in spec:
            allowed = self.enums.get(spec["enum"], [])
            if value not in allowed:
                out.append(Violation("ENUM", f"{name}={value!r} not in {spec['enum']} {allowed}", ref))
        if typ in ("int", "float"):
            if "min" in spec and value < spec["min"]:
                out.append(Violation("RANGE", f"{name}={value} < min {spec['min']}", ref))
            if "max" in spec and value > spec["max"]:
                out.append(Violation("RANGE", f"{name}={value} > max {spec['max']}", ref))
        return out

    def _curie_prefix(self, curie: str) -> str | None:
        m = CURIE_RE.match(curie or "")
        return m.group(1) if m else None

    def validate_node(self, node: Mapping[str, Any]) -> list[Violation]:
        out: list[Violation] = []
        nid = node.get("id", "")
        ref = nid or "<no id>"
        cls = node.get("class")

        if cls not in self.node_classes:
            return [Violation("UNKNOWN_NODE_CLASS", f"{cls!r} is not a declared node class", ref)]

        spec = self.node_classes[cls]
        prefix = self._curie_prefix(nid)
        if prefix is None:
            out.append(Violation("CURIE", f"id {nid!r} is not a CURIE", ref))
        elif prefix not in self.prefixes:
            out.append(Violation("PREFIX", f"prefix {prefix!r} is not declared", ref))
        elif prefix not in spec.get("id_prefixes", []):
            out.append(Violation("PREFIX",
                                 f"prefix {prefix!r} not permitted for {cls} "
                                 f"(allowed: {spec.get('id_prefixes')})", ref))

        props = node.get("properties") or {}
        declared = spec.get("properties") or {}
        for pname, pspec in declared.items():
            if pspec.get("required") and props.get(pname) is None:
                out.append(Violation("MISSING_PROPERTY", f"{cls}.{pname} is required", ref))
        for pname, value in props.items():
            if pname not in declared:
                out.append(Violation("UNDECLARED_PROPERTY", f"{pname!r} is not declared on {cls}", ref))
            elif value is not None:
                out.extend(self._check_value(pname, declared[pname], value, ref))
        return out

    def validate_edge(self,
                      edge: Mapping[str, Any],
                      node_index: Mapping[str, Mapping[str, Any]] | None = None) -> list[Violation]:
        out: list[Violation] = []
        subj, pred, obj = edge.get("subject"), edge.get("predicate"), edge.get("object")
        ref = f"{subj} -{pred}-> {obj}"

        if pred not in self._by_predicate:
            return [Violation("UNKNOWN_PREDICATE", f"{pred!r} is not a declared predicate", ref)]

        ec = None
        if node_index is not None:
            snode, onode = node_index.get(subj), node_index.get(obj)
            if snode is None:
                out.append(Violation("DANGLING", f"subject {subj!r} not in node index", ref))
            if onode is None:
                out.append(Violation("DANGLING", f"object {obj!r} not in node index", ref))
            if snode is None or onode is None:
                return out

            ec = self.find_edge_class(snode.get("class"), pred, onode.get("class"))
            if ec is None:
                allowed = [f"({c.subject}, {c.predicate}, {c.object})" for c in self._by_predicate[pred]]
                return [Violation(
                    "UNDECLARED_TRIPLE",
                    f"({snode.get('class')}, {pred}, {onode.get('class')}) is not a declared "
                    f"edge class; permitted: {allowed}", ref)]

            for endpoint, node, constraints in (("subject", snode, ec.subject_constraints),
                                                ("object", onode, ec.object_constraints)):
                for cname, cval in constraints.items():
                    actual = (node.get("properties") or {}).get(cname)
                    if actual != cval:
                        out.append(Violation(
                            "CONSTRAINT",
                            f"{endpoint} must have {cname}={cval!r} for {pred}, got {actual!r}", ref))
        else:
            candidates = self._by_predicate[pred]
            ec = candidates[0] if len(candidates) == 1 else None

        quals = edge.get("qualifiers") or {}
        if ec is not None:
            for qname, value in quals.items():
                if qname not in ec.qualifiers:
                    out.append(Violation("UNDECLARED_QUALIFIER",
                                         f"{qname!r} is not declared on {pred}", ref))
                elif value is not None:
                    out.extend(self._check_value(qname, ec.qualifiers[qname], value, ref))
            for qname in ec.required_qualifiers:
                if quals.get(qname) is None:
                    out.append(Violation("MISSING_QUALIFIER", f"{pred} requires {qname}", ref))

        for pname in self.required_provenance:
            if edge.get(pname) in (None, "", []):
                out.append(Violation("MISSING_PROVENANCE", f"{pname} is required on every edge", ref))

        tier = edge.get("evidence_tier")
        if tier is not None and tier not in self.enums["evidence_tier"]:
            out.append(Violation("ENUM", f"evidence_tier={tier!r} not in {self.enums['evidence_tier']}", ref))

        date = edge.get("first_asserted_date")
        if date is not None and not _is_iso_date(date) and not isinstance(date, (_dt.date, _dt.datetime)):
            out.append(Violation("TYPE", f"first_asserted_date={date!r} is not an ISO date", ref))

        if ec is not None and ec.is_model_output and tier in self.enums["evidence_tier"]:
            out.append(Violation("MODEL_OUTPUT_IN_EVIDENCE",
                                 f"{pred} is model output and must not carry an evidence tier", ref))
        return out

    def validate_batch(self,
                       nodes: Iterable[Mapping[str, Any]],
                       edges: Iterable[Mapping[str, Any]]) -> list[Violation]:
        nodes = list(nodes)
        out: list[Violation] = []
        index: dict[str, Mapping[str, Any]] = {}
        for n in nodes:
            out.extend(self.validate_node(n))
            nid = n.get("id")
            if nid in index:
                out.append(Violation("DUPLICATE_NODE", f"{nid} appears more than once", str(nid)))
            index[nid] = n
        for e in edges:
            out.extend(self.validate_edge(e, index))
        return out


def load_schema(path: str | Path | None = None) -> Schema:
    if path is None:
        path = Path(__file__).resolve().parents[2] / "schema" / "kg_schema.yaml"
    return Schema.load(path)
