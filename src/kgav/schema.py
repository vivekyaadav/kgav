"""Schema loading and validation for the antiviral repurposing KG.

The single job of this module: make it impossible for a triple that is not
declared in kg_schema.yaml to enter the graph, and impossible for an edge to
enter without provenance and a date.

Every validate_* function returns a list of Violation. Empty list == valid.
Nothing here raises on invalid data; callers decide whether to reject the row
or fail the batch.
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

# A metapath hop prefixed with this is traversed object-to-subject.
REVERSE = "<"

# Qualifier-spec key declaring that a bool is DERIVED from the absence of
# other qualifiers, rather than ingested. One definition, in the schema.
DERIVED_FROM_ABSENCE = "derived_from_absence_of"


def parse_hop(token: str) -> tuple[str, bool]:
    """A metapath hop token -> (predicate, traversed_in_reverse).

    "<PARTICIPATES_IN" means walk the predicate backwards, Pathway to Protein.
    Direction is stated in the metapath rather than inferred from whether a
    node happens to have outgoing edges -- that inference is what silently
    truncated every symmetric predicate. See the metapaths block in the schema.
    """
    if token.startswith(REVERSE):
        return token[len(REVERSE):], True
    return token, False


@dataclass(frozen=True)
class Violation:
    code: str
    message: str
    ref: str = ""

    def __str__(self) -> str:  # pragma: no cover - display only
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
    # A label is a measurement an evaluation partitions into positives and
    # negatives. DISTINCT FROM held_out, which is about traversal: see the
    # edge_classes header in kg_schema.yaml. Every held-out predicate was also
    # a label until MEASURED_INACTIVE_AGAINST, which is held out and is not.
    evaluation_label: bool = False
    is_model_output: bool = False
    identity_qualifiers: list = field(default_factory=list)
    # When true, two edges from DIFFERENT knowledge sources are different
    # facts even with identical subject, predicate, object and qualifiers.
    # Needed for measurements: a Vero E6 panel reading and a ChEMBL-sourced
    # EC50 are two experiments, and collapsing them lets layer precedence
    # silently discard one. Qualifier-based identity cannot express this,
    # because the selectivity layer DEPENDS on sharing the bare triple with
    # the chemistry layer to supply its CC50 -- both are infores:chembl, so
    # keying on the source keeps that merge and splits only the panel.
    identity_includes_source: bool = False
    # None means UNDECLARED, which is different from False. A traversal may
    # not guess, so validate_metapaths rejects any metapath hop over a
    # predicate that left this unset.
    symmetric: bool | None = None

    @property
    def key(self) -> tuple:
        return (self.subject, self.predicate, self.object)


_TYPE_CHECKS = {
    "string": lambda v: isinstance(v, str),
    "int": lambda v: isinstance(v, int) and not isinstance(v, bool),
    "float": lambda v: isinstance(v, (int, float)) and not isinstance(v, bool),
    "bool": lambda v: isinstance(v, bool),
    "date": lambda v: isinstance(v, (_dt.date, _dt.datetime)) or _is_iso_date(v),
    "list[string]": lambda v: isinstance(v, (list, tuple)) and all(isinstance(x, str) for x in v),
}


def _is_iso_date(v: Any) -> bool:
    if not isinstance(v, str):
        return False
    try:
        _dt.date.fromisoformat(v[:10])
        return True
    except ValueError:
        return False


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
                evaluation_label=bool(e.get("evaluation_label", False)),
                is_model_output=bool(e.get("is_model_output", False)),
                identity_qualifiers=list(e.get("identity_qualifiers") or []),
                identity_includes_source=bool(
                    e.get("identity_includes_source", False)),
                symmetric=e.get("symmetric"),
            )
            for e in doc["edge_classes"]
        ]
        self._by_key = {ec.key: ec for ec in self.edge_classes}
        self._by_predicate: dict[str, list[EdgeClass]] = {}
        for ec in self.edge_classes:
            self._by_predicate.setdefault(ec.predicate, []).append(ec)

    # -- loading ---------------------------------------------------------
    @classmethod
    def load(cls, path: str | Path) -> Schema:
        with open(path, "r", encoding="utf-8") as fh:
            return cls(yaml.safe_load(fh))

    # -- lookups ---------------------------------------------------------
    def identity_qualifiers(self, predicate: str) -> list[str]:
        """Qualifiers that make two edges DIFFERENT FACTS rather than two
        witnesses to one. The assembler includes them in the edge key."""
        for ec in self._by_predicate.get(predicate, []):
            if ec.identity_qualifiers:
                return ec.identity_qualifiers
        return []

    def source_identity_predicates(self) -> set[str]:
        """Predicates whose facts are identified partly by their source."""
        return {ec.predicate for ec in self.edge_classes
                if ec.identity_includes_source}

    def symmetry(self) -> dict[str, bool]:
        """predicate -> whether the relation holds in both directions.

        Predicates that do not declare it are ABSENT from this mapping, not
        defaulted. A traversal that meets an absent predicate must fail loudly
        rather than fall back to guessing from local topology.
        """
        return {ec.predicate: bool(ec.symmetric)
                for ec in self.edge_classes if ec.symmetric is not None}

    def is_symmetric(self, predicate: str) -> bool | None:
        return self.symmetry().get(predicate)

    def derived_absence_qualifiers(self, predicate: str) -> dict[str, tuple[str, ...]]:
        """qualifier -> the qualifiers whose ABSENCE makes it true.

        A derived qualifier is a statement about the rest of the edge, so it
        is computed from the edge and never ingested. `unquantified` was the
        opposite of that: three modules each asserted it from their own local
        view, and the assembler's precedence picked a winner. See the schema
        note on INHIBITS.unquantified.

        A source qualifier that is not declared on the same edge class raises.
        Silently ignoring it would make the derivation read a field that can
        never be populated, so the qualifier would be true for every edge and
        look like data.
        """
        out: dict[str, tuple[str, ...]] = {}
        for ec in self._by_predicate.get(predicate, []):
            for qname, spec in ec.qualifiers.items():
                if not isinstance(spec, Mapping):
                    continue
                sources = spec.get(DERIVED_FROM_ABSENCE)
                if not sources:
                    continue
                undeclared = [s for s in sources if s not in ec.qualifiers]
                if undeclared:
                    raise ValueError(
                        f"{predicate}.{qname} derives from {undeclared}, which "
                        f"{'is' if len(undeclared) == 1 else 'are'} not declared "
                        f"on {predicate}")
                out[qname] = tuple(sources)
        return out

    def derive_qualifiers(self, predicate: str, quals: Mapping[str, Any]) -> dict:
        """`quals` with every derived qualifier recomputed from the evidence.

        Called by the assembler AFTER merging, because the layer holding the
        EC50 and the layer holding the CC50 are different layers and neither
        can see the other. Returns a new dict; the input is not mutated.
        """
        out = dict(quals)
        for qname, sources in self.derived_absence_qualifiers(predicate).items():
            out[qname] = not any(out.get(s) is not None for s in sources)
        return out

    def validate_metapaths(self) -> list[Violation]:
        """Every hop must name a declared predicate with a declared symmetry.

        This is the check that stops H1 regressing. Without it a new metapath
        over an undeclared predicate would silently resume the old behaviour
        of inferring direction from whether a node had outgoing edges.
        """
        out: list[Violation] = []
        sym = self.symmetry()
        for name, mp in self.metapaths.items():
            path = list(mp.get("path") or [])
            if len(path) < 3 or len(path) % 2 == 0:
                out.append(Violation(
                    "METAPATH_MALFORMED",
                    f"path must alternate class, predicate, class; got {len(path)} entries",
                    name))
                continue
            for i in range(1, len(path) - 1, 2):
                predicate, reverse = parse_hop(path[i])
                if predicate not in self._by_predicate:
                    out.append(Violation(
                        "METAPATH_UNKNOWN_PREDICATE",
                        f"hop {i // 2 + 1} traverses {predicate!r}, which is not a "
                        f"declared predicate", name))
                    continue
                if predicate not in sym:
                    out.append(Violation(
                        "METAPATH_UNDECLARED_SYMMETRY",
                        f"hop {i // 2 + 1} traverses {predicate!r}, which does not "
                        f"declare `symmetric`; traversal direction would have to be "
                        f"guessed from local topology", name))
                elif reverse and sym[predicate]:
                    out.append(Violation(
                        "METAPATH_REDUNDANT_REVERSE",
                        f"hop {i // 2 + 1} marks {predicate!r} reverse, but it is "
                        f"symmetric and is already traversed both ways; the marker "
                        f"signals a misunderstanding of the relation", name))
        return out

    def held_out_predicates(self) -> set[str]:
        """Predicates that must be stripped from the TRAINING graph."""
        return {ec.predicate for ec in self.edge_classes if ec.held_out}

    def evaluation_label_predicates(self) -> set[str]:
        """Predicates whose edges ARE labels, keyed by their object.

        Not the same set as held_out_predicates, and the difference is the
        point. A split that partitions every held-out edge into positives and
        negatives buckets MEASURED_INACTIVE_AGAINST by viral protein, putting
        protein nodes where viruses belong. Every predicate here terminates on
        an OrganismTaxon or a Disease; validate_labels() enforces that, so a
        predicate cannot be declared a label and then key the label sets on
        something that is not a label's subject.
        """
        return {ec.predicate for ec in self.edge_classes if ec.evaluation_label}

    def validate_labels(self) -> list[Violation]:
        """A label must terminate on a virus or a disease, and be held out.

        The first rule stops the bug this flag exists for from reappearing by
        a different route: a label keyed on a protein silently fills the label
        sets with protein nodes. The second is a consistency check -- a label
        left in the training graph is the leak the whole temporal protocol is
        built to prevent, so declaring one without holding it out is an error
        rather than a choice.
        """
        out: list[Violation] = []
        for ec in self.edge_classes:
            if not ec.evaluation_label:
                continue
            if ec.object not in ("OrganismTaxon", "Disease"):
                out.append(Violation(
                    "LABEL_BAD_OBJECT",
                    f"{ec.predicate} is declared evaluation_label but terminates "
                    f"on {ec.object}; a label is a measurement about a (compound, "
                    f"virus) or (compound, disease) pair, and keying label sets "
                    f"on anything else fills them with the wrong node class",
                    ec.predicate))
            if not ec.held_out:
                out.append(Violation(
                    "LABEL_NOT_HELD_OUT",
                    f"{ec.predicate} is a label but is not held_out, so it stays "
                    f"in the training graph and can be traversed to predict "
                    f"itself", ec.predicate))
        return out

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

    # -- internal helpers ------------------------------------------------
    def _check_value(self, name: str, spec: Mapping, value: Any, ref: str) -> list[Violation]:
        out: list[Violation] = []
        typ = spec.get("type", "string")
        check = _TYPE_CHECKS.get(typ)
        if check is None:
            out.append(Violation("SCHEMA_BAD_TYPE", f"unknown declared type {typ!r} for {name}", ref))
            return out
        if not check(value):
            out.append(Violation("TYPE", f"{name}={value!r} is not {typ}", ref))
            return out
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

    # -- node validation -------------------------------------------------
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
        else:
            if prefix not in self.prefixes:
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
                continue
            if value is not None:
                out.extend(self._check_value(pname, declared[pname], value, ref))
        return out

    # -- edge validation -------------------------------------------------
    def validate_edge(self,
                      edge: Mapping[str, Any],
                      node_index: Mapping[str, Mapping[str, Any]] | None = None) -> list[Violation]:
        """Validate one edge.

        node_index maps node id -> node dict. When supplied, endpoint classes
        and endpoint property constraints are enforced; without it only the
        predicate, qualifiers and provenance can be checked.
        """
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

        # qualifiers
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

            # A derived qualifier may be absent -- a single layer cannot
            # compute one -- but it may never DISAGREE with the evidence on
            # its own edge. 12,756 edges claimed unquantified=true while
            # carrying a measured value, 1,264 of them alongside a
            # selectivity index, and nothing rejected them.
            for qname, sources in self.derived_absence_qualifiers(pred).items():
                if quals.get(qname) is None:
                    continue
                expected = not any(quals.get(s) is not None for s in sources)
                if bool(quals[qname]) is not expected:
                    present = [s for s in sources if quals.get(s) is not None]
                    out.append(Violation(
                        "DERIVED_QUALIFIER",
                        f"{qname}={quals[qname]!r} contradicts the edge: "
                        f"{'carries ' + ', '.join(present) if present else 'carries no measured value'}"
                        f", so it must be {expected}", ref))

        # A SYMMETRIC EDGE FROM A NODE TO ITSELF ASSERTS NOTHING, and it
        # validates perfectly: the triple is declared, the qualifiers are
        # legal, the provenance is present. ingest_folds wrote 36 of them --
        # FOLD_SIMILAR_TO from a shared class node to itself -- and M8 produced
        # zero paths for every run since, because walk_metapath discards a
        # neighbour already in `visited`. Nothing rejected them and nothing
        # reported the consequence.
        #
        # Scoped to symmetric predicates. A directed self-reference can be
        # meaningful (a protein interacting with itself is a homodimer); "X is
        # similar to X" is not.
        if subj is not None and subj == obj and ec is not None and ec.symmetric:
            out.append(Violation(
                "SYMMETRIC_SELF_LOOP",
                f"{pred} from a node to itself asserts nothing: a symmetric "
                f"relation between one node and itself is not a fact, and a "
                f"metapath hop over it is discarded as already-visited", ref))

        # provenance — absolute
        for pname in self.required_provenance:
            if edge.get(pname) in (None, "", []):
                out.append(Violation("MISSING_PROVENANCE", f"{pname} is required on every edge", ref))

        tier = edge.get("evidence_tier")
        if tier is not None and tier not in self.enums["evidence_tier"]:
            out.append(Violation("ENUM", f"evidence_tier={tier!r} not in {self.enums['evidence_tier']}", ref))

        date = edge.get("first_asserted_date")
        if date is not None and not _is_iso_date(date) and not isinstance(date, (_dt.date, _dt.datetime)):
            out.append(Violation("TYPE", f"first_asserted_date={date!r} is not an ISO date", ref))

        # model output must never be written into an evidence tier
        if ec is not None and ec.is_model_output and tier in self.enums["evidence_tier"]:
            out.append(Violation("MODEL_OUTPUT_IN_EVIDENCE",
                                 f"{pred} is model output and must not carry an evidence tier", ref))
        return out

    # -- batch -----------------------------------------------------------
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
