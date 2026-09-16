"""Parameterised graph tools.

WHY NOT LET THE MODEL WRITE QUERIES. A language model asked to generate
traversal code against a 21-predicate schema fails quietly: a malformed
traversal returns an empty result, and an empty result is indistinguishable
from "no evidence exists". The model then reports the latter. Every tool here
is a fixed function with typed parameters, so a failure is a failure rather
than a confident absence.

WHAT EACH TOOL RETURNS. A structured record plus a verbalisation — a sentence
rendering with provenance inline. The model reads text, so leaving
verbalisation to the model means it invents phrasing for facts it half
understands. Rendering it here keeps the wording tied to the data.

BOUNDED BY CONSTRUCTION. Every tool caps its output. An unbounded neighbour
expansion from a well-connected protein returns thousands of rows and fills the
context window with ubiquitin.
"""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from kgav.agent import guards

MAX_NEIGHBOURS = 50
MAX_PATHS = 20
MAX_SUBGRAPH_EDGES = 100


@dataclass
class ToolResult:
    """Every tool returns this shape.

    `warnings` are not decoration: they carry the cell-context and
    host-directed caveats that make a path interpretable, and the verification
    stage checks they survived into the answer.
    """
    ok: bool
    kind: str
    data: Any = None
    verbalised: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    edge_ids: list[str] = field(default_factory=list)
    note: str = ""

    def to_json(self) -> str:
        return json.dumps({"ok": self.ok, "kind": self.kind, "data": self.data,
                           "text": self.verbalised, "warnings": self.warnings,
                           "edge_ids": self.edge_ids, "note": self.note},
                          indent=2, default=str)


class GraphTools:
    """The graph, loaded once, exposed as a fixed set of callable tools."""

    def __init__(self, release: Path, name_index: dict[str, str] | None = None):
        self.nodes: dict[str, dict] = {}
        self.edges: list[dict] = []
        self.out: dict[tuple[str, str], list[dict]] = defaultdict(list)
        self.inv: dict[tuple[str, str], list[dict]] = defaultdict(list)
        self.by_id: dict[str, dict] = {}
        # lowercase display name -> node id, for entity resolution
        self.names: dict[str, str] = dict(name_index or {})

        release = Path(release)

        # Compound names are not in the graph: nodes are keyed on InChIKey and
        # the ChEMBL layer stores only chembl_id, because path scoring never
        # needs a name. An agent has a human on the other end who types
        # "remdesivir", so the release carries a prebuilt name index rather
        # than the agent opening a 25 GB database at runtime.
        self.preferred_names: dict[str, str] = {}
        names_file = release / "NAMES.json"
        if names_file.exists():
            self.names.update(json.loads(names_file.read_text()))
        # Assembly merges the per-study selectivity edges into one, keeping the
        # worst value. That is right for filtering and wrong for display, so
        # the spread is read from the pre-assembly layer.
        self.si_range: dict[tuple[str, str], list[float]] = {}
        self.si_cell: dict[tuple[str, str], str] = {}
        si_layer = release.parent / (release.name + "-selectivity") / "edges.jsonl"
        if si_layer.exists():
            for line in si_layer.read_text().splitlines():
                if not line.strip():
                    continue
                e = json.loads(line)
                si = (e.get("qualifiers") or {}).get("selectivity_index")
                if si is not None:
                    self.si_range.setdefault((e["subject"], e["object"]), []).append(si)
                    cl = (e.get("qualifiers") or {}).get("cell_line")
                    if cl:
                        self.si_cell[(e["subject"], e["object"])] = cl

        pref_file = release / "PREF_NAMES.json"
        if pref_file.exists():
            self.preferred_names.update(json.loads(pref_file.read_text()))

        for line in (release / "nodes.jsonl").read_text().splitlines():
            if line.strip():
                n = json.loads(line)
                self.nodes[n["id"]] = n
                label = (n["properties"].get("gene_symbol")
                         or n["properties"].get("label"))
                if label:
                    self.names.setdefault(str(label).lower(), n["id"])

        for i, line in enumerate(
                (release / "edges.jsonl").read_text().splitlines()):
            if not line.strip():
                continue
            e = json.loads(line)
            e["_id"] = e.get("edge_id") or f"E{i}"
            self.edges.append(e)
            self.by_id[e["_id"]] = e
            self.out[(e["predicate"], e["subject"])].append(e)
            self.inv[(e["predicate"], e["object"])].append(e)

    # ------------------------------------------------------------- helpers
    def _label(self, node_id: str) -> str:
        """Display name, most specific first.

        protein_family precedes gene_symbol deliberately: every one of the 16
        mature peptides in the coronavirus replicase shares the gene symbol
        "rep", so preferring the symbol discards exactly the target
        specificity the polyprotein expansion existed to recover. A route
        reading "inhibits rep" is useless; "inhibits nsp5" is the finding.
        """
        n = self.nodes.get(node_id)
        if not n:
            return node_id
        p = n["properties"]
        if n["class"] == "SmallMolecule":
            # Reverse the name index: a human reads "nirmatrelvir", not an
            # accession.
            # Prefer the canonical name. A reverse lookup over the synonym
            # index returns whichever alias happened to be inserted last --
            # "pf07321332" for nirmatrelvir, "nsc-187208" for chloroquine.
            name = self.preferred_names.get(node_id)
            if name:
                return name
            # Fall back to a label on the node before the bare accession. A
            # release without PREF_NAMES.json should still render readable
            # names rather than "CHEMBL1".
            return str(p.get("label") or p.get("chembl_id") or node_id)
        if p.get("mature_peptide") and p.get("protein_family"):
            return f"{p['protein_family']} ({p['mature_peptide']})"
        return str(p.get("protein_family") or p.get("gene_symbol")
                   or p.get("label") or node_id)

    @staticmethod
    def _potency_note(nm: float) -> str:
        """Qualify a potency value in words.

        Without this, "IC50 75 nM" and "IC50 7280 nM" read as equivalent
        claims. The second is 7.3 uM -- above the threshold this project uses
        for "not pursued" -- and in the case that prompted this, it is
        chloroquine reported against a protease it does not inhibit, from a
        broad 2020 screening panel.
        """
        if nm <= 100:
            return ", potent"
        if nm <= 1000:
            return ", moderate"
        return ", WEAK — above 1 uM and probably not a real mechanism"

    def _verbalise(self, e: dict) -> str:
        """One edge as a sentence with its provenance inline.

        Templates live here rather than in the model, so a fact is always
        described the same way and always with its source attached.
        """
        s, o = self._label(e["subject"]), self._label(e["object"])
        q = e.get("qualifiers") or {}
        src = (e.get("primary_knowledge_source") or "?").replace("infores:", "")
        # Cap citations: one nirmatrelvir edge carries 51 PMIDs, which is
        # unreadable and would consume the context window.
        all_pubs = e.get("publications") or []
        pubs = ", ".join(all_pubs[:3])
        if len(all_pubs) > 3:
            pubs += f" +{len(all_pubs) - 3} more"
        support = e.get("support_count", 1)

        templates = {
            "INHIBITS": lambda: (
                f"{s} inhibits the viral protein {o}"
                + (f" (IC50 {q['ic50_nm']:.0f} nM{self._potency_note(q['ic50_nm'])})"
                   if q.get("ic50_nm") else "")
                + (f", assayed against the {q['domain']} site" if q.get("domain") else "")),
            "TARGETS": lambda: (
                f"{s} engages the host protein {o}"
                + (f" as an {q['direction']}" if q.get("direction") not in (None, "unknown")
                   else " (mechanism not stated: potency measured, not direction)")),
            "HOST_FACTOR_FOR": lambda: (
                f"{s} is a host {q.get('direction', 'factor')} factor for {o}"
                + (f", found in {support} independent CRISPR screens" if support > 1
                   else ", found in 1 CRISPR screen")
                + (f" ({q['cell_line']})" if q.get("cell_line") else "")),
            "PHYSICALLY_INTERACTS_WITH": lambda: (
                f"{s} and {o} were detected in the same complex or neighbourhood"),
            "ENCODED_BY": lambda: f"{s} is encoded by {o}",
            "BELONGS_TO": lambda: f"{o} carries the gene {s}",
            "PARTICIPATES_IN": lambda: f"{s} participates in {o}",
            "HAS_ANTIVIRAL_ACTIVITY_AGAINST": lambda: (
                f"{s} shows activity against {o}"
                + (f" (EC50 {q['ec50_nm']:.0f} nM)" if q.get("ec50_nm") else "")),
            "CHEMICALLY_SIMILAR_TO": lambda: (
                f"{s} is structurally similar to {o} "
                f"(Tanimoto {q.get('tanimoto_ecfp4', 0):.2f})"),
            "CAUSES": lambda: f"{s} causes {o}",
        }
        body = templates.get(e["predicate"], lambda: f"{s} {e['predicate']} {o}")()
        cite = f" [{src}"
        if pubs:
            cite += f", {pubs}"
        cite += f", {e.get('first_asserted_date', '?')[:4]}; {e['_id']}]"
        return body + cite

    # --------------------------------------------------------------- tools
    def resolve(self, text: str, limit: int = 5) -> ToolResult:
        """Free text to node identifiers.

        Ambiguity is returned, not resolved. Entity resolution is where a
        pipeline silently answers a question about the wrong molecule.
        """
        t = (text or "").strip().lower()
        if not t:
            return ToolResult(False, "resolve", note="empty query")
        exact = self.names.get(t)
        hits = [exact] if exact else []
        if not exact:
            hits = [nid for name, nid in self.names.items()
                    if t in name][:limit]
        if not hits:
            return ToolResult(False, "resolve",
                              note=f"'{text}' does not appear in this graph. It "
                                   "may be absent rather than misspelled: the "
                                   "graph covers compounds with published "
                                   "coronavirus or human-target activity.")
        return ToolResult(
            True, "resolve",
            data=[{"id": h, "label": self._label(h),
                   "class": self.nodes[h]["class"]} for h in hits],
            note="ambiguous — ask which was meant" if len(hits) > 1 else "")

    def explain(self, compound: str, virus: str,
                max_paths: int = MAX_PATHS) -> ToolResult:
        """Mechanistic paths from a compound to a virus, with provenance.

        The supported use of this graph. Given a compound with observed
        activity, this returns why it might act, as chains of individually
        sourced facts.
        """
        verdict = guards.check(guards.QueryKind.EXPLAIN, virus)
        paths = self._find_paths(compound, virus, max_paths)
        if not paths:
            return ToolResult(
                False, "explain", warnings=verdict.warnings,
                note=f"No mechanistic path connects {self._label(compound)} to "
                     f"{self._label(virus)} in this graph. That is an absence of "
                     "evidence, not evidence of absence: most compounds have an "
                     "activity measurement and no target annotation.")

        lines, warns, ids = [], list(verdict.warnings), []
        for i, (route, edges) in enumerate(paths, 1):
            lines.append(f"Route {i} ({route}):")
            for e in edges:
                lines.append("  " + self._verbalise(e))
                ids.append(e["_id"])
            nodes = [self.nodes.get(e["subject"]) or {} for e in edges]
            nodes += [self.nodes.get(e["object"]) or {} for e in edges]
            # Attach each warning to the route that triggered it. Collecting
            # them at the compound level produced a cell-context caveat naming
            # SIGMAR1 printed against a route through Spike glycoprotein --
            # a mismatch that discredits every other warning on the page.
            for w in guards.path_warnings(nodes, route):
                lines.append(f"  ! {w}")

        warns.append(self._selectivity_for(compound, virus))
        return ToolResult(True, "explain",
                          data=[{"route": r, "n_edges": len(es)} for r, es in paths],
                          verbalised=lines,
                          warnings=list(dict.fromkeys(w for w in warns if w)),
                          edge_ids=ids)

    def _best_direct_potency(self, compound: str) -> float | None:
        """Most potent INHIBITS measurement, in nM."""
        vals = [q[k]
                for e in self.out.get(("INHIBITS", compound), [])
                for q in [e.get("qualifiers") or {}]
                for k in ("ic50_nm", "ki_nm", "kd_nm")
                if q.get(k) is not None]
        return min(vals) if vals else None

    def _selectivity_for(self, compound: str, virus: str) -> str:
        """Selectivity as a RANGE, not a single value.

        The ingest keeps the worst selectivity across studies, which is right
        for filtering a positive set and wrong for display: reporting
        nirmatrelvir as "cytotoxic, SI 1.8" from one unflattering paper, with
        no sign that others disagree, states something false about a licensed
        drug. Studies disagree because cell lines and protocols differ, and
        the disagreement is information.
        """
        vals = list(self.si_range.get((compound, virus), []))
        if not vals:
            vals = [q["selectivity_index"]
                    for e in self.out.get(("HAS_ANTIVIRAL_ACTIVITY_AGAINST", compound), [])
                    if e["object"] == virus
                    and (q := e.get("qualifiers") or {}).get("selectivity_index") is not None]
        if not vals:
            return guards.selectivity_note(None, None)
        note = guards.selectivity_note(min(vals), True,
                                       n_studies=len(vals), si_max=max(vals))
        cl = guards.cell_line_note(self.si_cell.get((compound, virus)))
        return (note + " " + cl).strip() if cl else note

    def _find_paths(self, compound: str, virus: str, limit: int
                    ) -> list[tuple[str, list[dict]]]:
        """Enumerate metapath instances. Bounded and non-revisiting."""
        found: list[tuple[str, list[dict]]] = []

        # M1: compound -> viral protein -> gene -> virus
        for e1 in self.out.get(("INHIBITS", compound), []):
            for e2 in self.out.get(("ENCODED_BY", e1["object"]), []):
                for e3 in self.out.get(("BELONGS_TO", e2["object"]), []):
                    if e3["object"] == virus:
                        found.append(("M1 direct-acting", [e1, e2, e3]))
                        if len(found) >= limit:
                            return found

        # M2: compound -> host protein -> virus (dependency only)
        for e1 in self.out.get(("TARGETS", compound), []):
            for e2 in self.out.get(("HOST_FACTOR_FOR", e1["object"]), []):
                q = e2.get("qualifiers") or {}
                if e2["object"] == virus and q.get("direction") == "dependency":
                    found.append(("M2 host-directed", [e1, e2]))
                    if len(found) >= limit:
                        return found
        return found

    def triage(self, compounds: list[str], virus: str) -> ToolResult:
        """Rank compounds you already have.

        Supported because the recall ceiling does not apply: these compounds
        are in the graph by construction.
        """
        verdict = guards.check(guards.QueryKind.TRIAGE, virus)
        # Report anything that could not be scored. Silently returning three
        # of four requested compounds is a failure a user cannot detect.
        missing = [c for c in compounds if c not in self.nodes]
        scored = []
        for c in compounds:
            if c in missing:
                continue
            paths = self._find_paths(c, virus, MAX_PATHS)
            direct = sum(1 for r, _ in paths if r.startswith("M1"))
            host = len(paths) - direct
            sel = self._selectivity_for(c, virus)
            scored.append({"compound": c, "label": self._label(c),
                           "direct_acting_routes": direct,
                           "host_directed_routes": host,
                           "selectivity": sel})
        # Counting routes is wrong and inverted the controls: chloroquine
        # ranked ABOVE nirmatrelvir because it has a 160 nM Spike hit and a
        # 7,280 nM nsp5 hit -- two routes, one of them a screening artifact --
        # and then won the tiebreak on host-directed routes, which score BELOW
        # chance and should count for nothing.
        #
        # Rank on what the evaluation actually established: potent
        # direct-acting evidence, then selectivity. Host-directed routes are
        # reported for context and contribute zero.
        for r in scored:
            r["best_potency_nm"] = self._best_direct_potency(r["compound"])
            r["selectivity_index"] = self.si_range.get(
                (r["compound"], virus), [None])[0] \
                if self.si_range.get((r["compound"], virus)) else None
        def key(r):
            si = r["selectivity_index"]
            potent = (r["best_potency_nm"] is not None
                      and r["best_potency_nm"] <= 1000)
            # SELECTIVITY BEFORE POTENCY. Ordering on potency first ranked
            # chloroquine (160 nM, SI 15) above remdesivir (1,560 nM, SI 107)
            # -- and potency without selectivity is precisely the signal that
            # made chloroquine look promising in 2020. Selectivity is what
            # took M1's AUC from 0.729 to 0.826.
            #
            # Compounds with an unverified selectivity index rank below those
            # with one: unknown is not a passing grade.
            # NO POTENCY GATE. A 1 uM threshold has no basis in this
            # project's evaluation and put remdesivir (1,560 nM, SI 107) below
            # chloroquine (160 nM, SI 15) purely by sitting the wrong side of
            # an arbitrary line. The ranking may only use quantities the
            # evaluation measured: direct-acting evidence, then selectivity,
            # with potency as a tiebreak.
            del potent
            return (0 if r["direct_acting_routes"] else 1,
                    0 if si is not None else 1,
                    -(si or 0),
                    r["best_potency_nm"] if r["best_potency_nm"] is not None else 1e12)
        scored.sort(key=key)
        lines = []
        for i, r in enumerate(scored, 1):
            pot = (f"best direct IC50 {r['best_potency_nm']:.0f} nM"
                   if r["best_potency_nm"] is not None else "no direct-acting data")
            si = (f", SI {r['selectivity_index']:.0f}"
                  if r["selectivity_index"] is not None else ", SI unknown")
            lines.append(f"{i}. {r['label']}: {pot}{si} "
                         f"({r['direct_acting_routes']} direct-acting, "
                         f"{r['host_directed_routes']} host-directed routes)")
        warns = list(verdict.warnings)
        warns.append(
            "Ranking is by direct-acting route count. Direct-acting paths "
            "separate selective antivirals from measured inactives at AUC 0.826; "
            "host-directed paths score below chance and are shown for context "
            "only.")
        if missing:
            warns.append(f"{len(missing)} requested compound(s) are not in this "
                         f"graph and were not ranked: {', '.join(missing)}")
        return ToolResult(True, "triage", data=scored, verbalised=lines,
                          warnings=warns, note=f"{len(missing)} unranked" if missing else "")

    def evidence(self, subject: str, object_: str | None = None,
                 predicate: str | None = None, limit: int = 25) -> ToolResult:
        """Every edge supporting a specific claim, with provenance."""
        hits = [e for e in self.edges
                if e["subject"] == subject
                and (object_ is None or e["object"] == object_)
                and (predicate is None or e["predicate"] == predicate)][:limit]
        if not hits:
            return ToolResult(False, "evidence",
                              note="No edge in this graph makes that claim.")
        lines = [self._verbalise(e) for e in hits]
        warns = [guards.evidence_note(e.get("evidence_tier"),
                                      (e.get("qualifiers") or {}).get("detection_method"))
                 for e in hits]
        return ToolResult(True, "evidence",
                          data=[{k: v for k, v in e.items() if not k.startswith("_")}
                                for e in hits],
                          verbalised=lines,
                          warnings=list(dict.fromkeys(warns)),
                          edge_ids=[e["_id"] for e in hits])

    def neighbours(self, node_id: str, predicate: str | None = None,
                   limit: int = MAX_NEIGHBOURS) -> ToolResult:
        """Adjacent nodes, capped."""
        if node_id not in self.nodes:
            return ToolResult(False, "neighbours", note=f"{node_id} not in graph")
        hits = [e for e in self.edges
                if (e["subject"] == node_id or e["object"] == node_id)
                and (predicate is None or e["predicate"] == predicate)]
        truncated = len(hits) > limit
        hits = hits[:limit]
        return ToolResult(
            True, "neighbours",
            data=[{"edge": e["predicate"],
                   "other": e["object"] if e["subject"] == node_id else e["subject"]}
                  for e in hits],
            verbalised=[self._verbalise(e) for e in hits],
            edge_ids=[e["_id"] for e in hits],
            note=f"truncated to {limit}" if truncated else "")

    def profile(self, node_id: str) -> ToolResult:
        """What the graph holds about one entity."""
        n = self.nodes.get(node_id)
        if not n:
            return ToolResult(False, "profile", note=f"{node_id} not in graph")
        counts: dict[str, int] = defaultdict(int)
        for e in self.edges:
            if e["subject"] == node_id or e["object"] == node_id:
                counts[e["predicate"]] += 1
        lines = [f"{self._label(node_id)} ({n['class']}, {node_id})"]
        lines += [f"  {p}: {c}" for p, c in sorted(counts.items(),
                                                   key=lambda kv: -kv[1])]
        return ToolResult(True, "profile",
                          data={"id": node_id, "class": n["class"],
                                "properties": n["properties"],
                                "edge_counts": dict(counts)},
                          verbalised=lines)
