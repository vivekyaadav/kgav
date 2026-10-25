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
from kgav.baselines import qualifier_ok
from kgav.schema import load_schema

MAX_NEIGHBOURS = 50
MAX_PATHS = 20
MAX_SUBGRAPH_EDGES = 100
# Pathways LISTED by pathways(), not pathways considered. Every one found is
# counted and the overflow is reported with how many of it carry a dependency
# factor, so the cap bounds the output without hiding the census.
MAX_PATHWAYS = 25


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

    def __init__(self, release: Path, name_index: dict[str, str] | None = None,
                 calibration=None, landscape=None):
        # ONE LANDSCAPE FOR THE WHOLE RELEASE, queried per virus.
        #
        # Unlike calibration there is no per-virus mapping to mis-key: the
        # landscape is a single measurement over every taxon and
        # Landscape.for_virus returns that taxon's row or None, never a
        # neighbour's. None here is not neutral -- guards.pathway_caveat
        # renders it as "not measured on this release", which is the honest
        # branch and still appears in the answer. The one thing an answer may
        # not do is stay silent about it.
        self.landscape = landscape
        # ONE CALIBRATION PER VIRUS, NOT ONE PER INSTANCE.
        #
        # This was a single `self.calibration`, handed to every caveat
        # regardless of which virus was asked about. explain() and triage()
        # both TAKE a virus argument and neither used it to pick the
        # calibration, so an agent constructed with SARS-CoV-2's numbers and
        # asked about Zika would state Zika's direct-acting route "separates
        # measured actives from measured inactives at AUC 0.675" -- SARS-CoV-2's
        # figure, on the discriminates branch, which is the most confident
        # sentence this module can emit. Zika's own M1 is 0.511, at chance.
        #
        # Latent while every virus was at chance anyway. Live the moment a
        # second family arrived with asymmetric calibration, and the mismatch
        # was always detectable: Calibration carries its own `.virus`, and
        # nothing compared it to the virus being answered about.
        #
        # Accepts a single Calibration (for one virus) or a {virus:
        # Calibration} mapping. A virus with no calibration resolves to None,
        # which routes the caveats to their honest "uncalibrated" branch --
        # never to another virus's numbers.
        self._calibrations: dict[str, object] = {}
        self._pending_calibration = calibration
        # METAPATH CONSTRAINTS, READ FROM THE SCHEMA RATHER THAN RESTATED.
        #
        # _find_paths enumerates M1 and M2 as its own loops. M1 applied NO
        # constraints and M2 a hand-copied one, so a constraint added to the
        # schema reached dwpc_scores -- which passes mp.get("constraints") --
        # and never reached here. When M1 gained
        # INHIBITS.target_resolution=chain_from_assay_description, the scorer
        # stopped walking parent_unresolved edges and the agent carried on
        # explaining routes through them: it could narrate a mechanism the
        # ranking refuses. Same shape as the agent holding one calibration
        # for every virus -- a fact declared in one place and consumed in
        # one of two.
        _schema = load_schema()
        self._constraints: dict[str, dict] = {
            name: (mp.get("constraints") or {})
            for name, mp in _schema.metapaths.items()}
        # THE SAME CAP THE SCORER AND THE LANDSCAPE USE, read from the schema
        # rather than restated. pathways() must exclude the same oversized
        # pathways pathway_landscape.py excludes, or the agent would describe
        # a route through a 1,389-member pathway that the measured null it
        # quotes never considered.
        self._max_pathway_members: int | None = (
            _schema.retrieval_limits.get("max_pathway_members"))
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

        # After the nodes, because reconciling a calibration's key needs
        # _label to translate a taxon CURIE into the display label the results
        # file is keyed by. calibration.load stores whatever key the CALLER
        # passed -- "SARS-CoV-2" from kgav_ask, "NCBITaxon:2697049" from a
        # caller holding the release -- so neither form can be assumed.
        self._index_calibrations(self._pending_calibration)
        del self._pending_calibration

    def family_of(self, virus: str) -> str | None:
        """The viral family of a taxon, from the graph's own OrganismTaxon.

        ingest_spine writes `family` from config/viruses.yaml, so this is the
        register's answer rather than an inference from the name.
        """
        n = self.nodes.get(virus) or {}
        fam = (n.get("properties") or {}).get("family")
        return str(fam) if fam else None

    # --------------------------------------------------- calibration lookup
    def _same_virus(self, a: str, b: str) -> bool:
        """Whether two virus keys name one virus, across CURIE and label.

        "NCBITaxon:2697049" and "SARS-CoV-2" share no characters, so this
        needs the release's own node labels rather than string surgery --
        which is the same reason calibration.load takes a node_labels map.
        """
        if a == b:
            return True
        for x, y in ((a, b), (b, a)):
            if x in self.nodes and self._label(x) == y:
                return True
        return False

    def _index_calibrations(self, calibration) -> None:
        if calibration is None:
            return
        if hasattr(calibration, "channels"):          # one Calibration
            key = getattr(calibration, "virus", "") or ""
            self._calibrations[key] = calibration
            return
        for key, cal in dict(calibration).items():
            # A MIS-KEYED MAP IS A SILENT WRONG-FIGURE GENERATOR, which is the
            # bug this method exists to end, so it raises rather than being
            # tolerated. Keying SARS-CoV-2's calibration under Zika would
            # reintroduce the exact failure by a different route.
            own = getattr(cal, "virus", "") or ""
            if own and not self._same_virus(key, own):
                raise ValueError(
                    f"calibration keyed under {key!r} was measured on {own!r}; "
                    f"a calibration must be keyed by the virus it measured")
            self._calibrations[key] = cal

    def calibration_for(self, virus: str):
        """The calibration measured ON this virus, or None.

        None rather than a fallback: every caveat has an honest uncalibrated
        branch, and "its reliability is unknown" is weaker than a figure and
        true, where another virus's figure is neither.
        """
        for key, cal in self._calibrations.items():
            if self._same_virus(key, virus):
                return cal
        return None

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
        verdict = guards.check(guards.QueryKind.EXPLAIN, virus,
                               calibration=self.calibration_for(virus))
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
            #
            # THE SAME WARNINGS GO INTO `warns`, not only into the printed
            # route. They used to go only here, so cell_context and
            # host_directed_route never reached ToolResult.warnings, never
            # reached the brief's required_warnings, and verify()'s
            # warning-survival check never required them: the model could drop
            # the SIGMAR1 caveat entirely and the answer still passed. That is
            # the one caveat this project exists to attach -- chloroquine's
            # route is real, well-replicated, and led the field astray in 2020.
            # Both halves were tested and the join between them was not:
            # test_agent_tools asserts the text lands in the right route block,
            # test_agent_verify hands verify() a hand-built caveat directly.
            route_warnings = guards.path_warnings(
                nodes, route, calibration=self.calibration_for(virus),
                landscape=self.landscape, virus=virus)
            # A ROUTE WHOSE OWN EDGE IS CONTRADICTED SAYS SO. Attached per
            # route, like the cell-context caveat above and for the same
            # reason: a contradiction naming one protein printed against a
            # route through another discredits every warning on the page.
            route_warnings = list(route_warnings) + self._contradictions(edges)
            for w in route_warnings:
                lines.append(f"  ! {w}")
            warns.extend(route_warnings)

        warns.append(self._selectivity_for(compound, virus))
        return ToolResult(True, "explain",
                          data=[{"route": r, "n_edges": len(es)} for r, es in paths],
                          verbalised=lines,
                          warnings=list(dict.fromkeys(w for w in warns if w)),
                          edge_ids=ids)

    def _contradictions(self, edges: list[dict]) -> list[str]:
        """Caveats for any INHIBITS edge on this route with an inactive twin.

        The agent loads EVERY edge, held-out predicates included, so the
        contradicting measurement is already in hand -- this reads it rather
        than fetching anything. Scoring correctly never walks it; explaining
        correctly cannot ignore it.
        """
        out: list[str] = []
        for e in edges:
            if e.get("predicate") != "INHIBITS":
                continue
            twins = [t for t in self.out.get(
                ("MEASURED_INACTIVE_AGAINST", e["subject"]), [])
                if t.get("object") == e["object"]]
            if not twins:
                continue
            q = twins[0].get("qualifiers") or {}
            bits = [f"{k}={q[k]:g}" for k in ("ic50_nm", "ec50_nm", "ki_nm")
                    if q.get(k) is not None]
            rel = q.get("relation")
            detail = ((rel or "") + " " + ", ".join(bits)).strip()
            out.append(guards.caveat(
                "contradicted",
                guards.contradicted_caveat(self._label(e["subject"]),
                                           self._label(e["object"]), detail)))
        return out

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
        cl = guards.cell_line_note(self.si_cell.get((compound, virus)),
                                   family=self.family_of(virus))
        return (note + " " + cl).strip() if cl else note

    def _find_paths(self, compound: str, virus: str, limit: int
                    ) -> list[tuple[str, list[dict]]]:
        """Enumerate metapath instances. Bounded and non-revisiting."""
        found: list[tuple[str, list[dict]]] = []

        # The constraints come from the SCHEMA, via baselines.qualifier_ok --
        # the same function dwpc_scores applies -- so a constraint declared
        # once is honoured by both the scorer and the agent. M2's dependency
        # check used to be written out here by hand; it is the schema's
        # HOST_FACTOR_FOR.direction=dependency and is now read from there.
        c1, c2 = self._constraints.get("M1"), self._constraints.get("M2")

        # M1: compound -> viral protein -> gene -> virus
        for e1 in self.out.get(("INHIBITS", compound), []):
            if not qualifier_ok(e1, "INHIBITS", c1):
                continue
            for e2 in self.out.get(("ENCODED_BY", e1["object"]), []):
                if not qualifier_ok(e2, "ENCODED_BY", c1):
                    continue
                for e3 in self.out.get(("BELONGS_TO", e2["object"]), []):
                    if e3["object"] != virus:
                        continue
                    if not qualifier_ok(e3, "BELONGS_TO", c1):
                        continue
                    found.append(("M1 direct-acting", [e1, e2, e3]))
                    if len(found) >= limit:
                        return found

        # M2: compound -> host protein -> virus (dependency only)
        for e1 in self.out.get(("TARGETS", compound), []):
            if not qualifier_ok(e1, "TARGETS", c2):
                continue
            for e2 in self.out.get(("HOST_FACTOR_FOR", e1["object"]), []):
                if e2["object"] != virus:
                    continue
                if not qualifier_ok(e2, "HOST_FACTOR_FOR", c2):
                    continue
                found.append(("M2 host-directed", [e1, e2]))
                if len(found) >= limit:
                    return found
        return found

    def triage(self, compounds: list[str], virus: str) -> ToolResult:
        """Rank compounds you already have.

        Supported because the recall ceiling does not apply: these compounds
        are in the graph by construction.
        """
        verdict = guards.check(guards.QueryKind.TRIAGE, virus,
                               calibration=self.calibration_for(virus))
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
            # made chloroquine look promising in 2020. That reasoning stands
            # on its own; the evidence once cited for it does not.
            # "Selectivity took M1's AUC from 0.729 to 0.826" was measured by
            # filtering positives to compounds with verified selectivity data
            # while leaving negatives whole, which also moved cov_pos 49% ->
            # 69% and left cov_neg at 4%. Most of that gain was annotation
            # density: symmetric, M1 reads 0.660, and 0.499 once same-paper
            # evidence is withheld. Selectivity still comes first because
            # potency alone is the chloroquine failure, not because of 0.826.
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
        warns.append(guards.caveat(
            "ranking",
            # Describes the sort key above, in its actual order. The previous
            # wording said "by direct-acting route count", which is what the
            # ranking deliberately STOPPED doing -- counting routes is what
            # put chloroquine above nirmatrelvir on a screening artifact.
            "compounds are ordered by whether they have any direct-acting "
            "route, then by whether a selectivity index exists, then by higher "
            "selectivity, with the most potent direct-acting measurement as a "
            "final tiebreak. The NUMBER of routes does not affect the order, "
            "host-directed routes contribute nothing to it, and a compound "
            "whose selectivity is unknown ranks below one with a measured "
            "index. "
            + guards.direct_acting_caveat(self.calibration_for(virus))
            + ". Host-directed routes contribute nothing to this order and "
            "are shown for context only."))
        # UNCONDITIONAL, AND THE RANKING CAVEAT ABOVE IS NOT A SUBSTITUTE.
        # That one describes the sort key: which quantities move a compound up
        # the list. This one states what the landscape measured about the
        # route itself. A reader can accept "host-directed routes contribute
        # nothing to this order" and still read the host-directed count beside
        # each compound as weak supporting evidence, which is the claim the
        # degree-matched null refuses. Appended even with no landscape, where
        # it reads "not measured on this release": an unstated null and a
        # measured one must not look the same from outside.
        warns.append(guards.caveat(
            "pathway_route", guards.pathway_caveat(self.landscape, virus)))
        if missing:
            warns.append(f"{len(missing)} requested compound(s) are not in this "
                         f"graph and were not ranked: {', '.join(missing)}")
        return ToolResult(True, "triage", data=scored, verbalised=lines,
                          warnings=warns, note=f"{len(missing)} unranked" if missing else "")

    def _is_viral(self, node_id: str) -> bool:
        """Whether a node is a viral protein, from its own properties.

        A missing node is not viral: an absent id cannot be asserted to be
        anything, and treating it as viral would silently drop it from a
        human-protein census.
        """
        n = self.nodes.get(node_id)
        if not n:
            return False
        return bool((n.get("properties") or {}).get("is_viral"))

    def _dependency_factors(self, virus: str) -> set[str]:
        """Proteins this virus DEPENDS on, by CRISPR screen.

        The direction qualifier is required rather than assumed: a
        restriction factor is the opposite claim, and pooling the two makes
        every count below mean nothing. Same definition as
        pathway_landscape.dependency_factors, because the agent must not
        count a different set from the one the measured null was computed
        over.
        """
        out = set()
        for e in self.inv.get(("HOST_FACTOR_FOR", virus), []):
            if (e.get("qualifiers") or {}).get("direction") == "dependency":
                out.add(e["subject"])
        return out

    def pathways(self, compound: str, virus: str,
                 limit: int = MAX_PATHWAYS) -> ToolResult:
        """Which pathways this compound's targets sit in, and whether the
        virus depends on any member.

        WHY THIS IS NOT explain(). _find_paths walks M1 and M2 only, so the
        pathway metapath the schema declares as M4 is never traversed by the
        agent at all -- a question about pathway membership had no answer
        before this, only a silence that read like a negative. And a path
        sampler would be the wrong answer even once M4 were walked: it
        returns MAX_PATHS instances chosen by iteration order, which tells a
        reader nothing about how many pathways there are, how many were
        skipped, or which of the compound's targets contributed none. The
        question here is a census, so the answer is counted, and every
        absence in it is reported as an absence.

        WHAT IT REFUSES TO IMPLY. Membership is not benefit. Most TARGETS
        edges record potency with no direction, so the compound may activate
        a protein the virus depends on as easily as inhibit it, and the
        caveat says so from the edges actually used. The measured
        degree-matched null rides along too: on this release the drugged
        proteome is no likelier to hit a dependency factor than chance, so a
        pathway hit is context, not evidence.
        """
        verdict = guards.check(guards.QueryKind.EXPLAIN, virus,
                               calibration=self.calibration_for(virus))
        base_warnings = list(verdict.warnings) + [guards.caveat(
            "pathway_route", guards.pathway_caveat(self.landscape, virus))]

        if compound not in self.nodes:
            return ToolResult(False, "pathways", warnings=base_warnings,
                              note=f"{compound} is not in this graph, so it "
                                   f"has no targets to place in a pathway.")

        # Human targets only. A viral protein can be hit by a compound and is
        # never a host pathway member, so including one would put a viral
        # protein in a census of the human proteome.
        target_edges = [e for e in self.out.get(("TARGETS", compound), [])
                        if not self._is_viral(e["object"])]
        if not target_edges:
            return ToolResult(
                False, "pathways", warnings=base_warnings,
                note=f"{self._label(compound)} has no human target in this "
                     f"graph, so there is no pathway to place it in. That is "
                     f"a property of the source and not a finding about the "
                     f"compound: most small molecules here carry an activity "
                     f"measurement and no target annotation.")

        n_directed = sum(
            1 for e in target_edges
            if (e.get("qualifiers") or {}).get("direction"))
        deps = self._dependency_factors(virus)
        cap = self._max_pathway_members

        # pathway id -> {members, this compound's targets in it}
        hits: dict[str, dict] = {}
        oversized: set[str] = set()
        no_pathway: list[str] = []
        edge_ids = [e["_id"] for e in target_edges]

        for te in target_edges:
            prot = te["object"]
            memberships = self.out.get(("PARTICIPATES_IN", prot), [])
            placed = False
            for pe in memberships:
                pw = pe["object"]
                if pw not in hits and pw not in oversized:
                    members = {m["subject"]
                               for m in self.inv.get(("PARTICIPATES_IN", pw), [])}
                    # THE SAME CAP THE LANDSCAPE APPLIED. Quoting a measured
                    # null computed over capped pathways while reporting an
                    # uncapped one would be two different questions under one
                    # answer.
                    if cap is not None and len(members) > cap:
                        oversized.add(pw)
                        continue
                    hits[pw] = {"members": members, "targets": set(),
                                "edges": []}
                if pw in oversized:
                    continue
                hits[pw]["targets"].add(prot)
                hits[pw]["edges"].append(pe["_id"])
                placed = True
            if not placed:
                no_pathway.append(prot)

        rows = []
        for pw, h in hits.items():
            dep_members = sorted(h["members"] & deps)
            rows.append({
                "pathway": pw, "label": self._label(pw),
                "n_members": len(h["members"]),
                "targets_in_pathway": sorted(h["targets"]),
                "dependency_members": dep_members,
                "n_dependency_members": len(dep_members),
            })
        # Pathways the virus depends on something in, first; then by how much
        # of the pathway is a dependency factor. NOT a ranking of relevance:
        # ordering decides what survives the limit, and a reader is better
        # served by seeing the hits than by seeing the first N by node id.
        rows.sort(key=lambda r: (-r["n_dependency_members"], -len(
            r["targets_in_pathway"]), r["pathway"]))
        with_dep = [r for r in rows if r["n_dependency_members"]]
        shown = rows[:limit]

        lines = [
            (f"{self._label(compound)} has {len(target_edges)} human "
             f"target(s) in this graph."),
            (f"{len(rows)} pathway(s) at or under the {cap}-member cap "
             f"contain at least one of them; {len(with_dep)} of those "
             f"contain at least one protein {self._label(virus)} depends "
             f"on."),
        ]
        if oversized:
            lines.append(
                f"{len(oversized)} further pathway(s) were EXCLUDED for "
                f"exceeding the {cap}-member cap. They are not absent from "
                f"the graph, only from this answer and from the measured "
                f"null it quotes.")
        if no_pathway:
            lines.append(
                f"{len(no_pathway)} target(s) sit in no pathway at or under "
                f"the cap: {', '.join(self._label(p) for p in no_pathway)}.")
        if not with_dep:
            lines.append(
                f"No pathway containing a target of this compound contains a "
                f"protein {self._label(virus)} is known to depend on. That is "
                f"an absence in the CRISPR layer as much as a statement about "
                f"this compound.")
        for r in shown:
            dep = (", ".join(self._label(d) for d in r["dependency_members"])
                   if r["dependency_members"] else "none")
            lines.append(
                f"  {r['label']} ({r['n_members']} members) — via "
                f"{', '.join(self._label(t) for t in r['targets_in_pathway'])}"
                f"; dependency factors in it: {dep}")
        if len(rows) > limit:
            lines.append(f"  … {len(rows) - limit} further pathway(s) not "
                         f"shown, of which "
                         f"{sum(1 for r in rows[limit:] if r['n_dependency_members'])} "
                         f"contain a dependency factor.")

        for h in hits.values():
            edge_ids.extend(h["edges"])
        warns = base_warnings + [guards.caveat(
            "target_direction",
            guards.target_direction_caveat(n_directed, len(target_edges)))]
        return ToolResult(True, "pathways", data=rows, verbalised=lines,
                          warnings=warns, edge_ids=edge_ids,
                          note=f"{len(oversized)} pathway(s) over the cap"
                               if oversized else "")

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
