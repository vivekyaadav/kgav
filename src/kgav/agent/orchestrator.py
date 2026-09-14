"""The orchestrator: a fixed pipeline, not an agent loop.

WHY NOT A LOOP. A model that decides its own next action is non-deterministic
by construction, and a researcher who runs the same query twice and gets two
different answers cannot cite the tool. Reproducibility is not a nice property
here -- it is the difference between something that can appear in a methods
section and something that cannot.

Per-step reliability also compounds. Five chained decisions at 90% each leave
59%, and 90% is generous for a 14B model reasoning over a 21-predicate schema.

THE FIVE STAGES, and which are deterministic:

  1. classify      deterministic  keyword intent classification
  2. guard         deterministic  refuse what the evaluation cannot support
  3. resolve       deterministic  free text -> node ids, ambiguity returned
  4. retrieve      deterministic  fixed tool calls, bounded output
  5. synthesise    MODEL          narrate the retrieved subgraph

Only stage 5 is generative. The model never chooses what to fetch, never
decides whether a question may be answered, and never resolves an entity --
those are the decisions where a plausible wrong answer is worse than an error.

WHAT THE MODEL RECEIVES. A structured brief: the retrieved facts already
verbalised with provenance, the warnings that must survive into the answer,
and the list of edge identifiers it is permitted to cite. It is not given the
graph, a query language, or discretion about scope.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from kgav.agent import guards
from kgav.agent.tools import GraphTools, ToolResult

MAX_BRIEF_EDGES = 100


@dataclass
class Brief:
    """Everything the model is given, and nothing else."""
    question: str
    intent: str
    allowed: bool
    refusal: str = ""
    facts: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    citable_edge_ids: list[str] = field(default_factory=list)
    resolved: dict[str, str] = field(default_factory=dict)
    clarification: str = ""
    trace: list[str] = field(default_factory=list)

    def to_prompt(self) -> str:
        """The brief as text for a language model.

        The instruction block is deliberately restrictive. The model's job is
        to narrate retrieved facts, not to add domain knowledge: anything it
        contributes from its own weights is unciteable and unverifiable, and
        the verification stage will reject it.
        """
        if not self.allowed:
            return self.refusal

        parts = [
            "You are reporting findings from a curated knowledge graph.",
            "",
            "RULES:",
            ("- Use ONLY the facts listed below. Add nothing from your own "
            "knowledge, however confident you are."),
            ("- Cite the edge identifier in square brackets after every factual "
            "claim, e.g. [E335139]."),
            ("- Every warning below MUST appear in your answer. They are not "
            "optional context; they are what makes the facts interpretable."),
            "- If the facts do not answer the question, say so plainly.",
            "",
            f"QUESTION: {self.question}",
            "",
            "FACTS:",
        ]
        parts += [f"  {f}" for f in self.facts] or ["  (none retrieved)"]
        if self.warnings:
            parts += ["", "WARNINGS THAT MUST APPEAR IN YOUR ANSWER:"]
            parts += [f"  - {w}" for w in self.warnings]
        parts += ["", f"CITABLE EDGE IDS: {', '.join(self.citable_edge_ids) or 'none'}"]
        return "\n".join(parts)


class Orchestrator:
    def __init__(self, tools: GraphTools, default_virus: str = guards.SUPPORTED_VIRUS):
        self.tools = tools
        self.default_virus = default_virus

    # ------------------------------------------------------------- stage 3
    def _resolve_mentions(self, question: str) -> tuple[dict[str, str], str]:
        """Longest-match entity resolution over the name index.

        Longest first so "hydroxychloroquine" is not resolved as
        "chloroquine" -- the same substring hazard that required Vero E6 to be
        matched before Vero and SARS-CoV-2 before SARS-CoV.
        """
        text = (question or "").lower()
        found: dict[str, str] = {}
        for name in sorted(self.tools.names, key=len, reverse=True):
            if len(name) < 4:
                continue
            # Suppress a shorter name only where every occurrence sits inside
            # a longer match. "chloroquine" is a substring of
            # "hydroxychloroquine", but a question naming both must resolve
            # both -- the earlier rule silently dropped one of four requested
            # compounds from a ranking.
            if name not in text:
                continue
            covered = text
            for seen in found:
                covered = covered.replace(seen, " " * len(seen))
            if name in covered:
                found[name] = self.tools.names[name]
            if len(found) >= 6:
                break
        if not found:
            return {}, ("No entity in this question matches a node in the graph. "
                        "Name a compound, protein or virus the graph contains.")
        return found, ""

    def _virus_in(self, resolved: dict[str, str]) -> str | None:
        for nid in resolved.values():
            if self.tools.nodes.get(nid, {}).get("class") == "OrganismTaxon":
                return nid
        return None

    def _compounds_in(self, resolved: dict[str, str]) -> list[str]:
        return [nid for nid in resolved.values()
                if self.tools.nodes.get(nid, {}).get("class") == "SmallMolecule"]

    # ------------------------------------------------------------ pipeline
    def build_brief(self, question: str) -> Brief:
        trace: list[str] = []

        intent = guards.classify_intent(question)
        trace.append(f"intent={intent.value}")

        resolved, problem = self._resolve_mentions(question)
        trace.append(f"resolved={len(resolved)}")
        virus = self._virus_in(resolved) or self.default_virus

        verdict = guards.check(intent, virus)
        trace.append(f"allowed={verdict.allowed}")
        if not verdict.allowed:
            # The refusal is produced BEFORE retrieval, so no partial answer
            # can leak out alongside it.
            return Brief(question=question, intent=intent.value, allowed=False,
                         refusal=verdict.reason, trace=trace)

        if problem:
            return Brief(question=question, intent=intent.value, allowed=True,
                         clarification=problem, trace=trace)

        compounds = self._compounds_in(resolved)
        results: list[ToolResult] = []

        if intent is guards.QueryKind.TRIAGE and compounds:
            results.append(self.tools.triage(compounds, virus))
        elif intent is guards.QueryKind.PROFILE:
            results += [self.tools.profile(nid) for nid in resolved.values()]
        elif intent is guards.QueryKind.EVIDENCE and len(resolved) >= 2:
            ids = list(resolved.values())
            results.append(self.tools.evidence(ids[0], ids[1]))
        elif compounds:
            results += [self.tools.explain(c, virus) for c in compounds]
        else:
            results += [self.tools.profile(nid) for nid in resolved.values()]

        facts, warns, ids = [], list(verdict.warnings), []
        for r in results:
            trace.append(f"tool={r.kind} ok={r.ok}")
            facts.extend(r.verbalised)
            warns.extend(r.warnings)
            ids.extend(r.edge_ids)
            if not r.ok and r.note:
                facts.append(f"(no result: {r.note})")

        # Bound the brief. An unbounded subgraph fills the context window and
        # the model starts summarising rather than reporting.
        if len(facts) > MAX_BRIEF_EDGES:
            facts = facts[:MAX_BRIEF_EDGES]
            facts.append(f"(truncated to {MAX_BRIEF_EDGES} facts)")

        return Brief(question=question, intent=intent.value, allowed=True,
                     facts=facts,
                     warnings=list(dict.fromkeys(w for w in warns if w)),
                     citable_edge_ids=list(dict.fromkeys(ids)),
                     resolved={k: v for k, v in resolved.items()},
                     trace=trace)

    # ------------------------------------------------------------ stage 5
    def answer(self, question: str, model: Any = None) -> dict:
        """Full pipeline. `model` is any callable taking a prompt.

        With no model supplied the brief is returned unsynthesised, which is
        the mode to use when the facts matter more than the prose -- and is
        what the tests exercise, since a test of the pipeline should not
        depend on a model's wording.
        """
        brief = self.build_brief(question)
        out: dict[str, Any] = {"question": question, "intent": brief.intent,
                               "allowed": brief.allowed, "trace": brief.trace,
                               "facts": brief.facts, "warnings": brief.warnings,
                               "citable_edge_ids": brief.citable_edge_ids}
        if not brief.allowed:
            out["answer"] = brief.refusal
            return out
        if brief.clarification:
            out["answer"] = brief.clarification
            return out
        if model is None:
            out["answer"] = None
            out["prompt"] = brief.to_prompt()
            return out
        out["answer"] = model(brief.to_prompt())
        return out
