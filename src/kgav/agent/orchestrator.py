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
from kgav.agent.verify import caveat_labels, dropped_warnings
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

    def to_prompt(self, must_include: list | None = None) -> str:
        """The brief as text for a language model.

        `must_include` names caveats a previous attempt dropped. They are
        quoted at the TOP of the prompt rather than appended, because the
        rules below already say every warning must appear and a 7B model
        dropped one anyway -- repeating the same instruction further down is
        not a different instruction. Naming the specific label it omitted is.

        The instruction block is deliberately restrictive. The model's job is
        to narrate retrieved facts, not to add domain knowledge: anything it
        contributes from its own weights is unciteable and unverifiable, and
        the verification stage will reject it.
        """
        if not self.allowed:
            return self.refusal

        parts = []
        if must_include:
            labels = caveat_labels(must_include)
            parts += [
                "YOUR PREVIOUS ANSWER WAS REJECTED. It omitted "
                + ("these required warnings: " if len(labels) > 1
                   else "this required warning: ")
                + "; ".join(labels) + ".",
                "",
                "Write the answer again and include each of them verbatim as "
                "a labelled sentence. An answer missing any of these labels "
                "is discarded, so a shorter answer is not a better one.",
                "",
            ]
        parts += [
            "You are reporting findings from a curated knowledge graph.",
            "",
            "RULES:",
            ("- Use ONLY the facts listed below. Add nothing from your own "
            "knowledge, however confident you are."),
            ("- Cite the edge identifier in square brackets after every "
            "factual claim, using ONLY identifiers from the list at the end. "
            "If that list says none, cite NOTHING: do not invent an "
            "identifier and do not reuse one from these instructions."),
            ("- Every warning below MUST appear in your answer. They are not "
            "optional context; they are what makes the facts interpretable."),
            ("- Keep each warning's LABEL -- the capitalised words before its "
            "colon -- exactly as written. You may reword the rest of the "
            "warning; the label is how it is checked."),
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
        if self.citable_edge_ids:
            parts += ["", f"CITABLE EDGE IDS: {', '.join(self.citable_edge_ids)}"]
        else:
            # Printing "none" beside an example identifier invited the model to
            # reuse the example: given an empty list it cited E335139 -- the id
            # that used to appear in the instruction above -- seven times, for
            # facts that have no edges at all. The example is now gone and the
            # empty case says explicitly what to do.
            parts += ["", ("CITABLE EDGE IDS: none. These facts are summary "
                           "counts, not individual assertions. Write the "
                           "answer with NO square-bracket citations.")]
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

        # Bound the brief. An unbounded subgraph fills the context window and
        # the model starts summarising rather than reporting.
        #
        # TRUNCATION TAKES THE IDS WITH THE FACTS. It used to cut `facts` and
        # leave `citable_edge_ids` whole, so past the limit the model held
        # identifiers for facts it had never seen -- and verify()'s
        # invented-citation check passed them, because they WERE in the
        # allowed list. The result was a claim with no retrieved evidence
        # behind it, carrying a real identifier, reported as verified: the
        # failure ab86858 set out to close, arriving by a different route.
        #
        # Whole results are kept or dropped rather than cutting mid-result,
        # because a tool's facts and its edge ids correspond as a set and not
        # line by line: a route header and a warning line are facts with no id
        # of their own, so no index into one list addresses the other.
        facts, warns, ids = [], list(verdict.warnings), []
        dropped = 0
        for r in results:
            trace.append(f"tool={r.kind} ok={r.ok}")
            warns.extend(r.warnings)          # warnings are never truncated
            lines = list(r.verbalised)
            if not r.ok and r.note:
                lines.append(f"(no result: {r.note})")
            if facts and len(facts) + len(lines) > MAX_BRIEF_EDGES:
                dropped += 1
                continue
            facts.extend(lines)
            ids.extend(r.edge_ids)
        if dropped:
            facts.append(f"({dropped} further result(s) omitted to bound the "
                         f"brief at {MAX_BRIEF_EDGES} facts; their evidence is "
                         f"not cited below)")
            trace.append(f"truncated_results={dropped}")

        return Brief(question=question, intent=intent.value, allowed=True,
                     facts=facts,
                     warnings=list(dict.fromkeys(w for w in warns if w)),
                     citable_edge_ids=list(dict.fromkeys(ids)),
                     resolved={k: v for k, v in resolved.items()},
                     trace=trace)

    # ------------------------------------------------------------ stage 5
    def answer(self, question: str, model: Any = None,
               max_attempts: int = 2) -> dict:
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
            out["clarification"] = True
            return out
        if model is None:
            out["answer"] = None
            out["prompt"] = brief.to_prompt()
            return out

        # ONE RETRY, AND ONLY FOR A DROPPED CAVEAT. Printing the answer with
        # the caveat bolted underneath is the behaviour kgav_ask's docstring
        # rejects -- the fluent prose is what sticks and the caveat is read
        # second, if at all. Withholding is right, but withholding on the
        # first attempt made it the normal outcome for every host-directed
        # question, so the model gets one chance with the omitted label named.
        # Nothing else is retried: an invented citation or a fabricated number
        # is not an oversight to nudge, and re-rolling until a check passes is
        # how a verifier becomes a sampler.
        prompt = brief.to_prompt()
        out["prompt"] = prompt
        text = model(prompt)
        out["attempts"] = 1
        missing = dropped_warnings(text, brief.warnings)
        if missing and max_attempts > 1:
            retry_prompt = brief.to_prompt(must_include=missing)
            out["retry_prompt"] = retry_prompt
            out["dropped_on_first_attempt"] = [str(w) for w in missing]
            text = model(retry_prompt)
            out["attempts"] = 2
        out["answer"] = text
        return out
