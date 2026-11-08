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

import re
from dataclasses import dataclass, field
from typing import Any

from kgav.agent import guards
from kgav.agent.verify import (
    EDGE_ID_RE,
    FACT_ID_RE,
    caveat_labels,
    dropped_warnings,
)
from kgav.agent.tools import GraphTools, ToolResult
from kgav.virus_register import MIN_SEARCH_NAME

MAX_BRIEF_EDGES = 100

#: Entities resolved out of one question. It bounds DOWNSTREAM work -- one
#: explain() per compound, one triage row per compound -- not the scan,
#: which runs over the whole name index either way because _mentions is
#: evaluated before any cap is consulted.
#:
#: IT USED TO BE A BARE `if len(found) >= 6: break` INSIDE the
#: longest-first loop, and that is a different thing entirely. Drug names
#: are long and virus names are short, so seven named compounds filled
#: the quota before "zika" (four characters) was ever reached:
#:
#:   "rank hydroxychloroquine, dexamethasone, remdesivir, nitazoxanide,
#:    camostat mesylate, ivermectin, favipiravir for Zika"
#:
#: resolved six compounds, dropped ivermectin AND Zika, fell back to the
#: default virus, and refused -- about SARS-CoV-2. The same question with
#: five compounds answers correctly about Zika. A cap applied during a
#: length-ordered scan does not bound the work, it decides which KIND of
#: entity survives.
MAX_RESOLVED = 6

#: Citable edge ids a brief may carry. SEPARATE FROM THE FACT BOUND, and
#: it had no bound at all: `MAX_BRIEF_EDGES` capped `facts` while
#: `citable_edge_ids` grew with whatever the tools returned. On a graph
#: with 300 dependency factors and 400 pathways, pathways() returns 29
#: census lines and 1,500 PARTICIPATES_IN ids, and Brief.to_prompt joins
#: every one of them into a single line:
#:
#:     prompt     14,579 chars
#:     the ids     9,762 chars   -- 67% of it
#:     the facts    3,329 chars
#:
#: Two thirds of the prompt bare identifiers, which is the
#: context-window problem MAX_BRIEF_EDGES exists to prevent, in the half
#: of the brief it did not cover.
#:
#: Higher than the fact bound on purpose. One shown fact can legitimately
#: carry several ids -- assembly merges per-study edges and a path is a
#: chain -- so a cap at the fact bound would discard ids whose facts are
#: in front of the model.
MAX_BRIEF_IDS = 400


def shown_edge_ids(facts: list[str]) -> set[str]:
    """Edge ids the model was actually shown a fact for.

    TWO CITATION SHAPES, and reading only one is a bug I shipped into a
    test run. tools._verbalise renders a fact's provenance as
    "[chembl, PMID:1, 2021; E335139]", while EDGE_ID_RE matches a
    bracket whose whole content is ids -- the shape a MODEL writes. Using
    EDGE_ID_RE alone on the facts finds nothing in any real fact.
    """
    joined = " ".join(facts or [])
    out = set(FACT_ID_RE.findall(joined))
    for group in EDGE_ID_RE.findall(joined):
        out |= {c.strip() for c in group.split(",")}
    return out


def bounded_ids(facts: list[str], ids: list[str]) -> tuple[list[str], int]:
    """`(citable ids for the brief, how many were dropped)`, deduplicated.

    PREFERRING THE IDS WHOSE FACT THE MODEL CAN SEE. An id with no fact
    behind it is provenance for a COUNTED claim: real, and useless to a
    model deciding what to cite -- and verify() no longer lets one exempt
    a sentence from the uncited-claim check.

    THIS PREFERENCE IS DEFENSIVE AND NOT REACHABLE TODAY, which is worth
    saying rather than implying. build_brief takes one branch per
    question, so a brief's results are the same shape: a pathway brief
    carries census lines with no inline ids, and an evidence or explain
    brief carries inline ids and stays far under the cap. A mixed brief
    -- one id-poor result supplying the overflow and an id-rich one
    behind it -- is what the ordering protects against, and no current
    branch builds one. A mutation that replaced this with a plain
    `ids[:MAX_BRIEF_IDS]` survived every test written through
    build_brief for exactly that reason, so it is tested here directly
    instead of through a fixture that cannot produce the condition.
    """
    ids = list(dict.fromkeys(ids or []))
    if len(ids) <= MAX_BRIEF_IDS:
        return ids, 0
    shown = shown_edge_ids(facts)
    keep = [i for i in ids if i in shown][:MAX_BRIEF_IDS]
    keep += [i for i in ids if i not in shown][:MAX_BRIEF_IDS - len(keep)]
    return keep, len(ids) - len(keep)


#: A name must sit on a word boundary to count as a mention.
#:
#: THE BUG THIS CLOSES. Resolution was `if name not in text` -- a bare
#: substring test over every node name of four characters or more. The
#: index holds brand names, and "Depen" is a brand name for
#: penicillamine, so the question "which pathways does SARS-CoV-2
#: DEPENd on?" resolved penicillamine and the agent answered about a
#: compound nobody had mentioned, with every caveat attached as though
#: it were the question asked. A node profile for a misread question is
#: confusing; a confident, fully-caveated answer about the wrong
#: compound is the kind of thing that reaches a methods section.
#:
#: The boundary is "not alphanumeric" rather than \b, because \b is
#: defined against \w and the names here end in digits and punctuation:
#: \bsars-cov-2\b fails to match "sars-cov-2." on the trailing "2"
#: against "." -- it is a boundary, so \b is satisfied -- but
#: \bchloroquine\b inside "hydroxychloroquine" is NOT suppressed by \b
#: at the front, since "y" and "c" are both word characters and there
#: is no boundary to assert. Testing the neighbouring character
#: directly handles both, and keeps hyphens and digits inside a name
#: from being read as separators.
def _mentions(text: str, name: str) -> bool:
    """Whether `name` occurs in `text` other than inside a longer word."""
    return re.search(rf"(?<![a-z0-9]){re.escape(name)}(?![a-z0-9])",
                     text) is not None


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
                "a labelled sentence: the label in capitals, then a COLON, "
                "then the warning. An answer missing any of these labels "
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
            # THE COLON IS SAID OUT LOUD because it is enforced (fix-list
            # A5). verify._label_is_present requires "LABEL:", not the bare
            # word, so an answer that works the label into a sentence --
            # "the evidence here is thin" for EVIDENCE -- is discarded. The
            # warnings below are already in that form, so copying one is
            # enough; asking for it explicitly costs a line and saves a
            # retry, which is the same reason dropped_warnings exists.
            parts += ["", "WARNINGS THAT MUST APPEAR IN YOUR ANSWER, each as "
                          "its LABEL then a colon then the text:"]
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
    def _resolve_mentions(
            self, question: str) -> tuple[dict[str, str], str, list[str]]:
        """Longest-match entity resolution over the name index.

        `(resolved, problem, omitted)`. `problem` means the question
        cannot be answered at all; `omitted` names entities that matched
        and did not fit MAX_RESOLVED, which is a caveat and not a
        refusal.

        Longest first so "hydroxychloroquine" is not resolved as
        "chloroquine" -- the same substring hazard that required Vero E6 to be
        matched before Vero and SARS-CoV-2 before SARS-CoV.

        A TAXON IS NEVER THE THING DROPPED. The cap used to break out of
        this loop, which is ordered by name LENGTH, so the virus -- the
        shortest name in a question full of drug names, and the one
        entity that decides the calibration, the routes, every caveat
        and the refusal itself -- was the first casualty. The cap now
        applies after matching and only to non-taxon entities.
        """
        text = (question or "").lower()
        found: dict[str, str] = {}
        for name in sorted(self.tools.names, key=len, reverse=True):
            if len(name) < MIN_SEARCH_NAME:
                continue
            # Suppress a shorter name only where every occurrence sits inside
            # a longer match. "chloroquine" is a substring of
            # "hydroxychloroquine", but a question naming both must resolve
            # both -- the earlier rule silently dropped one of four requested
            # compounds from a ranking.
            if not _mentions(text, name):
                continue
            covered = text
            for seen in found:
                covered = covered.replace(seen, " " * len(seen))
            if _mentions(covered, name):
                found[name] = self.tools.names[name]
        if not found:
            return {}, ("No entity in this question matches a node in the graph. "
                        "Name a compound, protein or virus the graph contains."), []

        def _is_taxon(nid: str) -> bool:
            return self.tools.nodes.get(nid, {}).get("class") == "OrganismTaxon"

        taxa = {n: i for n, i in found.items() if _is_taxon(i)}
        # THE REMAINDER IN THE ORDER THE QUESTION ASKED, not the order
        # this loop matched in. Matching is by descending name LENGTH,
        # which is the right order for resolving overlapping names and an
        # arbitrary one for deciding whose compound gets dropped: it kept
        # "camostat mesylate" over "remdesivir" for no reason a user
        # could predict. First-named, first-kept at least matches what
        # they wrote.
        rest = sorted((n for n in found if n not in taxa),
                      key=lambda n: text.find(n))
        room = max(0, MAX_RESOLVED - len(taxa))
        keep = dict(taxa)
        for n in rest[:room]:
            keep[n] = found[n]
        # REPORTED, NOT SILENT. Dropping two of seven requested compounds
        # from a ranking is a bounded harm; doing it without saying so is
        # not, because the answer then looks like a complete comparison
        # of everything the user asked about.
        omitted = [n for n in rest[room:]]
        return keep, "", omitted

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

        resolved, problem, omitted = self._resolve_mentions(question)
        trace.append(f"resolved={len(resolved)}")
        if omitted:
            trace.append(f"omitted={len(omitted)}")
        named = self._virus_in(resolved)
        virus = named or self.default_virus
        trace.append(f"virus={virus}{'' if named else ' (assumed)'}")

        # A DEFAULT VIRUS THAT IS NOT IN THIS GRAPH IS NOT AN ANSWER.
        # guards.check would go on to refuse about it, and the refusal
        # states figures -- "it holds 0 antiviral measurement(s)" -- and
        # names the virus through _label, which returns the raw CURIE for
        # a node that is absent. So a user who named a virus the cap then
        # dropped could be told "This release cannot answer that about
        # NCBITaxon:2697049", a sentence that is wrong about a virus they
        # never mentioned, in a vocabulary no reader shares.
        if named is None and virus not in self.tools.nodes:
            return Brief(
                question=question, intent=intent.value, allowed=True,
                clarification=(
                    "This question does not name a virus that this release "
                    "holds, and the virus decides everything about the "
                    "answer: which routes exist, what has been measured, and "
                    "which caveats apply. Name one of "
                    + ", ".join(sorted(
                        self.tools._label(nid)
                        for nid, d in self.tools.nodes.items()
                        if d.get("class") == "OrganismTaxon")) + "."),
                trace=trace + ["default_virus_absent"])

        # This virus's own calibration, so the caveat quotes what was
        # measured on it rather than a sentence asserting things about
        # other viruses.
        verdict = guards.check(
            intent, virus, calibration=self.tools.calibration_for(virus),
            routes=self.tools.routes_for(virus))

        # AN ASSUMED VIRUS IS SAID OUT LOUD. The fallback above is silent,
        # and the virus decides the calibration, the caveats and the
        # cell-context biology -- so a question meaning Zika answered about
        # SARS-CoV-2 is wrong in every one of those at once. Latent while the
        # graph held one virus; live since it holds eight across two families.
        #
        # AND ON A REFUSAL ABOVE ALL. The condition was `named is None and
        # verdict.allowed`, so the one path where the substitution does the
        # most damage -- the agent refusing, in detail and with figures,
        # about a virus the user never named -- was the path that said
        # nothing about it. The caveat is built either way now, and the
        # refusal below carries it in its TEXT, because a refusal reaches
        # the user as brief.refusal and nothing else does.
        if named is None:
            others = sorted(
                self.tools._label(nid) for nid, d in self.tools.nodes.items()
                if d.get("class") == "OrganismTaxon" and nid != virus)
            verdict.warnings.append(guards.caveat(
                "virus_assumed",
                guards.virus_assumed_caveat(self.tools._label(virus), others)))
        if omitted:
            verdict.warnings.append(guards.caveat(
                "entities_omitted",
                guards.omitted_entities_caveat(omitted, len(resolved))))
        trace.append(f"allowed={verdict.allowed}")
        if not verdict.allowed:
            # The refusal is produced BEFORE retrieval, so no partial answer
            # can leak out alongside it.
            #
            # The caveats go in the TEXT and in `warnings`. A refusal is
            # delivered as brief.refusal verbatim -- answer() sets
            # out["answer"] = brief.refusal and returns -- so a warning
            # carried only in the field would reach nobody.
            reason = verdict.reason
            if verdict.warnings:
                reason = "\n\n".join([*map(str, verdict.warnings), reason])
            return Brief(question=question, intent=intent.value, allowed=False,
                         refusal=reason,
                         warnings=list(dict.fromkeys(
                             str(w) for w in verdict.warnings)),
                         trace=trace)

        if problem:
            return Brief(question=question, intent=intent.value, allowed=True,
                         clarification=problem, trace=trace)

        compounds = self._compounds_in(resolved)
        results: list[ToolResult] = []

        if intent is guards.QueryKind.PATHWAY and compounds:
            results += [self.tools.pathways(c, virus) for c in compounds]
        elif intent is guards.QueryKind.PATHWAY and virus:
            # THE OTHER DIRECTION. Without this, "which pathways does
            # Zika depend on" fell through every branch below to a node
            # profile -- an answer about the taxon node, to a question
            # about pathways, with nothing saying the question had not
            # been understood.
            results.append(self.tools.virus_pathways(virus))
        elif intent is guards.QueryKind.TRIAGE and compounds:
            results.append(self.tools.triage(compounds, virus))
        elif intent is guards.QueryKind.PROFILE:
            results += [self.tools.profile(nid) for nid in resolved.values()]
        elif intent is guards.QueryKind.EVIDENCE and len(resolved) >= 2:
            ids = list(resolved.values())
            # Taxa are inserted into `resolved` first, but edges run
            # compound -> virus, so the question's order says nothing about
            # edge direction. Try both before reporting "no edge".
            hit = self.tools.evidence(ids[0], ids[1])
            if not hit.ok:
                hit = self.tools.evidence(ids[1], ids[0])
            results.append(hit)
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

        # THE FIRST RESULT IS ADMITTED WHOLE, which the `facts and` above
        # says and the external audit read as unbounded. It is bounded in
        # fact, by arithmetic rather than by this guard: the largest
        # single result any tool can return is explain()'s MAX_PATHS (20)
        # paths of one header plus M1's three edges = 80 lines, and
        # _find_paths walks only M1 and M2. Every other tool caps lower.
        # Dropping the first result instead would leave an EMPTY brief,
        # which is worse than an over-long one -- so it is admitted and
        # the overflow is RECORDED, where a raised per-tool cap becomes
        # visible instead of silently doubling the prompt.
        if len(facts) > MAX_BRIEF_EDGES:
            trace.append(f"brief_over_fact_bound={len(facts)}")

        # IDS BOUNDED TOO. See bounded_ids.
        kept, n_dropped = bounded_ids(facts, ids)
        if n_dropped:
            trace.append(f"truncated_ids={n_dropped}")
        ids = kept

        return Brief(question=question, intent=intent.value, allowed=True,
                     facts=facts,
                     warnings=list(dict.fromkeys(w for w in warns if w)),
                     citable_edge_ids=ids,
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
