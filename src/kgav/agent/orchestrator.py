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



#: How far a typed word may sit from a curated virus name and still be
#: offered as a suggestion. TWO, and never applied automatically: the four
#: dengue serotypes are all ONE edit from each other, and DENV-2 is the
#: only cell in this project that survives a power analysis while DENV-1
#: cannot be powered at any arm size. Correcting silently between them
#: would invert the conclusion with nothing on screen to show it.
MAX_SUGGEST_EDITS = 2


def _or_list(items) -> str:
    """a, b, c or d -- so four candidates do not read as three "or"s."""
    items = list(items)
    if len(items) <= 2:
        return " or ".join(items)
    return ", ".join(items[:-1]) + " or " + items[-1]


def _edits(a: str, b: str) -> int:
    """Levenshtein distance. Small inputs; no dependency worth adding."""
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1,
                           prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def near_names(question: str, candidates: dict[str, str],
               max_edits: int = MAX_SUGGEST_EDITS) -> dict[str, list[str]]:
    """{typed phrase: [curated names within max_edits]}, never applied.

    SUGGESTED, NOT CORRECTED, and the distinction is the whole point.
    Measured on this release: 6,499 pairs of DIFFERENT compounds have
    names one character apart -- "penicillin g" and "penicillin v",
    "pentobarbital" and "phenobarbital", "folic acid" and "oleic acid" --
    and all four dengue serotypes are one edit from each other. An
    autocorrect is a silent substitution of the subject, which is the
    defect class this function exists to replace, not to re-enter by a
    friendlier route.

    So every candidate within range is returned, and a caller that finds
    more than one must say so rather than choose.
    """
    words = re.findall(r"[a-z0-9][a-z0-9-]*", (question or "").lower())
    phrases = set(words)
    for n in (2, 3):
        phrases |= {" ".join(words[i:i + n])
                    for i in range(len(words) - n + 1)}
    out: dict[str, list[str]] = {}
    for phrase in phrases:
        if len(phrase) < 4:
            continue
        # TWO EDITS ONLY FOR SOMETHING LONG ENOUGH TO HAVE MEANT IT.
        # A flat budget of two let "the" suggest TBEV -- two edits is
        # two thirds of a three-letter word -- so an ordinary English
        # sentence produced a clarification about a virus nobody typed.
        # Every real near-miss measured here is ONE edit: yellew/yellow,
        # zikaa/zika, denv-5/denv-1.
        budget = max_edits if len(phrase) >= 6 else 1
        hits = {name for name in candidates
                if abs(len(name) - len(phrase)) <= budget
                and 0 < _edits(phrase, name) <= budget}
        # A WHOLE WORD INSIDE A CURATED NAME. "nile" returned nothing
        # while "west nile" worked, because the two differ by five
        # characters and no edit budget covers that. Matching a whole
        # word -- not a prefix, not a substring -- suggests West Nile
        # from "nile" and both coronaviruses from "sars", which is the
        # right answer in each case because it lists every candidate
        # rather than choosing.
        hits |= {name for name in candidates
                 if name != phrase
                 and re.search(rf"(?<![a-z0-9]){re.escape(phrase)}"
                               rf"(?![a-z0-9])", name)}
        hits = sorted(hits)
        if hits:
            out[phrase] = hits
    return out


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
            # THE FLOOR, EXCEPT FOR A CURATED VIRUS NAME. _mentions below
            # matches on word boundaries, so "jev" cannot hit inside
            # "jevons"; the floor exists for graph labels, which are not
            # curated and are numerous. Without this exemption JEV -- the
            # name every table in this project prints -- resolved to
            # nothing and the question was answered about SARS-CoV-2.
            if (len(name) < MIN_SEARCH_NAME
                    and name not in self.tools._register_names):
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

    def _viruses_in(self, resolved: dict[str, str]) -> list[str]:
        """EVERY taxon the question named, in the order it named them.

        THIS RETURNED THE FIRST ONE AND DROPPED THE REST, SILENTLY.
        "explain remdesivir for Zika and DENV-2" resolved both, answered
        about DENV-2 and said nothing about Zika -- no caveat, nothing in
        the trace, and not even the "(assumed)" marker, because a virus
        HAD been named, just not the one the reader might have meant.

        The winner was whichever name sorted first in a dict built by
        descending name LENGTH, so "denv-2" beat "zika" for no reason a
        user could predict. This module's own docstring states the rule
        it was breaking: "REPORTED, NOT SILENT ... dropping two of seven
        requested compounds from a ranking is a bounded harm; doing it
        without saying so is not." That rule was honoured for compounds
        and not for viruses, which decide the calibration, the caveats
        and the biology.

        Order is the order of the QUESTION, so a clarification lists them
        the way they were written rather than the way they were matched.
        """
        return [nid for nid in resolved.values()
                if self.tools.nodes.get(nid, {}).get("class")
                == "OrganismTaxon"]

    def _disease_taxa(self, resolved: dict[str, str]) -> list[str]:
        """Taxa reachable from a resolved Disease by its CAUSES edge.

        A question naming a disease rather than a virus -- "what is known
        about covid-19" -- names a virus as far as a reader is concerned.
        The graph records which: sixteen OrganismTaxon -CAUSES-> Disease
        edges. Derived from them rather than from a table here, so a
        disease added later needs no edit.
        """
        out = []
        for nid in resolved.values():
            if self.tools.nodes.get(nid, {}).get("class") != "Disease":
                continue
            for e in self.tools.inv.get(("CAUSES", nid), []):
                if e["subject"] not in out:
                    out.append(e["subject"])
        return out

    def _virus_in(self, resolved: dict[str, str]) -> str | None:
        """The single named virus, or None when there is not exactly one.

        None for "no virus" AND for "more than one", because the caller's
        job in both cases is to stop and ask rather than to pick.
        """
        taxa = self._viruses_in(resolved) or self._disease_taxa(resolved)
        return taxa[0] if len(taxa) == 1 else None

    def _compounds_in(self, resolved: dict[str, str]) -> list[str]:
        return [nid for nid in resolved.values()
                if self.tools.nodes.get(nid, {}).get("class") == "SmallMolecule"]

    # ------------------------------------------------------------ pipeline
    #: Questions about the agent itself rather than about the graph.
    #: Matched before entity resolution, because they name no entity and
    #: the resolver can only report that.
    #:
    #: THE TWO COMMONEST FIRST QUESTIONS BOTH FAILED. "hi introduce
    #: yourself" and "list all the viruses i can talk to you about" both
    #: returned "No entity in this question matches a node in the graph"
    #: -- the second one especially badly, because the agent lists those
    #: sixteen viruses in three other replies and simply would not do it
    #: when asked. Found by the owner using the frontend, 2026-11-21.
    ABOUT = ("introduce yourself", "who are you", "what are you",
             "what can you do", "what can i ask", "what do you do",
             "how do you work", "help", "hello", "hi there")
    LIST_VIRUSES = ("list all the viruses", "list the viruses",
                    "which viruses", "what viruses", "all the viruses",
                    "viruses you", "viruses do you", "viruses can i")

    def _about(self) -> str:
        held = sorted(self.tools._label(nid)
                      for nid, d in self.tools.nodes.items()
                      if d.get("class") == "OrganismTaxon")
        n_cmpd = sum(1 for d in self.tools.nodes.values()
                     if d.get("class") == "SmallMolecule")
        return (
            "I answer questions from a curated antiviral knowledge graph: "
            f"{len(self.tools.nodes):,} nodes covering {n_cmpd:,} small "
            f"molecules and {len(held)} viruses. Every answer is built "
            "from retrieved edges, each carrying the source it came "
            "from, and a written answer that cannot be checked against "
            "those edges is withheld rather than shown.\n\n"
            "What I can do:\n"
            "  - explain why a compound with known activity might work "
            "against a virus, route by route\n"
            "  - rank a set of compounds YOU name\n"
            "  - say which pathway a compound acts through\n"
            "  - retrieve provenanced evidence for one specific claim\n"
            "  - say what the graph holds about a compound, protein or "
            "virus\n\n"
            "What I refuse: generating a shortlist of NEW candidates. "
            "That is measured, not modest — the refusal quotes the "
            "measurement for the virus you name.\n\n"
            "Viruses I hold: " + ", ".join(held) + ".")

    def build_brief(self, question: str) -> Brief:
        trace: list[str] = []

        # ABOUT THE AGENT, NOT ABOUT THE GRAPH. Answered from the
        # release's own census, so the counts cannot drift from it.
        low = " ".join((question or "").lower().split())
        if any(p in low for p in self.LIST_VIRUSES):
            held = sorted(self.tools._label(nid)
                          for nid, d in self.tools.nodes.items()
                          if d.get("class") == "OrganismTaxon")
            return Brief(
                question=question, intent="about", allowed=True,
                clarification=(
                    f"This release holds {len(held)} viruses:\n\n  "
                    + "\n  ".join(held)
                    + "\n\nName one and ask about it. What has been "
                      "MEASURED differs sharply between them, and every "
                      "answer says which applies."),
                trace=trace + ["about_viruses"])
        if any(p == low or p in low for p in self.ABOUT):
            return Brief(
                question=question, intent="about", allowed=True,
                clarification=self._about(),
                trace=trace + ["about_agent"])

        intent = guards.classify_intent(question)
        trace.append(f"intent={intent.value}")

        resolved, problem, omitted = self._resolve_mentions(question)
        trace.append(f"resolved={len(resolved)}")
        if omitted:
            trace.append(f"omitted={len(omitted)}")
        # MORE THAN ONE VIRUS IS A QUESTION, NOT A CHOICE. Answering about
        # one of them and saying nothing is what this used to do.
        all_taxa = self._viruses_in(resolved) or self._disease_taxa(resolved)
        if len(all_taxa) > 1:
            labels = [self.tools._label(t) for t in all_taxa]
            trace.append("viruses=" + ",".join(labels))
            return Brief(
                question=question, intent=intent.value, allowed=True,
                clarification=(
                    f"This question names {len(labels)} viruses — "
                    + ", ".join(labels)
                    + ". The virus decides which routes exist, what has been "
                      "measured and which caveats apply, and they differ: ask "
                      "about one at a time, and the answers are not "
                      "interchangeable."),
                resolved=resolved,
                trace=trace + ["ambiguous_virus"])

        named = self._virus_in(resolved)

        # DOES THE ANSWER DEPEND ON THE VIRUS? Computed once, and used by
        # both the suggestion below and the "name one" path after it.
        #
        # Gating only the second was wrong: "what is the evidence that
        # remdesivir inhibits nsp12?" names no virus and does not need
        # one -- evidence(a, b) takes two node ids -- yet the suggester
        # fired and asked about TBEV. A clarification on a question the
        # agent could have answered is its own kind of wrong answer.
        compounds_named = [
            nid for nid in resolved.values()
            if self.tools.nodes.get(nid, {}).get("class") == "SmallMolecule"]
        # The EVIDENCE carve-out holds only when evidence() can actually
        # run, which needs two resolved nodes. With one it falls through
        # to `elif compounds` and calls explain(compound, virus) -- so
        # "what is the evidence that remdesivir inhibits nsp12?", where
        # nsp12 does not resolve, was answered about an ASSUMED virus
        # under an exemption granted on the grounds that it would not use
        # one.
        evidence_pair = (intent is guards.QueryKind.EVIDENCE
                         and len(resolved) >= 2)
        needs_virus = (intent in guards.ROUTED_KINDS
                       or intent is guards.QueryKind.PROSPECT
                       or (compounds_named and not evidence_pair))

        # NO VIRUS NAMED, BUT SOMETHING CLOSE TO ONE WAS TYPED. Ask.
        #
        # A typo used to reach the default silently-but-for-a-caveat:
        # "yellew fever virus" was answered about SARS-CoV-2. Suggesting
        # is safe where correcting is not -- every candidate within range
        # is listed and none is applied, because the four dengue
        # serotypes are one edit apart and their conclusions are
        # opposite.
        # OR NOTHING RESOLVED AT ALL. "what is known about zikaa?" is a
        # PROFILE question, so needs_virus is false -- but it resolved
        # nothing, and "did you mean Zika?" is a better reply than "no
        # entity in this question matches a node in the graph".
        if named is None and (needs_virus or not resolved):
            near = near_names(question, {
                k: v for k, v in self.tools._register_names.items()
                if v in self.tools.nodes})
            if near:
                # ONE LINE PER SET OF CANDIDATES, keyed by what is being
                # suggested rather than by what was typed: "yellew fever"
                # and "yellew fever virus" are the same suggestion twice,
                # and a reader counting lines would think two things were
                # wrong with their question.
                by_target: dict[tuple, list[str]] = {}
                for phrase, hits in sorted(near.items()):
                    labels = tuple(sorted({self.tools._label(
                        self.tools._register_names[h]) for h in hits}))
                    by_target.setdefault(labels, []).append(phrase)
                lines = []
                for labels, phrases in by_target.items():
                    typed = max(phrases, key=len)
                    lines.append(f"  \u201c{typed}\u201d — did you mean "
                                 + _or_list(labels) + "?")
                trace.append("suggested")
                return Brief(
                    question=question, intent=intent.value, allowed=True,
                    clarification=(
                        "No virus in this question matched one this release "
                        "holds, but something close did:\n\n"
                        + "\n".join(lines)
                        + "\n\nName it and ask again. Nothing is corrected "
                          "for you: the four dengue serotypes are one "
                          "character apart and their results are not "
                          "interchangeable, so a silent fix would change "
                          "the answer without showing you that it had."),
                    resolved=resolved,
                trace=trace + ["near_virus_name"])


        # NO VIRUS NAMED, AND THE VIRUS DECIDES THE ANSWER: ASK.
        #
        # The fallback to self.default_virus was defended by
        # guards.virus_assumed_caveat's docstring: "The default itself is
        # not wrong: SARS-CoV-2 is the only virus any protocol here
        # evaluated. Defaulting without saying so is." That was true when
        # written and is not true now -- fourteen ChEMBL cells, three
        # plate cells and ten prospective cells have been evaluated since,
        # six of the ten clearing chance. The premise that made a default
        # defensible is gone, and what remained was a caveat apologising
        # for a guess.
        #
        # THE CONDITION IS WHAT RETRIEVAL WILL ACTUALLY USE, not the
        # intent alone. ROUTED_KINDS and PROSPECT are virus-decided by
        # definition; so is any question that resolved a COMPOUND,
        # because retrieval then calls explain(compound, virus) or
        # triage(compounds, virus) whatever the intent classifier made of
        # the wording.
        #
        # Keying on intent alone missed exactly the case
        # virus_assumed_caveat's own docstring cites: "is chloroquine
        # useful?" classifies as UNKNOWN, falls to `elif compounds` in
        # retrieval, and was answered about SARS-CoV-2.
        #
        # PROFILE and EVIDENCE with no compound are about a node that may
        # not be a virus at all, so they keep the default and the caveat
        # that names it.
        # AND ONLY WHEN SOMETHING DID RESOLVE. With an empty `resolved`
        # the question named nothing the graph holds -- "explain the
        # mechanism of aspirin" -- and the useful reply names what kinds
        # of thing it does hold, which the `problem` path below already
        # does. Asking for a virus there answers a question the user did
        # not hit yet.
        # PROSPECT IS EXEMPT FROM THE "something resolved" REQUIREMENT.
        # "which drugs should we test?" is a COMPLETE question missing
        # only a virus -- there is no entity for it to have named -- and
        # the refusal it would otherwise get is virus-specific: the AUCs
        # it quotes, and whether the closing offer about explaining a
        # mechanism is even true for that virus. Refusing in detail about
        # SARS-CoV-2 is the case guards' own comment calls the one where
        # "the substitution does the most damage".
        if named is None and needs_virus and (
                resolved or intent is guards.QueryKind.PROSPECT):
            held = sorted(self.tools._label(nid)
                          for nid, d in self.tools.nodes.items()
                          if d.get("class") == "OrganismTaxon")
            trace.append("no_virus_named")
            return Brief(
                question=question, intent=intent.value, allowed=True,
                clarification=(
                    "This question does not name a virus, and for a question "
                    "of this kind the virus decides the answer: which routes "
                    "exist, what has been measured, which caveats apply, and "
                    "what a refusal is entitled to say.\n\n"
                    "This release holds " + ", ".join(held) + ".\n\n"
                    "Nothing is assumed for you. Until 2026-11-20 this "
                    "answered about SARS-CoV-2 under a caveat, on the "
                    "grounds that it was the only virus any protocol here "
                    "had evaluated — which stopped being true."),
                resolved=resolved,
                trace=trace + ["virus_not_named"])

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
        # A COMPOUND QUESTION THAT NAMED NO COMPOUND THIS GRAPH HOLDS.
        #
        # "why might frobnicatinib work against SARS-CoV-2?" resolved the
        # virus, found no compound, fell through to the final `else` in
        # retrieval and returned a PROFILE OF THE VIRUS -- edge counts --
        # with nothing anywhere saying the compound had not been found. A
        # reader takes that as an answer about frobnicatinib. It is an
        # answer about SARS-CoV-2.
        #
        # EXPLAIN and TRIAGE are the two intents that are ABOUT a
        # compound; PROFILE and EVIDENCE legitimately take any node.
        # AFTER the "nothing resolved at all" and "default virus absent"
        # paths, and only when something DID resolve. Placed before them
        # at first, it pre-empted both: a question resolving nothing was
        # told "no compound matched" instead of "no entity matched", and
        # "rank these compounds" against an absent default virus lost the
        # clarification that names the viruses the release holds. Caught
        # by two existing tests that had the orderings right.
        if resolved and intent in (guards.QueryKind.EXPLAIN,
                                   guards.QueryKind.TRIAGE):
            if not [nid for nid in resolved.values()
                    if self.tools.nodes.get(nid, {}).get("class")
                    == "SmallMolecule"]:
                trace.append("no_compound")
                return Brief(
                    question=question, intent=intent.value, allowed=True,
                    clarification=(
                        "No compound in this question matches one this "
                        "release holds, and the question asks about a "
                        "compound. Nothing is substituted for it: answering "
                        "about the virus instead would read as an answer "
                        "about the compound you named.\n\n"
                        "Name it as the graph does — a generic name or a "
                        "ChEMBL id — or ask what the graph stores about the "
                        "virus directly."),
                    resolved=resolved,
                trace=trace + ["no_compound_resolved"])

        verdict = guards.check(
            intent, virus, calibration=self.tools.calibration_for(virus),
            temporal=self.tools.temporal_for(virus),
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
