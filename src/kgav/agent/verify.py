"""Verification: does the answer say only what the facts support?

WHY A SEPARATE STAGE. Self-critique within one context is weak, because the
model is checking work it has already committed to. This stage runs with a
fresh context containing only the claims and the evidence, and — more
importantly — most of it is not a model at all.

WHAT IS DETERMINISTIC HERE, and it is nearly everything:

  citation validity   every [Exxx] resolves to an edge that was retrieved
  citation coverage   every factual sentence carries one
  warning survival    every mandatory warning appears in the answer
  number provenance   every figure in the answer appears in the facts
  refusal integrity   a refusal is not accompanied by an answer

Only the optional `model_check` is generative, and it is advisory. The
deterministic checks decide.

WHY NUMBER PROVENANCE MATTERS MOST. A language model narrating "IC50 75 nM"
may render it as "approximately 100 nM" or transpose a digit, and a plausible
wrong number is harder to catch than a plausible wrong sentence. Every numeric
token in the answer is therefore matched against the retrieved facts.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

from kgav.agent.guards import CAVEAT_LABELS

# Hex, not just digits: edge ids are now blake2b digests of the fact key
# (assemble.fact_id). [0-9a-f]+ still matches the old E<digits> line-index
# form, so a saved answer from before the change still parses.
EDGE_ID_RE = re.compile(
    r"\[([A-Za-z]?E[0-9a-f]+(?:\s*,\s*[A-Za-z]?E[0-9a-f]+)*)\]")

#: The id inside a FACT's provenance tail, which EDGE_ID_RE cannot read.
#: tools._verbalise renders a citation as "[chembl, PMID:1, 2021; E335139]",
#: and EDGE_ID_RE matches only a bracket whose whole content is ids -- the
#: form a MODEL writes. The first version of shown_ids below used
#: EDGE_ID_RE on the facts, found nothing in any real fact, and stripped
#: the exemption from every correctly cited answer. Three pre-existing
#: tests failed on it, which is the only reason it is not in the release.
FACT_ID_RE = re.compile(r"\[[^\]]*;\s*([A-Za-z]?E[0-9a-f]+)\s*\]")
NUMBER_RE = re.compile(r"\d+(?:[.,]\d+)*")

# A NUMBER CARRYING A UNIT IS A MAGNITUDE CLAIM, whatever its size. The two
# exclusions below it -- years, and integers at or under ten -- exist to stop
# "the 2021 study" and "Route 1" being read as claims about quantity, and they
# were applied to EVERY number. So a fabricated "Ki of 3 nM" verified clean
# against facts reading 190 nM, and so would "1960 nM", dropped as a year.
# Measured: claimed Ki of 3 and 9 nM both passed; 11 and 750 were caught.
#
# Single-digit potencies are not hypothetical here. nirmatrelvir's own Ki on
# v0.1 is 3 nM, so the exempted range is exactly where the most potent
# compounds live -- and a 100-fold error in a potency is the kind of plausible
# wrong number this module's docstring says is hardest to catch.
# THE SEPARATOR MAY BE A HYPHEN. "19-fold" and "19-day" are unit-bearing
# and were read as a bare 19, which matched any 19 in the facts -- so a
# potency reread as a ratio verified clean. Added 2026-11-02.
# LONGEST ALTERNATIVE FIRST. Python's alternation is ordered, so with
# `nm` ahead of `nmol` the regex captured "nmol" as "nm" and left "ol"
# dangling -- which happened to be harmless here, and meant three
# entries in _UNIT_ALIASES could never fire. Dead code that looks like
# it is normalising something is worse than no normalising. Found by a
# surviving mutant that deleted the alias map with no test noticing.
UNIT_RE = re.compile(
    r"[\s-]*(?:nmol|µmol|umol|molar|percent|hours?|hrs?|days?|kda|fold"
    r"|nm|µm|um|mm|pm|mg|µg|ug|kg|ml|µl|ul|da\b|min"
    r"|m\b|l\b|h\b|x\b|%|×)",
    re.IGNORECASE)
SENTENCE_RE = re.compile(r"(?<=[.!?])\s+")

# Openings that carry no claim of their own. A sentence is judged on
# what FOLLOWS one of these, not excused by starting with one.
#
# THEY USED TO EXEMPT THE WHOLE SENTENCE, which made every one of these
# verify clean against a fact about an IC50:
#
#   "Note that nirmatrelvir cures every coronavirus infection in humans."
#   "This compound cures every coronavirus infection in humans everywhere."
#
# Found by external audit, 2026-11-02.
HEDGE_STARTS = (
    "this", "these", "however", "note", "in summary", "overall", "the graph",
    "no ", "there is no", "i cannot", "it is not", "caution", "warning",
    "the cell line", "note that", "that",
)

# THE CAVEAT LABELS ARE NOT LISTED HERE ANY MORE. They were, and the
# audit's own example is what that cost:
#
#   "CALIBRATION: scoring for this virus is fully validated and the
#    ranking can be trusted."
#
# -- a fabrication that inverts the one thing the calibration caveat
# exists to say, exempted because it began with the label.
#
# _is_supplied already covers a restated caveat, by CONTENT: its
# docstring says in so many words that it "replaces prefix-matching a
# caveat's label", and the labels stayed behind anyway. A Caveat is a
# string carrying its own label, and _supplied_corpus holds the
# warnings verbatim, so a genuine restatement is still exempt and a
# claim wearing the label is not.
_CAVEAT_LABEL_RE = re.compile(
    r"^\s*(?:" + "|".join(re.escape(v) for v in CAVEAT_LABELS.values())
    + r")\s*:\s*", re.IGNORECASE)


def _strip_openings(sentence: str) -> str:
    """The claim a sentence makes, with its hedges and labels removed.

    Applied repeatedly and boundedly: "CELL CONTEXT: this was measured
    in Vero cells" carries two, and one pass would leave the other
    doing the exempting.
    """
    cur = _CAVEAT_LABEL_RE.sub("", sentence or "").strip()
    for _ in range(4):
        low = cur.lower()
        for h in sorted(HEDGE_STARTS, key=len, reverse=True):
            if low.startswith(h):
                cur = cur[len(h):].lstrip(" ,:;-\u2014")
                break
        else:
            break
    return cur


@dataclass
class Finding:
    severity: str        # "fail" | "warn"
    code: str
    detail: str


@dataclass
class Verification:
    passed: bool
    findings: list[Finding] = field(default_factory=list)
    cited: list[str] = field(default_factory=list)
    uncited_sentences: list[str] = field(default_factory=list)

    @property
    def failures(self) -> list[Finding]:
        return [f for f in self.findings if f.severity == "fail"]

    def report(self) -> str:
        if self.passed and not self.findings:
            return "verified: every claim is cited and grounded"
        lines = ["verification failed" if not self.passed else "verified with warnings"]
        lines += [f"  [{f.severity}] {f.code}: {f.detail}" for f in self.findings]
        return "\n".join(lines)


def _sentences(text: str) -> list[str]:
    """Split into sentences, after neutralising markdown.

    Models format. Numbered list markers ("1.", "2.") are split as sentence
    ends by any regex keyed on a period, and bold markers put "**" at the
    start of a line so a hedge-prefix check never matches. Both produced
    uncited-claim warnings on text that made no uncited claim.
    """
    t = text or ""
    t = re.sub(r"^\s*\d+[.)]\s+", "", t, flags=re.MULTILINE)   # list markers
    t = re.sub(r"^\s*[-*•]\s+", "", t, flags=re.MULTILINE)      # bullets
    t = t.replace("**", "").replace("__", "")           # emphasis
    return [s.strip() for s in SENTENCE_RE.split(t) if s.strip()]


#: Punctuation that differs between a caveat as supplied and the same
#: words as a model restates them: it ends a borrowed clause with a full
#: stop where the brief had a comma. Dropped on BOTH sides, so
#: containment compares words.
_PUNCT_RE = re.compile(r"[.,;:!?()\[\]\u2014\u2013-]+")


def _normalise(text: str, drop_punct: bool = False) -> str:
    """Lowercase, markdown-free, citation-free, single-spaced.

    `drop_punct` is for CONTAINMENT only. A model that restates half a
    caveat ends it with a full stop where the brief had a comma --
    "... measured in Vero cells." against "... measured in Vero cells,
    which lack TMPRSS2." -- and the substring test then failed on the
    punctuation alone. Harmless while an uncited claim was a warning
    and the caveat labels exempted the sentence anyway; both of those
    changed on 2026-11-02, and this surfaced immediately.
    """
    s = EDGE_ID_RE.sub(" ", text or "")
    s = s.replace("**", "").replace("__", "").replace("*", "")
    if drop_punct:
        s = _PUNCT_RE.sub(" ", s)
    return " ".join(s.lower().split())


def _supplied_corpus(facts: list, warnings: list) -> str:
    """Everything the brief handed the model, normalised for containment."""
    return _normalise(" ".join([str(x) for x in (facts or [])]
                               + [str(x) for x in (warnings or [])]),
                      drop_punct=True)


def _is_supplied(sentence: str, corpus: str, min_words: int = 6) -> bool:
    """Whether this sentence was COPIED from the brief rather than asserted.

    Replaces prefix-matching a caveat's label, which only ever exempted a
    caveat's FIRST sentence. Each caveat is several sentences long -- the
    cell-context one is three -- so its continuations were reported as
    uncited claims, on text the brief itself supplied. Forcing the caveats in
    verbatim made that worse: every answer then carried the warning.

    Containment rather than equality, because a model legitimately merges a
    caveat sentence into surrounding prose. Short sentences are excluded from
    the shortcut: "it does not transfer." appears inside almost any corpus by
    accident, and exempting it would let a real uncited claim through.
    """
    norm = _normalise(sentence, drop_punct=True)
    if len(norm.split()) < min_words:
        return False
    return norm in corpus


def _caveat_marker(warning) -> str:
    """What must appear in an answer for this caveat to count as surviving.

    A guards.Caveat carries an explicit label, so the answer must keep the
    LABEL and may paraphrase everything after it however it likes. That is the
    whole point: enforcement is tied to the caveat's identity, not its prose,
    so rewriting the wording cannot silently change what is required.

    A plain string -- one not built by guards.caveat(), or one that has been
    through JSON and lost its attributes -- falls back to the older
    prose-derived key. Enforcing something imperfect beats enforcing nothing,
    and for a labelled caveat the fallback recovers the same label anyway,
    since the label is the text before the first colon.
    """
    label = getattr(warning, "label", None)
    if label:
        return label.lower()
    head = warning.split(":")[0].strip().lower()
    return head if head and len(head) < 40 else warning[:40].lower()


def caveat_labels(warnings: list) -> list[str]:
    """The LABEL of each caveat, uppercased -- what an answer must keep.

    Public because the synthesis stage names these back to the model on a
    retry, and reaching across modules for _caveat_marker made that an
    accident waiting to be renamed.
    """
    out = []
    for w in warnings or []:
        marker = _caveat_marker(w)
        if marker:
            out.append(marker.rstrip(": ").upper())
    return out


def dropped_warnings(answer: str, required: list) -> list:
    """Required caveats whose label is missing from `answer`.

    Extracted so the synthesis stage can ask the same question the verifier
    will, and retry before an answer is thrown away. A 7B model drops a
    caveat for brevity reliably enough that one withheld answer per
    host-directed question was the normal outcome.
    """
    low = (answer or "").lower()
    out = []
    for w in required or []:
        marker = _caveat_marker(w)
        if marker and marker not in low:
            out.append(w)
    return out


def _caveat_name(warning) -> str:
    """How a dropped caveat is named in the report."""
    return getattr(warning, "key", None) or str(warning)[:70]


#: Spellings that mean the same magnitude, so "19 nM" and "19 nm" and
#: "19 nmol" compare equal and "19 uM" does not compare equal to any of
#: them. Anything not listed keeps its own lowercased spelling.
_UNIT_ALIASES = {
    "nmol": "nm", "\u00b5m": "um", "umol": "um", "\u00b5mol": "um",
    "molar": "m", "percent": "%", "x": "fold", "\u00d7": "fold",
    "\u00b5g": "ug", "\u00b5l": "ul", "hrs": "h", "hr": "h",
    "hours": "h", "hour": "h", "days": "day", "da": "da",
}

#: A fact may carry the unit in the KEY rather than after the number --
#: `ic50_nm=19`, `ki_nm: 3`. Matched so the facts side can pair those too;
#: without it every unit-bearing answer would read as ungrounded.
_KEYED_UNIT_RE = re.compile(
    r"[a-z0-9_]*_(nm|um|mm|pm|molar|pct|percent|fold|mg|ug|kg|ml|ul|kda|da)"
    r"\s*[=:]\s*$", re.IGNORECASE)


def _unit_token(raw: str) -> str:
    """The unit, stripped of its separator and normalised.

    The separator is part of the match -- UNIT_RE starts `[\\s-]*` so that
    "19-fold" is unit-bearing -- and `.strip()` alone leaves the hyphen
    on, so "19-fold" and "19 fold" did not compare equal. Introduced by
    the hyphen change in the same commit that needed it.
    """
    u = (raw or "").strip(" \t\r\n-").lower()
    return _UNIT_ALIASES.get(u, u)


def _numbers(text: str, exclude_ordinals: bool = True,
             facts_side: bool = False) -> set[str]:
    """Numeric tokens, normalised so 1,560 and 1560 compare equal.

    Years and small ordinals are excluded: "the 2021 study" and "Route 1" are
    not claims about magnitude. A number carrying a unit is always kept.

    EXCLUDE_ORDINALS=FALSE FOR THE FACTS SIDE. The check is "every number in
    the ANSWER appears in the facts", so filtering the facts can only invent
    findings. It was symmetric while both sides dropped the same tokens; once
    the answer side became unit-aware it was not, and an answer reading
    "5 nM" against a fact reading "ic50_nm=5" was reported as an ungrounded
    figure. The facts are the ground truth here and are taken whole.

    A NUMBER IS PAIRED WITH ITS UNIT, which it was not until 2026-11-02.
    The set held bare numeric tokens, so any number in the facts could be
    reattached to any claim with any unit. Against a fact reading
    `IC50 19 nM`, every one of these verified clean:

        "IC50 19 uM"              a 1000-fold error
        "19 mM"                   a million-fold error
        "inhibits nsp5 by 19%"    a potency reread as an efficacy
        "inhibits nsp5 19-fold"   the same, as a ratio

    The earlier UNIT_RE work made the exclusions unit-aware -- it stopped
    "3 nM" being dropped as a small ordinal -- and never checked the
    PAIRING. Found by external audit.

    `facts_side` emits the bare number as well as the pair, because a fact
    is ground truth for both a united and an un-united mention of it. The
    answer side emits the pair alone, which is what makes a wrong unit a
    finding.
    """
    out = set()
    t = text or ""
    for m in NUMBER_RE.finditer(t):
        norm = m.group(0).replace(",", "")
        unit_m = UNIT_RE.match(t, m.end())
        unit = _unit_token(unit_m.group(0)) if unit_m else ""
        if not unit and facts_side:
            keyed = _KEYED_UNIT_RE.search(t, 0, m.start())
            if keyed:
                unit = _unit_token(keyed.group(1))
        if exclude_ordinals and not unit:
            if norm.isdigit() and len(norm) == 4 and norm.startswith(("19", "20")):
                continue      # a year
            if norm.isdigit() and int(norm) <= 10:
                continue      # an ordinal or list index
        out.add(f"{norm} {unit}" if unit else norm)
        if facts_side and unit:
            # A fact grounds a mention of the same number without a unit.
            out.add(norm)
    return out


def verify(answer: str, facts: list[str], citable_edge_ids: list[str],
           required_warnings: list[str] | None = None,
           allowed: bool = True, refusal: str = "",
           model_check=None) -> Verification:
    """Check an answer against the brief it was generated from."""
    findings: list[Finding] = []
    required_warnings = required_warnings or []

    # A refusal must stand alone. An answer accompanying one means the
    # pipeline produced the very output the guard exists to prevent.
    if not allowed:
        if answer and answer.strip() != (refusal or "").strip():
            findings.append(Finding(
                "fail", "refusal_contaminated",
                "a refusal was issued but the output contains other content"))
        return Verification(passed=not findings, findings=findings)

    cited: list[str] = []
    for group in EDGE_ID_RE.findall(answer or ""):
        cited += [c.strip() for c in group.split(",")]

    # 1. Every citation must resolve to an edge that was actually retrieved.
    allowed_ids = set(citable_edge_ids)
    invented = [c for c in cited if c not in allowed_ids]
    if invented:
        findings.append(Finding(
            "fail", "invented_citation",
            f"cites edge ids that were not retrieved: {', '.join(sorted(set(invented)))}"))

    # 2. Every factual sentence must carry a citation -- UNLESS nothing was
    #    citable. Triage and profile return summary counts rather than
    #    individual edges, and demanding citations where none exist is what
    #    drove the model to invent them. The invented-citation check above
    #    still applies: with an empty allowed-list, ANY citation is invented.
    #
    #    ONLY THIS CHECK IS SKIPPED. An earlier version RETURNED here, which
    #    silently disabled the warning and number checks below for every
    #    triage and profile query -- the two shapes that never carry edge ids.
    #    A fabricated potency and a dropped cell-context caveat both passed,
    #    reported as "verified: every claim is cited and grounded". Having no
    #    citable edges means citations cannot be checked; it says nothing
    #    about whether the numbers are real or the caveats survived.
    #    AND A CITATION ONLY EXEMPTS A SENTENCE IF ITS FACT WAS SHOWN.
    #    `if EDGE_ID_RE.search(s): continue` exempted any sentence
    #    carrying ANY citable id, and a tool may hand over far more ids
    #    than it shows facts for: pathways() returns 29 census lines and
    #    1,500 PARTICIPATES_IN ids, none of which appears in a line. So
    #
    #      "Promiscuousdrug eradicates SARS-CoV-2 in every patient
    #       tested and is the standard of care worldwide [E1001]."
    #
    #    verified CLEAN with every caveat kept: E1001 is citable, so the
    #    citation is not invented, and the citation exempted the sentence
    #    from this check. That is the hole the uncited_claim promotion
    #    was meant to close, reached through the ids instead of the
    #    hedges.
    #
    #    Those ids are provenance for a COUNTED claim, not support for a
    #    sentence. They stay citable -- the count they back is real -- and
    #    they no longer buy an exemption.
    joined = " ".join(facts or [])
    shown_ids: set[str] = set(FACT_ID_RE.findall(joined))
    for group in EDGE_ID_RE.findall(joined):
        shown_ids |= {c.strip() for c in group.split(",")}

    uncited: list[str] = []
    if citable_edge_ids:
        supplied = _supplied_corpus(facts, required_warnings)
        for s in _sentences(answer):
            here = set()
            for group in EDGE_ID_RE.findall(s):
                here |= {c.strip() for c in group.split(",")}
            if here & shown_ids:
                continue
            # Copied from the brief, so it is grounded by construction --
            # the model did not assert it, the retrieval did. Tested
            # FIRST and on the whole sentence, because the corpus holds
            # each caveat with its label attached.
            if _is_supplied(s, supplied):
                continue
            core = _strip_openings(s)
            if _is_supplied(core, supplied):
                continue
            # The LENGTH TEST APPLIES TO THE CLAIM, not to the sentence
            # carrying it: "However, no." is a remark and "Note that
            # nirmatrelvir cures every coronavirus infection" is not.
            if len(core.split()) < 6:
                continue
            uncited.append(s)
        if uncited:
            findings.append(Finding(
                # FAIL, NOT WARN. As a warning it never touched `passed`,
                # so an answer carrying a fabricated claim with no
                # citation was reported as verified. The whole purpose of
                # this module is that an unsupported sentence does not
                # reach a reader.
                "fail", "uncited_claim",
                f"{len(uncited)} sentence(s) make a claim without a "
                f"citation: " + " | ".join(s[:60] for s in uncited[:3])))

    # 3. Mandatory warnings must survive into the answer. This is the check
    #    that stops a cell-context caveat being dropped for brevity.
    for w in dropped_warnings(answer, required_warnings):
        findings.append(Finding(
            "fail", "warning_dropped",
            f"required warning not present in the answer: {_caveat_name(w)}"))

    # 3b. A TRUNCATED ANSWER IS NOT AN ANSWER. num_predict cut the chloroquine
    #     answer at "measured in Vero E", mid-word and mid-caveat. The three
    #     checks above cannot see it: the caveat's label had already appeared,
    #     every citation so far resolved, and a half-written number passes a
    #     containment test. So a cut-off answer could lose the second half of
    #     a caveat, or report "160" from "1,600", and verify clean. Demanding
    #     terminal punctuation is crude and catches exactly this.
    stripped = (answer or "").rstrip()
    if stripped and stripped[-1] not in ".!?\"')]:":
        findings.append(Finding(
            "fail", "truncated_answer",
            f"the answer ends mid-sentence ({stripped[-40:]!r}), so content "
            f"after the cut -- possibly the rest of a required caveat -- is "
            f"missing"))

    # 4. Every number in the answer must appear in the facts. A transposed
    #    digit in a potency value is harder to notice than a wrong sentence.
    fact_numbers = _numbers(" ".join(facts) + " " + " ".join(required_warnings),
                            exclude_ordinals=False, facts_side=True)
    # Strip citation brackets first: the digits inside [E999999] are already
    # reported by the citation check, and repeating them as a numeric finding
    # adds noise to a report a human has to read.
    answer_numbers = _numbers(EDGE_ID_RE.sub(" ", answer or ""))
    stray = sorted(answer_numbers - fact_numbers)
    if stray:
        findings.append(Finding(
            "fail", "ungrounded_number",
            f"figures not present in the retrieved facts: {', '.join(stray[:6])}"))

    # 5. Optional model-based check, advisory only.
    if model_check is not None:
        verdict = model_check(answer, facts)
        if verdict and not verdict.get("supported", True):
            findings.append(Finding(
                "warn", "model_dispute",
                str(verdict.get("reason", "second model disputes the answer"))))

    return Verification(passed=not any(f.severity == "fail" for f in findings),
                        findings=findings, cited=sorted(set(cited)),
                        uncited_sentences=uncited)


def verify_response(response: dict, answer: str | None = None,
                    model_check=None) -> Verification:
    """Verify an Orchestrator.answer() result."""
    text = answer if answer is not None else (response.get("answer") or "")
    return verify(text,
                  facts=response.get("facts") or [],
                  citable_edge_ids=response.get("citable_edge_ids") or [],
                  required_warnings=response.get("warnings") or [],
                  allowed=response.get("allowed", True),
                  refusal=response.get("answer") or "",
                  model_check=model_check)


CHECK_PROMPT = """You are checking whether a summary is supported by evidence.

You will be given FACTS and a SUMMARY. Answer only whether every claim in the
summary follows from the facts. Do not use your own knowledge; a claim that is
true in general but absent from the facts is NOT supported.

Reply with JSON only: {"supported": true|false, "reason": "..."}

FACTS:
%s

SUMMARY:
%s
"""
