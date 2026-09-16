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

EDGE_ID_RE = re.compile(r"\[([A-Za-z]?E\d+(?:\s*,\s*[A-Za-z]?E\d+)*)\]")
NUMBER_RE = re.compile(r"\d+(?:[.,]\d+)*")
SENTENCE_RE = re.compile(r"(?<=[.!?])\s+")

# Sentences that make no factual claim and need no citation.
HEDGE_STARTS = (
    "this", "these", "however", "note", "in summary", "overall", "the graph",
    "no ", "there is no", "i cannot", "it is not", "caution", "warning",
    # Restating a required caveat is not an uncited claim: the caveat is
    # supplied by the brief, not asserted by the model.
    "selective", "cytotoxic", "selectivity unknown", "cell context",
    "cell line", "host-directed route", "the cell line",
)


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


def _numbers(text: str) -> set[str]:
    """Numeric tokens, normalised so 1,560 and 1560 compare equal.

    Years and small ordinals are excluded: "the 2021 study" and "Route 1" are
    not claims about magnitude.
    """
    out = set()
    for raw in NUMBER_RE.findall(text or ""):
        norm = raw.replace(",", "")
        if norm.isdigit() and (len(norm) == 4 and norm.startswith(("19", "20"))):
            continue          # a year
        if norm.isdigit() and int(norm) <= 10:
            continue          # an ordinal or list index
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
    uncited: list[str] = []
    if not citable_edge_ids:
        return Verification(
            passed=not any(f.severity == "fail" for f in findings),
            findings=findings, cited=sorted(set(cited)))
    for s in _sentences(answer):
        if EDGE_ID_RE.search(s):
            continue
        low = s.lower()
        if any(low.startswith(h) for h in HEDGE_STARTS):
            continue
        if len(s.split()) < 6:
            continue
        uncited.append(s)
    if uncited:
        findings.append(Finding(
            "warn", "uncited_claim",
            f"{len(uncited)} sentence(s) make a claim without a citation"))

    # 3. Mandatory warnings must survive into the answer. This is the check
    #    that stops a cell-context caveat being dropped for brevity.
    " ".join(facts).lower()
    ans_low = (answer or "").lower()
    for w in required_warnings:
        head = w.split(":")[0].strip().lower()
        key = head if head and len(head) < 40 else w[:40].lower()
        if key and key not in ans_low:
            findings.append(Finding(
                "fail", "warning_dropped",
                f"required warning not present in the answer: {w[:70]}"))

    # 4. Every number in the answer must appear in the facts. A transposed
    #    digit in a potency value is harder to notice than a wrong sentence.
    fact_numbers = _numbers(" ".join(facts) + " " + " ".join(required_warnings))
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
