#!/usr/bin/env python3
"""kgav ask — query the graph from the command line.

    kgav_ask.py "why might nirmatrelvir work against SARS-CoV-2?"
    kgav_ask.py --no-model "rank nirmatrelvir and remdesivir for SARS-CoV-2"
    kgav_ask.py --json "what is known about SIGMAR1"

THE ANSWER IS PRINTED ONLY IF IT VERIFIES. A tool that shows a hallucinated
answer with a warning underneath is a tool that shows hallucinated answers: the
caveat is read second, if at all, and the fluent prose is what sticks. On
failure the verification report is printed and the brief is offered instead.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.agent.llm import DEFAULT_MODEL, ModelUnavailable, OllamaModel
from kgav.agent.orchestrator import Orchestrator
from kgav.agent.tools import GraphTools
from kgav.agent.verify import verify_response

RULE = "─" * 74


def render(out: dict, verification=None, brief_only: bool = False) -> str:
    lines: list[str] = []

    if not out["allowed"]:
        lines += ["REFUSED", RULE, out["answer"]]
        return "\n".join(lines)

    if out.get("facts"):
        lines += ["EVIDENCE", RULE]
        lines += [f"  {f}" for f in out["facts"]]
        lines.append("")

    if out.get("warnings"):
        lines += ["CAVEATS", RULE]
        for w in out["warnings"]:
            lines += [f"  • {w}", ""]

    if brief_only or out.get("answer") is None:
        lines += [("(no model — evidence shown above; run without --no-model "
                  "for a synthesised answer)")]
        return "\n".join(lines)

    if verification is not None and not verification.passed:
        lines += ["ANSWER WITHHELD — FAILED VERIFICATION", RULE]
        for f in verification.findings:
            lines.append(f"  [{f.severity}] {f.code}: {f.detail}")
        lines += ["",
                  ("  The model produced text that is not supported by the "
                  "retrieved evidence."),
                  "  The evidence above is what the graph actually contains."]
        return "\n".join(lines)

    lines += ["ANSWER", RULE, out["answer"]]
    if verification is not None and verification.findings:
        lines += ["", "verification notes:"]
        lines += [f"  [{f.severity}] {f.code}: {f.detail}"
                  for f in verification.findings]
    return "\n".join(lines)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description="Query the kgav graph.")
    ap.add_argument("question", nargs="+")
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--host", default="http://localhost:11434")
    ap.add_argument("--no-model", action="store_true",
                    help="show retrieved evidence without synthesis")
    ap.add_argument("--json", action="store_true",
                    help="machine-readable output including the verification")
    ap.add_argument("--show-prompt", action="store_true")
    args = ap.parse_args()
    question = " ".join(args.question)

    if not (args.release / "nodes.jsonl").exists():
        print(f"release not found at {args.release}", file=sys.stderr)
        return 2

    orch = Orchestrator(GraphTools(args.release))

    model = None
    if not args.no_model:
        m = OllamaModel(model=args.model, host=args.host)
        ok, note = m.available()
        if ok:
            model = m
        else:
            # Degrade to the brief rather than failing: the evidence is useful
            # without prose, and a missing model is not a missing answer.
            print(f"note: {note} — showing evidence without synthesis\n",
                  file=sys.stderr)

    try:
        out = orch.answer(question, model)
    except ModelUnavailable as e:
        print(f"note: {e} — showing evidence without synthesis\n",
              file=sys.stderr)
        out = orch.answer(question, None)

    if args.show_prompt and out.get("prompt"):
        print("PROMPT\n" + RULE)
        print(out["prompt"])
        print()

    verification = None
    if out["allowed"] and out.get("answer"):
        verification = verify_response(out)

    if args.json:
        print(json.dumps({
            **{k: v for k, v in out.items() if k != "prompt"},
            "verification": None if verification is None else {
                "passed": verification.passed,
                "findings": [{"severity": f.severity, "code": f.code,
                              "detail": f.detail}
                             for f in verification.findings],
                "cited": verification.cited},
        }, indent=2))
        return 0 if (verification is None or verification.passed) else 1

    print(render(out, verification, brief_only=args.no_model))
    return 0 if (verification is None or verification.passed) else 1


if __name__ == "__main__":
    raise SystemExit(main())
