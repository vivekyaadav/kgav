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
from kgav import calibration as kgav_calibration
from kgav.labels import measured_reach
from kgav import landscape as kgav_landscape
from kgav.pipeline import release_dir

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


def json_payload(out: dict, verification=None) -> dict:
    """The --json surface, as a value, so it can be compared with `render`.

    THE TWO SURFACES DISAGREED ABOUT WHETHER AN UNVERIFIED ANSWER MAY BE
    SHOWN (fix-list A6). `render` withholds it outright -- "ANSWER WITHHELD
    -- FAILED VERIFICATION", the findings, and a line saying the evidence
    is what the graph actually contains -- while the --json branch printed
    `out` whole, answer included. The exit code was 1 either way, which
    protects a shell `&&` and nothing that reads the JSON: `jq -r .answer`
    returned model text that had failed verification, with the failure in a
    sibling key.

    The answer is withheld here too, on the same condition. The findings
    stay -- they are the diagnostic value of the failure -- and so do the
    facts and warnings, because those are the graph's and not the model's.

    Extracted to module level for the same reason `render` is: the test
    that matters asks BOTH for the same case and fails if they disagree,
    and it could not reach this code while it lived inside main().

    PROMPTS ARE EXCLUDED BY SHAPE, NOT BY NAME. This read
    `if k != "prompt"`, naming one key, and the orchestrator sets two:
    `prompt` and, when the first answer dropped a required warning,
    `retry_prompt`. So every --json reply to a question that needed a
    retry carried ~2,900 characters of prompt -- the system rules, the
    whole brief, and the "YOUR PREVIOUS ANSWER WAS REJECTED" preamble --
    in a surface whose author had deliberately removed the other one.
    `--show-prompt` is the opt-in for seeing a prompt, and it prints
    `out["prompt"]` only, so the retry prompt was readable in the
    machine surface and nowhere else.

    This is the defect class the repository already has a name for: a
    filter that lists what it excludes goes stale the moment the thing
    it filters grows a sibling. Matching `*prompt*` means the next one
    is covered on the day it is added, and
    tests/test_agent_json_surface.py fails if a prompt-shaped key ever
    reaches a caller again.
    """
    payload = {k: v for k, v in out.items() if "prompt" not in k}
    if verification is not None and not verification.passed:
        payload["answer"] = None
        payload["answer_withheld"] = (
            "failed verification; see verification.findings. The model "
            "produced text that is not supported by the retrieved "
            "evidence, and `facts` is what the graph actually contains.")
    payload["verification"] = None if verification is None else {
        "passed": verification.passed,
        "findings": [{"severity": f.severity, "code": f.code,
                      "detail": f.detail}
                     for f in verification.findings],
        "cited": verification.cited}
    return payload


def build_parser(root: Path) -> argparse.ArgumentParser:
    """The CLI surface, separate from main() so the DEFAULTS are testable.

    Built inside main() until 2026-11-07, which is why --results could sit
    for weeks pointing a directory-shaped argument at a file: every unit
    test of the consumer passed it a directory, because the author of the
    test knew what it wanted. Nothing could assert what the CLI passes
    without running the whole agent. tests/test_agent_cli_wiring.py now
    does, against these defaults rather than against a copy of them.
    """
    ap = argparse.ArgumentParser(description="Query the kgav graph.")
    ap.add_argument("question", nargs="+")
    ap.add_argument("--release", type=Path, default=release_dir(root))
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--host", default="http://localhost:11434")
    ap.add_argument("--no-model", action="store_true",
                    help="show retrieved evidence without synthesis")
    ap.add_argument("--json", action="store_true",
                    help="machine-readable output including the verification")
    ap.add_argument("--show-prompt", action="store_true")
    ap.add_argument("--results", type=Path,
                    default=root / "data/results/hard_negatives.json",
                    help="evaluation run whose measured AUCs the caveats "
                         "quote. Absent, a host-directed route is reported as "
                         "uncalibrated rather than given a number.")
    ap.add_argument("--temporal-results", type=Path,
                    default=root / "data/results/"
                                   "temporal_structure_scorer.json",
                    help="the PROSPECTIVE run whose AUCs the "
                         "candidate-generation refusal quotes. A separate "
                         "file from --results because it is a separate "
                         "study with its own registration and its own "
                         "shape; --results' `temporal` section is empty on "
                         "this release, which is why the refusal used to "
                         "say no prospective result had been measured.")
    ap.add_argument("--landscape", type=Path,
                    default=root / "data/results/pathway_landscape.json",
                    help="pathway_landscape.py --out run whose degree-matched "
                         "null the caveats quote. Absent, the pathway route "
                         "is reported as unmeasured rather than assumed null.")
    ap.add_argument("--virus-label", default="SARS-CoV-2",
                    help="with --no-all-viruses, the single results-file key "
                         "to load a calibration for")
    ap.add_argument("--all-viruses", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="load a calibration for every virus the results file "
                         "holds (default). --no-all-viruses loads only "
                         "--virus-label, which is how this behaved when an "
                         "agent asked about Zika could quote SARS-CoV-2's AUC.")
    return ap


def build_orchestrator(args, note=None) -> Orchestrator:
    """Everything between a parsed namespace and a ready Orchestrator.

    EXTRACTED SO THE SERVER CANNOT LOAD IT DIFFERENTLY. This was the body
    of main(), which meant scripts/agent_server.py had two options: import
    main() and get a CLI, or copy sixty lines of loading and let the two
    drift. The second is how the agent came to quote SARS-CoV-2's AUC for
    a question about Zika -- one caller passed a single calibration where
    another passed a dict -- and how --results was handed to
    measured_reach as a file for weeks. The project has a name for this:
    a contract a caller can forget is not a contract.

    So there is one loader, it takes the namespace `build_parser` produces,
    and the server builds its namespace from that same parser rather than
    from a dict of its own defaults. A default that moves, moves for both.

    `note` receives the degradation messages -- an unreadable results file,
    a refused landscape, an unmeasurable reach. None sends them to stderr,
    which is what the CLI wants; the server collects them and serves them
    as data, because a note printed into a server log is a note the user
    asking the question never sees.
    """
    if note is None:
        def note(msg: str) -> None:
            # "note: " is the prefix the CLI printed before this was
            # extracted, kept so stderr reads the same as it always has.
            print(f"note: {msg}", file=sys.stderr)

    # The caveats quote measured AUCs, so the run that measured them is
    # loaded here and named in the text. A missing or unreadable results file
    # is not fatal: the route is then reported as uncalibrated, which is
    # weaker than a number and true.
    # EVERY VIRUS THE FILE HOLDS, not just one.
    #
    # This loaded a single calibration for --virus-label and handed it to
    # GraphTools, which then used it for whatever virus the question named.
    # Asked about Zika, the agent quoted SARS-CoV-2's AUC. Now each virus gets
    # its own, keyed by the label the results file uses, and a virus the file
    # does not cover resolves to None -- reported as uncalibrated rather than
    # borrowing another virus's numbers.
    cals: dict[str, object] = {}

    # THE PARAGRAPH THAT STOOD HERE DESCRIBED LOADING THIS FILE NO LONGER
    # DOES. It said the prospective calibration was "read from the same
    # document, by the same loader, with section='temporal'" -- true
    # until 2026-11-20, left in place when that loop was replaced by the
    # dedicated file and loader below. A comment that describes the code
    # it sits above, inaccurately, is how a maintainer "re-fixes" a thing
    # back to the defect. The history is in the commit; the current
    # arrangement is documented where it happens.
    temporal: dict[str, object] = {}
    if args.results.exists():
        try:
            labels = (kgav_calibration.available_viruses(args.results)
                      if args.all_viruses else [args.virus_label])
            for label in labels:
                cals[label] = kgav_calibration.load(args.results, label)
        except (OSError, ValueError, KeyError) as exc:
            note(f"could not read {args.results.name} ({exc}); "
                 f"routes will be reported as uncalibrated")

    # THE PROSPECTIVE CALIBRATION, FROM ITS OWN FILE.
    #
    # This used to read --results' `temporal` section with the
    # cross-sectional loader. That section is `{}` on this release --
    # hard_negatives.py skipped the temporal protocol, and records that it
    # did under `temporal_skipped` -- so every calibration it produced was
    # unusable and the refusal took its no-measurement branch: "this
    # release carries no measured prospective result". Meanwhile
    # temporal_structure_scorer.json held ten evaluated cells with six
    # clearing chance.
    #
    # So the agent was understating its own evidence. That is the safe
    # direction and it is still wrong: guards._temporal_sentence has a
    # branch for a scorer that clears prospectively, written precisely so
    # the reader is told before being refused, and nothing could reach it.
    #
    # A SEPARATE LOADER BECAUSE IT IS A SEPARATE SHAPE. The prospective
    # study keys `cells` by taxon CURIE with each scorer beside per-cell
    # counts; `load` expects {virus: {channel: row}}. Pointed at this file
    # it returns an empty calibration -- no error, no note, just silence,
    # which is the failure mode that hid this for eleven days.
    #
    # THE REFUSAL ITSELF DOES NOT MOVE. guards.check returns
    # allowed=False for every PROSPECT query whatever the calibration
    # says; `temporal` shapes the reason and nothing else, and
    # temporal_for feeds that one call site. Checked before wiring it.
    if args.temporal_results.exists():
        try:
            taxa = (kgav_calibration.available_temporal_viruses(
                args.temporal_results)
                if args.all_viruses else [args.virus_label])
            for taxon in taxa:
                # THE RELEASE, so a file measured on another graph is
                # refused rather than quoted. landscape.load and
                # measured_reach are both given it for the same reason.
                temporal[taxon] = kgav_calibration.load_temporal(
                    args.temporal_results, taxon, release=args.release)
            # A FILE THAT YIELDS NOTHING MUST SAY SO. A refusal record
            # carries no `cells`, so available_temporal_viruses returns
            # [] and this loop runs zero times -- silently, which is the
            # shape of the defect this whole change was fixing.
            if not temporal:
                note(f"{args.temporal_results.name} yielded no prospective "
                     f"cell (a refusal record, or a run that evaluated "
                     f"none); the candidate-generation refusal will "
                     f"report no measured prospective result")
            elif not any(getattr(c, "usable", False)
                         for c in temporal.values()):
                reasons = {getattr(c, "refused_reason", "")
                           for c in temporal.values()}
                note(f"every prospective cell in "
                     f"{args.temporal_results.name} is unusable: "
                     + "; ".join(sorted(r for r in reasons if r)))
        except (OSError, ValueError, KeyError) as exc:
            note(f"could not read {args.temporal_results.name} ({exc}); "
                 f"the candidate-generation refusal will report no "
                 f"measured prospective result")
    else:
        note(f"{args.temporal_results.name} is absent; the "
             f"candidate-generation refusal will report no measured "
             f"prospective result. `make temporal-structure-scorer` "
             f"produces it.")
    # THE PATHWAY NULL, loaded the same way and for the same reason. Passing
    # the release name lets landscape.load refuse a run measured on a
    # different one rather than quote it. An absent or refused file is not
    # fatal: the caveat then says the overlap has not been measured on this
    # release, which is weaker than the finding and true.
    land = kgav_landscape.load(args.landscape, release=args.release.name)
    if not land.usable:
        note(f"{land.refused_reason}; the pathway route will be "
             f"reported as unmeasured")
    # MEASURED REACH, the fact that separates `evaluable` from ROUTES
    # REACH NO LABEL. Loaded like the two above and refused the same
    # way: a metapath_comparison.json that is absent, stale or stamped
    # for another release gives None, and the agent then says reach is
    # unmeasured instead of asserting it is zero.
    #
    # .parent BECAUSE --results IS A FILE AND THIS WANTS THE DIRECTORY.
    # hard_negatives.json is one run in data/results; metapath_comparison
    # .json is another, and measured_reach reads the second by name. Passed
    # args.results, it looked for
    # data/results/hard_negatives.json/metapath_comparison.json, which
    # cannot exist -- so the agent reported reach as unmeasured on every
    # question ever asked of it, including the five viruses whose reach is
    # measured and is zero, which is the exact caveat measured_reach was
    # written to stop. Returning None for both "absent" and "you handed me
    # the wrong kind of path" is what hid it; labels.measured_reach now
    # refuses the second rather than reporting it as the first.
    reach = measured_reach(args.results.parent, args.release)
    if reach is None:
        note("metapath_comparison.json is absent, stale or stamped "
             "for another release; whether a route reaches a label will "
             "be reported as unmeasured. `make check-results` says which.")
    orch = Orchestrator(GraphTools(args.release, calibration=cals or None,
                                   temporal=temporal or None,
                                   landscape=land, reach=reach))

    return orch


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    args = build_parser(root).parse_args()
    question = " ".join(args.question)

    if not (args.release / "nodes.jsonl").exists():
        print(f"release not found at {args.release}", file=sys.stderr)
        return 2

    orch = build_orchestrator(args)

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

    # A clarification is a system message, not a synthesised claim, so it is
    # not verified. Keying that on "has no edge ids" was WRONG and dangerous:
    # triage and profile legitimately produce no edge ids, so the model's
    # output went unverified and it invented E335139, E335140 and E335141,
    # which printed. A clarification is identified by the orchestrator having
    # set it, not by the shape of the result.
    verification = None
    if out["allowed"] and out.get("answer") and not out.get("clarification"):
        verification = verify_response(out)

    if args.json:
        print(json.dumps(json_payload(out, verification), indent=2))
        return 0 if (verification is None or verification.passed) else 1

    print(render(out, verification, brief_only=args.no_model))
    return 0 if (verification is None or verification.passed) else 1


if __name__ == "__main__":
    raise SystemExit(main())
