"""What a caller of the agent's machine-readable surface can receive.

Two properties, both of which have been broken in this repository before.

ONE: AN ANSWER THAT FAILED VERIFICATION IS NOT IN THE PAYLOAD. Not
flagged in it, not beside a failure key -- absent. `render` has withheld
it since fix-list A6; `--json` printed `out` whole until the same fix, so
`jq -r .answer` returned model text the terminal refused to show. The
test for that lives with the two surfaces, and this file adds the case
that surface now has a second consumer: scripts/agent_server.py serves
this payload to a browser, where an answer present in JSON is an answer
present in the DOM.

TWO: NO PROMPT REACHES A CALLER. `json_payload` filtered `k != "prompt"`
-- one key, named by hand -- and the orchestrator sets two, the second
being `retry_prompt` whenever the first answer dropped a required
warning. Every --json reply to such a question carried the system rules,
the whole brief and the rejection preamble, in the surface whose author
had deliberately stripped the other prompt. Measured on
"why might nirmatrelvir work against SARS-CoV-2?": 2,881 characters.

The filter is now by shape, so this file asserts the property rather
than the spelling: whatever prompt-shaped keys the orchestrator grows,
none of them may arrive.
"""
from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _kgav_ask():
    """The script as a module; it is not importable as a package member."""
    spec = importlib.util.spec_from_file_location(
        "kgav_ask_json_surface", ROOT / "scripts" / "kgav_ask.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@dataclass
class FakeFinding:
    severity: str = "fail"
    code: str = "uncited_claim"
    detail: str = "a sentence with no edge id behind it"


@dataclass
class FakeVerification:
    passed: bool
    findings: list = field(default_factory=list)
    cited: list = field(default_factory=list)


#: An orchestrator result shaped like the real one, including BOTH prompt
#: keys. Spelled here rather than produced by running the agent so the
#: test needs no release and no model -- the property is about the
#: filter, not about any particular question.
def _out(answer: str = "NIRMATRELVIR inhibits nsp5 [E83e08310cff11af7]."):
    return {
        "question": "why might nirmatrelvir work against SARS-CoV-2?",
        "intent": "explain",
        "allowed": True,
        "trace": ["intent=explain"],
        "facts": ["Route 1 (M1 direct-acting):"],
        "warnings": ["VIRUS NOT VALIDATED: ..."],
        "citable_edge_ids": ["E83e08310cff11af7"],
        "answer": answer,
        "attempts": 2,
        "dropped_on_first_attempt": ["VIRUS NOT VALIDATED: ..."],
        "prompt": "You are reporting findings from a curated knowledge graph.",
        "retry_prompt": "YOUR PREVIOUS ANSWER WAS REJECTED. It omitted ...",
    }


# ------------------------------------------------------------ property one
def test_a_failed_verification_leaves_no_answer_in_the_payload():
    payload = _kgav_ask().json_payload(
        _out(), FakeVerification(passed=False, findings=[FakeFinding()]))
    assert payload["answer"] is None
    assert "answer_withheld" in payload
    # The findings survive -- they are the diagnostic value of the failure.
    assert payload["verification"]["findings"][0]["code"] == "uncited_claim"


def test_the_withheld_answer_text_is_nowhere_in_the_payload():
    """Absent, not merely nulled in one key.

    The model's sentence must not survive in a sibling -- which is how
    the original defect read to anyone checking only `.answer`.
    """
    # Named `canary` on purpose. Every diff in this repository is swept
    # for credential-shaped words before it is pushed, and a fixture that
    # trips that sweep on every future commit trains the reader to skim
    # it. The variable holds a sentence no graph supports; that is its
    # whole job.
    canary = "nirmatrelvir cures covid in every patient"
    payload = _kgav_ask().json_payload(
        _out(answer=canary),
        FakeVerification(passed=False, findings=[FakeFinding()]))

    def strings(obj):
        if isinstance(obj, str):
            yield obj
        elif isinstance(obj, dict):
            for v in obj.values():
                yield from strings(v)
        elif isinstance(obj, (list, tuple)):
            for v in obj:
                yield from strings(v)

    assert not [s for s in strings(payload) if canary in s]


def test_a_passing_verification_keeps_the_answer():
    """The guard above must not be satisfied by withholding everything."""
    payload = _kgav_ask().json_payload(
        _out(), FakeVerification(passed=True, cited=["E83e08310cff11af7"]))
    assert payload["answer"]
    assert "answer_withheld" not in payload
    assert payload["verification"]["passed"] is True


# ------------------------------------------------------------ property two
def test_no_prompt_shaped_key_reaches_the_caller():
    """Both prompt keys, and any future sibling, by shape."""
    for verification in (None,
                         FakeVerification(passed=True),
                         FakeVerification(passed=False,
                                          findings=[FakeFinding()])):
        payload = _kgav_ask().json_payload(_out(), verification)
        leaked = sorted(k for k in payload if "prompt" in k)
        assert not leaked, (
            f"{leaked} reached the caller. The orchestrator's prompts carry "
            f"the system rules and the whole brief, and --show-prompt is the "
            f"opt-in for reading one.")


def test_the_fixture_really_carries_two_prompt_keys():
    """THE GUARD ON THE TEST ABOVE.

    A payload filter cannot be shown to strip prompts by feeding it a
    result that has none -- that test passes on a filter that does
    nothing at all, which is this project's "a check that cannot fail is
    worse than no check". The orchestrator sets `prompt` at two call
    sites and `retry_prompt` at one; both spellings are in the fixture,
    and this asserts they are still there.
    """
    assert sorted(k for k in _out() if "prompt" in k) == [
        "prompt", "retry_prompt"]


def test_the_orchestrator_sets_no_prompt_key_this_file_does_not_know_about():
    """So a third prompt key is a failure here, not a silent leak.

    The filter is by shape and would catch it anyway; this names the
    source of truth so the fixture above cannot quietly fall behind the
    code it stands in for.
    """
    import ast

    src = (ROOT / "src" / "kgav" / "agent" / "orchestrator.py").read_text()
    found = set()
    for node in ast.walk(ast.parse(src)):
        # out["<something>prompt<something>"] = ...
        if (isinstance(node, ast.Subscript)
                and isinstance(node.value, ast.Name)
                and node.value.id == "out"
                and isinstance(node.slice, ast.Constant)
                and isinstance(node.slice.value, str)
                and "prompt" in node.slice.value):
            found.add(node.slice.value)
    assert found == {"prompt", "retry_prompt"}, (
        f"orchestrator.py sets prompt-shaped keys {sorted(found)}; the "
        f"fixture in this file covers {['prompt', 'retry_prompt']}. Add the "
        f"new one to _out() so the filter is tested against it.")
