"""The HTTP surface, and the one claim the frontend is allowed to make.

WHY THIS FILE EXISTS. scripts/agent_server.py is a second consumer of the
agent, and the first one that puts its output in front of someone who is
not reading a terminal. Every property the CLI protects has to survive
the move: an answer that failed verification must not arrive, no prompt
must arrive, and the figures on the page must be the release's own.

THE PREVIOUS FRONTEND IS WHY THE LAST ONE IS CHECKED. It typed five
release counts into web/client/src/pages/Home.tsx under the sentence "the
layer counts below are real figures from release v0.1". Four were wrong,
by up to elevenfold. The test that caught it had to find the numbers by
running a regex over TSX, and it could only ever check the five that
happened to be there -- a sixth class could go unshown with every one of
them correct, which it did.

web/chat.html is checked a different way, because it is built a different
way: it hand-types NO figure, so the test is that there is nothing in it
to check. A page that fetches /meta cannot disagree with the release.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
UI = ROOT / "web" / "chat.html"
RELEASE = ROOT / "data" / "releases" / "v0.1"

needs_release = pytest.mark.skipif(
    not (RELEASE / "nodes.jsonl").exists(),
    reason="no assembled release on disk")


def _server_module():
    spec = importlib.util.spec_from_file_location(
        "agent_server_under_test", ROOT / "scripts" / "agent_server.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ===================================================== the frontend file
def test_the_ui_exists_where_the_server_looks_for_it():
    """The server's UI constant and the file on disk, not two opinions."""
    assert UI.exists(), f"{UI} is missing"
    assert _server_module().UI == UI


@needs_release
def test_the_ui_hand_types_no_figure_from_the_release():
    """THE WHOLE REASON THE OLD UI's COUNT TEST COULD BE DELETED.

    Asked against the release rather than against a digit pattern. A
    "four or more digits" rule has to be taught about hex colours and
    millisecond divisors, and a rule with exemptions drifts into
    exempting the thing it was written for. This asks the question
    directly: does any number written into the page equal a figure the
    release reports?

    It is the exact failure the old frontend had -- five counts typed
    under "real figures from release v0.1", four of them wrong -- and it
    cannot be passed by a page that happens to have no four-digit
    literals today.
    """
    manifest = json.loads((RELEASE / "MANIFEST.json").read_text())
    census = dict(manifest["nodes_by_class"])
    census["__total__"] = sum(manifest["nodes_by_class"].values())
    for key in ("edges", "nodes", "edge_count", "node_count"):
        if isinstance(manifest.get(key), int):
            census[key] = manifest[key]

    # FOUR FIGURES OR MORE, and the small ones are named rather than
    # silently skipped. A stylesheet contains 10px, 16px and 14px; this
    # release contains Disease 10 and OrganismTaxon 16. Those collide by
    # coincidence, and a test that reports a coincidence as a defect gets
    # relaxed until it reports nothing. So the guarantee is bounded and
    # stated: a census figure of 1,000 or more cannot be written into this
    # page. The classes below that floor are listed in the failure text so
    # the next reader knows what is NOT covered here.
    floor = 1000
    small = {k: v for k, v in census.items() if v < floor}
    figures = {str(v): k for k, v in census.items() if v >= floor}

    literals = set(re.findall(r"\b\d+\b", UI.read_text()))
    typed = sorted(figures[n] + "=" + n for n in figures.keys() & literals)
    assert not typed, (
        f"web/chat.html contains {typed}, which are figures this release "
        f"reports. Every number the page shows must arrive from /meta or "
        f"from an /ask payload.")
    assert small, (
        "no census figure is under the floor, so the exemption above is "
        "dead and the floor can be removed")


def test_the_ui_names_no_sha():
    """A 16-hex-digit literal is a release fingerprint, wherever it sits."""
    found = re.findall(r"\b[0-9a-f]{16}\b", UI.read_text())
    assert not found, (
        f"web/chat.html names {found}; the fingerprint belongs to the "
        f"release and reaches the page through /meta.")


def test_the_ui_asks_the_server_for_its_figures():
    """The other half: a page that types nothing AND fetches nothing shows
    nothing, and would pass the two tests above."""
    src = UI.read_text()
    assert '"/meta"' in src, "the page never fetches /meta"
    assert '"/ask"' in src, "the page never posts to /ask"


def test_the_ui_renders_the_withheld_case():
    """A page that only renders `answer` would show an empty reply where
    the CLI shows why the answer is absent. The withheld payload is the
    one case the server can return that has no answer AND is not an
    error, so the page must have a branch for it."""
    src = UI.read_text()
    assert "answer_withheld" in src
    assert "failed verification" in src.lower()


def test_the_ui_loads_nothing_over_the_network_but_this_server():
    """No CDN, no Google Fonts, no analytics. This runs on a machine with
    unpublished results on it and must work with the network off."""
    src = UI.read_text()
    hosts = re.findall(r"https?://[^\s\"'<>)]+", src)
    outside = [h for h in hosts
               if not re.match(r"https?://(127\.0\.0\.1|localhost)", h)]
    assert not outside, f"web/chat.html reaches out to {outside}"


# ======================================================== the server API
def test_the_server_binds_loopback_by_default():
    """The payloads name compounds, targets and unpublished AUCs."""
    args = _server_module().build_parser().parse_args([])
    assert args.bind == "127.0.0.1"


def test_the_server_takes_the_clis_defaults_not_a_copy_of_them():
    """Agent builds its namespace from kgav_ask.build_parser.

    Spelled as an assertion about the parser rather than about the
    values, because the point is the shared source: a default that moves
    must move for both surfaces. The CLI's own wiring test
    (tests/test_agent_cli_wiring.py) is what pins the values.
    """
    src = (ROOT / "scripts" / "agent_server.py").read_text()
    assert "ask.build_parser(ROOT).parse_args" in src, (
        "the server stopped building its arguments with the CLI's parser; "
        "it now has defaults of its own that nothing keeps in step")
    assert "build_orchestrator" in src, (
        "the server stopped using the CLI's loader")


def test_the_server_does_not_reimplement_the_payload():
    """It must call json_payload, not assemble a reply of its own.

    CODE ONLY, NOT PROSE. The first version of this test searched the
    file as text and failed on the module docstring, which describes the
    withholding in order to explain why the server does not implement
    it. A test that cannot tell an explanation from an implementation
    would be fixed by deleting the explanation, which is the wrong
    direction.
    """
    import ast

    src = (ROOT / "scripts" / "agent_server.py").read_text()
    assert "json_payload" in src

    tree = ast.parse(src)
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef,
                             ast.AsyncFunctionDef, ast.ClassDef)):
            body = getattr(node, "body", [])
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                docstrings.add(id(body[0].value))

    spelled = [n.value for n in ast.walk(tree)
               if isinstance(n, ast.Constant) and isinstance(n.value, str)
               and id(n) not in docstrings and "answer_withheld" in n.value]
    assert not spelled, (
        "the server spells the withholding itself; that decision belongs "
        "to json_payload, which both surfaces share")


# ============================================ against the real release
@pytest.mark.slow
@needs_release
def test_the_served_payload_is_the_clis_payload():
    """End to end, on the release, with no model involved.

    Evidence-only on both sides: the model is the one part of this that
    is neither deterministic nor local to the repository, and the
    property under test is that the two surfaces agree about the graph.
    """
    import subprocess

    question = "why might nirmatrelvir work against SARS-CoV-2?"

    srv = _server_module()
    agent = srv.Agent([])
    served = agent.answer(question, use_model=False)

    cli = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "kgav_ask.py"),
         "--no-model", "--json", question],
        capture_output=True, text=True, cwd=ROOT, timeout=300)
    assert cli.returncode == 0, cli.stderr
    printed = json.loads(cli.stdout)

    # The server adds three fields about the request itself and shares
    # everything else. Anything the CLI has that the server lacks, or any
    # shared key whose value differs, is a divergence.
    extra = {"elapsed_seconds", "model_used", "model_note",
             "synthesis_requested"}
    assert set(served) - extra == set(printed), (
        f"server-only keys {sorted(set(served) - extra - set(printed))}, "
        f"CLI-only keys {sorted(set(printed) - set(served))}")
    for key in printed:
        assert served[key] == printed[key], f"{key} differs"


@pytest.mark.slow
@needs_release
def test_meta_reports_the_release_on_disk():
    """Not a remembered one. Every figure the page shows passes through
    here, so this is where a stale number would enter."""
    from kgav.provenance import release_fingerprint

    meta = _server_module().Agent([]).meta()
    fp = release_fingerprint(RELEASE)
    assert meta["nodes_sha"] == fp["nodes_sha"]
    assert meta["edges_sha"] == fp["edges_sha"]

    manifest = json.loads((RELEASE / "MANIFEST.json").read_text())
    assert meta["nodes_by_class"] == manifest["nodes_by_class"]
    assert meta["node_total"] == sum(manifest["nodes_by_class"].values())


@pytest.mark.slow
@needs_release
def test_a_prospective_question_is_still_refused_over_http():
    """The refusal is the safety property most worth re-checking on a new
    surface: it is the one a friendlier interface invites someone to
    route around."""
    served = _server_module().Agent([]).answer(
        "which drugs should we test against Zika?", use_model=False)
    assert served["allowed"] is False
    assert served["answer"], "a refusal with no reason is not a refusal"
