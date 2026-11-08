"""Gate for the CLI layer: a failing answer is withheld, not annotated."""
import pytest

from kgav.agent.llm import DEFAULT_MODEL, ModelUnavailable, OllamaModel


def test_adapter_is_deterministic_by_default():
    """A researcher who runs the same query twice and gets two answers cannot
    cite the tool, and everything below this point is already deterministic."""
    m = OllamaModel()
    assert m.temperature == 0.0
    assert m.model == DEFAULT_MODEL


def test_unreachable_host_raises_a_distinct_error():
    """Distinct from a bad answer: the caller should fall back to printing the
    brief, which is still useful, not report an empty answer."""
    m = OllamaModel(host="http://127.0.0.1:9", timeout=2)
    with pytest.raises(ModelUnavailable):
        m("hello")


def test_availability_check_does_not_raise():
    ok, note = OllamaModel(host="http://127.0.0.1:9").available()
    assert ok is False and isinstance(note, str)


# --------------------------------------------------------------- rendering
def _out(answer, facts, ids, warnings, allowed=True):
    return {"allowed": allowed, "answer": answer, "facts": facts,
            "citable_edge_ids": ids, "warnings": warnings, "intent": "explain",
            "trace": []}


def test_failed_verification_withholds_the_answer():
    """A tool that shows a hallucinated answer with a warning underneath is a
    tool that shows hallucinated answers."""
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from kgav_ask import render

    from kgav.agent.verify import verify_response

    out = _out("Compound X inhibits nsp5 [E999].",
               ["real fact [E1]"], ["E1"], [])
    v = verify_response(out)
    text = render(out, v)
    assert "ANSWER WITHHELD" in text
    assert "Compound X inhibits nsp5" not in text


def test_passing_answer_is_shown():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from kgav_ask import render

    from kgav.agent.verify import verify_response

    out = _out("Nirmatrelvir inhibits nsp5 [E1].",
               ["nirmatrelvir inhibits nsp5 [chembl, 2021; E1]"], ["E1"], [])
    v = verify_response(out)
    assert v.passed, v.report()
    assert "ANSWER" in render(out, v)


def test_refusal_is_rendered_without_evidence():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from kgav_ask import render

    out = _out("This graph cannot support that.", [], [], [], allowed=False)
    text = render(out, None)
    assert "REFUSED" in text and "EVIDENCE" not in text


# ------------------------------------- A6: the two surfaces must agree
# `render` withheld a failing answer; the --json branch printed `out`
# whole, answer included. The exit code was 1 either way, which protects a
# shell `&&` and nothing that reads the JSON: `jq -r .answer` returned
# model text that had failed verification, with the failure in a sibling
# key. The file docstring above says "a failing answer is withheld, not
# annotated" and had only ever been true of one of the two surfaces.


def _surfaces():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from kgav_ask import json_payload, render
    return render, json_payload


def test_the_json_surface_withholds_what_the_text_surface_withholds():
    """Asked of BOTH surfaces for one case, so neither can be fixed alone.

    The hallucinated citation is the same fixture
    test_failed_verification_withholds_the_answer uses, which is the point:
    that test passed throughout while this one would have failed.
    """
    from kgav.agent.verify import verify_response
    render, json_payload = _surfaces()

    out = _out("Compound X inhibits nsp5 [E999].",
               ["real fact [E1]"], ["E1"], [])
    v = verify_response(out)
    assert not v.passed, "the fixture no longer fails verification"

    text = render(out, v)
    payload = json_payload(out, v)
    assert "Compound X inhibits nsp5" not in text
    assert payload["answer"] is None, (
        "the JSON surface still carries an answer that failed verification")
    assert "Compound X inhibits nsp5" not in json_dumps(payload), (
        "the unverified text survives somewhere else in the payload")
    assert payload["answer_withheld"]
    # The graph's own output is NOT withheld: the facts are what it
    # contains and the failure is the model's.
    assert payload["facts"] == ["real fact [E1]"]
    assert payload["verification"]["passed"] is False
    assert payload["verification"]["findings"], "the diagnosis is gone too"


def test_the_json_surface_shows_an_answer_that_passes():
    """The other half: a check that withholds everything is a switch that
    is always off."""
    from kgav.agent.verify import verify_response
    _render, json_payload = _surfaces()

    out = _out("Real fact [E1].", ["real fact [E1]"], ["E1"], [])
    v = verify_response(out)
    assert v.passed, v.findings
    payload = json_payload(out, v)
    assert payload["answer"] == "Real fact [E1]."
    assert "answer_withheld" not in payload


def test_the_payload_carries_no_prompt():
    """It never did, and it is the one key that would leak the model's
    instructions into a machine-readable surface."""
    _render, json_payload = _surfaces()
    out = _out("x", [], [], [])
    out["prompt"] = "SYSTEM PROMPT"
    assert "prompt" not in json_payload(out, None)


def json_dumps(obj) -> str:
    import json
    return json.dumps(obj, default=str)


# ------------------------------------- A8: every failure is ModelUnavailable
# Marked *(suspected)* in the audit; measured, with only urllib.error.URLError
# caught, against four stub servers:
#
#   non-JSON body            -> JSONDecodeError ESCAPED
#   read timeout             -> TimeoutError ESCAPED
#   HTTP 500                 -> handled, but told the user to check whether a
#                               service that had just answered was running
#   JSON with no "response"  -> returned "" as if it were an answer
#
# The timeout is the sharpest: the class SETS timeout=180 and catches
# URLError, but a read timeout raises socket.timeout -- TimeoutError since
# 3.10 -- straight out of urlopen unwrapped. The one failure it configures
# for was the one it did not handle.

import http.server
import threading
import time


def _stub(handler_cls):
    srv = http.server.HTTPServer(("127.0.0.1", 0), handler_cls)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, f"http://127.0.0.1:{srv.server_address[1]}"


class _Quiet(http.server.BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass


class _NotJson(_Quiet):
    def do_POST(self):
        self.send_response(200)
        self.send_header("Content-Type", "text/html")
        self.end_headers()
        self.wfile.write(b"<html>a proxy error page</html>")


class _ServerError(_Quiet):
    def do_POST(self):
        self.send_response(500)
        self.end_headers()
        self.wfile.write(b"boom")


class _NeverAnswers(_Quiet):
    def do_POST(self):
        # Three seconds, against a one-second client timeout. Thirty made
        # the suite 30s slower, because shutdown() waits for the handler
        # and the point is only to outlast the client.
        time.sleep(3)


class _NoResponseField(_Quiet):
    def do_POST(self):
        import json as _json
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        self.wfile.write(_json.dumps({"error": "model not found"}).encode())


class _Good(_Quiet):
    def do_POST(self):
        import json as _json
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        self.wfile.write(_json.dumps({"response": "  the answer  "}).encode())


@pytest.mark.parametrize("handler,timeout,expect", [
    (_NotJson, 5, "not"),
    (_ServerError, 5, "HTTP 500"),
    (_NeverAnswers, 1, "did not reply within 1s"),
    (_NoResponseField, 5, "model not found"),
])
def test_every_adapter_failure_is_model_unavailable(handler, timeout, expect):
    """ModelUnavailable is the right class for all of them: its own
    docstring says the caller should fall back to printing the brief
    unsynthesised, which is still useful, and that is as true of a garbled
    reply as of an unreachable socket. A crash loses the brief too."""
    srv, host = _stub(handler)
    try:
        with pytest.raises(ModelUnavailable) as got:
            OllamaModel(host=host, timeout=timeout)("hello")
        assert expect in str(got.value), str(got.value)
    finally:
        srv.shutdown()


def test_an_http_error_does_not_blame_the_service_for_being_down():
    """It answered. "Is the service running? systemctl status ollama" sends
    the reader after the wrong thing."""
    srv, host = _stub(_ServerError)
    try:
        with pytest.raises(ModelUnavailable) as got:
            OllamaModel(host=host, timeout=5)("hello")
        assert "Is the service running" not in str(got.value)
    finally:
        srv.shutdown()


def test_a_good_reply_is_returned_stripped():
    """The other half: a check that rejects everything is an adapter that
    never answers."""
    srv, host = _stub(_Good)
    try:
        assert OllamaModel(host=host, timeout=5)("hello") == "the answer"
    finally:
        srv.shutdown()
