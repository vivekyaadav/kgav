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
