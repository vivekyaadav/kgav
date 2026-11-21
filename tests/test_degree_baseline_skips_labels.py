"""The null hypothesis must not be able to read the answer sheet.

`baselines.degree_ranking` has said since fix-list L7 that it ranks over
the TRAVERSABLE graph, not the whole one. It delivered that only when its
caller passed `skip_predicates` to `Graph.load`, which defaults to the
empty set. `scripts/structure_scorer.py` never passed it, so the three
registered scorer studies compared their structure statistic against a
degree baseline built from the label edges those studies are scored
against -- and nothing failed, in twelve hundred tests, because the
docstring described behaviour the default did not provide.

That was found on 2026-11-17, and only because 577 clinical-label edges
landed and moved a number documented as unable to move. A defect found by
coincidence is a defect with no test, so these are the tests.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from collections import Counter  # noqa: E402
from kgav.baselines import Graph, degree_ranking  # noqa: E402
from kgav.schema import load_schema  # noqa: E402


def test_degree_ranking_refuses_a_graph_that_skipped_nothing():
    """The contract enforced where it matters, not where it is written."""
    g = Graph()
    g.degree["INCHIKEY:A"] = 3
    # the PRIVATE counter: the public one is a property that refuses on
    # exactly this graph, which is the guard added beneath this test.
    g._degree_unlabelled["INCHIKEY:A"] = 3
    with pytest.raises(ValueError, match="answer sheet"):
        degree_ranking(g)


def test_degree_ranking_accepts_a_graph_that_skipped_something():
    """So the refusal above is about the missing argument and not about
    an empty graph."""
    g = Graph()
    g.skipped_predicates = frozenset({"HAS_ANTIVIRAL_ACTIVITY_AGAINST"})
    assert degree_ranking(g) == {}


def test_load_records_what_it_was_asked_to_skip(tmp_path):
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text("")
    (rel / "edges.jsonl").write_text("")
    g = Graph.load(rel, skip_predicates={"TREATS"}, symmetry={})
    assert g.skipped_predicates == frozenset({"TREATS"})
    g2 = Graph.load(rel, skip_predicates=None, symmetry={})
    assert g2.skipped_predicates == frozenset()


def test_every_script_that_ranks_by_degree_loads_a_skipping_graph():
    """THE GENERAL FORM, because the specific one was a single missing
    keyword in a single line.

    Any script calling `degree_ranking` must also name `skip_predicates`
    on its `Graph.load`. The refusal above catches it at runtime; this
    catches it at commit time and names the file, which a stack trace
    from a thirty-minute protocol does not do cheaply.
    """
    offenders = []
    for path in sorted((ROOT / "scripts").glob("*.py")):
        text = path.read_text()
        if "degree_ranking" not in text:
            continue
        loads = re.findall(r"Graph\.load\((.*?)\)", text, re.S)
        if not loads:
            continue
        if not all("skip_predicates" in call for call in loads):
            offenders.append(path.name)
    assert not offenders, (
        f"{offenders} rank by degree over a graph that was not told to skip "
        "the held-out predicates, so the null hypothesis sees the labels")


def test_the_held_out_set_still_contains_both_label_families():
    """The skip is only as good as the declaration behind it. If a label
    predicate stopped being held_out, every caller above would pass the
    argument and skip nothing that matters."""
    held = load_schema().held_out_predicates()
    for pred in ("HAS_ANTIVIRAL_ACTIVITY_AGAINST", "TREATS", "IN_TRIAL_FOR"):
        assert pred in held, f"{pred} is no longer held out of traversal"


def test_the_amendment_states_the_invariant_and_the_invariant_holds():
    """A figure in a registration, derived rather than typed.

    The amendment first said "Twelve scripts in this repository pass
    `schema.held_out_predicates()`", and this test derived that count and
    fired the same day, because `scripts/clinical_labels.py` was written
    and made it thirteen. A TALLY OF CALLERS IS NOT THE CLAIM: what makes
    the omission an omission rather than a convention is that every other
    caller passes it. So the amendment states the invariant and this
    checks the invariant.
    """
    import re as _re

    passing = sorted(
        p.name for p in (ROOT / "scripts").glob("*.py")
        if "held_out_predicates()" in p.read_text()
        and _re.search(r"Graph\.load\(", p.read_text(), _re.S))
    prereg = (ROOT / "doc2" / "registrations" / "prereg_structure_scorer.md").read_text()
    marker = "# Amendment, 2026-11-17: the comparator was reading"
    assert marker in prereg, "the comparator amendment is gone"
    amendment = " ".join(prereg.split(marker)[1].split())
    assert ("Every other `Graph.load` caller in `scripts/` that touches the "
            "held-out set passes `schema.held_out_predicates()`") in amendment
    assert "structure_scorer.py" in passing, (
        "the script the amendment is about no longer passes the argument")
    # And the invariant itself: nothing in scripts/ reads degree_ranking
    # over a graph it did not ask to skip. That is the test above; this
    # asserts the amendment is not describing an empty set.
    assert len(passing) > 1, (
        "only one script passes held_out_predicates(), so 'every other "
        "caller' describes nothing and the amendment's argument is vacuous")


# ---------------------------------------------- the 2026-11-19 class guards
def test_degree_unlabelled_refuses_when_nothing_was_skipped(tmp_path):
    """THE LIE, GUARDED AT THE POINT IT IS TOLD.

    On a graph loaded with no skip, this counter IS `degree`, so the name
    is false. The refusal added on 2026-11-17 was put on `degree_ranking`
    alone, and two days later `repurposing_candidates.py` was found
    reading the counter one level down and publishing a shortlist whose
    graph-evidence column counted the held-out labels.
    """
    rel = tmp_path / "r"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text("")
    (rel / "edges.jsonl").write_text("")
    g = Graph.load(rel, skip_predicates=None, symmetry={})
    with pytest.raises(ValueError, match="identical to degree"):
        _ = g.degree_unlabelled
    # and it is readable the moment the caller has actually skipped
    ok = Graph.load(rel, skip_predicates={"TREATS"}, symmetry={})
    assert ok.degree_unlabelled == Counter()


def test_graph_load_refuses_to_guess_skip_predicates(tmp_path):
    """A CONTRACT A CALLER CAN FORGET IS NOT A CONTRACT.

    The old signature defaulted to None, so a caller who never thought
    about held-out edges got a graph where they are traversable. Omitting
    it is now an error; passing None is a decision and is allowed.
    """
    rel = tmp_path / "r"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text("")
    (rel / "edges.jsonl").write_text("")
    with pytest.raises(TypeError, match="requires skip_predicates"):
        Graph.load(rel, symmetry={})
    assert Graph.load(rel, skip_predicates=None,
                      symmetry={}).skipped_predicates == frozenset()


def test_every_graph_load_in_the_repository_states_its_choice():
    """THE GENERAL FORM. The guards above fail a forgotten argument at
    runtime; this fails one at review time, including in a script no
    test happens to execute.

    AST, NOT A REGEX. The first draft scanned the source text and
    reported six hits, every one of them the string "Graph.load(" inside
    a comment or a docstring -- including the docstrings of the two
    guards themselves. A checker that cannot tell code from prose about
    code fails on the documentation that explains it.
    """
    import ast
    bad = []
    for d in ("scripts", "src", "tests"):
        for p in sorted((ROOT / d).rglob("*.py")):
            tree = ast.parse(p.read_text(), filename=str(p))
            # ALIASES COUNT. tests/test_baselines.py imports
            # `Graph as G` and calls `G.load(...)`; the first version of
            # this check looked for the literal name `Graph` and walked
            # straight past it, so nine tests failed at runtime on a
            # guard this was supposed to have caught at review time.
            aliases = {"Graph"}
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom):
                    aliases |= {a.asname for a in node.names
                                if a.name == "Graph" and a.asname}
            # A CALL WRITTEN TO PROVE THE REFUSAL IS NOT A FORGOTTEN
            # ARGUMENT. The test directly above omits skip_predicates on
            # purpose, inside pytest.raises, so those line ranges are
            # exempt -- and only those. Exempting the file instead would
            # have switched the check off for every other call in it.
            proving = set()
            for node in ast.walk(tree):
                if not isinstance(node, ast.With):
                    continue
                if not any("raises" in ast.dump(i.context_expr)
                           for i in node.items):
                    continue
                proving.update(range(node.lineno, (node.end_lineno or
                                                   node.lineno) + 1))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                f = node.func
                if not (isinstance(f, ast.Attribute) and f.attr == "load"
                        and isinstance(f.value, ast.Name)
                        and f.value.id in aliases):
                    continue
                if any(k.arg == "skip_predicates" for k in node.keywords):
                    continue
                if node.lineno in proving:
                    continue
                bad.append(f"{p.relative_to(ROOT)}:{node.lineno}")
    assert not bad, (
        "these Graph.load calls do not say what to do with the held-out "
        "predicates:\n  " + "\n  ".join(bad))
