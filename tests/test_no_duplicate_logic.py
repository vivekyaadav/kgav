"""One calculation, one implementation -- enforced instead of rediscovered.

This project has found the same defect ten times: percentile arithmetic in
combine and diagnose_controls; `trustworthy` in baselines and temporal_split;
hits_range; the metapath constraints; the Emit class; `_p`; the "484-odd"
figure; find_chain; the vhppi gate order; and `_vs`. Every one was found by
tripping over a wrong answer, and each time the fix reached one site.

So the sweep that found the last two runs on every commit. It compares
FUNCTION BODIES by AST -- names, docstrings and string contents normalised
away -- because the hazard is textual duplication that drifts, and a copy
renamed or reworded is the same hazard.

ALLOWED lists the pairs that are deliberate. A duplicate belongs there with
a reason, not because the test was noisy.
"""
import ast
import hashlib
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
MIN_STATEMENTS = 3

#: Deliberate duplicates. Key is the sorted tuple of "path:function".
#:
#: EMPTY, and it started with an entry. phase0_verify_folds.find_chain and
#: folds.find_chain were byte-identical until folds.find_chain learned to
#: keep the PRO_ id, so they have already drifted -- and the second test
#: below refused the exemption on exactly that ground. phase0 is a frozen
#: verification script whose independence is the point of one, so the drift
#: is correct; what would not have been correct is an exemption sitting here
#: afterwards, quietly covering those two paths for whatever lands next.
ALLOWED: dict[tuple[str, ...], str] = {}


class _Normalise(ast.NodeTransformer):
    """Collapse string contents so two copies that differ only in wording,
    a message or a docstring still register as the same body."""

    def visit_Constant(self, node):
        if isinstance(node.value, str):
            return ast.copy_location(ast.Constant(value="S"), node)
        return node


def _body_signature(fn: ast.AST) -> tuple[str | None, int]:
    f = ast.parse(ast.unparse(fn)).body[0]
    if (f.body and isinstance(f.body[0], ast.Expr)
            and isinstance(f.body[0].value, ast.Constant)
            and isinstance(f.body[0].value.value, str)):
        f.body = f.body[1:]                      # drop the docstring
    if len(f.body) < MIN_STATEMENTS:
        return None, 0
    f = _Normalise().visit(f)
    ast.fix_missing_locations(f)
    return hashlib.sha1(ast.unparse(f).encode()).hexdigest()[:12], len(f.body)


def _collect() -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for path in sorted(ROOT.glob("src/**/*.py")) + sorted(ROOT.glob("scripts/*.py")):
        rel = path.relative_to(ROOT).as_posix()
        try:
            tree = ast.parse(path.read_text())
        except SyntaxError:                      # pragma: no cover
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                sig, _n = _body_signature(node)
                if sig:
                    out.setdefault(sig, []).append(f"{rel}:{node.name}")
    return out


def test_no_function_body_is_implemented_twice():
    offenders = []
    for sig, sites in _collect().items():
        files = {s.rsplit(":", 1)[0] for s in sites}
        if len(files) < 2:
            continue                              # overloads in one file are fine
        key = tuple(sorted(sites))
        if key in ALLOWED:
            continue
        offenders.append(f"{sig}: " + ", ".join(sorted(sites)))
    assert not offenders, (
        "the same function body exists in more than one file:\n  "
        + "\n  ".join(offenders)
        + "\n\nMove it to one module and import it, or add the pair to "
          "ALLOWED with the reason it is deliberate. This project has found "
          "this defect ten times; every one was found by hitting a wrong "
          "answer first.")


def test_the_allowed_list_does_not_outlive_its_entries():
    """A stale exemption is how a real duplicate gets waved through later."""
    present = {tuple(sorted(sites)) for sites in _collect().values()
               if len({s.rsplit(":", 1)[0] for s in sites}) > 1}
    stale = [k for k in ALLOWED if k not in present]
    assert not stale, (
        f"ALLOWED exempts pairs that are no longer duplicated: {stale}. "
        "Remove them, or the next real duplicate at those paths passes.")
