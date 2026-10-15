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


# ------------------------------------------- schema failure had nine copies
def test_no_script_prints_its_own_schema_failure():
    """NINE near-identical copies of the SCHEMA FAIL reporter existed: one
    per ingest script plus assemble_release. Eight printed the head of the
    list; the ninth -- the one not shared -- grouped by violation code
    first and was strictly better.

    kgav.schema.report_violations is now the only one. A new ingest script
    that rolls its own would lose the per-code breakdown silently, so the
    sweep is the thing that keeps it from happening a tenth time.
    """
    # BY AST, NOT BY TEXT. The first version matched the literal string
    # anywhere in the file and fired on dryrun_promote.py, where "SCHEMA
    # FAIL" appears in a comment explaining a past incident. A sweep that
    # cannot tell prose from a print call is a sweep someone will disable.
    offenders = []
    for p in sorted((ROOT / "scripts").glob("*.py")):
        tree = ast.parse(p.read_text())
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == "print"):
                continue
            emitted = "".join(
                v.value for a in node.args for v in ast.walk(a)
                if isinstance(v, ast.Constant) and isinstance(v.value, str))
            if "SCHEMA FAIL" in emitted:
                offenders.append(f"{p.name}:{node.lineno}")
    assert not offenders, (
        f"{offenders} print SCHEMA FAIL directly instead of calling "
        f"kgav.schema.report_violations")


def test_a_rare_violation_code_is_not_buried_under_a_common_one():
    """The eight copies took violations[:8] off the head of one list. The
    folds layer reports 127 violations and showed 8, all [DANGLING]; a
    second code among the other 119 was unanswerable from the output."""
    import io
    from contextlib import redirect_stdout
    from dataclasses import dataclass

    from kgav.schema import report_violations

    @dataclass
    class V:
        code: str
        ref: str = ""
        message: str = "m"

        def __str__(self):
            return f"[{self.code}] {self.ref}: {self.message}"

    vs = [V("DANGLING", f"n{i}") for i in range(125)] + [V("ENUM", "x")]
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = report_violations(vs)
    out = buf.getvalue()
    assert rc == 1
    assert "126 violations, 2 distinct codes" in out
    assert "[ENUM]" in out                       # the rare one is shown
    assert "and 122 more [DANGLING]" in out      # the remainder is counted

    # and a single kind does not read as "1 distinct code(s)"
    buf = io.StringIO()
    with redirect_stdout(buf):
        report_violations([V("DANGLING", "n")])
    assert "1 distinct code\n" in buf.getvalue()


def test_a_clean_batch_prints_pass_and_returns_zero():
    import io
    from contextlib import redirect_stdout

    from kgav.schema import report_violations
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = report_violations([])
    assert rc == 0 and "schema: PASS" in buf.getvalue()


# --------------------------------- a printed list must say what it left out
def test_no_printing_loop_silently_drops_the_rest():
    """FOURTEEN loops printed a constant-capped slice. Exactly ONE said how
    many it had left out -- diagnose_flavivirus_join's "... N more". Nine
    were copies of the schema reporter; the other four printed a head and
    stopped, so ten conflicts and four hundred read identically.

    It is the same failure as the join list that counted fourteen dangling
    chains and named ten, in the report that gates LAYER_PRECEDENCE, and as
    the allowed-transfer list cut at 14 of 36 while a comment beside it
    drew a conclusion the hidden 22 refuted.

    kgav.schema.print_capped is the honest one's shape. The sweep is what
    stops a fifteenth.
    """
    offenders = []
    for p in sorted((ROOT / "scripts").glob("*.py")):
        tree = ast.parse(p.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.For):
                continue
            it = node.iter
            if not (isinstance(it, ast.Subscript)
                    and isinstance(it.slice, ast.Slice)
                    and it.slice.lower is None
                    and isinstance(it.slice.upper, ast.Constant)):
                continue
            if any(isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                   and n.func.id == "print" for n in ast.walk(node)):
                offenders.append(f"{p.name}:{node.lineno}")
    assert not offenders, (
        f"{offenders} print a capped slice directly; use "
        f"kgav.schema.print_capped so the remainder is counted")


def test_print_capped_states_the_remainder_and_stays_quiet_otherwise():
    import io
    from contextlib import redirect_stdout

    from kgav.schema import print_capped

    buf = io.StringIO()
    with redirect_stdout(buf):
        print_capped(list(range(400)), 10, unit="more conflict(s)")
    out = buf.getvalue()
    assert out.count("\n") == 11                      # 10 rows + the tail
    assert "... and 390 more conflict(s)" in out

    buf = io.StringIO()
    with redirect_stdout(buf):
        print_capped([1, 2], 10)
    assert buf.getvalue() == "    1\n    2\n"          # nothing left out


def test_print_capped_counts_with_a_thousands_separator():
    """124 and 1240 must not read alike in a report people skim."""
    import io
    from contextlib import redirect_stdout

    from kgav.schema import print_capped
    buf = io.StringIO()
    with redirect_stdout(buf):
        print_capped(list(range(1250)), 10)
    assert "... and 1,240 more" in buf.getvalue()
