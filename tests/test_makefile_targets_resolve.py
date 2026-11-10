"""Every Makefile target must be runnable with no arguments, or say so.

WHY THIS FILE EXISTS. `release-validate` expanded to
`python scripts/validate_release.py data/releases/$(REL)` with no default
for REL anywhere in the Makefile, so it ran `data/releases/` and exited 2
on every invocation it ever had. Nothing was unvalidated -- the same check
runs inside `make rebuild` for every layer -- but the one command a reader
is told to use to check a release by hand had never worked.

tests/test_kg_audit_report.py already asserts that every `make` target the
audit report CITES exists. Existing and running are different claims, and
only the first was tested: a target is a line in a file, and a target that
cannot run is a line in a file. This closes the other half.

THE RULE. A recipe may expand an uppercase variable only if one of these
holds, because each one means a bare `make <target>` still does something
defined:

  1. the Makefile gives it a default (`VAR ?= x` or `VAR = x`)
  2. the expansion supplies one inline (`$(or $(VAR),default)`, `$(if ...)`)
  3. the recipe refuses without it (`test -n "$(VAR)" || { ...; exit; }`)

Case 3 is not a loophole -- `diff-releases` and `compare-results` use it,
and a target that names what it needs and exits 2 is honest. What is
forbidden is the fourth case: silent expansion to the empty string.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAKEFILE = ROOT / "Makefile"

#: `$(NAME)` with an all-caps name. Lowercase `$(shell ...)`, `$(or ...)`
#: and friends are make FUNCTIONS, not variables, and are skipped by the
#: case of the first letter rather than by a list that would go stale.
VAR = re.compile(r"\$\((?P<name>[A-Z][A-Z0-9_]*)\)")

#: Make's own automatic variables, which are always defined.
AUTOMATIC = {"MAKE", "MAKEFLAGS", "CURDIR", "SHELL", "MAKECMDGOALS", "PWD"}


def _recipes() -> dict[str, list[str]]:
    """target name -> its recipe lines, from the Makefile as written."""
    out: dict[str, list[str]] = {}
    current: str | None = None
    for line in MAKEFILE.read_text().splitlines():
        if line.startswith("\t"):
            if current is not None:
                out[current].append(line)
            continue
        m = re.match(r"^([A-Za-z0-9_.-]+):(?!=)", line)
        current = m.group(1) if m else None
        if current is not None:
            out.setdefault(current, [])
    return out


def _defaulted() -> set[str]:
    """Variables the Makefile itself assigns, so a bare `make` has a value."""
    return set(re.findall(r"^([A-Z][A-Z0-9_]*)\s*[?:+]?=",
                          MAKEFILE.read_text(), re.M))


def test_no_recipe_expands_an_undefined_variable():
    """THE TEST THAT WOULD HAVE CAUGHT `REL`.

    Proven to fail against the Makefile as it stood on 2026-11-10: with
    both assignments removed it reports `release-validate: $(REL)` -- the
    target that genuinely could not run -- and `diff-releases: $(ARGS)`,
    which could, and whose emptiness is now declared rather than inferred.
    Two findings, one real, and the test cannot tell them apart; that is
    why the fix for the second was to state the intent in the Makefile
    instead of teaching this test a name to skip.
    """
    defaulted = _defaulted() | AUTOMATIC
    bad: list[str] = []
    for target, recipe in _recipes().items():
        body = "\n".join(recipe)
        # Variables this recipe defaults inline or refuses without.
        # `-n "$(X)"` is matched wherever it appears, not only directly
        # after `test`, because diff-releases guards two at once with
        # `test -n "$(A)" -a -n "$(B)"` and anchoring on `test` saw only
        # the first -- which is how this test first reported $(B).
        inline: set[str] = set()
        inline |= set(re.findall(r"\$\(or\s+\$\(([A-Z][A-Z0-9_]*)\)", body))
        inline |= set(re.findall(r"\$\(if\s+\$\(([A-Z][A-Z0-9_]*)\)", body))
        inline |= set(re.findall(r'-n\s+"\$\(([A-Z][A-Z0-9_]*)\)"', body))
        for m in VAR.finditer(body):
            name = m.group("name")
            if name in defaulted or name in inline:
                continue
            bad.append(f"{target}: $({name})")
    assert not bad, (
        "these recipes expand a variable with no default, no inline "
        "fallback and no guard, so a bare `make <target>` runs with it "
        "empty:\n  " + "\n  ".join(sorted(bad)))


def test_REL_specifically_has_a_default():
    """The instance, pinned separately from the rule above.

    The general test passes the moment `REL` gains a default OR gains a
    guard; this one records which was chosen, because a guard here would
    mean `make release-validate` still needs an argument and the docs that
    print it bare would still be wrong.
    """
    assert re.search(r"^REL\s*\?=\s*v0\.1\s*$", MAKEFILE.read_text(), re.M), (
        "REL must default to v0.1 so `make release-validate` runs bare")
