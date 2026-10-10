.PHONY: install lint-fast test test-all test-serial test-slowest lint gate release-validate check-results clean

install:
	pip install -e ".[dev]"

# --dist loadfile keeps a file's tests on one worker, which matters because the
# release fixtures are session-scoped PER WORKER: the default per-test
# distribution can land test_agent_tools' two real-graph tests on two workers
# and parse the release twice. By file, the three release-reading files parse
# once each, in parallel.
# F821 is undefined-name, F811 redefinition. They cost milliseconds and catch
# the class of bug the conftest refactor introduced: two tests had their
# BODIES rewritten to use a fixture while their SIGNATURES were not, so they
# referenced an undefined name. A NameError only fires when the test actually
# RUNS, so it was invisible on a machine where the release is absent and those
# tests skip. `test` depends on this so the combination cannot pass again.
#
# Deliberately NOT the full ruff default set: that currently reports 29
# findings, nearly all import ordering, and a backlog of style warnings in
# front of a correctness gate is a gate nobody runs. `make lint` still runs
# everything.
lint-fast:
	ruff check --select F821,F811 src tests scripts

test: lint-fast
	pytest -q -n auto --dist loadfile -m "not slow"

# Everything, including the raw-source integration checks.
test-all:
	pytest -q -n auto --dist loadfile

# Serial. Use when a failure needs a readable traceback, or to check that a
# failure is real rather than two workers colliding -- xdist runs each test in
# a separate PROCESS, so the code under test needs no thread-safety, but a
# test that writes outside tmp_path would still race.
test-serial:
	pytest -q

# Where the time actually goes. The suite went 8s without a release on disk
# and 85s with one, because four files each re-read the 102 MB edges.jsonl.
test-slowest:
	pytest -q --durations=15

lint:
	ruff check src tests scripts

gate:
	@if python scripts/validate_release.py data/releases/smoke; then \
		echo "FAIL: validator accepted planted errors"; exit 1; \
	else \
		echo "OK: validator rejected planted errors"; \
	fi

release-validate:
	python scripts/validate_release.py data/releases/$(REL)

# Flags results computed from a release that has since changed. The stale
# "before" column that made a falling positive count read as a rise was
# undetectable because results carried no provenance.
# AUDIT WHAT THE SCRIPTS WRITE. This pointed at data/results-corrected, a
# snapshot from 2026-09-17 that nothing writes to and only this target reads,
# while every script defaults to data/results. So the audit reported 11 stale
# files in a directory nobody uses and said nothing about the one in use --
# which then held one current file, two stale, and eight with no provenance
# block at all. An integrity check aimed at the wrong directory is worse than
# none: it returns a verdict about neither.
#
# data/results-corrected is superseded. Two directories with identical
# filenames, one audited and one used, is the hazard itself.
check-results:
	python scripts/check_results.py $(or $(RESULTS),data/results)

clean:
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
	rm -rf .pytest_cache .ruff_cache *.egg-info build dist
