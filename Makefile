.PHONY: install test test-serial test-slowest lint gate release-validate check-results clean

install:
	pip install -e ".[dev]"

# --dist loadfile keeps a file's tests on one worker, which matters because the
# release fixtures are session-scoped PER WORKER: the default per-test
# distribution can land test_agent_tools' two real-graph tests on two workers
# and parse the release twice. By file, the three release-reading files parse
# once each, in parallel.
test:
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
check-results:
	python scripts/check_results.py $(or $(RESULTS),data/results-corrected)

clean:
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
	rm -rf .pytest_cache .ruff_cache *.egg-info build dist
