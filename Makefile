.PHONY: install lint-fast test test-all test-serial test-slowest lint gate release-validate check-results check-raw fetch-raw pathway-landscape clean

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

# FETCH THE RAW INPUTS, AND AUDIT WHAT IS ALREADY THERE.
#
# data/raw/* is gitignored and was populated by hand, so a fresh clone could
# not rebuild v0.1 and nothing recorded which UniProt release the proteins
# came from. check-raw is the audit and needs no network: it reports a file
# that disagrees with config/viruses.yaml, a built virus whose proteome is
# absent, and a file placed by hand with no URL, date or checksum.
# Two halves, because the inputs arrived two ways. fetch_raw --check
# validates the proteomes it fetched against the register. The host layer
# was placed by hand before provenance existed -- 329 MB feeding over
# three-quarters of the graph -- so record_raw_provenance re-hashes it
# against the sidecars written for it. Both run; the target fails if
# either does.
check-raw:
	python scripts/fetch_raw.py --check
	python scripts/record_raw_provenance.py --check

# Stamp a raw directory that was populated by hand. Writes no URL and no
# retrieval date -- neither is recoverable and a guess would be a
# fabricated record -- but records a checksum, so the next change to one
# of these files is visible.
record-provenance:
	python scripts/record_raw_provenance.py

fetch-raw:
	python scripts/fetch_raw.py --all-built

# REBUILD THE RELEASE FROM THE RAW INPUTS, IN DEPENDENCY ORDER.
#
# Not a list of commands here, because the order is not a straight line:
# selectivity, similarity and gtop read the ASSEMBLED release rather than the
# per-layer ones, so assemble_release runs in the middle and again at the end.
# Typed out by hand that goes wrong silently -- the layer still builds, from
# the PREVIOUS release, and nothing says so. src/kgav/pipeline.py holds the
# order as data; scripts/rebuild_release.py executes it.
#
# A layer whose raw inputs are absent is skipped and named. A layer whose
# inputs are present and which then fails stops the run, because the next
# assemble would write a release missing it while every count still read
# correctly.
rebuild:
	python scripts/rebuild_release.py

rebuild-plan:
	python scripts/rebuild_release.py --plan

# WHAT THE GRAPH IS BUILT FROM, in two timelines that are easy to conflate:
# when each raw input landed here, and what period the evidence covers. A
# file downloaded yesterday can carry evidence from 1998, so one date
# without the other misleads in opposite directions.
data-timeline:
	python scripts/data_timeline.py

# WHAT THE SCHEMA DECLARES AGAINST WHAT THE RELEASE HOLDS. assemble_release
# lists what IS there, so a declared class with no nodes prints no line and
# a reader counting lines sees a complete graph. This prints the absences,
# derived from the schema rather than a list kept here, and carries the
# recorded remedy for each gap so it is not worked out twice.
data-audit:
	python scripts/data_audit.py

# THE PATHWAY NULL, AS AN ARTIFACT THE AGENT CAN CITE. Without the results
# file the agent reports the pathway route as unmeasured, which is honest and
# weaker than the finding; with it, every ranking carries the measured
# degree-matched null. Re-run this whenever the release changes --
# src/kgav/landscape.py refuses a file measured on a different one rather
# than quoting it.
pathway-landscape:
	python scripts/pathway_landscape.py \
	  --out data/results/pathway_landscape.json

clean:
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
	rm -rf .pytest_cache .ruff_cache *.egg-info build dist
