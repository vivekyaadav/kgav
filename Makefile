.PHONY: kg-report kg-report-check install lint-fast test test-all test-serial test-slowest lint gate release-validate check-results check-raw fetch-raw pathway-landscape label-overlap clean

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
# Deliberately a two-rule subset, so `make test` pays milliseconds for it.
# `make lint` runs the whole declared set and is what CI runs.
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

# The rule set lives in pyproject.toml [tool.ruff.lint], with the reason for
# every rule left out. ruff is pinned in the dev extra: unpinned, the set CI
# enforced was whatever the newest ruff defaulted to, which is how this
# command went from 29 findings to 281 with no commit in between. CI runs
# this target rather than its own copy of the command.
lint:
	ruff check src tests scripts

# NEGATIVE CONTROL. data/releases/smoke carries four deliberately planted
# violations; this proves the validator finds THOSE, by code and by count.
#
# It used to be `if validator succeeds then fail else OK`, which accepted any
# non-zero exit -- an ImportError, a missing directory, a schema that would
# not load -- as proof the planted errors were caught. The one thing it could
# not distinguish is a working validator from a broken one.
#
# If this target fails after a schema change, read the diff it prints before
# editing the expectation: a code that stopped firing is the gate working.
gate:
	python scripts/validate_release.py data/releases/smoke \
		--expect "CONSTRAINT=1,MISSING_PROVENANCE=1,MISSING_QUALIFIER=2"

# REL HAD NO DEFAULT AND THIS TARGET THEREFORE NEVER RAN (fix-list L14's
# family). `make release-validate` expanded to `data/releases/` and exited 2
# with "no nodes.parquet or nodes.jsonl in data/releases" on every invocation
# since the target was written. The schema check it performs was real -- it
# runs inside `make rebuild` on every layer -- so nothing was unvalidated;
# what was broken was the one command a reader is told to use to check the
# release by hand. Found 2026-11-10 by running it.
#
# tests/test_makefile_targets_resolve.py now asserts that every variable a
# recipe expands either carries a default here, is wrapped in $(or)/$(if), or
# is guarded by a `test -n` in its own recipe. A target that cannot run is
# the same defect class as a check that cannot fail.
REL ?= v0.1

# DECLARED OPTIONAL, not exempted in the test. `diff-releases` passes
# $(ARGS) straight through to the script and an empty expansion is the
# intended default, but an undeclared variable and a deliberately empty one
# look identical to a reader and to a checker. Assigning it here states
# which it is, so the rule in tests/test_makefile_targets_resolve.py stays
# uniform and needs no list of names to ignore -- a list that would be the
# same defect as a check written with the literals it is checking.
ARGS ?=

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

# WHAT THE max_phase GATE COSTS (fix-list B1). Needs the 30 GB chembl_37
# sqlite, so it is a diagnostic and not a test; ~20s.
diagnose-phase-gate:
	python scripts/diagnose_phase_gate.py

# WHAT RESTATING THE TASK WOULD COST (fix-list B3). Runs on the release
# alone in ~25s; add DB=... for the drug_mechanism definition, which is the
# literal reading of "known mechanism" and needs the gitignored sqlite.
# L19's unrun query: does ONE layer file two tautomers under one InChIKey?
# The 32 conflations L19 found came from a CROSS-layer diagnostic, and the
# only build with two layers per compound was B1, which was declined. So on
# the shipped release that diagnostic reports nothing, and whether a single
# layer does the same thing decides whether L19 is a limit to STATE or a
# defect to FIX. LAYER= selects ncats (default), chembl or hosttargets;
# the chembl layers need DB=.
diagnose-tautomer-conflation:
	python scripts/diagnose_tautomer_conflation.py \
		--layer $(or $(LAYER),ncats) $(if $(DB),--db $(DB),)

diagnose-label-reachability:
	python scripts/diagnose_label_reachability.py $(if $(DB),--db $(DB),)

# The design facts docs/prereg_structure_scorer.md is registered against.
# SOURCES defaults to the registered population, ChEMBL only; pass
# SOURCES=infores:ncats-opendata for the matched secondary cell, or
# SOURCES=- for the pooled task the registration deliberately excludes.
# DEFAULTING TO THE REGISTERED POPULATION, NOT TO POOLED. The first version
# of this target produced an empty --sources when SOURCES was unset, so the
# default silently ran the pooled task -- the one population the
# registration excludes by name, because a structure scorer can separate
# ChEMBL antivirals from NCATS library compounds on collection alone.
# The registered run. Writes data/results/structure_scorer.json.
structure-scorer:
	python scripts/structure_scorer.py $(ARGS)

# The symmetric-control follow-up, docs/prereg_structure_scorer_symmetric.md.
# A SEPARATE TARGET AND A SEPARATE FILE: the publication-disjoint run is its
# comparator, its figures are quoted in a verdict and pinned by tests, and
# overwriting it would destroy the only thing it is being compared against.
structure-scorer-symmetric:
	python scripts/structure_scorer.py --control knn $(ARGS)

# The strength-matched follow-up,
# docs/prereg_structure_scorer_strength_matched.md. It DERIVES every k from
# data/results/structure_scorer.json and refuses if that file is absent,
# because there is then no strength to match.
structure-scorer-strength:
	python scripts/structure_scorer.py --control strength $(ARGS)

# The read-out: what the graph suggests, by name, for each virus. Not a
# protocol -- it points the validated scorer at the question.
repurposing-candidates:
	python scripts/repurposing_candidates.py $(ARGS)

STRUCTURE_SOURCES = $(or $(SOURCES),infores:chembl)
diagnose-structure-design:
	python scripts/diagnose_structure_scorer_design.py \
	  $(if $(filter-out -,$(STRUCTURE_SOURCES)),--sources $(STRUCTURE_SOURCES),)

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

# WHAT HAS BEEN INGESTED, BUILT, AND IS STILL MISSING, per layer and per
# virus. docs/kg_audit_report.txt has cited this as `make audit` in three
# places since it was written, and there was no such target: the script
# existed, the name did not, so the one line telling a reader how to
# check the per-virus tiers named a command that would not run. Found on
# 2026-11-02 by tests/test_kg_audit_report.py, which now checks every
# make target that report cites against this file.
audit:
	python scripts/audit_build.py

# THE PATHWAY NULL, AS AN ARTIFACT THE AGENT CAN CITE. Without the results
# file the agent reports the pathway route as unmeasured, which is honest and
# weaker than the finding; with it, every ranking carries the measured
# degree-matched null. Re-run this whenever the release changes --
# src/kgav/landscape.py refuses a file measured on a different one rather
# than quoting it.
# THE REPORT'S RELEASE CENSUS. Section 4 of docs/kg_audit_report.txt is
# generated from the assembled release and stamped with its content
# hashes, so a rebuild without a regenerate is detectable rather than
# merely wrong. kg-report-check is the read-only half, for the gate.
kg-report:
	python scripts/write_kg_report.py --placeholder

kg-report-check:
	python scripts/write_kg_report.py --check

pathway-landscape:
	python scripts/pathway_landscape.py \
	  --out data/results/pathway_landscape.json

# THE SHARED-PMID CONFOUND, AS AN ARTIFACT RATHER THAN CONSOLE OUTPUT. The
# 92.2% figure is why docs/evaluation.md calls its own M1 result partly
# self-scoring, and it lived only in this script's stdout while being
# quoted in seven places across two documents. Written here, it is stamped
# against the release and audited by check-results like everything else.
label-overlap:
	python scripts/diagnose_label_overlap.py \
	  --out data/results/label_overlap.json

# WHAT THE AGENT SAYS, for one virus from EACH tier in section 9 of the
# audit report -- not just the one with the most data. Prints for a
# person to read; there is nothing to assert. The first run of it found
# that the refusal written for Yellow fever was never reached, which
# every unit test had missed because a unit test asks the question its
# author thought to ask.
tier-check:
	bash scripts/tier_spotcheck.sh

# REGENERATE EVERY RESULTS FILE, from the block in docs/evaluation.md
# rather than from a copy of it. The order is load bearing and was
# being transcribed by hand after each rebuild, where a dropped flag
# changes a protocol without failing.
reproduce:
	python scripts/run_reproduce_block.py

reproduce-print:
	python scripts/run_reproduce_block.py --print

# TWO RELEASES, EDGE BY EDGE. edges_sha says they differ; this says
# whether the graph differs or only the byte order.
#   make diff-releases A=~/kgav-build5 B=data/releases/v0.1
diff-releases:
	@test -n "$(A)" -a -n "$(B)" || { echo "set A=<dir> B=<dir>"; exit 2; }
	python scripts/diff_releases.py "$(A)" "$(B)" $(ARGS)

# DID A REBUILD MOVE ANY FIGURE, as opposed to just the release stamp.
# BEFORE= is a copy of data/results taken before reassembling; make one
# OUTSIDE the repo first, or an untracked tree of results sits where
# `git add -A` will find it.
#   make compare-results BEFORE=~/kgav-results.pre-dedup
compare-results:
	@test -n "$(BEFORE)" || { echo "set BEFORE=<dir>"; exit 2; }
	python scripts/compare_results.py "$(BEFORE)" data/results --fields

clean:
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
	rm -rf .pytest_cache .ruff_cache *.egg-info build dist
