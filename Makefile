.PHONY: agent-ui agent-ask auc-contribution hit-rate identity-fixed-points auc-delta kg-report kg-report-check scorer-table scorer-table-check clinical-labels degree-sensitivity landscape-count-move temporal-structure-scorer install lint-fast test test-all test-serial test-slowest lint gate release-validate check-results check-raw fetch-raw pathway-landscape label-overlap clean

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

# DID THE SCORES GO DOWN. check-results says whether a results file is
# current; this says whether the rebuild that made it current cost anything.
# Defaults to the previous commit and to a floor of one AUC point on cells
# both releases REPORT -- the scorer's own trustworthy / hits_determined /
# readable flags decide that, because the material falls a rebuild produces
# are in cells with six negatives, which the project already refuses.
#   make auc-delta ARGS="--against 5151de8"
auc-delta:
	python scripts/compare_results_aucs.py --fail-below -0.01 $(ARGS)

# FIX-LIST L27's closing evidence: is every node id what the normaliser
# derives from that node's own stored SMILES? The whole release is an RDKit
# pass per molecule -- minutes, not seconds -- so it lives here rather than
# in the suite, which samples it. See the script's header for the trade.
identity-fixed-points:
	python scripts/diagnose_identity_fixed_points.py $(ARGS)

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
# THE CHAT FRONTEND, and the agent behind it. Loads the release once and
# holds it, so a question costs what the question costs rather than the
# ~5s reload `kgav_ask.py` pays every time it starts.
#
# LOOPBACK ONLY. scripts/agent_server.py defaults --bind to 127.0.0.1
# because the payloads name compounds, targets and unpublished AUCs, and
# the model is reached at localhost:11434. Pass BIND=0.0.0.0 if you mean
# to expose it.
#
# STDLIB ONLY: http.server plus one HTML file. There is no npm install
# step to forget, which is the whole reason the previous frontend --
# 472MB of node_modules for a page that made no network call -- is gone.
PORT ?= 8765
BIND ?= 127.0.0.1
agent-ui:
	python scripts/agent_server.py --port $(PORT) --bind $(BIND)

# The same agent from the command line, for a one-off question. Q is the
# question; everything after it goes to the script, so
# `make agent-ask Q="..." ARGS=--no-model` skips synthesis.
agent-ask:
	@test -n "$(Q)" || { echo 'usage: make agent-ask Q="your question"'; exit 2; }
	python scripts/kgav_ask.py $(ARGS) "$(Q)"

# WHERE THE AUC COMES FROM, AS A SHARE PER APPROACH. Reads every figure
# from data/results and computes no AUC of its own. Answers the question
# "percentage-wise, what did each approach contribute" without pretending
# the approaches stack: section 2 is the one genuine decomposition,
# because degree and structure rank the same candidates in the same cell.
#
# NOT IN THE REPRODUCE BLOCK, and not written to data/results: it derives
# nothing new, so there is nothing for `make check-results` to find stale.
# Pass --out to keep a copy.
auc-contribution:
	python scripts/auc_contribution.py

# THE WET-LAB HIT RATE: of the top k ranked compounds, how many came back
# active. AUC does not answer this -- it is a statement about ordering,
# and a screen only ever orders the top of the list.
#
# POST-HOC AND NOT A REGISTERED STATISTIC. The temporal registration
# fixed AUC before any AUC existed; precision@k is computed here after
# those AUCs were read, on the same held-out labels. Descriptive only,
# p-values uncorrected for selection, and deliberately NOT written to
# data/results, where it would read as a claim-carrying result.
#
# Reports both registered controls prospectively and all three exclusion
# strengths on the plate cells, because the figure that matters is how
# far the precision falls when a candidate's nearest known actives are
# withheld.
hit-rate:
	python scripts/hit_rate_at_k.py

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

# EVERY SCORER AUC, EVERY VIRUS, SIDE BY SIDE. The three studies write
# three results files with three shapes, so the comparison that decides
# whether any of this is about chemistry -- each protocol's AUC against
# its own degree baseline, per virus -- was being reassembled by hand in
# a terminal after every rebuild and written down nowhere. Generated and
# stamped, like the census: a table built from a mixture of releases is
# detectable rather than merely wrong. Runs as the last step of
# `make reproduce`, so the numbers are printed every time they move.
# DOES THE DWPC WEIGHTING READ THE ANSWER SHEET. Three arms over every
# channel: as published, with the held-out labels out of the step weight
# only, and out of the hub cut as well. docs/evaluation.md's claim
# section quotes the middle arm beside its own M7 figure, so this ships
# as a stamped artifact rather than a terminal printout -- a figure in
# that document has to be derivable by a test. It REGISTERS NOTHING:
# switching the weighting is a decision that needs its own registration.
#
# MEASURED 5m18s, and nearly all of it is M4 walking ~3,950 drugs three
# times -- M1+M6+M7 together are 8.6s. Running only the channels that
# CAN move was considered and rejected: the five that never move are the
# control, and they are what licenses the claim "4 of the 61 cells move
# and none crosses alpha". Dropping them would save five minutes and the
# ability to say the finding is small.
degree-sensitivity:
	python scripts/diagnose_label_degree_weighting.py \
	  --out data/results/degree_weighting_sensitivity.json $(ARGS)

# WHY TWO LANDSCAPE CELLS MOVED, closing the question left open in
# docs/prereg_pathway_overlap.md's note of 2026-11-17. MEASURED 6.7s: it
# reads two recorded artifacts and builds the factor sets once. It needs
# GIT, not a rebuild -- the earlier artifact comes from a revision, which
# is the whole reason the controlled `ncats` re-run that note named is
# unnecessary.
#
# IT IS IN THE REPRODUCE BLOCK, and the first draft of this comment
# argued it should not be -- that comparing against a historical release
# makes it a one-off. That was wrong, and the guard caught it: anything
# under data/results is globbed by check-results and will be reported
# stale with nothing able to refresh it, which is the state
# tests/test_reproduce_block.py exists to prevent. The set-difference
# mechanism also works across any two artifacts however many releases
# apart, so the comparison stays valid as the before-side ages.
#
# IT EXITS 1 when the operand difference does not account for every move,
# or when any cell's factor set changed size. That breaks the block on
# purpose: it means an ingest moved a landscape count and nothing has
# explained it yet.
landscape-count-move:
	python scripts/diagnose_landscape_count_move.py \
	  --out data/results/landscape_count_move.json $(ARGS)

# THE ONE INDEPENDENT LABEL AXIS, registered in advance by
# docs/prereg_clinical_labels.md. One cell, family of one, alpha 0.05.
# It reads the release and the three scorer files read nothing from it;
# it is in the reproduce block before scorer-table because
# check-results audits data/results and nothing else regenerates this.
# DOES THE SCORER PREDICT FORWARD. Registered in advance by
# docs/prereg_temporal_structure_scorer.md, signed 2026-11-19 before any
# AUC against a post-2021 label existed. Pre-2021 actives are the
# reference set, post-2021 labels the evaluated population; the
# statistic, fingerprint and interval are imported from
# structure_scorer.py rather than copied.
#
# It reads data/releases/v0.1-train2021, which `make reproduce` builds in
# step 1, so it must come after the temporal splits and does not depend
# on the scorer protocols. It REFUSES -- writing a refusal record, not
# numbers -- if the split is stale against the release, if the evaluable
# cell count is not the registered 10, or if fingerprint coverage drops
# below 95%.
temporal-structure-scorer:
	python scripts/temporal_structure_scorer.py \
	  --out data/results/temporal_structure_scorer.json $(ARGS)

clinical-labels:
	python scripts/clinical_labels.py $(ARGS)

scorer-table:
	python scripts/write_scorer_table.py

scorer-table-check:
	python scripts/write_scorer_table.py --check

# FIX-LIST P10: the marked tables, written from the derivations their own
# tests already compute. --check reports drift and exits 1 without writing.
rewrite-tables:
	python scripts/rewrite_derived_tables.py $(ARGS)

rewrite-tables-check:
	python scripts/rewrite_derived_tables.py --check

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
