# 12 — Reproducing everything

```
make install       # pip install -e ".[dev]"
make test          # the gate: lint subset + the suite, slow tests deselected
make lint          # the full declared ruff set
make gate          # the validator's NEGATIVE CONTROL
make check-results # is every artifact current for the release on disk
make reproduce     # regenerate every results file, in order
```

## The order is read from the document, not transcribed

`make reproduce` runs the Reproducing block in `doc2/08-results.md` **as the
document holds it**. The order is load-bearing and was being transcribed by
hand after each rebuild, where a dropped flag changes a protocol without
failing. `make reproduce-print` prints the steps without running them, and
`tests/test_reproduce_block.py` fails if an audited results file has no step
that regenerates it — because anything under `data/results` is globbed by
`check-results` and would be reported stale with nothing able to refresh it.

## Every target, and the command it runs

57 targets. The **recipe** is shown rather than a
description, because a Makefile comment block is not delimited per target and
the first version of this table put a sentence about one target beside
another.

| target | the command it runs |
|:---|:---|
| `make install` | `pip install -e ".[dev]"` |
| `make lint-fast` | `ruff check --select F821,F811 src tests scripts` |
| `make test` | `pytest -q -n auto --dist loadfile -m "not slow"` |
| `make test-all` | `pytest -q -n auto --dist loadfile` |
| `make test-serial` | `pytest -q` |
| `make test-slowest` | `pytest -q --durations=15` |
| `make lint` | `ruff check src tests scripts` |
| `make gate` | `python scripts/validate_release.py data/releases/smoke \`<br>`--expect "CONSTRAINT=1,MISSING_PROVENANCE=1,MISSING_QUALIFIER=2"` |
| `make release-validate` | `python scripts/validate_release.py data/releases/$(REL)` |
| `make check-results` | `python scripts/check_results.py $(or $(RESULTS),data/results)` |
| `make auc-delta` | `python scripts/compare_results_aucs.py --fail-below -0.01 $(ARGS)` |
| `make identity-fixed-points` | `python scripts/diagnose_identity_fixed_points.py $(ARGS)` |
| `make diagnose-phase-gate` | `python scripts/diagnose_phase_gate.py` |
| `make diagnose-tautomer-conflation` | `python scripts/diagnose_tautomer_conflation.py \`<br>`--layer $(or $(LAYER),ncats) $(if $(DB),--db $(DB),)` |
| `make diagnose-label-reachability` | `python scripts/diagnose_label_reachability.py $(if $(DB),--db $(DB),)` |
| `make structure-scorer` | `python scripts/structure_scorer.py $(ARGS)` |
| `make structure-scorer-symmetric` | `python scripts/structure_scorer.py --control knn $(ARGS)` |
| `make structure-scorer-strength` | `python scripts/structure_scorer.py --control strength $(ARGS)` |
| `make repurposing-candidates` | `python scripts/repurposing_candidates.py $(ARGS)` |
| `make diagnose-structure-design` | `python scripts/diagnose_structure_scorer_design.py \`<br>`$(if $(filter-out -,$(STRUCTURE_SOURCES)),--sources $(STRUCTURE_SOURCES),)` |
| `make check-raw` | `python scripts/fetch_raw.py --check`<br>`python scripts/record_raw_provenance.py --check` |
| `make record-provenance` | `python scripts/record_raw_provenance.py` |
| `make fetch-raw` | `python scripts/fetch_raw.py --all-built` |
| `make rebuild` | `python scripts/rebuild_release.py` |
| `make rebuild-plan` | `python scripts/rebuild_release.py --plan` |
| `make data-timeline` | `python scripts/data_timeline.py` |
| `make agent-ui` | `python scripts/agent_server.py --port $(PORT) --bind $(BIND)` |
| `make agent-ask` | `python scripts/kgav_ask.py $(ARGS) "$(Q)"` |
| `make auc-contribution` | `python scripts/auc_contribution.py` |
| `make hit-rate` | `python scripts/hit_rate_at_k.py` |
| `make screen-power` | `python scripts/screen_power.py` |
| `make metagraph-figure` | `python scripts/write_metagraph_figure.py --render` |
| `make metagraph-figure-check` | `python scripts/write_metagraph_figure.py --check` |
| `make arch-figure` | `python scripts/write_architecture_figure.py --render` |
| `make arch-figure-check` | `python scripts/write_architecture_figure.py --check` |
| `make protocol-disagreement` | `python scripts/diagnose_protocol_disagreement.py` |
| `make prospective-analogy` | `python scripts/diagnose_prospective_analogy.py` |
| `make data-audit` | `python scripts/data_audit.py` |
| `make doc2` | `python scripts/write_doc2.py` |
| `make doc2-check` | `python scripts/write_doc2.py --check` |
| `make audit` | `python scripts/audit_build.py` |
| `make kg-report` | `python scripts/write_kg_report.py --placeholder` |
| `make kg-report-check` | `python scripts/write_kg_report.py --check` |
| `make degree-sensitivity` | `python scripts/diagnose_label_degree_weighting.py \`<br>`--out data/results/degree_weighting_sensitivity.json $(ARGS)` |
| `make landscape-count-move` | `python scripts/diagnose_landscape_count_move.py \`<br>`--out data/results/landscape_count_move.json $(ARGS)` |
| `make temporal-structure-scorer` | `python scripts/temporal_structure_scorer.py \`<br>`--out data/results/temporal_structure_scorer.json $(ARGS)` |
| `make clinical-labels` | `python scripts/clinical_labels.py $(ARGS)` |
| `make scorer-table` | `python scripts/write_scorer_table.py` |
| `make scorer-table-check` | `python scripts/write_scorer_table.py --check` |
| `make pathway-landscape` | `python scripts/pathway_landscape.py \`<br>`--out data/results/pathway_landscape.json` |
| `make label-overlap` | `python scripts/diagnose_label_overlap.py \`<br>`--out data/results/label_overlap.json` |
| `make tier-check` | `bash scripts/tier_spotcheck.sh` |
| `make reproduce` | `python scripts/run_reproduce_block.py` |
| `make reproduce-print` | `python scripts/run_reproduce_block.py --print` |
| `make diff-releases` | `python scripts/diff_releases.py "$(A)" "$(B)" $(ARGS)` |
| `make compare-results` | `python scripts/compare_results.py "$(BEFORE)" data/results --fields` |
| `make clean` | `find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null \|\| true`<br>`rm -rf .pytest_cache .ruff_cache *.egg-info build dist` |

## What each artifact is, and what regenerates it

| file | bytes | kind | edges sha | generated |
|:---|---:|:---|:---|:---|
| `baselines_NCBITaxon_11053.json` | 22,298 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_11060.json` | 66,018 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_11069.json` | 22,297 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_11070.json` | 21,542 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_11072.json` | 22,180 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_11082.json` | 66,111 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_11084.json` | 27,471 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_11089.json` | 21,608 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_11137.json` | 66,081 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_1335626.json` | 66,085 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_2697049.json` | 66,026 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_277944.json` | 65,992 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_290028.json` | 38,294 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_31631.json` | 66,129 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_64320.json` | 65,995 | result | 873623c9 | 2026-11-17 |
| `baselines_NCBITaxon_694009.json` | 43,377 | result | 873623c9 | 2026-11-17 |
| `clinical_labels.json` | 3,428 | result | 873623c9 | 2026-11-17 |
| `degree_weighting_sensitivity.json` | 102,617 | result | 873623c9 | 2026-11-17 |
| `hard_negatives.json` | 13,676 | result | 873623c9 | 2026-11-17 |
| `hard_negatives_row1.json` | 67,443 | result | 873623c9 | 2026-11-17 |
| `hard_negatives_row2.json` | 70,334 | result | 873623c9 | 2026-11-17 |
| `hard_negatives_row3.json` | 52,538 | result | 873623c9 | 2026-11-17 |
| `hard_negatives_row4.json` | 52,204 | result | 873623c9 | 2026-11-17 |
| `label_overlap.json` | 699 | result | 873623c9 | 2026-11-17 |
| `landscape_count_move.json` | 9,843 | result | 873623c9 | 2026-11-18 |
| `metapath_comparison.json` | 68,897 | result | 873623c9 | 2026-11-17 |
| `pathway_landscape.json` | 64,350 | result | 873623c9 | 2026-11-17 |
| `repurposing_candidates.json` | 394,933 | result | 873623c9 | 2026-11-21 |
| `structure_scorer.json` | 83,992 | result | 873623c9 | 2026-11-17 |
| `structure_scorer_strength_matched.json` | 63,243 | result | 873623c9 | 2026-11-17 |
| `structure_scorer_symmetric.json` | 131,970 | result | 873623c9 | 2026-11-17 |
| `temporal_2021_computed_only.json` | 17,319 | result | 873623c9 | 2026-11-17 |
| `temporal_2021_exclude.json` | 17,221 | result | 873623c9 | 2026-11-17 |
| `temporal_2021_include.json` | 39,232 | result | 873623c9 | 2026-11-17 |
| `temporal_2021_include_overlap-keep.json` | 2,798 | REFUSED | 873623c9 | 2026-11-17 |
| `temporal_structure_scorer.json` | 36,538 | result | 873623c9 | 2026-11-19 |

1 of 36 is a **recorded refusal** rather
than a result. A run that refuses still has to leave a record: the protocol
aborts, writes a stamped file carrying the audit failures and **no results
block**, and `check_results` reports it as `refused` rather than `stale`.
Before that, a refusal left the previous run's numbers on disk — computed from
a graph three revisions old — with nothing saying the protocol now refuses to
run, and `check_results` advised "regenerate it", which was impossible.

## The four statuses, and why they are four

| status | meaning | what to do |
|:---|:---|:---|
| `ok` | the stamp matches the release on disk | nothing |
| `stale` | computed from a graph that has since moved | re-run it |
| `unstamped` | no provenance block at all | its numbers cannot be attributed to any graph |
| `refused` | the protocol declined to run, and said why | resolve the audit failure |
| `absent` | (derived releases and layers only) never built | build it |

A **missing** release is not a match. `_sha256` returns `None` for a file that
does not exist, and two `None`s compare equal — so checking against a typo'd
or deleted release path reported every file as current. That is a false
"current" from the one tool whose job is to catch false "current" readings,
which is worse than no check at all.

## Timings, measured on this workstation

| command | time | note |
|:---|---:|:---|
| `make test` | ~2 min | 3,054 tests, slow deselected, `-n auto --dist loadfile` |
| `make lint` | seconds | ruff, pinned in the dev extra |
| `make check-results` | ~10 s | hashes every release file once |
| `make protocol-disagreement` | ~3 s | |
| `make prospective-analogy` | ~10 s | |
| `make hit-rate` | ~12 s | |
| `make screen-power` | ~2 min | power enumerated exactly, not approximated |
| `make degree-sensitivity` | ~5 min | nearly all of it M4 walking ~3,950 drugs three times |
| `make repurposing-candidates` | ~8 min | all cells; one virus is ~100 s |
| `make rebuild` | ~25 min | ChEMBL and NCATS are the bulk of it |

`--dist loadfile` keeps a file's tests on one worker, which matters because
the release fixtures are session-scoped **per worker**: the default per-test
distribution lands the real-graph tests on several workers and parses the
102 MB `edges.jsonl` repeatedly.

## The gates that are not tests

| gate | what it would catch |
|:---|:---|
| `make gate` | the validator failing to find four planted violations, **by code and by count** |
| `make auc-delta` | a rebuild dropping any claim-carrying AUC by more than one point |
| `make check-results` | a figure quoted from a graph that has moved |
| `make check-raw` | a raw input that disagrees with the register, or was placed by hand with no checksum |
| `make kg-report-check` | the release census drifting from the release |
| `make scorer-table-check` | the published AUC table drifting from the results files |
| `make doc2-check` | any document in `doc2/` drifting from the code or the artifacts |

| `make metagraph-figure-check` | the metagraph figure losing or inventing a predicate |
| `make doc2-check` | **this folder** drifting from the code |

`make auc-delta` defaults to the previous commit and to a floor of one AUC
point, on cells both releases **report** — the scorer's own
`trustworthy` / `hits_determined` / `readable` flags decide that, because the
material falls a rebuild produces are in cells with six negatives, which the
project already refuses.

## Rebuilding the graph

```
make check-raw      # audit what is on disk, no network
make fetch-raw      # fetch the proteomes
make rebuild-plan   # what would run, with each step's verdict
make rebuild        # run it
```

A layer whose raw inputs are absent is **skipped and named**. A layer whose
inputs are present and which then fails **stops the run**, because the next
assemble would write a release missing it while every count still read
correctly. `rebuild-plan` reports a third state — a skipped layer whose
directory already exists is **reused**, not absent, and that is the more
dangerous of the two.

One environment variable is required by one layer: `ORCS_KEY`, exactly 32
characters, for BioGRID ORCS. It lives in `.env`, is never committed, and
`ingest_orcs` exits before any request without it — so a run without the key
looks like a layer that produced nothing, which is why `pipeline.availability`
reports the reason by name.

---

*The block below was `doc2/08-results.md`'s and moved here on
2026-11-21 when that document was retired. It is read and run by
`scripts/run_reproduce_block.py`, which extracts the fenced `bash`
section under the heading — so this is executable text, not an
illustration, and editing it changes what `make reproduce` does.*

## Reproducing

**Order matters, for two reasons.** `temporal_split.py` cuts the training
graph that `hard_negatives.py` reads, so the split comes first — and the
three `--undated` protocols write to *one* output directory, so the last one
run is the split that stays on disk. Every `hard_negatives.py` protocol
likewise writes to one `hard_negatives.json`, so the last one run is the row
that stays. Both lists below therefore end on the run whose numbers the rest
of doc2/08-results.md quotes.

```bash
# ---- after any rebuild of data/releases/v0.1, everything below is required.
# Each block regenerates files that `make check-results` audits; a block
# skipped leaves a file describing a graph that no longer exists.

# 1. the temporal split, which is also an input to step 3
python scripts/temporal_split.py --cutoff 2021 --undated exclude
python scripts/temporal_split.py --cutoff 2021 --undated computed_only
python scripts/temporal_split.py --cutoff 2021 --undated include   # LAST: the
                                                                   # split the
                                                                   # rest reads

# control: must refuse, naming the already-known labels -- 51 of them on
# this release, 12 positive and 38 negative for SARS-CoV-2 and 1 negative
# for SARS-CoV. It said "the 19 already-known compounds", which was the
# figure before the NCATS entry assays enlarged the SARS-CoV-2 label set
# and gave SARS-CoV its first overlap; nothing could fail on it, because
# the refusal record it describes was never read by a test. One is read
# now. Writes a
# refusal record and exits 1 -- expected, and it writes no split, so it
# cannot disturb the one above.
python scripts/temporal_split.py --cutoff 2021 --undated include \
    --overlap-policy keep || true

# 2. per-virus baselines -- one file each, every taxon in
#    config/viruses.yaml. Six of the eight flaviviruses promoted on
#    2026-10-23 have no HOST_FACTOR_FOR edge, so their host-mediated
#    channels are empty by construction; the file is still written, and an
#    empty channel recorded is what lets a later ingest show up as a change.
#
#    ONE INVOCATION, not sixteen. dwpc_scores computes
#    {metapath: {drug: {virus}}} for EVERY virus at once and the script then
#    used one slice, so a loop recomputed the same traversal over 25,383
#    drugs sixteen times -- about 24 of this block's 60 minutes. --all-viruses
#    reads the register, so a virus promoted later is evaluated without
#    anyone editing this list, which went stale once already when it was
#    eight taxa long.
python scripts/run_baselines.py --all-viruses

# 3. metapath comparison, and the hard-negative protocols
python scripts/compare_metapaths.py

# rows 1-5, IN ORDER, in ONE process.
#
# Each --protocol string is parsed by hard_negatives.py's own parser, so
# nothing here could not be run on its own: copy any one of them onto the
# command line and it behaves exactly as it did when these were five
# commands. The results file is rewritten after each protocol, so the LAST
# one is still the row that stays on disk -- row 5, the protocol the rest
# of doc2/08-results.md quotes. It used to run first, which meant the no-op
# control overwrote it and hard_negatives.json held row 2 while being read
# as row 5.
#
# ONE PROCESS because dwpc_scores depends on the GRAPH alone. These five
# differ only in their LABELS -- which compounds count as active, which as
# negative, which screen they came from -- so run separately they recompute
# the same traversal: five protocols over the full graph and the training
# graph is ten DWPC computations of which four are distinct. One invocation
# measured 2m41s on v0.1, so the five were about thirteen and a half
# minutes of which roughly six were duplicated work.
#
# The no-op control is NOT in this list. It used to be a sixth run, and it
# verified nothing: every protocol writes to one hard_negatives.json, so
# row 2's file had been overwritten two steps before the control ran and
# the control's own output was overwritten one step after. It cost ~2.5
# minutes a rebuild and nobody could read the comparison it was supposed to
# make. It runs in tests/test_hard_negatives.py instead, on every commit,
# with a second arm proving the fixture could have shown the two settings
# disagreeing.
#
# EACH ROW WRITES ITS OWN FILE. They all wrote hard_negatives.json, so only
# row 5 survived and rows 1-4 left no record at all -- which is how the
# table in "## The sequence" came to report four rows the graph had stopped
# producing, with check-results certifying the one file it could see. Row 5
# keeps the canonical name, because calibrate.py and the agent read it.
#
# ROW 5 NAMES A PROTOCOL AS WELL AS A SCREEN, since 2026-11-17.
# --label-source alone stopped meaning "one screen" when
# infores:ncats-opendata began publishing pseudotyped-entry assays
# alongside the cytopathic-effect one: build_labels groups by virus, so
# the SARS-CoV-2 cell pooled both protocols and this row's headline
# figures were computed over 149 positives from two concentration
# ranges rather than 89 from one. --label-assay-type restores what the
# row was always described as measuring.
python scripts/hard_negatives.py \
    --protocol "--selectivity selective-only --results-name hard_negatives_row1.json" \
    --protocol "--selectivity all --results-name hard_negatives_row2.json" \
    --protocol "--selectivity verified-only --selectivity-applies both --results-name hard_negatives_row3.json" \
    --protocol "--selectivity verified-only --selectivity-applies both --publication-disjoint --results-name hard_negatives_row4.json" \
    --protocol "--selectivity all --publication-disjoint --label-source infores:ncats-opendata --label-assay-type cell_based_antiviral"
python scripts/calibrate.py          # which channels were evaluable, and why not

# 4. the pathway landscape. It depends on nothing above -- it reads the
#    release directly and scores no compounds -- but check-results globs
#    data/results and audits it like everything else, so leaving it out
#    meant a reader could run this block top to bottom and still be told
#    a file was stale. It was missing here from the day the file was added
#    until 2026-10-28; it looked fine only because it was being rebuilt by
#    hand. tests/test_reproduce_block.py now fails on a results file that
#    no step here regenerates.
make pathway-landscape

#    WHY ANY LANDSCAPE CELL MOVED since the pre-entry-assay release. It
#    reads the artifact just written and the same artifact at an earlier
#    revision out of GIT, so it needs a repository rather than a rebuild.
#    `n_observed` is `len(targets & factors)`, so a cell can move with no
#    edge in its own channel changing -- which is how two moved cells came
#    to be recorded with a cause that could not have been theirs
#    (doc2/registrations/prereg_pathway_overlap.md, amendment 2026-11-18). IT EXITS 1
#    when the operand difference does not account for every move, or when
#    any cell's factor set changed size, and that failure is the point:
#    it means an ingest moved a count and nothing has explained it. Here
#    rather than left to a reader for the usual reason -- check-results
#    globs data/results and would report it stale while no step refreshed
#    it.
make landscape-count-move

#    the shared-PMID confound, for the same reason. The 88.5% figure is why
#    doc2/08-results.md calls its own M1 result partly self-scoring, and until
#    2026-11-02 it existed nowhere but that script's console output while
#    being quoted in six places here and one in kg_audit_report.txt.
#    Nothing could check any of them. It is a stamped results file now, so
#    it is audited below and tests/test_label_overlap_figures.py holds both
#    documents to it.
make label-overlap

#    the structure-derived scorer, registered in advance by
#    doc2/registrations/prereg_structure_scorer.md and run for the first time
#    2026-11-11. It reads the release directly and depends on nothing
#    above. It belongs here for the reason the pathway landscape does:
#    check-results globs data/results and would report it stale while no
#    step regenerated it, and tests/test_reproduce_block.py fails on a
#    results file no step here refreshes.
make structure-scorer

#    and the symmetric-control follow-up, registered by
#    doc2/registrations/prereg_structure_scorer_symmetric.md. It writes its OWN file:
#    the run above is its comparator and the verdict beneath that
#    registration quotes both, so one target overwriting the other would
#    destroy the comparison. Order between the two does not matter; each
#    reads the release and neither reads the other's output, except that
#    this one prints a 'prior' column from it when it is present.
make structure-scorer-symmetric

#    and the strength-matched follow-up. ORDER MATTERS HERE, unlike the
#    two above: it derives every k_cell from
#    data/results/structure_scorer.json's counters and refuses outright if
#    that file is absent, so the publication-control run must precede it.
make structure-scorer-strength

#    and the read-out, which is what the scorer is for: the top candidates
#    per virus with the known active each one resembles, named. It is not a
#    protocol and has no registration of its own -- it applies the one the
#    three studies above validated. It also prints how many suggestions
#    were dropped as the same compound as something already known, which
#    is fix-list L27.
make repurposing-candidates

#    the DWPC weighting sensitivity: three arms over every channel,
#    answering whether the step weight reads the held-out labels the
#    channels are scored against. The claim section above quotes its
#    middle arm beside M7's published p, so it is an audited artifact
#    rather than a terminal printout. It registers nothing.
make degree-sensitivity

#    the clinical-label cell, registered by
#    doc2/registrations/prereg_clinical_labels.md and run for the first time
#    2026-11-17. It reads the release and nothing above it; it is here
#    because check-results globs data/results and would report it stale
#    while no step regenerated it, and tests/test_reproduce_block.py
#    fails on a results file no step here refreshes. It is the only
#    protocol in this block measured against a label that is not an
#    in-vitro assay.
make clinical-labels

#    DOES IT PREDICT FORWARD. Registered in advance by
#    doc2/registrations/prereg_temporal_structure_scorer.md, signed 2026-11-19. It
#    reads the 2021 split built in step 1 -- pre-2021 actives as the
#    reference set, post-2021 labels as the evaluated population -- so
#    it must follow the splits and depends on none of the scorer
#    protocols above. Everything else in this block is retrospective:
#    known actives separated from known inactives measured at the same
#    time. This is the only cell that asks whether the separation
#    transfers. It REFUSES rather than reporting numbers if the split is
#    stale, if the family is not the registered 10, or if fingerprint
#    coverage falls below 95%, and a refusal exits 1.
make temporal-structure-scorer

#    EVERY SCORER AUC, EVERY VIRUS, SIDE BY SIDE -- each protocol against
#    its own degree baseline, in doc2/scorer_aucs.txt. It computes nothing
#    and reads the three files written just above, so it goes last among
#    them. It is here rather than left to a reader because the comparison
#    that decides whether any of this is about chemistry was being
#    reassembled by hand after every rebuild: three results files with
#    three shapes, and no written-down answer to "what does each protocol
#    say about each virus". `make scorer-table-check` fails if the file on
#    disk stops describing the results files, and it runs in the test
#    suite rather than here: `make check-results` has to stay the last
#    command in this block, because it audits everything above it.
make scorer-table

# 5. every results file, AND the derived release they depend on, must
#    describe the graph on disk
make check-results
```

`hard_negatives.py` will **skip** its temporal column, rather than compute
one, if `data/releases/v0.1-train2021` was cut from a different graph or
carries no stamp — which is why step 1 precedes step 3. The skip and its
reason are recorded in `hard_negatives.json` under `temporal_skipped`, and
the split's own policy under `train_split`.

Diagnostics behind the causes: `diagnose_hop_attrition.py` (hop-by-hop
survivors), `diagnose_label_overlap.py` (population overlap, shared PMIDs),
`diagnose_assay_type.py` (cellular readouts filed against proteins),
`diagnose_temporal_leak.py` (which source put a label on which side of the
cutoff).
