# 10 — Reproducing everything

```
make install       # pip install -e ".[dev]"
make test          # the gate: lint subset + the suite, slow tests deselected
make lint          # the full declared ruff set
make gate          # the validator's NEGATIVE CONTROL
make check-results # is every artifact current for the release on disk
make reproduce     # regenerate every results file, in order
```

## The order is read from the document, not transcribed

`make reproduce` runs the Reproducing block in `docs/evaluation.md` **as the
document holds it**. The order is load-bearing and was being transcribed by
hand after each rebuild, where a dropped flag changes a protocol without
failing. `make reproduce-print` prints the steps without running them, and
`tests/test_reproduce_block.py` fails if an audited results file has no step
that regenerates it — because anything under `data/results` is globbed by
`check-results` and would be reported stale with nothing able to refresh it.

## Every target, and the command it runs

60 targets. The **recipe** is shown rather than a
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
| `make arch-figures` | `dot -Tpng docs/figures/arch.dot    -o docs/figures/kgav_architecture.png`<br>`dot -Tsvg docs/figures/arch.dot    -o docs/figures/kgav_architecture.svg`<br>`dot -Tpdf -Gsize=8.07,11.49            docs/figures/arch_a4.dot \`<br>`-o docs/figures/kgav_architecture_a4.pdf`<br>`dot -Tpng -Gsize=8.07,11.49 -Gdpi=300  docs/figures/arch_a4.dot \`<br>`-o docs/figures/kgav_architecture_a4.png` |
| `make results-doc` | `python scripts/write_results_doc.py` |
| `make results-doc-check` | `python scripts/write_results_doc.py --check` |
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
| `make rewrite-tables` | `python scripts/rewrite_derived_tables.py $(ARGS)` |
| `make rewrite-tables-check` | `python scripts/rewrite_derived_tables.py --check` |
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
| `clinical_labels.json` | 3,414 | result | 873623c9 | 2026-11-17 |
| `degree_weighting_sensitivity.json` | 102,617 | result | 873623c9 | 2026-11-17 |
| `hard_negatives.json` | 13,676 | result | 873623c9 | 2026-11-17 |
| `hard_negatives_row1.json` | 67,443 | result | 873623c9 | 2026-11-17 |
| `hard_negatives_row2.json` | 70,334 | result | 873623c9 | 2026-11-17 |
| `hard_negatives_row3.json` | 52,538 | result | 873623c9 | 2026-11-17 |
| `hard_negatives_row4.json` | 52,204 | result | 873623c9 | 2026-11-17 |
| `label_overlap.json` | 699 | result | 873623c9 | 2026-11-17 |
| `landscape_count_move.json` | 9,829 | result | 873623c9 | 2026-11-18 |
| `metapath_comparison.json` | 68,897 | result | 873623c9 | 2026-11-17 |
| `pathway_landscape.json` | 64,350 | result | 873623c9 | 2026-11-17 |
| `repurposing_candidates.json` | 394,919 | result | 873623c9 | 2026-11-21 |
| `structure_scorer.json` | 83,978 | result | 873623c9 | 2026-11-17 |
| `structure_scorer_strength_matched.json` | 63,229 | result | 873623c9 | 2026-11-17 |
| `structure_scorer_symmetric.json` | 131,956 | result | 873623c9 | 2026-11-17 |
| `temporal_2021_computed_only.json` | 17,319 | result | 873623c9 | 2026-11-17 |
| `temporal_2021_exclude.json` | 17,221 | result | 873623c9 | 2026-11-17 |
| `temporal_2021_include.json` | 39,232 | result | 873623c9 | 2026-11-17 |
| `temporal_2021_include_overlap-keep.json` | 2,798 | REFUSED | 873623c9 | 2026-11-17 |
| `temporal_structure_scorer.json` | 36,524 | result | 873623c9 | 2026-11-19 |

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
| `make test` | ~2 min | 3,061 tests, slow deselected, `-n auto --dist loadfile` |
| `make lint` | seconds | ruff, pinned in the dev extra |
| `make check-results` | ~10 s | hashes every release file once |
| `make results-doc-check` | ~25 s | re-runs the hit-rate script so the qualification is measured, not copied |
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
| `make results-doc-check` | `docs/results.md` drifting from the artifacts |
| `make rewrite-tables-check` | any marked table in a document drifting from the test that derives it |
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
