# Changes needed — audit of 2026-11-03

Source: a line-by-line read of `src/`, `scripts/`, `tests/`, `web/`, `config/`, and
measurements taken directly on `data/releases/v0.1` (78,789 nodes, 392,907 edges).
Every item was confirmed against code or data unless marked *(suspected)*.

Diagrams (drawn from the code, defects in red): `docs/figures/kgav_architecture.svg`,
`docs/figures/kgav_metagraph.svg` (PNG versions alongside; `.dot` sources to redraw).

![architecture](figures/kgav_architecture.png)

![metagraph](figures/kgav_metagraph.png)

How to run the agent: `docs/agent_guide.md`.

Order of work: **C → B → re-run → docs**. Estimated 3–4 working days in total;
compute is minutes (a results regeneration is about 11 min).

---

## C. Evidence dates (temporal split is invalid until fixed) — ~2–3 days

| # | Problem | Where | Fix |
|---|---|---|---|
| C1 | Every virus–host PPI edge (12.5k) is dated `2024-01-01`; `--date` help says "fallback when no PMID year", but it is applied to every row | `scripts/ingest_vhppi.py:156`, `src/kgav/vhppi.py:424` | Date each edge by its PMID's publication year (PMID already parsed in `parse_row`); fall back only when there is no PMID |
| C2 | ~~All 141,019 STRING edges are dated `2023-08-28` (STRING v12 release date) → `v0.1-train2021` holds **0** `PHYSICALLY_INTERACTS_WITH` edges, so M3/M5 vanish from the temporal evaluation~~ **MEASURED AND DECLARED 2026-11-09, not repaired.** The stated symptom was wrong in one direction and worse in another: the 2021 graph holds **11,744** `PHYSICALLY_INTERACTS_WITH` edges, and **every one is virus–host** (VirHostNet, IntAct, MINT, DIP). All 141,019 host–host edges are STRING's and **0** survive, because STRING publishes no per-interaction date. So M3 (virus–host middle hop) runs and M5 (host–host) cannot — which is why M5 was absent from every SARS-CoV-2 temporal table ever published. `temporal_split.py` now writes `training_graph_predicates` and `training_graph_predicate_sources`, so this is a recorded fact rather than one that had to be walked by hand. **Not repairable offline:** STRING v11.0 (2019) is a 139 MB download and `data/raw` is hand-populated |
| C3 | ~~All GtoP `TARGETS` edges dated `1970-01-01` (= "undated") → always kept in the training graph~~ **MEASURED AND DECLARED 2026-11-09, not repaired.** Confirmed: 5,306 undated `TARGETS` edges enter the 2021 graph under the default `include` — 3,601 GtoP, 1,564 ChEMBL, 141 both — and `undated_census` in the temporal output now names them by predicate and source. **The proposed fix has no offline source.** GtoP edges carry PMIDs but no years, and a PMID→year index derived from the release is UNSOUND: 42.9% of PMIDs get different years from different dated edges, with spreads up to 41 years, because `first_asserted_date` on a merged edge belongs to no single publication in its list. Dating them by the GtoP release date (2026-10-12) would drop all 3,601 from every pre-2023 training graph, which is a protocol change after seeing the numbers and belongs in a prereg |
| C4 | ~~All 3,384 similarity edges dated 1970 and computed over the *full* compound set → post-cutoff compounds leak into the 2021 graph (M6/M7)~~ **FIXED 2026-11-09, and it was holding up every temporal result.** `kgav.temporal.derived_year` dates an undated COMPUTED edge by the later of its endpoints' first appearances. Measured: of 3,173 computed edges, 2,864 are datable and **1,412 (44% of the layer) were post-cutoff**. Effect: SARS-CoV-2 reachable test compounds **231 (14.7%) → 64 (4.1%)**, evaluable viruses **7 → 3**, COMBINED AUC 0.526 [0.503, 0.549] → 0.509 [0.486, 0.532], degree 0.555 → 0.516, and **all five temporal cells whose interval excluded chance are gone**. `exclude` gained two evaluable viruses, because 1,452 computed edges are now genuinely dated pre-cutoff rather than discarded as dateless |
| C5 | Fold edges dated 1970 with source `infores:interpro` are not treated as "computed" → dropped under `computed_only` | `src/kgav/folds.py:712,766`, `temporal.COMPUTED_SOURCES` | Add the folds source to `COMPUTED_SOURCES`, or date them |
| C6 | ~~5,306 undated `TARGETS` kept under the default `--undated include`~~ **MEASURED 2026-11-09; the policy change is a prereg decision.** The 5,306 are now counted by source in `undated_census`. Making `computed_only` the default is defensible and would be chosen after seeing that `include`'s one chance-clearing result was a leak, which is exactly the ordering this project forbids. See C3 |
| C7 | ~~`temporal_split.py` never calls `metapath_reach`; a channel with zero paths is silently skipped (`if sc:`)~~ **FIXED 2026-11-09.** A channel that scores nothing now records a `not_evaluated` block with its metapath, the training graph's edge count for each of its predicates, the per-source breakdown, and which predicates are missing — so "never ran" is distinguishable from "ran and found nothing". The reason text says only what the census supports: M5's predicates are all present and its host–host subset is not, which the source breakdown shows and a predicate-level count cannot |
| C8 | `hard_negatives.py` temporal column applies no overlap policy (compounds measured on both sides of the cutoff) | `scripts/hard_negatives.py:406-408` | Use `split.test_labels` from `temporal.build_split` |

## B. Label reachability (why every channel sits at AUC ≈ 0.5) — ~1 day

Measured: **13,916 of 19,506 labelled compounds (71%) have no edge besides their label**.
Of those, only **1,409** have any human target at pChEMBL ≥ 6 in ChEMBL 37 (any phase).
The rest are screening compounds with no known target; no data source fixes them.

| # | Problem | Fix |
|---|---|---|
| B1 | `host_targets` keeps only `max_phase ≥ 2` drugs | Drop the phase filter (keep pChEMBL ≥ 6) → recovers ≤1,409 compounds |
| B2 | ~~Unreached compounds are imputed to the bottom of every ranking, so AUC mostly measures coverage~~ **FIXED 2026-11-09.** `kgav.labels.auc_floor` computes what a given coverage scores with no ranking skill (`0.5 + (cov_pos - cov_neg) / 2`, derived from the pair counting), `auc_skill` the fraction of the floor-to-ceiling span a scorer used, and both are emitted beside `auc_ceiling`. Both protocols are now decomposed in `docs/evaluation.md`, pinned by `tests/test_auc_decomposition.py`. **Cross-sectional:** 86% of COMBINED's margin above chance is coverage asymmetry and 0% of degree's is, since degree ranks every compound; normalised, every scorer uses 13.0–26.3% of its available span and none separates from the others. **Temporal:** 53 of 69 AUC cells equal their floor to four decimal places, median deviation exactly zero, largest 0.0067, largest floor-to-ceiling span 0.0320 — so the tie-safe AUC on that protocol **cannot** separate, whatever scorer is applied. On 2026-11-08 the decomposition showed COMBINED's chance-clearing 0.5263 to be exactly its own floor (ceiling 0.5266); the C4 leak fix put it below the bar the next day, so two independent findings retired the same result. An interim note on 2026-11-08 said the temporal coverages were not recorded — wrong, the fields are `auc_coverage_pos`/`auc_coverage_neg` |
| B3 | **Decision needed:** restate the task as "compounds with a known mechanism" (≈30% of labels) | Owner's call before step 2 |

## L. Label quality

| # | Problem | Where | Fix |
|---|---|---|---|
| L1 | 69% of `HAS_ANTIVIRAL_ACTIVITY_AGAINST` edges are measured *inactive* (censored `>`; 7,333 NCATS curve class 4). `labels.py` re-classifies, but the agent and exports read them as active | `src/kgav/chembl.py`, `src/kgav/ncats.py` | Emit inactives under `MEASURED_INACTIVE_AGAINST` (or rename the predicate to "MEASURED_AGAINST") |
| L2 | ~~645 contradictory ChEMBL (compound, virus) pairs are resolved by SQL row order, not dropped (`labels.py:23` claims they are dropped)~~ **NOT REPRODUCIBLE AS WRITTEN, checked 2026-11-09.** `build_labels` does drop them: 150 (compound, virus) pairs on v0.1 have edges that disagree, all of them excluded and counted in `ambiguous_dropped`, and no `ORDER BY`/`GROUP BY` in `chembl.py`'s activity query collapses a contradiction before it gets there. The real defect was one level down and is recorded as L3: the contradiction inside a single MERGED edge, decided by `POTENCY_FIELDS` order. Fixing L3 took `ambiguous_dropped` from 150 to 114 |
| L3 | ~~1,568 merged edges hold several potency fields and one `relation`; 385 change verdict depending on the field read~~ **FIXED 2026-11-09.** Re-measured on v0.1: 1,578 label edges carry two potency fields and **173** disagree on the verdict — `_potency` returned the first field present in `POTENCY_FIELDS`, so the verdict came from the declaration order of a module constant, and reversing the tuple flipped 118 positives to negatives and 55 the other way. `classify` now returns `None` for a contested edge, which is the rule this module already states for ambiguity and that `build_labels` already applies to compounds whose SEPARATE edges disagree. Effect: active 5,643→5,525, inactive 20,293→20,238, undecidable 342→515, and **ambiguous_dropped 150→114** — 36 compounds were being discarded as self-contradictory on the strength of a verdict chosen by tuple order, so the fix recovers them. The canonical 89-positive protocol is unchanged (the NCATS CPE screen carries one potency field per row). A rebuild confirmed the release byte-identical: `chembl.py` classifies single pre-merge rows, so the contradiction is created by `merge_edges` and only the label predicate is read after it — pinned by `test_the_protein_layer_has_no_contested_edges` |
| L4 | NULL ChEMBL relation becomes `=` (bound read as exact) | `chembl.py:564`, `host_targets.py:370` | Count and drop NULL relations |
| L5 | Censored CC50 written as exact with `selectivity_verified: true`; MIN/MIN SI is not conservative; tox CTE unrestricted | `src/kgav/selectivity.py:233-237, 404-462` | Carry the relation; pair max EC50 with min CC50; restrict to the same assay/cell line |
| L6 | 434 of 3,384 similarity edges join two keys of the same molecule (same 14-char InChIKey skeleton) | `src/kgav/similarity.py` | Exclude same-skeleton pairs |
| L7 | Degree baseline counts held-out label edges | `src/kgav/baselines.py:178` | Count degree after the held-out skip |

## E. Statistics and reporting

| # | Problem | Fix |
|---|---|---|
| E1 | ~~No multiple-testing correction across ~150 tests in `compare_metapaths.py` / `temporal_split.py`~~ **FIXED 2026-11-09.** `kgav.baselines.correction_families` derives four candidate families for any run and both scripts emit a `correction` block; `tests/test_correction_families.py` asserts `docs/evaluation.md`'s sensitivity table against it, verdict column included. The four sizes had been typed by hand and nothing could recompute one — three were right, and the fourth, `cells_computed` at 61, had a definition recorded nowhere (cells the run PRODUCED: not cells with a p-value, 34, nor cells planned, 120). Cross-sectional: 10 / 61 / 120 / 150, and only SARS-CoV-2 M7 survives the smallest. Temporal: 2 / 17 / 24 / 30, and no readable cell survives any — while the four significant rows that DO clear the strictest family all fail `trustworthy`, so multiplicity is the weaker of the two gates. **The family is reported, not chosen**: fixing one now, with the p-values on screen, is the ordering this project forbids, so the decision stays in the prereg |
| E2 | `baselines_*.json` Hits@k/MRR depend on tie order; no range flag | Use `hits_at_k_range` there too |
| E3 | ~~`vs_degree_readable` is True exactly when degree is significant (inverted)~~ **NOT A DEFECT, checked 2026-11-07.** It is true exactly then and that is correct: the flag asks whether the RATIO can be read, and a quotient over a denominator at chance is large because the divisor is small (SARS-CoV-2 M4 reads 2.04x a degree lift of 0.75, p=0.8167, while M4's own lift is 1.53 at p=0.2630). Whether the CHANNEL's result is trustworthy is a different question and is `confound_warning`'s BASELINE AT CHANCE verdict. The two share the words "degree at chance" and one comment in `baselines.py` answered both at once, which is what this item read as an inversion; that comment is now split and three tests pin the semantics. | `baselines.py:924-925` |
| E4 | `pathway_landscape.json` predates the per-row seeding fix; `check_results` still calls it current | Re-run; make `check_results` compare the code commit too |
| E5 | `network_proximity` averages observed vs null over different denominators *(suspected)* | Penalise unreachable targets symmetrically |

## D. Docs that contradict the project's own results

| # | Claim | Truth (from `data/results`) |
|---|---|---|
| D1 | ~~`docs/evaluation.md:55` "no scorer excludes chance… degree only cross-sectionally"~~ **FIXED 2026-11-08, and the sentence it replaced became true on 2026-11-09.** The claim block prints the derived census of every cell whose AUC interval excludes chance, asserted both ways by `tests/test_evaluation_table.py::test_the_claim_block_enumerates_every_cell_that_clears_chance`. It was seven cells (five temporal) when derived, which is what made the original wording false; fixing the similarity leak removed all five, so it is now two cells, both cross-sectional, both SARS-CoV-2 — the original words, true for the first time. A sentence that is right today and was wrong yesterday with the same words is the whole argument for deriving it |
| D2 | ~~`evaluation.md:59` "4 of 4, no exception"~~ **FIXED 2026-11-08.** Re-measured: five cross-sectional viruses have a channel at p<0.05 — SARS-CoV-2 (M7), DENV-2 (M1, M4, M5), Zika (M5), West Nile (M4), MERS-CoV (M4, M5) — and degree clears it in one, MERS-CoV. The sentence had counted the viruses where DEGREE was significant. Corrected in `evaluation.md` and in the amendment to `docs/prereg_temporal_auc.md`, which cited it |
| D3 | ~~`evaluation.md:49-69` temporal AUCs~~ **FIXED 2026-11-08.** Every cell of the temporal table was stale by two rebuilds; it is re-derived and now pinned by `tests/test_temporal_table.py`. Two further defects fell out of doing it: the pre-registered AUC bar of 0.586 was computed on the reach count (48) instead of the AUC's own n (`auc_n_pos` 1,580 / `auc_n_neg` 982, where the bar is 0.523), and "degree is at chance on the temporal split" was false — degree read 0.555 [0.533, 0.578] there. Both are recorded in the prereg amendment. **The second one became true again on 2026-11-09**: with the similarity leak fixed degree reads 0.516 [0.493, 0.539], containing chance, so the original claim was right and both the figure and the mechanism offered for it were wrong |
| D4 | `docs/v1_release_plan.md:16-40` node/edge counts and overlap ratios | Stale; see `pathway_landscape.json` |
| D5 | `scripts/diagnose_doc_figures.py` matches figures by value, not by virus/channel → false passes | Match on (virus, channel, metric) |

**State the real headline:** ~~degree beats every reasoning channel on every
scored virus.~~ **Re-measured 2026-11-08 and withdrawn.** Degree's apparent
margin was substantially the baseline counting the held-out label edges
themselves (L7). With them excluded, cross-sectional SARS-CoV-2 reads degree
0.577 [0.514, 0.639] against COMBINED 0.575 [0.512, 0.637] — 0.002 apart and
not distinguishable. The real headline is that an eight-channel ensemble and
counting a compound's edges perform identically, and both barely clear
chance.

## A. Agent (`src/kgav/agent`, `scripts/kgav_ask.py`)

| # | Problem | Status |
|---|---|---|
| A1 | Evidence question returned "no edge" because taxa are listed first and edges run compound → virus | **Fixed 2026-11-03** (`orchestrator.py`, tries both orders; agent tests pass) |
| A2 | Reach caveat never loads: `--results` is a file, `labels.measured_reach` expects a directory | `kgav_ask.py:147` → pass `args.results.parent` |
| A3 | Triage ranks by selectivity taken from held-out label edges, without the ANSWER SHEET caveat | `tools.py:220-231, 766-769, 859-891` |
| A4 | "Best direct IC50" mixes Ki/Kd and is not virus-scoped | `tools.py:745-752, 895` |
| A5 | Caveat check is a substring match (e.g. "evidence" satisfies EVIDENCE) | `verify.py:266-271` |
| A6 | `--json` prints answers that failed verification | `kgav_ask.py:191-192` |
| A7 | Refusal says "every scorer performs at chance" without a measurement (`temporal` never passed) | `orchestrator.py:366-368` |
| A8 | Ollama timeout / bad JSON is uncaught → crash *(suspected)* | `llm.py:80-86` |
| A9 | 7B model's prose often fails verification → answer withheld | Better prompt, or a bigger model |

## P. Pipeline, ingest and process

| # | Problem | Fix |
|---|---|---|
| P1 | CI runs full `ruff check` → 281 findings → always red (Makefile says 29) | Pin ruff; run `F821,F811` + fix or select a ruleset |
| P2 | `make gate` / CI treat **any** non-zero exit (even an ImportError) as "validator rejected planted errors" | Assert the expected violation codes |
| P3 | Only `vhppi` resolves secondary UniProt accessions; `normalize_bio.UniProtResolver` is unused | Use the alias index in chembl, host_targets, gtop |
| P4 | ORCS stores a calibrated `direction` beside an uncalibrated `effect_size` (signs disagree) | Multiply `effect_size` by the same sign |
| P5 | Shared Disease nodes get the last virus's label ("DENV-4 disease") | Use the MONDO label |
| P6 | M8 gate admits unknown-vs-unknown catalytic types as a match | Refuse `unknown` |
| P7 | `ncats` assumes molar `conc_*` columns | Infer/verify the unit |
| P8 | `host_targets` dedups on (drug, target, type), so the kept potency depends on row order | Include value/assay in the key |
| P9 | Stale results dirs `data/results-corrected`, `data/results-prechange` | Delete or move out of the repo |
| P10 | 22 one-off `diagnose_*` scripts, about 190 tests that only check doc phrases | Archive the diagnostics; replace phrase tests with generated tables |
| P11 | `web/` is an unused Manus template (dead components, dev endpoint writing logs, Forge key proxy) | Remove until a frontend is needed |
| P12 | `.gitignore` has `!data/releases/**/MANIFEST.json`, but `data/releases/*` ignores the directories first, so git can never re-include them: no release manifest (layer hashes, counts) is tracked | Ignore `data/releases/*/*` and un-ignore `!data/releases/*/` + `!data/releases/*/MANIFEST.json`, then add them |
