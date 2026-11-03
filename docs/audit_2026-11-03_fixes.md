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
| C2 | All 141,019 STRING edges are dated `2023-08-28` (STRING v12 release date) → `v0.1-train2021` holds **0** `PHYSICALLY_INTERACTS_WITH` edges, so M3/M5 vanish from the temporal evaluation | `src/kgav/host.py` | Build the ≤2021 training graph from STRING v11.0 (2019), or exclude STRING from the temporal protocol and say so |
| C3 | All GtoP `TARGETS` edges dated `1970-01-01` (= "undated") → always kept in the training graph | `src/kgav/gtop.py:375,486` | Use the interaction's reference year |
| C4 | All 3,384 similarity edges dated 1970 and computed over the *full* compound set → post-cutoff compounds leak into the 2021 graph (M6/M7) | `src/kgav/similarity.py`, `src/kgav/temporal.py:52-87` | Recompute similarity per split using only pre-cutoff compounds |
| C5 | Fold edges dated 1970 with source `infores:interpro` are not treated as "computed" → dropped under `computed_only` | `src/kgav/folds.py:712,766`, `temporal.COMPUTED_SOURCES` | Add the folds source to `COMPUTED_SOURCES`, or date them |
| C6 | 5,306 undated `TARGETS` kept under the default `--undated include` | `src/kgav/temporal.py` | Make `computed_only` the default once C3 is fixed |
| C7 | `temporal_split.py` never calls `metapath_reach`; a channel with zero paths is silently skipped (`if sc:`) | `scripts/temporal_split.py:193-196` | Report reach for every channel, including zero |
| C8 | `hard_negatives.py` temporal column applies no overlap policy (compounds measured on both sides of the cutoff) | `scripts/hard_negatives.py:406-408` | Use `split.test_labels` from `temporal.build_split` |

## B. Label reachability (why every channel sits at AUC ≈ 0.5) — ~1 day

Measured: **13,916 of 19,506 labelled compounds (71%) have no edge besides their label**.
Of those, only **1,409** have any human target at pChEMBL ≥ 6 in ChEMBL 37 (any phase).
The rest are screening compounds with no known target; no data source fixes them.

| # | Problem | Fix |
|---|---|---|
| B1 | `host_targets` keeps only `max_phase ≥ 2` drugs | Drop the phase filter (keep pChEMBL ≥ 6) → recovers ≤1,409 compounds |
| B2 | Unreached compounds are imputed to the bottom of every ranking, so AUC mostly measures coverage | Report AUC on the reachable subset **plus** coverage, as the primary number |
| B3 | **Decision needed:** restate the task as "compounds with a known mechanism" (≈30% of labels) | Owner's call before step 2 |

## L. Label quality

| # | Problem | Where | Fix |
|---|---|---|---|
| L1 | 69% of `HAS_ANTIVIRAL_ACTIVITY_AGAINST` edges are measured *inactive* (censored `>`; 7,333 NCATS curve class 4). `labels.py` re-classifies, but the agent and exports read them as active | `src/kgav/chembl.py`, `src/kgav/ncats.py` | Emit inactives under `MEASURED_INACTIVE_AGAINST` (or rename the predicate to "MEASURED_AGAINST") |
| L2 | 645 contradictory ChEMBL (compound, virus) pairs are resolved by SQL row order, not dropped (`labels.py:23` claims they are dropped) | `src/kgav/assemble.py:84-109, 247` | Keep each measurement as its own edge (assay id in the key) or mark conflicts and drop in `build_labels` |
| L3 | 1,568 merged edges hold several potency fields and one `relation`; 385 change verdict depending on the field read | `src/kgav/labels.py:69-95` | Fixed by L2 (one measurement per edge) |
| L4 | NULL ChEMBL relation becomes `=` (bound read as exact) | `chembl.py:564`, `host_targets.py:370` | Count and drop NULL relations |
| L5 | Censored CC50 written as exact with `selectivity_verified: true`; MIN/MIN SI is not conservative; tox CTE unrestricted | `src/kgav/selectivity.py:233-237, 404-462` | Carry the relation; pair max EC50 with min CC50; restrict to the same assay/cell line |
| L6 | 434 of 3,384 similarity edges join two keys of the same molecule (same 14-char InChIKey skeleton) | `src/kgav/similarity.py` | Exclude same-skeleton pairs |
| L7 | Degree baseline counts held-out label edges | `src/kgav/baselines.py:178` | Count degree after the held-out skip |

## E. Statistics and reporting

| # | Problem | Fix |
|---|---|---|
| E1 | No multiple-testing correction across ~150 tests in `compare_metapaths.py` / `temporal_split.py` | Bonferroni or BH, as `landscape.py` already does |
| E2 | `baselines_*.json` Hits@k/MRR depend on tie order; no range flag | Use `hits_at_k_range` there too |
| E3 | `vs_degree_readable` is True exactly when degree is significant (inverted) | `baselines.py:924-925` |
| E4 | `pathway_landscape.json` predates the per-row seeding fix; `check_results` still calls it current | Re-run; make `check_results` compare the code commit too |
| E5 | `network_proximity` averages observed vs null over different denominators *(suspected)* | Penalise unreachable targets symmetrically |

## D. Docs that contradict the project's own results

| # | Claim | Truth (from `data/results`) |
|---|---|---|
| D1 | `docs/evaluation.md:55` "no scorer excludes chance… degree only cross-sectionally" | SARS-CoV-2 degree 0.624 [0.561, 0.686] and COMBINED 0.581 [0.518, 0.643] exclude 0.5; degree also excludes chance under the temporal split |
| D2 | `evaluation.md:59` "4 of 4, no exception" | 5 viruses with p<0.05; degree cleared it in 1 |
| D3 | `evaluation.md:49-69` temporal AUCs | Stale (pre-GtoP release) |
| D4 | `docs/v1_release_plan.md:16-40` node/edge counts and overlap ratios | Stale; see `pathway_landscape.json` |
| D5 | `scripts/diagnose_doc_figures.py` matches figures by value, not by virus/channel → false passes | Match on (virus, channel, metric) |

**State the real headline:** degree beats every reasoning channel on every scored virus.

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
