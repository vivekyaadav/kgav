# 15 — Appendix: every script, and every artifact

## Scripts

86 scripts, 31,730 lines. Four of them exit non-zero on
`--help` **on purpose** — two print a usage block because they take positional
directories, one needs a release, and one **refuses to run at all** because it
is a phase-1 schema patch and would downgrade the schema from
0.16.0. The test that covers the sweep greps for a traceback
rather than checking exit codes, for exactly that reason.

| script | lines | what it does |
|:---|---:|:---|
| `agent_server.py` | 358 | kgav agent server — the chat frontend's backend. |
| `assemble_release.py` | 309 | Day 8 driver: assemble layers, validate the whole graph, write a manifest. |
| `auc_contribution.py` | 576 | Where the AUC actually comes from, as a share per approach. |
| `audit_build.py` | 512 | What has been ingested, what has been built, and what is missing. |
| `calibrate.py` | 79 | What each reasoning channel is measured to be worth, and whether it earns |
| `check_results.py` | 181 | Flag results files that no longer describe the release on disk. |
| `clinical_labels.py` | 274 | The registered clinical-label analysis: doc2/registrations/prereg_clinical_labels.md. |
| `compare_metapaths.py` | 209 | Compare metapaths across all seven viruses. |
| `compare_results.py` | 163 | Two results directories: did a rebuild change any NUMBER, or only the stamp? |
| `compare_results_aucs.py` | 356 | Every AUC in data/results, this release against an earlier commit. |
| `data_audit.py` | 325 | What the schema declares against what the release actually holds. |
| `data_timeline.py` | 327 | What data this graph is built from, when it arrived, and what it covers. |
| `diagnose_assay_type.py` | 180 | How many INHIBITS edges claim `biochemical` on a CELLULAR measurement? |
| `diagnose_comparator_cost.py` | 188 | What building each missing comparator virus would buy M8, and what it costs. |
| `diagnose_controls.py` | 309 | Why each biology control ranks where it does. |
| `diagnose_doc_figures.py` | 196 | Which figures in doc2/08-results.md still appear in the results files. |
| `diagnose_fingerprint_blindspots.py` | 162 | What the ECFP4 fingerprint cannot tell apart, and what that is worth. |
| `diagnose_flavivirus_join.py` | 540 | Would a flavivirus build produce a graph with labels, or an empty shell? |
| `diagnose_flavivirus_targets.py` | 416 | Which flavivirus proteins do ChEMBL's assay descriptions actually name? |
| `diagnose_frontier_size.py` | 102 | How large does the frontier actually get, per metapath? |
| `diagnose_hop_attrition.py` | 152 | Where do the host-directed metapaths lose their candidates? |
| `diagnose_identity_fixed_points.py` | 199 | Is every node id what the normaliser derives from that node's own SMILES? |
| `diagnose_label_degree_weighting.py` | 297 | Does DWPC weighting read the answer sheet, and by how much. |
| `diagnose_label_overlap.py` | 236 | Two questions the hop-attrition run raised but could not answer. |
| `diagnose_label_reachability.py` | 398 | WHAT RESTATING THE TASK WOULD COST (fix-list B3), measured per virus. |
| `diagnose_landscape_count_move.py` | 331 | Which cells of the pathway landscape moved between two releases, and why. |
| `diagnose_m1_potency.py` | 309 | Is there enough measured potency to break M1's tie, or is there not? |
| `diagnose_orphans.py` | 239 | Name every VIRAL orphan in a release, and say which failure it is. |
| `diagnose_phase_gate.py` | 247 | What host_targets' max_phase gate costs, as three options rather than one. |
| `diagnose_property_conflicts.py` | 256 | Which layers disagree about a node property, and whether it is a real |
| `diagnose_prospective_analogy.py` | 200 | Is the prospective result prediction, or retrieval of a known series? |
| `diagnose_protocol_disagreement.py` | 260 | Why one cell can anti-rank cross-sectionally and clear prospectively. |
| `diagnose_reach_ceiling.py` | 202 | Is a channel's reach limited by missing data or by a retrieval parameter. |
| `diagnose_scorer.py` | 285 | How much of each channel's score is the compound's own degree. |
| `diagnose_serotype_attribution.py` | 331 | Does the assay description name the virus the row is filed under? |
| `diagnose_source_scope.py` | 243 | What kind of literature is ChEMBL, and is it the kind this project needs? |
| `diagnose_structure_scorer_design.py` | 279 | The design facts a structure-derived scorer has to be registered against. |
| `diagnose_tautomer_conflation.py` | 314 | Does ONE layer file two tautomers under one InChIKey? (fix-list L19) |
| `diagnose_temporal_leak.py` | 180 | Which compounds sit on both sides of the temporal cutoff, and what put them there. |
| `diagnose_unresolved_ids.py` | 329 | What the VirHostNet rows counted as `unresolved_id` actually are. |
| `diff_releases.py` | 117 | Are two releases the same GRAPH, or the same graph in a different ORDER? |
| `dryrun_promote.py` | 649 | What would change if these surveyed viruses were promoted to `viruses:`. |
| `fetch_raw.py` | 477 | Fetch the UniProt files ingest_spine reads, and record where they came from. |
| `fetch_vhppi.py` | 368 | Fetch virus-host interactions per virus from a PSICQUIC service. |
| `filter_intact.py` | 283 | Cut the viral rows out of IntAct's full release, without unpacking it. |
| `hard_negatives.py` | 651 | Evaluate against MEASURED negatives, cross-sectionally and temporally. |
| `hit_rate_at_k.py` | 447 | The wet-lab hit rate: of the top-ranked compounds, how many were active. |
| `ingest_chembl.py` | 235 | Day 7 driver: chemistry from ChEMBL 37. |
| `ingest_curated.py` | 198 | Driver: host factors from published screens ORCS does not carry. |
| `ingest_folds.py` | 1020 | Phase 2 driver — build the fold layer across three virus families. |
| `ingest_gtop.py` | 282 | Driver: drug -> human protein TARGETS from the Guide to PHARMACOLOGY. |
| `ingest_host.py` | 154 | Day 4 driver: build the host layer release. |
| `ingest_host_targets.py` | 234 | Driver: drug -> human protein TARGETS edges. Unblocks M2 through M6. |
| `ingest_indications.py` | 195 | Driver: drug -> Disease TREATS / IN_TRIAL_FOR from ChEMBL indications. |
| `ingest_ncats.py` | 409 | Ingest the NCATS OpenData SARS-CoV-2 CPE screen and its counterscreen. |
| `ingest_orcs.py` | 278 | Day 6 driver: host dependency/restriction factors from BioGRID ORCS. |
| `ingest_selectivity.py` | 230 | Build the selectivity layer from same-document EC50/CC50 pairs. |
| `ingest_similarity.py` | 183 | Build CHEMICALLY_SIMILAR_TO edges so M6 (cold start) has something to walk. |
| `ingest_spine.py` | 615 | Day 3: build the taxonomy and viral proteome spine. |
| `ingest_vhppi.py` | 364 | Day 5 driver: virus-host PPI bridge layer. |
| `kgav_ask.py` | 431 | kgav ask — query the graph from the command line. |
| `label_headroom.py` | 376 | How many more labels would it take, and are there compounds to label. |
| `partial_degree.py` | 416 | Does any channel carry signal beyond degree. |
| `pathway_landscape.py` | 999 | The druggable host-factor landscape, without a scorer. |
| `phase0_verify_folds.py` | 361 | Phase 0 — verify fold classifications AT CHAIN LEVEL. |
| `phase1_schema_patch.py` | 242 | Phase 1 — schema 0.8.0: TargetClass, fold edges, and metapath M8. |
| `phase4_hcv_retrospective.py` | 270 | Phase 4: the broad set, measured exactly as doc2/registrations/prereg_phase4_hcv.md |
| `phase4_survey_chemistry.py` | 1096 | Phase 4, step 0 — can the HCV retrospective be measured at all? |
| `rebuild_release.py` | 334 | Run the whole ingest chain in dependency order and stop on a real failure. |
| `record_raw_provenance.py` | 219 | Stamp a raw input that was placed by hand, and re-verify it later. |
| `repurposing_candidates.py` | 1509 | What the graph suggests, by name, for each virus. |
| `run_baselines.py` | 466 | Run the three non-learned baselines and the biology controls. |
| `run_reproduce_block.py` | 136 | Run the Reproducing block in doc2/12-reproducing.md, from the document. |
| `screen_power.py` | 289 | How many compounds a prospective wet-lab screen would need, per cell. |
| `snapshot_release_operands.py` | 174 | What this release IS, in the figures a later rebuild will move. |
| `structure_scorer.py` | 887 | The structure-derived scorer, exactly as registered. |
| `temporal_split.py` | 375 | Build the temporal split, run the leakage audits, and evaluate on it. |
| `temporal_structure_scorer.py` | 321 | Does the structure scorer predict forward? Pre-2021 references, post-2021 labels. |
| `validate_release.py` | 131 | Validate a KG release against the frozen schema. |
| `verify_mondo.py` | 152 | Check every MONDO term the register declares against the ontology. |
| `verify_polymerase_folds.py` | 220 | Verify the polymerase result, which contradicted the plan. |
| `write_architecture_figure.py` | 576 | Generate the architecture figures from the code that defines the pipeline. |
| `write_doc2.py` | 2180 | Generate doc2/: the whole project described from the code and the artifacts. |
| `write_kg_report.py` | 280 | Generate doc2/kg_release_census.txt from the assembled release. |
| `write_metagraph_figure.py` | 583 | Generate doc2/figures/meta.dot from the assembled release. |
| `write_scorer_table.py` | 269 | Generate doc2/scorer_aucs.txt: every scorer AUC, every virus, side by side. |

## Artifacts

36 files in `data/results/`, every one stamped with the
content hashes of the graph it was computed from and the commit that computed
it. 1 is a recorded refusal.

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

## The constants, and where each one is decided

| constant | value | decided by |
|:---|---:|:---|
| `inactive_above_nm` | 10,000 nM | **measured** from the censored-bound census |
| `active_below_nm` | 10,000 nM | a decision about the task, not derived |
| `dose_series_policy` | `include` | fix-list L10, recorded with both branches runnable |
| ECFP4 | radius 2, 2,048 bits | `kgav.similarity`, one definition |
| `CHANCE` | 0.5 | `kgav.calibration` — imported by the agent, never respelled |
| `MIN_REACHED_POSITIVES` | 10 | matches the minimum class size the hard-negative script requires |
| landscape `ALPHA` | 0.05 | `kgav.landscape`, imported by the agent so one answer cannot carry two alphas |
| `hub_degree_cut` | 460 | the `Protein` p99.9 on this release |
| `THIN_ROUTE` | 25 | not a cutoff: the route is used either way and the number is stated |
| `SI_THRESHOLD` | 10.0 | the conventional selectivity floor, applied as a report rather than a filter |

## Correction families, and the labelled populations

| run | fingerprinted | family | corrected alpha |
|:---|---:|---:|---:|
| publication-disjoint | 19,168 | 14 | 0.003571 |
| symmetric k=10 | 19,168 | 14 | 0.003571 |
| strength-matched | 19,168 | 14 | 0.003571 |
| prospective | — | 10 | 0.005 |
| clinical labels | — | 1 | 0.05 |
| pathway landscape | — | 44 | 0.05 (uncorrected; Bonferroni applied per table) |

## The floors that make an imputed AUC readable

For the plate-matched cell, the degree baseline reaches
100.0% of the positives so its floor is 0.5000 —
exactly chance — and its skill is 14.0%. The ensemble reaches
42.7%, so its floor is 0.5578 and its skill is
13.3%. Those two floors are the whole difference between the
two AUCs, which sit 0.5666 and 0.5700.

## Every registration

12 documents in `registrations/`, each signed before the
result it governs existed, each amended by appending. The subject is read out
of each file's own title line, so a registration added to the directory
appears here without anyone editing this table — the predecessor of this
index was hand-written, and a registration added without an edit to it was
simply absent.

| registration | subject |
|:---|:---|
| [`prereg_clinical_labels.md`](registrations/prereg_clinical_labels.md) | a label that is not an in-vitro assay |
| [`prereg_denv2_prospective_screen.md`](registrations/prereg_denv2_prospective_screen.md) | a prospective wet-lab screen, DENV-2 |
| [`prereg_m1_potency.md`](registrations/prereg_m1_potency.md) | letting potency into the direct-acting channel |
| [`prereg_ncats_entry_assays.md`](registrations/prereg_ncats_entry_assays.md) | the NCATS pseudotyped-entry assays as plate-matched labels |
| [`prereg_partial_degree.md`](registrations/prereg_partial_degree.md) | does any channel carry signal beyond degree |
| [`prereg_pathway_overlap.md`](registrations/prereg_pathway_overlap.md) | the pathway and interactome overlap, re-run |
| [`prereg_phase4_hcv.md`](registrations/prereg_phase4_hcv.md) | the HCV retrospective, as the evidence permits it |
| [`prereg_structure_scorer.md`](registrations/prereg_structure_scorer.md) | a structure-derived scorer, and the leak it has to survive |
| [`prereg_structure_scorer_strength_matched.md`](registrations/prereg_structure_scorer_strength_matched.md) | the same strength, symmetrically |
| [`prereg_structure_scorer_symmetric.md`](registrations/prereg_structure_scorer_symmetric.md) | a control that handicaps both classes equally |
| [`prereg_temporal_auc.md`](registrations/prereg_temporal_auc.md) | a tie-safe metric for the temporal protocol |
| [`prereg_temporal_structure_scorer.md`](registrations/prereg_temporal_structure_scorer.md) | does the structure scorer predict forward, or only separate what is already known? |
