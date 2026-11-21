# 04 — The build pipeline

`src/kgav/pipeline.py` holds the order **as data**, and
`scripts/rebuild_release.py` executes it. A hand-typed sequence gets it wrong
silently: a layer built out of order still builds, from stale inputs, and
nothing says so.

```
make rebuild-plan     # what a rebuild would run, with each step's verdict
make rebuild          # run it
```

## 13 declared layers, 12 in this release

| layer | script | phase | status | reads | raw inputs | source version |
|:---|:---|:---|:---|:---|:---|:---|
| spine | `scripts/ingest_spine.py` | 1 | required | raw | `data/raw/proteomes/*.json.gz` | — |
| host | `scripts/ingest_host.py` | 1 | required | raw | `data/raw/host/*` | string `v12.0` |
| vhppi | `scripts/ingest_vhppi.py` | 1 | optional | raw | `data/raw/vhppi/vhn_*.tab27 or data/raw/vhppi/intact_*.tab27` | — |
| orcs | `scripts/ingest_orcs.py` | 1 | optional | raw | `data/raw/orcs/screens_all.json` | — |
| curated | `scripts/ingest_curated.py` | 1 | optional | raw | `data/raw/curated_screens/*.tsv` | — |
| chembl | `scripts/ingest_chembl.py` | 1 | optional | raw | `data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db` | chembl `chembl_37` |
| hosttargets | `scripts/ingest_host_targets.py` | 1 | optional | raw | `data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db` | chembl `chembl_37` |
| ncats | `scripts/ingest_ncats.py` | 1 | optional | raw | `data/raw/ncats/cpe_assay.csv`, `data/raw/ncats/cpe_tox.csv`, `data/raw/ncats/pp_sars_cov_2.csv`, `data/raw/ncats/pp_sars_cov_2_tox.csv`, `data/raw/ncats/pp_sars_cov.csv`, `data/raw/ncats/pp_sars_cov_tox.csv`, `data/raw/ncats/pp_mers.csv`, `data/raw/ncats/pp_mers_tox.csv` | — |
| folds | `scripts/ingest_folds.py` | 1 | optional | raw | `data/raw/folds/phase2_chains.json` | — |
| selectivity | `scripts/ingest_selectivity.py` | 2 | optional | peers | `data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db` | chembl `chembl_37` |
| gtop | `scripts/ingest_gtop.py` | 2 | optional | peers | `data/raw/gtop/interactions.csv`, `data/raw/gtop/ligands.csv or data/raw/gtop/ligand_id_mapping.csv` | guidetopharmacology `2026.3.1` |
| similarity | `scripts/ingest_similarity.py` | 2 | optional | peers | — | rdkit `2026.03.6` |
| indications | `scripts/ingest_indications.py` | 2 | optional | peers | `data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db`, `data/raw/mondo/mondo.json` | chembl_indications: no release string (declares db, evidence_tier, ref_type, refs_table, table) |

**Absent from this release: `curated` — nothing matches data/raw/curated_screens/*.tsv.** That
is the documented behaviour rather than a failure: a layer whose raw inputs
are missing is skipped and named. It is also the reason the pathway landscape
has no cell for five of the flaviviruses, so the skip is visible in the
results rather than only in the plan.

**Phase 1** reads raw inputs only. **Phase 2** reads the layer *directories*
built before it. The assembly runs **once, last**.

### Why one assembly is enough now, and was not before

Until 2026-11-11 three layers read the **assembled release**, so an assemble
had to run in the middle — and *that* middle assemble consumed the previous
run's output of those same three layers. Two defects followed, both measured
rather than argued:

- **Non-convergence.** `ingest_similarity` partitions compounds by their
  `TARGETS` / `INHIBITS` / `MEASURED_INACTIVE_AGAINST` subjects and `gtop`
  contributes `TARGETS`, so the two were mutually dependent through the
  assembly between them. One revert walked `edges_sha` through three values
  before settling, and the output said "the release CHANGED this run" every
  time.
- **A hard stop.** Renaming a phase-2 layer's qualifier made the middle
  assemble validate 4,264 edges carrying the now-undeclared name and exit 1 —
  *before* the ingest that would have replaced them. Not avoidable by running
  twice: every pass began by assembling the same stale directory.

The three layers now read `peer_nodes` / `peer_edges`, which serve the layer
directories named by `layers_before(name)` **in LAYER_PRECEDENCE order**, so
property resolution matches the assembler's. Both substitutions were measured
before the change and are exact rather than close: the union of those
directories holds the assembled node set precisely, and phase 1 + `gtop` gives
`similarity` the identical 9,791-compound partition the assembly gave it.
Reading phase 1 **alone** would not have been neutral — 1,082 compounds are
annotated only by `gtop` — which is why the order inside phase 2 is a declared
dependency and not a convention.

### Required, optional, and the third state

`spine` and `host` are the graph; without them there is nothing to assemble
and the run stops. Every other layer may legitimately have absent inputs, and
a missing input is a **skip reported by name**, not a failure.

`pipeline.plan()` reports a third state, and leaving it out made the module's
own report wrong on the day it shipped: `assemble_release` builds its layer set
from whatever directory is **on disk**, regardless of whether the ingest ran
this time. So a skipped layer whose directory already exists is not absent —
it is **reused**. That is the more dangerous skip. An absent layer makes a
metapath reach nothing, which shows up as a zero; a reused one produces
plausible numbers from a graph that no longer matches the register.

## Assembly

`src/kgav/assemble.py`. Per-layer validation cannot see cross-layer problems;
this pass exists for the ones that only appear once the layers sit together.

**Precedence.** `spine > host > vhppi > orcs > curated > chembl > selectivity > gtop > hosttargets > similarity > ncats > indications > folds` — earlier wins a qualifier conflict.
`folds` is last because it is **computed** (a fold assignment and a gate
verdict, not an experiment). `gtop` sits above `hosttargets` because both
write `TARGETS` and only one of them knows the direction. `curated` sits below
`orcs` because it carries the same claim from weaker provenance.

**One fact, two witnesses.** A merged edge takes the union of publications,
the **strongest** evidence tier, the **earliest real** assertion date, and a
`support_count`. Two details were bugs:

- `_earliest` falls back to the `1970` placeholder only if that is all there
  is. A plain `min()` over ISO strings makes `1970-01-01` win every
  comparison, so merging a dated edge with an undated witness turned a dated
  fact into an undated one — which under `undated_policy="include"` keeps it
  in **every** training graph regardless of its real date. A temporal leak
  created by the merge, in the module whose whole purpose is the earliest
  date.
- `primary_knowledge_source` splits before unioning. It used to union the two
  *strings*, so an already-merged side carried a compound one and every
  further witness appended the whole previous string again. Merging was
  neither associative nor idempotent, a fact with many witnesses grew its
  source list quadratically, and the release held 30 distinct "sources" that
  were permutations of six real ones.

**The qualifier ladder.** Precedence cannot decide a *graded* qualifier.
Taking the higher layer's value unconditionally would have replaced ChEMBL's
`inhibitor` or `agonist` with GtoP's `modulator` on 49 edges — a stated
direction overwritten by a vaguer one, by a layer placed above it precisely
because it was assumed to know more. The ladder's three rules are arithmetic,
not policy: a claim beats no claim; specific beats graded; **specific against
opposite is neither** and becomes `contested`. That last rule costs 40 edges
their direction, and two curated sources contradicting each other is not a
fact.

`contested` is **absorbing and needs no special case**: it is absent from the
ladder, so it is an unranked specific claim, and two differing unranked claims
already fall to the irreconcilable branch. An explicit
`if CONTESTED in (hi, lo)` guard was written first and **survived mutation** —
dead code that read as load-bearing.

**Merge statistics for this release.** 30,904 edges merged;
64,424 qualifier conflicts settled by precedence alone;
14 irreconcilable.

| merge statistic | count |
|:---|---:|
| `edge_merged` | 30,904 |
| `qualifier_precedence_decided:assay_id` | 10,361 |
| `qualifier_precedence_decided:potency_max_nm` | 6,122 |
| `qualifier_precedence_decided:potency_min_nm` | 6,098 |
| `qualifier_precedence_decided:potency_spread_field` | 4,724 |
| `qualifier_precedence_decided:detection_method` | 4,676 |
| `qualifier_precedence_decided:ec50_nm` | 4,550 |
| `qualifier_precedence_decided:assay_type` | 4,197 |
| `qualifier_precedence_decided:potency_rows` | 4,037 |
| `qualifier_precedence_decided:ic50_nm` | 3,853 |
| `qualifier_precedence_decided:cell_line` | 3,846 |
| `qualifier_precedence_decided:document_id` | 3,277 |
| `qualifier_precedence_decided:screen_id` | 2,849 |
| `qualifier_precedence_decided:effect_size` | 1,816 |
| `qualifier_precedence_decided:source_efficacy_pct` | 1,632 |
| `qualifier_precedence_decided:relation` | 1,490 |
| `qualifier_precedence_decided:ki_nm` | 454 |
| `qualifier_precedence_decided:source_curve_class` | 272 |
| `qualifier_kept_more_specific:outside_dose_series` | 100 |
| `qualifier_precedence_decided:kd_nm` | 88 |
| `property_conflict:smiles` | 66 |
| `property_conflict:chembl_id` | 66 |
| `qualifier_precedence_decided:screen_type` | 63 |
| `qualifier_kept_more_specific:direction` | 49 |
| `qualifier_irreconcilable:direction` | 14 |
| `qualifier_precedence_decided:nominal_protein_target` | 9 |
| `qualifier_precedence_decided:cc50_nm` | 5 |
| `qualifier_precedence_decided:selectivity_index` | 5 |
| `property_conflict:is_approved` | 4 |
| `property_conflict:max_phase` | 2 |

A qualifier with no ladder and two different claims is settled by precedence
**silently**, and the counters above are what makes the size of that visible.
Only `direction` and `outside_dose_series` have ladders, so every measurement
qualifier — `ic50_nm`, `cell_line`, `relation`, `assay_type`, `effect_size` —
takes the higher layer's value and the other witness disappears. It is not
reconciled here: a ladder over `ic50_nm` would have to rank two numbers, and
"which assay is better" is not something the assembler can know.

**Stable bytes.** `Assembly.write` sorts records by a content-derived id and
keys *within* each record, because `merge_edges` does `out = dict(a)` and so
the arriving order decided key order in a merged edge — the same fact, the
same values, different bytes. Before that, four builds of this release
produced four `edges_sha` values and no tool in the project could say whether
the graph had changed. `nodes_sha` was stable across all four, which read as
evidence that nodes were fine; it was luck.

**Stable fact ids.** `edge_id` is `blake2b` over the edge key, so it inherits
exactly the assembler's identity semantics. It replaced a **line index**: the
agent cited `f"E{i}"`, the edge's position in `edges.jsonl`, and nothing had
ever written the field — so every citation the agent had ever emitted was a
line number, and those renumber on every rebuild.

## Three provenance stamps, one vocabulary

| stamp | file | answers |
|:---|:---|:---|
| results | the `provenance` block in every `data/results/*.json` | which graph and which code produced this number |
| derived release | `DERIVED_FROM.json` | which release this temporal split was cut from, and under which policy |
| layer | `BUILT_BY.json` | which code wrote this layer, hashing the script's whole transitive `kgav.*` import closure |
| source | `SOURCE.json` | which *release of the database* the layer read |

All four use the same four statuses — `ok`, `stale`, `unstamped`, `absent` —
because a third vocabulary for one idea is its own defect. They exist because
a results file recorded numbers and nothing about where they came from: a
comparison's "before" column was five months stale, and SARS-CoV-2's M1
positive count appeared to **rise** 736 → 889 when the true movement on one
release was 1,011 → 889, a fall. That error was undetectable by inspection.

**The `dirty` flag reported on its own output**, twice. It was computed from a
bare `git status` over the whole tree, and `data/results/*.json` is tracked —
so the first file a run wrote dirtied the tree and every stamp after it
recorded `dirty: True`. 13 of 14 files from one regeneration said so whatever
happened. One directory over, the release `MANIFEST.json` is the only tracked
file under `data/releases`, written in one phase and stamped in another, so it
was **always** true by construction: a manifest records HEAD's hash and is
committed afterwards, so its committed version necessarily names the previous
commit. Measured: 7 of 7 committed versions said `dirty: True`, each naming
the commit before its own. This release's stamp reads
`e2bc1305fb19`, dirty `False`.

## Validation

```
make release-validate     # the assembled release against the frozen schema
make gate                 # the NEGATIVE CONTROL
```

`make gate` runs the validator against `data/releases/smoke`, which carries
four deliberately planted violations, and asserts the exact codes and counts:
`CONSTRAINT=1,MISSING_PROVENANCE=1,MISSING_QUALIFIER=2`. It used to be *"if
the validator fails, OK"* — which accepted an `ImportError`, a missing
directory or an unloadable schema as proof the planted errors were caught. The
one thing it could not distinguish was a working validator from a broken one.
