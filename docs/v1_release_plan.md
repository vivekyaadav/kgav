# v1: the pathway-analysis release

Written 2026-10-24, from the state of release v0.1 on that date. Everything
below is either measured by a script in this repository or a decision
recorded with its reason. Where a figure appears, the command that produces
it is named, so this document can be checked rather than believed.

## What v1 is

The pathway analysis over every small molecule the graph holds, shipped
with the agent in `src/kgav/agent/`, and the paper written from it. A
scoring method for drug repurposing comes after, not in this release.

## Where the graph stands

Release v0.1: 78,789 nodes, 392,375 edges, schema 0.14.0, `schema: PASS`.
Sixteen viruses — seven coronaviruses and nine flaviviruses — after the
2026-10-23 promotion took the built set from eight.

    make data-audit        what the schema declares against what exists
    make data-timeline     when each input arrived, and what period it covers
    make rebuild           the ingest chain in dependency order
    make check-results     which results files describe a graph that is gone

### The finding the paper is about

`scripts/pathway_landscape.py`, on v0.1: drugged human proteins overlap
CRISPR dependency factors at close to the rate a degree-matched random
protein set does, in all eight taxa that have dependency factors. **No
taxon clears the corrected threshold in either direction.**

<!-- TABLE: landscape_by_taxon — asserted cell-by-cell by
     tests/test_release_plan_figures.py against
     data/results/pathway_landscape.json. ratio is n_observed/expected,
     derived the way kgav.landscape.TaxonOverlap.ratio derives it. Do not
     hand-edit; re-derive with `make pathway-landscape`. -->

| virus | dependency factors | observed | matched expectation | ratio | p naive | p matched |
|---|---:|---:|---:|---:|---:|---:|
| HCoV-NL63 | 57 | 10 | 7.4 | **1.35** | 0.04148 | 0.1988 |
| SARS-CoV-2 | 3,460 | 384 | 358.4 | **1.07** | 0.0003959 | 0.05964 |
| HCoV-229E | 1,367 | 128 | 135.3 | **0.95** | 0.6012 | 0.7675 |
| HCoV-OC43 | 657 | 62 | 65.2 | **0.95** | 0.5537 | 0.6843 |
| Zika | 815 | 76 | 80.1 | **0.95** | 0.5994 | 0.7093 |
| MERS-CoV | 678 | 60 | 66.4 | **0.90** | 0.7494 | 0.8206 |
| DENV-2 | 50 | 2 | 6.3 | **0.32** | 0.9583 | 0.9898 |
| West Nile | 18 | 0 | 2.5 | **0.00** | 1 | 1 |

**THIS PARAGRAPH CLAIMED A DEPLETION IN EVERY TAXON AND THAT WAS FALSE**
(fix-list D4). It read "0.86–0.94 times the rate ... in every one of the
eight taxa ... with matched p from 0.73 to 1.00", and on this release the
ratios run **0.00 to 1.35** with **two taxa above 1.0** — HCoV-NL63 at 1.35
and SARS-CoV-2 at 1.07 are ENRICHED, not depleted — and matched p runs from
**0.05964**, not 0.73. The naive p for SARS-CoV-2 is 3.96e-04, not 1.28e-04.
Every figure came from the 2026-10-28 run and nothing re-derived them
across four rebuilds; this table is now machine-checked.

The conclusion survives the correction and gets weaker, which is the
honest direction. **SARS-CoV-2's naive p of 3.96e-04 does disappear under
degree matching** — 0.05964, against a Bonferroni-corrected alpha the
landscape computes from its own family size — so the headline control still
works. What cannot be said any more is that the overlap is *depleted*: it
is indistinguishable from the degree-matched background, in both
directions, in all eight taxa.

**DENV-2 at 0.32 and West Nile at 0.00 are still NOT evidence of
depletion.** They hold 50 and 18 dependency factors against matched
expectations of 6.3 and 2.5, so a ratio there is one or two Poisson draws;
West Nile observed zero. They contribute scope, not statistical weight, and
the paper must say so. **HCoV-NL63 at 1.35 is the same situation in the
other direction** — 10 observed against 7.3 expected — and the old
paragraph had no room for that case because it asserted every taxon fell
in one narrow band.

---

## A. Blocking — a false claim ships without these

### A1. Decide what `docs/evaluation.md` claims

`make check-results` reports 14 stale files: every baseline, every temporal
split, the hard negatives and the metapath comparison were computed on the
eight-virus graph. `docs/evaluation.md` publishes figures derived from
them.

This is a decision, not a job. Either regenerate — the Reproducing block
now loops over all sixteen taxa, and six of the new files will record empty
host-mediated channels by construction — or scope the document explicitly
to the scoring work and state that v1 is the pathway release. What cannot
stand is a document stating numbers that no longer describe the graph.

### A2. `host` layer provenance — **DONE 2026-11-11, and the report was the last part to be fixed**

The item as written: `make data-timeline` reports the `host` layer as five
files, 328,846,356 bytes, `mtime only`. That is the human proteome and
STRING v12.0, feeding 157,211 `PARTICIPATES_IN` and 141,019
`PHYSICALLY_INTERACTS_WITH` edges — over three-quarters of the graph — with
no recorded URL, release or checksum.

**The checksums were written on 2026-10-25** by
`scripts/record_raw_provenance.py` (`make record-provenance`), which stamps
a hand-placed file with its sha256 and states `url: null`,
`retrieved: null`, `provenance: "reconstructed"` rather than inventing a
download it cannot verify.

**AND `make data-timeline` WENT ON PRINTING `mtime only` FOR ANOTHER
SEVENTEEN DAYS.** The host layer declares the glob `data/raw/host/*`, so
each `.source.json` this project wrote was counted as an input: ten files
instead of five, 328,848,975 bytes instead of 328,846,356, and the newest
of the ten was a sidecar — they were written that day, the data files are
dated 2023 and September — and a sidecar has no sidecar of its own. **The
report whose job is to say whether provenance exists was driven negative by
the act of recording it.** The figures quoted in this item were always the
right ones; the report moved away from them. Fixed in
`data_timeline.inputs_for`, where a `.source.json` is no longer an input,
and the verdict is now taken over every input rather than whichever is
newest. `tests/test_data_timeline_provenance.py` holds both halves.

**Every layer the pipeline declares now reads `sidecar`**, including the
three that share `chembl_37.db` — 30,480,314,368 bytes, hashed in 113s.
`curated` and `similarity` read `-` because they have no input on disk.

What is still NOT done, and is a different job: a fetcher that lets a fresh
clone download these with a real URL and a real retrieval date, the way
`fetch_raw.py` does for proteomes. The sidecars make what is here
auditable; they do not make it reproducible elsewhere.

### A3. The scope statement, written down

All of it is measured, and the honesty is the contribution:

- Sixteen viruses built; **eight have host-factor data**, eight do not.
  ORCS carries screens for five coronaviruses, Zika, DENV-2 and West Nile.
  DENV-1, DENV-3, DENV-4 and yellow fever have none in this release. **JEV
  and TBEV have no genome-wide screen in the literature at all** — a limit
  of the field, not of the build.
- Each flavivirus is **one strain's polyprotein** reassigned to the species
  taxon. "DENV-1 NS3" in this graph is strain 408685's NS3. This is what
  `chain_source` has always meant.
- **CLOSED BY gtop.** 12,155 of 14,381 `TARGETS` edges stated no
  direction (84.5%) before that layer landed; on the assembled v0.1 it is
  10,395 of 17,982 (57.8%), and the directed count rose 2,226 -> 7,587.
  The original statement, kept because the reasoning still holds: so
  "this drug targets a dependency factor" is not "inhibiting it would
  help", and the landscape says so in its own output.
- The landscape covers **11.8% of molecules**; the rest have no human
  target in this graph. The ceiling a curated drug-target source could add
  is 1,679, not the full 22,381.
- HCoV-HKU1 is built and carries **eight edges and nothing else** — no host
  factor, no antiviral activity. It answers no question and should be named
  rather than left for a reader to find in a node count.

## B. Blocking for the agent

### A4. The agent must not imply a ranking is predictive

`GraphTools.triage` ranks compounds and `explain` finds paths. The measured
result is that the pathway route carries no enrichment in any taxon. An
agent that ranks without saying so is making the claim the landscape
refuses.

This repository has been here before — `7ac7090` ("The agent could explain
a route the scorer refuses") and `c7ebf26` ("The agent explained a route
its own source contradicts") — so the guard pattern exists in
`src/kgav/agent/guards.py`. It needs a test, not a sentence in prose: a
caveat printed under a fluent answer is read second, if at all.

### A5. A pathway capability

`explain` finds paths between a compound and a virus. Nothing answers the
v1 question for an arbitrary compound: which pathways do this drug's
targets sit in, and does the virus depend on any member of them. That is
`pathway_landscape.py`'s computation at the level of one drug, over the
3,002 molecules that carry a human target.

## C. Improves v1, does not block it

- **C1. Ingest IntAct.** `data/raw/vhppi/intact_viral.tab27` is written
  (139 MB, 12,496 rows) and read by `ingest_vhppi` as `infores:intact`. It
  lands at the next `make rebuild`. Gives M3 six flaviviruses instead of
  zero. Thin — 121 flavivirus rows, and Zika has three against Shah et
  al.'s 368, so IMEx has not curated that paper.
- **C2. Guide to PHARMACOLOGY. DONE.** `src/kgav/gtop.py` and
  `scripts/ingest_gtop.py` are built, tested AND BUILT INTO v0.1:
  `docs/kg_release_census.txt` lists `gtop` among the 11 assembled layers
  with 6,765 edges from `infores:guidetopharmacology`. This entry said
  "the files are not downloaded" after they had been, which also means
  the `direction: unknown` fraction it promised to move has already
  moved -- see scope_and_limits.md for the current figure. OPEN: the
  GtoPdb licence question is use-versus-access, and the edges are
  already in the release, so it needs settling rather than planning.
- **C3. MONDO terms for 13 of 16 viruses.** Three Disease nodes today, so
  `TREATS` and `IN_TRIAL_FOR` — both marked `evaluation_label` in the
  schema — are empty, and every evaluation rests on
  `HAS_ANTIVIRAL_ACTIVITY_AGAINST` alone. ChEMBL's `drug_indication` table
  fills both and is already on disk in `chembl_37.db`; the layer queries
  `activities` and `assays` and never touches it. Blocked behind the
  Disease axis. Not needed for the pathway paper.

## D. Out of scope for v1

- **The temporal split and everything downstream.** Decided 2026-10-24: the
  pathway analysis uses the whole graph.

  **THIS NOTE WAS WRONG IN FOUR WAYS UNTIL 2026-11-11, AND C1/C2/V4 HAD
  ALREADY FIXED OR MEASURED EVERY ONE OF THEM.** It read: *"all 150,432
  `PHYSICALLY_INTERACTS_WITH` edges are dropped from any pre-2021 training
  cut, because `ingest_host` stamps every STRING edge with
  `string_release_date` (2023-08-28) and `ingest_vhppi` stamps every edge
  with its `--date` default (2024-01-01). M3 and M5 therefore have zero
  edges in the temporal protocol."* Against the release and the temporal
  results:

      153,533   PHYSICALLY_INTERACTS_WITH edges in v0.1
      141,019   STRING, host-host, dropped from the 2021 cut
       11,744   SURVIVE the cut, and every one is virus-host
            0   of the survivors are STRING's

  The quotation above is what was claimed; the block is what is true. No
  figure from the superseded sentence is restated outside that quotation,
  which is the convention `tests/test_release_plan_temporal_note.py`
  enforces -- it reads the plan with quoted material removed, so a wrong
  number recorded as history passes and the same number asserted does not.

  **"All" and "zero" were the load-bearing words and both were false.** The
  11,744 come from VirHostNet (5,222), IntAct and its overlaps (6,362),
  MINT (147), DIP (12) and HPIDB (1), and they sum exactly. So **M3, the
  virus-host middle hop, RUNS in the temporal protocol and only M5,
  host-host, cannot** -- which is the reason M5 was absent from every
  SARS-CoV-2 temporal table ever published here, a fact this note's "M3 and
  M5 therefore have zero edges" actively concealed.

  The two mechanisms it blamed are also no longer there. `ingest_vhppi`
  dates each edge from its own PMID year (`vhppi.PUB_YEAR_RE`, fix-list
  C1), so the `--date` default is a fallback and not the stamp; and
  `string_release_date` was renamed `string_file_mtime` (V4) precisely
  because it is a download mtime and not a release date.

  **What survives of the note is the recommendation, and it is unchanged:**
  STRING has no per-edge date and should carry `PLACEHOLDER_DATE`, which
  routes it through `undated_policy` rather than silently outside the
  cutoff. That is a decision, not a correction -- it would move all 141,019
  edges from one side of the cut to "undated" -- and it is fix-list C2,
  which also records that STRING v11.0 (2019) is a 139 MB download and
  `data/raw` is hand-populated. None of this removes an edge from v0.1.

  The figures above are re-derived by
  `tests/test_release_plan_temporal_note.py` from the release and from
  `data/results/temporal_2021_include.json`'s
  `training_graph_predicate_sources`, which `temporal_split.py` writes for
  this purpose. They are not retyped here again.
- **The structure-derived scorer — REGISTERED 2026-11-11, not run.**
  `docs/prereg_structure_scorer.md`, signed before the code and before the
  measurement. It is out of scope for v1 because v1 is the pathway analysis
  and this is a scorer, but it is the one candidate measured to have room to
  rank: every evaluated channel sits within +0.0050 of its imputed floor
  because 74.8% of labelled compounds are unreachable, while all 13,417
  labelled compounds in the registered population carry a SMILES. Fourteen
  readable cells, alpha = 0.05/14. Nothing is computed yet, and the
  registration can still be amended until the first score exists.
- The resistance axis: `VirusStrain`, `Variant`, and `HAS_STRAIN`,
  `OCCURS_IN`, `HAS_VARIANT`, `CONFERS_RESISTANCE_TO`.
- `COMBINED_WITH`, `SUBCLASS_OF`, `SUBCLASS_OF_DISEASE`.
- Curating published flavivirus screens into `src/kgav/curated.py`. The
  only route to more rows in the overlap table, and a large job. JEV and
  TBEV stay out of reach regardless.

## Order

    A2  →  A3  →  A4  →  A5  →  release
    DONE

with C1 folded into whichever rebuild comes next. A1 is a decision to take
before writing, not work to schedule.

**A2 closed 2026-11-11.** C1 also landed: `infores:intact` contributes
6,617 edges on the assembled release, so it is no longer waiting for a
rebuild. The three remaining items are the agent and scope ones, and A3's
figures have moved since it was written — the drugged-proteome share it
quotes as 11.8% is 16.3% on this release, which
`docs/scope_and_limits.md` now derives by test rather than restating.
