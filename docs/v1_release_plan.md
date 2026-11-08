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

Release v0.1: 78,789 nodes, 392,363 edges, schema 0.13.0, `schema: PASS`.
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
| HCoV-NL63 | 57 | 10 | 7.3 | **1.36** | 0.04148 | 0.1869 |
| SARS-CoV-2 | 3,460 | 384 | 359.2 | **1.07** | 0.0003959 | 0.07096 |
| HCoV-OC43 | 657 | 62 | 65.3 | **0.95** | 0.5537 | 0.6887 |
| HCoV-229E | 1,367 | 128 | 134.9 | **0.95** | 0.6012 | 0.7596 |
| Zika | 815 | 76 | 80.4 | **0.95** | 0.5994 | 0.7241 |
| MERS-CoV | 678 | 60 | 66.4 | **0.90** | 0.7494 | 0.8256 |
| DENV-2 | 50 | 2 | 6.4 | **0.31** | 0.9583 | 0.989 |
| West Nile | 18 | 0 | 2.5 | **0.00** | 1 | 1 |

**THIS PARAGRAPH CLAIMED A DEPLETION IN EVERY TAXON AND THAT WAS FALSE**
(fix-list D4). It read "0.86–0.94 times the rate ... in every one of the
eight taxa ... with matched p from 0.73 to 1.00", and on this release the
ratios run **0.00 to 1.36** with **two taxa above 1.0** — HCoV-NL63 at 1.36
and SARS-CoV-2 at 1.07 are ENRICHED, not depleted — and matched p runs from
**0.071**, not 0.73. The naive p for SARS-CoV-2 is 3.96e-04, not 1.28e-04.
Every figure came from the 2026-10-28 run and nothing re-derived them
across four rebuilds; this table is now machine-checked.

The conclusion survives the correction and gets weaker, which is the
honest direction. **SARS-CoV-2's naive p of 3.96e-04 does disappear under
degree matching** — 0.071, against a Bonferroni-corrected alpha the
landscape computes from its own family size — so the headline control still
works. What cannot be said any more is that the overlap is *depleted*: it
is indistinguishable from the degree-matched background, in both
directions, in all eight taxa.

**DENV-2 at 0.31 and West Nile at 0.00 are still NOT evidence of
depletion.** They hold 50 and 18 dependency factors against matched
expectations of 6.4 and 2.5, so a ratio there is one or two Poisson draws;
West Nile observed zero. They contribute scope, not statistical weight, and
the paper must say so. **HCoV-NL63 at 1.36 is the same situation in the
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

### A2. `host` layer provenance

`make data-timeline` reports the `host` layer as five files, 328,846,356
bytes, `mtime only`. That is the human proteome and STRING v12.0, feeding
157,211 `PARTICIPATES_IN` and 141,019 `PHYSICALLY_INTERACTS_WITH` edges —
over three-quarters of the graph — with no recorded URL, release or
checksum.

`fetch_raw.py` writes a `.source.json` sidecar for every proteome it
fetches. Nothing does this for the host layer, which was placed by hand.
This is the first thing a reviewer asks about.

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
  pathway analysis uses the whole graph. Note for whenever it is picked up
  again — all 150,432 `PHYSICALLY_INTERACTS_WITH` edges are dropped from
  any pre-2021 training cut, because `ingest_host` stamps every STRING edge
  with `string_release_date` (2023-08-28) and `ingest_vhppi` stamps every
  edge with its `--date` default (2024-01-01). M3 and M5 therefore have
  zero edges in the temporal protocol, and `edge_dropped_future` is a
  single scalar that cannot show it. `vhppi.parse_row` already extracts a
  per-edge PMID and `ingest_mitab` discards it. STRING has no per-edge date
  and should carry `PLACEHOLDER_DATE`, which routes it through
  `undated_policy` rather than silently outside the cutoff. None of this
  removes an edge from v0.1.
- The resistance axis: `VirusStrain`, `Variant`, and `HAS_STRAIN`,
  `OCCURS_IN`, `HAS_VARIANT`, `CONFERS_RESISTANCE_TO`.
- `COMBINED_WITH`, `SUBCLASS_OF`, `SUBCLASS_OF_DISEASE`.
- Curating published flavivirus screens into `src/kgav/curated.py`. The
  only route to more rows in the overlap table, and a large job. JEV and
  TBEV stay out of reach regardless.

## Order

    A2  →  A3  →  A4  →  A5  →  release

with C1 folded into whichever rebuild comes next. A1 is a decision to take
before writing, not work to schedule.
