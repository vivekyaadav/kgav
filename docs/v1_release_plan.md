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

Release v0.1: 78,782 nodes, 386,192 edges, schema 0.13.0, `schema: PASS`.
Sixteen viruses — seven coronaviruses and nine flaviviruses — after the
2026-10-23 promotion took the built set from eight.

    make data-audit        what the schema declares against what exists
    make data-timeline     when each input arrived, and what period it covers
    make rebuild           the ingest chain in dependency order
    make check-results     which results files describe a graph that is gone

### The finding the paper is about

`scripts/pathway_landscape.py`, on v0.1: drugged human proteins overlap
CRISPR dependency factors at 0.86–0.94 times the rate a degree-matched
random protein set does, in every one of the eight taxa that have
dependency factors, with matched p from 0.73 to 1.00. The naive p of
1.28e-04 for SARS-CoV-2 disappears under degree matching.

Eight taxa, two families. Zika at 0.88 is not a repetition of the
coronavirus result; it is a different family, a different screen
literature, different labs, landing in the same band.

Two taxa — DENV-2 at 0.26 and West Nile at 0.00 — are NOT evidence of
depletion. They have 50 and 18 dependency factors against expectations of
7.5 and 3.1, so a ratio there is one or two Poisson draws. They contribute
scope, not statistical weight, and the paper must say so.

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
- **12,155 of 14,381 `TARGETS` edges state no direction** (84.5%). So
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
- **C2. Guide to PHARMACOLOGY.** `src/kgav/gtop.py` and
  `scripts/ingest_gtop.py` are built and tested; the files are not
  downloaded. The only thing that moves the 84.5% unknown direction.
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
