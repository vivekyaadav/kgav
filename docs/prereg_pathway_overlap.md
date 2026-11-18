# Pre-registration: the pathway and interactome overlap, re-run

**SIGNED 2026-11-11, RETROSPECTIVELY, AND THE ORDERING IS RECORDED BELOW
RATHER THAN SMOOTHED OVER.** The values under *Fixed in advance* are now
fixed, and the three decisions flagged **DECISION** -- the permutation
count, which compounds the drug-target query admits, and whether the six
empty-rejection-region tests stay in the correction family -- are decided as
this document proposes them, unchanged. A fourth candidate -- whether
`auc_ceiling` replaces `MIN_REACHED_POSITIVES` as the readability gate --
belongs to the hard-negatives analysis and not to this one, and is left out
rather than listed here for completeness.

**THE SIGNATURE POSTDATES THE RUN, WHICH THIS DOCUMENT'S OWN DRAFT
FORBADE.** The draft read *"this is not a registration until it is signed
and committed, and the re-run must not be executed before then"*. The re-run
was executed on **2026-11-11 at 04:22:49 UTC** and signed hours later, so
that sentence was broken rather than satisfied, and no wording here repairs
it. A weaker claim is true and is checkable, and it is the whole of what this
signature asserts:

> **Every registered value predates the run. The signature does not.**

**AND THE RUN HAS SINCE BEEN REPEATED TWICE, ON TWO LATER GRAPHS.** The
artifact on disk was generated at **2026-11-17 at 15:00:19 UTC**
against the clinical-label release (`edges_sha` 873623c9a9df2049,
78,649 nodes **396,223 edges**), which is the run every derived figure
below is now read against. The run before it, at 11:10:53 UTC, was
against the entry-assay release and is superseded; its figures are not
quoted anywhere below. That does not move the signature or any registered value:
the protocol, the family, the permutation count, the alpha and the
empty-rejection decision were all fixed before 2026-11-11 and are
unchanged. What moved is the expectations the graph produces, and the
figures that are functions of them are re-derived by
`tests/test_landscape.py` rather than maintained here. **The 04:22:49 run's
ordering problem stands exactly as recorded**; repeating the run on a new
graph does not retroactively make the signature precede it.

That is a narrower guarantee than the draft asked for and it is worth
stating what it still rules out: no parameter was chosen, and no test
dropped from the correction family, after a p-value was in hand. What it
does NOT rule out is the one thing only an ante-dated signature could --
that signing was itself conditional on having seen the result. Nothing in
this repository can settle that, so it is named here and left standing.

| registered value | fixed by | relative to the run |
|---|---|---|
| Permutations: 100,000 | `b443f4b`, 2026-11-02 | 9 days before |
| The six empty-rejection-region tests stay in the family | `b443f4b`, 2026-11-02 | 9 days before |
| Correction family: every test, both channels, both directions | `b443f4b`, 2026-11-02 | 9 days before |
| Drug-target admission: `floor` | `0c2d3e4`, 2026-11-10 | 1 day before |
| **This signature** | 2026-11-11, after the run | **after** |
| The re-run on the L27-repaired release | 2026-11-11 11:28:33 UTC | **after the signature** |
| The L26 re-run, recording the operand | 2026-11-16 23:04:38 UTC | **after the signature** |

**AND THE RESULT NOW ON DISK IS THE ONE RUN IN THE RIGHT ORDER.** The run
this signature argues about was executed at **04:22:49 UTC** and signed
afterwards. The L27 identity repair then rebuilt the release, and
`data/results/pathway_landscape.json` was regenerated at **11:28:33 UTC** --
*after* the signature, which is the order the draft asked for and the first
time this document has had it. Every registered value predates that run by
construction, since they were all fixed before the signature, so the
re-run's figures are the ones the registration can claim without
qualification. **The file on disk is now a third run, at 2026-11-16 23:04:38
UTC**, made to close fix-list L26 by recording the operand; it is later
still, so the same ordering argument covers it, and `compare_results.py`
shows it moved no figure -- only `measured_on` and the new `operands`
key. The 04:22:49 run's ordering problem is kept recorded above
because it happened, not because it still governs the figures below. The
operand it moved is recorded under *A registered operand moved*.

**ONE EDIT LANDED AFTER THE RUN AND IT IS NOT A REGISTERED VALUE.** Commit
`d43c316` (2026-11-11 10:01:55 +05:30, nine minutes after the run) changed
this file by exactly one line, `edges_sha e26239a4e87865b3` ->
`ad7121bc2195a93a` in the release identification below. So at the moment the
run executed, this document named a release the run did not read -- the
release had already moved and the identification had not caught up. It is a release
identification rather than a parameter, and the whole of it is visible in
`git show d43c316 -- docs/prereg_pathway_overlap.md`: one insertion, one
deletion.

**THAT PAIR OF HASHES IS A QUOTATION FROM THE HISTORY, AND I OVERWROTE IT.**
Commit `5ab4af1` moved every live release identification in this directory to
the L27-repaired release, and it also replaced the `ad7121bc2195a93a` in the
sentence above with `545fb7d87b0aec08` -- which turned a true statement about
what `d43c316` did into a false one, because `d43c316` wrote
`ad7121bc2195a93a` and `5ab4af1` itself wrote `545fb7d87b0aec08`. The same
commit's message claimed the opposite in as many words: *"87 occurrences of
the five old values sit across 18 files and most are deliberate history ... a
global replace would have destroyed that, so each edit was driven by a
failing test."* No test covered this sentence, so no failing test drove it and
nothing objected. `tests/test_prereg_is_obeyed.py::
test_the_quoted_post_run_edit_is_the_diff_git_shows` now reads the two hashes
out of this sentence and checks them against `git show d43c316`, so the next
such replacement fails.

**THE RUN'S CODE STAMP SAYS `dirty: true`, AND IT DID NOT REACH THIS
ANALYSIS.** `pathway_landscape.json` records code commit `21c0f92` with a
dirty working tree, because the run happened mid-session with the
`potency_rows` rename uncommitted. Both files this analysis consists of are
byte-identical between that commit and the signed tree --
`scripts/pathway_landscape.py` at sha256 `39d20e56b370`,
`src/kgav/landscape.py` at `887f0405372f` -- and no commit from `21c0f92` to
the signature touches either. The uncommitted content was elsewhere.

**ONE REGISTERED PARAMETER WAS NOT RECOVERABLE FROM THE ARTIFACT, AND IS
NOW -- CLOSED 2026-11-17.** The paragraph below is kept as it was written,
because it is what this signature recorded about itself; what follows it is
the closure. The
operand this document fixes -- the drugged proteome under `floor`, 1,939
human proteins -- is recorded in **neither** `pathway_landscape.json` nor its
provenance block, and it is the one parameter this document says the verdict
depends on (0 of 22 rows clearing Bonferroni under `floor`, 3 under
`relaxed`). A reader cannot confirm from the result file which operand
produced it. It is re-derivable from the release, and re-deriving it on the
current release gives **1,944 drugged of a 20,416-protein human universe**,
which matches what is registered below; `tests/test_prereg_is_obeyed.py`
pins that so the agreement is checked rather than claimed. **At signature it
was 1,946 and the L27 repair moved it**, which is recorded below under *A
registered operand moved*. **The artifact
should state it, and that is the open item this signature leaves behind:**
recording it requires regenerating the file, which is another 100,000
permutations, and doing so AFTER this signature is the correct order for the
first time in this document's history.

**CLOSED 2026-11-17, AND THE ORDER WAS KEPT.** `pathway_landscape.json` now
carries an `operands` block among its study parameters -- not in its
provenance block, because these are parameters of the analysis and not of
the build. It records the drugged proteome as **n = 1,944, a sha256 over
its sorted ids, and the 1,944 ids themselves**, and the human universe as
**n = 20,416 with its own sha256**. Three values are checked against each
other rather than asserted: `tests/test_prereg_is_obeyed.py` reads the
figure out of the DECISION paragraph below and out of the JSON and compares
them, and its slow companion re-derives both sets from the release, so
document, artifact and graph have to agree and no two of them can agree
while being stale together.

**THE MEMBERS ARE STORED BECAUSE A COUNT WAS MEASURABLY NOT ENOUGH.** *A
registered operand moved*, below, records that the L27 rebuild took this
operand from 1,946 to 1,939 and that **which seven proteins left could not
be answered**, the superseded release having been rebuilt in place. With
the ids on disk that question is a set difference. The universe's 20,416
ids are deliberately NOT stored -- 400 KB against the drugged set's 38 KB,
in a file that was 9 KB -- and the omission is named in the file itself
rather than left to be discovered, with the fingerprint retained so a move
in it is still loud.

**THE ADMISSION POLICY IS STILL NOT NAMED, and that is deliberate rather
than overlooked.** Neither this script nor the release manifest records
whether `floor` or `relaxed` built the graph, so a file asserting `"floor"`
would be stating something it cannot check -- the defect this document
exists to avoid. The operand identifies it instead: 1,944 under `floor`
against 2,489 under `relaxed`, which is exactly the discrimination the
reader was missing.

**AND NOTHING ELSE IN THE FILE MOVED.** The re-run was compared leaf by
leaf against the superseded artifact with `scripts/compare_results.py`,
which drops the provenance stamp and reports every changed value: the only
differences are `measured_on` and the new `operands` key. Every p-value,
every taxon row and every interactome figure is identical, so this closure
adds a parameter and does not restate a result.

Date drafted: 2026-11-02. Release: v0.1 (78,649 nodes / 395,646 edges,
`nodes_sha` dc6d76a2b3cfb2a1, `edges_sha` c0c4daccd48e9979). Scripts:
`scripts/pathway_landscape.py`, read by `src/kgav/landscape.py`.

**The release identification has now been updated thirteen times -- once on
2026-11-07, three times on 2026-11-08, five times on 2026-11-09, three
times on 2026-11-11 and once on 2026-11-17 -- and the analysis was not
re-registered any of them**,
because no rebuild has changed a measured p-value. One rebuild has now
changed a registered OPERAND, which is recorded below under *A registered
operand moved* and is a different thing from a measured figure. **That count was wrong until
2026-11-11 and so was its breakdown.** It read "eight times -- 2026-11-08
and seven times on 2026-11-09", which was maintained by hand and matched
neither the total nor any date: git puts three of them on 2026-11-08, not
one, and five on 2026-11-09, not seven. It is now taken from the history,
derived by `tests/test_prereg_is_obeyed.py` rather than maintained, and
anyone can re-take it --

    git log -G'nodes_sha. [0-9a-f]{16}' -- docs/prereg_pathway_overlap.md

-- where the thirteenth and oldest commit it lists is the drafting rather
than an update. **The pattern has to pair the key with an actual hash.** Written
as `nodes_sha.*edges_sha` it matched the line above that DOCUMENTS it, so the
commit adding this paragraph counted itself as a twelfth update and the test
failed on its own commit. **And the obvious command for this is the wrong
one.** `git log -L
/^Date drafted/,+1:` tracks the line by content and silently omits
`d43c316`, the one commit this header is about, because it re-resolves the
regex as it walks; the line-number form `-L 14,15:` agrees with the pickaxe
above on all twelve but needs the line number the file had at the revision
being read, which this signature's own insertion has already changed. The narrative below numbers only the first eight and stops; the
ninth and tenth are the 2026-11-11 rebuilds (the identity-policy fix and the
`potency_rows` rename), the eleventh is the `edges_sha` correction described
in the header, which followed the run by nine minutes, and the twelfth is
`5ab4af1`, which moved this identification onto the L27-repaired release
(78,500 nodes / 391,727 edges). The twelfth is the first of them that moved a
registered operand, and it is the one described under *A registered operand
moved*. The fifth and
sixth are P13, which added four qualifiers to `TARGETS` edges and removed
them again from 551 merged ones, and P14, which gave the gtop layer a range
of its own so that 16,436 of 17,982 edges keep one; this analysis reads
`TARGETS` as a SET of (compound, protein) pairs and no qualifier of it, so a
qualifier-only change cannot reach it for the same reason the date-only
change of P8 could not. The seventh is E6, which records `fold_mode` and
`fold_strength` on the 24 `FOLD_SIMILAR_TO` edges, and the eighth is L9,
which records what a ChEMBL curator wrote beside the number on the three
label predicates -- neither predicate is one this analysis reads, and both
are qualifier-only changes.
`data/results/pathway_landscape.json` is regenerated against each new
release and differs from the previously committed version only in its
**provenance and its own run date**: the code commit, the timestamp,
`edges_sha` and `measured_on`. Every p-value, every overlap ratio and every
degree-matched null is byte-identical, checked field by field rather than
asserted.

The fourth rebuild is the first one that touched a predicate this analysis
DOES read. P8 moved `first_asserted_date` on 453 `TARGETS` edges -- 350
earlier, 103 off the 1970 placeholder -- and the figures did not move,
because this analysis reads the drugged proteome as a SET and dates nothing:
it asks which human proteins carry a `TARGETS` edge at all, not when the
edge was asserted. The earlier rebuilds moved only
`HAS_ANTIVIRAL_ACTIVITY_AGAINST` (two edges added by L4, 261 removed by L5),
which it does not read. A registration that quoted a release no longer on disk would be
unverifiable, which is what `tests/test_landscape.py::
test_the_quoted_release_matches_the_census` exists to prevent; it caught
this.

## A registered operand moved

**2026-11-11, AFTER THE SIGNATURE.** Fix-list L27 repaired the compound
identity policy -- node ids were ChEMBL's `standard_inchi_key` verbatim while
both SQL layers declared `salt_collapsed=True`, so a salt and its parent were
two nodes -- and the release was rebuilt: 78,789 nodes / 392,375 edges ->
78,500 / 391,727. **One registered operand changed with it. The drugged
proteome under `floor` is 1,939 human proteins, not 1,946.** The human
universe is unchanged at 20,416.

**THE REGISTERED RULE DID NOT CHANGE AND IS NOT BEING RE-REGISTERED.** What
this document fixes is the admission policy -- `--phase-gate floor` -- and
1,939 is that rule's value on this release exactly as 1,946 was its value on
the signed one. Re-registering would be the wrong move for the same reason
the eleven release identifications before this one were not re-registered: a
registration that had to be re-signed whenever the graph was repaired would
make every repair a protocol change. **What has to be true for that reasoning
to hold is that the operand is derived and not chosen, and that is checked:**
`tests/test_prereg_is_obeyed.py::test_the_registered_operand_is_the_one_the_release_gives`
loads the release and recomputes it, so this figure cannot drift from the
release by hand, in either direction.

**THE VERDICT DID NOT MOVE.** The re-run on the repaired release gives **0 of
22 rows clearing Bonferroni**, the same as the signed run, at the same
alpha = 0.05/44. The largest effect is still interactome/NCBITaxon:11082 at
`p_matched` 0.006210 against a threshold of 0.001136, and SARS-CoV-2's
interactome row moved 0.235 -> 0.2237. So the direction of this analysis's
conclusion does not depend on the seven proteins.

**AND THE SEVEN PROTEINS CANNOT BE NAMED, WHICH IS THIS SIGNATURE'S OWN OPEN
ITEM COMING DUE.** The header records that `pathway_landscape.json` stores the
permutation count, the seed scheme, the degree basis and the family size and
**not** the drugged proteome, and that the operand is therefore recoverable
only by re-deriving it from the release. The superseded release was rebuilt in
place, so it is not on disk and the old operand is no longer derivable at all:
`1,946` survives only as the number written in this document. Which seven
proteins left, and by what mechanism, is therefore **not answerable from this
repository** -- and had the artifact recorded the operand as that open item
says it should, it would have been a set difference. A gap that was an
abstract criticism on the day it was written cost a concrete answer nine
hours later, and that is the whole argument for closing it.

**2026-11-17: THE OPERAND MOVED AGAIN, AND THIS TIME THE PROTEINS CAN BE
NAMED.** The entry-assay ingest added three NCATS pseudotyped-particle
screens and rebuilt the release (78,500 / 391,727 -> 78,649 / 395,646).
**The drugged proteome under `floor` is now 1,944 human proteins, not
1,939.** Five joined and **none left**:

    UniProtKB:Q16280   UniProtKB:Q4FZB7   UniProtKB:Q86Y97
    UniProtKB:Q92993   UniProtKB:Q9H7Z6

The human universe is unchanged, and provably so rather than apparently:
its sha256 is the same `f32b6e78f99a96a5...` the previous artifact
recorded. **This is the closure above earning its cost on the first
rebuild after it landed.** Nine hours of work bought the ability to answer
"which proteins" with a set difference instead of "not answerable from this
repository" -- and the answer came out five-joined-none-left, which is the
shape a widening-by-addition should have and a repair-by-collapse should
not.

**THE RULE IS STILL NOT RE-REGISTERED**, for the reason given above: 1,944
is `--phase-gate floor`'s value on this release exactly as 1,939 was its
value on the previous one and 1,946 on the signed one. **THE VERDICT STILL
DID NOT MOVE** -- 0 of 22 rows clear Bonferroni at alpha = 0.05/44, for the
third release running. The largest effect is now interactome/HCoV-229E at
`p_matched` 0.005170 against the 0.001136 threshold, where it was
interactome/NCBITaxon:11082 at 0.006210; that row now reads 0.005820. So
the conclusion does not depend on the five proteins either.

**TWO MECHANISMS ARE CONSISTENT WITH A DROP AND NEITHER IS ESTABLISHED.**
Collapsing a salt onto its parent merges their edges and cannot by itself
remove a protein from the set, since the protein endpoint survives the merge;
so the drop has to come from an edge that stopped being emitted rather than
from one that was merged. The candidates are the 980 edges that sat on
duplicate nodes and the 31 compounds that had been in both label classes.
Distinguishing them needs the old operand, which is the paragraph above.
Neither is asserted here.

## Why this is a re-registration, and what that costs

The other four pre-registrations in this directory open by saying they were
written **before** the code, "because a prediction made afterwards is a
description". This one cannot say that. The analysis ran on 2026-10-28 and
its figures are quoted in `docs/scope_and_limits.md` and in the paper.

**So the 2026-10-28 run was unregistered and its figures are exploratory.**
Writing registered predictions for a measurement already in hand would be
the exact failure this directory's format exists to prevent, and no
paragraph here retroactively fixes that. What was genuinely still ahead, when this
was drafted on 2026-11-02, was a **re-run** -- forced by two defects found
that day that change every p-value in the file, and whose parameters were
then undecided. This document registers that re-run. **It has since been
executed, on 2026-11-11, before this signature**; the header says what that
costs and what it does not.

The published null does not change status because of this. It was
exploratory before and it is exploratory now; what changes is that the
re-run will not be.

## The defects forcing the re-run

**One generator consumed in row order.** `matched_null` already documents
three bugs of this kind -- unsorted pools, unsorted bin order, degree ties
broken by set iteration -- and each was the same sentence: the seed stopped
determining the draw. All three were fixed and a single run-wide generator
survived them. Measured on the real function with `Random(0)` and identical
data, one taxon's null total:

    rows A, B        A = 1209
    rows C, A, B     A = 1165      a virus added BEFORE it
    rows B, A        A = 1165      the same two, other order

Adding a taxon to `config/viruses.yaml`, reordering the register, or running
one channel alone moves every downstream p-value. The curated layer will add
five viruses with dependency factors, so this is not hypothetical.

**Two resolvable steps below the threshold -- UNTIL 2026-11-11.** At 2,000
permutations an empirical p takes values k/2001, so below the corrected
0.001136 there were exactly two attainable values -- 0.000500 and
0.001000 -- and the third, 0.001499, is above it. The test could
distinguish "significant" from "not" and nothing finer. **The run now on
disk uses 100,000 and has 113 attainable values below alpha.**

Both are fixed in code, and the file on disk now comes from the fixed
behaviour: it records `seed_scheme` as `per-row: Random(f'{seed}|{channel}|
{virus}')` and `permutations: 100000`. The reporting path is kept for the
older artifacts -- `landscape.load` describes a file with NO `seed_scheme`
key as shared-generator in its protocol string, so a stale result cannot
present itself as this one.

## The measurement

Unchanged, and stated here because a re-run that quietly changes the
estimand is a different study. For each (channel, taxon): the overlap
between the drugged proteome and that taxon's factor set, against a
**degree-matched permutation null** that draws, from each degree bin, as
many proteins as the drugged set has in that bin.

    pathway channel       factors = CRISPR dependency factors
    interactome channel   factors = host proteins a viral protein touches

**The estimand is the overlap count**, compared against its matched
expectation. Both tails from the same null.

## Fixed in advance

- **Permutations: 100,000.** **DECISION**. The hard floor is one attainable
  value below alpha, which needs N+1 > 1/alpha = 880. Ten steps needs
  N >= 8,799. 100,000 gives 113 steps and costs 50x the current runtime.
  The proposal is 100,000 rather than 8,799 because the family size will
  change when the curated layer lands, and a count chosen to just clear
  today's alpha has to be revisited then.

- **Drug-target admission: `--phase-gate floor`.** **DECISION**, and it is
  here because 2026-11-10 proved this analysis's verdict depends on it. The
  ChEMBL drug-target query keeps compounds at `max_phase >= 2`; fix-list B1
  built a `relaxed` policy that also keeps any compound the graph already
  holds. **The drugged proteome -- this analysis's first operand -- is 1,944
  human proteins under `floor` and 2,489 under `relaxed` (the `relaxed`
  figure measured on the superseded release; see *A registered operand
  moved*), and the verdict follows it:** interactome/NCBITaxon:2697049 reads
  `p_matched` 0.2451 under `floor` on the current release -- it read 0.235 at
  signature and 0.2237 after L27 -- and 0.000500 under `relaxed` on the
  signed release, and the run as a whole goes from 0 of 22 rows clearing
  Bonferroni to 3. A registration that leaves this free
  fixes the permutation count to three digits while leaving the operand
  undefined, which is the same incompleteness the seeding scheme had.

  **AND THE TWO POLICIES ARE NOT EQUALLY ADMISSIBLE HERE, WHICH IS WHY THE
  REGISTERED ONE IS `floor`.** `relaxed` admits compounds *the graph already
  holds*, which means compounds someone screened against a virus. Their
  human targets therefore over-represent virus-interacting proteins by
  construction: of the 543 proteins `relaxed` newly drugs, **96.5% appear in
  a PPI against 85.4% of the 1,946 then already drugged, and 45.7% are dependency
  factors against 37.2%**. The degree-matched null controls for how well
  studied a protein is, not for whether its drug was chosen for antiviral
  screening.

  **THERE IS NO INTERNAL CONTROL AVAILABLE UNDER `relaxed`.** Every compound
  the policy can admit is one that was screened, so the change contains no
  unscreened comparison group. The only admission rule whose compounds are
  independent of this hypothesis is the one that drops the gate entirely --
  B1's option C, 676,997 new compounds -- which was rejected as too large to
  evaluate. B1 was then declined outright on a sharper measurement of the
  same fact: of the compounds `relaxed` gives a first host-side route, 84.6%
  carry positive antiviral labels against a 24.7% base rate for DENV-2. So a
  re-run under `relaxed` is not a sensitivity analysis of this one. It is a
  different hypothesis, on an operand selected on something close to the
  outcome, and it would have to be registered as such.

- **Correction family: every test the run performs, both channels, both
  directions.** On v0.1 that is 22 rows x 2 = 44, so alpha = 0.05/44 =
  0.001136. Computed by the script, which is the only thing that knows how
  many channels ran; `landscape.load` takes the declared count as a floor
  and may raise it, never lower it.

- **The six empty-rejection-region tests stay in the family.** **DECISION**.
  Six of the 22 downward tests cannot reject at any observation: the most
  extreme downward value is zero hits with probability exp(-m), and where
  that exceeds alpha -- m < 6.78 -- no observation clears the bar.

      interactome  JEV          0.83 expected
      interactome  West Nile    0.91
      pathway      West Nile    2.51
      interactome  DENV-3       4.43
      interactome  TBEV         4.62
      pathway      DENV-2       6.28

  Dropping them would give 0.05/38 = 0.001316 and is defensible: emptiness
  depends on the matched expectation, which the degree-matched design fixes
  before any outcome is looked at, so it is not choosing tests after seeing
  them. The proposal is to KEEP them, for two reasons. The threshold moves
  16%, from 0.001136 to 0.001316, and the smallest p in the family is
  0.005170 -- 3.9x the looser bar -- so nothing turns on it. And a family
  defined as "the tests that could have rejected" has to be recomputed
  whenever an expectation moves, which makes the correction depend on the
  release in a way a reader cannot check by counting.

- **Both procedures reported: Bonferroni and Benjamini-Hochberg at
  q = 0.05.** On the 2026-11-11 figures both reject nothing, and the
  smallest p is 0.005170 against a Bonferroni bar of 0.001136 and a BH
  critical value of 0.001136 at k=1. Bonferroni over 22 correlated rows is
  over-conservative and that objection is worth answering with a number
  rather than a sentence. Primary is Bonferroni; BH is reported beside it.

- **One generator per row, seeded `Random(f"{seed}|{channel}|{virus}")`,
  default seed 0.** A string seed because `Random` hashes a `str` through
  sha512 of its bytes rather than through `hash()`, so `PYTHONHASHSEED`
  cannot move it -- the property the three earlier fixes were about,
  asserted in subprocesses rather than reasoned about. The scheme is
  recorded in the output file, because a file written under either scheme
  is structurally identical and differs in every p-value.

- **Downward readability gate: `expected > -ln(alpha)`, derived from the
  run's own alpha.** It was a constant 3.0, from `exp(-3) = 0.0498 < 0.05`
  -- the uncorrected alpha -- while every verdict in the module is judged
  against alpha/n_tests. Derived per row, so a change in family size moves
  it automatically.

- **Primary direction: enrichment. Secondary: depletion.** Both are
  computed from the same null and both are in the correction family. The
  asymmetry is in interpretation only: the project's claim is about
  pathway-mediated routes existing, and depletion is reported because
  without it a low ratio cannot be told from an ordinary low draw.

- **Degree basis: `--degree-from study`**, the default, not the assembled
  graph. The graph basis is circular -- it derives degree from the same
  TARGETS and virus-host edges that define the sets being compared -- and
  biases every ratio downward.

- **No floor on taxon size for being tested.** A taxon with a small
  expectation can still return a significant excess, and excluding it would
  be a power decision disguised as a coverage decision. The minimum
  detectable effect is reported per row instead, in
  `docs/scope_and_limits.md`.

## Registered predictions

**P1. No test clears the corrected alpha, under Bonferroni or under BH.**
The null this document points at, stated so that failing to reject it is a
result rather than an absence.

**P2. Raising the permutation count will not produce a rejection.** The
smallest p in the 2026-11-17 family is 0.005170, which is 4.5x the
Bonferroni threshold. Resolution was never what blocked a finding, and
this prediction says so in advance: if a rejection appears at 100,000
permutations it did not come from finer resolution and something else
changed.

**AND THE SMALLEST p IS A ONE-SIDED ONE, which is worth saying because
it was misread while re-deriving these figures on 2026-11-17.** The
family holds both directions: `p_matched` for enrichment and
`p_matched_less` for depletion, and
`tests/test_landscape.py::test_the_smallest_p_and_its_multiple_of_alpha_are_right`
takes the minimum over both. 0.005170 is `NCBITaxon:11137`'s interactome
DEPLETION p. The smallest enrichment p is 0.005820
(`NCBITaxon:11082`), and reading only that column makes the family look
0.00065 further from rejection than it is. The test already derived the
right number; the attempted "correction" was the error.

**P2 WAS TESTED ON 2026-11-11 AND HELD.** The landscape ran at 100,000
permutations: 113 attainable values below alpha against two at 2,000, and
**zero rows clear the threshold in either direction**, unchanged. Every
p-value moved by less than 0.015 -- SARS-CoV-2's taxa row went 0.07096 to
0.05964 -- so the 2,000-permutation estimates were accurate and simply
could not support the alpha this document fixes. The negative result is
now powered rather than merely unrejected.

**P3. Per-row seeding will move individual p-values and not the verdict.**
Every row's draws change, so every p changes. The ordering of ratios and
the set of rejections should not. A verdict that flips on the seeding scheme
alone would mean the result was a draw rather than a measurement.

**P4. The six empty-rejection-region rows stay empty.** Emptiness is a
property of alpha and the expectation, not of N, so no permutation count
fixes it. If the curated layer later gives DENV-3, TBEV or JEV real
dependency-factor sets, their expectations rise and some of the six leave
that state -- which is the only thing that can.

## What would change the conclusion

A test clearing the corrected alpha under both procedures, on a row whose
minimum detectable effect is below the observed ratio. That would be the
first pathway-level finding in this project and would justify a per-pathway
test rather than a pooled overlap.

## What this cannot do

- It does not test pathway ENRICHMENT in the gene-set sense. The statistic
  is an overlap between two flat protein sets against a degree-matched
  background. `scripts/pathway_landscape.py` computes no per-pathway
  quantity, and the claim in the paper is being renamed accordingly (P1 in
  the audit fix list) rather than defended here.
- It does not establish mechanism for any row it cannot reject, and eight of
  22 rows need a 2.4x excess or more before they could. (Seven until
  2026-11-10, eight after that rebuild, seven again after 2026-11-11 --
  `interactome` DENV-1 crossed from 2.52x to 2.39x -- and eight again after
  2026-11-17, when the entry-assay rebuild moved that same row back to
  2.53x. Two rows sit within 0.05 of the line, so this count moves with the
  graph. It has been derived
  from the results file since 2026-11-10 rather than typed in two
  documents, which is the only reason the move is visible.)
- It does not settle `MIN_REACHED_POSITIVES` or `MIN_CLASS_SIZE`. Both are
  10; `MIN_REACHED_POSITIVES` justifies itself by matching
  `MIN_CLASS_SIZE`, and `MIN_CLASS_SIZE` at `scripts/hard_negatives.py:29`
  carries no comment at all. The pair is circular and neither is used by
  this analysis, so it is named here and left to its own registration.
- It does not license pooling the two channels. They ask different
  questions of the same virus -- what it NEEDS against what it TOUCHES --
  and are reported side by side.

## Implementation constraints

No scipy. The Poisson tail behind the minimum-detectable-effect table is
hand-written with an iterative pmf, because m = 425 overflows a factorial,
and is checked against a direct summation at sizes where the factorial
still fits.

100,000 permutations over 22 rows is 50x the current run. If that is not
affordable the fallback is 8,799, which clears ten steps at today's family
size and must be revisited when the family changes -- and the fallback has
to be chosen HERE, not after seeing how long the first attempt takes.


---

# Note, 2026-11-17: two cells' observed sets moved and the cause is not settled

**NOT AN AMENDMENT TO THE PROTOCOL.** Nothing registered changes. This
records an observation made while re-deriving the figures above, so it is
on the record before anyone reads the two numbers as noise.

Between the pre-entry-assay release and the one on disk, `n_observed`
moved in two interactome cells — `NCBITaxon:31631` 62 -> 63 and
`NCBITaxon:11084` 3 -> 4 — and every p_matched in the family shifted in
the fourth decimal. The drugged proteome growing 1,939 -> 1,944 (recorded
above) accounts for a shift in every p, because it is the observed set
the null is built around.

**IT DOES NOT OBVIOUSLY ACCOUNT FOR A NEW MEMBER OF AN INTERACTOME SET.**
`interactome_partners` excludes hub neighbours by `g.degree[nbr] >
hub_cut`, and `Graph.degree` counts **every** predicate including the
held-out labels — which is the same property that made the degree
comparator in `docs/prereg_structure_scorer.md` defective, found the same
day. `study_degree` was written precisely to keep that out of the
degree-MATCHING ("WHY THE DEFAULT DEGREE IS NOT SAFE HERE"), and the hub
cut still reads it.

**SETTLED 2026-11-17, AND THE ANSWER IS NO.** All 577 clinical edges
run `SmallMolecule -> Disease`; not one has a `Protein` endpoint, so no
protein's degree moves and the hub cut above cannot see them. The
observed set is the drugged proteome, derived from compound-to-protein
edges, which they also do not touch. **The clinical ingest cannot move
this analysis at all**, which makes the controlled re-run unnecessary
rather than merely unaffordable.

`scripts/diagnose_label_degree_weighting.py` reports the same thing from
the other side: the hub set is **62 nodes with the held-out degree in
and 62 with it out**, identical membership. Held-out degree reaches 22
Protein nodes and every one of them is VIRAL, while the hub cut is
class-scoped to host proteins acting as path connectors.

**SO THE TWO MOVED CELLS BELONG TO THE ENTRY-ASSAY INGEST, and one of
them is now named.** `NCBITaxon:11084`'s interactome set went 3 -> 4
because **`UniProtKB:Q92993`** -- one of the five proteins that joined
the drugged proteome in that ingest, named in the amendment above --
physically interacts with a TBEV protein. `NCBITaxon:31631`'s 62 -> 63
is **not** explained by any of the five: none of them partners a
SARS-CoV protein over `PHYSICALLY_INTERACTS_WITH`. That one remains
open, and it is open against the entry-assay ingest rather than against
the clinical one.

---

# Amendment, 2026-11-18: the one open cell is closed, and the note above looked in the wrong channel

**NOT AN AMENDMENT TO THE PROTOCOL.** Nothing registered changes, no
figure in `data/results` moves, and the verdict is untouched: still 0 of
22 rows clearing Bonferroni. This closes the question the note above
left open, and corrects two errors in how that note stated it.

Measured by `scripts/diagnose_landscape_count_move.py`, artifact
`data/results/landscape_count_move.json`, `make landscape-count-move`
(6.7s). Derived by `tests/test_landscape_count_move.py`.

## The note named the wrong channel, and the plan named the wrong virus

Both are errors in the record rather than in the analysis, and both
sent a reader looking where the answer could not be:

- The note calls `NCBITaxon:31631` 62 -> 63 an **interactome** cell. It
  is a **dependency** cell. OC43's interactome cell is 28 and did not
  move. So the note's search for a new `PHYSICALLY_INTERACTS_WITH`
  partner was conducted in a channel that cell never reads, and **no
  result that search returned could have been the right one** — which
  is why it came back with nothing and the cell stayed open.
- `docs/v1_release_plan.md`'s open-work table called `NCBITaxon:31631`
  **SARS-CoV**. `config/viruses.yaml` gives 31631 as **HCoV-OC43**;
  SARS-CoV is 694009, whose cells are 167 and 63 and neither moved.

The second error is the more instructive one. A taxon id is not
self-describing, and this repository has already paid for that four
times over in `virus_register.py`'s docstring. The test added here
resolves every moved cell's name through `config/viruses.yaml` rather
than accepting a name written in prose.

## The mechanism, which the note did not state

`overlap_row` computes

    n_observed = len(targets & factors)

where `targets` is the drugged proteome and `factors` is the cell's
dependency or interactome set. **A cell can therefore move without any
edge in its own channel changing**: it is enough that a protein already
inside `factors` joins `targets`. The note reasoned about new edges,
and both moves are instead intersection growth against a set that did
not move at all.

The drugged proteome went 1,939 -> 1,944, five joined, none left — and
`diagnose_landscape_count_move.py` re-derives those five **by set
difference over the two artifacts' recorded `members`**, obtaining
exactly the five this document names above. That is an independent
confirmation of a figure nobody could re-measure once the earlier
release was overwritten, and it is what the member list in `operands`
was added for.

## Both cells, with a cause each

    channel      taxon             virus       before  after   via
    -----------------------------------------------------------------------
    dependency   NCBITaxon:31631   HCoV-OC43       62     63   UniProtKB:Q16280
    interactome  NCBITaxon:11084   TBEV             3      4   UniProtKB:Q92993

Q92993 is the protein the note already named for TBEV. Its reasoning —
that it "physically interacts with a TBEV protein" — is true, and is
also not the mechanism: the interaction was there before, and what
changed was Q92993 becoming drugged.

## What makes this a measurement rather than a story that fits

Two numbers can be explained by almost anything. The script states both
falsifiers before it runs and reports them over **all 22 cells**:

    moves with no operand cause               none
    predicted moves that did not happen       none

The second is the one that does the work. The operand difference
predicts a nonzero delta in **exactly 2 of 22 cells**, and they are
exactly the two that moved. Twenty cells predict zero and moved zero.
Had the five proteins landed in a third cell's factor set, that cell
would have had to move too, and this would have failed.

An independent check on the same claim: if no channel edge changed,
every cell's `n_dependency` — the size of `factors` — is identical
across the two artifacts. **It is, in all 22**, including the 657 of
OC43's dependency set and the 35 of TBEV's interactome set.

## The controlled re-run is not needed

The note named a controlled `ncats` re-run as the way to close this,
and `docs/v1_release_plan.md` carried it as remaining work. **It is
unnecessary.** The earlier artifact is in git (`57602d9`), it records
the operand's members, and the factor sets are rebuilt from the release
on disk — so the comparison the re-run would have produced is already
available, at 6.7 seconds instead of a nine-minute layer rebuild on a
release that no longer exists to be compared against.

**SETTLED 2026-11-18. Both moved cells have a named cause, the cause is
the same in both, and no edge in either channel changed.**
