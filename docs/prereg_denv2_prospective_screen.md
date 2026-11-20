# Pre-registration: a prospective wet-lab screen, DENV-2

**Signed 2026-11-20, before any candidate list for this screen exists.**

Nothing in this document names a compound. The list is drawn by a single
recorded command **after** this is signed, and the command is fixed below.
A shortlist drawn first and judged afterwards is not a prospective test of
anything, and this project has the whole of `docs/evaluation.md` as a
record of what that costs.

---

## 1. What is being tested

Whether the structure scorer's ranking of **never-screened** compounds
against DENV-2 hits at a higher rate in the laboratory than compounds
drawn at random from the same pool.

Not whether the scorer is "good". Not whether any particular compound
works. One number: the difference in hit rate between two arms.

## 2. Why this is worth doing at all, and why it is DENV-2 alone

Every retrospective result in this repository is consistent with the
scorer retrieving analogues rather than predicting activity. The
measurement that says so most sharply is in the `claim_hit_rate` block of
`docs/evaluation.md`: on the one population matched by assay rather than
by source, the top of this ranking holds **no actives at all** once each
candidate's nearest known actives are withheld, while the same cell's AUC
clears its registered alpha.

DENV-2 is the one cell where that does not happen. Prospectively — trained
on activity published at or before 2021, tested on activity published
after — its precision@10 is **70% against a base rate of 6.51%**, a lift
of **10.7x**, and it holds under the registered k10 control. Full table:
`make hit-rate`.

**AND THE OTHER NINE CELLS ARE NOT REGISTRABLE, WHICH IS A MEASUREMENT
AND NOT A PREFERENCE.** `make screen-power` sizes a two-arm screen for
every prospective cell:

    virus            base   prec   lift     k @a=0.05    k @a=0.005
    DENV-1         46.4%    60%   1.3x  out of reach  out of reach
    DENV-2          6.5%    70%  10.7x            10            16
    DENV-3         47.0%    60%   1.3x  out of reach  out of reach
    DENV-4         43.8%    80%   1.8x            26            46
    West Nile      48.9%    20%   0.4x  out of reach  out of reach
    Yellow fever   37.0%    10%   0.3x  out of reach  out of reach
    HCoV-229E      70.4%   100%   1.4x            22            34
    SARS-CoV-2     56.7%    90%   1.6x            26            44
    HCoV-OC43      67.7%    55%   0.8x  out of reach  out of reach
    Zika           18.5%    10%   0.5x  out of reach  out of reach

Six of ten cannot reach 80% power at **any** arm size up to 120, for two
different reasons that the lift column separates:

- **No enrichment.** West Nile, Yellow fever, Zika and HCoV-OC43 rank at
  or below their own base rate. There is nothing to confirm.
- **No headroom.** DENV-1 and DENV-3 rank at 60% against a 47% base. That
  is a lift of 1.3, and separating 1.3 from 1.0 takes hundreds of
  compounds however high the AUC is. Their AUCs are 0.8025 and 0.7300 —
  among the highest in the study — which is the clearest statement
  available that **an AUC does not tell you whether a screen is worth
  running.**

Three more cells can be powered and are deliberately not registered here:
DENV-4 (46 per arm), SARS-CoV-2 (44) and HCoV-229E (34). Each needs more
than twice DENV-2's compounds to establish a lift under 2x, and
HCoV-229E's pool is 70.4% active, so its result would mostly say the
library was already good.

**THE FAMILY IS ONE CELL, AND REGISTERING MORE WOULD COST THIS ONE.** At
alpha 0.05 DENV-2 needs 10 compounds per arm; at the 0.05/10 a ten-cell
family would impose it needs 16. Registering cells that cannot be powered
does not merely waste compounds, it raises the bar on the cell that can.
So: **family size 1, alpha 0.05**, and this document claims nothing about
any other virus.

## 3. The design

Two arms, **20 compounds each, 40 in total**, from the same candidate pool,
plated together, read blind.

    arm A    the top 20 by max-Tanimoto to known DENV-2 actives
    arm B    20 drawn uniformly at random from the same pool

**WHY A RANDOM ARM AND NOT A THRESHOLD.** A one-arm screen must be
compared against an assumed base rate, and the base rate of an unscreened
pool is exactly the thing nobody knows —
`scripts/hit_rate_at_k.py` is explicit that its own base rates are an
upper bound, because they come from compounds somebody chose to screen.
With a random arm on the same plates, the comparator is measured in the
same experiment, under the same assay, cell line, passage and batch. It
costs twice the compounds and it is the only version whose answer does not
rest on a number supplied from outside.

**Blinding.** Arm membership is not disclosed to whoever runs or reads the
assay. The plate map is generated with the two arms interleaved and the
mapping held by someone who does not read the plates.

## 4. The candidate pool, fixed now

A compound is eligible if and only if, on release `v0.1`
(`edges_sha 873623c9a9df2049`):

1. it is a `SmallMolecule` with a parseable SMILES and an ECFP4
   fingerprint (radius 2, 2048 bits);
2. it carries **no** `HAS_ANTIVIRAL_ACTIVITY_AGAINST` and **no**
   `MEASURED_INACTIVE_AGAINST` edge to `NCBITaxon:11060` — never screened
   against DENV-2 in either direction;
3. it is not the same structure as any known DENV-2 active, by
   `kgav.normalize_chem.same_structure` (fix-list L27). A compound that
   differs from a known active only in salt form or tautomer is not a
   prediction;
4. it is orderable — see below.

**ORDERABILITY IS A FILTER ON THE POOL, NOT ON ARM A.** It must be applied
to the whole eligible set *before* either arm is drawn, because applying
it after ranking would let arm A be "the top compounds we could buy" while
arm B stayed "anything at all", and the difference between the arms would
then include a purchasability difference. Whoever determines orderability
does so **without seeing the ranking**.

## 5. The statistic

One-sided Fisher exact test on the 2x2 of hits against misses, arm A
versus arm B. Alpha **0.05**, family size **1**, as justified in section 2.

**Power, derived by `make screen-power`** at the measured
precision (70%) and base rate (6.51%):

    k per arm    total    power at a=0.05    power at a=0.005
           10       20              88.7%               48.4%
           16       32              98.6%               88.6%
           20       40              99.8%               96.4%

20 per arm is chosen over the 10 that alpha 0.05 requires, for one stated
reason: at 20 the screen still has **96.4%** power at the alpha a ten-cell
family would impose, so if this design is later extended to other viruses
this cell's result survives the correction rather than having to be rerun.

At 20 per arm, with arm B hitting at its base rate, the expected counts
are **14.0 hits in arm A and 1.30 in arm B**. The test rejects at
alpha 0.05 with **5 or more** hits in arm A against 0 in arm B.

## 6. A hit, defined before anything is plated

A compound is a hit if it reduces DENV-2 viral signal with an
**EC50 at or below 10 micromolar** in a cell-based assay, measured in
concentration-response over at least five concentrations, with the
plate's own positive and negative controls passing.

**Cytotoxicity is recorded and is not part of the primary endpoint.** CC50
is measured on the same compounds, and a selectivity index is reported
beside every hit. A compound that reduces signal by killing cells is a hit
by this definition and the registration says so now, rather than being
reclassified afterwards. The secondary endpoint below is where
selectivity enters.

**Secondary, declared now and not substituted for the primary:** the same
test on hits with a selectivity index of 10 or greater. If the primary
clears and the secondary does not, the finding is "the ranking enriches
for antiviral signal, and that signal is substantially cytotoxic", which
is a different and weaker claim. It will be reported in those words.

## 7. Predictions

- **P1.** Arm A hits at **70% or above**. This is the conservative figure:
  the k10-controlled retrospective precision. The *uncontrolled*
  retrospective precision@10 for this cell is 100%, so P1 is deliberately
  set below what the retrospective data suggests.
- **P2.** Arm B hits at **between 2% and 15%**, bracketing the measured
  6.51%. Arm B is the one quantity this screen measures that nothing in
  the repository can predict, so the interval is wide on purpose.
- **P3.** Arm A's mean max-Tanimoto to the known actives exceeds arm B's.
  A ranking that did not separate the arms on the statistic it is built
  from would mean the draw was wrong, not that the scorer failed.
- **P4.** At least 15 of arm A's 20 compounds carry **no** graph edge
  other than the ones used to rank them. The shortlist's ranking uses no
  graph information, and this records that the screen is not quietly a
  test of graph connectivity.

## 8. Refusal conditions — the screen reports nothing if any holds

Stated now so none of them can be weighed after the numbers exist.

1. **Fewer than 20 eligible orderable compounds** in either arm. A short
   arm changes the test's power, and a power computed after the fact is
   not a power.
2. **Any plate control fails.** No partial read, no dropping a plate and
   keeping the rest.
3. **Arm B's hit rate falls outside P2's 2%–15%.** Then the pool that was
   screened is not the pool that was ranked, and the comparison is
   between two different populations. This is a refusal and not a result
   in either direction — including if arm B hits at 0%, which would make
   the pool's base rate unmeasurable from this experiment.
4. **The release moves.** Both arms are drawn from `edges_sha
   873623c9a9df2049`. If the release is rebuilt before the draw, this
   registration is void and must be re-signed against the new one.
5. **Arm membership is disclosed** to anyone who runs or reads the assay
   before the plates are read.

## 9. The decision rule — exactly one of these is announced

- **CONFIRMED.** Fisher exact p <= 0.05, arm A above arm B, no refusal
  condition. Announced as: *the structure scorer's ranking enriched for
  DENV-2 activity in compounds never screened against it, at a measured
  rate, prospectively.* This is the strongest claim this project would
  have. It licenses nothing about any other virus, and section 2 is the
  reason.
- **NOT CONFIRMED.** p > 0.05 with both arms readable. Announced as: *the
  ranking did not enrich for activity in a screen with 99.8% power to
  detect the retrospective effect.* At this power, this outcome is
  evidence the retrospective 10.7x was itself analogue retrieval, and
  `docs/evaluation.md`'s claim block is amended to say so.
- **INVERTED.** Arm B above arm A at p <= 0.05. Announced in those words.
  Nothing here licenses inverting the ranking as a strategy; see the
  prospective cells that clear below chance.
- **REFUSED.** Any condition in section 8. The reason is published with
  the refusal and no hit rate is reported.

## 10. What a CONFIRMED result would and would not license

**Would:** that for DENV-2, on this release, ranking never-screened
compounds by chemical similarity to known actives enriches for activity
at a measured rate.

**Would not:** anything about the other fifteen viruses — six of their
cells cannot be powered at all, and two rank *below* their own base rate.
Anything about the knowledge graph: the ranking uses no graph edges, and
`make auc-contribution` measures the graph's contribution to every scorer
in this project as zero or negative against its own degree baseline.
Anything about clinical use: the endpoint is an in-vitro EC50 in one
assay, and the one independent label axis available here anti-correlates
with this statistic (0.4087 against a degree baseline of 0.7751).

## 11. How the list is drawn, after this is signed

One command, recorded with its output hash:

    python scripts/repurposing_candidates.py \
        --virus NCBITaxon:11060 \
        --exclude-labelled \
        --top 20 \
        --random-arm 20 --seed <recorded at draw time> \
        --out data/results/denv2_screen_arms.json

`--exclude-labelled`, `--random-arm` and `--seed` **do not exist yet**.
They are named here so the command cannot be adjusted after the list is
seen, and implementing them is the next commit. The seed is recorded in
the artifact and in the commit message at draw time.

---

## Verdict

*Unsigned. Nothing has been drawn and nothing has been screened.*

This section is where the outcome of section 9 is written, once, with the
counts. If this document is read with this line still present, the screen
has not been run.
