# 08 — What this graph cannot do

Every limit here has a figure behind it, and the figure is read out of an
artifact rather than typed. A limit restated with its own copy of the numbers
is a second thing to go stale — which `docs/scope_and_limits.md` has had
happen to it twice.

## One strain stands for each flavivirus

The register names, per virus, the single accession whose mature-peptide
chains the graph carries. Seven of nine surveyed flaviviruses are curated
under a **strain** taxon rather than the species, so a lineage-walking lookup
would have "worked" for all of them and left no record of which strain it
picked. Declaring it makes the choice auditable; it does not make one strain
representative of a serotype.

## Most drug–target edges do not say which direction

16 predicates carry edges and `TARGETS` is the
host-directed channel's entry point. On this release 10,177 of 17,635
`TARGETS` edges state no direction — 57.8%, down from 84.5% before the Guide
to PHARMACOLOGY layer. "This drug targets a dependency factor" therefore
cannot always become "inhibiting it should block infection", and M2's
`direction = dependency` constraint silently excludes the undirected ones.

## The analysable set is a fraction of the molecules

78,649 nodes include 25,243 small molecules, and the
host-directed channels can only reach the ones carrying a human target. The
cold-start similarity layer exists because 953 of 1,846 SARS-CoV-2 positives
have **no metapath at all** and 947 of those have no target edges whatsoever —
they carry a measured activity and nothing else.

## An AUC above chance is not an expected hit rate

AUC is the probability that a randomly chosen active outranks a randomly
chosen inactive: a statement about the **whole ordering**. A screen never sees
the whole ordering; it sees the top of a list. The two can disagree
completely, and here they do. On the one population matched by assay rather
than by source, the head of the ranking holds **no actives at all** once each
candidate's nearest known actives are withheld, under either registered
control — while the cell's AUC clears its alpha either way.

So a row on a shortlist is the claim that a compound **resembles** a known
active, with the known active named. `make hit-rate` measures the head at
three exclusion strengths and labels itself post-hoc.

## A Tanimoto is not readable until you know what a failure scores

The median resemblance a **screened-and-inactive** compound of the same virus
reaches against its known actives spans a factor of three across the cells of
this release. The same 0.52 is therefore near the top of one cell's failures
and under the median of another's, so a flat display threshold is wrong in
both directions at once. Each row of the read-out carries its position in that
cell's own failed-resemblance distribution; the floor is the 50% line on it.

## A similarity score does not know what a compound is for

The ranking is max-Tanimoto and cannot tell a therapy from a contrast agent.
42 compounds in the pool are diagnostic, contrast or
radiotracer agents by ATC and are dropped; 205 antineoplastics
and 219 withdrawn compounds are labelled and kept, because the
confound is real and so is the class. **2,238 of
25,243 molecules resolve to no ChEMBL record at all** by id or by
structure, so their role could not be checked — and unchecked is not clean.

## The graph's own channels are at chance

0 of 8 metapath channels has a 95%
interval excluding 0.5 on the plate-matched population
(none). The ensemble clears, and
0.5578 of its 0.5666 is its coverage asymmetry
rather than its ranking. A channel reaching fewer than
10 measured actives is reported as **not evaluable**
rather than with a number, because its AUC describes the imputed floor.

## Prospectively, two cells are measured getting it backwards

2 of 10 cells clear **below** chance
prospectively: West Nile, Zika. That is not a null — it is the scorer,
trained on what was known and asked to rank what came next, putting the
actives last. Both cells clear **above** chance cross-sectionally, so a
cross-sectional verdict alone does not license a shortlist, and the read-out
refuses them by name.

## No prospective cell has a usable comparator

0 of 10 cells have a degree
baseline reaching at least 20% of the compounds
being ranked — the pre-cutoff training graph reaches 0.3%–6.5% of them. The
registered decision rule's second clause is therefore satisfied by nothing,
and the outcome is reported as *measurable, not comparable*, which is the
outcome the registration said to expect.

## The clinical axis contradicts the chemistry

0.4087 [0.3734, 0.4440], clearing
**below** chance, against degree at
0.7751. Whether a compound reached a trial is predicted by
how well studied it is. This is the only independent label axis in the
repository.

## The pathway route is a null, with its power stated

44 tests at 100,000 degree-matched permutations, and
**8 viruses have no cell at all**: DENV-1, DENV-3, DENV-4, HCoV-HKU1, JEV, SARS-CoV, TBEV, Yellow fever. A null
result is only a finding if the reader knows what would have survived, so the
minimum detectable effect is published per cell. Failing to detect a 2×
excess is not evidence against a 2× excess when the test needed 7×.

## The labels are one kind of measurement at one threshold

Until the indications layer landed, every AUC in the repository rested on
`HAS_ANTIVIRAL_ACTIVITY_AGAINST` alone — one predicate derived from in-vitro
concentration-response at one cut. A scorer that recovered the idiosyncrasies
of that screening and nothing else would pass all of them. There are now
271 `IN_TRIAL_FOR` and 1 `TREATS` edges, and
they go the other way.

Within the in-vitro axis: 4,839 label edges carry a potency ChEMBL flagged as
read **outside** its tested dose series (23 of them positives), 157 verdicts
rest on a relation nobody stated, 169 edges carry potency fields that
disagree, and 56.5% of positives at 10 µM are not positive at 1 µM.

## The circularity confound is measured, not removed

779 of 943 compounds carrying both an `INHIBITS` edge
and an antiviral label cite the same publication on both.
58 pairs cannot be withheld at all, because one side
carries no PMID. `--publication-disjoint` removes what it can reach; those
58 it cannot.

## Nine nodes hold two compounds' screening results

The identity policy collapses salts, charges and tautomers, which is the right
default and is not free. A Tanimoto of 1.0 means the same ECFP4 bits, **not**
the same compound: radius 2 folded into 2,048 bits cannot
see a saturated ring's size, a peptide's length, or a residue's position, and
rows at 1.0 against a different skeleton are marked with their heavy-atom
difference rather than read as identities.

## What the agent will not do, as a consequence

- refuse to rank **novel** chemical matter for any virus — compounds enter
  screening because they are novel, which is what a prior-evidence graph
  cannot represent;
- refuse any ranking question for a virus with no route at all, because an
  empty ranking reads as "no compound looks promising" and the difference
  between finding nothing and being unable to look is the whole distinction;
- attach a cell-context caveat to every host-directed path, because
  chloroquine's route through SIGMAR1 is real, well replicated, and led the
  field astray in 2020.
