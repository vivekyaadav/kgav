#!/usr/bin/env bash
# EXTRACTED VERBATIM from docs/evaluation.md '## Reproducing'.
# Not set -e: one step is REQUIRED to exit 1 (the overlap-policy
# control), and aborting there would lose every step after it.
set -u
cd ~/projects/kgav

# ---- after any rebuild of data/releases/v0.1, everything below is required.
# Each block regenerates files that `make check-results` audits; a block
# skipped leaves a file describing a graph that no longer exists.

# 1. the temporal split, which is also an input to step 3
python scripts/temporal_split.py --cutoff 2021 --undated exclude
python scripts/temporal_split.py --cutoff 2021 --undated computed_only
python scripts/temporal_split.py --cutoff 2021 --undated include   # LAST: the
                                                                   # split the
                                                                   # rest reads

# control: must refuse, naming the 19 already-known compounds. Writes a
# refusal record and exits 1 -- expected, and it writes no split, so it
# cannot disturb the one above.
python scripts/temporal_split.py --cutoff 2021 --undated include \
    --overlap-policy keep || true

# 2. per-virus baselines -- one file each, all SIXTEEN taxa in
#    config/viruses.yaml. It was eight until 2026-10-23, when the remaining
#    flaviviruses were promoted. Six of the eight added here have no
#    HOST_FACTOR_FOR edge, so their host-mediated channels are empty by
#    construction; the file is still written, and an empty channel recorded
#    is what lets a later ingest show up as a change.
for t in 2697049 694009 1335626 11137 277944 31631 290028 64320 \
         11060 11053 11069 11070 11082 11072 11089 11084; do
    python scripts/run_baselines.py --virus "NCBITaxon:$t"
done

# 3. metapath comparison, and the hard-negative protocols
python scripts/compare_metapaths.py

# rows 1-4, in order
python scripts/hard_negatives.py --selectivity selective-only
python scripts/hard_negatives.py --selectivity all
python scripts/hard_negatives.py --selectivity verified-only --selectivity-applies both
python scripts/hard_negatives.py --selectivity verified-only \
    --selectivity-applies both --publication-disjoint

# no-op control: must reproduce row 2 exactly
python scripts/hard_negatives.py --selectivity all --selectivity-applies both

# row 5 -- the current protocol, run LAST so its numbers are the ones left
# on disk. Running it first, as this block used to, meant the no-op control
# overwrote it and hard_negatives.json held row 2 while being read as row 5.
python scripts/hard_negatives.py --selectivity all --publication-disjoint \
    --label-source infores:ncats-opendata
python scripts/calibrate.py          # which channels were evaluable, and why not

# 4. the pathway landscape. It depends on nothing above -- it reads the
#    release directly and scores no compounds -- but check-results globs
#    data/results and audits it like everything else, so leaving it out
#    meant a reader could run this block top to bottom and still be told
#    a file was stale. It was missing here from the day the file was added
#    until 2026-10-28; it looked fine only because it was being rebuilt by
#    hand. tests/test_reproduce_block.py now fails on a results file that
#    no step here regenerates.
make pathway-landscape

# 5. every results file, AND the derived release they depend on, must
#    describe the graph on disk
make check-results
