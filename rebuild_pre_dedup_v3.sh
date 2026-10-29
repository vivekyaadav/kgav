#!/usr/bin/env bash
# Rebuild the PRE-DEDUP release in a throwaway worktree.
#
# v1 hardcoded a commit SHA from Claude's clone, which `git am` had
# rewritten, so it did not exist here.
#
# v2 found the commit correctly and then broke on the symlink:
#   ln -s SRC DIR, when DIR already exists, creates DIR/<basename SRC>.
# The worktree checkout creates a real data/raw (.gitignore keeps
# data/**/.gitkeep tracked), so the link landed at data/raw/raw and
# every phase-1 layer saw an empty data/raw. The build then REFUSED to
# start, correctly, rather than assembling a graph from nothing.
#
# v3 removes the directory first and VERIFIES the inputs are visible
# before running anything.
#
# Nothing here writes to ~/projects/kgav/data: data/raw is a read-only
# symlink, everything else is the worktree's own.
set -u
MAIN="$HOME/projects/kgav"
WT="$HOME/kgav-pre-dedup"
NEEDLE='sources.update(p for p in str(src).split'

cd "$MAIN" || exit 1

FIX=$(git log --format=%H -S"$NEEDLE" -- src/kgav/assemble.py | tail -1)
if [ -z "${FIX:-}" ]; then
    echo "FAILED: no commit here introduces the dedup line into"
    echo "src/kgav/assemble.py. Searched for: $NEEDLE"
    exit 1
fi
echo "dedup fix commit : $FIX"
git log -1 --format='                   %s' "$FIX"

HAVE=$(git show "$FIX:src/kgav/assemble.py" | grep -c "sources.update" || true)
PREV=$(git show "$FIX^:src/kgav/assemble.py" | grep -c "sources.update" || true)
echo "occurrences      : $HAVE at the commit, $PREV at its parent"
if [ "$HAVE" -lt 1 ] || [ "$PREV" -ne 0 ]; then
    echo "FAILED: that commit does not cleanly introduce the fix."
    exit 1
fi
PARENT=$(git rev-parse "$FIX^") || exit 1
echo "pre-dedup parent : $PARENT"

if [ -e "$WT" ]; then
    echo "removing stale worktree at $WT"
    git worktree remove --force "$WT" 2>/dev/null || rm -rf "$WT"
    git worktree prune
fi
git worktree add "$WT" "$PARENT" || exit 1

# THE DIRECTORY MUST GO FIRST. See the v2 note above.
rm -rf "$WT/data/raw"
ln -s "$MAIN/data/raw" "$WT/data/raw" || exit 1
mkdir -p "$WT/data/interim" "$WT/data/releases"

# AND THE LINK MUST BE CHECKED, not assumed. v2's failure mode was a
# link that existed, pointed somewhere real, and exposed nothing.
echo
echo "=== raw inputs visible through the link ==="
MISSING=0
for want in proteomes host orcs chembl ncats; do
    if [ -e "$WT/data/raw/$want" ]; then
        echo "  ok       data/raw/$want"
    else
        echo "  MISSING  data/raw/$want"
        MISSING=$((MISSING + 1))
    fi
done
if [ "$MISSING" -gt 0 ]; then
    echo
    echo "FAILED: the symlink does not expose the raw inputs."
    echo "What is actually at $WT/data/raw:"
    ls -la "$WT/data/raw" | head
    echo "And at $MAIN/data/raw:"
    ls "$MAIN/data/raw" | head
    exit 1
fi

cd "$WT" || exit 1
echo
echo "=== plan (expect 10 layers to RUN, 2 assemble passes) ==="
python scripts/rebuild_release.py --plan || exit 1
echo
echo "If any layer above says REQUIRED or absent for a reason other"
echo "than curated/gtop, stop and send this output rather than"
echo "building: the 'before' release would not match the one that was"
echo "overwritten, and comparing against it would be worse than not"
echo "comparing at all."
echo
read -r -p "press enter to build (~16 min), or ctrl-c to stop: " _
time python scripts/rebuild_release.py 2>&1 | tail -15
echo
echo "Then:  cd $MAIN && python compare_releases.py \\"
echo "         $WT/data/releases/v0.1 data/releases/v0.1"
