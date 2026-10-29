#!/usr/bin/env bash
# Is assemble_release.py deterministic?
#
# Rebuilding the pre-dedup code produced edges_sha 2f51efb4f897538e
# where the release it was meant to reproduce had a9a68fe6a41380c6 --
# same node set, same 389,306 edges, different content. And the count
# of edges carrying a repeated curator came out 568 where the original
# release had 542.
#
# Same code, same inputs, different graph. If that is assembly rather
# than the layers, then the before/after comparison run today
# attributes nothing: two runs of identical code differ by themselves,
# so no difference between them can be credited to a code change.
#
# This re-runs ONLY the assemble step, on interim data that does not
# change, four times. 11s each.
#
#   A, B  default environment
#   C, D  PYTHONHASHSEED=0
#
# A != B          -> assembly is nondeterministic
# A != B, C == D  -> and it is string-hash randomisation, which is the
#                    usual cause: a set or dict iterated in a
#                    different order per process, deciding which of
#                    two witnesses is merged into which.
# A == B          -> assembly is stable and the difference came from a
#                    phase-1 or phase-3 LAYER instead; rerun this
#                    against those.
set -u
WT="$HOME/kgav-pre-dedup"
cd "$WT" || { echo "no worktree at $WT"; exit 1; }
[ -f data/interim/.keep ] || true

fp () {
  python - <<'PY'
import sys
sys.path.insert(0, "src")
from pathlib import Path
from kgav.provenance import release_fingerprint
f = release_fingerprint(Path("data/releases/v0.1"))
n = sum(1 for l in Path("data/releases/v0.1/edges.jsonl")
        .read_text().splitlines() if l.strip())
import json
dup = 0
for l in Path("data/releases/v0.1/edges.jsonl").read_text().splitlines():
    if l.strip():
        s = str(json.loads(l).get("primary_knowledge_source") or "")
        p = [q for q in s.split("|") if q]
        if len(p) != len(set(p)):
            dup += 1
print(f"{f['edges_sha']}  edges {n:,}  repeated-curator {dup:,}")
PY
}

run () {
  echo "--- run $1 ${2:-(default env)}"
  if [ -n "${2:-}" ]; then
    env PYTHONHASHSEED=0 python scripts/assemble_release.py >/dev/null 2>&1 \
      || { echo "    assemble FAILED"; return 1; }
  else
    python scripts/assemble_release.py >/dev/null 2>&1 \
      || { echo "    assemble FAILED"; return 1; }
  fi
  printf '    '; fp
}

echo "=== assemble_release.py, same interim data, four runs ==="
run A
run B
run C PYTHONHASHSEED=0
run D PYTHONHASHSEED=0
echo
echo "A and B equal?  -> assembly is deterministic, look at the layers."
echo "A and B differ, C and D equal -> string-hash randomisation."
echo "C and D differ too -> something else varies per process."
