#!/usr/bin/env python3
"""What building each missing comparator virus would buy M8, and what it costs.

THE QUESTION THIS ANSWERS. The fold layer finds 36 allowed transfers; only
2 have both endpoints in the spine, so M8 can walk 2. Unblocking the rest
means the spine building comparator chains it does not build today, and
"which comparator next" is a scope decision that should be made against
numbers rather than against whichever virus name came to mind.

MARGINAL, NOT BLAME-COUNTED. ingest_folds' blocking_chains answers "how
many blocked transfers does this chain appear in", which double-counts a
pair whose BOTH comparators are unbuilt: HCV NS5B -> Poliovirus 3D is
charged to each. Building one of those two unblocks nothing. So this
script reports the marginal gain -- transfers that become walkable if that
virus ALONE is added -- and the pairs needing two viruses separately.

Run after scripts/ingest_folds.py has written its cache.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import yaml

from ingest_folds import (SCORE_ELIGIBLE_ROLES, allowed_pairs, load_cache,  # noqa: E402
                          missing_inputs, walkable_transfers)
from kgav.folds import gate  # noqa: E402


def register_state(cfg: dict) -> dict[str, tuple[str, str]]:
    """taxon -> (state, detail). Three states the config already encodes."""
    out: dict[str, tuple[str, str]] = {}
    for v in cfg["viruses"]:
        out[str(v["taxon"])] = ("BUILT", v.get("proteome_id", "-"))
    for s in cfg.get("surveyed") or []:
        tax = str(s["taxon"])
        if tax in out:
            continue
        has = s.get("proteome_id") and s.get("chain_source")
        out[tax] = ("SURVEYED" if has else "SURVEYED (no chain_source)",
                    s.get("proteome_id") or "-")
    return out


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path,
                    default=root / "data/raw/folds/phase2_chains.json")
    ap.add_argument("--spine", type=Path,
                    default=root / "data/releases/v0.1-spine")
    ap.add_argument("--config", type=Path, default=root / "config/viruses.yaml")
    args = ap.parse_args()

    if not args.cache.exists():
        print(f"no chain cache at {args.cache}; run scripts/ingest_folds.py",
              file=sys.stderr)
        return 2
    raw = json.loads(args.cache.read_text())
    missing = missing_inputs(raw)
    if missing:
        print(f"{args.cache} predates {', '.join(missing)}; re-run "
              f"scripts/ingest_folds.py without --offline", file=sys.stderr)
        return 2
    chains = load_cache(raw)
    spine_ids = {json.loads(line)["id"]
                 for line in (args.spine / "nodes.jsonl").read_text().splitlines()
                 if line.strip()}
    cfg = yaml.safe_load(args.config.read_text())
    state = register_state(cfg)

    eligible = [c for c in chains
                if c.ok and c.spec.role in SCORE_ELIGIBLE_ROLES]
    results = [(a, b, gate(a, b))
               for i, a in enumerate(eligible) for b in eligible[i + 1:]
               if a.spec.family != b.spec.family]
    walk, blocked = walkable_transfers(results, spine_ids)

    print(f"{len(allowed_pairs(results))} allowed transfers | "
          f"{len(walk)} walkable today | {len(blocked)} blocked\n")

    # taxon of each unbuilt chain, and the pairs it blocks
    def unbuilt(c):
        return c.node_id not in spine_ids

    solo: dict[str, list] = defaultdict(list)   # one missing end
    duo: dict[tuple, list] = defaultdict(list)  # two missing ends
    for a, b, g in blocked:
        gaps = sorted({c.spec.virus_taxon for c in (a, b) if unbuilt(c)})
        (solo[gaps[0]] if len(gaps) == 1 else duo[tuple(gaps)]).append((a, b, g))

    def tax(t):
        return t.removeprefix("NCBITaxon:")

    def label(t):
        st, pid = state.get(tax(t), ("NOT IN REGISTER", "-"))
        return f"{st:<27} {pid}"

    names = {}
    for c in chains:
        names.setdefault(c.spec.virus_taxon, c.spec.label.split()[0])

    print("building ONE virus unblocks (marginal — the other end is already "
          "in the spine):")
    rows = sorted(solo.items(), key=lambda kv: -len(kv[1]))
    for t, pairs in rows:
        s = Counter(g.strength for _a, _b, g in pairs)
        print(f"  {names.get(t, tax(t)):<14} {tax(t):<8} "
              f"+{len(pairs):>2} transfer(s)  "
              f"({', '.join(f'{k} {v}' for k, v in sorted(s.items()))})"
              f"   {label(t)}")
    if not rows:
        print("  (none — every blocked transfer is missing BOTH ends)")

    if duo:
        n = sum(len(v) for v in duo.values())
        print(f"\nneeding TWO viruses before a single transfer appears "
              f"({n} transfer(s)):")
        for (t1, t2), pairs in sorted(duo.items(), key=lambda kv: -len(kv[1])):
            print(f"  {names.get(t1, tax(t1))} + {names.get(t2, tax(t2))}"
                  f"   {len(pairs)} transfer(s)")

    print("\nthe same viruses, by what the register already knows about them:")
    allgaps = set(solo) | {t for pair in duo for t in pair}
    by_state: dict[str, list[str]] = defaultdict(list)
    for t in allgaps:
        by_state[state.get(tax(t), ("NOT IN REGISTER", "-"))[0]].append(
            names.get(t, tax(t)))
    for st in sorted(by_state):
        print(f"  {st:<27} {', '.join(sorted(by_state[st]))}")

    print("\nA transfer is only worth unblocking if its strength is worth "
          "walking.\nCount the strong ones: those are the conservation "
          "hypothesis's own case.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
