#!/usr/bin/env python3
"""Day 6 driver: host dependency/restriction factors from BioGRID ORCS.

Downloads per-screen gene results (cached to disk, so re-runs cost nothing),
determines polarity, and emits HOST_FACTOR_FOR edges.

    export ORCS_KEY=...
    python scripts/ingest_orcs.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.aliases import AliasIndex
from kgav.emit import Emit
from kgav.orcs import (
    build_gene_to_proteins,
    ingest_screens,
    is_virus_resistance_screen,
    polarity_mode,
    virus_taxon,
)
from kgav.schema import load_schema

API = "https://orcsws.thebiogrid.org"


def fetch_screen(sid: str, key: str, cache: Path, retries: int = 3) -> list[dict] | None:
    f = cache / f"screen_{sid}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            f.unlink()
    url = f"{API}/screen/{sid}/?accesskey={key}&format=json"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                body = r.read()
            data = json.loads(body)
            if isinstance(data, dict) and data.get("STATUS") == "ERROR":
                print(f"    screen {sid}: {data.get('MESSAGE','error')[:60]}")
                return None
            f.write_bytes(body)
            return data
        except (urllib.error.URLError, json.JSONDecodeError, TimeoutError) as exc:
            if attempt == retries - 1:
                print(f"    screen {sid}: {exc}")
                return None
            time.sleep(2 ** attempt)
    return None


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--orcs-dir", type=Path, default=root / "data" / "raw" / "orcs")
    ap.add_argument("--host", type=Path, default=root / "data" / "releases" / "v0.1-host")
    ap.add_argument("--spine", type=Path, default=root / "data" / "releases" / "v0.1-spine")
    ap.add_argument("--out", type=Path, default=root / "data" / "releases" / "v0.1-orcs")
    ap.add_argument("--overrides", type=Path,
                    default=root / "config" / "orcs_polarity_overrides.json")
    args = ap.parse_args()

    key = os.environ.get("ORCS_KEY", "")
    if len(key) != 32:
        raise SystemExit("ORCS_KEY must be a 32-character key; export it first")

    screens = json.loads((args.orcs_dir / "screens_all.json").read_text())
    candidates = [s for s in screens if virus_taxon(s)]
    usable = [s for s in candidates if is_virus_resistance_screen(s)]
    print(f"{len(screens):,} human screens | {len(candidates)} coronavirus | "
          f"{len(usable)} virus-resistance")

    modes = {"screen_level": 0, "signed": 0, "unresolvable": 0}
    for s in usable:
        modes[polarity_mode(s)[0]] += 1
    print(f"  polarity: {modes['screen_level']} screen-level, {modes['signed']} signed, "
          f"{modes['unresolvable']} unresolvable")

    cache = args.orcs_dir / "screens"
    cache.mkdir(exist_ok=True)
    fetchable = [s for s in usable if polarity_mode(s)[0] != "unresolvable"]
    print(f"\nfetching gene results for {len(fetchable)} screens")
    gene_rows: dict[str, list[dict]] = {}
    for i, s in enumerate(fetchable, 1):
        sid = str(s["SCREEN_ID"])
        rows = fetch_screen(sid, key, cache)
        if rows:
            gene_rows[sid] = rows
        if i % 10 == 0:
            print(f"  {i}/{len(fetchable)}")
    print(f"  {len(gene_rows)} screens with gene data")

    index = AliasIndex.from_releases(args.spine, args.host)
    gene_to_protein = build_gene_to_proteins(
        [args.host / "edges.jsonl", args.spine / "edges.jsonl"])
    print(f"\n{len(gene_to_protein):,} genes map to host proteins")

    overrides = json.loads(args.overrides.read_text()) if args.overrides.exists() else {}
    em = Emit()
    s = ingest_screens(em, usable, gene_rows, gene_to_protein,
                       "infores:biogrid-orcs", overrides)

    print(f"\n{s['edges']:,} edges from {s['screens_used']} screens")
    print(f"  {s['dependency']:,} dependency, {s['restriction']:,} restriction")
    print("\n  skipped:")
    for k in ("not_a_resistance_screen", "polarity_unresolvable", "no_gene_data",
              "no_hits", "uncalibrated_sign", "gene_not_in_graph",
              "no_direction", "duplicate", "override_applied",
              "positive_is_dependency", "positive_is_restriction",
              "rationale_contradicted_by_anchors"):
        if s[k]:
            print(f"    {k:<28} {s[k]:>7,}")
    if em.notes:
        print("\n  sign tripwire:")
        for k, v in em.notes.most_common():
            print(f"    {k}")

    violations = load_schema().validate_batch(index.nodes.values(), em.edges)
    if violations:
        print(f"\nSCHEMA FAIL — {len(violations)} violations")
        for v in violations[:8]:
            print(f"  {v}")
        return 1
    print("\nschema: PASS")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("")
    (args.out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in em.edges))
    print(f"{len(em.edges):,} edges -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
