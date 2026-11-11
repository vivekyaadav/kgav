"""The nine conflated nodes, and how far into the evaluation they reach.

docs/scope_and_limits.md states this limit, and the first draft of that
section asserted the nine were unreachable by any evaluated channel and so
could not affect a published figure. MEASURED, 8 of 9 carry an evaluation
label and 5 carry a traversable edge; one holds 16 TARGETS edges. The claim
was written from the assumption that nine tautomer oddities would be
peripheral, and they are ordinary screened drugs.

The comfortable version was the wrong one, which is the reason these three
counts are asserted against the release instead of stated in prose.
"""
from __future__ import annotations

import collections
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "data/releases/v0.1"
DOC = ROOT / "docs/scope_and_limits.md"
LABELS = {"HAS_ANTIVIRAL_ACTIVITY_AGAINST", "MEASURED_INACTIVE_AGAINST"}


def conflated_ids() -> set[str]:
    """The nine, recomputed from the raw ncats input via the one owner."""
    import csv
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.ncats import substance_id
    from kgav.normalize_chem import normalize_chemical, tautomer_conflations
    cache: dict[str, object] = {}
    derivations: dict[str, dict[str, set[str]]] = {}
    for name in ("cpe_assay.csv", "cpe_tox.csv", "cpe.tsv"):
        f = ROOT / "data/raw/ncats" / name
        if not f.exists():
            continue
        delim = "\t" if f.suffix == ".tsv" else ","
        with f.open(newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh, delimiter=delim):
                keys = {k.lower().replace("_", ""): k for k in row if k}
                sk, ik = keys.get("smiles"), keys.get("sampleid")
                if not sk:
                    continue
                smi = (row.get(sk) or "").strip()
                if not smi:
                    continue
                if smi not in cache:
                    try:
                        cache[smi] = normalize_chemical(smi)
                    except Exception:
                        cache[smi] = None
                cid = cache[smi]
                if cid is None:
                    continue
                derivations.setdefault(f"INCHIKEY:{cid.inchikey}", {}) \
                    .setdefault(smi, set()).add(substance_id((row.get(ik) or "").strip()))
    return set(tautomer_conflations(derivations))


@pytest.mark.slow
@pytest.mark.skipif(not (RELEASE / "edges.jsonl").exists()
                    or not (ROOT / "data/raw/ncats/cpe_assay.csv").exists(),
                    reason="needs the release and data/raw")
def test_the_documents_three_counts_are_the_releases():
    ids = conflated_ids()
    by = collections.defaultdict(collections.Counter)
    for line in (RELEASE / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        for end in (e["subject"], e["object"]):
            if end in ids:
                by[end][e["predicate"]] += 1
    labelled = sum(1 for n in ids
                   if any(p in LABELS for p in by.get(n, ())))
    reachable = sum(1 for n in ids
                    if any(p not in LABELS for p in by.get(n, ())))
    text = " ".join(DOC.read_text().split())

    assert f"**{len(ids)} nodes that policy pools" in text or \
           re.search(rf"\*\*{len(ids)} nodes that policy pools", text), (
        f"the release has {len(ids)} conflated nodes and the document says "
        f"otherwise")
    assert f"**{labelled} of the {len(ids)} carry a" in text, (
        f"{labelled} of {len(ids)} carry an evaluation label; the document "
        f"states a different count")
    assert f"{reachable} carry a traversable edge" in text, (
        f"{reachable} of {len(ids)} carry a traversable edge; the document "
        f"states a different count")
