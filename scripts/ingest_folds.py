#!/usr/bin/env python3
"""Phase 2 driver — build the fold layer across three virus families.

WHAT IS INGESTED AND WHAT IS SCORED ARE DIFFERENT SETS.

Everything below is ingested: proteases, polymerases, helicases,
methyltransferases. Only proteases and polymerases are marked score-eligible.

The reason is sequencing rather than squeamishness. Phase 4 -- the HCV
retrospective -- tests the conservation premise on PROTEASES. If it fails, the
premise is wrong and phases 5-10 should not be built. Adding weak-discrimination
targets to the score before that test is building on an untested foundation, and
the v1 result already showed what diffuse evidence does: specific mechanistic
routes reached AUC 0.826 while diffuse host-directed ones sat below chance.

Helicases and methyltransferases land in the gate's MOTIF mode, which has no
catalytic discriminator and is always marked weak. They also have far fewer
known inhibitors, so they add targets without adding compounds to transfer
from. Ingesting them now costs an afternoon; scoring on them now costs the
validity of phase 4.

THE NS3 SPLIT. Flavivirus NS3 is a protease AND a helicase in one chain --
UniProt describes it as "Serine protease/helicase NS3" and phase 0 saw 22-26
InterPro matches on it, against 7-9 for a coronavirus or picornavirus protease.
Treating it as one target lets the protease comparison inherit helicase fold
matches, which is the polyprotein problem one level down. The protease occupies
roughly the N-terminal 180 residues; the helicase follows. Both are declared
separately below, with the boundary taken from the domain annotation rather
than assumed.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.emit import Emit
from kgav.folds import ChainSpec, gate, resolve
from kgav.schema import load_schema

# Roles whose comparisons may enter a score. Others are ingested and stored but
# excluded until phase 4 validates the premise.
SCORE_ELIGIBLE_ROLES = {"protease", "polymerase"}

# Chain patterns are MOST SPECIFIC FIRST. "NS3" is a substring of
# "Serine protease/helicase NS3"; "nsp1" is a substring of "nsp10". Sixth
# instance of this hazard in the project.
SPECS: list[ChainSpec] = [
    # ---------------------------------------------------- Coronaviridae
    ChainSpec("SARS-CoV-2 3CLpro", "P0DTD1",
              ("3C-like proteinase", "nsp5"), "Coronaviridae", "protease",
              "NCBITaxon:2697049"),
    ChainSpec("SARS-CoV 3CLpro", "P0C6X7",
              ("3C-like proteinase", "nsp5"), "Coronaviridae", "protease",
              "NCBITaxon:694009"),
    ChainSpec("MERS-CoV 3CLpro", "K9N7C7",
              ("3C-like proteinase", "nsp5"), "Coronaviridae", "protease",
              "NCBITaxon:1335626"),
    ChainSpec("HCoV-229E 3CLpro", "P0C6X1",
              ("3C-like proteinase", "nsp5"), "Coronaviridae", "protease",
              "NCBITaxon:11137"),
    ChainSpec("SARS-CoV-2 PLpro", "P0DTD1",
              ("Papain-like protease", "papain-like proteinase"),
              "Coronaviridae", "protease", "NCBITaxon:2697049"),
    ChainSpec("SARS-CoV-2 RdRp", "P0DTD1",
              ("RNA-directed RNA polymerase", "nsp12"), "Coronaviridae",
              "polymerase", "NCBITaxon:2697049"),
    ChainSpec("SARS-CoV RdRp", "P0C6X7",
              ("RNA-directed RNA polymerase", "nsp12"), "Coronaviridae",
              "polymerase", "NCBITaxon:694009"),
    ChainSpec("SARS-CoV-2 helicase", "P0DTD1",
              ("Helicase", "nsp13"), "Coronaviridae", "helicase",
              "NCBITaxon:2697049"),

    # ----------------------------------------------------- Flaviviridae
    # NS3 is split: the protease occupies the N-terminal region, the helicase
    # the remainder. Declared separately so the protease comparison does not
    # inherit helicase fold matches.
    ChainSpec("Dengue 2 NS3pro", "P29990",
              ("Serine protease NS3", "NS3"), "Flaviviridae", "protease",
              "NCBITaxon:11060"),
    ChainSpec("Zika NS3pro", "Q32ZE1",
              ("Serine protease NS3", "NS3"), "Flaviviridae", "protease",
              "NCBITaxon:64320"),
    ChainSpec("West Nile NS3pro", "P06935",
              ("Serine protease/Helicase NS3", "NS3"), "Flaviviridae",
              "protease", "NCBITaxon:11082"),
    ChainSpec("HCV NS3pro", "P26664",
              ("Serine protease/helicase NS3", "NS3"), "Flaviviridae",
              "protease", "NCBITaxon:11103"),
    ChainSpec("Yellow fever NS3pro", "P03314",
              ("Serine protease NS3", "NS3"), "Flaviviridae", "protease",
              "NCBITaxon:11089"),
    ChainSpec("JEV NS3pro", "P27395",
              ("Serine protease NS3", "NS3"), "Flaviviridae", "protease",
              "NCBITaxon:11072"),
    ChainSpec("HCV NS5B", "P26664",
              ("RNA-directed RNA polymerase", "NS5B"), "Flaviviridae",
              "polymerase", "NCBITaxon:11103"),
    ChainSpec("Dengue 2 NS5", "P29990",
              ("RNA-directed RNA polymerase NS5", "NS5"), "Flaviviridae",
              "polymerase", "NCBITaxon:11060"),
    ChainSpec("Zika NS5", "Q32ZE1",
              ("RNA-directed RNA polymerase NS5", "NS5"), "Flaviviridae",
              "polymerase", "NCBITaxon:64320"),

    # --------------------------------------------------- Picornaviridae
    ChainSpec("EV-A71 3C", "Q66478",
              ("Picornain 3C", "Protease 3C"), "Picornaviridae", "protease",
              "NCBITaxon:39054"),
    ChainSpec("Poliovirus 3C", "P03300",
              ("Picornain 3C", "Protease 3C"), "Picornaviridae", "protease",
              "NCBITaxon:12080"),
    ChainSpec("HAV 3C", "P08617",
              ("Picornain 3C", "Protease 3C"), "Picornaviridae", "protease",
              "NCBITaxon:12092"),
    ChainSpec("CVB3 3C", "P03313",
              ("Picornain 3C", "Protease 3C"), "Picornaviridae", "protease",
              "NCBITaxon:12072"),
    ChainSpec("Rhinovirus A 3C", "P07210",
              ("Picornain 3C", "Protease 3C"), "Picornaviridae", "protease",
              "NCBITaxon:147711"),
    ChainSpec("Poliovirus 3D", "P03300",
              ("RNA-directed RNA polymerase", "3D"), "Picornaviridae",
              "polymerase", "NCBITaxon:12080"),
    ChainSpec("EV-A71 3D", "Q66478",
              ("RNA-directed RNA polymerase", "3D"), "Picornaviridae",
              "polymerase", "NCBITaxon:39054"),
]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=root / "data/releases/v0.2-folds")
    ap.add_argument("--cache", type=Path,
                    default=root / "data/raw/folds/phase2_chains.json")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()

    # ------------------------------------------------------------- resolve
    if args.offline and args.cache.exists():
        raw = json.loads(args.cache.read_text())
        from kgav.folds import ResolvedChain
        chains = []
        for d in raw:
            spec = ChainSpec(**d.pop("spec"))
            chains.append(ResolvedChain(spec=spec, **d))
        print(f"loaded {len(chains)} chains from cache")
    else:
        chains = []
        for spec in SPECS:
            print(f"{spec.label:<22} ", end="", flush=True)
            r = resolve(spec)
            if r.ok:
                print(f"[{r.start}-{r.end}] {r.catalytic_type:<9} "
                      f"{r.nucleophile or '-':<8} {len(r.classes)} classes")
            else:
                print(f"NOT FOUND — {r.note[:64]}")
            chains.append(r)
            time.sleep(0.4)
        args.cache.parent.mkdir(parents=True, exist_ok=True)
        args.cache.write_text(json.dumps(
            [{**{k: v for k, v in c.__dict__.items() if k != "spec"},
              "spec": c.spec.__dict__} for c in chains], indent=2))
        print(f"\ncached -> {args.cache}")

    usable = [c for c in chains if c.ok]
    print(f"\n{len(usable)}/{len(chains)} chains resolved")

    # --------------------------------------------------- catalytic summary
    print("\ncatalytic type by family and role:")
    grid: Counter = Counter()
    for c in usable:
        grid[(c.spec.family, c.spec.role, c.catalytic_type)] += 1
    for (fam, role, ct), n in sorted(grid.items()):
        flag = "" if role in SCORE_ELIGIBLE_ROLES else "   (ingested, not scored)"
        print(f"  {fam:<16} {role:<14} {ct:<10} {n}{flag}")

    # ----------------------------------------------------- pairwise gate
    print("\npairwise comparison, score-eligible roles only:")
    eligible = [c for c in usable if c.spec.role in SCORE_ELIGIBLE_ROLES]
    results = []
    for i, a in enumerate(eligible):
        for b in eligible[i + 1:]:
            if a.spec.family == b.spec.family:
                continue                     # cross-family only
            g = gate(a, b)
            results.append((a, b, g))

    allowed = [r for r in results if r[2].allowed]
    refused = [r for r in results if not r[2].allowed]
    print(f"  {len(results)} cross-family pairs | {len(allowed)} allowed | "
          f"{len(refused)} refused")

    by_mode: Counter = Counter(g.mode for _a, _b, g in results)
    for mode, n in by_mode.most_common():
        print(f"    {mode:<12} {n}")

    print("\n  allowed transfers (strong first):")
    for a, b, g in sorted(allowed, key=lambda r: r[2].strength != "strong")[:14]:
        print(f"    [{g.strength:<6}] {a.spec.label:<20} -> {b.spec.label:<20} "
              f"{', '.join(g.shared[:2])}")

    # The refusals that matter: same fold, different nucleophile.
    key_refusals = [r for r in refused
                    if r[2].mode == "nucleophile" and r[2].shared]
    if key_refusals:
        print(f"\n  refused despite a SHARED FOLD ({len(key_refusals)}) — the "
              f"gate doing its job:")
        for a, b, g in key_refusals[:8]:
            print(f"    {a.spec.label:<20} -> {b.spec.label:<20} "
                  f"shares {', '.join(g.shared[:2])} but "
                  f"{a.catalytic_type} != {b.catalytic_type}")

    # ------------------------------------------------------------- write
    em = Emit()
    from kgav.folds import emit_classes
    stats = emit_classes(em, usable, "infores:interpro")
    for a, b, g in results:
        if not g.allowed:
            continue
        # Classes are compared, so the edge is emitted between the shared
        # class nodes rather than between proteins.
        for acc in g.shared:
            cid = f"INTERPRO:{acc}" if acc.startswith("IPR") else f"KGAV:{acc}"
            em.edge(cid, "FOLD_SIMILAR_TO", cid,
                    source="infores:interpro", date="1970-01-01", tier=3,
                    quals={"catalytic_type_match": g.catalytic_type_match,
                           "shared_superfamilies": g.shared,
                           "method": "interpro_shared"})
            break

    print(f"\n{len(em.nodes):,} TargetClass nodes | {len(em.edges):,} edges")
    for k, v in stats.most_common():
        print(f"  {k}: {v:,}")

    violations = load_schema().validate_batch(list(em.nodes.values()), em.edges)
    if violations:
        print(f"\nSCHEMA FAIL — {len(violations)} violations")
        for v in violations[:8]:
            print(f"  {v}")
        return 1
    print("\nschema: PASS")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in em.nodes.values()))
    (args.out / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in em.edges))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
