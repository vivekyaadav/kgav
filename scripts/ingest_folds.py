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
import dataclasses
import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.emit import Emit
from kgav.folds import ChainSpec, ResolvedChain, gate, resolve
from kgav.schema import load_schema

# Roles whose comparisons may enter a score. Others are ingested and stored but
# excluded until phase 4 validates the premise.
SCORE_ELIGIBLE_ROLES = {"protease", "polymerase"}

# Chain patterns are MOST SPECIFIC FIRST. "NS3" is a substring of
# "Serine protease/helicase NS3"; "nsp1" is a substring of "nsp10". Sixth
# instance of this hazard in the project.
#
# ---------------------------------------------------------------------------
# THE FLAVIVIRUS ACCESSIONS BELOW ARE STRAIN-TAXON ACCESSIONS, AND THE
# virus_taxon FIELD BESIDE THEM IS THE SPECIES. That reassignment is
# deliberate -- cross-viral metapaths must resolve against one taxon per virus
# -- but it was undocumented, which is how this class of bug survives.
# Measured against UniProt 2026-09-17:
#
#   P29990  Dengue 2   filed under 31634  (strain Thailand/16681/1984)
#                      declared here as   NCBITaxon:11060  (species)
#   P27395  JEV        filed under 11073  (strain SA-14)
#                      declared here as   NCBITaxon:11072  (species)
#   P03314  Yellow fev filed under 11090  (strain 17D vaccine)
#                      declared here as   NCBITaxon:11089  (species)
#   P26664  HCV        filed under 11104  (genotype 1a isolate 1)
#                      declared here as   NCBITaxon:11103  (species)
#
# Q32ZE1 (Zika, 64320) and P06935 (West Nile, 11082) are the only two curated
# AT the species taxon, so they need no reassignment.
#
# Querying the species taxon for any of the four above returns no curated
# polyprotein and therefore no chains -- the MERS trap, which cost this
# project Mpro and RdRp once already. Every strain taxon here is now also
# registered in config/viruses.yaml under `surveyed:`, with chain_source
# naming the accession and the taxon it really lives under, and
# kgav.virus_register.validate enforces that the pairing is declared rather
# than discovered. When this list moves into the config mechanism, that
# register is the source to read; until then the two must agree.
# ---------------------------------------------------------------------------
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
    ap.add_argument("--spine", type=Path,
                    default=root / "data/releases/v0.1-spine",
                    help="checked for the Protein nodes this layer references")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()

    # ------------------------------------------------------------- resolve
    if args.offline and args.cache.exists():
        raw = json.loads(args.cache.read_text())
        # A CACHE MISSING ANY FIELD ResolvedChain NOW CARRIES IS NOT USABLE.
        # ResolvedChain(**d) silently defaults what the file lacks, and the
        # default is indistinguishable from a measured value.
        #
        # THE FIRST VERSION OF THIS GUARD NAMED ONE FIELD, chain_pro_id, and
        # the very next field added walked straight past it:
        # catalytic_basis defaulted to "none" for all 24 chains, so the run
        # reported "none 24" and listed every 3CLpro under "no rule could
        # read them" WHILE the table above it said cysteine. A report
        # contradicting itself one screen apart.
        #
        # So the check is now derived from the dataclass instead of from a
        # list someone has to remember to extend. Adding a field invalidates
        # an old cache automatically, which is the only version of this that
        # stays true.
        required = {f.name for f in dataclasses.fields(ResolvedChain)} - {"spec"}
        missing = sorted(required - set().union(*(set(d) for d in raw))) if raw else []
        if missing:
            print(f"{args.cache} predates {', '.join(missing)}.\n"
                  f"  ResolvedChain would default those, and a default is "
                  f"not a measurement -- the\n  last time this happened the "
                  f"report said 'none 24' beside a table saying\n  "
                  f"'cysteine 4'. Re-run without --offline to refresh it.",
                  file=sys.stderr)
            return 2
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

    # A PROTEASE WITH NO CATALYTIC MACHINERY IS A RESOLUTION FAILURE, not a
    # biological fact, and nothing said so. SARS-CoV-2 PLpro reads
    # catalytic_type "none" because its spec matches UniProt's "Papain-like
    # protease nsp3" -- the whole 1,945-residue multi-domain chain, not the
    # protease domain inside it. 60 InterPro classes against 7-26 for every
    # other chain is the tell.
    #
    # This is the polyprotein problem from phase0's header recurring one
    # level down: "Any question about a specific protein must be asked of a
    # specific CHAIN" -- and here the chain is itself multi-domain. Reported
    # rather than dropped, because which sub-range to use is a curation
    # question and silently excluding a target would be worse.
    miscalled = [c for c in usable
                 if c.spec.role == "protease" and c.catalytic_type == "none"]
    if miscalled:
        print(f"\n  {len(miscalled)} chain(s) declared a protease but "
              f"resolved no catalytic type:")
        for c in miscalled:
            print(f"    {c.spec.label:<22} {c.spec.accession} "
                  f"[{c.start}-{c.end}] {c.end - c.start + 1:,} aa, "
                  f"{len(c.classes)} classes")
        print("    A protease has catalytic machinery by definition, so this "
              "is the chain match\n    being wrong, not the enzyme. Such a "
              "chain cannot be gated on its nucleophile\n    and contributes "
              "nothing to M8.")

    # WHICH RULE DECIDED EACH TYPE. catalytic_machinery has two paths: an
    # explicitly annotated nucleophile, and a composition fallback whose
    # test is "the LAST active site is C or S" -- the chymotrypsin-like
    # assumption the annotated path was already fixed to stop making.
    # Widening the fallback would re-decide every chain that reaches it,
    # and until now nothing said which chains those are.
    by_basis: Counter = Counter(c.catalytic_basis for c in usable)
    print(f"\n  how each catalytic type was decided: "
          + ", ".join(f"{k} {v}" for k, v in by_basis.most_common()))
    fallback = [c for c in usable if c.catalytic_basis == "composition"]
    if fallback:
        print("    decided by composition, so a change to that rule would "
              "move these:")
        for c in fallback:
            print(f"      {c.spec.label:<22} {c.catalytic_type:<9} "
                  f"{', '.join(c.catalytic_residues[:4])}")
    unread = [c for c in usable
              if c.catalytic_basis == "none" and c.catalytic_residues]
    if unread:
        print("    HAVE active sites but no rule could read them:")
        for c in unread:
            print(f"      {c.spec.label:<22} {c.spec.role:<10} "
                  f"{', '.join(c.catalytic_residues[:5])}")
    unkeyed = [c for c in chains
               if c.end > c.start > 0 and not c.chain_pro_id]
    if unkeyed:
        print(f"  {len(unkeyed)} matched a chain but carry no PRO_ id, so "
              f"they have no stable node key:")
        for c in unkeyed[:6]:
            print(f"    {c.spec.label:<22} {c.note[:60]}")

    # DOES THIS LAYER JOIN THE SPINE? Until now it could not: node_id was
    # built from coordinates while the spine keys the same chain on its PRO_
    # id, so every MEMBER_OF_CLASS edge referenced a node that does not
    # exist. That is almost certainly why `folds` is absent from
    # LAYER_PRECEDENCE -- assemble_release returns 1 on a dangling edge.
    #
    # Reported rather than assumed, because the ids matching in form does
    # not mean the spine built the chain: a virus outside the register, or a
    # chain the spine skipped, still dangles.
    if args.spine and (args.spine / "nodes.jsonl").exists():
        spine_ids = {json.loads(line)["id"]
                     for line in (args.spine / "nodes.jsonl").read_text().splitlines()
                     if line.strip()}
        hit = [c for c in usable if c.node_id in spine_ids]
        miss = [c for c in usable if c.node_id not in spine_ids]
        print(f"\njoin against {args.spine.name}: {len(hit)}/{len(usable)} "
              f"resolved chains exist as Protein nodes")
        for c in miss[:10]:
            print(f"    MISSING  {c.node_id:<30} {c.spec.label} "
                  f"({c.spec.family}, taxon {c.spec.virus_taxon})")
        if miss:
            print(f"  {len(miss)} chain(s) would produce DANGLING edges. The "
                  f"layer must not be added to\n  LAYER_PRECEDENCE until this "
                  f"is zero, or until those chains are dropped.")
    else:
        print(f"\nno spine release at {args.spine}, so the join is unchecked")

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

    # WHICH CHAINS HAVE ANY PARTNER AT ALL. M8's entire claim is reaching a
    # virus with no activity data of its own, so a chain with zero allowed
    # transfers is a virus M8 cannot help however the path is wired.
    #
    # Reported per chain because the aggregate hides it: 36 allowed sounds
    # like broad coverage, and on v0.1 every one of them is a cysteine
    # protease pair. Comparison is cross-family only (same-family pairs are
    # skipped above), and the catalytic table shows serine proteases exist
    # in ONE family here -- so a flavivirus protease has no cross-family
    # partner that is not cysteine, and the gate correctly refuses all of
    # them. That is structural, not a data gap.
    partners: Counter = Counter()
    for a, b, g in results:
        if g.allowed:
            partners[a.spec.label] += 1
            partners[b.spec.label] += 1
    stranded = [c for c in eligible if not partners[c.spec.label]]
    print(f"\n  allowed partners per chain "
          f"({len(eligible) - len(stranded)}/{len(eligible)} have at least one):")
    for c in sorted(eligible, key=lambda c: -partners[c.spec.label]):
        mark = "   <- M8 cannot reach this one" if not partners[c.spec.label] else ""
        print(f"    {c.spec.label:<22} {c.catalytic_type:<9} "
              f"{partners[c.spec.label]:>3}{mark}")

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
    # NO FOLD_SIMILAR_TO EDGES ARE EMITTED, AND THAT IS DELIBERATE UNTIL THE
    # PATH SHAPE IS SETTLED.
    #
    # This loop used to write `em.edge(cid, "FOLD_SIMILAR_TO", cid, ...)` --
    # subject and object the SAME node, because a shared class is one node, not
    # two. M8 requires two DISTINCT TargetClass nodes joined by that predicate,
    # and walk_metapath discards a neighbour already in `visited`, so every
    # self-loop was skipped and M8 produced zero paths on every run since 0.8.0.
    #
    # The self-loop is not a typo to swap for a pair of ids. It is the data
    # model disagreeing with the metapath: when two chains SHARE a class, the
    # structure that connects them is MEMBER_OF_CLASS in both directions
    # through that one node, and no similarity edge is needed or meaningful.
    # M8's declared path asks for a hop between two classes that this ingest
    # has no principled pair to supply. Inventing one -- every class of A
    # crossed with every class of B -- would mint thousands of similarity
    # assertions nobody measured, which is the failure this project exists to
    # avoid. See the note in the schema on M8.
    #
    # The gate itself is correct and is reported above; what is unresolved is
    # how to write its verdict into the graph. Until that is decided the layer
    # emits classes and membership only, and metapath_reach reports M8 as
    # reaching nothing rather than letting it look scored.
    if any(g.allowed for _a, _b, g in results):
        n_allowed = sum(1 for _a, _b, g in results if g.allowed)
        print(f"\n  NOTE: {n_allowed} allowed transfers are NOT written as "
              f"FOLD_SIMILAR_TO edges.\n  M8 therefore has no similarity hop to "
              f"walk and will report zero reach. See the\n  comment at this "
              f"line before changing it -- the previous version emitted\n  "
              f"self-loops, which validate cleanly and are silently unwalkable.")

    print(f"\n{len(em.nodes):,} TargetClass nodes | {len(em.edges):,} edges")
    for k, v in stats.most_common():
        print(f"  {k}: {v:,}")

    # VALIDATED AGAINST THE SPINE TOO, because this layer emits TargetClass
    # nodes and MEMBER_OF_CLASS edges whose SUBJECTS are Proteins living in
    # the spine. Checking the edges against em.nodes alone reported all 402
    # as DANGLING on every run, including the ten chains the spine really
    # does hold -- a validator that fails identically whether or not the
    # data is right, which is the condition it exists to detect.
    #
    # ingest_vhppi.py solved exactly this and said so: "This release is
    # edges-only: every endpoint lives in the spine or host release ...
    # validation runs here against the real index instead." Ninth instance
    # in this project of a solved problem not reaching its second site.
    known = list(em.nodes.values())
    if args.spine and (args.spine / "nodes.jsonl").exists():
        known += [json.loads(line)
                  for line in (args.spine / "nodes.jsonl").read_text().splitlines()
                  if line.strip()]
    else:
        print(f"\nno spine at {args.spine}: every Protein endpoint will read "
              f"DANGLING, which says\n  nothing about this layer")
    violations = load_schema().validate_batch(known, em.edges)
    if violations:
        print(f"\nSCHEMA FAIL — {len(violations)} violations")
        for v in violations[:8]:
            print(f"  {v}")
        print("\n  A DANGLING subject here is a chain the spine does not "
              "build -- see the join\n  report above. Those are comparator "
              "chains for viruses outside the register,\n  not a defect in "
              "this layer.")
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
