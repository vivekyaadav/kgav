#!/usr/bin/env python3
"""What would change if these surveyed viruses were promoted to `viruses:`.

Writes nothing. Promotes in memory and reports the STRUCTURAL consequences --
what validates, what breaks, which layers would run, and what M8 would be
expected to do. The data volumes are a separate question already answered by
diagnose_flavivirus_join.py, and this script deliberately does not recompute
them: duplicating one measurement in two places is the failure mode this
repository has paid for five times (AUC 0.826, the coverage floor, lift's
trustworthy rule, the taxon map in four scripts, and the gate this report
itself used to assert).

    python scripts/dryrun_promote.py --promote DENV-2 Zika
    python scripts/dryrun_promote.py --promote DENV-2 Zika "West Nile"

WHY A DRY RUN AT ALL. Promotion rebuilds the release, which restates every
figure in docs/evaluation.md and changes what the agent answers. A wrong
assumption should cost a read, not a rebuild.
"""
from __future__ import annotations

import argparse
import ast
import copy
import gzip
import importlib.util
import json
import os
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.assemble import LAYER_PRECEDENCE
from kgav.channels import layer_predicates, reachable_for_virus
from kgav.virus_register import declared_taxa, validate

_CFG_CACHE: dict = {}
SCHEMA: dict = yaml.safe_load((ROOT / "schema" / "kg_schema.yaml").read_text())


def _phase0_targets() -> dict:
    """phase0_verify_folds.TARGETS, loaded by path because scripts/ is not a
    package. Returns {} if it cannot be loaded, rather than guessing its
    contents."""
    path = ROOT / "scripts" / "phase0_verify_folds.py"
    if not path.exists():
        return {}
    try:
        spec = importlib.util.spec_from_file_location("_p0", path)
        mod = importlib.util.module_from_spec(spec)
        # REGISTER BEFORE EXEC. phase0_verify_folds defines a @dataclass, and
        # dataclasses resolves sys.modules[cls.__module__] while processing
        # it. Without this line that lookup returns None and the import dies
        # with "'NoneType' object has no attribute '__dict__'" -- which says
        # nothing about the real cause until you read the traceback.
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        return dict(getattr(mod, "TARGETS", {}) or {})
    except Exception as exc:                       # noqa: BLE001
        print(f"  (could not load phase0_verify_folds: {exc})")
        return {}


def _hardcoded_counts() -> list[tuple[str, int, str]]:
    """Assertions that pin the number of built viruses.

    test_surveyed_viruses_are_not_built asserts len(doc["viruses"]) == 7. A
    promotion that does not update it fails the suite, which is the gate
    working -- but it should be found here rather than by a red test.
    """
    out = []
    pat = re.compile(r'len\(\s*doc\[\s*["\']viruses["\']\s*\]\s*\)\s*==\s*(\d+)')
    for path in sorted((ROOT / "tests").glob("*.py")):
        for n, line in enumerate(path.read_text().splitlines(), 1):
            m = pat.search(line)
            if m:
                out.append((f"{path.name}:{n}", int(m.group(1)), line.strip()))

    # AND THE NAMED SET THAT REPLACED THE COUNT. test_virus_register.BUILT
    # lists every virus the build is meant to contain. Reporting only the
    # count made this section print "no test pins the number of built
    # viruses" the moment the count was replaced by something better --
    # a dry run that goes quiet when the guard improves is worse than no
    # dry run, because quiet reads as safe.
    # ANY comparison against the current built count, anywhere in tests.
    # The narrow scan above found one breakage; the promotion of DENV-2 and
    # Zika broke THREE, because test_ingest_spine.py asserts `len(taxa) == 7`
    # and `len(labels) == 7` -- same magic number, different expression. A
    # scan that only recognises one spelling under-reports, and a dry run
    # that under-reports is the thing it exists to prevent.
    #
    # Reported as CANDIDATES, not failures: `== 7` also appears for the seven
    # metapaths in test_calibration.py, which has nothing to do with viruses.
    # Over-reporting costs a glance; under-reporting cost a broken build.
    # TWO COUNTS, NOT ONE. This scanned for the built-virus count only, and
    # the 2026-10-12 Zika promotion broke four tests it had not named. One of
    # them, test_declared_taxa_still_defaults_to_the_built_section, asserted
    # `len(built) == 11` -- the size of the TAXON MAP, which grows with a
    # promotion too and is a different number from the virus count because of
    # isolate_taxa and chain_source. The other three asserted output ORDER or
    # string membership, which no count scan can see; those are what
    # `make test` is for, and it must be run after the register edit and
    # BEFORE the ingest chain. That ordering is the whole lesson of the first
    # promotion attempt.
    counts = set()
    if _CFG_CACHE:
        counts.add(len(_CFG_CACHE["viruses"]))
        counts.add(len(declared_taxa(_CFG_CACHE)))
    for built_now in sorted(counts):
        loose = re.compile(rf"==\s*{built_now}\b")
        for path in sorted((ROOT / "tests").glob("*.py")):
            for n, line in enumerate(path.read_text().splitlines(), 1):
                if not loose.search(line) or "assert" not in line:
                    continue
                ref = f"{path.name}:{n}"
                if any(ref == e[0] for e in out):
                    continue
                out.append((ref, None, line.strip()))

    tv = ROOT / "tests" / "test_virus_register.py"
    if tv.exists():
        tree = ast.parse(tv.read_text())
        for node in ast.walk(tree):
            if (isinstance(node, ast.Assign)
                    and any(isinstance(x, ast.Name) and x.id == "BUILT"
                            for x in node.targets)
                    and isinstance(node.value, (ast.Tuple, ast.List))):
                names = [e.value for e in node.value.elts
                         if isinstance(e, ast.Constant)]
                out.append((
                    f"test_virus_register.py:{node.lineno} BUILT",
                    len(names),
                    f"BUILT names {len(names)}: {', '.join(names)}"))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "config/viruses.yaml")
    ap.add_argument("--promote", nargs="+", required=True,
                    help="names exactly as they appear under `surveyed:`")
    ap.add_argument("--proteomes", type=Path,
                    default=ROOT / "data/raw/proteomes")
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data/raw/chembl/chembl_37"
                                   "/chembl_37_sqlite/chembl_37.db")
    args = ap.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    _CFG_CACHE.update(cfg)
    surveyed = {e["name"]: e for e in cfg.get("surveyed") or []}

    unknown = [n for n in args.promote if n not in surveyed]
    if unknown:
        sys.exit(f"not in `surveyed:`: {unknown}\n"
                 f"available: {sorted(surveyed)}")

    promoted = copy.deepcopy(cfg)
    moving = [surveyed[n] for n in args.promote]
    promoted["viruses"] = list(promoted["viruses"]) + moving
    promoted["surveyed"] = [e for e in promoted["surveyed"]
                            if e["name"] not in set(args.promote)]

    print("DRY RUN -- nothing is written\n")
    print(f"promoting {len(moving)}: {', '.join(args.promote)}")
    print(f"  viruses:  {len(cfg['viruses'])} -> {len(promoted['viruses'])}")
    print(f"  surveyed: {len(cfg['surveyed'])} -> {len(promoted['surveyed'])}")

    # --------------------------------------------------------- preconditions
    # THE SECTION THIS SCRIPT DID NOT HAVE, AND THE OMISSION COST A BUILD.
    # ingest_spine reads a LOCAL CACHE: data/raw/proteomes/<proteome_id>.json.gz,
    # and there is no fetcher in this repository. When the file is absent it
    # records proteome_missing:<virus>, returns early, and leaves the virus
    # with an OrganismTaxon node and nothing else.
    #
    # The damage is downstream, not there. declared_taxa now contains the new
    # taxa, so ingest_chembl queries them, finds thousands of activities and
    # emits INHIBITS edges to flavivirus proteins that no layer created:
    # 1,881 DANGLING, SCHEMA FAIL, layer not written. assemble_release then
    # runs happily on the OLD chembl layer plus the new orphan taxon nodes,
    # producing a release whose edges are unchanged, whose nodes_sha is not,
    # and which invalidates every results file in exchange for nothing.
    #
    # One stat() per virus would have caught it.
    print("\n" + "=" * 78)
    print("0. PRECONDITIONS")
    print("=" * 78)
    blocked = []
    for e in moving:
        pid = e.get("proteome_id")
        pfile = args.proteomes / f"{pid}.json.gz"
        ok = pfile.exists()
        if not ok:
            blocked.append(e["name"])
        print(f"  {'ok ' if ok else 'ABSENT':<7} {e['name']:<14} "
              f"{pid}.json.gz")
        cs = e.get("chain_source") or {}
        acc = cs.get("accession")
        extras = e.get("extra_accessions") or []

        # THE EXTRA FILE, WHEN THE REGISTER DECLARES ONE. ingest_spine reads
        # EXTRA_<species taxon>.json.gz and its curated entries REPLACE the
        # proteome -- but the pfile.exists() guard returns early first, so
        # both files are required.
        if extras:
            xfile = args.proteomes / f"EXTRA_{e['taxon']}.json.gz"
            xok = xfile.exists()
            if not xok:
                blocked.append(e["name"])
            print(f"  {'ok ' if xok else 'ABSENT':<7} {'':<14} "
                  f"EXTRA_{e['taxon']}.json.gz  for {', '.join(extras)}")

        if ok:
            try:
                entries = (json.loads(gzip.decompress(pfile.read_bytes()))
                           .get("results", []))
                print(f"          {len(entries)} entries; register claims "
                      f"{cs.get('chains')} chains on {acc}")
                # DOES THE FILE ACTUALLY CONTAIN THE CHAIN SOURCE? This is
                # the check that would have caught DENV-2 locally. P29990
                # has no proteome cross-reference at all and is filed under
                # strain 31634, so UP000095876 does not contain it and the
                # build produces a virus with zero chains -- which the
                # register's own comment predicted and nothing enforced.
                accs = {x.get("primaryAccession") for x in entries}
                if acc and acc not in accs and acc not in extras:
                    print(f"          MISSING: {acc} is not in this proteome "
                          f"and no extra_accessions declares it, so the "
                          f"chains\n          this virus needs are "
                          f"unreachable")
                    blocked.append(e["name"])
                chains = [x for x in entries
                          if x.get("primaryAccession") == acc
                          for f in x.get("features", [])
                          if f.get("type") == "Chain"]
                if chains:
                    print(f"          {len(chains)} CHAIN features on {acc} "
                          f"in the file")

                # CAN ANYTHING REACH THE TAXON? ingest_virus gates the whole
                # Gene block on `if symbol`, so an entry with no `genes`
                # field and no declared fallback yields no Gene node, hence
                # no ENCODED_BY and no BELONGS_TO. M1 is INHIBITS ->
                # ENCODED_BY -> BELONGS_TO, so the virus is then absent from
                # the direct-acting channel entirely.
                #
                # The 2026-10-12 Zika build is why this check exists. Q32ZE1
                # has no `genes` field, the build printed "genes=0", schema
                # PASSed, 3,653 rows resolved to NS3 -- and nothing could
                # walk from a compound to the Zika taxon. One read of the
                # cached file would have said so beforehand, which is what
                # this now does.
                declared = e.get("gene_symbol")
                with_genes = [x.get("primaryAccession") for x in entries
                              if x.get("genes")]
                if with_genes:
                    print(f"          gene symbols from UniProt on "
                          f"{len(with_genes)}/{len(entries)} entries")
                elif declared:
                    print(f"          no UniProt gene names; register "
                          f"declares gene_symbol: {declared}")
                else:
                    print(f"          NO GENE SYMBOL: no entry has a `genes` "
                          f"field and the register declares no\n          "
                          f"gene_symbol, so ingest_spine emits no Gene node "
                          f"-- no ENCODED_BY, no\n          BELONGS_TO, and "
                          f"M1 cannot reach this taxon at all")
                    blocked.append(e["name"])
            except (OSError, ValueError) as exc:
                print(f"          UNREADABLE: {exc}")
                blocked.append(e["name"])

    if blocked:
        print(f"\n  BLOCKED: {', '.join(sorted(set(blocked)))}")
        print("  ingest_spine would record proteome_missing and build an "
              "OrganismTaxon with no\n  proteins. ingest_chembl would then "
              "query the new taxa, emit INHIBITS edges\n  to proteins no "
              "layer created, and fail its schema check -- after which\n"
              "  assemble_release still writes a release, invalidating every "
              "results file for\n  no new facts. Fetch the proteome(s) "
              "first.")
        # WAS: "There is no fetcher in this repo, so whatever produced the
        # existing files should produce these". True when written, false since
        # scripts/fetch_raw.py landed, and a precondition section telling the
        # reader there is no way to satisfy the precondition is worse than
        # silence.
        print("\n  FETCH THEM WITH:")
        for e in moving:
            print(f"    python scripts/fetch_raw.py --virus {e['name']}")
        print("\n  fetch_raw validates the response shape, checks the entry's "
              "Chain count\n  against this register's chain_source, and "
              "refuses to write a file that\n  disagrees -- so a bad download "
              "cannot reach ingest_spine. `make check-raw`\n  re-validates "
              "what is already on disk without network.")
    else:
        print("\n  every promoted virus has its proteome cached")

    other = []
    unavailable: set[str] = set()
    if len(os.environ.get("ORCS_KEY", "")) != 32:
        other.append("ORCS_KEY is not a 32-character key: ingest_orcs.py "
                     "exits without writing,\n    so the flavivirus "
                     "screens would not be ingested even once the patterns "
                     "match")
        unavailable.add("orcs")
    if not list((ROOT / "data" / "raw" / "vhppi").glob("vhn_*.tab27")):
        other.append("no vhn_*.tab27 files: ingest_vhppi has nothing to read")
        unavailable.add("vhppi")
    if not args.db.exists():
        other.append(f"no ChEMBL db at {args.db}: ingest_chembl and "
                     f"ingest_selectivity cannot run")
    for o in other:
        print(f"\n  NOTE: {o}")

    # WHICH CHANNELS EACH MISSING INPUT ACTUALLY BLOCKS.
    #
    # The notes above used to be an undifferentiated list, which made
    # ORCS_KEY and the VirHostNet files read as blockers on the whole
    # flavivirus effort. They are not: M1 is
    # SmallMolecule -INHIBITS-> Protein -ENCODED_BY-> Gene -BELONGS_TO->
    # OrganismTaxon, so it needs chembl and spine and nothing host-mediated.
    # Derived from the schema's own metapath declarations and from the
    # predicates each layer release actually emitted -- not a list anyone
    # maintains, because this is precisely the text that goes stale.
    print("\n  WHICH CHANNELS WOULD REACH THE PROMOTED VIRUS")
    # LAYER_PRECEDENCE is the project's own statement of what a layer is, so
    # the temporal training splits (v0.1-train2021*) are excluded by it rather
    # than by a list kept here. Whatever is skipped gets printed.
    layer_preds, skipped = layer_predicates(ROOT / "data" / "releases", "v0.1",
                                            include=set(LAYER_PRECEDENCE))
    if skipped:
        print(f"    not layers, skipped: {', '.join(skipped)}")
    if not layer_preds:
        print("    no v0.1-* layer releases on disk; cannot derive")
    else:
        # WHICH LAYERS WILL HOLD DATA ABOUT *THIS* VIRUS.
        #
        # The first version of this report asked only whether some layer
        # anywhere emits a predicate, and so announced "M3 evaluable" for
        # Zika: PHYSICALLY_INTERACTS_WITH exists in the host layer with
        # 150,432 edges. But M3 reaches the taxon through a VIRAL protein, and
        # host's PPIs are human-human -- only vhppi bridges host to virus, and
        # it has no flavivirus files. Wired is not the same as reachable.
        covering = {"spine"} if not blocked else set()
        vhn = sorted((ROOT / "data" / "raw" / "vhppi").glob("vhn_*.tab27"))
        promoted_taxa = {t for e in moving
                         for t in ([str(e["taxon"])]
                                   + [str(x) for x in (e.get("isolate_taxa") or [])]
                                   + [str((e.get("chain_source") or {}).get("taxon")
                                          or e["taxon"])])}
        if vhn and any(f"taxid:{t}" in f.read_text() for f in vhn
                       for t in promoted_taxa):
            covering.add("vhppi")
        if "orcs" not in unavailable:
            covering.add("orcs")
        if args.db.exists():
            try:
                from kgav.chembl import query_activities
                if query_activities(args.db, sorted(promoted_taxa)):
                    covering.add("chembl")
            except Exception as e:                   # noqa: BLE001 - reported
                print(f"    could not check ChEMBL coverage: {e}")
        print(f"    layers that will hold data for "
              f"{', '.join(e['name'] for e in moving)}: "
              f"{', '.join(sorted(covering)) or 'none'}")

        reach = reachable_for_virus(SCHEMA.get("metapaths") or {},
                                    layer_preds, covering)
        for code in sorted(reach):
            i = reach[code]
            src = ", ".join(sorted(i["sources"])) or "no layer emits it"
            if i["reachable"]:
                print(f"    {code:<4} REACHES IT   enters via {i['entry']} "
                      f"({src})")
            else:
                print(f"    {code:<4} cannot reach it: enters via "
                      f"{i['entry']} ({src})")
        yes = sorted(c for c, i in reach.items() if i["reachable"])
        if yes:
            print(f"\n    {', '.join(yes)} would reach it. M1 is the "
                  f"direct-acting channel and the\n    only one carrying a "
                  f"defensible result here, so a first evaluation is not\n"
                  f"    waiting on ORCS or VirHostNet -- those gate the "
                  f"host-mediated channels.")
        else:
            print("\n    No channel would reach it. Promoting now adds a "
                  "taxon and no evidence.")

    # ------------------------------------------------------------- register
    print("\n" + "=" * 78)
    print("1. REGISTER VALIDATION")
    print("=" * 78)
    violations = validate(promoted)
    if not violations:
        print("  validate(): no violations")
    for v in violations:
        print(f"  {v.code:<26} {v.ref}")
        print(f"      {v.message}")

    # the families each entry declares, against the defaults block
    print(f"\n  defaults.family: {(cfg.get('defaults') or {}).get('family')}")
    for e in moving:
        fam = e.get("family") or "(inherits defaults)"
        print(f"    {e['name']:<14} family: {fam}")
    print("  An entry declaring its own family overrides the default for that\n"
          "  entry only. Nothing in src/ or scripts/ branches on family, so\n"
          "  this is metadata rather than behaviour -- checked, not assumed.")

    # ------------------------------------------------------------ taxon map
    print("\n" + "=" * 78)
    print("2. TAXON MAP  (declared_taxa -> every ingest that resolves a taxon)")
    print("=" * 78)
    before, after = declared_taxa(cfg), declared_taxa(promoted)
    print(f"  {len(before)} -> {len(after)} entries")
    for tax in sorted(set(after) - set(before)):
        print(f"    + {tax:<10} -> {after[tax]}")

    # ------------------------------------------------------- what breaks
    print("\n" + "=" * 78)
    print("3. WHAT BREAKS")
    print("=" * 78)
    counts = _hardcoded_counts()
    if not counts:
        print("  no test pins the number of built viruses")
    for ref, n, line in counts:
        if n is None:
            print(f"  {'CANDIDATE':<11} {ref}  compares against "
                  f"{len(cfg['viruses'])}, the current built count")
        else:
            verdict = "OK" if n == len(promoted["viruses"]) else "WOULD FAIL"
            print(f"  {verdict:<11} {ref}  pins {n}, would be "
                  f"{len(promoted['viruses'])}")
        print(f"              {line}")
    print("\n  test_surveyed_viruses_are_not_built is the gate that keeps a\n"
          "  surveyed virus out of a build. A promotion updates it; it does\n"
          "  not get deleted.")

    # ----------------------------------------------------------- the layers
    print("\n" + "=" * 78)
    print("4. LAYERS THAT WOULD RUN")
    print("=" * 78)
    scripts = sorted(p.stem.removeprefix("ingest_")
                     for p in (ROOT / "scripts").glob("ingest_*.py"))
    for s in scripts:
        key = s.replace("_", "")
        listed = key in LAYER_PRECEDENCE
        print(f"  {'listed  ' if listed else 'UNLISTED'}  ingest_{s}.py")
    missing = [s for s in scripts if s.replace("_", "") not in LAYER_PRECEDENCE]
    if missing:
        print("\n  assemble_release.py builds its layer list FROM "
              "LAYER_PRECEDENCE and filters\n  by existence, so an unlisted "
              "layer is not merely last -- it is silently\n  absent while "
              "every count still looks right. "
              "test_every_ingest_script_has_a\n  _precedence_entry is a "
              "strict xfail on exactly this.")

    # --------------------------------------------------------------- M8
    print("\n" + "=" * 78)
    print("5. WHAT M8 WOULD BE EXPECTED TO DO")
    print("=" * 78)
    targets = _phase0_targets()
    # ONE ACCESSION, MANY ROLES. A flavivirus genome is a SINGLE polyprotein,
    # so chain_source.accession carries the protease and the polymerase and
    # everything else: P29990 is pinned in phase0 twice, as "Dengue 2 NS3"
    # (serine protease) and "Dengue 2 NS5" (polymerase). Keying a dict on the
    # accession kept whichever came last and reported DENV-2's expected
    # nucleophile as None, hiding the protease -- the half this whole section
    # is about. Coronaviruses do not collide this way, because ORF1ab and the
    # structural proteins are separate accessions, which is why a dict looked
    # sufficient when I wrote it.
    by_acc: dict[str, list[tuple[str, dict]]] = {}
    for _name, _spec in targets.items():
        by_acc.setdefault(_spec["acc"], []).append((_name, _spec))
    print("  folds.py records that the protease fold is SHARED across families")
    print("  (SSF50494, G3DSA:2.40.10.10 for Coronaviridae 3CLpro,")
    print("  Picornaviridae 3C and Flaviviridae NS3 alike) and that what")
    print("  separates them is the catalytic nucleophile. The gate's")
    print("  NUCLEOPHILE mode requires an identical residue.\n")
    print(f"  {'virus':<14}{'accession':<11}{'role':<11}"
          f"{'nucleophile':<13}pinned as")
    for e in moving:
        cs = e.get("chain_source") or {}
        acc = cs.get("accession") or "?"
        hits = sorted(by_acc.get(acc) or [], key=lambda h: h[0])
        if not hits:
            print(f"  {e['name']:<14}{acc:<11}{'-':<11}{'not pinned':<13}-")
            continue
        for i, (name, spec_) in enumerate(hits):
            print(f"  {(e['name'] if i == 0 else ''):<14}"
                  f"{(acc if i == 0 else ''):<11}"
                  f"{spec_.get('role')!s:<11}"
                  f"{(spec_.get('expect') or '-- none --')!s:<13}{name}")
        if len(hits) > 1:
            print(f"  {'':<14}{'':<11}-> {len(hits)} roles on ONE accession: "
                  f"a flavivirus polyprotein, {cs.get('chains')} chains")
    corona = [t for t in targets.values()
              if t.get("family") == "Coronaviridae" and t.get("role") == "protease"]
    if corona and any(by_acc.get((e.get("chain_source") or {}).get("accession"))
                      for e in moving):
        print("\n  PRE-REGISTERED PREDICTION, and it is a NEGATIVE one:")
        print("  coronavirus 3CLpro is cysteine, flavivirus NS3 is serine, so")
        print("  NUCLEOPHILE mode must emit NO protease edge between them.")
        print("  Polymerases have no equivalent discriminator (motif C is GDD")
        print("  vs SDD), so MOTIF mode may emit edges and must mark them WEAK.")
        print("  phase0_verify_folds.py already labels the flavivirus")
        print("  proteases 'NEGATIVE CONTROL'.")
        fams_now = sorted({(v.get("family")
                            or (cfg.get("defaults") or {}).get("family"))
                           for v in cfg["viruses"]})
        fams_after = sorted({(v.get("family")
                              or (cfg.get("defaults") or {}).get("family"))
                             for v in promoted["viruses"]})
        print(f"\n  FAMILIES IN THE BUILD: {', '.join(fams_now)}")
        if fams_after != fams_now:
            print(f"  after promotion:       {', '.join(fams_after)}")
        if len(fams_now) < 2:
            print("  A cross-family channel needs a SECOND family to reach, so"
                  " M8 emits\n  nothing today for that reason as well as the "
                  "unlisted layer -- two\n  faults wearing one label.")
        else:
            print("  Two or more families are already built, so the layer "
                  "being unlisted in\n  LAYER_PRECEDENCE is now the only "
                  "thing keeping M8 empty.")

    print("\n" + "=" * 78)
    print("WHAT THIS DOES NOT DO")
    print("=" * 78)
    print("""  Nothing is written and no virus is promoted. Volumes -- compounds,
  protein rows, publications, ORCS screens, VirHostNet coverage -- come from
  diagnose_flavivirus_join.py and are not recomputed here. Whether the
  chain_source accessions still carry the chain counts the register claims is
  a UniProt question this script does not ask; ingest_spine does, and
  chain_source.verified records when it was last measured.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
