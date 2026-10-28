#!/usr/bin/env python
"""The druggable host-factor landscape, without a scorer.

WHY THIS IS A DIFFERENT QUESTION FROM THE ONE THAT FAILED. Every protocol in
docs/evaluation.md asked which drug works against a virus. That needs the
antiviral labels, 92.2% of which share a publication with the evidence edge
that was supposed to predict them, and no channel survived the control.

This asks what the drugged human proteome overlaps with what each virus
depends on. It reads TARGETS, PARTICIPATES_IN and HOST_FACTOR_FOR, and it
never touches HAS_ANTIVIRAL_ACTIVITY_AGAINST or any other label predicate.
The shared-PMID confound has nothing to attach to here.

THE CONFOUND MOVES RATHER THAN VANISHING, and that is the whole reason for
the matched test below. Well-studied proteins carry more pathway
annotations, are likelier to be drug targets, AND are likelier to come out
of a CRISPR screen with enough coverage to be called a dependency factor.
An over-representation test against the whole proteome measures that, not
biology. So every enrichment is reported twice: naive, and against a
background sampled to match the observed degree distribution. The gap
between them is a measurement, not a footnote.

EVERY MOLECULE IS ACCOUNTED FOR. The graph holds far more small molecules
than carry a human target, and a landscape over the subset that silently
called itself "all drugs" would be the defect this repository keeps
finding. The census prints the whole population and names what each part
can and cannot contribute.

WHAT THIS CANNOT SAY. ChEMBL's bulk activity table records potency, not
whether a compound agonises or antagonises its target. Inhibiting a
dependency factor should block infection; activating one should not. So
"drug targets a dependency factor" is not "drug would help", and the census
reports how many edges carry a known direction so the gap is visible.

USAGE

  python scripts/pathway_landscape.py
  python scripts/pathway_landscape.py --virus NCBITaxon:2697049
  python scripts/pathway_landscape.py --permutations 10000
  python scripts/pathway_landscape.py --out data/results/pathway_landscape.json

Reads the release only, ranks no drug, recommends no compound. With --out it
writes the overlap table to a results file so the agent can state the measured
null instead of implying a ranking the landscape refuses -- see
src/kgav/landscape.py, which loads it and will not quote a run from a
different release. Without --out it writes nothing, as before.

A --virus run writes a PARTIAL file and marks it so: a landscape covering one
taxon cannot support "no taxon is enriched", and the loader must not be handed
a file that looks complete because it was filtered.
"""
from __future__ import annotations

import argparse
import random
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.baselines import Graph  # noqa: E402
from kgav.schema import load_schema  # noqa: E402
from kgav.provenance import write_results
from kgav.vhppi import BINARY_METHODS

# THE ONE IMPLEMENTATION, imported. A second copy of an exact hypergeometric
# tail would drift from the one Phase 4 validated against scipy over 400
# random tables.
from phase4_hcv_retrospective import fisher_one_sided  # noqa: E402

TARGETS = "TARGETS"
PARTICIPATES_IN = "PARTICIPATES_IN"
HOST_FACTOR_FOR = "HOST_FACTOR_FOR"
PHYSICALLY_INTERACTS_WITH = "PHYSICALLY_INTERACTS_WITH"
DEPENDENCY = "dependency"
DEGREE_BINS = 10
MIN_PATHWAY_TARGETS = 5
DEFAULT_PERMUTATIONS = 2000


# ------------------------------------------------------------------ census
def census(g: Graph) -> dict:
    """Every small molecule, by what it can contribute.

    Printed before anything else because the analysable set is a minority
    of the graph's molecules, and a landscape that reported only the
    minority while calling it "all drugs" would be false in the first line.
    """
    molecules = {n for n in g.nodes if g.klass(n) == "SmallMolecule"}
    has_targets, has_inhibits, directed = set(), set(), 0
    for (pred, node), edges in g.out.items():
        if node not in molecules:
            continue
        if pred == TARGETS:
            has_targets.add(node)
            for _nbr, e in edges:
                q = e.get("qualifiers") or {}
                if q.get("direction") and q["direction"] != "unknown":
                    directed += 1
        elif pred == "INHIBITS":
            has_inhibits.add(node)
    both = has_targets & has_inhibits
    neither = molecules - has_targets - has_inhibits
    # HOW MUCH OF THE UNREACHABLE SET A DRUG-TARGET SOURCE COULD REACH.
    # 86% of the molecules carry no human target, and "add another source"
    # is only worth doing if that source covers them. DrugCentral, IUPHAR
    # and DrugBank curate APPROVED and clinical-stage drugs, not bulk
    # screening compounds, so the clinical phase of the unreachable set is
    # the ceiling on what they could add -- measured here rather than
    # guessed at.
    phases: Counter = Counter()
    for m in (molecules - has_targets):
        ph = (g.nodes[m].get("properties") or {}).get("max_phase")
        if ph is None:
            phases["unknown"] += 1
        elif ph >= 4:
            phases["approved"] += 1
        elif ph >= 2:
            phases["clinical"] += 1
        elif ph >= 1:
            phases["phase 1"] += 1
        else:
            phases["preclinical"] += 1
    # A PARTITION, and `both` is NOT a part of it. The first version
    # returned four counts that read as four buckets; on v0.1 they summed
    # to 21,714 against 21,667 molecules, because the 47 compounds with
    # both edge kinds are already inside `analysable`. The test asserted
    # the four summed and passed, because the fixture had no compound with
    # both -- the defect was in the fixture, not only in the arithmetic.
    return {
        "molecules": len(molecules),
        "analysable": len(has_targets),
        "inhibits_only": len(has_inhibits - has_targets),
        "both": len(both),
        "neither": len(neither),
        "parts_sum": len(has_targets) + len(has_inhibits - has_targets)
        + len(neither),
        "target_edges": sum(len(g.out.get((TARGETS, m), ())) for m in has_targets),
        "directed_edges": directed,
        "untargeted_by_phase": phases,
    }


# ------------------------------------------------------------- the surface
def drug_targets(g: Graph) -> set[str]:
    """Human proteins reached by at least one TARGETS edge."""
    out = set()
    for (pred, node), edges in g.out.items():
        if pred != TARGETS:
            continue
        for nbr, _e in edges:
            if g.klass(nbr) == "Protein":
                out.add(nbr)
    return out


def human_proteins(g: Graph) -> set[str]:
    """Proteins that could in principle be drugged or be a host factor.

    The first version sampled the null from EVERY Protein node, viral ones
    included. A viral protein can never be a human drug target nor a host
    dependency factor, so each one added to a degree bin is a guaranteed
    miss that lowers the matched expectation and flatters the observed
    overlap. The universe is the human proteome.
    """
    return {n for n in g.nodes
            if g.klass(n) == "Protein" and not g.is_viral(n)}


def dependency_factors(g: Graph, virus: str) -> set[str]:
    """Proteins the virus DEPENDS on, by CRISPR screen.

    The direction qualifier is required, not assumed: a restriction factor
    is the opposite claim, and pooling the two would make every count mean
    nothing.
    """
    out = set()
    for (pred, node), edges in g.out.items():
        if pred != HOST_FACTOR_FOR:
            continue
        for nbr, e in edges:
            if nbr != virus:
                continue
            if (e.get("qualifiers") or {}).get("direction") == DEPENDENCY:
                out.add(node)
    return out


def viral_proteins(g: Graph) -> dict[str, set[str]]:
    """{virus taxon: its protein nodes}, for the interactome channel."""
    out: dict[str, set[str]] = defaultdict(set)
    for nid, node in g.nodes.items():
        if node.get("class") != "Protein":
            continue
        props = node.get("properties") or {}
        if not props.get("is_viral"):
            continue
        taxon = props.get("taxon_id")
        if taxon:
            out[str(taxon)].add(nid)
    return dict(out)


def interactome_partners(g: Graph, proteins: set[str],
                         hub_cut: int | None = None
                         ) -> tuple[set[str], set[str]]:
    """Human proteins physically interacting with any of `proteins`.

    A SECOND EVIDENCE CHANNEL, NOT A SECOND SOURCE FOR THE FIRST. A
    dependency factor is a gene the virus NEEDS: knock it out and
    replication drops. A physical interactor is a protein the viral
    protein TOUCHES. Neither implies the other -- a virus binds plenty it
    does not depend on, and depends on plenty it never touches -- so the
    two sets are measured separately and never pooled. Pooling them would
    make the count mean nothing, the same reason dependency_factors
    refuses to pool restriction factors with dependency ones.

    It exists because the CRISPR layer reaches six of nine flaviviruses
    not at all: DENV-1, DENV-3, DENV-4, JEV, TBEV and Yellow fever have
    no screen in this ORCS release, so the dependency channel is silent
    about them by construction. Four of those now carry virus-host PPI.

    THE HOST SIDE IS WHAT MATTERS HERE, which is why edges attributed by
    the register's ppi_accession_map are usable despite resolving only to
    a polyprotein: the human partner is identified exactly either way.

    Returns (partners, hubs). Hubs are returned rather than silently
    dropped: ubiquitin and ribosomal proteins are the classic sticky
    partners of pull-down data, and how much of an interactome they
    account for is a property of the assay worth reporting. The release
    keeps them -- assemble_release reports the cut rather than applying
    it -- so the decision is the consumer's and this one states it.
    """
    partners: set[str] = set()
    hubs: set[str] = set()
    for vp in sorted(proteins):
        for nbr, _e in g.neighbours(PHYSICALLY_INTERACTS_WITH, vp):
            if g.klass(nbr) != "Protein" or g.is_viral(nbr):
                continue
            if hub_cut and g.degree[nbr] > hub_cut:
                hubs.add(nbr)
                continue
            partners.add(nbr)
    return partners, hubs


def binary_share(g: Graph, proteins: set[str]) -> tuple[int, int]:
    """(edges from a binary assay, edges in total) into `proteins`.

    WHAT THE EVIDENCE SAYS, carried forward from vhppi.py: VirHostNet's
    SARS-CoV-2 rows are dominated by proximity labelling, which marks
    anything within about 10 nm. "Physically interacts with" in this
    graph mostly means "sits in the same complex neighbourhood", and a
    pathway result built on it should not be read as binding.
    """
    binary = total = 0
    for vp in sorted(proteins):
        for nbr, e in g.neighbours(PHYSICALLY_INTERACTS_WITH, vp):
            if g.klass(nbr) != "Protein" or g.is_viral(nbr):
                continue
            total += 1
            if (e.get("qualifiers") or {}).get("detection_method") in BINARY_METHODS:
                binary += 1
    return binary, total


def host_factor_coverage(g: Graph) -> dict[str, Counter]:
    """{taxon: Counter of HOST_FACTOR_FOR directions}.

    WHY THIS IS PRINTED RATHER THAN INFERRED. The overlap table drops a
    taxon with no DEPENDENCY factor, and on v0.1 three dropped: HCoV-HKU1,
    SARS-CoV and Zika. "No dependency factor" has at least two causes that
    look identical from outside -- no screen ever reached the taxon, or
    screens reached it and every hit was a RESTRICTION factor, which is the
    opposite claim. Those call for completely different work, so the report
    separates them instead of leaving the reader to guess.

    Zika in particular is matched by src/kgav/orcs.py's VIRUS_PATTERNS,
    added 2026-10-11, so a zero there is a finding about the build rather
    than about the virus.
    """
    out: dict[str, Counter] = defaultdict(Counter)
    for (pred, node), edges in g.out.items():
        if pred != HOST_FACTOR_FOR:
            continue
        for nbr, e in edges:
            d = (e.get("qualifiers") or {}).get("direction") or "unstated"
            out[nbr][d] += 1
    return dict(out)


def pathway_members(g: Graph, cap: int | None) -> dict[str, set[str]]:
    """{pathway: proteins in it}, dropping pathways over the schema's cap.

    The cap is the schema's own `max_pathway_members`, for the reason
    baselines._oversized_pathways gives: a GO term with 1,389 members makes
    every drug adjacent to every host factor, and an enrichment over it
    measures the annotation, not the biology.
    """
    members: dict[str, set[str]] = defaultdict(set)
    for (pred, node), edges in g.out.items():
        if pred != PARTICIPATES_IN:
            continue
        for nbr, _e in edges:
            members[nbr].add(node)
    if cap:
        members = {p: m for p, m in members.items() if len(m) <= cap}
    return dict(members)


# --------------------------------------------------- the matched test
def study_degree(g: Graph) -> Counter:
    """Host-host PPI degree: a study-volume proxy independent of the test.

    WHY THE DEFAULT DEGREE IS NOT SAFE HERE. Graph.degree counts every
    edge of every predicate, so it includes the edges that DEFINE both
    sides of this comparison:

      TARGETS                     ~8.8 per drugged protein, and the
                                  drugged set is the observed set
      PHYSICALLY_INTERACTS_WITH   several per virus partner, and the
                                  partners are the hit set

    Matching on that degree is circular twice over. Being drugged pushes
    a protein into a high bin for a reason unrelated to touching a
    virus; being a partner pushes it into a high bin by definition, so
    the high bins the null draws from are enriched for hits. The
    expectation inflates, the observed does not, and every ratio lands
    under 1 -- which is what the 2026-10-28 run shows for all nine taxa
    in the interactome channel, six of them at the 1/(N+1) floor.

    This counts only PHYSICALLY_INTERACTS_WITH edges with NO viral
    endpoint: the STRING host-host interactome. Well-studied proteins
    have more of them, which is the confound worth controlling, and a
    host-host edge is in neither the observed set nor the hit set.

    Decided from the NODES rather than from the edge's source label.
    STRING and VirHostNet write the same predicate, and tonight's
    provenance work is the reminder that a label can be wrong while a
    node property cannot be.
    """
    deg: Counter = Counter()
    for (pred, node), edges in g.out.items():
        if pred != PHYSICALLY_INTERACTS_WITH or g.is_viral(node):
            continue
        for nbr, _e in edges:
            if g.is_viral(nbr):
                continue
            deg[node] += 1
            deg[nbr] += 1
    return deg


def degree_bins(g: Graph, universe: set[str], bins: int = DEGREE_BINS,
                degree: Counter | None = None
                ) -> dict[str, int]:
    """Assign each protein to a degree quantile bin.

    Quantiles rather than equal-width: node degree is heavy-tailed, so
    equal-width bins would put almost everything in the first one and match
    nothing.
    """
    # THE NODE ID IS PART OF THE KEY, not decoration. Degree is heavily
    # tied -- thousands of proteins share a low degree -- and sorting on
    # degree alone leaves ties in whatever order the SET yielded, which
    # python varies between processes through string hash randomisation.
    # Three runs of this script on one release and one seed returned
    # SARS-CoV-2 expectations of 353.3, 352.6 and 352.4.
    deg = g.degree if degree is None else degree
    ordered = sorted(universe, key=lambda n: (deg.get(n, 0), n))
    out = {}
    n = len(ordered)
    for i, node in enumerate(ordered):
        out[node] = min(bins - 1, i * bins // n) if n else 0
    return out


def matched_null(observed: set[str], universe: set[str],
                 bins: dict[str, int], hits: set[str], draws: int,
                 rng: random.Random) -> list[int]:
    """Overlap with `hits` for `draws` degree-matched random sets.

    Each draw takes, from every degree bin, as many proteins as `observed`
    has in that bin. So the null set looks like the drug-target set in
    study volume and differs only in which proteins it picked.
    """
    # SORTED, for the same reason degree_bins sorts: rng.sample draws by
    # position, so a pool in set order makes the seed stop determining the
    # draw. A figure in a paper has to come back the same on a rerun.
    by_bin: dict[int, list[str]] = defaultdict(list)
    for node in sorted(universe):
        by_bin[bins[node]].append(node)
    want = Counter(bins[n] for n in observed if n in bins)
    # SORTED AGAIN, and this one is the subtlest of the three. A Counter
    # iterates in INSERTION order, and it was built by walking `observed`,
    # which is a set -- so the bins came out in a different sequence per
    # process, rng.sample was called in a different order, and the draws
    # differed even after the pools themselves were sorted.
    wanted = sorted(want.items())
    out = []
    for _ in range(draws):
        total = 0
        for b, k in wanted:
            pool = by_bin[b]
            if k >= len(pool):
                total += sum(1 for p in pool if p in hits)
                continue
            total += sum(1 for p in rng.sample(pool, k) if p in hits)
        out.append(total)
    return out


def overlap_row(targets: set[str], universe: set[str], bins: dict[str, int],
                factors: set[str], permutations: int,
                rng: random.Random) -> tuple[dict, float]:
    """One channel's overlap between the drugged proteome and `factors`.

    ONE IMPLEMENTATION FOR BOTH CHANNELS. The interactome channel asks
    the identical question of a different set -- how much of the drugged
    proteome lands in it, against a degree-matched background -- so it
    takes this rather than a copy. Two copies of a statistic drift, and
    this project has already paid that twice over: two trustworthy
    rules, five taxon maps.
    """
    a = len(targets & factors)
    b = len(targets) - a
    cc = len(factors - targets)
    d = len(universe) - a - b - cc
    p_naive = fisher_one_sided(a, b, cc, d, greater=True)
    null = matched_null(targets, universe, bins, factors, permutations, rng)
    exp = sum(null) / len(null) if null else 0.0
    # BOTH TAILS, from the same null. Six of the eight taxa on v0.1 come
    # back under ratio 1.0 and two far under, and with only the greater
    # tail there was no way to tell a depletion from an ordinary low draw.
    row = {"n_dependency": len(factors), "n_observed": a, "expected": exp,
           "p_naive": p_naive,
           "p_matched": empirical_p(a, null),
           "p_matched_less": empirical_p(a, null, greater=False)}
    # The ratio is NOT stored: a third number that can disagree with the
    # two it is derived from, and infinite whenever the matched
    # background expects nothing, which JSON cannot write.
    return row, (a / exp if exp else float("inf"))


def empirical_p(observed: int, null: list[int], greater: bool = True) -> float:
    """One-sided, with the +1 correction.

    Without it a p of 0 is reportable, which claims more precision than the
    number of permutations can support.

    `greater` picks the tail. The enrichment tail was the only one computed
    for a long time, and the agent built on this file reported every sub-1
    ratio as "not a measured depletion" because there was no lower tail to
    measure one with. Default True so every existing caller keeps the
    behaviour it was written against.

    BOTH TAILS FROM ONE NULL is deliberate: the same permutation draw answers
    both questions, so a second run would only add sampling noise between
    two numbers that must be consistent with each other.
    """
    if greater:
        extreme = sum(1 for v in null if v >= observed)
    else:
        extreme = sum(1 for v in null if v <= observed)
    return (extreme + 1) / (len(null) + 1)


# ----------------------------------------------------------------- report
def _write(out: Path, args, rows: dict[str, dict],
           uncovered: list[str], **channels) -> None:
    """Record the overlap table with the run that produced it.

    `partial` is the field that matters. A --virus run measures one taxon, and
    a file that did not say so would let the agent state "no taxon is
    enriched" on the strength of a single row. The loader refuses a partial
    file for exactly that reason.
    """
    # THROUGH write_results, LIKE EVERY OTHER RESULTS WRITER. This wrote
    # its own stamp -- measured_on, measured_by, and `release` as the
    # directory NAME -- which is strictly weaker than the shared one and
    # failed in exactly the way kgav.provenance was written to prevent:
    # v0.1 was reassembled in place on 2026-10-27, keeping its name while
    # becoming a different graph, and a file recording "release": "v0.1"
    # cannot tell the two apart. The shared stamp fingerprints the CONTENT
    # of nodes.jsonl and edges.jsonl, so it can.
    #
    # measured_on and measured_by are kept: the loader reads the first and
    # the second names the script, which the provenance block's `code`
    # (a commit) does not.
    doc = {
        "measured_on": date.today().isoformat(),
        "measured_by": "scripts/pathway_landscape.py",
        "release": args.release.name,
        "permutations": args.permutations,
        "seed": args.seed,
        # WHICH YARDSTICK. Two files measured on different bases are not
        # comparable, and a file that does not say which it used cannot
        # be checked -- the same failure as a results file recording a
        # release by name.
        "degree_basis": args.degree_from,
        "partial": bool(args.virus),
        "partial_reason": (f"run with --virus {args.virus}, so it covers that "
                           f"taxon only" if args.virus else ""),
        "taxa": rows,
        "uncovered": uncovered,
        # SEPARATE KEYS, NEVER MERGED INTO `taxa`. A reader or loader that
        # knows only the dependency channel must not silently pick up
        # interactome numbers under the keys it already reads: the two
        # answer different questions and a mixed table means nothing.
        **channels,
    }
    write_results(out, doc, args.release)
    print(f"\n  wrote {out} — {len(rows)} taxon(s) measured, "
          f"{len(uncovered)} with no dependency factor"
          + ("  [PARTIAL: one taxon only]" if args.virus else ""))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path, default=ROOT / "data/releases/v0.1")
    ap.add_argument("--degree-from", choices=("all", "study"), default="all",
                    help="what the degree-matched null matches ON. 'all' is "
                         "total graph degree, which includes the edges "
                         "defining both the observed set and the hits. "
                         "'study' uses host-host PPI only, which is in "
                         "neither. See study_degree.")
    ap.add_argument("--virus", default=None)
    ap.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=None,
                    help="write the overlap table here for src/kgav/"
                         "landscape.py to load")
    args = ap.parse_args(argv)

    if not (args.release / "nodes.jsonl").exists():
        print(f"no release at {args.release}")
        return 2

    schema = load_schema()
    g = Graph.load(args.release, symmetry=schema.symmetry())
    rng = random.Random(args.seed)

    print("PATHWAY LANDSCAPE. No scorer, no ranking, no antiviral label.")
    print("Reads TARGETS, PARTICIPATES_IN and HOST_FACTOR_FOR only.\n")
    print(f"release {args.release.name} | {len(g.nodes):,} nodes\n")

    c = census(g)
    print("  CENSUS — every small molecule in the graph")
    print(f"    {c['molecules']:8,}  small molecules in total")
    print(f"    {c['analysable']:8,}  carry a human TARGETS edge  <- the "
          f"analysable set")
    print(f"    {c['inhibits_only']:8,}  carry INHIBITS only, no human target")
    print(f"    {c['neither']:8,}  carry neither")
    print(f"    {'':8}  ---")
    print(f"    {c['parts_sum']:8,}  those three partition the population")
    print(f"    {c['both']:8,}  of the analysable set ALSO carry INHIBITS "
          f"(a subset, not a fourth part)")
    assert c["parts_sum"] == c["molecules"], "the census must partition"
    share = c["analysable"] / c["molecules"] if c["molecules"] else 0
    print(f"\n    The landscape covers {share:.1%} of the molecules. The rest")
    print("    have no human target in this graph, so they have no pathway to")
    print("    place them in. That is a property of the source, not a bug.")
    ph = c["untargeted_by_phase"]
    reachable = ph["approved"] + ph["clinical"]
    print(f"    Of the {sum(ph.values()):,} with no human target:")
    for k in ("approved", "clinical", "phase 1", "preclinical", "unknown"):
        if ph[k]:
            print(f"      {ph[k]:8,}  {k}")
    print(f"\n    A curated drug-target source (DrugCentral, IUPHAR,")
    print(f"    DrugBank) covers approved and clinical-stage drugs, so")
    print(f"    {reachable:,} is the CEILING on what one could add here,")
    print(f"    not {sum(ph.values()):,}. The preclinical and unknown rest")
    print(f"    are bulk screening compounds no such source curates.\n")
    print(f"    {c['target_edges']:,} TARGETS edges, of which "
          f"{c['directed_edges']:,} state a direction.")
    print("    The rest are potency without agonist/antagonist, so targeting")
    print("    a dependency factor is NOT a claim that the drug would help.\n")

    targets = drug_targets(g)
    universe = human_proteins(g)
    members = pathway_members(g, schema.retrieval_limits.get(
        "max_pathway_members"))
    in_pathways = {p for ms in members.values() for p in ms}
    print(f"  DRUGGED PROTEOME")
    print(f"    {len(targets):,} distinct human proteins carry a TARGETS edge")
    print(f"    {len(members):,} pathways at or under the schema's "
          f"{schema.retrieval_limits.get('max_pathway_members')}-member cap")
    print(f"    {len(targets & in_pathways):,} of the targets sit in at least "
          f"one such pathway\n")

    viruses = [args.virus] if args.virus else sorted(
        n for n in g.nodes if g.klass(n) == "OrganismTaxon")

    cov = host_factor_coverage(g)
    kinds = sorted({k for c in cov.values() for k in c} | {DEPENDENCY})
    print("  HOST-FACTOR COVERAGE, per taxon")
    print(f"    {'taxon':<18} " + " ".join(f"{k:>12}" for k in kinds)
          + f" {'total':>7}")
    for v in viruses:
        c = cov.get(v, Counter())
        total = sum(c.values())
        flag = "" if total else "   <- NO screen reached this taxon"
        print(f"    {v:<18} " + " ".join(f"{c.get(k, 0):12,}" for k in kinds)
              + f" {total:7,}{flag}")
    print("    A taxon absent from the overlap table below has no DEPENDENCY")
    print("    factor. That is not the same as having no screen: a taxon with")
    print("    restriction hits only has data and the opposite finding, and")
    print("    the two need different work.\n")
    print("  OVERLAP WITH DEPENDENCY FACTORS, naive against degree-matched")
    print(f"  {args.permutations:,} permutations, seed {args.seed}\n")
    print(f"    {'virus':<18} {'dep':>6} {'drugged':>8} {'exp':>7} "
          f"{'ratio':>6} {'p naive':>9} {'p up':>10} {'p down':>10}")
    any_virus = False
    uncovered = []
    rows: dict[str, dict] = {}
    # HOISTED. This was rebuilt inside the loop, once per virus, from the
    # same two arguments. Identical output, 16x less work.
    basis = study_degree(g) if args.degree_from == "study" else None
    bins = degree_bins(g, universe, degree=basis)
    if args.degree_from == "study":
        print(f"  degree matched on HOST-HOST PPI only "
              f"({sum(basis.values()) // 2:,} edges), not on total graph "
              f"degree:\n    total degree includes the TARGETS edges that "
              f"define the drugged set and the\n    virus-host edges that "
              f"define the interactome, so matching on it is circular.\n")
    for virus in viruses:
        deps = dependency_factors(g, virus)
        if not deps:
            uncovered.append(virus)
            continue
        any_virus = True
        row, ratio = overlap_row(targets, universe, bins, deps,
                                 args.permutations, rng)
        rows[virus] = row
        print(f"    {virus:<18} {len(deps):6,} {row['n_observed']:8,} "
              f"{row['expected']:7.1f} {ratio:6.2f} {row['p_naive']:9.2e} "
              f"{row['p_matched']:10.4f} {row['p_matched_less']:10.4f}")
    # ---------------------------------------------- the interactome channel
    # A SECOND CHANNEL, REPORTED BESIDE THE FIRST AND NEVER MERGED INTO IT.
    # Six of nine flaviviruses have no CRISPR screen in this ORCS release,
    # so the dependency channel is silent about them by construction. Four
    # of those six carry virus-host PPI as of 2026-10-27, which is a
    # different question about the same virus: what it TOUCHES rather than
    # what it NEEDS.
    hub_cut = int((schema.retrieval_limits or {}).get("hub_degree_cut") or 0)
    vprot = viral_proteins(g)
    irows: dict[str, dict] = {}
    iuncovered = []
    icov: dict[str, dict] = {}
    for virus in viruses:
        partners, hubs = interactome_partners(g, vprot.get(virus, set()),
                                              hub_cut or None)
        binary, total = binary_share(g, vprot.get(virus, set()))
        icov[virus] = {"n_partners": len(partners), "n_hubs": len(hubs),
                       "binary_edges": binary, "total_edges": total}
        if not partners:
            iuncovered.append(virus)
            continue
        row, ratio = overlap_row(targets, universe, bins, partners,
                                 args.permutations, rng)
        row["n_hubs_excluded"] = len(hubs)
        row["binary_edges"] = binary
        row["total_edges"] = total
        irows[virus] = row

    print("  INTERACTOME COVERAGE, per taxon  (host proteins a viral protein")
    print(f"  physically interacts with; hubs over degree {hub_cut or 'n/a'} "
          f"excluded)")
    print(f"    {'taxon':<18} {'partners':>9} {'hubs cut':>9} "
          f"{'edges':>8} {'binary':>8}")
    for virus in viruses:
        c = icov[virus]
        flag = "" if c["n_partners"] else "   <- no virus-host PPI for this taxon"
        print(f"    {virus:<18} {c['n_partners']:9,} {c['n_hubs']:9,} "
              f"{c['total_edges']:8,} {c['binary_edges']:8,}{flag}")
    print("    A BINARY edge is two-hybrid or pull-down: it reads out a")
    print("    contact. The rest are proximity or co-complex, which say the")
    print("    two sit in the same neighbourhood. This channel does not")
    print("    separate them, so nothing here means 'binds'.\n")

    if irows:
        print("  OVERLAP WITH INTERACTION PARTNERS, naive against "
              "degree-matched")
        print(f"  {args.permutations:,} permutations, seed {args.seed}\n")
        print(f"    {'virus':<18} {'partners':>8} {'drugged':>8} {'exp':>7} "
              f"{'ratio':>6} {'p naive':>9} {'p up':>10} {'p down':>10}")
        for virus, row in irows.items():
            exp = row["expected"]
            ratio = row["n_observed"] / exp if exp else float("inf")
            print(f"    {virus:<18} {row['n_dependency']:8,} "
                  f"{row['n_observed']:8,} {exp:7.1f} {ratio:6.2f} "
                  f"{row['p_naive']:9.2e} {row['p_matched']:10.4f} "
                  f"{row['p_matched_less']:10.4f}")
        print()
        print("    THIS IS NOT THE TABLE ABOVE WITH MORE VIRUSES. A gene the")
        print("    virus needs and a protein it touches are different claims,")
        print("    measured from different experiments, and a virus may")
        print("    appear in one table and not the other for reasons that")
        print("    have nothing to do with its biology.")
    if iuncovered:
        print(f"\n    {len(iuncovered)} of {len(viruses)} taxa have no "
              f"virus-host PPI and are absent above:")
        for v in iuncovered:
            print(f"      {v}")

    if args.out:
        _write(args.out, args, rows, uncovered,
               interactome=irows, interactome_uncovered=iuncovered,
               interactome_coverage=icov)
    if not any_virus:
        print("    NO VIRUS HAS DEPENDENCY FACTORS in this release, so there")
        print("    is no overlap to measure. That is the finding, not an")
        print("    empty table.")
        return 1

    # A VIRUS THAT DROPS OUT OF THE TABLE IS A SCOPE LIMIT, NOT A BLANK.
    # The first version skipped these silently, so a reader counted five
    # rows against eight taxa and had nothing to tell them why.
    if uncovered:
        print()
        print(f"    {len(uncovered)} of {len(viruses)} taxa have NO dependency "
              f"factor at all and are absent above:")
        for v in uncovered:
            print(f"      {v}")
        print("    The CRISPR layer does not cover them, so the landscape")
        print("    has nothing to say about those taxa. That is a scope")
        print("    limit of the source and belongs in the paper's methods.")

    print()
    print("  'exp' is how many dependency factors a DEGREE-MATCHED random")
    print("  protein set of the same size hits. Where the matched p is far")
    print("  weaker than the naive one, the naive test was measuring study")
    print("  volume: well-studied proteins are both likelier to be drugged")
    print("  and likelier to be called dependency factors.")
    print()
    print("  Nothing here ranks a drug or claims one would work. Direction is")
    print("  unknown on most TARGETS edges, and a CRISPR dependency factor is")
    print("  a target worth inhibiting, not evidence that any given compound")
    print("  does inhibit it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
