#!/usr/bin/env python
"""WHAT RESTATING THE TASK WOULD COST (fix-list B3), measured per virus.

B3 asks whether the evaluated task should be restated as "compounds with a
known mechanism (~30% of labels)". That is a decision about what the project
claims, not a refactor, and it was filed with an estimate nobody could check.
This measures it.

"KNOWN MECHANISM" IS AMBIGUOUS BY AN ORDER OF MAGNITUDE, which is the first
thing this reports and the reason the row's "~30%" cannot be read as one
number. Four definitions, all defensible, all different:

  mechanism    ChEMBL drug_mechanism names an action_type for the compound.
               The literal reading of "known mechanism". Needs --db.
  targets      the compound carries a TARGETS edge in this release.
  mech_or_tgt  either of the above.
  any_edge     the compound has ANY traversable edge besides its own label,
               which is the 71% figure section B of the fix list quotes.

WHAT TO READ FIRST. Not the kept-label counts -- the two columns beside
them. `cells` is how many (virus, task) cells still hold at least
MIN_REACHED_POSITIVES positives AND as many negatives, because a cell below
that cannot describe a scorer rather than its ceiling. `pos/neg ratio` is
whether the restriction is correlated with the LABEL: a definition that
keeps 30% of positives and 21% of negatives has selected on the answer, and
every AUC computed under it inherits that. Fix-list B1 was declined for
exactly this, measured the same way, and the lesson generalises: an
admission rule drawn from annotation is drawn from what somebody chose to
study.

THE SIGN IS NOT CONSTANT ACROSS VIRUSES, which is worse than a bias. Most
viruses' labels are ChEMBL-only, where the annotated compounds are the
antiviral research compounds -- so restricting raises the positive rate.
SARS-CoV-2 carries the NCATS panel, whose negatives are library compounds
that DO carry host-target and mechanism annotations -- so restricting
lowers it. `build_labels` already warns about this in as many words; this
shows what it costs.

AND PRECISION IS THE OTHER HALF OF THE DECISION. A smaller task is not a
worse task if it can still answer the question, so the last section reports,
for each surviving cell, the AUC it would take to exclude chance at all.

Cross-checks, because a measurement nothing contradicts is a measurement
nothing checks: the per-virus label counts are compared against the
n_positives / n_negatives independently recorded in
data/results/baselines_*.json, and the TARGETS objects are counted by taxon
rather than assumed to be human.

Needs the release. --db additionally needs the 30 GB gitignored chembl
sqlite; without it the `mechanism` columns are reported as unavailable
rather than silently folded into the others.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import Graph                       # noqa: E402
from kgav.calibration import MIN_REACHED_POSITIVES     # noqa: E402
from kgav.labels import auc_interval, build_labels     # noqa: E402
from kgav.schema import load_schema                     # noqa: E402
from kgav.pipeline import release_dir

LABEL_PREDICATES = ("HAS_ANTIVIRAL_ACTIVITY_AGAINST", "MEASURED_INACTIVE_AGAINST")
#: The order they are reported in, loosest last.
DEFINITIONS = ("mechanism", "targets", "mech_or_tgt", "any_edge")


def build_parser(root: Path | None = None) -> argparse.ArgumentParser:
    root = ROOT if root is None else root
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--release", type=Path, default=release_dir(root))
    ap.add_argument("--results", type=Path, default=root / "data/results",
                    help="for the independent label-count cross-check")
    ap.add_argument("--db", type=Path, default=None,
                    help="chembl sqlite, for the drug_mechanism definition. "
                         "Omitted: the mechanism columns are reported as "
                         "unavailable, never merged into another definition.")
    return ap


def mechanism_nodes(db: Path, nodes_by_chembl: dict[str, str]) -> set[str]:
    """Release node ids whose compound has a curated ChEMBL mechanism.

    Keyed through the node's own `chembl_id` property rather than by name or
    InChIKey: drug_mechanism is keyed on molregno/chembl_id and the graph is
    keyed on InChIKey, and the join has to be the one the ingest used or the
    answer is about a different set of compounds.
    """
    from kgav.host_targets import load_mechanisms
    mech_ids = {c for c, _ in load_mechanisms(db)}
    return {nid for cid, nid in nodes_by_chembl.items() if cid in mech_ids}


def read_nodes(release: Path) -> tuple[dict[str, str], dict[str, str], set[str]]:
    """(chembl_id -> node id, virus id -> label, human protein ids)."""
    by_chembl: dict[str, str] = {}
    names: dict[str, str] = {}
    human: set[str] = set()
    for line in (release / "nodes.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        n = json.loads(line)
        props = n.get("properties") or {}
        if n["class"] == "SmallMolecule":
            cid = props.get("chembl_id")
            if cid:
                by_chembl.setdefault(cid, n["id"])
        elif n["class"] == "OrganismTaxon":
            names[n["id"]] = props.get("label", n["id"])
        elif n["class"] == "Protein" and props.get("taxon_id") == "NCBITaxon:9606":
            human.add(n["id"])
    return by_chembl, names, human


def header_figures(release: Path, traversable: set[str]) -> dict:
    """Section B's own three figures, re-derived.

    It quotes "13,916 of 19,506 labelled compounds (71%) have no edge besides
    their label". Both operands are counted here from the release rather than
    carried: the subjects of an activity-label edge, and how many of those
    appear nowhere else in the graph.
    """
    labelled: set[str] = set()
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] in LABEL_PREDICATES:
            labelled.add(e["subject"])
    bare = {c for c in labelled if c not in traversable}
    return {"labelled": len(labelled), "bare": len(bare),
            "share_bare": len(bare) / len(labelled) if labelled else 0.0}


def label_sources(release: Path) -> set[str]:
    """Knowledge sources that assert a label edge in this release.

    Read rather than listed, because `build_labels(sources=...)` takes
    whatever is here and a hand-written pair would silently stop covering a
    screen the moment one is added.
    """
    from kgav.labels import edge_sources
    out: set[str] = set()
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] in LABEL_PREDICATES:
            out |= edge_sources(e)
    return out


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    schema = load_schema()
    held = schema.held_out_predicates()
    g = Graph.load(args.release, skip_predicates=held, symmetry=schema.symmetry())
    labels = build_labels(args.release)
    by_chembl, names, human = read_nodes(args.release)

    # TARGETS BY TAXON, COUNTED. "has a human target" and "has a target" are
    # the same restriction on v0.1 and that is a fact about this release, not
    # a property of the predicate -- so it is counted, and reported when it
    # stops being true.
    targets: set[str] = set()
    by_taxon = {"human": 0, "other": 0}
    for (pred, node), adj in g.out.items():
        if pred != "TARGETS":
            continue
        if adj:
            targets.add(node)
        for nb, _ in adj:
            by_taxon["human" if nb in human else "other"] += 1
    print(f"release: {args.release}")
    print(f"TARGETS edges by object taxon: {by_taxon['human']:,} human, "
          f"{by_taxon['other']:,} other")
    if by_taxon["other"]:
        print("  NOTE: `targets` and `targets_human` are no longer the same "
              "restriction; this script reports the first.")

    traversable = {c for c in by_chembl.values() if g.degree_unlabelled[c] > 0}
    traversable |= {n for n in g.degree_unlabelled if g.degree_unlabelled[n] > 0}
    hdr = header_figures(args.release, traversable)
    print(f"\nSECTION B's HEADER, re-derived: {hdr['bare']:,} of "
          f"{hdr['labelled']:,} labelled compounds ({hdr['share_bare']:.1%}) "
          f"have no edge besides their label")

    mech: set[str] | None = None
    if args.db:
        if not args.db.exists():
            print(f"\n--db {args.db} does not exist", file=sys.stderr)
            return 1
        mech = mechanism_nodes(args.db, by_chembl)
        print(f"curated-mechanism compounds in the graph: {len(mech):,}")
    else:
        print("\nmechanism columns: UNAVAILABLE (no --db); not folded into "
              "another definition")

    def admits(name: str):
        if name == "mechanism":
            return None if mech is None else (lambda c: c in mech)
        if name == "targets":
            return lambda c: c in targets
        if name == "mech_or_tgt":
            if mech is None:
                return None
            return lambda c: c in mech or c in targets
        return lambda c: g.degree_unlabelled[c] > 0

    viruses = sorted(set(labels.positives) | set(labels.negatives),
                     key=lambda v: names.get(v, v))
    active = ["all"] + [d for d in DEFINITIONS if admits(d) is not None]
    tests = {"all": (lambda c: True)} | {d: admits(d) for d in active[1:]}

    print(f"\nLABELS KEPT, per virus. pos/neg and the positive share; "
          f"* marks a cell with >= {MIN_REACHED_POSITIVES} of each")
    head = f"{'virus':13s}" + "".join(f"{d:>23s}" for d in active)
    print(head)
    print("-" * len(head))
    totals = {d: [0, 0] for d in active}
    cells = {d: 0 for d in active}
    surviving: dict[str, list[tuple[str, int, int]]] = {d: [] for d in active}
    for v in viruses:
        pos = labels.positives.get(v, set())
        neg = labels.negatives.get(v, set())
        if not pos and not neg:
            continue
        row = f"{names.get(v, v)[:12]:13s}"
        for d in active:
            f = tests[d]
            p = sum(1 for c in pos if f(c))
            n = sum(1 for c in neg if f(c))
            totals[d][0] += p
            totals[d][1] += n
            ok = p >= MIN_REACHED_POSITIVES and n >= MIN_REACHED_POSITIVES
            cells[d] += ok
            if ok:
                surviving[d].append((names.get(v, v), p, n))
            rate = p / (p + n) if p + n else 0.0
            row += f"{p:6,}/{n:6,}{'*' if ok else ' '}{rate:6.1%} "
        print(row)
    print("-" * len(head))
    row = f"{'TOTAL':13s}"
    for d in active:
        p, n = totals[d]
        row += f"{p:6,}/{n:6,} {p / (p + n):6.1%} " if p + n else f"{'-':>22s} "
    print(row)

    print("\nWHAT EACH DEFINITION COSTS. `ratio` is kept-positives over "
          "kept-negatives:\n  1.00x is label-neutral, above it the "
          "restriction prefers positives, below it negatives.")
    p0, n0 = totals["all"]
    print(f"  {'definition':13s} {'cells':>6s} {'pos kept':>10s} "
          f"{'neg kept':>10s} {'ratio':>8s}")
    for d in active:
        p, n = totals[d]
        kp, kn = (p / p0 if p0 else 0.0), (n / n0 if n0 else 0.0)
        ratio = f"{kp / kn:.2f}x" if kn else "n/a"
        print(f"  {d:13s} {cells[d]:6d} {kp:10.1%} {kn:10.1%} {ratio:>8s}")

    # THE SIGN, PER VIRUS, COUNTED. The pooled ratio above is one number and
    # it hides the thing that makes this restriction dangerous: it enriches
    # positives in most viruses and negatives in the one holding most of the
    # labels, so a pooled reading and a per-virus reading move in opposite
    # directions. A cell with no admitted label at all has no sign and is
    # reported as such rather than counted either way.
    print("\nTHE SIGN IS NOT CONSTANT ACROSS VIRUSES. Positive share under "
          "each definition\nagainst the same virus's share on the full task:")
    for d in active[1:]:
        up = down = flat = undef = 0
        f = tests[d]
        for v in viruses:
            pos = labels.positives.get(v, set())
            neg = labels.negatives.get(v, set())
            if not (pos or neg):
                continue
            p = sum(1 for c in pos if f(c))
            n = sum(1 for c in neg if f(c))
            if p + n == 0:
                undef += 1
                continue
            before = len(pos) / (len(pos) + len(neg))
            after = p / (p + n)
            if after > before:
                up += 1
            elif after < before:
                down += 1
            else:
                flat += 1
        print(f"  {d:13s} raises it for {up:2d} virus(es), lowers it for "
              f"{down:2d}, unchanged {flat:2d}, no admitted label {undef:2d}")

    print("\nWHAT THE SURVIVING CELLS COULD DETECT. The AUC each would need "
          "before its\ninterval excludes chance at all -- a task that cannot "
          "reject is not a task.")
    for d in active[1:]:
        if not surviving[d]:
            print(f"  {d:13s} no cell survives")
            continue
        for name, p, n in surviving[d]:
            need = 0.5
            while need < 1.0:
                lo, _ = auc_interval(need, p, n)
                if lo > 0.5:
                    break
                need += 0.001
            lo, hi = auc_interval(0.5, p, n)
            print(f"  {d:13s} {name:12s} {p:5,}/{n:6,}  "
                  f"chance interval [{lo:.3f}, {hi:.3f}]  needs AUC >= {need:.3f}")

    # ---- THE ONE CONSTRUCTION THAT CAN BE LABEL-NEUTRAL, and what it costs.
    #
    # Pooling screens is what makes the restriction select on the label: the
    # positives are ChEMBL antiviral research compounds and the panel's
    # negatives are library compounds, and the two populations carry
    # annotation at different rates for reasons that have nothing to do with
    # antiviral activity. build_labels(sources=...) exists to remove that by
    # construction -- both classes off the same plates under one protocol --
    # so the restriction is measured again WITHIN each screen. A screen
    # covering one virus yields one cell, which is the cost.
    print("\nWITHIN A SINGLE SCREEN, which is the only construction that "
          "matches the two\nclasses by design rather than by a filter "
          "(build_labels(sources=...)):")
    for screen in sorted(label_sources(args.release)):
        one = build_labels(args.release, sources={screen})
        covered = [v for v in viruses
                   if one.positives.get(v) or one.negatives.get(v)]
        print(f"  {screen}  ({len(covered)} virus(es) covered)")
        for v in covered:
            pos = one.positives.get(v, set())
            neg = one.negatives.get(v, set())
            for d in active[1:]:
                f = tests[d]
                p = sum(1 for c in pos if f(c))
                n = sum(1 for c in neg if f(c))
                if not (pos and neg):
                    continue
                kp, kn = p / len(pos), n / len(neg)
                ratio = f"{kp / kn:.2f}x" if kn else "n/a"
                mark = ("*" if p >= MIN_REACHED_POSITIVES
                        and n >= MIN_REACHED_POSITIVES else " ")
                need = 0.5
                while need < 1.0 and p and n:
                    lo, _ = auc_interval(need, p, n)
                    if lo > 0.5:
                        break
                    need += 0.001
                print(f"    {names.get(v, v)[:12]:13s} {d:12s} "
                      f"{p:5,}/{n:6,}{mark} kept {kp:5.1%}/{kn:5.1%} "
                      f"ratio {ratio:>6s}"
                      + (f"  needs AUC >= {need:.3f}" if mark == "*" else ""))

    # ---- CROSS-CHECK against figures this script did not compute.
    print("\ncross-check: per-virus label counts against the n_positives / "
          "n_negatives\nindependently recorded in each baselines result")
    checked = mismatched = 0
    for f in sorted(args.results.glob("baselines_NCBITaxon_*.json")):
        d = json.loads(f.read_text())
        v = d.get("virus")
        if v is None:
            continue
        checked += 1
        p = len(labels.positives.get(v, set()))
        n = len(labels.negatives.get(v, set()))
        if (p, n) != (d["n_positives"], d["n_negatives"]):
            mismatched += 1
            print(f"  MISMATCH {v}: here {p:,}/{n:,}, "
                  f"result {d['n_positives']:,}/{d['n_negatives']:,}")
    print(f"  {checked} result file(s), {mismatched} mismatch(es)")
    return 1 if mismatched else 0


if __name__ == "__main__":
    raise SystemExit(main())
