#!/usr/bin/env python3
"""Run the three non-learned baselines and the biology controls.

The metric is not the point. The point is whether known antivirals surface and
whether known artifacts do not.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav import calibration as C
from kgav.baselines import (
    Graph,
    combine,
    degree_ranking,
    dwpc_scores,
    hits_at_k,
    mrr,
    network_proximity,
    rank,
)
from kgav.labels import build_labels
from kgav.provenance import write_results
from kgav.schema import load_schema

# Compounds whose behaviour tells us whether the graph works.
#   nirmatrelvir  nsp5 inhibitor  -> must rank high via M1
#   remdesivir    nsp12 inhibitor -> must rank high via M1
#   hydroxychloroquine / chloroquine: abundant SARS-CoV-2 activity generated
#   largely in Vero E6, where entry is cathepsin-dependent. Both failed in
#   TMPRSS2-expressing airway cells. If these rank high the graph is
#   reproducing a cell-line artifact, and with no CC50 data we cannot filter
#   them on selectivity.
CONTROLS = {
    "nirmatrelvir": "positive",
    "remdesivir": "positive",
    "GS-441524": "positive",
    "hydroxychloroquine": "artifact",
    "chloroquine": "artifact",
}


def name_index(g: Graph, db: Path | None) -> dict[str, str]:
    """chembl_id -> node id, plus a lowercase pref_name lookup if the ChEMBL
    database is available."""
    by_chembl = {n["properties"].get("chembl_id"): nid
                 for nid, n in g.nodes.items()
                 if n["class"] == "SmallMolecule" and n["properties"].get("chembl_id")}
    out: dict[str, str] = {}
    if db and db.exists():
        import sqlite3
        con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        for cid, pref in con.execute(
                "SELECT chembl_id, lower(pref_name) FROM molecule_dictionary "
                "WHERE pref_name IS NOT NULL"):
            if cid in by_chembl:
                out[pref] = by_chembl[cid]
        con.close()
    return out


# WHETHER THE DWPC COLUMN IS A VALIDATED RANKING IS PART OF THE RESULT, so
# the two decisions behind it live here as plain functions rather than inside
# main(): a script whose behaviour can only be checked by running it against
# a 74,872-node release is a script whose behaviour is not checked.
DWPC_MODES = ("uncalibrated", "calibrated", "refused")


def virus_labels(g) -> dict[str, str]:
    """taxon CURIE -> the display label hard_negatives keys its output by.

    Derived exactly as hard_negatives.run derives it, from the node's own
    `label` property, because that is the key its results file holds. A
    second rule for producing the same key would be a second chance to get
    it wrong.
    """
    return {v: ((g.nodes.get(v, {}).get("properties") or {}).get("label") or v)
            for v in g.viruses()}


def load_calibration(calibrated: bool, path: Path, virus: str, section: str,
                     node_labels: dict[str, str] | None = None):
    """The Calibration to score with, or None to score uncalibrated.

    Refuses rather than falling back, in two different ways.

    A MISSING FILE is fatal: asking for calibration and being handed an
    uncalibrated order under the same column heading is the exact failure the
    flag exists to prevent, and the file describes one specific evaluation
    run with no default to stand in for it.

    A KEY MISS is fatal too, and this is the one the first version of this
    function got wrong. hard_negatives keys its results by DISPLAY LABEL,
    not by taxon CURIE, and calibration.load says plainly that a CURIE
    "cannot be translated by string surgery" so the lookup "simply misses"
    unless node_labels is supplied -- available_viruses() existing "so the
    miss can be reported with what the file does hold, rather than as an
    empty calibration that reads like an empty evaluation."

    Passing the CURIE with no node_labels did exactly that. The calibration
    came back with ZERO channels, can_rank() was false because there was
    nothing to rank with, and run_baselines printed

        DWPC RANKING REFUSED. NO channel was evaluable for NCBITaxon:2697049
        ... every one reaches fewer than 10 of the measured actives
        ... Not evaluated: .

    -- a sentence about every channel having been measured and failed, with
    an empty list where the channels should be, from a file that had never
    been read for that key. The refusal happened to agree with the truth on
    this release, which is what made it look right. So node_labels is
    supplied, and an empty channel set is now a reported miss rather than a
    verdict.
    """
    if not calibrated:
        return None
    if not path.exists():
        raise FileNotFoundError(
            f"--calibrated needs {path}, which does not exist. Run "
            f"scripts/hard_negatives.py first; this file describes one "
            f"specific evaluation run and there is no default to fall back on.")
    cal = C.load(path, virus, section=section, node_labels=node_labels)
    if not cal.channels and not cal.refused_reason:
        have = C.available_viruses(path, section)
        raise LookupError(
            f"no channels recorded for {virus!r} in {section} of "
            f"{path.name}.\n  that file holds: "
            f"{', '.join(have) if have else '(nothing)'}\n  protocol it "
            f"records: {cal.protocol}\n  This is a KEY MISS, not a refusal: "
            f"nothing was measured for this virus in this section, which is "
            f"a different statement from every channel having been measured "
            f"and failed.")
    return cal


def dwpc_mode(cal) -> str:
    """Which of the three claims this run makes about the DWPC column.

    A calibrated ranking, an uncalibrated one and a refusal are three
    different statements. combine's own comment is why the uncalibrated case
    needs a name at all: max-over-channels is biased upward by the channel
    count, so a drug at the top of a CHANCE-LEVEL pool scores 1.0 and
    outranks one with a real route.
    """
    if cal is None:
        return "uncalibrated"
    return "calibrated" if cal.can_rank() else "refused"


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--db", type=Path,
                    default=root / "data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db")
    ap.add_argument("--virus", default="NCBITaxon:2697049")
    ap.add_argument("--out", type=Path, default=root / "data/results")
    ap.add_argument("--top", type=int, default=25)
    # calibrate.py has been printing "run_baselines.py  # with --calibrated
    # once wired" since the calibration module landed. This is that wiring.
    ap.add_argument("--calibrated", action="store_true",
                    help="restrict the DWPC ranking to channels whose measured "
                         "interval excludes chance, and refuse to rank if none do")
    ap.add_argument("--calibration", type=Path,
                    default=root / "data/results/hard_negatives.json",
                    help="the evaluation run --calibrated reads its channels from")
    ap.add_argument("--section", default="cross_sectional",
                    choices=["cross_sectional", "temporal",
                             "cross-sectional", "post-2021"])
    args = ap.parse_args()

    schema = load_schema()
    held = schema.held_out_predicates()
    print(f"held out of traversal (evaluation labels): {sorted(held)}")

    g = Graph.load(args.release, skip_predicates=held, symmetry=schema.symmetry())
    print(f"{len(g.nodes):,} nodes | {len(g.drugs()):,} drugs | "
          f"{len(g.viruses())} viruses")

    # Ground truth: the cell-based antiviral edges we withheld, with polarity
    # decided by labels.classify rather than by the relation alone. Reading
    # only the relation counted 578 compounds measured above 10 uM as
    # positives -- they are measurements of INACTIVITY -- and a positive set
    # containing measured inactives makes Hits@k and MRR uninterpretable.
    labels = build_labels(args.release)
    positives = labels.positives.get(args.virus, set())
    negatives = labels.negatives.get(args.virus, set())
    print(f"{len(positives):,} compounds measured ACTIVE against {args.virus} "
          f"(ground truth), {len(negatives):,} measured inactive")

    print("\ncomputing DWPC over M1-M6")
    dwpc = dwpc_scores(g, schema.metapaths, schema.retrieval_limits)
    for name in sorted(dwpc):
        reach = sum(1 for by_virus in dwpc[name].values() if args.virus in by_virus)
        print(f"  {name}: {reach:,} drugs reach this virus")

    # WHETHER THE DWPC COLUMN IS A VALIDATED RANKING, SAID OUT LOUD.
    #
    # This script used to call combine(dwpc, args.virus) and print the result
    # in a column headed "dwpc", with nothing to distinguish it from a
    # measured ranking. combine's own comment explains why that is not a
    # neutral default: max-over-channels is biased upward by the channel
    # count, so a drug sitting at the top of a CHANCE-LEVEL pool scores 1.0
    # and outranks one with a real route. Reading the column without knowing
    # which channels cleared chance is reading noise as an order.
    #
    # --calibrated restricts the max to channels whose measured 95% interval
    # excludes 0.5 and refuses to rank when none do, which on v0.1 is every
    # one of them. Refusing is the result in that case, not a failure to
    # produce one. Without the flag the uncalibrated ranking still prints --
    # every published number reproduces -- but it is now LABELLED, and the
    # calibration file is consulted so the label can say what calibration
    # would have said.
    try:
        cal = load_calibration(args.calibrated, args.calibration,
                               args.virus, args.section, virus_labels(g))
    except (FileNotFoundError, LookupError) as exc:
        sys.exit(str(exc))
    mode = dwpc_mode(cal)

    scorers = {"degree": degree_ranking(g)}
    if mode == "uncalibrated":
        scorers["dwpc"] = combine(dwpc, args.virus)
        print("\n  dwpc column is UNCALIBRATED: the max runs over every "
              "channel that\n  reaches this virus, whether or not it was "
              "measured above chance.")
        if args.calibration.exists():
            peek = C.load(args.calibration, args.virus, section=args.section,
                          node_labels=virus_labels(g))
            if not peek.channels:
                print(f"  calibration on {args.calibration.name} records "
                      f"nothing for this virus in {args.section}, so it says "
                      f"nothing either way.")
            elif not peek.can_rank():
                print(f"  calibration on {args.calibration.name} would refuse "
                      f"to rank: {peek.why_not()}")
            else:
                print(f"  calibration would keep: "
                      f"{', '.join(sorted(peek.weights()))}  (--calibrated)")
    elif mode == "calibrated":
        scorers["dwpc"] = combine(dwpc, args.virus, calibration=cal)
        print(f"\n  dwpc column is CALIBRATED on {args.section}, over "
              f"{', '.join(sorted(cal.weights()))}")
    else:
        print(f"\n  DWPC RANKING REFUSED. {cal.why_not()}")
        print("  No channel's measured interval excludes chance, so there is "
              "no order to\n  report. The degree baseline and the biology "
              "controls below still stand:\n  evidence for a named compound "
              "does not require a ranking.")
    print("\ncomputing network proximity")
    prox = network_proximity(g, args.virus, schema.retrieval_limits)
    if prox:
        scorers["proximity"] = prox
    print(f"  {len(prox):,} drugs scored")

    print(f"\n{'scorer':<12} {'n':>7} {'Hits@10':>8} {'Hits@50':>8} "
          f"{'Hits@100':>9} {'MRR':>7}")
    ranked_all = {}
    for name, scores in scorers.items():
        r = rank(scores)
        ranked_all[name] = r
        print(f"{name:<12} {len(r):>7,} {hits_at_k(r, positives, 10):>8} "
              f"{hits_at_k(r, positives, 50):>8} {hits_at_k(r, positives, 100):>9} "
              f"{mrr(r, positives):>7.4f}")

    names = name_index(g, args.db)
    if names:
        print("\nbiology controls (rank out of n; lower is higher-ranked):")
        print(f"{'compound':<22} {'expect':<10} " +
              " ".join(f"{s:>12}" for s in ranked_all))
        for cname, kind in CONTROLS.items():
            nid = names.get(cname)
            if not nid:
                print(f"{cname:<22} {kind:<10} " +
                      " ".join(f"{'not in graph':>12}" for _ in ranked_all))
                continue
            cells = []
            for r in ranked_all.values():
                pos = next((i for i, (d, _v) in enumerate(r, 1) if d == nid), None)
                cells.append(f"{pos:>12,}" if pos else f"{'unscored':>12}")
            print(f"{cname:<22} {kind:<10} " + " ".join(cells))

    # ranked_all has no "dwpc" key when the ranking was refused, and
    # indexing it there would end the run with a KeyError after the controls
    # table -- losing the results file for a run that worked exactly as
    # intended.
    if "dwpc" not in ranked_all:
        print("\nno top-N by DWPC: the ranking was refused above.")
    else:
        print(f"\ntop {args.top} by DWPC:")

    for i, (drug, score) in enumerate(ranked_all.get("dwpc", [])[:args.top], 1):
        contrib = sorted(((n, dwpc[n][drug][args.virus]) for n in dwpc
                          if drug in dwpc[n] and args.virus in dwpc[n][drug]),
                         key=lambda kv: -kv[1])
        routes = ",".join(f"{n}:{v:.3f}" for n, v in contrib[:3])
        label = next((k for k, v in names.items() if v == drug), None)
        mark = " *" if drug in positives else "  "
        print(f"{i:>3}{mark} {(label or drug)[:34]:<36} {score:.4f}  {routes}")
    # Inside the branch: on the refusal path there is no starred list for the
    # legend to describe, and printing it under "no top-N by DWPC" claims a
    # table that is not there.
    if "dwpc" in ranked_all:
        print("  * = has measured activity (was held out of traversal)")

    args.out.mkdir(parents=True, exist_ok=True)
    # HOW THE DWPC COLUMN WAS PRODUCED IS PART OF THE RESULT. A calibrated
    # ranking, an uncalibrated one and a refusal are three different claims,
    # and a file recording only the numbers cannot tell them apart.
    # check_results tracks staleness against the release and the schema, not
    # against the flags a script was invoked with, so the flags go in here.
    payload = {
        "virus": args.virus,
        "dwpc_ranking": mode,
        "dwpc_channels": (sorted(cal.weights()) if mode == "calibrated"
                          else None),
        "calibration_section": args.section if cal is not None else None,
        "calibration_protocol": cal.protocol if cal is not None else None,
        "held_out": sorted(held),
        "n_positives": len(positives),
        "n_negatives": len(negatives),
        "metrics": {n: {"n": len(r),
                        "hits@10": hits_at_k(r, positives, 10),
                        "hits@50": hits_at_k(r, positives, 50),
                        "hits@100": hits_at_k(r, positives, 100),
                        "mrr": mrr(r, positives)}
                    for n, r in ranked_all.items()},
        "top": {n: [{"drug": d, "score": s} for d, s in r[:200]]
                for n, r in ranked_all.items()},
    }
    write_results(args.out / f"baselines_{args.virus.replace(':', '_')}.json",
                  payload, args.release, schema.version)
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
