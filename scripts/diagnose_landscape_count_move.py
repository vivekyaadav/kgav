#!/usr/bin/env python3
"""Which cells of the pathway landscape moved between two releases, and why.

    python scripts/diagnose_landscape_count_move.py
    python scripts/diagnose_landscape_count_move.py --out data/results/x.json

READ-ONLY, AND NOT A PROTOCOL CHANGE. It re-reads two recorded
artifacts and the release on disk. It changes no registered figure and
decides nothing; it closes an open question by measuring it.

THE OPEN QUESTION. doc2/registrations/prereg_pathway_overlap.md's note of 2026-11-17
records that `n_observed` moved in two cells between the pre-entry-assay
release and the one on disk, names a cause for one of them, and leaves
the other open -- "it is open against the entry-assay ingest rather
than against the clinical one" -- with a controlled re-run named as the
way to close it. The re-run is not needed. Both artifacts record the
drugged proteome's MEMBERS, so the operand's change is a set difference,
and the mechanism is arithmetic on sets already on disk.

THE MECHANISM, which the note did not state. `overlap_row` computes

    n_observed = len(targets & factors)

where `targets` is the drugged proteome and `factors` is the cell's
dependency or interactome set. A cell can therefore move WITHOUT ANY
EDGE IN ITS OWN CHANNEL CHANGING: it is enough that a protein already
inside `factors` joined `targets`. The note instead looked for a new
`PHYSICALLY_INTERACTS_WITH` partner, which is why one cell came back
explained and one did not.

TWO ERRORS IN THAT NOTE, both of which this script's output contradicts
rather than asserts:

  * `NCBITaxon:31631` is **HCoV-OC43**, not SARS-CoV (config/viruses.yaml).
    doc2/plan.md's open-work table named the wrong virus.
  * the cell is in the **dependency** channel, not the interactome one.
    OC43's interactome cell is 28 and did not move. So the note's search
    through PPI partners was looking in a channel the cell never used,
    and no answer it returned could have been the right one.

WHAT WOULD FALSIFY THE EXPLANATION, stated before it is run: a cell that
moved but whose predicted delta is zero, or a cell whose predicted delta
is nonzero but which did not move. Either means the operand difference
is not the whole cause and something in a channel did change. The script
reports both classes and `explains_every_move` is false if either is
non-empty.

A SECOND, INDEPENDENT CHECK on the same claim. If no channel edge
changed, then every cell's `n_dependency` -- the size of `factors`, which
this analysis stores -- is identical across the two artifacts. That is
checked for all cells, not just the two that moved, and reported as
`factor_sets_unchanged`.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.baselines import Graph  # noqa: E402
from kgav.pipeline import release_dir  # noqa: E402
from kgav.provenance import write_results  # noqa: E402
from kgav.schema import load_schema, print_capped  # noqa: E402

import pathway_landscape as PL  # noqa: E402

#: The artifact this analysis writes and re-reads.
ARTIFACT = "data/results/pathway_landscape.json"

#: The commit holding the PRE-ENTRY-ASSAY artifact. This is a historical
#: fact, not a tunable: it is the last revision of ARTIFACT measured on
#: the release before the three NCATS pseudotyped screens landed. It is
#: not trusted blindly -- the script requires the two artifacts to
#: disagree about the drugged proteome's sha256, so pointing this at the
#: wrong revision produces a refusal rather than a wrong answer.
BEFORE_REV = "57602d9"

#: The two channels the landscape reports, and the artifact key holding
#: each. Never pooled -- see interactome_partners' docstring.
CHANNELS = (("dependency", "taxa"), ("interactome", "interactome"))


class MoveError(RuntimeError):
    """The inputs cannot answer the question asked of them."""


def read_artifact(rev: str | None, path: Path) -> dict:
    """The artifact as of a git revision, or from the working tree."""
    if rev is None:
        return json.loads(path.read_text())
    rel = path.relative_to(ROOT) if path.is_absolute() else path
    out = subprocess.run(["git", "show", f"{rev}:{rel.as_posix()}"],
                         cwd=ROOT, capture_output=True, text=True)
    if out.returncode != 0:
        raise MoveError(f"cannot read {rel} at {rev}: "
                        f"{out.stderr.strip() or 'git show failed'}")
    return json.loads(out.stdout)


def drugged(artifact: dict) -> tuple[set[str], str, int]:
    """(members, sha256, n) of an artifact's drugged proteome operand.

    The members list is what makes the re-run unnecessary, so its absence
    is an error with a name rather than an empty set that quietly
    explains nothing.
    """
    op = ((artifact.get("operands") or {}).get("drugged_proteome") or {})
    members = op.get("members")
    if not members:
        raise MoveError("this artifact records no drugged_proteome.members, "
                        "so the operand difference cannot be derived from "
                        "it; only a controlled re-run could answer")
    return set(members), str(op.get("sha256") or ""), int(op.get("n") or 0)


def factor_sets(g: Graph, schema) -> dict[tuple[str, str], set[str]]:
    """{(channel, taxon): factors} exactly as pathway_landscape builds them.

    IMPORTED, NOT REIMPLEMENTED. `dependency_factors` and
    `interactome_partners` are the functions that produced the numbers
    being explained; a copy of either could agree with itself and
    disagree with the analysis.
    """
    hub_cut = int((schema.retrieval_limits or {}).get("hub_degree_cut") or 0)
    vprot = PL.viral_proteins(g)
    out: dict[tuple[str, str], set[str]] = {}
    for taxon in sorted(vprot):
        partners, _hubs = PL.interactome_partners(g, vprot[taxon],
                                                  hub_cut or None)
        out[("interactome", taxon)] = partners
    for taxon in sorted({t for t in vprot}):
        out[("dependency", taxon)] = PL.dependency_factors(g, taxon)
    return out


def compare(before: dict, after: dict,
            factors: dict[tuple[str, str], set[str]]) -> dict:
    """Observed moves against the moves the operand difference predicts."""
    b_members, b_sha, b_n = drugged(before)
    a_members, a_sha, a_n = drugged(after)
    if b_sha and a_sha and b_sha == a_sha:
        raise MoveError(
            "both artifacts record the same drugged proteome "
            f"(sha256 {a_sha[:16]}), so there is no operand change to "
            "explain -- check --before-rev points at the earlier release")
    joined = sorted(a_members - b_members)
    left = sorted(b_members - a_members)

    cells: dict[str, dict] = {}
    for channel, key in CHANNELS:
        b_block, a_block = before.get(key) or {}, after.get(key) or {}
        for taxon in sorted(set(b_block) | set(a_block)):
            b_row, a_row = b_block.get(taxon), a_block.get(taxon)
            fac = factors.get((channel, taxon), set())
            pred = (len([p for p in joined if p in fac])
                    - len([p for p in left if p in fac]))
            obs = (None if (b_row is None or a_row is None)
                   else a_row["n_observed"] - b_row["n_observed"])
            cells[f"{channel}|{taxon}"] = {
                "channel": channel, "taxon": taxon,
                "n_observed_before": b_row and b_row.get("n_observed"),
                "n_observed_after": a_row and a_row.get("n_observed"),
                "observed_delta": obs,
                "predicted_delta": pred,
                "n_factors": len(fac),
                "n_dependency_before": b_row and b_row.get("n_dependency"),
                "n_dependency_after": a_row and a_row.get("n_dependency"),
                "joined_inside": [p for p in joined if p in fac],
                "left_inside": [p for p in left if p in fac],
            }

    moved = {k: c for k, c in cells.items() if c["observed_delta"]}
    predicted = {k: c for k, c in cells.items() if c["predicted_delta"]}
    unexplained = sorted(k for k, c in moved.items()
                         if c["predicted_delta"] != c["observed_delta"])
    overpredicted = sorted(k for k, c in predicted.items()
                           if c["predicted_delta"] != (c["observed_delta"] or 0))
    changed_factors = sorted(
        k for k, c in cells.items()
        if c["n_dependency_before"] is not None
        and c["n_dependency_after"] is not None
        and c["n_dependency_before"] != c["n_dependency_after"])

    return {
        "drugged_proteome": {
            "before": {"n": b_n, "sha256": b_sha},
            "after": {"n": a_n, "sha256": a_sha},
            "joined": joined, "left": left,
        },
        "cells": cells,
        "n_cells": len(cells),
        "moved": sorted(moved),
        "predicted_to_move": sorted(predicted),
        "unexplained_moves": unexplained,
        "moves_predicted_that_did_not_happen": overpredicted,
        "explains_every_move": not unexplained and not overpredicted,
        "cells_with_changed_factor_set": changed_factors,
        "factor_sets_unchanged": not changed_factors,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release", type=Path, default=release_dir(ROOT))
    ap.add_argument("--artifact", type=Path, default=ROOT / ARTIFACT)
    ap.add_argument("--before-rev", default=BEFORE_REV,
                    help=f"git revision holding the earlier artifact "
                         f"(default {BEFORE_REV})")
    ap.add_argument("--before", type=Path, default=None,
                    help="read the earlier artifact from this path instead "
                         "of from a git revision")
    ap.add_argument("--out", type=Path, default=None,
                    help="write a STAMPED results file here")
    args = ap.parse_args(argv)

    s = load_schema()
    after = read_artifact(None, args.artifact)
    if args.before is not None:
        before = read_artifact(None, args.before)
        before_from = str(args.before)
    else:
        before = read_artifact(args.before_rev, args.artifact)
        before_from = f"{args.before_rev}:{ARTIFACT}"

    g = Graph.load(args.release, skip_predicates=s.held_out_predicates(),
                   symmetry=s.symmetry())
    try:
        r = compare(before, after, factor_sets(g, s))
    except MoveError as exc:
        print(f"REFUSED: {exc}")
        return 2

    dp = r["drugged_proteome"]
    print(f"release {args.release.name}   schema {s.version}")
    print(f"  before  {before_from}")
    print(f"          measured_on {before.get('measured_on')}   "
          f"drugged proteome {dp['before']['n']:,}  "
          f"sha {dp['before']['sha256'][:16]}")
    print(f"  after   {args.artifact}")
    print(f"          measured_on {after.get('measured_on')}   "
          f"drugged proteome {dp['after']['n']:,}  "
          f"sha {dp['after']['sha256'][:16]}")
    print(f"\nTHE OPERAND DIFFERENCE: {len(dp['joined'])} joined, "
          f"{len(dp['left'])} left")
    print_capped(dp["joined"], 20, unit="more joined")
    print_capped(dp["left"], 20, unit="more left")

    print("\n" + "=" * 78)
    print("EVERY CELL. 'pred' is the delta the operand difference alone")
    print("predicts: the joined proteins already inside that cell's factor")
    print("set, minus the ones that left. 'obs' is what the two artifacts")
    print("record. They must agree in BOTH directions.")
    print("=" * 78)
    print(f"  {'channel':12} {'taxon':20} {'factors':>8} {'before':>7} "
          f"{'after':>6} {'obs':>5} {'pred':>5}  why")
    print("  " + "-" * 84)
    for key in sorted(r["cells"]):
        c = r["cells"][key]
        why = ", ".join(p.split(":")[-1] for p in c["joined_inside"]) or "-"
        if c["left_inside"]:
            why += " / -" + ", -".join(p.split(":")[-1]
                                       for p in c["left_inside"])
        flag = "  <--" if c["observed_delta"] else ""
        print(f"  {c['channel']:12} {c['taxon']:20} {c['n_factors']:8,} "
              f"{_n(c['n_observed_before']):>7} {_n(c['n_observed_after']):>6} "
              f"{_d(c['observed_delta']):>5} {_d(c['predicted_delta']):>5}  "
              f"{why}{flag}")

    print(f"\n{len(r['moved'])} of {r['n_cells']} cell(s) moved; "
          f"{len(r['predicted_to_move'])} predicted to move.")
    for key in r["moved"]:
        c = r["cells"][key]
        print(f"    {c['channel']} {c['taxon']}: "
              f"{c['n_observed_before']} -> {c['n_observed_after']} "
              f"via {', '.join(c['joined_inside']) or 'NOTHING IN THE OPERAND'}")

    print("\nFALSIFIERS, both of which must be empty:")
    print(f"    moves with no operand cause:      "
          f"{r['unexplained_moves'] or 'none'}")
    print(f"    predicted moves that did not happen: "
          f"{r['moves_predicted_that_did_not_happen'] or 'none'}")
    print("\nINDEPENDENT CHECK -- no channel edge changed, so every cell's")
    print("factor set should be the same size in both artifacts:")
    print(f"    cells with a changed factor set:  "
          f"{r['cells_with_changed_factor_set'] or 'none'}")

    ok = r["explains_every_move"] and r["factor_sets_unchanged"]
    print("\n" + ("THE OPERAND DIFFERENCE EXPLAINS EVERY MOVE, and no "
                  "channel\nedge changed. No controlled re-run is needed."
                  if ok else
                  "NOT FULLY EXPLAINED. Something in a channel changed; the "
                  "cells\nnamed above need the controlled re-run after all."))

    if args.out:
        payload = dict(r)
        payload |= {
            "measured_by": "scripts/diagnose_landscape_count_move.py",
            "question_from": "doc2/registrations/prereg_pathway_overlap.md",
            "before_artifact": before_from,
            "after_artifact": str(args.artifact),
            "before_measured_on": before.get("measured_on"),
            "after_measured_on": after.get("measured_on"),
            "mechanism": "n_observed = len(targets & factors); a cell moves "
                         "when a protein already inside factors joins the "
                         "drugged proteome, with no edge in its own channel "
                         "changing",
        }
        write_results(args.out, payload, args.release, s.version)
        print(f"\nwrote {args.out}")

    print("\nNOT A PROTOCOL CHANGE. Every registered figure stands; this "
          "names the\ncause of two that moved.")
    return 0 if ok else 1


def _n(v) -> str:
    return "--" if v is None else f"{v:,}"


def _d(v) -> str:
    return "--" if v is None else f"{v:+d}" if v else "0"


if __name__ == "__main__":
    raise SystemExit(main())
