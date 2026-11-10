#!/usr/bin/env python3
"""Run the whole ingest chain in dependency order and stop on a real failure.

    make rebuild              # the chain
    make rebuild-plan         # what it would run, no execution

WHY NOT A LIST OF COMMANDS IN THE MAKEFILE. Three layers read the ASSEMBLED
release rather than the per-layer ones, so assemble_release has to run in the
middle and again at the end. src/kgav/pipeline.py holds that order as data and
this script executes it; neither restates the other.

WHAT COUNTS AS A FAILURE. A layer whose raw inputs are absent is SKIPPED and
named -- data/raw is gitignored and populated by hand, so an absent input is
an ordinary state of a checkout. A layer whose inputs are present and which
then exits non-zero STOPS the run, because the next assemble would otherwise
write a release missing that layer while every count still looked right.

spine and host are required: without them there is no graph, and the run
refuses to start rather than assembling an empty one.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.pipeline import blocking, layer_dir, plan, release_dir
from kgav.provenance import release_fingerprint
from kgav.similarity import annotating_subjects_sha
from kgav.provenance import write_layer


# THE CLOSING LINE USED TO ASSERT STALENESS RATHER THAN MEASURE IT (L13).
# Every run ended "Results in data/results were computed from the previous
# release and are now stale", which is false of the one run a reader most
# needs it right about: a rebuild that changed nothing. The 2026-11-10
# re-stamp rebuild came back byte-identical -- edges_sha f8a1ae20bfd89b37
# unmoved -- and check-results reported 27 current, 0 stale.
#
# A pure function because the alternative is a sentence only a 7-minute
# rebuild can exercise, and that is how the wrong sentence survived.
#
# NOTE THE SCOPE, WHICH IS NARROWER THAN IT LOOKS. An unchanged release
# proves this RUN staled nothing. It does not prove the results were current
# beforehand -- they may have been stale already, from an earlier run. Only
# check_results.py compares each file's own stamp, so the pointer to it
# stays in all three branches rather than being presented as redundant.
def staleness_verdict(before: dict | None, after: dict | None,
                     annotating_before: str | None = None,
                     annotating_after: str | None = None) -> str:
    """What this run did to data/results, measured from the two fingerprints.

    `before` and `after` are release_fingerprint() dicts taken around the
    run, or None/all-None when no release was on disk.

    THE FIFTH BRANCH IS L20's (option (c), taken 2026-11-11). The four
    original branches cover no release / none after / unchanged / changed,
    and "changed because it has not converged yet" reads to any operator as
    "changed because the inputs did". It is not a hypothetical: reverting B1
    went af01e12ce8a86891 -> 8402bbaf9ed5c531, differing from the intended
    release by exactly one CHEMICALLY_SIMILAR_TO edge, and a second run went
    -> f8a1ae20bfd89b37 and stopped. Both printed "The release CHANGED this
    run".

    `annotating_before`/`annotating_after` are similarity.
    annotating_subjects_sha() taken around the run. When they differ, what
    `ingest_similarity` reads at phase 3 is not what it read this run, so
    phase 3 would produce something different if run again and this release
    is NOT the fixed point.

    WHAT THIS DOES NOT DO. It does not iterate, and it does not rebuild the
    phase-3 layers from the phase-1 directories; those are L20's options (a)
    and (b) and both are still open. This branch only stops the run from
    reporting convergence it has not measured. On this repository's history
    nine of ten rebuilds would not have triggered it.
    """
    def shas(fp):
        if not fp:
            return None
        pair = (fp.get("nodes_sha"), fp.get("edges_sha"))
        return None if pair == (None, None) else pair

    b, a = shas(before), shas(after)
    if b is None:
        return ("There was no release on disk before this run, so nothing in "
                "data/results\ncan have been computed from it. "
                "`make check-results` reports what each file\nwas stamped "
                "against.")
    if a is None:
        return ("This run left no release on disk, which should not happen "
                "after a\nsuccessful rebuild. Do not trust data/results "
                "until `make check-results`\nsays what they match.")
    if a == b:
        return (f"This run did not change the release (edges_sha "
                f"{after['edges_sha']}), so it\nmade nothing in "
                "data/results stale. Whether those files were already\n"
                "current is a separate question: `make check-results`.")
    changed = (f"The release CHANGED this run (edges_sha "
               f"{before['edges_sha']} -> {after['edges_sha']}),\nso every "
               f"result computed from the previous one is now stale. "
               f"`make check-results`\nreports which.")
    if (annotating_before is not None and annotating_after is not None
            and annotating_before != annotating_after):
        return (changed + "\n\nAND THIS RELEASE MAY NOT BE THE FIXED POINT "
                "-- RUN THE REBUILD AGAIN.\nThe set of compounds carrying "
                f"TARGETS / INHIBITS / MEASURED_INACTIVE_AGAINST\nmoved this "
                f"run ({annotating_before} -> {annotating_after}). That set is "
                "what\ningest_similarity partitions on, and phase 2 assembles "
                "the PREVIOUS run's\nphase-3 output, so similarity has not "
                "yet seen the graph this run\nproduced. A second rebuild "
                "either converges or says this again.\n(fix-list L20)")
    return changed


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--plan", action="store_true",
                    help="print the steps and exit, running nothing")
    ap.add_argument("--root", type=Path, default=ROOT)
    args = ap.parse_args(argv)

    # Taken before any step runs, so the comparison at the end is a
    # measurement and not a recollection.
    before = release_fingerprint(release_dir(args.root))
    # And the one set that decides whether phase 3 has settled (L20).
    annotating_before = annotating_subjects_sha(release_dir(args.root))

    steps = plan(args.root)
    layers = [s for s in steps if s["kind"] == "layer"]
    runnable = [s for s in layers if s["run"]]
    skipped = [s for s in layers if not s["run"]]

    print(f"{len(runnable)} layer(s) to run, {len(skipped)} skipped, "
          f"2 assemble passes\n")
    for s in steps:
        if s["kind"] == "assemble":
            print(f"  phase {s['phase']}  ASSEMBLE   {s['script']}")
        elif s["run"]:
            # THE FLAGS, IN THE PLAN. A plan that does not show them cannot
            # be used to check what the rebuild will do, and the one flag
            # that exists is the difference between reading a cache and
            # calling a third-party API mid-rebuild.
            flags = f"  {' '.join(s['args'])}" if s.get("args") else ""
            print(f"  phase {s['phase']}  run        {s['name']}{flags}")
        elif s["required"]:
            print(f"  phase {s['phase']}  REQUIRED   {s['name']}"
                  f"  -- {s['reason']}")
        elif s["existing"]:
            when = time.strftime("%Y-%m-%d %H:%M",
                                 time.localtime(s["written"]))
            print(f"  phase {s['phase']}  REUSE      {s['name']}"
                  f"  -- {s['reason']}; keeping the layer written {when}")
        else:
            print(f"  phase {s['phase']}  absent     {s['name']}"
                  f"  -- {s['reason']}")
    print()

    # WHICH REUSED LAYERS WERE WRITTEN BY CODE THAT HAS SINCE CHANGED.
    # The REUSED note further down says reuse is correct only if nothing a
    # layer depends on has changed, and then leaves the reader to know. For
    # the code half of "anything", this answers it: provenance.check_layer
    # hashes the layer's whole import closure.
    #
    # PRINTED BEFORE THE REQUIRED-LAYER REFUSAL, deliberately. A run that
    # cannot start is exactly when someone is reading this output, and a
    # diagnosis withheld until the build is otherwise healthy is withheld
    # when it is most wanted.
    #
    # It is a WARNING, not a refusal. A stale phase-1 layer is reused and
    # merely out of date; a stale PHASE 3 layer is the shape that killed
    # this script twice on 2026-11-02, because the phase-2 assemble
    # consumes last run's phase-3 output, validates it, and exits 1 before
    # the ingest that would replace it can run. So the phase is named and
    # the remedy differs by phase.
    stale = [s for s in layers if s["code"] == "stale"]
    if stale:
        print("WRITTEN BY CODE THAT HAS SINCE CHANGED:")
        for s in stale:
            print(f"  phase {s['phase']}  {s['code_reason']}")
        # PHASE, NOT RUNNABILITY. The first version of this said "phase 3
        # and will NOT be re-run", which is the wrong condition and would
        # have stayed silent through the failure it was written for: on
        # 2026-11-02 gtop WAS runnable and DID re-run. It never got the
        # chance, because phase 2 assembles before phase 3 runs, so the
        # stale output is validated first and the run exits 1 there. A
        # stale phase-3 layer is dangerous whether or not its ingest is
        # queued behind the assemble that will reject it.
        p3 = [s for s in stale if s["phase"] == 3]
        if p3:
            names = ", ".join(s["name"] for s in p3)
            queued = [s for s in p3 if s["run"]]
            print(f"\n  {names} {'is' if len(p3) == 1 else 'are'} PHASE 3, "
                  f"and the phase-2 assemble runs first.")
            print("  It reads the stale output and validates it, so if those "
                  "qualifiers no longer\n  validate this run exits there -- "
                  + ("before the phase-3 ingest that would\n  replace them, "
                     "even though it is queued to run this pass.\n"
                     if queued else "and the ingest is not queued either.\n")
                  + "  Delete the directory instead: an absent layer is left "
                    "out, a stale one is not.")
        print()
    unstamped = [s for s in layers if s["code"] == "unstamped"]
    if unstamped:
        print("UNSTAMPED (built before layer stamps existed, so what wrote "
              "them is unknown):")
        print(f"  {', '.join(s['name'] for s in unstamped)}")
        print("  Rebuilding any of them records it. Not an error.\n")

    stop = blocking(args.root)
    if stop:
        print("REFUSING TO START. Required layer(s) cannot run:")
        for line in stop:
            print(f"  {line}")
        print("\nWithout these there is no graph to assemble, and "
              "assemble_release\nwould write a release anyway.")
        return 2

    # TWO KINDS OF SKIP, AND THIS SAID ONLY ONE OF THEM. assemble_release
    # builds its layer set from whatever v0.1-* directories exist, so a
    # skipped ingest whose directory is already on disk does not drop the
    # layer -- it carries the old one forward. The first version of this
    # message claimed every skipped layer was absent from the release, which
    # was wrong for exactly the case that matters most.
    reused = [s for s in skipped if s["existing"]]
    absent = [s for s in skipped if not s["existing"]]
    if absent:
        print("ABSENT from the release this run produces: "
              f"{', '.join(s['name'] for s in absent)}")
        print("  No directory on disk, so assemble_release skips them. Any "
              "metapath that\n  needs one reaches nothing, and the counts "
              "still look right.\n")
    if reused:
        print("REUSED, not rebuilt: "
              f"{', '.join(s['name'] for s in reused)}")
        print("  assemble_release reads every v0.1-* directory that exists, "
              "so these are\n  carried forward from when they were last "
              "written. That is correct only if\n  nothing they depend on "
              "has changed. If you have edited config/viruses.yaml --\n"
              "  promoted a virus, say -- a reused layer has no data for it "
              "and the build\n  will look complete anyway.\n")

    if args.plan:
        print("--plan: nothing was run.")
        return 0

    t0 = time.time()
    for s in steps:
        if s["kind"] == "layer" and not s["run"]:
            continue
        label = s["script"]
        # THE FLAGS THE LAYER DECLARES, from pipeline.LAYERS. Only folds
        # has any (--offline), and until 2026-11-07 there was no way to
        # pass them: this called every script bare, so the one layer that
        # fetches from a live third-party API did so on every rebuild and
        # an InterPro 408 ended the run after seven layers.
        extra = list(s.get("args") or ())
        shown = f"{label} {' '.join(extra)}".strip()
        print(f"\n{'=' * 70}\nphase {s['phase']}  {shown}\n{'=' * 70}",
              flush=True)
        started = time.time()
        rc = subprocess.call([sys.executable,
                              str(args.root / "scripts" / label), *extra])
        took = time.time() - started
        if rc != 0:
            print(f"\nFAILED after {took:.0f}s: {label} exited {rc}.")
            print("Stopping here. The release on disk is from before this "
                  "run, so it is\nconsistent -- rerunning after the fix "
                  "picks up where this stopped,\nsince every layer already "
                  "written is reused.")
            return 1
        print(f"  ok  ({took:.0f}s)")
        # STAMP IT HERE, not inside the ingest. Twelve scripts write their
        # own output and a thirteenth would forget; this is the one place
        # that knows a layer ingest just succeeded.
        if s["kind"] == "layer":
            write_layer(layer_dir(args.root, s["name"]),
                        args.root / "scripts" / s["script"], args.root)

    print(f"\n{'=' * 70}")
    print(f"rebuild complete in {time.time() - t0:.0f}s")
    if absent:
        print(f"WITHOUT: {', '.join(s['name'] for s in absent)}")
    if reused:
        print(f"REUSED (not rebuilt this run): "
              f"{', '.join(s['name'] for s in reused)}")
    after = release_fingerprint(release_dir(args.root))
    print(staleness_verdict(
        before, after,
        annotating_before, annotating_subjects_sha(release_dir(args.root))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
