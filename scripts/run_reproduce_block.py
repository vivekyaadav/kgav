#!/usr/bin/env python3
"""Run the Reproducing block in docs/evaluation.md, from the document.

ONE SOURCE OF TRUTH. The block is 21 commands whose ORDER is load
bearing -- three temporal_split.py runs write to one directory so the
last one is the split everything downstream reads, and six
hard_negatives.py variants write to one file so the last one is the
row the document quotes. tests/test_reproduce_block.py already fails
if a results file exists that no step here regenerates.

What was missing was a way to RUN it. It was being transcribed by
hand after every rebuild, and a dropped flag there silently changes a
protocol rather than failing.

NOT `set -e`. One step is required to exit non-zero: the
--overlap-policy keep control must refuse, naming the compounds that
appear on both sides of the cutoff. Aborting there would skip every
step after it, which is most of the block.

    make reproduce           run it
    make reproduce-print     print it and run nothing

ONE AT A TIME. Two concurrent runs destroy each other: the three
temporal_split.py protocols write to ONE directory and the six
hard_negatives.py variants to ONE file, so the block is correct only
as a sequence. On 2026-10-30 a second run was started 42 minutes into
a first; its step 1 rewrote data/releases/v0.1-train2021 with the
`exclude` split while the first run's hard_negatives.py was reading
the `include` one, 40 seconds before it wrote its output.

Nothing detected that. check-results stamps each file against the
graph, and both runs used the same graph -- so the mixture would have
been certified current. It was caught by reading timestamps.

So this takes an exclusive lock and refuses to start while another
run holds it. The lock lives on a file descriptor, so it is released
when the process exits however it exits, including SIGKILL.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import fcntl
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "evaluation.md"
HEADING = "## Reproducing"

LOCK = ROOT / ".reproduce.lock"
# A FIXED NAME, not tempfile. A killed run never reached the cleanup
# in `finally`, leaving /tmp/tmp*.sh behind every time. One path,
# overwritten, gitignored, and inspectable while a run is going.
SCRIPT = ROOT / ".reproduce.sh"

PREAMBLE = (
    "#!/usr/bin/env bash\n"
    "# EXTRACTED FROM docs/evaluation.md -- edit the document, not this.\n"
    "# Not set -e: the --overlap-policy keep control is REQUIRED to exit\n"
    "# non-zero, and aborting there would skip the rest of the block.\n"
    "set -u\n"
)


def block(doc_text: str) -> str:
    """The bash block under '## Reproducing'.

    Raises rather than returning empty: a silently empty block would
    run nothing and report success, which is the shape of every
    vacuous check this project has had to remove.
    """
    try:
        i = doc_text.index(HEADING)
    except ValueError as e:
        raise LookupError(f"{HEADING!r} not found in the document") from e
    m = re.search(r"```bash\n(.*?)\n```", doc_text[i:], re.S)
    if not m:
        raise LookupError(f"no bash block under {HEADING!r}")
    body = m.group(1)
    if not [ln for ln in body.splitlines()
            if ln.strip() and not ln.strip().startswith("#")]:
        raise LookupError(f"the block under {HEADING!r} holds no commands")
    return body


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--doc", type=Path, default=DOC)
    ap.add_argument("--print", dest="show", action="store_true",
                    help="print the block and run nothing")
    args = ap.parse_args(argv)

    try:
        body = block(args.doc.read_text())
    except (OSError, LookupError) as e:
        print(f"cannot read the Reproducing block: {e}", file=sys.stderr)
        return 2

    n = len([ln for ln in body.splitlines()
             if ln.strip() and not ln.strip().startswith("#")])
    if args.show:
        print(body)
        return 0

    fd = os.open(LOCK, os.O_RDWR | os.O_CREAT, 0o644)
    lock = os.fdopen(fd, "r+")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        lock.seek(0)
        who = lock.read().strip() or "(unknown)"
        print(f"ANOTHER REPRODUCE RUN HOLDS THE LOCK:\n    {who}\n\n"
              f"Two runs destroy each other: three temporal_split.py "
              f"protocols write to\none directory and six "
              f"hard_negatives.py variants to one file, so the block "
              f"is\ncorrect only as a sequence. Wait for it, or stop "
              f"it, then start again.\n\nLock: {LOCK}",
              file=sys.stderr)
        return 3
    lock.seek(0)
    lock.truncate()
    lock.write(f"pid {os.getpid()} started "
               f"{_dt.datetime.now().isoformat(timespec='seconds')}\n")
    lock.flush()

    print(f"running {n} command(s) from {args.doc} '{HEADING}'\n", flush=True)
    SCRIPT.write_text(PREAMBLE + body + "\n")
    return subprocess.call(["bash", str(SCRIPT)], cwd=ROOT)


if __name__ == "__main__":
    raise SystemExit(main())
