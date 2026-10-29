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
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "evaluation.md"
HEADING = "## Reproducing"

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

    print(f"running {n} command(s) from {args.doc} '{HEADING}'\n", flush=True)
    with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as fh:
        fh.write(PREAMBLE + body + "\n")
        path = fh.name
    try:
        return subprocess.call(["bash", path], cwd=ROOT)
    finally:
        Path(path).unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
