#!/usr/bin/env python3
"""kgav agent server — the chat frontend's backend.

    make agent-ui                 # then open http://127.0.0.1:8765
    python scripts/agent_server.py --port 8765

WHY A SERVER AND NOT JUST THE CLI. `kgav_ask.py` spends ~8 seconds loading
the release before it can answer anything, and it reloads it for every
question. Here the release, the calibrations, the landscape and the reach
table are loaded once at startup and held, so a question costs what the
question costs. That is the whole reason this file exists; it adds no
capability the CLI does not have.

WHAT IT DOES NOT DO. It does not reimplement the agent. It calls
`kgav_ask.build_orchestrator` for loading and `kgav_ask.json_payload` for
the response, which is the surface the CLI's --json flag serves. Those two
functions live at module level in that script precisely so a second caller
cannot quietly grow a second opinion -- the project already paid for that
lesson twice, once when the agent quoted SARS-CoV-2's AUC for a question
about Zika, and once when --json printed an answer that `render` withheld.

THE WITHHOLDING IS THE POINT, AND IT IS NOT REIMPLEMENTED HERE.
`json_payload` sets `answer` to None and fills `answer_withheld` when
verification fails. This server serves that payload as it comes back, so
an answer the CLI refuses to print is an answer the browser never
receives -- not one it receives and declines to render, which would put
the model's unverified text one devtools tab away from a reader.

BINDS LOCALHOST ONLY, by default and on purpose. The payloads name
compounds, targets and unpublished AUCs, and the model is reached at
localhost:11434. --bind exists for anyone who means to expose it and says
so; it is not the default, because a research tool that listens on every
interface the first time it is run is a tool that leaks on a conference
wifi.

ONE QUESTION AT A TIME. The orchestrator holds the loaded graph and the
model is serial anyway, so /ask takes a lock. Concurrent callers queue
rather than racing a shared GraphTools.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import socket
import sys
import threading
import time
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.agent.llm import DEFAULT_MODEL, ModelUnavailable, OllamaModel
from kgav.agent.verify import verify_response
from kgav.provenance import release_fingerprint

#: The chat UI. One file, no build step, no node_modules, no CDN: it is
#: served from disk on every request so editing it needs no restart.
UI = ROOT / "web" / "chat.html"

#: Largest request body accepted. A question is a sentence; anything of
#: this size is a mistake or an attack, and reading it into memory first
#: to find out is how a server falls over.
MAX_BODY = 16 * 1024


def _load_kgav_ask():
    """scripts/kgav_ask.py as a module.

    Spelled the way tests/test_agent_cli_wiring.py spells it, and for the
    same reason: the script is not a package member, and the alternative
    to importing it is duplicating it.
    """
    spec = importlib.util.spec_from_file_location(
        "kgav_ask_for_server", ROOT / "scripts" / "kgav_ask.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


class Agent:
    """The loaded agent, plus the facts the UI displays about it.

    EVERY FIGURE THE PAGE SHOWS COMES FROM HERE. The UI hand-types no
    count, no sha and no model name; it asks /meta. The previous frontend
    typed five release counts into a source file, four were wrong, and the
    test that caught it had to parse TSX with a regex to do so. A number
    that exists in one place cannot drift from itself.
    """

    def __init__(self, argv: list[str]) -> None:
        ask = _load_kgav_ask()
        self.ask = ask
        # THE CLI'S OWN PARSER, not a dict of defaults copied from it. A
        # default that moves, moves for both surfaces. The positional
        # question is required by the parser and unused here.
        self.args = ask.build_parser(ROOT).parse_args([*argv, "placeholder"])
        self.notes: list[str] = []

        if not (self.args.release / "nodes.jsonl").exists():
            raise SystemExit(f"release not found at {self.args.release}")

        t0 = time.monotonic()
        self.orch = ask.build_orchestrator(self.args, note=self.notes.append)
        self.load_seconds = round(time.monotonic() - t0, 2)
        self.lock = threading.Lock()
        self.asked = 0

    # ------------------------------------------------------------- meta
    def model(self) -> OllamaModel:
        return OllamaModel(model=self.args.model, host=self.args.host)

    def meta(self) -> dict:
        """What the page needs to describe the system it is talking to."""
        ok, note = self.model().available()
        fp = release_fingerprint(self.args.release)
        manifest = self.args.release / "MANIFEST.json"
        census: dict = {}
        if manifest.exists():
            try:
                census = json.loads(manifest.read_text()).get(
                    "nodes_by_class", {})
            except (OSError, ValueError):
                census = {}
        return {
            "release": self.args.release.name,
            "release_path": str(self.args.release),
            "nodes_sha": fp.get("nodes_sha"),
            "edges_sha": fp.get("edges_sha"),
            "nodes_by_class": census,
            "node_total": sum(census.values()) if census else None,
            "model": self.args.model,
            "model_host": self.args.host,
            "model_available": ok,
            "model_note": None if ok else note,
            "results_file": str(self.args.results),
            "landscape_file": str(self.args.landscape),
            "load_seconds": self.load_seconds,
            "asked": self.asked,
            # Degradations, as data. The CLI prints these to stderr, where
            # a browser user would never see them, and every one of them
            # changes what the caveats are allowed to claim.
            "notes": list(self.notes),
        }

    # -------------------------------------------------------------- ask
    def answer(self, question: str, use_model: bool) -> dict:
        """One question. Returns exactly the CLI's --json payload, plus timing.

        The model fallback is the CLI's: an unreachable model degrades to
        the unsynthesised brief, because the evidence is useful without
        prose and a missing model is not a missing answer.
        """
        with self.lock:
            t0 = time.monotonic()
            model = None
            model_note = None
            if use_model:
                m = self.model()
                ok, note = m.available()
                if ok:
                    model = m
                else:
                    model_note = note

            try:
                out = self.orch.answer(question, model)
            except ModelUnavailable as exc:
                model_note = str(exc)
                out = self.orch.answer(question, None)

            # Identical condition to the CLI's. A clarification is a system
            # message rather than a synthesised claim, and keying this on
            # "has no edge ids" instead of on the orchestrator's own flag is
            # how the model once got to invent three edge ids that printed.
            verification = None
            if out["allowed"] and out.get("answer") and not out.get(
                    "clarification"):
                verification = verify_response(out)

            payload = self.ask.json_payload(out, verification)
            self.asked += 1

        payload["elapsed_seconds"] = round(time.monotonic() - t0, 2)
        payload["model_used"] = model is not None
        payload["model_note"] = model_note
        payload["synthesis_requested"] = use_model
        return payload


class Handler(BaseHTTPRequestHandler):
    server_version = "kgav-agent"
    agent: Agent                      # set on the server instance below

    # -------------------------------------------------------- plumbing
    def log_message(self, fmt, *a):   # noqa: A003 - BaseHTTPRequestHandler API
        sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % a))

    def _send(self, code: int, body: bytes, ctype: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        # This serves a local research tool over http on loopback. No
        # caching, so an edited chat.html is the one that loads.
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, code: int, obj: dict) -> None:
        self._send(code, json.dumps(obj, indent=2).encode(),
                   "application/json; charset=utf-8")

    # ------------------------------------------------------------- GET
    def do_GET(self) -> None:                       # noqa: N802 - API
        path = self.path.split("?", 1)[0]
        if path in ("/", "/index.html"):
            if not UI.exists():
                self._json(500, {"error": f"UI missing at {UI}"})
                return
            self._send(200, UI.read_bytes(), "text/html; charset=utf-8")
        elif path == "/meta":
            self._json(200, self.agent.meta())
        else:
            self._json(404, {"error": f"no route {path}"})

    # ------------------------------------------------------------ POST
    def do_POST(self) -> None:                      # noqa: N802 - API
        if self.path.split("?", 1)[0] != "/ask":
            self._json(404, {"error": f"no route {self.path}"})
            return
        try:
            length = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            self._json(400, {"error": "bad Content-Length"})
            return
        if length > MAX_BODY:
            self._json(413, {"error": f"body over {MAX_BODY} bytes"})
            return
        try:
            req = json.loads(self.rfile.read(length) or b"{}")
        except ValueError as exc:
            self._json(400, {"error": f"body is not JSON: {exc}"})
            return

        question = (req.get("question") or "").strip()
        if not question:
            self._json(400, {"error": "no question"})
            return
        # Default TRUE, matching the CLI, where synthesis is on unless
        # --no-model is passed. The UI's toggle is the --no-model flag.
        use_model = bool(req.get("use_model", True))

        try:
            self._json(200, self.agent.answer(question, use_model))
        except Exception as exc:                    # noqa: BLE001
            # A traceback to the terminal and a short message to the page.
            # Serving the traceback would put file paths in a browser; not
            # serving anything would leave the UI spinning forever.
            traceback.print_exc()
            self._json(500, {"error": f"{type(exc).__name__}: {exc}"})


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Serve the kgav chat frontend and the agent behind it.")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--bind", default="127.0.0.1",
                    help="interface to listen on. LOOPBACK BY DEFAULT: the "
                         "payloads name compounds, targets and unpublished "
                         "AUCs. Pass 0.0.0.0 only if you mean it.")
    ap.add_argument("--model", default=DEFAULT_MODEL,
                    help="passed through to the agent")
    ap.add_argument("--no-open", action="store_true",
                    help="do not print the ready-to-click URL banner")
    return ap


def main() -> int:
    args, passthrough = build_parser().parse_known_args()
    # --model is accepted by both parsers; hand it to the agent's.
    agent = Agent([*passthrough, "--model", args.model])

    Handler.agent = agent
    try:
        httpd = ThreadingHTTPServer((args.bind, args.port), Handler)
    except OSError as exc:
        print(f"cannot listen on {args.bind}:{args.port} — {exc}",
              file=sys.stderr)
        return 2
    httpd.daemon_threads = True

    meta = agent.meta()
    if not args.no_open:
        rule = "─" * 62
        print(rule)
        print(f"  kgav agent  →  http://{args.bind}:{args.port}")
        print(rule)
        print(f"  release   {meta['release']}  "
              f"nodes {meta['nodes_sha']} / edges {meta['edges_sha']}")
        print(f"  graph     loaded in {meta['load_seconds']}s"
              + (f", {meta['node_total']:,} nodes"
                 if meta["node_total"] else ""))
        if meta["model_available"]:
            print(f"  model     {meta['model']} at {meta['model_host']}")
        else:
            print(f"  model     UNAVAILABLE ({meta['model_note']}) — "
                  f"evidence will be served without synthesis")
        for n in meta["notes"]:
            print(f"  note      {n}")
        print(rule)
        print("  ctrl-c to stop")
        print()
        # FLUSHED, because stdout is block-buffered the moment it is not a
        # terminal. `make agent-ui | tee log` and `nohup ... > log` both
        # held the whole banner -- the URL included -- in a buffer until
        # the process exited, so the one line a reader needs in order to
        # use the server appeared only after they killed it.
        sys.stdout.flush()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping")
    finally:
        httpd.server_close()
    return 0


def free_port() -> int:
    """An unused port, for tests that need to bind a real server."""
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


if __name__ == "__main__":
    raise SystemExit(main())
