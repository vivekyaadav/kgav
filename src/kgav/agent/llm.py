"""Local model adapter and command-line entry point.

WHAT THE MODEL IS FOR, AND WHAT IT IS NOT FOR.

The model narrates a brief that has already been assembled deterministically.
It does not choose what to fetch, decide whether a question may be answered,
resolve an entity, or rank anything. Those are the decisions where a plausible
wrong answer is worse than an error, and they are fixed code.

This matters for model selection: a 7B model is sufficient here precisely
because the hard reasoning is already done. The task is "restate these fifteen
facts in prose, keep every citation, keep every warning" -- constrained
narration, not open reasoning.

VERIFICATION IS NOT OPTIONAL AND NOT ADVISORY.

Synthesis is followed automatically by verification against the brief, and an
answer that fails is NOT PRINTED. The failure report is printed instead. A tool
that shows a hallucinated answer with a warning underneath is a tool that shows
hallucinated answers -- the caveat is read second, if at all, and the fluent
text is what sticks.

WHY OLLAMA RATHER THAN AN API. The model sees retrieved facts about compounds
and targets, which for academic work is unremarkable but for anything
commercial would leave the machine. Local also means no per-query cost, so the
verifier can retry without a budget conversation.
"""
from __future__ import annotations

import json
import urllib.error
import urllib.request
from dataclasses import dataclass

DEFAULT_MODEL = "qwen2.5:7b-instruct-q5_K_M"
DEFAULT_HOST = "http://localhost:11434"


class ModelUnavailable(RuntimeError):
    """Raised when the model cannot be reached.

    Distinct from a model that answers badly: the caller should fall back to
    printing the brief unsynthesised, which is still useful, rather than
    reporting an empty answer.
    """


@dataclass
class OllamaModel:
    model: str = DEFAULT_MODEL
    host: str = DEFAULT_HOST
    temperature: float = 0.0
    num_predict: int = 700
    timeout: int = 180

    def __call__(self, prompt: str) -> str:
        payload = {
            "model": self.model,
            "prompt": prompt,
            "stream": False,
            "options": {
                # Deterministic. A researcher who runs the same query twice and
                # gets two different answers cannot cite the tool, and the
                # whole pipeline below this point is deterministic already.
                "temperature": self.temperature,
                "top_p": 1.0,
                "seed": 0,
                "num_predict": self.num_predict,
            },
        }
        req = urllib.request.Request(
            f"{self.host}/api/generate",
            data=json.dumps(payload).encode(),
            headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as r:
                return json.loads(r.read().decode()).get("response", "").strip()
        except urllib.error.URLError as e:
            raise ModelUnavailable(
                f"cannot reach ollama at {self.host}: {e}. Is the service "
                f"running? `systemctl status ollama`") from e

    def available(self) -> tuple[bool, str]:
        try:
            req = urllib.request.Request(f"{self.host}/api/tags")
            with urllib.request.urlopen(req, timeout=10) as r:
                tags = json.loads(r.read().decode()).get("models", [])
            names = [t.get("name", "") for t in tags]
            if self.model in names:
                return True, f"{self.model} ready"
            return False, (f"{self.model} not pulled. Available: "
                           f"{', '.join(names) or 'none'}")
        except Exception as e:                      # noqa: BLE001
            return False, f"ollama unreachable: {e}"
