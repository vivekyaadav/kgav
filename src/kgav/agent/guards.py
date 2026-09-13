"""What the agent is permitted to claim.

This module exists because the evaluation produced a specific, measured account
of what this graph can and cannot support, and an agent that ignores it would
be worse than no agent. The failures are not random: they are systematic,
plausible-looking, and in one case actively harmful.

WHAT THE MEASUREMENTS SAY

  Direct-acting paths separate SELECTIVE antivirals from measured-inactive
  compounds at AUC 0.823, cross-sectionally, for SARS-CoV-2.

  The same paths have NO discrimination prospectively: AUC 0.500 under a
  temporal split, because 77% of compounds screened after the cutoff are absent
  from the graph built before it.

  Host-directed paths are BELOW chance against measured negatives (AUC
  0.391-0.482) and get worse when cytotoxic compounds are removed.

  Nothing works on any virus other than SARS-CoV-2; on four of them a
  degree-only baseline outperforms every metapath.

THE THREE RULES THAT FOLLOW

  1. Refuse prospective questions. "What should we test against X?" is the
     question the graph cannot answer, and answering it anyway produces
     confident, mechanistically coherent, wrong suggestions.
  2. Attach a cell-context warning to every host-directed path. Chloroquine's
     route through SIGMAR1 is real, well-replicated, and led the field astray
     in 2020. The graph cannot encode that the route fails in airway
     epithelium.
  3. State selectivity status on every compound claim. 38% of apparently
     active compounds are cytotoxic, and a compound with no verified index is
     not the same as one that passed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

# Measured on release v0.1. Quoted by the refusal messages so a user is given
# the reason rather than a bare no.
TEMPORAL_AUC = 0.500
CROSS_SECTIONAL_AUC = 0.823
ABSENT_FRACTION = 0.77
SUPPORTED_VIRUS = "NCBITaxon:2697049"
SI_THRESHOLD = 10.0

# Host factors whose route is known to be cell-context dependent. Chloroquine
# and hydroxychloroquine act through the first two; the endosomal machinery is
# the route that made them look effective in Vero E6 and useless in airway
# cells.
CELL_CONTEXT_SENSITIVE = {
    "SIGMAR1", "TMEM97",                      # sigma receptors
    "ATP6AP1", "ATP6AP2", "ATP6V1A", "ATP6V0C", "ATP6V1B2",   # v-ATPase
    "CTSL", "CTSB",                           # endosomal proteases
    "RAB7A", "RAB5A", "NPC1",                 # endosomal trafficking
    "PIK3C3",                                 # autophagy/endosome
}


class QueryKind(str, Enum):
    """What a user is asking for. Determines whether the graph may answer."""
    EXPLAIN = "explain"          # why might this compound work? -> supported
    TRIAGE = "triage"            # rank these compounds -> supported
    EVIDENCE = "evidence"        # what supports this claim? -> supported
    PROFILE = "profile"          # what is known about X? -> supported
    PROSPECT = "prospect"        # what should we test? -> REFUSED
    UNKNOWN = "unknown"


@dataclass
class Verdict:
    """Whether a request may proceed, and what must accompany the answer."""
    allowed: bool
    reason: str = ""
    warnings: list[str] = field(default_factory=list)


def classify_intent(text: str) -> QueryKind:
    """Coarse intent classification.

    Deliberately keyword-based rather than model-based. A language model asked
    to classify intent will occasionally route a prospective question into a
    permitted category, and the whole point of this layer is that the refusal
    must not be probabilistic.
    """
    t = (text or "").lower()

    prospective = (
        "what should we test", "what should i test", "which drug should",
        "suggest a drug", "suggest drugs", "recommend a drug", "recommend drugs",
        "find a cure", "find a treatment", "what would work against",
        "propose candidates", "novel candidates", "new candidates",
        "predict which", "best drug for",
    )
    if any(p in t for p in prospective):
        return QueryKind.PROSPECT

    if any(p in t for p in ("why might", "why does", "how does", "mechanism",
                            "explain", "how might")):
        return QueryKind.EXPLAIN
    if any(p in t for p in ("rank", "prioritise", "prioritize", "triage",
                            "which of these", "order these")):
        return QueryKind.TRIAGE
    if any(p in t for p in ("evidence", "citation", "source", "which paper",
                            "how do we know", "provenance")):
        return QueryKind.EVIDENCE
    if any(p in t for p in ("what is known", "profile", "tell me about",
                            "what do we have")):
        return QueryKind.PROFILE
    return QueryKind.UNKNOWN


def check(kind: QueryKind, virus: str | None = None) -> Verdict:
    """Whether a query of this kind, about this virus, may be answered."""
    if kind is QueryKind.PROSPECT:
        return Verdict(
            allowed=False,
            reason=(
                "This graph cannot support prospective candidate generation, and "
                "the limit is measured rather than assumed. Under a temporal "
                f"split — training on evidence published before 2022 and testing "
                f"on compounds measured after — every scorer performs at chance "
                f"(AUC {TEMPORAL_AUC:.3f}). The reason is that "
                f"{ABSENT_FRACTION:.0%} of compounds screened after the cutoff do "
                "not appear anywhere in the graph built before it: compounds "
                "enter screening because they are novel chemical matter, which "
                "is exactly what a prior-evidence graph cannot represent.\n\n"
                "What this graph can do: explain the mechanism of a compound "
                "with known activity, rank a set of compounds you already have, "
                "and retrieve provenanced evidence for a specific claim."
            ),
        )

    warnings: list[str] = []
    if virus and virus != SUPPORTED_VIRUS:
        warnings.append(
            "Scoring is not validated for this virus. Cross-sectional AUC is "
            "approximately 0.50 for every metapath on all coronaviruses except "
            "SARS-CoV-2, and on four of them a degree-only baseline outperforms "
            "every mechanistic route. Paths retrieved here are evidence to read, "
            "not a ranking to trust."
        )
    return Verdict(allowed=True, warnings=warnings)


def path_warnings(path_nodes: list[dict], metapath: str | None = None) -> list[str]:
    """Warnings a retrieved path must carry into the answer.

    Attached at retrieval time rather than left to documentation, because the
    failure mode is that a plausible-looking host-directed path is presented
    without the caveat that makes it interpretable.
    """
    out: list[str] = []

    symbols = {
        (n.get("properties") or {}).get("gene_symbol", "").upper()
        for n in path_nodes
    }
    hits = sorted(symbols & CELL_CONTEXT_SENSITIVE)
    if hits:
        out.append(
            f"CELL CONTEXT: this path runs through {', '.join(hits)}, whose "
            "contribution depends on the cell type used. The endosomal and "
            "sigma-receptor routes are prominent in Vero E6 cells and largely "
            "absent in TMPRSS2-expressing airway cells. Chloroquine's apparent "
            "SARS-CoV-2 activity in 2020 came from exactly this route and did "
            "not translate to patients."
        )

    # Route labels carry a description ("M2 host-directed"), so match the
    # identifier prefix rather than the whole string.
    code = (metapath or "").split()[0] if metapath else ""
    if code in {"M2", "M3", "M4", "M5", "M6"}:
        out.append(
            "HOST-DIRECTED ROUTE: measured below chance against known-inactive "
            "compounds (AUC 0.391-0.482 depending on route), and worse once "
            "cytotoxic compounds are excluded. Read this path as a mechanistic "
            "hypothesis to evaluate, not as evidence of activity."
        )
    return out


# Cell lines where an antiviral result is known not to transfer. Vero E6
# expresses little TMPRSS2, so SARS-CoV-2 enters it almost entirely by the
# endosomal route; compounds blocking that route score well there and fail in
# airway epithelium. Chloroquine has a defensible selectivity index of ~15 in
# Vero E6 and no clinical benefit. The same compound's EC50 ranges 330 nM to
# 10,900 nM across cell lines in this graph -- a thirtyfold spread driven by
# cell type alone.
NON_TRANSFERRING_CELL_LINES = ("vero", "bhk", "huh-7", "huh7", "hel 299", "mrc5")


def cell_line_note(cell_line: str | None) -> str:
    if not cell_line:
        return ""
    c = cell_line.lower()
    if any(x in c for x in NON_TRANSFERRING_CELL_LINES):
        return (
            f"CELL LINE: measured in {cell_line}, which lacks meaningful TMPRSS2 "
            "expression. SARS-CoV-2 enters such cells by the endosomal route, so "
            "compounds acting on endosomal acidification score well here and fail "
            "in TMPRSS2-expressing airway cells. A selectivity index from this "
            "system does not transfer."
        )
    return f"CELL LINE: measured in {cell_line}."


def selectivity_note(selectivity_index: float | None,
                     verified: bool | None,
                     n_studies: int = 1,
                     si_max: float | None = None) -> str:
    """One sentence on a compound's cytotoxicity status.

    'No data' and 'passed' are different states and must not read the same.
    38% of apparently active compounds with matched data are cytotoxic.
    """
    if not verified or selectivity_index is None:
        return (
            "SELECTIVITY UNKNOWN: no cytotoxicity measurement paired with an "
            "activity measurement in the same study, so it is not known whether "
            "the observed activity reflects antiviral effect or cell death. Of "
            "compounds where this could be checked, 38% were cytotoxic."
        )
    # How much weight the number carries. 96% of compounds with a verified
    # index have exactly one paired measurement, so quoting a bare value
    # implies a consensus that does not exist.
    if n_studies > 1 and si_max is not None and si_max > selectivity_index:
        basis = (f" Across {n_studies} studies the index ranges "
                 f"{selectivity_index:.1f} to {si_max:.1f}; the lowest is quoted, "
                 f"and cell line and protocol differ between them.")
    else:
        basis = (" This rests on a SINGLE paired measurement, so it indicates "
                 "rather than establishes the compound's selectivity.")

    if selectivity_index < SI_THRESHOLD:
        return (
            f"CYTOTOXIC: selectivity index {selectivity_index:.1f}, below the "
            f"threshold of {SI_THRESHOLD:g}. The apparent activity may reflect "
            f"cell death rather than antiviral effect." + basis
        )
    return (
        f"SELECTIVE: selectivity index {selectivity_index:.1f} "
        f"(CC50/EC50, both measured in the same study)." + basis
    )


def evidence_note(tier: int | None, detection_method: str | None = None) -> str:
    """How the underlying fact was established."""
    if detection_method and detection_method.startswith(("MI:1314", "MI:0400",
                                                         "MI:0676", "MI:0007")):
        return (
            "PROXIMITY EVIDENCE: this interaction was detected by proximity "
            "labelling or affinity co-purification, which shows the proteins are "
            "in the same complex neighbourhood rather than that they bind. "
            "Only 1.3% of virus-host edges in this graph come from a binary "
            "binding assay."
        )
    return {
        1: "Curated experimental result.",
        2: "Curated but inferred rather than directly measured.",
        3: "Computationally derived, not experimentally asserted.",
        4: "Text-mined; not curated.",
    }.get(tier or 0, "Provenance not recorded.")
