"""What the agent is permitted to claim.

This module exists because the evaluation produced a specific, measured account
of what this graph can and cannot support, and an agent that ignores it would
be worse than no agent. The failures are not random: they are systematic,
plausible-looking, and in one case actively harmful.

WHAT THE MEASUREMENTS SAY

  Direct-acting paths separate SELECTIVE antivirals from measured-inactive
  compounds cross-sectionally for SARS-CoV-2. The FIGURE is not written
  here: it depends on the protocol, and four protocols gave four answers.

  The same paths have NO discrimination prospectively: AUC 0.500 under a
  temporal split, because 77% of compounds screened after the cutoff are absent
  from the graph built before it.

  Host-directed paths are BELOW chance against measured negatives (AUC
  below chance on an early run) are NOT hardcoded here: _channel_caveat
  renders the figure from the calibration it is handed, and says the
  route is uncalibrated when handed none.

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
from typing import Self

# TEMPORAL_AUC = 0.500 AND ABSENT_FRACTION = 0.77 WERE HERE, both deleted.
#
# They were measured figures hardcoded into user-facing prose, which is the
# exact bug _channel_caveat and direct_acting_caveat were each rewritten to
# remove -- and their docstrings say so. The refusal message kept doing it.
#
# ABSENT_FRACTION was worse than stale, it was the wrong SHAPE. Measured on
# this release the post-cutoff absent fraction is 85.3% for SARS-CoV-2
# (1,343/1,575) and 99.0% for Zika (575/581), so there is no single number for
# a module-level constant to hold. One global value was always a per-virus
# quantity in disguise; a second viral family just made that obvious.
#
# TEMPORAL_AUC 0.500 happens to still be right -- post-2021 M1 reads 0.502 for
# SARS-CoV-2 and 0.502 for Zika -- which is precisely why it was dangerous.
# A figure that is correct today and unlinked to any run is one release away
# from being quoted wrongly with full confidence.
#
# The refusal now renders from a temporal Calibration when one is supplied and
# states the limit without a figure when none is. The MECHANISM is the durable
# claim and needs no number.
# CROSS_SECTIONAL_AUC = 0.826 WAS HERE, dead: defined, never read, while the
# number users actually saw was a separate literal in tools.py's ranking
# caveat. Two copies of one figure with nothing keeping them equal, and the
# figure itself has since fallen to 0.726 (filter removed), 0.660 (filter made
# symmetric) and 0.499 (same-paper evidence withheld). Direct-acting
# performance is rendered by direct_acting_caveat() from a Calibration now.
SUPPORTED_VIRUS = "NCBITaxon:2697049"
SI_THRESHOLD = 10.0

# BIOLOGY, NOT MEASUREMENT, WHICH IS WHY IT STAYS IN PROSE.
#
# Everything else quoted to a user in this module is rendered from a
# Calibration, because a measured figure is release-dependent and a hardcoded
# one goes stale -- that lesson cost two rewrites and a deleted constant. The
# statements below are different in kind: published virology that does not
# change when v0.1 is rebuilt. They belong in prose. What they must NOT do is
# describe one viral family's biology while answering about another, which is
# what they did until this patch.
#
# AND THEY CARRY SOURCES NOW. Every performance figure in this module is
# traceable to the evaluation that produced it, while the biology was asserted
# in bare prose -- four specific, falsifiable claims about deleted gene loci,
# entry routes and interferon restriction, with nothing a reader could check.
# A graph whose whole argument is provenance should not make its most
# confident statements the only unsourced ones. BIOLOGY_SOURCES holds the
# reference; cite() renders the short form into the prose, so the two cannot
# drift apart and an unknown key raises rather than rendering nothing.
#
# key -> (short id rendered inline, full reference)
#
# DOIs throughout rather than a mix: a DOI was verified for all four, a PMID
# for one, and a citation style that varies by what happened to be findable
# invites the reader to wonder what the difference means.
# Fragments joined with an explicit + rather than by adjacency: inside a
# collection literal a MISSING COMMA silently concatenates two entries,
# and in a citation registry that would merge one reference into another
# without a word. With +, the same slip is a syntax error.
BIOLOGY_SOURCES: dict[str, tuple[str, str]] = {
    "vero_ifn_deletion": (
        "doi:10.1093/dnares/dsu029",
        "Osada N, Kohara A, Yamaji T, Hirayama N, Kasai F, Sekizuka T, "
        + "Kuroda M, Hanada K (2014). The genome landscape of the African "
        + "green monkey kidney-derived Vero cell line. DNA Research "
        + "21(6):673-683. doi:10.1093/dnares/dsu029. PMID:25267831. "
        + "Reports a homozygous ~9 Mb deletion on chromosome 12 removing the "
        + "type I interferon gene cluster, along with CDKN2A and CDKN2B."),
    "chloroquine_not_in_airway": (
        "doi:10.1038/s41586-020-2575-3",
        "Hoffmann M, Mosbauer K, Hofmann-Winkler H, Kaul A, Kleine-Weber H, "
        + "Kruger N, et al. (2020). Chloroquine does not inhibit infection of "
        + "human lung cells with SARS-CoV-2. Nature 585(7826):588-590. "
        + "doi:10.1038/s41586-020-2575-3. Chloroquine blocks entry in Vero "
        + "and not in TMPRSS2-expressing Calu-3 lung cells."),
    "flavivirus_entry": (
        "doi:10.1371/journal.ppat.1000244",
        "van der Schaar HM, Rust MJ, Chen C, van der Ende-Metselaar H, "
        + "Wilschut J, Zhuang X, Smit JM (2008). Dissecting the cell entry "
        + "pathway of dengue virus by single-particle tracking in living "
        + "cells. PLoS Pathogens 4(12):e1000244. "
        + "doi:10.1371/journal.ppat.1000244. Dengue enters by clathrin-"
        + "mediated endocytosis and fuses from within acidified endosomes."),
    "flavivirus_ifn_restriction": (
        "doi:10.1016/j.chom.2016.03.010",
        "Lazear HM, Govero J, Smith AM, Platt DJ, Fernandez E, Miner JJ, "
        + "Diamond MS (2016). A mouse model of Zika virus pathogenesis. Cell "
        + "Host & Microbe 19(5):720-730. doi:10.1016/j.chom.2016.03.010. "
        + "Zika produces lethal neurological disease in mice lacking the type "
        + "I interferon receptor (Ifnar1-/-), which is why an "
        + "interferon-deficient background is needed to model it. The claim "
        + "supported is that restriction, not a wild-type comparison the "
        + "abstract does not report."),
}


def cite(key: str) -> str:
    """The short reference to render inline. Raises on an unknown key.

    Raising is the point: a claim whose source has been removed must break
    the render rather than quietly lose its citation.
    """
    return f"({BIOLOGY_SOURCES[key][0]})"


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


# --------------------------------------------------------------- caveats
# Stable key -> the visible label an answer must carry. The LABEL is what
# verification enforces, so the prose below it can be rewritten freely without
# changing what is checked.
#
# The trap this closes: verification used to derive its key from the caveat's
# own text -- the part before the first colon, or the first 40 characters when
# there was no colon. So a caveat without a colon demanded near-verbatim
# reproduction, and editing any caveat's wording silently changed the
# requirement without touching the verifier. Editing a LABEL is now a visible,
# deliberate act; editing prose is free.
CAVEAT_LABELS: dict[str, str] = {
    "cell_context":          "CELL CONTEXT",
    "host_directed_route":   "HOST-DIRECTED ROUTE",
    "cell_line":             "CELL LINE",
    "selectivity_unknown":   "SELECTIVITY UNKNOWN",
    "selectivity_cytotoxic": "CYTOTOXIC",
    "selectivity_selective": "SELECTIVE",
    "proximity_evidence":    "PROXIMITY EVIDENCE",
    "evidence_tier":         "EVIDENCE",
    "virus_not_validated":   "VIRUS NOT VALIDATED",
    "virus_assumed":         "VIRUS ASSUMED",
    "ranking":               "RANKING",
}


class Caveat(str):
    """A mandatory caveat: prose, plus a stable key and a visible label.

    Subclasses str so it renders, joins, dedups and serialises exactly like the
    plain string it replaces -- every existing call site keeps working. What it
    adds is identity: `key` names the caveat and `label` is the marker the
    answer must carry, both independent of the wording.
    """

    key: str
    label: str

    def __new__(cls, key: str, text: str) -> Self:
        label = CAVEAT_LABELS[key]
        obj = super().__new__(cls, f"{label}: {text}")
        obj.key = key
        obj.label = label
        return obj


def caveat(key: str, text: str) -> Caveat:
    """Build a labelled caveat. Unknown keys raise rather than pass silently."""
    return Caveat(key, text)


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


def _measured_at_chance(temporal) -> str:
    """The temporal result, from the run that measured it, or no figure."""
    if temporal is None or not getattr(temporal, "usable", False):
        return "every scorer performs at chance"
    evaluable = [c for c in temporal.channels.values() if c.evaluable]
    if not evaluable:
        return ("no scorer was even evaluable — each reaches too few of the "
                "post-cutoff actives for a measured performance")
    best = max(evaluable, key=lambda c: c.auc)
    if best.discriminates:
        # Stated rather than suppressed: if a channel ever clears chance
        # prospectively, a refusal built on "every scorer is at chance" has
        # stopped being true and the person is entitled to know that before
        # being refused.
        return (f"the strongest scorer reaches AUC {best.auc:.3f} "
                f"[{best.ci_lo:.3f}, {best.ci_hi:.3f}] under "
                f"{temporal.protocol}, which clears chance")
    return (f"every scorer performs at chance — the strongest is "
            f"{best.name} at AUC {best.auc:.3f} "
            f"[{best.ci_lo:.3f}, {best.ci_hi:.3f}] under {temporal.protocol}")


def _virus_caveat(calibration) -> str:
    """What is measured for THIS virus, from its own calibration.

    This was a hardcoded sentence asserting things about other viruses --
    "approximately 0.50 for every metapath on all coronaviruses except
    SARS-CoV-2, and on four of them a degree-only baseline outperforms every
    mechanistic route". It is now five, not four (Zika 0.527, HCoV-OC43 0.552,
    SARS-CoV 0.696, MERS-CoV 0.723, HCoV-229E 0.646), and "all coronaviruses"
    does not describe a graph containing a flavivirus.

    check() had the virus in hand and never consulted its calibration, while
    _channel_caveat twenty lines below did exactly that. Same bug, same file.
    """
    if calibration is None or not getattr(calibration, "usable", False):
        return ("no channel performance has been measured for this virus on "
                "this release, so its reliability is unknown. Paths retrieved "
                "here are evidence to read, not a ranking to trust.")
    disc = calibration.discriminating()
    if disc:
        return (f"scoring for this virus is measured under "
                f"{calibration.protocol}: {', '.join(disc)} clear chance and "
                f"no other channel does. Paths outside those are evidence to "
                f"read, not a ranking to trust.")
    return (f"scoring is not validated for this virus: "
            f"{calibration.why_not()} Paths retrieved here are evidence to "
            f"read, not a ranking to trust.")


def check(kind: QueryKind, virus: str | None = None,
          calibration=None, temporal=None) -> Verdict:
    """Whether a query of this kind, about this virus, may be answered.

    `calibration` is this virus's cross-sectional Calibration and `temporal`
    its post-cutoff one. Both optional: absent, the messages state the limit
    without a figure, which is weaker than a measurement and true.
    """
    if kind is QueryKind.PROSPECT:
        return Verdict(
            allowed=False,
            reason=(
                "This graph cannot support prospective candidate generation, and "
                "the limit is measured rather than assumed. Under a temporal "
                "split — training on evidence published before 2022 and testing "
                f"on compounds measured after — {_measured_at_chance(temporal)}. "
                "The reason is that the large majority of compounds screened "
                "after the cutoff do not appear anywhere in the graph built "
                "before it, and the fraction differs by virus — which is why "
                "no figure is quoted here as though it were one number. "
                "Compounds enter screening because they are novel chemical "
                "matter, which is exactly what a prior-evidence graph cannot "
                "represent. temporal_split.py reports the per-virus ceiling "
                "for the release in hand.\n\n"
                "What this graph can do: explain the mechanism of a compound "
                "with known activity, rank a set of compounds you already have, "
                "and retrieve provenanced evidence for a specific claim."
            ),
        )

    warnings: list[str] = []
    if virus and virus != SUPPORTED_VIRUS:
        warnings.append(caveat("virus_not_validated", _virus_caveat(calibration)))
    return Verdict(allowed=True, warnings=warnings)


def virus_assumed_caveat(label: str, others: list[str]) -> str:
    """Said out loud when no virus was named and one was assumed.

    Orchestrator.build_brief reads `self._virus_in(resolved) or
    self.default_virus`, so a question naming only a compound -- "is
    chloroquine useful?" -- was answered about SARS-CoV-2 with nothing
    marking the choice as an assumption. Harmless while the graph held one
    virus. The graph now holds eight taxa across TWO FAMILIES, so a question
    a user meant about Zika was answered about a coronavirus silently, and
    the virus decides the calibration, the caveats and the cell-context
    biology.

    The default itself is not wrong: SARS-CoV-2 is the only virus any
    protocol here evaluated. Defaulting without saying so is.
    """
    tail = (f" The graph also holds {', '.join(others[:6])}"
            + (f" and {len(others) - 6} more" if len(others) > 6 else "")
            + ".") if others else ""
    return (f"the question named no virus, so this answer is about {label}, "
            f"which is the only virus any protocol in this project "
            f"evaluated.{tail} Name a virus to ask about another.")


HOST_DIRECTED_CODES = frozenset({"M2", "M3", "M4", "M5", "M6"})


def _channel_caveat(code: str, calibration) -> str:
    """What this route is MEASURED to be worth, from the run that measured it.

    A HARDCODED NUMBER IS THE BUG THIS REPLACES. This caveat used to read
    "measured below chance (AUC 0.391-0.482 depending on route)", a figure
    from a run nobody could identify, and by the time the single-screen
    protocol was in place the real values were 0.492-0.523 -- at chance, not
    below it. So the agent was stating a stale measurement in user-facing
    text with nothing linking it to a run, which is precisely what
    provenance.py and calibration.py exist to prevent, reappearing inside the
    agent's own output.

    With no calibration the route is reported as UNCALIBRATED rather than
    given a number. Saying "its performance has not been measured on this
    release" is weaker than a figure and true, and a figure that might be
    stale is worse than no figure.
    """
    if calibration is None or not getattr(calibration, "usable", False):
        return ("no measured performance is attached to this route in this "
                "release, so its reliability is unknown. Read the path as a "
                "mechanistic hypothesis to evaluate, not as evidence of "
                "activity.")
    ch = calibration.channels.get(code)
    if ch is None or ch.auc is None:
        return (f"{code} was not evaluable under {calibration.protocol}, so "
                f"its reliability is unknown. Read the path as a mechanistic "
                f"hypothesis to evaluate, not as evidence of activity.")
    if not ch.evaluable:
        # An AUC computed mostly over positives the channel never reached
        # describes the imputed floor. Quoting it here would be the hardcoded
        # figure's error in a new place: more confidence than the evidence.
        return (f"{code} has NOT been evaluated under "
                f"{calibration.protocol} — it reaches only "
                f"{ch.n_reached_pos} of {ch.n_pos:,} measured actives, too "
                f"few for a measured performance. Its reliability is "
                f"unknown. Read the path as a mechanistic hypothesis to "
                f"evaluate, not as evidence of activity.")
    if ch.discriminates:
        return (f"{code} separates measured actives from measured inactives "
                f"at AUC {ch.auc:.3f} [{ch.ci_lo:.3f}, {ch.ci_hi:.3f}] under "
                f"{calibration.protocol}. That is above chance, and still a "
                f"mechanistic hypothesis rather than a measured activity for "
                f"this compound.")
    return (f"{code} scores AUC {ch.auc:.3f} [{ch.ci_lo:.3f}, {ch.ci_hi:.3f}] "
            f"against {ch.n_neg:,} compounds measured INACTIVE under "
            f"{calibration.protocol} — indistinguishable from chance. Read "
            f"this path as a mechanistic hypothesis to evaluate, not as "
            f"evidence of activity.")


def direct_acting_caveat(calibration, code: str = "M1") -> str:
    """How direct-acting evidence performed, from the run that measured it.

    The ranking caveat used to state "Direct-acting paths separate selective
    antivirals from measured inactives at AUC 0.826; host-directed paths score
    below chance". Both halves are now wrong. 0.826 came from a protocol that
    filtered positives to compounds with verified selectivity data and left
    negatives whole -- a 17x coverage asymmetry -- and the same metapath reads
    0.499 once evidence sharing a publication with its own label is withheld.
    Host-directed routes are at chance, not below it.

    A ranking is shown to a reader deciding what to test next, so the figure
    attached to it has to be the measured one or absent.
    """
    if calibration is None or not getattr(calibration, "usable", False):
        return ("no measured performance is attached to the direct-acting "
                "route in this release, so this order reflects the evidence "
                "available rather than a demonstrated ability to rank")
    ch = calibration.channels.get(code)
    if ch is None or ch.auc is None:
        return (f"{code} was not evaluable under {calibration.protocol}, so "
                f"this order reflects the evidence available rather than a "
                f"demonstrated ability to rank")
    if ch.discriminates:   # implies evaluable
        return (f"direct-acting evidence ({code}) separates measured actives "
                f"from measured inactives at AUC {ch.auc:.3f} "
                f"[{ch.ci_lo:.3f}, {ch.ci_hi:.3f}] under "
                f"{calibration.protocol}")
    if not ch.evaluable:
        return (f"direct-acting evidence ({code}) has NOT been evaluated "
                f"under {calibration.protocol} — it reaches only "
                f"{ch.n_reached_pos} of {ch.n_pos:,} measured actives — so "
                f"this order reflects the evidence available rather than a "
                f"demonstrated ability to rank")
    return (f"direct-acting evidence ({code}) scores AUC {ch.auc:.3f} "
            f"[{ch.ci_lo:.3f}, {ch.ci_hi:.3f}] under {calibration.protocol} — "
            f"indistinguishable from chance, so this order reflects the "
            f"evidence available rather than a demonstrated ability to rank")


def _family_of(path_nodes: list[dict]) -> str | None:
    """The viral family this path ends at, from the path itself.

    No new plumbing: ingest_spine writes `family` onto every OrganismTaxon
    from config/viruses.yaml, and a path that ends at a virus contains that
    node. src/kgav/chembl.py reads the same property.
    """
    for n in path_nodes or ():
        if n.get("class") == "OrganismTaxon":
            fam = (n.get("properties") or {}).get("family")
            if fam:
                return str(fam)
    return None


def _cell_context_text(hits: list[str], family: str | None) -> str:
    """Why this route's contribution depends on the cell system, per family.

    THE BUG THIS FIXES. One text described SARS-CoV-2 entry -- TMPRSS2, Vero
    E6, chloroquine in 2020 -- and was attached to every path through these
    genes, including a flavivirus one. The genes are genuinely involved in
    both (flavivirus entry is clathrin-mediated and acidification-dependent),
    so the warning belongs on both; the EXPLANATION does not transfer, and an
    explanation about the wrong virus is worse than none because it reads as
    specific knowledge.
    """
    genes = ", ".join(hits)
    if family == "Coronaviridae":
        return (f"this path runs through {genes}, whose contribution depends "
                "on the cell type used. The endosomal and sigma-receptor "
                "routes are prominent in Vero E6 cells and largely absent in "
                "TMPRSS2-expressing airway cells. Chloroquine's apparent "
                "SARS-CoV-2 activity in 2020 came from exactly this route and "
                f"did not translate to patients "
                f"{cite('chloroquine_not_in_airway')}.")
    if family == "Flaviviridae":
        return (f"this path runs through {genes}. Flavivirus entry is "
                "clathrin-mediated and depends on endosomal acidification "
                f"{cite('flavivirus_entry')}, so "
                "this route is genuine rather than a cell-line artifact -- but "
                "whether acting on it helps depends on the system. Vero, the "
                "standard flavivirus line, carries a deletion of the type I "
                f"interferon locus {cite('vero_ifn_deletion')}, and Zika is "
                "restricted by type I interferon to the point that modelling "
                "it needs an interferon-deficient host "
                f"{cite('flavivirus_ifn_restriction')}, so neither an effect "
                "nor its absence in "
                "Vero transfers to interferon-competent cells. Vero is also "
                "monkey kidney epithelium rather than the neural progenitor or "
                "placental tissue that determines Zika pathology.")
    return (f"this path runs through {genes}, whose contribution depends on "
            "the cell type used. No family-specific account of that dependence "
            "is recorded here for this virus, so read the route as a "
            "hypothesis whose cell-system dependence is unassessed.")


def path_warnings(path_nodes: list[dict], metapath: str | None = None,
                  calibration=None) -> list[str]:
    """Warnings a retrieved path must carry into the answer.

    Attached at retrieval time rather than left to documentation, because the
    failure mode is that a plausible-looking host-directed path is presented
    without the caveat that makes it interpretable.

    `calibration` is a kgav.calibration.Calibration. Supplying it replaces a
    hardcoded AUC with the one actually measured, and naming the protocol
    alongside it means a reader can tell which run the figure came from --
    four protocols produced four different AUCs for M1.
    """
    out: list[str] = []

    symbols = {
        (n.get("properties") or {}).get("gene_symbol", "").upper()
        for n in path_nodes
    }
    hits = sorted(symbols & CELL_CONTEXT_SENSITIVE)
    if hits:
        out.append(caveat("cell_context",
                          _cell_context_text(hits, _family_of(path_nodes))))

    # Route labels carry a description ("M2 host-directed"), so match the
    # identifier prefix rather than the whole string.
    code = (metapath or "").split()[0] if metapath else ""
    if code in HOST_DIRECTED_CODES:
        out.append(caveat("host_directed_route",
                          _channel_caveat(code, calibration)))
    return out


# Cell lines where an antiviral result is known not to transfer. Vero E6
# expresses little TMPRSS2, so SARS-CoV-2 enters it almost entirely by the
# endosomal route; compounds blocking that route score well there and fail in
# airway epithelium. Chloroquine has a defensible selectivity index of ~15 in
# Vero E6 and no clinical benefit. The same compound's EC50 ranges 330 nM to
# 10,900 nM across cell lines in this graph -- a thirtyfold spread driven by
# cell type alone.
NON_TRANSFERRING_CELL_LINES = ("vero", "bhk", "huh-7", "huh7", "hel 299", "mrc5")


# Vero's type I interferon locus is DELETED, which is a property of the line
# rather than of any virus -- so the fact is stated once and its consequence
# is drawn per family. BHK, Huh-7, HEL 299 and MRC5 are on the
# non-transferring list for their own reasons, and this does not assert the
# interferon fact about them.
IFN_DEFICIENT_LINES = ("vero",)


def cell_line_note(cell_line: str | None, family: str | None = None) -> str:
    """Whether a measurement in this cell system transfers, per viral family.

    THE BUG THIS FIXES. This said "lacks meaningful TMPRSS2 expression.
    SARS-CoV-2 enters such cells by the endosomal route" for ANY measurement
    in these lines. Vero is also the standard flavivirus line, so a Zika
    selectivity index came back explained by SARS-CoV-2 entry biology -- a
    confident, specific, irrelevant sentence. The reason a Vero result does
    not transfer for a flavivirus is the line's interferon deletion, not
    TMPRSS2, which flaviviruses do not use.
    """
    if not cell_line:
        return ""
    c = cell_line.lower()
    if not any(x in c for x in NON_TRANSFERRING_CELL_LINES):
        return caveat("cell_line", f"measured in {cell_line}.")

    if family == "Coronaviridae":
        return caveat(
            "cell_line",
            f"measured in {cell_line}, which lacks meaningful TMPRSS2 "
            "expression. SARS-CoV-2 enters such cells by the endosomal route, so "
            "compounds acting on endosomal acidification score well here and fail "
            "in TMPRSS2-expressing airway cells "
            f"{cite('chloroquine_not_in_airway')}. A selectivity index from this "
            "system does not transfer."
        )
    if family == "Flaviviridae":
        if any(x in c for x in IFN_DEFICIENT_LINES):
            return caveat(
                "cell_line",
                f"measured in {cell_line}, which carries a deletion of the "
                f"type I interferon locus {cite('vero_ifn_deletion')}. Zika is "
                "restricted by type I interferon to the point that modelling "
                "it needs an interferon-deficient host "
                f"{cite('flavivirus_ifn_restriction')}, so a selectivity index "
                "from this system "
                "captures neither interferon-dependent activity nor its "
                "absence, and does not transfer to interferon-competent cells."
            )
        return caveat(
            "cell_line",
            f"measured in {cell_line}, a line on which antiviral results are "
            "known not to transfer reliably. The specific reason is not "
            "characterised here for this family, so treat the index as "
            "system-dependent."
        )
    return caveat(
        "cell_line",
        f"measured in {cell_line}, a line on which antiviral results are known "
        "not to transfer reliably. No family-specific reason is recorded for "
        "this virus, so treat the index as system-dependent."
    )


def selectivity_note(selectivity_index: float | None,
                     verified: bool | None,
                     n_studies: int = 1,
                     si_max: float | None = None) -> str:
    """One sentence on a compound's cytotoxicity status.

    'No data' and 'passed' are different states and must not read the same.
    38% of apparently active compounds with matched data are cytotoxic.
    """
    if not verified or selectivity_index is None:
        return caveat(
            "selectivity_unknown",
            "no cytotoxicity measurement paired with an "
            "activity measurement in the same study, so it is not known whether "
            "the observed activity reflects antiviral effect or cell death. Of "
            "compounds where this could be checked, 38% were cytotoxic."
        )
    # How much weight the number carries. 96% of compounds with a verified
    # index have exactly one paired measurement, so quoting a bare value
    # implies a consensus that does not exist.
    if n_studies > 1 and si_max is not None and si_max > selectivity_index:
        # "the lowest is quoted" was false: the ingest quotes the MEDIAN of
        # the exact measurements, which is less conservative. Taking the worst
        # was the earlier rule and it reported nirmatrelvir at SI 1.8 and
        # remdesivir at 0.9 by selecting an outlier; the median replaced it.
        # The sentence describing the number must describe the number.
        basis = (f" Across {n_studies} studies the index ranges "
                 f"{selectivity_index:.1f} to {si_max:.1f}; the median of the "
                 f"exact measurements is quoted, and cell line and protocol "
                 f"differ between them.")
    else:
        basis = (" This rests on a SINGLE paired measurement, so it indicates "
                 "rather than establishes the compound's selectivity.")

    if selectivity_index < SI_THRESHOLD:
        return caveat(
            "selectivity_cytotoxic",
            f"selectivity index {selectivity_index:.1f}, below the "
            f"threshold of {SI_THRESHOLD:g}. The apparent activity may reflect "
            f"cell death rather than antiviral effect." + basis
        )
    return caveat(
        "selectivity_selective",
        f"selectivity index {selectivity_index:.1f} "
        f"(CC50/EC50, both measured in the same study)." + basis
    )


def evidence_note(tier: int | None, detection_method: str | None = None) -> str:
    """How the underlying fact was established."""
    if detection_method and detection_method.startswith(("MI:1314", "MI:0400",
                                                         "MI:0676", "MI:0007")):
        return caveat(
            "proximity_evidence",
            "this interaction was detected by proximity "
            "labelling or affinity co-purification, which shows the proteins are "
            "in the same complex neighbourhood rather than that they bind. "
            "Only 1.3% of virus-host edges in this graph come from a binary "
            "binding assay."
        )
    return caveat("evidence_tier", {
        1: "Curated experimental result.",
        2: "Curated but inferred rather than directly measured.",
        3: "Computationally derived, not experimentally asserted.",
        4: "Text-mined; not curated.",
    }.get(tier or 0, "Provenance not recorded."))
