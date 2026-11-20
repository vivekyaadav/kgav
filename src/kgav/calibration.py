"""What each reasoning channel is measured to be worth.

THE PROBLEM THIS SOLVES. baselines.combine ranks a drug by its BEST
percentile across metapaths, which was right when the alternative was summing
incomparable scales, and is wrong once some channels carry no signal. Max over
k channels is biased upward by k: a drug sitting at the top of a pure-noise
pool scores 1.0 and outranks one with a real route. Measured on the
single-screen protocol, COMBINED scores 0.516 while the best single channel
scores 0.523 -- pooling made the ranking worse, which is the arithmetic of max
over eight channels whose intervals all span 0.5.

WHAT A CHANNEL HAS TO EARN. Not a point AUC above 0.5 -- half of pure noise
clears that -- but a 95% interval that EXCLUDES 0.5. On the current graph no
channel clears it:

    M1 0.513 [0.452, 0.574]    M5 0.516 [0.455, 0.577]
    M2 0.492 [0.432, 0.552]    M6 0.502 [0.442, 0.562]
    M3 0.523 [0.462, 0.584]    M7 0.499 [0.439, 0.559]
    M4 0.507 [0.446, 0.568]    degree 0.582 [0.520, 0.644]

So a calibrated combine over this release returns NOTHING, and that is the
correct output. An agent whose channels are all at chance should say it cannot
rank candidates, not emit a confident order derived from noise. The ranking
being empty is a finding; a plausible-looking ranking built from these numbers
would be the failure.

WHY IT IS LOADED, NOT HARDCODED. A weight written into the source is a fitted
parameter with no record of what fitted it. These come from a results file that
carries its own protocol, label source, n, and whether the publication-disjoint
control was applied -- so a channel's weight can always be traced to the run
that earned it, and a weight measured under a confounded protocol can be
refused. Four protocols produced four different AUCs for M1 tonight (0.823,
0.726, 0.660, 0.499); a number with no protocol attached is not a measurement.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

# A channel must beat chance with its interval, not its point estimate.
CHANCE = 0.5

# AN AUC OVER UNREACHED POSITIVES IS NOT A MEASUREMENT OF THE CHANNEL.
# evaluate_against_negatives imputes a floor score for every compound a
# scorer does not reach, which is the right choice -- a scorer evaluated only
# on what it happens to reach chooses its own test set. But it means a
# channel's AUC is dominated by that floor when its coverage is low, and the
# interval is computed from n_pos rather than from the number it actually
# reached. M1 reads 0.513 [0.453, 0.574] on 89 positives while reaching
# 3.4% of them: THREE compounds. The figure describes the imputation, and the
# interval implies far more evidence than three compounds provide.
#
# 10 matches hard_negatives.MIN_CLASS_SIZE, which is the minimum class size
# that script requires before evaluating a virus at all. The same number is
# the minimum for evaluating one channel.
MIN_REACHED_POSITIVES = 10

# Protocols whose numbers must never become weights, with the reason. Checked
# by name so a results file produced under one is refused rather than silently
# trusted -- the whole sequence in docs/evaluation.md is numbers that looked
# usable until the next control ran.
REFUSED_PROTOCOLS: dict[str, str] = {
    "selectivity-positives-only":
        "filters positives to compounds with verified selectivity data and "
        "leaves negatives whole: a 17x coverage asymmetry, which is what "
        "produced the 0.823 headline",
}


@dataclass
class ChannelCalibration:
    """One metapath's measured performance under one protocol."""
    name: str
    auc: float | None
    ci_lo: float | None
    ci_hi: float | None
    n_pos: int
    n_neg: int
    coverage_pos: float = 0.0
    coverage_neg: float = 0.0

    @property
    def n_reached_pos(self) -> int:
        """Positives this channel actually found a path to."""
        return round(self.coverage_pos * self.n_pos)

    @property
    def evaluable(self) -> bool:
        """Whether the AUC describes the channel or the imputation."""
        return (self.auc is not None
                and self.n_reached_pos >= MIN_REACHED_POSITIVES)

    @property
    def discriminates(self) -> bool:
        """Whether the interval clears chance.

        Deliberately strict. A point AUC above 0.5 is what half of pure noise
        does, and every channel on this release has a point estimate above or
        near 0.5 while none has an interval that excludes it.
        """
        return (self.evaluable and self.ci_lo is not None
                and self.ci_lo > CHANCE)

    @property
    def weight(self) -> float:
        """Discrimination above chance, or zero if it has not been shown.

        auc - 0.5 rather than the AUC itself: a channel at 0.52 is not worth
        96% of one at 0.54, it is worth a fifth as much, because the useful
        quantity is the margin over chance.
        """
        return max(0.0, self.auc - CHANCE) if self.discriminates else 0.0

    def describe(self) -> str:
        """One line an agent can put in front of a reader."""
        if self.auc is None:
            return f"{self.name}: not evaluable ({self.n_pos}+/{self.n_neg}-)"
        if not self.evaluable:
            # Say what is missing rather than a number. Reporting 0.513 for a
            # channel that reached three compounds is the same error as a
            # hardcoded AUC: more confidence than the evidence carries.
            return (f"{self.name}: NOT EVALUABLE — reaches {self.n_reached_pos} "
                    f"of {self.n_pos:,} measured actives "
                    f"({self.coverage_pos:.1%}), below the {MIN_REACHED_POSITIVES} "
                    f"needed for an AUC to describe the channel rather than "
                    f"the imputed floor")
        verdict = ("discriminates" if self.discriminates
                   else "indistinguishable from chance")
        return (f"{self.name}: AUC {self.auc:.3f} "
                f"[{self.ci_lo:.3f}, {self.ci_hi:.3f}] on {self.n_pos:,} "
                f"measured active / {self.n_neg:,} measured inactive, "
                f"reaching {self.n_reached_pos} of them — {verdict}")


@dataclass
class Calibration:
    """Every channel's measured performance for one virus, with provenance."""
    virus: str
    protocol: str
    channels: dict[str, ChannelCalibration] = field(default_factory=dict)
    source_file: str = ""
    refused_reason: str = ""

    @property
    def usable(self) -> bool:
        return not self.refused_reason

    def discriminating(self) -> list[str]:
        if not self.usable:
            return []
        return [n for n, c in sorted(self.channels.items()) if c.discriminates]

    def evaluable(self) -> list[str]:
        """Channels whose AUC describes the channel rather than the floor."""
        return [n for n, c in sorted(self.channels.items()) if c.evaluable]

    def not_evaluable(self) -> list[str]:
        return [n for n, c in sorted(self.channels.items()) if not c.evaluable]

    def weights(self) -> dict[str, float]:
        if not self.usable:
            return {}
        return {n: c.weight for n, c in self.channels.items() if c.weight > 0}

    def can_rank(self) -> bool:
        """Whether ANY channel has earned a place in a ranking."""
        return bool(self.weights())

    def why_not(self) -> str:
        """The sentence an agent should say when it cannot rank. Plain, and
        specific about what was measured rather than apologetic."""
        if self.refused_reason:
            return (f"calibration from {self.source_file} was refused: "
                    f"{self.refused_reason}")
        tested = [c for c in self.channels.values() if c.evaluable]
        untested = self.not_evaluable()
        if not tested:
            return (f"NO channel was evaluable for {self.virus} under "
                    f"{self.protocol}: every one reaches fewer than "
                    f"{MIN_REACHED_POSITIVES} of the measured actives, so "
                    f"their AUCs describe the imputed floor rather than the "
                    f"channels. Not evaluated: {', '.join(untested)}. The "
                    f"graph can report evidence for a named compound; it has "
                    f"no measured basis for a shortlist.")
        best = max(tested, key=lambda c: c.auc)
        tail = (f" {len(untested)} channel(s) were not evaluable at all "
                f"({', '.join(untested)}), reaching fewer than "
                f"{MIN_REACHED_POSITIVES} measured actives each."
                if untested else "")
        return (f"no reasoning channel beats chance for {self.virus} under "
                f"{self.protocol}: of {len(tested)} evaluable, the strongest "
                f"is {best.describe()}.{tail} Ranking candidates would order "
                f"them by noise, so the graph can report evidence for a named "
                f"compound but not a shortlist.")

    def report(self) -> list[str]:
        lines = [f"calibration: {self.protocol} / {self.virus}"
                 + (f"  (from {self.source_file})" if self.source_file else "")]
        if self.refused_reason:
            lines.append(f"  REFUSED: {self.refused_reason}")
            return lines
        for _n, c in sorted(self.channels.items()):
            lines.append(f"  {c.describe()}")
        if self.can_rank():
            w = self.weights()
            lines.append("  ranking weights: " + ", ".join(
                f"{k} {v:.3f}" for k, v in sorted(w.items())))
        else:
            lines.append(f"  {self.why_not()}")
        return lines


def _channel(name: str, row: dict) -> ChannelCalibration:
    return ChannelCalibration(
        name=name,
        auc=row.get("auc"),
        ci_lo=row.get("auc_ci_lo"),
        ci_hi=row.get("auc_ci_hi"),
        n_pos=int(row.get("n_pos") or 0),
        n_neg=int(row.get("n_neg") or 0),
        coverage_pos=float(row.get("coverage_pos") or 0.0),
        coverage_neg=float(row.get("coverage_neg") or 0.0),
    )


SECTION_ALIASES = {"cross-sectional": "cross_sectional",
                   "post-2021": "temporal", "temporal": "temporal"}


def available_viruses(results: Path, section: str = "cross_sectional") -> list[str]:
    """The keys a results file actually holds, for a lookup that missed."""
    doc = json.loads(Path(results).read_text())
    sec = doc.get(SECTION_ALIASES.get(section, section)) or {}
    return sorted(k for k in sec if isinstance(sec[k], dict))


def load(results: Path, virus: str, section: str = "cross_sectional",
         protocol: str | None = None,
         exclude: tuple[str, ...] = ("COMBINED", "degree"),
         node_labels: dict[str, str] | None = None) -> Calibration:
    """Read one virus's channel calibration out of a hard_negatives results file.

    COMBINED and degree are excluded by default and for different reasons.
    COMBINED is the thing being built, so weighting it by its own score would
    be circular. degree is the NULL HYPOTHESIS -- it is what a ranking has to
    beat, not a channel to rank with, and on the single-screen protocol it is
    the only scorer whose interval clears chance, which is precisely why it
    must not be smuggled in as a reasoning route.
    """
    doc = json.loads(Path(results).read_text())
    sec = doc.get(SECTION_ALIASES.get(section, section)) or {}
    # A STALE PARAGRAPH ONCE STOOD HERE saying "Accept either, because a
    # caller holding a CURIE is not wrong and a silent empty result is the
    # worst possible answer to a key mismatch." It described an intention,
    # not this function: the paragraph below is the behaviour, and the two
    # read as one until you traced the code. run_baselines then passed a
    # CURIE without node_labels, got an empty calibration, and reported it as
    # "NO channel was evaluable ... Not evaluated: ." -- a verdict on
    # measurements it had never read. Callers must translate; the miss is
    # theirs to report with available_viruses().
    # hard_negatives.run keys its output by the virus's DISPLAY LABEL, taken
    # from node properties -- not by the taxon CURIE. A CURIE cannot be
    # translated by string surgery, because "NCBITaxon:2697049" and
    # "SARS-CoV-2" share nothing: it needs the mapping, which only a caller
    # holding the release has. So `node_labels` translates when supplied and
    # the lookup simply misses when it is not -- and available_viruses() exists
    # so the miss can be reported with what the file does hold, rather than as
    # an empty calibration that reads like an empty evaluation.
    wanted = (node_labels or {}).get(virus, virus)
    rows = sec.get(wanted)
    if rows is None:
        for key, val in sec.items():
            if isinstance(val, dict) and key.lower() == wanted.lower():
                rows = val
                break
    rows = rows or {}
    proto = protocol or doc.get("protocol") or _infer_protocol(doc, section)
    cal = Calibration(virus=virus, protocol=proto,
                      source_file=str(Path(results).name))
    for key, reason in REFUSED_PROTOCOLS.items():
        if key in proto:
            cal.refused_reason = reason
    for name, row in rows.items():
        if name in exclude or not isinstance(row, dict):
            continue
        cal.channels[name] = _channel(name, row)
    return cal


# The keys hard_negatives.py writes at the top level of its results file.
# Read from there rather than an "args" block: the flags are recorded flat,
# and guessing a nesting that does not exist is how the first version of this
# loader returned an empty calibration without saying why.
PROTOCOL_KEYS = ("selectivity_filter", "selectivity_applies",
                 "publication_disjoint", "label_source")


def _infer_protocol(doc: dict, section: str) -> str:
    """Name the protocol from what the results file recorded about itself.

    A file written before those keys existed records none of them, and the
    right answer then is "unknown" rather than a guess: an unnamed protocol is
    exactly what REFUSED_PROTOCOLS exists to catch, and a confounded run
    silently promoted to "looks fine" is the failure mode.
    """
    if not any(k in doc for k in PROTOCOL_KEYS):
        return f"{section}/unknown"
    bits = [section]
    src = doc.get("label_source")
    if src:
        bits.append(",".join(src) if isinstance(src, list) else str(src))
    else:
        bits.append("labels-pooled")
    sel = doc.get("selectivity_filter")
    applies = doc.get("selectivity_applies")
    if sel and sel != "all":
        bits.append(f"{sel}-{applies or 'positives'}")
        if applies != "both":
            bits.append("selectivity-positives-only")
    if doc.get("publication_disjoint"):
        bits.append("publication-disjoint")
    else:
        # Named explicitly. Without it M1 reads 0.660 instead of 0.499, and a
        # protocol string that simply omits the control reads as though the
        # control had passed.
        bits.append("NO-publication-disjoint")
    return "/".join(bits)


# ====================================================================
# THE PROSPECTIVE CALIBRATION
#
# A SECOND LOADER, BECAUSE IT IS A SECOND FILE SHAPE AND NOT A SECOND
# SECTION. `load` above reads hard_negatives.json, whose `temporal`
# section is a {virus: {channel: row}} mapping. The prospective study
# writes temporal_structure_scorer.json, which keys `cells` by TAXON
# CURIE and puts each scorer under its own name beside per-cell counts.
# Pointing `load` at it returns an empty calibration, which is exactly
# what happened: hard_negatives.json's `temporal` section is `{}` on this
# release -- the temporal protocol was skipped there -- so the agent's
# PROSPECT refusal took its no-measurement branch and told every asker
# "this release carries no measured prospective result" while
# temporal_structure_scorer.json recorded six cells clearing chance.
#
# The refusal was weaker than the finding, which is the safe direction and
# still wrong: a refusal that understates what has been measured denies
# the reader a result they are entitled to before being refused, and
# guards._temporal_sentence has had a branch for precisely that case
# ("...which clears chance") that nothing could ever reach.

#: Scorers in the prospective file that must not become ranking weights.
#: `degree` is the NULL HYPOTHESIS -- the thing a scorer has to beat --
#: for the same reason `load` excludes it. On this release it is unusable
#: in all ten cells anyway (the training graph reaches 0.3%-6.5% of the
#: compounds being ranked), but it is excluded by role and not by its
#: having come out badly.
TEMPORAL_EXCLUDE = ("degree",)

#: File key -> what to call the scorer in front of a reader. The agent
#: prints these inside a refusal, and "primary" names this study's
#: position in its own registration rather than anything the scorer does.
#: An unlisted key passes through unchanged, so a new scorer is readable
#: as its key rather than silently renamed to something it is not.
TEMPORAL_SCORER_NAMES = {
    "primary": "max-Tanimoto to pre-cutoff actives",
    "symmetric_k10": "the same with each candidate's 10 nearest withheld",
}


def available_temporal_viruses(results: Path) -> list[str]:
    """The taxon keys temporal_structure_scorer.json holds.

    The counterpart of `available_viruses`, so a lookup that misses can be
    reported with what the file does contain rather than as an empty
    calibration -- which reads like an empty evaluation.
    """
    doc = json.loads(Path(results).read_text())
    cells = doc.get("cells") or {}
    return sorted(k for k in cells if isinstance(cells[k], dict))


def _temporal_channel(name: str, row: dict, n_pos: int, n_neg: int
                      ) -> ChannelCalibration:
    """One prospective scorer, carrying its CORRECTED interval.

    THE CORRECTED ONE, NOT THE 95% ONE. `discriminates` tests
    `ci_lo > CHANCE`, and the registration
    (docs/prereg_temporal_structure_scorer.md) decides whether a cell
    clears at alpha = 0.05/10 with z = 2.8070, not at z = 1.96. Handing
    it `auc_ci_lo` would let the agent announce cells as clearing chance
    that the registered family size says do not -- the agent claiming
    more than the study it is quoting, from the same file.

    Checked by tests/test_temporal_calibration.py, which requires the
    channels this produces to reproduce the file's own `clears` list
    exactly rather than approximately.
    """
    return ChannelCalibration(
        name=TEMPORAL_SCORER_NAMES.get(name, name),
        auc=row.get("auc"),
        ci_lo=row.get("auc_ci_lo_corrected"),
        ci_hi=row.get("auc_ci_hi_corrected"),
        n_pos=int(row.get("n_pos") or n_pos or 0),
        n_neg=int(row.get("n_neg") or n_neg or 0),
        coverage_pos=float(row.get("coverage_pos") or 0.0),
        coverage_neg=float(row.get("coverage_neg") or 0.0),
    )


def temporal_protocol(doc: dict) -> str:
    """The protocol string the agent quotes, built from the file's own record.

    Named from what the run recorded rather than from a constant, so a
    split built at a different cutoff cannot be described with this one's
    wording. The agent prints this, so it has to be readable as well as
    exact.
    """
    bits = ["temporal"]
    for key, prefix in (("split", ""), ("cutoff", "cutoff-"),
                        ("undated_policy", "undated-"),
                        ("overlap_policy", "")):
        val = doc.get(key)
        if val not in (None, ""):
            bits.append(f"{prefix}{val}")
    return "/".join(bits)


def load_temporal(results: Path, virus: str,
                  exclude: tuple[str, ...] = TEMPORAL_EXCLUDE) -> Calibration:
    """One virus's prospective calibration from temporal_structure_scorer.json.

    `virus` is matched against the file's `cells` keys, which are taxon
    CURIEs, and then case-insensitively against each cell's recorded
    `name` -- so a caller holding either resolves, and a caller holding
    neither misses rather than receiving somebody else's numbers.

    A cell the protocol refused comes back with `refused_reason` set, so
    `usable` is False and the agent takes its no-measurement branch. That
    is the distinction this file exists to preserve: a refused cell is
    not a cell that scored 0.5.
    """
    path = Path(results)
    doc = json.loads(path.read_text())
    cells = doc.get("cells") or {}
    proto = temporal_protocol(doc)
    cal = Calibration(virus=virus, protocol=proto, source_file=str(path.name))

    row = cells.get(virus)
    if row is None:
        for key, val in cells.items():
            if not isinstance(val, dict):
                continue
            if (key.lower() == virus.lower()
                    or str(val.get("name", "")).lower() == virus.lower()):
                row = val
                break
    if not isinstance(row, dict):
        # NOT A REFUSAL AND NOT A ZERO. A virus this study never evaluated
        # has no prospective calibration, and `usable` False with no
        # reason is how `load` reports the same thing.
        cal.refused_reason = (
            f"{path.name} holds no prospective cell for {virus}; it "
            f"evaluated {len(cells)} of "
            f"{doc.get('n_cells_total', len(cells))} candidate cells")
        return cal

    # THE WHOLE FILE'S REFUSAL, if the run itself refused. write_refusal
    # records the reasons under this key in place of results, and a
    # refusal record still carries the cells it had, so this has to be
    # checked before reading any of them.
    if doc.get("refused"):
        cal.refused_reason = ("the prospective run refused: "
                              + "; ".join(doc["refused"]))
        return cal

    n_pos = int(row.get("n_test_positives") or 0)
    n_neg = int(row.get("n_test_negatives") or 0)
    for name, sub in row.items():
        if name in exclude or not isinstance(sub, dict):
            continue
        if "auc" not in sub:
            continue                       # counts blocks, not scorers
        if sub.get("measurable") is False or sub.get("refusal"):
            continue                       # declined, not measured at 0.5
        cal.channels[name] = _temporal_channel(name, sub, n_pos, n_neg)
    if not cal.channels:
        cal.refused_reason = (
            f"no scorer in {path.name} produced a measurable AUC for "
            f"{virus}")
    return cal
