"""Parameterised graph tools.

WHY NOT LET THE MODEL WRITE QUERIES. A language model asked to generate
traversal code against a 21-predicate schema fails quietly: a malformed
traversal returns an empty result, and an empty result is indistinguishable
from "no evidence exists". The model then reports the latter. Every tool here
is a fixed function with typed parameters, so a failure is a failure rather
than a confident absence.

WHAT EACH TOOL RETURNS. A structured record plus a verbalisation — a sentence
rendering with provenance inline. The model reads text, so leaving
verbalisation to the model means it invents phrasing for facts it half
understands. Rendering it here keeps the wording tied to the data.

BOUNDED BY CONSTRUCTION. Every tool caps its output. An unbounded neighbour
expansion from a well-connected protein returns thousands of rows and fills the
context window with ubiquitin.
"""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from kgav import virus_register
from kgav.agent import guards
from kgav.baselines import qualifier_ok

# THE VERDICT COMES FROM ONE PLACE (fix-list L1). The agent used to
# render every HAS_ANTIVIRAL_ACTIVITY_AGAINST edge as "shows activity
# against", dropping `relation` entirely -- and 20,239 of the 26,017
# edges carrying that predicate, 77.8%, are measurements of INACTIVITY.
# A second copy of the threshold logic here would be the two-
# implementations failure this repository keeps finding, so classify()
# is imported and the predicate name is not trusted.
from kgav.labels import (
    DEFAULT_ACTIVE_BELOW_NM,
    DEFAULT_INACTIVE_ABOVE_NM,
    classify,
)

# The significance threshold, imported rather than restated: the
# landscape already declares one and two different alphas in one
# answer would be indefensible.
from kgav.landscape import ALPHA as LANDSCAPE_ALPHA
from kgav.schema import load_schema

MAX_NEIGHBOURS = 50
MAX_PATHS = 20
MAX_SUBGRAPH_EDGES = 100
# Pathways LISTED by pathways(), not pathways considered. Every one found is
# counted and the overflow is reported with how many of it carry a dependency
# factor, so the cap bounds the output without hiding the census.
MAX_PATHWAYS = 25
#: Rows in the by-density list. Shorter than MAX_PATHWAYS because it
#: is a second view of the same census and the brief is bounded.
MAX_RATIO_ROWS = 10
# Dependency-factor names printed per pathway. The count is always stated in
# full; only the naming is bounded.
MAX_DEPENDENCY_NAMES = 6
# Quantile bins for the pathway-size background. Quantiles, not equal width,
# for the reason degree_bins uses them: pathway size is heavy-tailed (median
# 4, p99 113 on v0.1) and equal-width bins put nearly everything in the first.
#: The smallest pooled sample a background rate may be estimated from.
#: A binomial rate has standard error sqrt(p(1-p)/n), worst at p = 0.5,
#: so n = 100 holds that error at or under 0.05 whatever the true rate.
#: Pooling stops as soon as the window reaches it, so a size with enough
#: pathways of its own is never smoothed against its neighbours.
#:
#: It replaces PATHWAY_SIZE_BINS = 10, which split tied sizes across
#: quantile bins -- see _pathway_background for what that cost.
MIN_BACKGROUND_SAMPLES = 100


@dataclass
class ToolResult:
    """Every tool returns this shape.

    `warnings` are not decoration: they carry the cell-context and
    host-directed caveats that make a path interpretable, and the verification
    stage checks they survived into the answer.
    """
    ok: bool
    kind: str
    data: Any = None
    verbalised: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    edge_ids: list[str] = field(default_factory=list)
    note: str = ""

    def to_json(self) -> str:
        return json.dumps({"ok": self.ok, "kind": self.kind, "data": self.data,
                           "text": self.verbalised, "warnings": self.warnings,
                           "edge_ids": self.edge_ids, "note": self.note},
                          indent=2, default=str)


class GraphTools:
    """The graph, loaded once, exposed as a fixed set of callable tools."""

    def __init__(self, release: Path, name_index: dict[str, str] | None = None,
                 calibration=None, landscape=None,
                 reach: dict[str, int] | None = None):
        # MEASURED REACH PER VIRUS, {display label: held-out positives the
        # metapaths reach}, or None for "not measured".
        #
        # The one fact in VirusRoutes that cannot come from the graph.
        # labels.measured_reach is its single owner and returns None
        # unless provenance says metapath_comparison.json is current for
        # this release, so a stale file makes reach UNKNOWN rather than
        # ZERO -- opposite claims, and the stale-figure mistake this
        # project has made twice already.
        #
        # A virus absent from the map also resolves to None. The map is
        # keyed by the label metapath_comparison.json uses, which is the
        # same display label _label() returns, so the lookup needs no
        # translation -- unlike calibration, which is keyed by whatever
        # the caller passed and needed a node_labels map to fix.
        self._reach = dict(reach or {})
        # ONE LANDSCAPE FOR THE WHOLE RELEASE, queried per virus.
        #
        # Unlike calibration there is no per-virus mapping to mis-key: the
        # landscape is a single measurement over every taxon and
        # Landscape.for_virus returns that taxon's row or None, never a
        # neighbour's. None here is not neutral -- guards.pathway_caveat
        # renders it as "not measured on this release", which is the honest
        # branch and still appears in the answer. The one thing an answer may
        # not do is stay silent about it.
        self.landscape = landscape
        # Memo for whole-graph scans that do not depend on the question.
        self._cache: dict[tuple, Any] = {}
        # ONE CALIBRATION PER VIRUS, NOT ONE PER INSTANCE.
        #
        # This was a single `self.calibration`, handed to every caveat
        # regardless of which virus was asked about. explain() and triage()
        # both TAKE a virus argument and neither used it to pick the
        # calibration, so an agent constructed with SARS-CoV-2's numbers and
        # asked about Zika would state Zika's direct-acting route "separates
        # measured actives from measured inactives at AUC 0.675" -- SARS-CoV-2's
        # figure, on the discriminates branch, which is the most confident
        # sentence this module can emit. Zika's own M1 is 0.511, at chance.
        #
        # Latent while every virus was at chance anyway. Live the moment a
        # second family arrived with asymmetric calibration, and the mismatch
        # was always detectable: Calibration carries its own `.virus`, and
        # nothing compared it to the virus being answered about.
        #
        # Accepts a single Calibration (for one virus) or a {virus:
        # Calibration} mapping. A virus with no calibration resolves to None,
        # which routes the caveats to their honest "uncalibrated" branch --
        # never to another virus's numbers.
        self._calibrations: dict[str, object] = {}
        self._pending_calibration = calibration
        # METAPATH CONSTRAINTS, READ FROM THE SCHEMA RATHER THAN RESTATED.
        #
        # _find_paths enumerates M1 and M2 as its own loops. M1 applied NO
        # constraints and M2 a hand-copied one, so a constraint added to the
        # schema reached dwpc_scores -- which passes mp.get("constraints") --
        # and never reached here. When M1 gained
        # INHIBITS.target_resolution=chain_from_assay_description, the scorer
        # stopped walking parent_unresolved edges and the agent carried on
        # explaining routes through them: it could narrate a mechanism the
        # ranking refuses. Same shape as the agent holding one calibration
        # for every virus -- a fact declared in one place and consumed in
        # one of two.
        _schema = load_schema()
        # THE PREDICATES THAT ARE THE ANSWER SHEET, read from the schema
        # and not restated here. There are FOUR -- the external audit
        # named one -- and a list kept in this file would drift the first
        # time a predicate was declared held_out in the schema alone,
        # which is the whole defect class this project keeps finding.
        self._held_out = frozenset(_schema.held_out_predicates())
        self._constraints: dict[str, dict] = {
            name: (mp.get("constraints") or {})
            for name, mp in _schema.metapaths.items()}
        # THE SAME CAP THE SCORER AND THE LANDSCAPE USE, read from the schema
        # rather than restated. pathways() must exclude the same oversized
        # pathways pathway_landscape.py excludes, or the agent would describe
        # a route through a 1,389-member pathway that the measured null it
        # quotes never considered.
        self._max_pathway_members: int | None = (
            _schema.retrieval_limits.get("max_pathway_members"))
        self.nodes: dict[str, dict] = {}
        self.edges: list[dict] = []
        self.out: dict[tuple[str, str], list[dict]] = defaultdict(list)
        self.inv: dict[tuple[str, str], list[dict]] = defaultdict(list)
        self.by_id: dict[str, dict] = {}
        # lowercase display name -> node id, for entity resolution
        self.names: dict[str, str] = dict(name_index or {})
        self._register_names: dict[str, str] = {}

        release = Path(release)

        # Compound names are not in the graph: nodes are keyed on InChIKey and
        # the ChEMBL layer stores only chembl_id, because path scoring never
        # needs a name. An agent has a human on the other end who types
        # "remdesivir", so the release carries a prebuilt name index rather
        # than the agent opening a 25 GB database at runtime.
        self.preferred_names: dict[str, str] = {}
        # EVERY VIRUS THE REGISTER DECLARES IS RESOLVABLE, whatever NAMES.json
        # happens to contain. That file carries the register's SHORT name and
        # nothing else, so "dengue virus 2" matched nothing, no taxon
        # resolved, and the question was answered about the default virus
        # under a caveat asserting that no virus had been named. Seeded
        # before NAMES.json so a release that does carry a better name still
        # wins; a name for a taxon absent from this graph is dropped, because
        # resolving to a node that is not here would trade a wrong answer for
        # a crash.
        # The doc is kept as well as the derived names, because the
        # label-matching map needs label_names() off the same document --
        # loading the register twice would let the two disagree.
        self._register_doc: dict = {}
        try:
            self._register_doc = virus_register.load()
            for key, curie in virus_register.search_names(
                    self._register_doc).items():
                self._register_names[key] = curie
        except (OSError, ValueError, KeyError):
            # A malformed register must not stop the agent answering about
            # the compounds it can still resolve.
            self._register_names = {}
            self._register_doc = {}
        names_file = release / "NAMES.json"
        if names_file.exists():
            self.names.update(json.loads(names_file.read_text()))
        # Assembly merges the per-study selectivity edges into one, keeping the
        # worst value. That is right for filtering and wrong for display, so
        # the spread is read from the pre-assembly layer.
        self.si_range: dict[tuple[str, str], list[float]] = {}
        self.si_cell: dict[tuple[str, str], str] = {}
        si_layer = release.parent / (release.name + "-selectivity") / "edges.jsonl"
        if si_layer.exists():
            for line in si_layer.read_text().splitlines():
                if not line.strip():
                    continue
                e = json.loads(line)
                si = (e.get("qualifiers") or {}).get("selectivity_index")
                if si is not None:
                    self.si_range.setdefault((e["subject"], e["object"]), []).append(si)
                    cl = (e.get("qualifiers") or {}).get("cell_line")
                    if cl:
                        self.si_cell[(e["subject"], e["object"])] = cl

        pref_file = release / "PREF_NAMES.json"
        if pref_file.exists():
            self.preferred_names.update(json.loads(pref_file.read_text()))

        for line in (release / "nodes.jsonl").read_text().splitlines():
            if line.strip():
                n = json.loads(line)
                self.nodes[n["id"]] = n
                label = (n["properties"].get("gene_symbol")
                         or n["properties"].get("label"))
                if label:
                    self.names.setdefault(str(label).lower(), n["id"])

        # Register names LAST and only for taxa this release holds. setdefault
        # so a name the release supplied keeps precedence, and the presence
        # check so a question naming a virus that was declared but not built
        # resolves to nothing -- which the assumed-virus caveat then reports
        # honestly -- rather than to a node id with no node behind it.
        for key, curie in self._register_names.items():
            if curie in self.nodes:
                self.names.setdefault(key, curie)

        for i, line in enumerate(
                (release / "edges.jsonl").read_text().splitlines()):
            if not line.strip():
                continue
            e = json.loads(line)
            e["_id"] = e.get("edge_id") or f"E{i}"
            self.edges.append(e)
            self.by_id[e["_id"]] = e
            self.out[(e["predicate"], e["subject"])].append(e)
            self.inv[(e["predicate"], e["object"])].append(e)

        # After the nodes, because reconciling a calibration's key needs
        # _label to translate a taxon CURIE into the display label the results
        # file is keyed by. calibration.load stores whatever key the CALLER
        # passed -- "SARS-CoV-2" from kgav_ask, "NCBITaxon:2697049" from a
        # caller holding the release -- so neither form can be assumed.
        self._index_calibrations(self._pending_calibration)
        del self._pending_calibration

    def family_of(self, virus: str) -> str | None:
        """The viral family of a taxon, from the graph's own OrganismTaxon.

        ingest_spine writes `family` from config/viruses.yaml, so this is the
        register's answer rather than an inference from the name.
        """
        n = self.nodes.get(virus) or {}
        fam = (n.get("properties") or {}).get("family")
        return str(fam) if fam else None

    # --------------------------------------------------- calibration lookup
    def _same_virus(self, a: str, b: str) -> bool:
        """Whether two virus keys name one virus, across CURIE and label.

        "NCBITaxon:2697049" and "SARS-CoV-2" share no characters, so this
        needs the release's own node labels rather than string surgery --
        which is the same reason calibration.load takes a node_labels map.
        """
        if a == b:
            return True
        for x, y in ((a, b), (b, a)):
            if x in self.nodes and self._label(x) == y:
                return True
        return False

    def _index_calibrations(self, calibration) -> None:
        if calibration is None:
            return
        if hasattr(calibration, "channels"):          # one Calibration
            key = getattr(calibration, "virus", "") or ""
            self._calibrations[key] = calibration
            return
        for key, cal in dict(calibration).items():
            # A MIS-KEYED MAP IS A SILENT WRONG-FIGURE GENERATOR, which is the
            # bug this method exists to end, so it raises rather than being
            # tolerated. Keying SARS-CoV-2's calibration under Zika would
            # reintroduce the exact failure by a different route.
            own = getattr(cal, "virus", "") or ""
            if own and not self._same_virus(key, own):
                raise ValueError(
                    f"calibration keyed under {key!r} was measured on {own!r}; "
                    f"a calibration must be keyed by the virus it measured")
            self._calibrations[key] = cal

    def routes_for(self, virus: str) -> guards.VirusRoutes:
        """Which feature channels this release holds for one virus.

        COUNTED FROM THE GRAPH BEING SERVED, not from a results file
        or a note. A virus whose routes all disappeared in a rebuild
        must stop being answerable in the same run, and a figure
        carried from elsewhere cannot do that.

        The label count is held-out positives. It is here so the
        refusal can say what exists, and is never an input to an
        answer: quoting a label as evidence is the leak the rest of
        this module exists to stop.

        ONE EXCEPTION TO THE PARAGRAPH ABOVE, added 2026-11-02:
        reached_labels comes from the reach map, not from the graph,
        because it cannot be counted from edges at all -- walking a
        metapath against a label set is what produces it. The guarantee
        the first paragraph protects still holds, by a different route:
        the map is refused unless provenance says it was measured on
        the release being served, so a figure from a previous graph is
        reported as unmeasured rather than quoted. Without it the agent
        cannot tell `evaluable` from ROUTES REACH NO LABEL and tells
        five viruses their reliability is unknown when it is zero.
        """
        inhibits = crispr = ppi = labels = 0
        for e in self.edges:
            pred, subj, obj = e["predicate"], e["subject"], e["object"]
            if pred == "HAS_ANTIVIRAL_ACTIVITY_AGAINST":
                labels += obj == virus
            elif pred == "HOST_FACTOR_FOR":
                crispr += (obj == virus
                           and (e.get("qualifiers") or {}).get("direction")
                           == "dependency")
            elif pred == "INHIBITS":
                inhibits += self._taxon_of(obj) == virus
            elif pred == "PHYSICALLY_INTERACTS_WITH":
                # EXACTLY ONE ENDPOINT NON-VIRAL, and the viral one is
                # this virus. The predicate is shared with STRING, so a
                # host-host edge has to be excluded too.
                #
                # The previous test was `(a == virus) != (b == virus)`
                # over `_taxon_of`, which returns None for a host
                # protein. That excludes host-host (None, None) and
                # same-virus viral pairs (virus, virus) -- and COUNTS a
                # cross-virus viral pair, because (virus, otherVirus)
                # gives True != False. The comment above it claimed
                # viral-viral was excluded "by the same check"; it was
                # excluded only when both proteins belonged to the virus
                # being asked about. Found by external audit,
                # 2026-11-02. The register records SARS-CoV's 73
                # VirHostNet interactions as viral-viral, so this
                # inflated a route count that decides whether the agent
                # refuses to rank at all.
                a, b = self._taxon_of(subj), self._taxon_of(obj)
                one_viral = (a is None) != (b is None)
                ppi += one_viral and virus in (a, b)
        name = self._label(virus)
        return guards.VirusRoutes(label=name, n_labels=labels,
                                  inhibits=inhibits, crispr=crispr, ppi=ppi,
                                  reached_labels=self._reach.get(name))

    def _taxon_of(self, node_id: str) -> str | None:
        """The virus a VIRAL protein belongs to, else None."""
        props = (self.nodes.get(node_id) or {}).get("properties") or {}
        return str(props["taxon_id"]) if props.get("is_viral") \
            and props.get("taxon_id") else None

    def calibration_for(self, virus: str):
        """The calibration measured ON this virus, or None.

        None rather than a fallback: every caveat has an honest uncalibrated
        branch, and "its reliability is unknown" is weaker than a figure and
        true, where another virus's figure is neither.
        """
        for key, cal in self._calibrations.items():
            if self._same_virus(key, virus):
                return cal
        return None

    # ------------------------------------------------------------- helpers
    def _label(self, node_id: str) -> str:
        """Display name, most specific first.

        protein_family precedes gene_symbol deliberately: every one of the 16
        mature peptides in the coronavirus replicase shares the gene symbol
        "rep", so preferring the symbol discards exactly the target
        specificity the polyprotein expansion existed to recover. A route
        reading "inhibits rep" is useless; "inhibits nsp5" is the finding.
        """
        n = self.nodes.get(node_id)
        if not n:
            return node_id
        p = n["properties"]
        if n["class"] == "SmallMolecule":
            # Reverse the name index: a human reads "nirmatrelvir", not an
            # accession.
            # Prefer the canonical name. A reverse lookup over the synonym
            # index returns whichever alias happened to be inserted last --
            # "pf07321332" for nirmatrelvir, "nsc-187208" for chloroquine.
            name = self.preferred_names.get(node_id)
            if name:
                return name
            # Fall back to a label on the node before the bare accession. A
            # release without PREF_NAMES.json should still render readable
            # names rather than "CHEMBL1".
            return str(p.get("label") or p.get("chembl_id") or node_id)
        if p.get("mature_peptide") and p.get("protein_family"):
            return f"{p['protein_family']} ({p['mature_peptide']})"
        return str(p.get("protein_family") or p.get("gene_symbol")
                   or p.get("label") or node_id)

    @staticmethod
    def _potency_note(nm: float) -> str:
        """Qualify a potency value in words.

        Without this, "IC50 75 nM" and "IC50 7280 nM" read as equivalent
        claims. The second is 7.3 uM -- above the threshold this project uses
        for "not pursued" -- and in the case that prompted this, it is
        chloroquine reported against a protease it does not inhibit, from a
        broad 2020 screening panel.
        """
        if nm <= 100:
            return ", potent"
        if nm <= 1000:
            return ", moderate"
        return ", WEAK — above 1 uM and probably not a real mechanism"

    def _activity_sentence(self, e: dict, s: str, o: str) -> str:
        """One HAS_ANTIVIRAL_ACTIVITY_AGAINST edge, saying what it measured.

        THE PREDICATE NAME IS WRONG FOR MOST OF ITS INSTANCES (fix-list
        L1). On release v0.1, `classify` over the 26,278 edges carrying
        it reads:

            inactive    20,239   77.8%
            active       5,525   21.0%
            undecidable    515    2.0%

        and this template rendered all of them as "shows activity
        against {o} (EC50 ... nM)", dropping `relation`. So an edge
        reading "EC50 > 50000 nM" -- a screen reporting that it tested
        to 50 uM and saw nothing -- was presented to the user as
        evidence of activity, with the bound's own number quoted as if
        it were a measured potency. That is the single most misleading
        sentence the agent could produce from this graph, and it could
        produce it about three quarters of the label set.

        WHY THE PREDICATE IS NOT RENAMED INSTEAD. The obvious fix is to
        route inactives to MEASURED_INACTIVE_AGAINST at ingest, as the
        protein layer does. At the ORGANISM level that would be wrong:
        `inactive_above_nm` is a protocol parameter, varied by
        `hard_negatives.py --inactive-above-nm` and recorded in its
        output, so routing by verdict would bake an evaluation threshold
        into the graph's topology and the predicate would stop matching
        the verdict the moment the threshold moved. The predicate means
        "this compound was measured against this virus"; what the
        measurement SAYS is `classify`'s job, and it stays there.

        The thresholds are the module defaults and are named in the
        sentence, because a verdict quoted without the cut it was taken
        at is not reproducible by the reader.
        """
        q = e.get("qualifiers") or {}
        verdict = classify(e, DEFAULT_INACTIVE_ABOVE_NM,
                           DEFAULT_ACTIVE_BELOW_NM)
        relation = q.get("relation") or "="
        value = next((float(q[f]) for f in ("ec50_nm", "ic50_nm")
                      if q.get(f) is not None), None)
        # The measurement, with its relation, never without.
        if value is None:
            shown = ""
        else:
            unit = "uM" if value >= 1000 else "nM"
            num = value / 1000 if unit == "uM" else value
            shown = f" ({relation} {num:,.4g} {unit})".replace("(= ", "(")

        if verdict == "inactive":
            return (f"{s} was MEASURED AGAINST {o} AND FOUND INACTIVE"
                    f"{shown}"
                    f" -- not evidence of activity")
        if verdict == "active":
            # EXPLICIT, not a conditional expression spanning a `+`:
            # "A + B if c else C" parses as "(A + B) if c else C", which
            # is right here and reads as though it were not.
            if value is None:
                return f"{s} shows activity against {o}"
            return (f"{s} shows activity against {o}{shown}"
                    + self._potency_note(value))
        # TWO WAYS TO DECIDE NOTHING, and they are not the same thing to
        # a reader. A bound below the inactivity floor is a measurement
        # that fails to settle the question; a missing value is no
        # measurement. Saying "a bound below the floor" about an edge
        # carrying no bound was the first version of this branch.
        if value is None:
            return (f"{s} was measured against {o} and no potency value "
                    f"was recorded, so the edge decides nothing")
        return (f"{s} was measured against {o} and the result decides "
                f"nothing{shown} -- a bound below the "
                f"{DEFAULT_INACTIVE_ABOVE_NM / 1000:,.0f} uM inactivity "
                f"floor proves neither activity nor its absence")

    def _verbalise(self, e: dict) -> str:
        """One edge as a sentence with its provenance inline.

        Templates live here rather than in the model, so a fact is always
        described the same way and always with its source attached.
        """
        s, o = self._label(e["subject"]), self._label(e["object"])
        q = e.get("qualifiers") or {}
        src = (e.get("primary_knowledge_source") or "?").replace("infores:", "")
        # Cap citations: one nirmatrelvir edge carries 51 PMIDs, which is
        # unreadable and would consume the context window.
        all_pubs = e.get("publications") or []
        pubs = ", ".join(all_pubs[:3])
        if len(all_pubs) > 3:
            pubs += f" +{len(all_pubs) - 3} more"
        support = e.get("support_count", 1)

        templates = {
            "INHIBITS": lambda: (
                f"{s} inhibits the viral protein {o}"
                + (f" (IC50 {q['ic50_nm']:.0f} nM{self._potency_note(q['ic50_nm'])})"
                   if q.get("ic50_nm") else "")
                + (f", assayed against the {q['domain']} site" if q.get("domain") else "")),
            "TARGETS": lambda: (
                f"{s} engages the host protein {o}"
                # THREE STATES, NOT TWO. "contested" means two curated
                # sources assert opposite directions, which is a different
                # thing from nobody stating one -- and rendering it through
                # the first branch produced "as an contested".
                + (" (sources disagree on whether it inhibits or activates "
                   "this target)" if q.get("direction") == "contested"
                   else f" as an {q['direction']}"
                   if q.get("direction") not in (None, "unknown")
                   else " (mechanism not stated: potency measured, not direction)")),
            "HOST_FACTOR_FOR": lambda: (
                f"{s} is a host {q.get('direction', 'factor')} factor for {o}"
                + (f", found in {support} independent CRISPR screens" if support > 1
                   else ", found in 1 CRISPR screen")
                + (f" ({q['cell_line']})" if q.get("cell_line") else "")),
            "PHYSICALLY_INTERACTS_WITH": lambda: (
                f"{s} and {o} were detected in the same complex or neighbourhood"),
            "ENCODED_BY": lambda: f"{s} is encoded by {o}",
            "BELONGS_TO": lambda: f"{o} carries the gene {s}",
            "PARTICIPATES_IN": lambda: f"{s} participates in {o}",
            "HAS_ANTIVIRAL_ACTIVITY_AGAINST": lambda: (
                self._activity_sentence(e, s, o)),
            "CHEMICALLY_SIMILAR_TO": lambda: (
                f"{s} is structurally similar to {o} "
                f"(Tanimoto {q.get('tanimoto_ecfp4', 0):.2f})"),
            "CAUSES": lambda: f"{s} causes {o}",
        }
        body = templates.get(e["predicate"], lambda: f"{s} {e['predicate']} {o}")()
        # MARKED HERE, ONCE, because this is the single place an edge
        # becomes a sentence. Marking in each tool would have left
        # whichever tool was written next unmarked -- and the audit
        # found two of the three that already were.
        if e["predicate"] in self._held_out:
            body = "[HELD-OUT LABEL] " + body
        cite = f" [{src}"
        if pubs:
            cite += f", {pubs}"
        cite += f", {e.get('first_asserted_date', '?')[:4]}; {e['_id']}]"
        return body + cite

    def _label_name_map(self) -> dict[str, set[str]]:
        """{name: taxa that name indicates} for matching against a LABEL.

        virus_register.label_names plus the display label of every taxon
        node the release holds, so a pathway naming a virus is recognised
        whatever the register calls it. Cached: it is the same for every
        row of every call.

        MIN_SEARCH_NAME IS IMPORTED, not re-spelled. This was a literal
        `len(n) >= 4` here, a second copy of a threshold whose owner is
        virus_register -- in the module whose family_of docstring is an
        account of what happened the last time one rule had five copies
        and three were wrong.
        """
        cached = self._cache.get(("label_names",))
        if cached is not None:
            return cached
        out: dict[str, set[str]] = defaultdict(set)
        for curie, names in virus_register.label_names(
                self._register_doc).items():
            if curie not in self.nodes:
                continue
            for n in names:
                out[n].add(curie)
        for nid, node in self.nodes.items():
            if node.get("class") != "OrganismTaxon":
                continue
            lab = str((node.get("properties") or {}).get("label") or "").lower()
            if len(lab) >= virus_register.MIN_SEARCH_NAME:
                out[lab].add(nid)
        out = dict(out)
        self._cache[("label_names",)] = out
        return out

    def _names_virus(self, label: str, virus: str) -> bool:
        """Whether this text label names THIS virus, by longest match.

        THE BUG THIS REPLACES. The test was `any(name in label)` over one
        virus's names, and "sars-cov" is a substring of "sars-cov-2" --
        as is "severe acute respiratory syndrome coronavirus" of the same
        name plus " 2". So on release v0.1 six Reactome pathways built
        from SARS-CoV-2 literature were flagged as naming SARS-CoV:

            Replication of the SARS-CoV-2 genome
            SARS-CoV-2 modulates autophagy
            SARS-CoV-2 modulates host translation machinery
            SARS-CoV-2 targets PDZ proteins in cell-cell junction
            SARS-CoV-2 targets host intracellular signalling and regulation
            SARS-CoV-2 activates/modulates innate and adaptive immune responses

        A circularity warning on the wrong virus is worse than none: it
        tells a reader the overlap for SARS-CoV was decided in advance
        when nothing about SARS-CoV was involved.

        SHADOWED BY SPAN, NOT GLOBALLY. A match is discarded only when a
        LONGER name covers the same starting offset, so

            "SARS-CoV and SARS-CoV-2 entry"

        still names both: the eight-character match at offset 0 is not
        shadowed by the ten-character one at offset 13. Suppressing
        whenever any longer name matched anywhere would lose the first.

        A name shared by several taxa keeps all of them -- for dengue it
        must, see virus_register.label_names.
        """
        low = (label or "").lower()
        spans: dict[tuple[int, int], set[str]] = defaultdict(set)
        for name, taxa in self._label_name_map().items():
            i = low.find(name)
            while i != -1:
                spans[(i, len(name))] |= taxa
                i = low.find(name, i + 1)
        if not spans:
            return False
        starts: dict[int, int] = {}
        for start, length in spans:
            starts[start] = max(starts.get(start, 0), length)
        return any(virus in taxa for (start, length), taxa in spans.items()
                   if length == starts[start])

    def _held_out_warnings(self, edges) -> list:
        """The HELD-OUT LABEL caveat for these edges, or nothing.

        A LIST so a caller can `warnings=[...] + self._held_out_warnings(h)`
        without a branch. The caveat's label is enforced by verify(), so
        this is the part a model cannot drop for brevity.
        """
        counts: dict[str, int] = defaultdict(int)
        for e in edges:
            if e["predicate"] in self._held_out:
                counts[e["predicate"]] += 1
        if not counts:
            return []
        return [guards.caveat("held_out_label",
                              guards.held_out_caveat(dict(counts)))]

    # --------------------------------------------------------------- tools
    def resolve(self, text: str, limit: int = 5) -> ToolResult:
        """Free text to node identifiers.

        Ambiguity is returned, not resolved. Entity resolution is where a
        pipeline silently answers a question about the wrong molecule.
        """
        t = (text or "").strip().lower()
        if not t:
            return ToolResult(False, "resolve", note="empty query")
        exact = self.names.get(t)
        hits = [exact] if exact else []
        if not exact:
            hits = [nid for name, nid in self.names.items()
                    if t in name][:limit]
        if not hits:
            return ToolResult(False, "resolve",
                              note=f"'{text}' does not appear in this graph. It "
                                   "may be absent rather than misspelled: the "
                                   "graph covers compounds with published "
                                   "coronavirus or human-target activity.")
        return ToolResult(
            True, "resolve",
            data=[{"id": h, "label": self._label(h),
                   "class": self.nodes[h]["class"]} for h in hits],
            note="ambiguous — ask which was meant" if len(hits) > 1 else "")

    def explain(self, compound: str, virus: str,
                max_paths: int = MAX_PATHS) -> ToolResult:
        """Mechanistic paths from a compound to a virus, with provenance.

        The supported use of this graph. Given a compound with observed
        activity, this returns why it might act, as chains of individually
        sourced facts.
        """
        # ROUTES ARE PASSED BY THE TOOL, not by whoever called it.
        # check() refuses a routed query only when `routes is not None`,
        # and three of the four tools omitted it -- so the no-route
        # refusal was a property of the orchestrator, which happens to
        # pass routes itself, rather than of the tool. Any second entry
        # point lost it silently, and tools.pathways on Yellow fever
        # (0 routes, the tier the refusal was written for) returned
        # ok=True. Found by external audit, 2026-11-02.
        verdict = guards.check(guards.QueryKind.EXPLAIN, virus,
                               calibration=self.calibration_for(virus),
                               routes=self.routes_for(virus))
        # THE VERDICT IS ACTED ON BEFORE ANY RETRIEVAL. This computed it
        # and then read only verdict.WARNINGS, so a zero-route virus
        # fell through to the no-paths branch below and was told "that
        # is an absence of evidence, not evidence of absence: most
        # compounds have an activity measurement and no target
        # annotation" -- a sentence about the COMPOUND, for a virus
        # nothing can reach. ok was False either way, which is why
        # a test asserting only `ok is False` could not tell them apart.
        if not verdict.allowed:
            return ToolResult(False, "explain",
                              warnings=list(verdict.warnings),
                              note=verdict.reason)
        paths = self._find_paths(compound, virus, max_paths)
        if not paths:
            return ToolResult(
                False, "explain", warnings=verdict.warnings,
                note=f"No mechanistic path connects {self._label(compound)} to "
                     f"{self._label(virus)} in this graph. That is an absence of "
                     "evidence, not evidence of absence: most compounds have an "
                     "activity measurement and no target annotation.")

        lines, warns, ids = [], list(verdict.warnings), []
        for i, (route, edges) in enumerate(paths, 1):
            lines.append(f"Route {i} ({route}):")
            for e in edges:
                lines.append("  " + self._verbalise(e))
                ids.append(e["_id"])
            nodes = [self.nodes.get(e["subject"]) or {} for e in edges]
            nodes += [self.nodes.get(e["object"]) or {} for e in edges]
            # Attach each warning to the route that triggered it. Collecting
            # them at the compound level produced a cell-context caveat naming
            # SIGMAR1 printed against a route through Spike glycoprotein --
            # a mismatch that discredits every other warning on the page.
            #
            # THE SAME WARNINGS GO INTO `warns`, not only into the printed
            # route. They used to go only here, so cell_context and
            # host_directed_route never reached ToolResult.warnings, never
            # reached the brief's required_warnings, and verify()'s
            # warning-survival check never required them: the model could drop
            # the SIGMAR1 caveat entirely and the answer still passed. That is
            # the one caveat this project exists to attach -- chloroquine's
            # route is real, well-replicated, and led the field astray in 2020.
            # Both halves were tested and the join between them was not:
            # test_agent_tools asserts the text lands in the right route block,
            # test_agent_verify hands verify() a hand-built caveat directly.
            route_warnings = guards.path_warnings(
                nodes, route, calibration=self.calibration_for(virus),
                landscape=self.landscape, virus=virus)
            # A ROUTE WHOSE OWN EDGE IS CONTRADICTED SAYS SO. Attached per
            # route, like the cell-context caveat above and for the same
            # reason: a contradiction naming one protein printed against a
            # route through another discredits every warning on the page.
            route_warnings = list(route_warnings) + self._contradictions(edges)
            for w in route_warnings:
                lines.append(f"  ! {w}")
            warns.extend(route_warnings)

        warns.append(self._selectivity_for(compound, virus))
        return ToolResult(True, "explain",
                          data=[{"route": r, "n_edges": len(es)} for r, es in paths],
                          verbalised=lines,
                          warnings=list(dict.fromkeys(w for w in warns if w)),
                          edge_ids=ids)

    def _contradictions(self, edges: list[dict]) -> list[str]:
        """Caveats for any INHIBITS edge on this route with an inactive twin.

        The agent loads EVERY edge, held-out predicates included, so the
        contradicting measurement is already in hand -- this reads it rather
        than fetching anything. Scoring correctly never walks it; explaining
        correctly cannot ignore it.
        """
        out: list[str] = []
        for e in edges:
            if e.get("predicate") != "INHIBITS":
                continue
            twins = [t for t in self.out.get(
                ("MEASURED_INACTIVE_AGAINST", e["subject"]), [])
                if t.get("object") == e["object"]]
            if not twins:
                continue
            q = twins[0].get("qualifiers") or {}
            bits = [f"{k}={q[k]:g}" for k in ("ic50_nm", "ec50_nm", "ki_nm")
                    if q.get(k) is not None]
            rel = q.get("relation")
            detail = ((rel or "") + " " + ", ".join(bits)).strip()
            out.append(guards.caveat(
                "contradicted",
                guards.contradicted_caveat(self._label(e["subject"]),
                                           self._label(e["object"]), detail)))
        return out

    def _best_direct_potency(self, compound: str) -> float | None:
        """Most potent INHIBITS measurement, in nM."""
        vals = [q[k]
                for e in self.out.get(("INHIBITS", compound), [])
                for q in [e.get("qualifiers") or {}]
                for k in ("ic50_nm", "ki_nm", "kd_nm")
                if q.get(k) is not None]
        return min(vals) if vals else None

    def _selectivity_for(self, compound: str, virus: str) -> str:
        """Selectivity as a RANGE, not a single value.

        The ingest keeps the worst selectivity across studies, which is right
        for filtering a positive set and wrong for display: reporting
        nirmatrelvir as "cytotoxic, SI 1.8" from one unflattering paper, with
        no sign that others disagree, states something false about a licensed
        drug. Studies disagree because cell lines and protocols differ, and
        the disagreement is information.
        """
        vals = list(self.si_range.get((compound, virus), []))
        if not vals:
            vals = [q["selectivity_index"]
                    for e in self.out.get(("HAS_ANTIVIRAL_ACTIVITY_AGAINST", compound), [])
                    if e["object"] == virus
                    and (q := e.get("qualifiers") or {}).get("selectivity_index") is not None]
        if not vals:
            return guards.selectivity_note(None, None)
        note = guards.selectivity_note(min(vals), True,
                                       n_studies=len(vals), si_max=max(vals))
        cl = guards.cell_line_note(self.si_cell.get((compound, virus)),
                                   family=self.family_of(virus))
        return (note + " " + cl).strip() if cl else note

    def _find_paths(self, compound: str, virus: str, limit: int
                    ) -> list[tuple[str, list[dict]]]:
        """Enumerate metapath instances. Bounded and non-revisiting."""
        found: list[tuple[str, list[dict]]] = []

        # The constraints come from the SCHEMA, via baselines.qualifier_ok --
        # the same function dwpc_scores applies -- so a constraint declared
        # once is honoured by both the scorer and the agent. M2's dependency
        # check used to be written out here by hand; it is the schema's
        # HOST_FACTOR_FOR.direction=dependency and is now read from there.
        c1, c2 = self._constraints.get("M1"), self._constraints.get("M2")

        # M1: compound -> viral protein -> gene -> virus
        for e1 in self.out.get(("INHIBITS", compound), []):
            if not qualifier_ok(e1, "INHIBITS", c1):
                continue
            for e2 in self.out.get(("ENCODED_BY", e1["object"]), []):
                if not qualifier_ok(e2, "ENCODED_BY", c1):
                    continue
                for e3 in self.out.get(("BELONGS_TO", e2["object"]), []):
                    if e3["object"] != virus:
                        continue
                    if not qualifier_ok(e3, "BELONGS_TO", c1):
                        continue
                    found.append(("M1 direct-acting", [e1, e2, e3]))
                    if len(found) >= limit:
                        return found

        # M2: compound -> host protein -> virus (dependency only)
        for e1 in self.out.get(("TARGETS", compound), []):
            if not qualifier_ok(e1, "TARGETS", c2):
                continue
            for e2 in self.out.get(("HOST_FACTOR_FOR", e1["object"]), []):
                if e2["object"] != virus:
                    continue
                if not qualifier_ok(e2, "HOST_FACTOR_FOR", c2):
                    continue
                found.append(("M2 host-directed", [e1, e2]))
                if len(found) >= limit:
                    return found
        return found

    def triage(self, compounds: list[str], virus: str) -> ToolResult:
        """Rank compounds you already have.

        Supported because the recall ceiling does not apply: these compounds
        are in the graph by construction.
        """
        verdict = guards.check(guards.QueryKind.TRIAGE, virus,
                               calibration=self.calibration_for(virus),
                               routes=self.routes_for(virus))
        if not verdict.allowed:
            return ToolResult(False, "triage",
                              warnings=list(verdict.warnings),
                              note=verdict.reason)
        # Report anything that could not be scored. Silently returning three
        # of four requested compounds is a failure a user cannot detect.
        missing = [c for c in compounds if c not in self.nodes]
        scored = []
        for c in compounds:
            if c in missing:
                continue
            paths = self._find_paths(c, virus, MAX_PATHS)
            direct = sum(1 for r, _ in paths if r.startswith("M1"))
            host = len(paths) - direct
            sel = self._selectivity_for(c, virus)
            scored.append({"compound": c, "label": self._label(c),
                           "direct_acting_routes": direct,
                           "host_directed_routes": host,
                           "selectivity": sel})
        # Counting routes is wrong and inverted the controls: chloroquine
        # ranked ABOVE nirmatrelvir because it has a 160 nM Spike hit and a
        # 7,280 nM nsp5 hit -- two routes, one of them a screening artifact --
        # and then won the tiebreak on host-directed routes, which score BELOW
        # chance and should count for nothing.
        #
        # Rank on what the evaluation actually established: potent
        # direct-acting evidence, then selectivity. Host-directed routes are
        # reported for context and contribute zero.
        for r in scored:
            r["best_potency_nm"] = self._best_direct_potency(r["compound"])
            r["selectivity_index"] = self.si_range.get(
                (r["compound"], virus), [None])[0] \
                if self.si_range.get((r["compound"], virus)) else None
        def key(r):
            si = r["selectivity_index"]
            potent = (r["best_potency_nm"] is not None
                      and r["best_potency_nm"] <= 1000)
            # SELECTIVITY BEFORE POTENCY. Ordering on potency first ranked
            # chloroquine (160 nM, SI 15) above remdesivir (1,560 nM, SI 107)
            # -- and potency without selectivity is precisely the signal that
            # made chloroquine look promising in 2020. That reasoning stands
            # on its own; the evidence once cited for it does not.
            # "Selectivity took M1's AUC from 0.729 to 0.826" was measured by
            # filtering positives to compounds with verified selectivity data
            # while leaving negatives whole, which also moved cov_pos 49% ->
            # 69% and left cov_neg at 4%. Most of that gain was annotation
            # density: symmetric, M1 reads 0.660, and 0.499 once same-paper
            # evidence is withheld. Selectivity still comes first because
            # potency alone is the chloroquine failure, not because of 0.826.
            #
            # Compounds with an unverified selectivity index rank below those
            # with one: unknown is not a passing grade.
            # NO POTENCY GATE. A 1 uM threshold has no basis in this
            # project's evaluation and put remdesivir (1,560 nM, SI 107) below
            # chloroquine (160 nM, SI 15) purely by sitting the wrong side of
            # an arbitrary line. The ranking may only use quantities the
            # evaluation measured: direct-acting evidence, then selectivity,
            # with potency as a tiebreak.
            del potent
            return (0 if r["direct_acting_routes"] else 1,
                    0 if si is not None else 1,
                    -(si or 0),
                    r["best_potency_nm"] if r["best_potency_nm"] is not None else 1e12)
        scored.sort(key=key)
        lines = []
        for i, r in enumerate(scored, 1):
            pot = (f"best direct IC50 {r['best_potency_nm']:.0f} nM"
                   if r["best_potency_nm"] is not None else "no direct-acting data")
            si = (f", SI {r['selectivity_index']:.0f}"
                  if r["selectivity_index"] is not None else ", SI unknown")
            lines.append(f"{i}. {r['label']}: {pot}{si} "
                         f"({r['direct_acting_routes']} direct-acting, "
                         f"{r['host_directed_routes']} host-directed routes)")
        warns = list(verdict.warnings)
        warns.append(guards.caveat(
            "ranking",
            # Describes the sort key above, in its actual order. The previous
            # wording said "by direct-acting route count", which is what the
            # ranking deliberately STOPPED doing -- counting routes is what
            # put chloroquine above nirmatrelvir on a screening artifact.
            "compounds are ordered by whether they have any direct-acting "
            "route, then by whether a selectivity index exists, then by higher "
            "selectivity, with the most potent direct-acting measurement as a "
            "final tiebreak. The NUMBER of routes does not affect the order, "
            "host-directed routes contribute nothing to it, and a compound "
            "whose selectivity is unknown ranks below one with a measured "
            "index. "
            + guards.direct_acting_caveat(self.calibration_for(virus))
            + ". Host-directed routes contribute nothing to this order and "
            "are shown for context only."))
        # UNCONDITIONAL, AND THE RANKING CAVEAT ABOVE IS NOT A SUBSTITUTE.
        # That one describes the sort key: which quantities move a compound up
        # the list. This one states what the landscape measured about the
        # route itself. A reader can accept "host-directed routes contribute
        # nothing to this order" and still read the host-directed count beside
        # each compound as weak supporting evidence, which is the claim the
        # degree-matched null refuses. Appended even with no landscape, where
        # it reads "not measured on this release": an unstated null and a
        # measured one must not look the same from outside.
        warns.append(guards.caveat(
            "pathway_route", guards.pathway_caveat(self.landscape, virus)))
        if missing:
            warns.append(f"{len(missing)} requested compound(s) are not in this "
                         f"graph and were not ranked: {', '.join(missing)}")
        return ToolResult(True, "triage", data=scored, verbalised=lines,
                          warnings=warns, note=f"{len(missing)} unranked" if missing else "")

    def _is_viral(self, node_id: str) -> bool:
        """Whether a node is a viral protein, from its own properties.

        A missing node is not viral: an absent id cannot be asserted to be
        anything, and treating it as viral would silently drop it from a
        human-protein census.
        """
        n = self.nodes.get(node_id)
        if not n:
            return False
        return bool((n.get("properties") or {}).get("is_viral"))

    def _pathway_members(self, pathway: str) -> set[str]:
        return {e["subject"]
                for e in self.inv.get(("PARTICIPATES_IN", pathway), [])}

    def _pathway_background(self, virus: str) -> dict:
        """The rate at which pathways contain a dependency factor, BY SIZE.

        WHY A BARE RATE WAS NOT ENOUGH, which the second live run showed as
        clearly as the first showed the need for any background at all. The
        unmatched figure came out at 58% against chloroquine's 90%, and the
        answer called that "+33% on the background". It is an artefact of
        size. The schema's own measurements put the median pathway at 4
        members and the 99th percentile at 113; a 4-member pathway contains
        one of SARS-CoV-2's 3,460 dependency factors about half the time and
        a 95-member pathway essentially always. Chloroquine's targets
        include mTOR, which sits in large pathways, so its set averages far
        above the median and a higher hit rate follows from the sizes alone.

        This is the SAME confound the landscape corrects at the protein
        level by degree-matching -- well-connected things are likelier to be
        hit -- reappearing one level up in the object being counted. So the
        background is binned by size, by quantile for the same reason
        degree_bins uses quantiles: pathway size is heavy-tailed and
        equal-width bins would put nearly everything in the first one.

        Returns the per-bin rates, the bin edges, and the unmatched totals,
        so an answer can report the matched expectation and still say what
        the raw background was.
        """
        key = ("pathway_background", virus)
        if key in self._cache:
            return self._cache[key]
        deps = self._dependency_factors(virus)
        cap = self._max_pathway_members
        sized: list[tuple[int, bool]] = []
        for pw, n in self.nodes.items():
            if n.get("class") != "Pathway":
                continue
            members = self._pathway_members(pw)
            if cap is not None and len(members) > cap:
                continue
            sized.append((len(members), bool(members & deps)))
        # PER EXACT SIZE, NOT QUANTILE BINS. Ties on size are the rule and
        # not the exception -- the median pathway has 4 members and
        # thousands share it -- and quantile binning over tied data is
        # unsound twice over.
        #
        # THE BUG IT REPLACES, found by external audit on 2026-11-02 and
        # reproduced here before this was written. `sized` held
        # (size, hit) tuples and the code sorted it with a bare
        # `sized.sort()`, which orders by size AND THEN by the hit flag --
        # and False < True. So inside every tied-size block each non-hit
        # sorted before each hit; the quantile cut then split that block
        # and handed the lower piece a rate of 0. `_size_matched_probs`
        # returns the FIRST bin whose upper edge reaches the size, which
        # for a tied size is exactly that zero-rate piece.
        #
        # Measured on a simulated distribution of 4,000 pathways, sizes
        # 1-200, with a UNIFORM true hit rate of 0.12:
        #
        #     size 1 (982 pathways, the commonest)  true 0.115 -> 0.000
        #     size 2 (599)                          true 0.109 -> 0.282
        #     size 5 (279)                          true 0.122 -> 0.085
        #
        # and where every pathway shares one size, the returned rate is
        # 0.000 whatever the truth. Too low an expectation means too high
        # a ratio and too small a p-value, so the error ran in the one
        # direction that manufactures findings -- in the headline number
        # of the feature. Sorting on the size alone fixes the ordering and
        # not the method: quantile cuts still slice one size across
        # several bins.
        #
        # So the rate is estimated for each size from the pathways of
        # THAT size, pooling outward over neighbouring sizes only while
        # the sample is too small to estimate a rate from.
        from collections import defaultdict

        per: dict[int, list[int]] = defaultdict(lambda: [0, 0])
        for size, hit in sized:
            per[size][1] += 1
            per[size][0] += int(hit)
        total = len(sized)
        out = {"by_size": {k: (v[0], v[1]) for k, v in sorted(per.items())},
               "total": total,
               "hits": sum(1 for _, h in sized if h),
               "n_dependency": len(deps)}
        self._cache[key] = out
        return out

    def _size_matched_probs(self, virus: str, sizes: list[int]) -> list[float]:
        """The background hit probability for a pathway of each given size.

        Estimated from the pathways of that exact size, widening to the
        nearest neighbouring sizes only until the pooled sample reaches
        MIN_BACKGROUND_SAMPLES. A size the graph holds no pathway of
        takes the nearest sizes that bracket it, by the same widening.
        """
        bg = self._pathway_background(virus)
        if not bg["total"]:
            return [0.0] * len(sizes)
        by_size = bg["by_size"]
        known = sorted(by_size)
        # THE FLOOR CANNOT EXCEED THE GRAPH. On the real release 100 is
        # 0.8% of 13,236 pathways, so pooling stays local. On a small
        # graph an unconditional 100 swallows the whole distribution and
        # every size gets the same rate -- the flat background this
        # method exists to replace, arrived at from the other side. A
        # window wider than a tenth of all pathways is not size-matched
        # in any useful sense, so that is the cap, and a noisier local
        # estimate is preferred to a precise global one.
        floor = min(MIN_BACKGROUND_SAMPLES, max(1, bg["total"] // 10))
        probs = []
        for s in sizes:
            # Widen outward from s by DISTANCE IN SIZE, nearest first, so
            # the estimate stays as local as the data allows. A pathway of
            # 4 members is informed by pathways of 3 and 5 long before one
            # of 200.
            hits = seen = 0
            for size in sorted(known, key=lambda k: (abs(k - s), k)):
                h, n = by_size[size]
                hits += h
                seen += n
                if seen >= floor:
                    break
            probs.append(hits / seen if seen else 0.0)
        return probs

    @staticmethod
    def _background_median(bg: dict) -> int:
        """The median member count over the background population.

        Read off `by_size`, which is a histogram, so the median is the
        size at which the running count passes half the population. The
        previous version took the lower edge of the middle QUANTILE BIN,
        which is a different quantity whenever the bins are uneven --
        and over tied sizes they always are.
        """
        total = bg.get("total") or 0
        if not total:
            return 0
        seen = 0
        for size, (_hits, n) in sorted(bg["by_size"].items()):
            seen += n
            if seen * 2 >= total:
                return size
        return 0

    def _size_matched_expectation(self, virus: str,
                                  sizes: list[int]) -> float:
        """How many of these pathways a size-matched background would hit."""
        return sum(self._size_matched_probs(virus, sizes))

    @staticmethod
    def _upper_tail(probs: list[float], k: int) -> float:
        """P(X >= k) for independent Bernoulli trials with these
        probabilities -- the Poisson-binomial upper tail, computed exactly.

        WHY A TEST AND NOT A THRESHOLD. This sentence used to turn on
        `ratio <= 1.05`, a number with nothing behind it deciding which of
        two things a reader was told. On the live run the ratio came out at
        1.04 and the verdict therefore hung on the second decimal place of
        a constant picked by hand -- exactly the fitted-parameter-with-no-
        record this project refuses everywhere else.

        Exact rather than sampled: a permutation null drawing pathways from
        matched size bins IS a sum of independent Bernoullis, so permuting
        would only add sampling noise to a quantity that can be computed.
        n is the number of pathways, so the O(n^2) recurrence is nothing.

        THE ASSUMPTION IT MAKES, which the caller must state: pathways
        overlap and share members, so "contains a dependency factor" is
        positively correlated between them and the trials are not
        independent. That makes this tail anti-conservative -- a small p
        here is a flag to examine, never a finding.
        """
        dist = [1.0]
        for p in probs:
            nxt = [0.0] * (len(dist) + 1)
            for i, v in enumerate(dist):
                nxt[i] += v * (1.0 - p)
                nxt[i + 1] += v * p
            dist = nxt
        return sum(dist[min(k, len(dist) - 1):]) if k <= len(dist) - 1 else 0.0

    def _dependency_factors(self, virus: str) -> set[str]:
        """Proteins this virus DEPENDS on, by CRISPR screen.

        The direction qualifier is required rather than assumed: a
        restriction factor is the opposite claim, and pooling the two makes
        every count below mean nothing. Same definition as
        pathway_landscape.dependency_factors, because the agent must not
        count a different set from the one the measured null was computed
        over.
        """
        out = set()
        for e in self.inv.get(("HOST_FACTOR_FOR", virus), []):
            if (e.get("qualifiers") or {}).get("direction") == "dependency":
                out.add(e["subject"])
        return out

    def pathways(self, compound: str, virus: str,
                 limit: int = MAX_PATHWAYS) -> ToolResult:
        """Which pathways this compound's targets sit in, and whether the
        virus depends on any member.

        WHY THIS IS NOT explain(). _find_paths walks M1 and M2 only, so the
        pathway metapath the schema declares as M4 is never traversed by the
        agent at all -- a question about pathway membership had no answer
        before this, only a silence that read like a negative. And a path
        sampler would be the wrong answer even once M4 were walked: it
        returns MAX_PATHS instances chosen by iteration order, which tells a
        reader nothing about how many pathways there are, how many were
        skipped, or which of the compound's targets contributed none. The
        question here is a census, so the answer is counted, and every
        absence in it is reported as an absence.

        WHAT IT REFUSES TO IMPLY. Membership is not benefit. Most TARGETS
        edges record potency with no direction, so the compound may activate
        a protein the virus depends on as easily as inhibit it, and the
        caveat says so from the edges actually used. The measured
        degree-matched null rides along too: on this release the drugged
        proteome is no likelier to hit a dependency factor than chance, so a
        pathway hit is context, not evidence.
        """
        verdict = guards.check(guards.QueryKind.EXPLAIN, virus,
                               calibration=self.calibration_for(virus),
                               routes=self.routes_for(virus))
        base_warnings = list(verdict.warnings) + [guards.caveat(
            "pathway_route", guards.pathway_caveat(self.landscape, virus))]
        # AND THE VERDICT HAS TO BE ACTED ON. Passing `routes` is
        # necessary and not sufficient: this took verdict.WARNINGS and
        # never read verdict.ALLOWED, so a refusal arrived as a caveat
        # attached to the answer it was supposed to prevent. Only
        # virus_pathways checked it. Found while fixing A6 of the
        # 2026-11-02 audit, which had identified the missing argument
        # and not this.
        if not verdict.allowed:
            return ToolResult(False, "pathways", warnings=base_warnings,
                              note=verdict.reason)

        if compound not in self.nodes:
            return ToolResult(False, "pathways", warnings=base_warnings,
                              note=f"{compound} is not in this graph, so it "
                                   f"has no targets to place in a pathway.")

        # Human targets only. A viral protein can be hit by a compound and is
        # never a host pathway member, so including one would put a viral
        # protein in a census of the human proteome.
        target_edges = [e for e in self.out.get(("TARGETS", compound), [])
                        if not self._is_viral(e["object"])]
        if not target_edges:
            # WHY IT HAS NO HUMAN TARGET MATTERS, and the graph knows which
            # of two very different reasons applies. pathway_landscape's
            # census already partitions the population this way: 3,002
            # molecules carry a human TARGETS edge, 3,689 carry INHIBITS
            # only, and 18,692 carry neither.
            #
            # Remdesivir is the second kind. Answering it with "most small
            # molecules here carry an activity measurement and no target
            # annotation" is true of the population and false about the
            # compound: its target is annotated, well-evidenced, and viral.
            # A direct-acting antiviral has no host pathway because of what
            # it is, not because the source is thin, and conflating the two
            # tells a reader the graph is missing something when it is not.
            inhibits = self.out.get(("INHIBITS", compound), [])
            if inhibits:
                viral = sorted({e["object"] for e in inhibits
                                if self._is_viral(e["object"])})
                # The schema puts object_constraints {is_viral: true} on
                # INHIBITS, so a non-viral object is a release that does not
                # satisfy its own schema. Dropping those silently would make
                # a broken release look like a thin one -- the same
                # conflation this branch exists to end -- so they are
                # counted and named instead.
                rogue = sorted({e["object"] for e in inhibits
                                if not self._is_viral(e["object"])})
                if viral:
                    bad = ("" if not rogue else
                           f" {len(rogue)} further INHIBITS edge(s) point at "
                           f"a NON-viral protein "
                           f"({', '.join(self._label(r) for r in rogue)}), "
                           f"which the schema's object_constraints forbid: "
                           f"this release does not satisfy its own schema, "
                           f"and those edges are reported rather than "
                           f"dropped.")
                    return ToolResult(
                        False, "pathways", warnings=base_warnings,
                        edge_ids=[e["_id"] for e in inhibits],
                        note=(f"{self._label(compound)} has no HUMAN target "
                              f"in this graph, so there is no host pathway "
                              f"to place it in — but it does inhibit "
                              f"{len(viral)} viral protein(s): "
                              f"{', '.join(self._label(v) for v in viral)}. "
                              f"It is direct-acting here, so a host-pathway "
                              f"census is not the question for it; ask for "
                              f"its mechanism instead. This is a property "
                              f"of the compound, not a gap in the "
                              f"source.{bad}"))
                if rogue:
                    return ToolResult(
                        False, "pathways", warnings=base_warnings,
                        edge_ids=[e["_id"] for e in inhibits],
                        note=(f"{self._label(compound)} has no TARGETS edge, "
                              f"and its {len(rogue)} INHIBITS edge(s) point "
                              f"at non-viral protein(s) "
                              f"({', '.join(self._label(r) for r in rogue)}) "
                              f"— which the schema's object_constraints on "
                              f"INHIBITS forbid. This release does not "
                              f"satisfy its own schema, so no pathway "
                              f"census is attempted; run "
                              f"scripts/validate_release.py."))
            return ToolResult(
                False, "pathways", warnings=base_warnings,
                note=f"{self._label(compound)} has no human target in this "
                     f"graph, and no viral one either, so there is no "
                     f"pathway to place it in. That is a property of the "
                     f"source and not a finding about the compound: most "
                     f"small molecules here carry an activity measurement "
                     f"and no target annotation at all.")

        # A CONTESTED DIRECTION IS NOT A STATED ONE. This counts edges
        # that say what the drug DOES to the target, and two sources
        # asserting opposite actions do not. Counting "contested" as
        # directed would report a mechanism where there is a
        # contradiction -- and "unknown" was already excluded only by
        # being falsy, which is luck rather than a rule.
        n_directed = sum(
            1 for e in target_edges
            if (e.get("qualifiers") or {}).get("direction")
            not in (None, "", "unknown", "contested"))
        deps = self._dependency_factors(virus)
        cap = self._max_pathway_members

        # pathway id -> {members, this compound's targets in it}
        hits: dict[str, dict] = {}
        oversized: set[str] = set()
        no_pathway: list[str] = []
        edge_ids = [e["_id"] for e in target_edges]

        for te in target_edges:
            prot = te["object"]
            memberships = self.out.get(("PARTICIPATES_IN", prot), [])
            placed = False
            for pe in memberships:
                pw = pe["object"]
                if pw not in hits and pw not in oversized:
                    members = {m["subject"]
                               for m in self.inv.get(("PARTICIPATES_IN", pw), [])}
                    # THE SAME CAP THE LANDSCAPE APPLIED. Quoting a measured
                    # null computed over capped pathways while reporting an
                    # uncapped one would be two different questions under one
                    # answer.
                    if cap is not None and len(members) > cap:
                        oversized.add(pw)
                        continue
                    hits[pw] = {"members": members, "targets": set(),
                                "edges": []}
                if pw in oversized:
                    continue
                hits[pw]["targets"].add(prot)
                hits[pw]["edges"].append(pe["_id"])
                placed = True
            if not placed:
                no_pathway.append(prot)

        rows = []
        for pw, h in hits.items():
            dep_members = sorted(h["members"] & deps)
            rows.append({
                "pathway": pw, "label": self._label(pw),
                "n_members": len(h["members"]),
                "targets_in_pathway": sorted(h["targets"]),
                "dependency_members": dep_members,
                "n_dependency_members": len(dep_members),
            })
        # Pathways the virus depends on something in, first; then by how much
        # of the pathway is a dependency factor. NOT a ranking of relevance:
        # ordering decides what survives the limit, and a reader is better
        # served by seeing the hits than by seeing the first N by node id.
        rows.sort(key=lambda r: (-r["n_dependency_members"], -len(
            r["targets_in_pathway"]), r["pathway"]))
        with_dep = [r for r in rows if r["n_dependency_members"]]
        shown = rows[:limit]

        bg = self._pathway_background(virus)
        sizes = [r["n_members"] for r in rows]
        probs = self._size_matched_probs(virus, sizes)
        expected = sum(probs)
        p_excess = self._upper_tail(probs, len(with_dep))
        here = len(with_dep) / len(rows) if rows else 0.0
        raw = bg["hits"] / bg["total"] if bg["total"] else 0.0
        med = sorted(sizes)[len(sizes) // 2] if sizes else 0
        # NO BACKGROUND IS NOT AN INFINITE RATIO. This read
        # `float("inf")` and rendered literally as `ratio inf` beside a
        # p-value of 0.0000 -- the strongest possible claim produced by
        # the weakest possible null. An expectation of zero means the
        # graph holds no comparable pathway that contains a dependency
        # factor, which is a statement about coverage and not about this
        # compound.
        ratio = (len(with_dep) / expected) if expected else None
        ratio_text = f"ratio {ratio:.2f}" if ratio is not None else (
            "no ratio: the size-matched background expects none, so there "
            "is nothing to compare against")
        lines = [
            (f"{self._label(compound)} has {len(target_edges)} human "
             f"target(s) in this graph."),
            (f"{len(rows)} pathway(s) at or under the {cap}-member cap "
             f"contain at least one of them; {len(with_dep)} of those "
             f"({here:.0%}) contain at least one protein "
             f"{self._label(virus)} depends on."),
            # NEVER THE COUNT WITHOUT A SIZE-MATCHED BACKGROUND. The raw
            # rate is reported too, and the gap between the two is the
            # measurement: a larger-than-typical pathway set hits more
            # often for no reason but its size.
            (f"A SIZE-MATCHED background expects {expected:.1f} of them to "
             f"contain one ({ratio_text}). The virus has "
             f"{bg['n_dependency']:,} dependency factors, so a pathway's "
             f"size largely decides whether it holds one: unmatched, "
             f"{bg['hits']:,} of {bg['total']:,} ({raw:.0%}) of all "
             f"pathways under the cap qualify, but that population has a "
             f"median of {self._background_median(bg)} members "
             f"against {med} here. "
             + (f"Against that expectation the excess is not significant "
                f"(p {p_excess:.2f}), so membership is no more than these "
                f"pathways' sizes already predict and is not evidence of "
                f"anything about this compound."
                if p_excess >= LANDSCAPE_ALPHA else
                f"The excess over the size-matched expectation reaches p "
                f"{p_excess:.4f}. Pathways overlap and share members, so "
                f"this tail assumes an independence that does not hold and "
                f"is anti-conservative: treat it as a flag to examine, not "
                f"a finding.")),
        ]
        if oversized:
            lines.append(
                f"{len(oversized)} further pathway(s) were EXCLUDED for "
                f"exceeding the {cap}-member cap. They are not absent from "
                f"the graph, only from this answer and from the measured "
                f"null it quotes.")
        if no_pathway:
            lines.append(
                f"{len(no_pathway)} target(s) sit in no pathway at or under "
                f"the cap: {', '.join(self._label(p) for p in no_pathway)}.")
        if not with_dep:
            lines.append(
                f"No pathway containing a target of this compound contains a "
                f"protein {self._label(virus)} is known to depend on. That is "
                f"an absence in the CRISPR layer as much as a statement about "
                f"this compound.")
        for r in shown:
            # BOUNDED, WITH THE REMAINDER COUNTED. The first live run put
            # forty-one protein names on a single line and did it for
            # twenty-five pathways. Nothing in it was wrong and none of it
            # was readable, and a brief built from it spends its budget on
            # name lists rather than on the facts the answer rests on. The
            # full set stays in `data` for a caller that wants it.
            names = [self._label(d) for d in
                     r["dependency_members"][:MAX_DEPENDENCY_NAMES]]
            extra = r["n_dependency_members"] - len(names)
            dep = (", ".join(names) + (f", +{extra} more" if extra else "")
                   if names else "none")
            lines.append(
                f"  {r['label']} ({r['n_members']} members) — via "
                f"{', '.join(self._label(t) for t in r['targets_in_pathway'])}"
                f"; dependency factors in it: {dep}")
        if len(rows) > limit:
            lines.append(f"  … {len(rows) - limit} further pathway(s) not "
                         f"shown, of which "
                         f"{sum(1 for r in rows[limit:] if r['n_dependency_members'])} "
                         f"contain a dependency factor.")

        for h in hits.values():
            edge_ids.extend(h["edges"])
        warns = base_warnings + [guards.caveat(
            "target_direction",
            guards.target_direction_caveat(n_directed, len(target_edges)))]
        return ToolResult(True, "pathways", data=rows, verbalised=lines,
                          warnings=warns, edge_ids=edge_ids,
                          note=f"{len(oversized)} pathway(s) over the cap"
                               if oversized else "")

    def virus_pathways(self, virus: str,
                       limit: int = MAX_PATHWAYS) -> ToolResult:
        """Which pathways this virus depends on, and which of them are
        druggable.

        THE QUESTION THE AGENT COULD NOT ANSWER. pathways() takes a
        COMPOUND and asks whether the virus depends on anything in the
        pathways that compound's targets sit in. "Which pathways does
        this virus depend on" is the other direction and is probably
        the first thing a reader asks of a graph built for pathway
        analysis. The orchestrator's PATHWAY branch required a
        compound, so the question fell through to a node profile.

        WHY THERE IS NO P-VALUE HERE, and why that is not timidity.
        pathways() reports an excess against a size-matched background
        because its pathway set is chosen by the COMPOUND and the
        outcome asked of it -- does it hold a dependency factor -- is a
        different event. Here the set is chosen BY holding a dependency
        factor, so selection and outcome are the same event and any
        enrichment tail would be circular by construction: every row
        would "contain a dependency factor" with probability one.

        What is not circular is DENSITY. A pathway of 90 members holds
        more of anything than a pathway of 4, so the count alone ranks
        by size. Each row therefore carries the count expected from its
        size and the virus's own dependency-factor density, and the
        ratio between them -- and the ratio is offered as a lead, not a
        finding, because the virus-level version of exactly this
        question is measured in pathway_landscape.py and nothing in it
        survived correction on this release.
        """
        verdict = guards.check(guards.QueryKind.PATHWAY, virus,
                               calibration=self.calibration_for(virus),
                               routes=self.routes_for(virus))
        if not verdict.allowed:
            return ToolResult(False, "virus_pathways",
                              warnings=list(verdict.warnings),
                              note=verdict.reason)
        base_warnings = list(verdict.warnings) + [guards.caveat(
            "pathway_route", guards.pathway_caveat(self.landscape, virus))]

        deps = self._dependency_factors(virus)
        if not deps:
            # AN ABSENCE IN THE CRISPR LAYER, NOT A FINDING ABOUT THE
            # VIRUS. Five of the sixteen viruses in this register have
            # no dependency factor at all, and answering with an empty
            # list would read as "this virus depends on no pathway".
            restriction = {e["subject"] for e
                           in self.inv.get(("HOST_FACTOR_FOR", virus), [])}
            extra = (f" It has {len(restriction)} HOST_FACTOR_FOR edge(s), "
                     f"none carrying direction=dependency, so what exists "
                     f"is restriction or undirected rather than absent."
                     if restriction else
                     " It has no HOST_FACTOR_FOR edge at all.")
            return ToolResult(
                False, "virus_pathways", warnings=base_warnings,
                note=(f"No CRISPR screen in this release reports a protein "
                      f"{self._label(virus)} depends on, so there is no "
                      f"pathway census to build for it.{extra} That is an "
                      f"absence in the screen layer, not a statement that "
                      f"this virus depends on nothing."))

        cap = self._max_pathway_members
        oversized: set[str] = set()
        rows = []
        edge_ids: list[str] = []
        placed: set[str] = set()
        for d in deps:
            for pe in self.out.get(("PARTICIPATES_IN", d), []):
                pw = pe["object"]
                if pw in oversized:
                    continue
                members = self._pathway_members(pw)
                # THE SAME CAP THE LANDSCAPE APPLIED, for the reason it
                # applied it: quoting a measured null computed over
                # capped pathways beside an uncapped census would be two
                # questions under one answer.
                if cap is not None and len(members) > cap:
                    oversized.add(pw)
                    continue
                placed.add(d)
                row = next((r for r in rows if r["pathway"] == pw), None)
                if row is None:
                    in_pw = sorted(members & deps)
                    row = {"pathway": pw, "label": self._label(pw),
                           "n_members": len(members),
                           "dependency_members": in_pw,
                           "n_dependency_members": len(in_pw),
                           "druggable_members": sorted(
                               m for m in in_pw
                               if self.out.get(("TARGETS", m)) or
                               self.inv.get(("TARGETS", m)))}
                    row["n_druggable"] = len(row["druggable_members"])
                    # CURATED FROM THE SAME LITERATURE AS THE SCREEN.
                    # Reactome carries virus-specific pathways -- on the
                    # live SARS-CoV-2 run "SARS-CoV-2 modulates
                    # autophagy" came fourth by density. A CRISPR screen
                    # for a virus's dependency factors overlapping a
                    # pathway Reactome built FROM that virus's
                    # literature is close to definitional, and it is the
                    # shared-publication confound this project chases
                    # elsewhere, reappearing in the pathway vocabulary
                    # rather than in the labels. Flagged and not
                    # dropped: the biology is real, and a reader has to
                    # know which part of the overlap was decided in
                    # advance.
                    row["names_the_virus"] = self._names_virus(
                        row["label"], virus)
                    rows.append(row)
                edge_ids.append(pe["_id"])

        if not rows:
            return ToolResult(
                False, "virus_pathways", warnings=base_warnings,
                note=(f"{self._label(virus)} has {len(deps):,} dependency "
                      f"factor(s) and none of them carries a "
                      f"PARTICIPATES_IN edge to a pathway at or under the "
                      f"{cap}-member cap, so there is no census to report. "
                      f"That is a gap in pathway annotation, not in the "
                      f"screen."))

        # Density, not count: a 90-member pathway holds more of anything
        # than a 4-member one.
        n_prot = sum(1 for n in self.nodes.values()
                     if n.get("class") == "Protein"
                     and not (n.get("properties") or {}).get("is_viral"))
        rate = len(deps) / n_prot if n_prot else 0.0
        for r in rows:
            r["expected_members"] = r["n_members"] * rate
            r["density_ratio"] = (r["n_dependency_members"]
                                  / r["expected_members"]
                                  if r["expected_members"] else float("inf"))
        rows.sort(key=lambda r: (-r["n_dependency_members"],
                                 -r["density_ratio"], r["pathway"]))
        shown = rows[:limit]
        # TWO ORDERINGS, BECAUSE ONE COLUMN DID NOT SELECT THE OTHER.
        # The list above is ordered by COUNT, which is ordered by size:
        # a 200-member pathway holds more of everything. Printing a
        # ratio beside it invites a reader to compare rows the ratio
        # never chose, and the pathways that ratio WOULD choose -- a
        # 9-member pathway with 6 dependency factors -- are not on the
        # page at all.
        #
        # The floor is expectation >= 1, which is this project's
        # existing rule for when a lift is readable (lift_row,
        # calibration.evaluable), not a threshold invented here. Below
        # it the denominator is a fraction of a hit and the ratio is
        # arithmetic rather than evidence: a 2-member pathway holding
        # both factors reads 5.9 and means nothing.
        eligible = [r for r in rows if r["expected_members"] >= 1.0]
        by_ratio = sorted(eligible,
                          key=lambda r: (-r["density_ratio"],
                                         -r["n_dependency_members"],
                                         r["pathway"]))[:MAX_RATIO_ROWS]
        thin_ratio = len(rows) - len(eligible)
        unplaced = sorted(deps - placed)

        lines = [
            (f"{self._label(virus)} depends on {len(deps):,} host protein(s) "
             f"by CRISPR screen."),
            (f"{len(rows):,} pathway(s) at or under the {cap}-member cap "
             f"contain at least one of them; "
             f"{sum(1 for r in rows if r['n_druggable']):,} of those contain "
             f"at least one that a compound in this graph targets."),
            (f"Each row's EXPECTED count is its size times this virus's "
             f"dependency-factor density ({len(deps):,} of {n_prot:,} host "
             f"proteins = {rate:.1%}), so the ratio says whether a pathway "
             f"holds more of them than its size alone predicts. No "
             f"significance is claimed for a row: pathways overlap heavily, "
             f"and the virus-level form of this question is measured in "
             f"pathway_landscape.py, where nothing survived correction on "
             f"this release. Read a high ratio as a lead to examine."),
        ]
        for r in shown:
            names = [self._label(d) for d
                     in r["dependency_members"][:MAX_DEPENDENCY_NAMES]]
            extra = r["n_dependency_members"] - len(names)
            dep = ", ".join(names) + (f", +{extra} more" if extra else "")
            # MARKED IN BOTH LISTS. It was only in the by-density one,
            # so a reader of the list printed FIRST -- the one most
            # likely to be read at all -- saw the pathway without the
            # thing that qualifies it.
            mark = ("  [NAMES THIS VIRUS]" if r.get("names_the_virus")
                    else "")
            lines.append(
                f"  {r['label']} ({r['n_members']} members) — "
                f"{r['n_dependency_members']} dependency factor(s), "
                f"{r['expected_members']:.1f} expected from size "
                f"(ratio {r['density_ratio']:.2f}); "
                f"{r['n_druggable']} druggable: {dep}{mark}")
        if len(rows) > limit:
            lines.append(f"  … {len(rows) - limit:,} further pathway(s) not "
                         f"shown in this ordering.")
        if by_ratio:
            lines.append(
                f"BY DENSITY instead — the same pathways reordered by how "
                f"much they exceed what their size predicts, among the "
                f"{len(eligible):,} that expect at least one dependency "
                f"factor. The list above is ordered by COUNT, which is "
                f"ordered by SIZE, so these are mostly different pathways "
                f"and are the ones the ratio column actually selects:")
            for r in by_ratio:
                mark = ("  [NAMES THIS VIRUS]" if r.get("names_the_virus")
                        else "")
                lines.append(
                    f"  {r['label']} ({r['n_members']} members) — "
                    f"{r['n_dependency_members']} of them, "
                    f"{r['expected_members']:.1f} expected "
                    f"(ratio {r['density_ratio']:.2f}); "
                    f"{r['n_druggable']} druggable{mark}")
            if thin_ratio:
                lines.append(
                    f"  {thin_ratio:,} pathway(s) are excluded from this "
                    f"second list for expecting under one dependency factor, "
                    f"where the ratio is arithmetic rather than evidence.")
        circular = [r for r in rows if r.get("names_the_virus")]
        if circular:
            lines.append(
                f"[NAMES THIS VIRUS] marks {len(circular):,} pathway(s) "
                f"whose own label names {self._label(virus)}. Reactome "
                f"curates those from this virus's literature, which is the "
                f"same literature the CRISPR screens come from, so their "
                f"overlap with its dependency factors is partly "
                f"definitional rather than a finding. They are shown "
                f"rather than dropped -- the biology in them is real -- "
                f"but the ratio on such a row is not evidence of the same "
                f"kind as the ratio on a row beside it.")
        if unplaced:
            lines.append(
                f"{len(unplaced):,} of the dependency factors sit in no "
                f"pathway under the cap and are absent from every row "
                f"above, so this census covers "
                f"{len(placed) / len(deps):.0%} of what the screens found.")
        if oversized:
            lines.append(f"{len(oversized):,} pathway(s) were excluded for "
                         f"exceeding the {cap}-member cap.")

        return ToolResult(True, "virus_pathways", data=rows, verbalised=lines,
                          warnings=base_warnings, edge_ids=edge_ids,
                          note=f"{len(oversized)} pathway(s) over the cap"
                               if oversized else "")

    def evidence(self, subject: str, object_: str | None = None,
                 predicate: str | None = None, limit: int = 25) -> ToolResult:
        """Every edge supporting a specific claim, with provenance."""
        hits = [e for e in self.edges
                if e["subject"] == subject
                and (object_ is None or e["object"] == object_)
                and (predicate is None or e["predicate"] == predicate)][:limit]
        if not hits:
            return ToolResult(False, "evidence",
                              note="No edge in this graph makes that claim.")
        lines = [self._verbalise(e) for e in hits]
        warns = [guards.evidence_note(e.get("evidence_tier"),
                                      (e.get("qualifiers") or {}).get("detection_method"))
                 for e in hits]
        # THE HELD-OUT CAVEAT GOES FIRST. The only warning this tool
        # used to attach for a label edge was evidence_note's "EVIDENCE:
        # Curated experimental result" -- so the answer sheet was
        # returned carrying the tier-1 reassurance and nothing else.
        return ToolResult(True, "evidence",
                          data=[{k: v for k, v in e.items() if not k.startswith("_")}
                                for e in hits],
                          verbalised=lines,
                          warnings=self._held_out_warnings(hits)
                          + list(dict.fromkeys(warns)),
                          edge_ids=[e["_id"] for e in hits])

    def neighbours(self, node_id: str, predicate: str | None = None,
                   limit: int = MAX_NEIGHBOURS) -> ToolResult:
        """Adjacent nodes, capped."""
        if node_id not in self.nodes:
            return ToolResult(False, "neighbours", note=f"{node_id} not in graph")
        hits = [e for e in self.edges
                if (e["subject"] == node_id or e["object"] == node_id)
                and (predicate is None or e["predicate"] == predicate)]
        truncated = len(hits) > limit
        hits = hits[:limit]
        return ToolResult(
            True, "neighbours",
            data=[{"edge": e["predicate"],
                   "other": e["object"] if e["subject"] == node_id else e["subject"]}
                  for e in hits],
            verbalised=[self._verbalise(e) for e in hits],
            # NOT NAMED BY THE AUDIT, and it had the same hole with no
            # warning at all. Found by probing every public tool against
            # one fixture instead of the two the audit listed.
            warnings=self._held_out_warnings(hits),
            edge_ids=[e["_id"] for e in hits],
            note=f"truncated to {limit}" if truncated else "")

    def profile(self, node_id: str) -> ToolResult:
        """What the graph holds about one entity."""
        n = self.nodes.get(node_id)
        if not n:
            return ToolResult(False, "profile", note=f"{node_id} not in graph")
        counts: dict[str, int] = defaultdict(int)
        for e in self.edges:
            if e["subject"] == node_id or e["object"] == node_id:
                counts[e["predicate"]] += 1
        # FEATURES AND LABELS, LISTED APART. One flat list of predicate
        # counts is how "which compounds should be considered for Yellow
        # fever" returned a taxon profile reading "1,030 antiviral
        # measurements" -- which looks like 1,030 supporting facts and is
        # 1,030 rows of the answer sheet. The refusal written for that
        # virus now catches the question (A4), but the profile itself
        # still reported it that way for any question that reaches here,
        # and the catch-all branch of build_brief reaches here often.
        held = {k: v for k, v in counts.items() if k in self._held_out}
        feat = {k: v for k, v in counts.items() if k not in self._held_out}

        def _rows(d):
            return [f"  {k}: {v:,}" for k, v in sorted(d.items(),
                                                       key=lambda kv: -kv[1])]

        lines = [f"{self._label(node_id)} ({n['class']}, {node_id})"]
        lines += _rows(feat)
        if held:
            lines.append("  held-out evaluation labels (the answer sheet, "
                         "not evidence):")
            lines += _rows(held)
        return ToolResult(True, "profile",
                          # edge_counts stays the FULL census, so a
                          # consumer summing it still gets the node's true
                          # degree; the split is additional, not a
                          # redefinition of a key other code reads.
                          data={"id": node_id, "class": n["class"],
                                "properties": n["properties"],
                                "edge_counts": dict(counts),
                                "feature_edge_counts": feat,
                                "held_out_label_counts": held},
                          warnings=[guards.caveat(
                              "held_out_label",
                              guards.held_out_caveat(held))] if held else [],
                          verbalised=lines)
