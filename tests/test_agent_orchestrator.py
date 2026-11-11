"""Gate for the orchestrator: the pipeline is deterministic, the refusal comes
before retrieval, and the model is given no discretion.
"""
import json

import pytest

from kgav.agent.verify import dropped_warnings
from kgav.agent.orchestrator import Orchestrator
from kgav.agent.tools import GraphTools

V = "NCBITaxon:2697049"
DRUG = "INCHIKEY:NIRMA"
CYTO = "INCHIKEY:CYTOX"


@pytest.fixture
def orch(tmp_path):
    def prot(pid, viral, sym=None, fam=None):
        p = {"taxon_id": V if viral else "NCBITaxon:9606", "is_viral": viral,
             "sequence_hash": "h", "reviewed": True}
        if sym:
            p["gene_symbol"] = sym
        if fam:
            p["protein_family"] = fam
        return {"id": pid, "class": "Protein", "properties": p}

    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        prot("UniProtKB:NSP5", True, fam="nsp5"),
        {"id": "KGAV:G", "class": "Gene",
         "properties": {"symbol": "rep", "taxon_id": V, "is_viral": True}},
        prot("UniProtKB:SIGMAR1", False, sym="SIGMAR1"),
        {"id": DRUG, "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "NIRMA", "label": "nirmatrelvir",
                        "chembl_id": "CHEMBL1", "is_approved": True,
                        "salt_collapsed": False, "stereo_collapsed": False, "tautomer_collapsed": True}},
        {"id": CYTO, "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "CYTOX", "label": "chloroquine",
                        "chembl_id": "CHEMBL2", "is_approved": True,
                        "salt_collapsed": False, "stereo_collapsed": False, "tautomer_collapsed": True}},
    ]

    def e(s, p, o, q=None):
        # chembl.py sets target_resolution on EVERY INHIBITS edge, both
        # branches, and the schema requires it -- so M1 constrains on it and
        # qualifier_ok reads a missing qualifier as a failed match. A fixture
        # omitting it is not a simpler fixture; it is an edge the
        # direct-acting channel refuses to walk, in the scorer and now in the
        # agent too.
        if p == "INHIBITS":
            q = dict(q or {})
            q.setdefault("target_resolution", "chain_from_assay_description")
        return {"subject": s, "predicate": p, "object": o, "qualifiers": q or {},
                "primary_knowledge_source": "infores:chembl", "evidence_tier": 1,
                "first_asserted_date": "2021-03-04"}

    edges = [
        e(DRUG, "INHIBITS", "UniProtKB:NSP5", {"assay_type": "biochemical", "ic50_nm": 25.0}),
        e("UniProtKB:NSP5", "ENCODED_BY", "KGAV:G"),
        e("KGAV:G", "BELONGS_TO", V),
        e(CYTO, "TARGETS", "UniProtKB:SIGMAR1",
          {"direction": "inhibitor", "assay_type": "binding"}),
        e("UniProtKB:SIGMAR1", "HOST_FACTOR_FOR", V,
          {"direction": "dependency", "screen_type": "CRISPRko", "cell_line": "CVCL_0574"}),
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))
    return Orchestrator(GraphTools(tmp_path))


# ------------------------------------------------------------------ refusal
def test_prospective_question_is_refused_before_retrieval(orch):
    """No partial answer may leak out alongside a refusal."""
    out = orch.answer("what should we test against SARS-CoV-2?")
    assert out["allowed"] is False
    assert out["facts"] == [] and out["citable_edge_ids"] == []
    # The refusal and its reason, not a figure. "0.50" was pinned here and in
    # guards' own test, so a stale number could not be removed from
    # user-facing prose without the suite failing.
    assert "cannot support prospective" in out["answer"]


def test_refusal_happens_even_when_entities_resolve(orch):
    out = orch.answer("recommend a drug for SARS-CoV-2 like nirmatrelvir")
    assert out["allowed"] is False


# ------------------------------------------------------------- determinism
def test_same_question_gives_the_same_brief(orch):
    """A researcher who runs a query twice and gets two answers cannot cite
    the tool."""
    a = orch.answer("why might nirmatrelvir work against SARS-CoV-2?")
    b = orch.answer("why might nirmatrelvir work against SARS-CoV-2?")
    assert a["facts"] == b["facts"]
    assert a["citable_edge_ids"] == b["citable_edge_ids"]


# --------------------------------------------------------------- resolution
def test_longest_match_wins(orch):
    """hydroxychloroquine must not resolve as chloroquine -- the same
    substring hazard as Vero E6 inside Vero."""
    orch.tools.names["hydroxychloroquine"] = "INCHIKEY:HCQ"
    orch.tools.nodes["INCHIKEY:HCQ"] = {
        "id": "INCHIKEY:HCQ", "class": "SmallMolecule",
        "properties": {"label": "hydroxychloroquine"}}
    brief = orch.build_brief("explain hydroxychloroquine")
    assert "INCHIKEY:HCQ" in brief.resolved.values()
    assert CYTO not in brief.resolved.values()


def test_unresolvable_question_asks_rather_than_guesses(orch):
    out = orch.answer("explain the mechanism of aspirin")
    assert out["allowed"] is True
    assert "does not match" in out["answer"] or "No entity" in out["answer"]


# ------------------------------------------------------------------ routing
def test_explain_intent_retrieves_a_mechanistic_path(orch):
    out = orch.answer("why might nirmatrelvir work against SARS-CoV-2?")
    assert out["intent"] == "explain"
    assert any("M1 direct-acting" in f for f in out["facts"])
    assert out["citable_edge_ids"]


def test_triage_intent_ranks_the_named_compounds(orch):
    out = orch.answer("rank nirmatrelvir and chloroquine for SARS-CoV-2")
    assert out["intent"] == "triage"
    assert any("nirmatrelvir" in f for f in out["facts"])


def test_profile_intent_summarises(orch):
    out = orch.answer("what is known about SIGMAR1")
    assert out["intent"] == "profile"
    assert out["facts"]


# ----------------------------------------------------------------- warnings
def test_warnings_survive_into_the_brief(orch):
    out = orch.answer("why might chloroquine work against SARS-CoV-2?")
    joined = " ".join(out["facts"]) + " ".join(out["warnings"])
    assert "CELL CONTEXT" in joined or "SELECTIVITY UNKNOWN" in joined


def test_prompt_requires_every_warning_to_appear(orch):
    brief = orch.build_brief("why might chloroquine work against SARS-CoV-2?")
    p = brief.to_prompt()
    assert "MUST appear in your answer" in p
    assert "Use ONLY the facts listed below" in p


def test_prompt_lists_only_retrieved_edge_ids(orch):
    brief = orch.build_brief("why might nirmatrelvir work against SARS-CoV-2?")
    p = brief.to_prompt()
    for eid in brief.citable_edge_ids:
        assert eid in p


# -------------------------------------------------------------------- model
def test_model_is_optional_and_receives_only_the_brief(orch):
    seen = {}

    def fake_model(prompt):
        seen["prompt"] = prompt
        return "answer text"

    out = orch.answer("why might nirmatrelvir work against SARS-CoV-2?", fake_model)
    assert out["answer"] == "answer text"
    assert "FACTS:" in seen["prompt"]
    # The model must not be handed a query language or the graph itself.
    # Checking for "SELECT" alone false-positives on "SELECTIVITY UNKNOWN" --
    # the fourth substring collision in this project, after SARS-CoV inside
    # SARS-CoV-2, Vero inside Vero E6, and chloroquine inside
    # hydroxychloroquine.
    p = seen["prompt"]
    assert "SELECT " not in p.upper().replace("SELECTIVITY", "")
    assert "jsonl" not in p and "MATCH (" not in p


def test_no_model_returns_the_brief_unsynthesised(orch):
    out = orch.answer("why might nirmatrelvir work against SARS-CoV-2?")
    assert out["answer"] is None and "prompt" in out


def test_trace_records_each_stage(orch):
    out = orch.answer("why might nirmatrelvir work against SARS-CoV-2?")
    joined = " ".join(out["trace"])
    assert "intent=" in joined and "resolved=" in joined and "allowed=" in joined


def test_dropping_the_cell_context_caveat_fails_verification(orch):
    """End to end: the caveat must survive from retrieval to the verdict.

    The gap this closes sat exactly here. test_agent_tools asserted the
    warning text reaches the printed route; test_agent_verify asserted
    verify() enforces a caveat handed to it directly. Neither asked whether
    explain() puts the caveat anywhere verify() reads, and it did not -- so a
    model that dropped the SIGMAR1 warning passed, for the one caveat this
    project exists to attach.
    """
    from kgav.agent.verify import verify_response

    out = orch.answer("why might chloroquine work against SARS-CoV-2?")
    assert any("CELL CONTEXT" in w for w in out["warnings"]), \
        "the caveat never reached the brief, so nothing downstream can require it"

    dropped = "Chloroquine engages SIGMAR1 [E0]. SELECTIVITY UNKNOWN: no paired data."
    v = verify_response(out, answer=dropped)
    assert not v.passed
    assert any(f.code == "warning_dropped" for f in v.failures)


def test_truncation_drops_edge_ids_with_their_facts(orch):
    """A citable id whose fact was cut is a licence to invent.

    verify()'s invented-citation check tests membership of
    citable_edge_ids, so an id left behind by truncation passes it -- the
    model can attach a real identifier to a claim with no retrieved evidence
    and the answer reports as verified.
    """
    from kgav.agent import orchestrator as orch_mod

    original = orch_mod.MAX_BRIEF_EDGES
    orch_mod.MAX_BRIEF_EDGES = 1
    try:
        # One result per compound, so there is a second result to drop. The
        # first is always kept whole: cutting inside a result would separate
        # facts from their ids again, and an empty brief is worse than a long
        # one.
        brief = orch.build_brief(
            "why might nirmatrelvir and chloroquine work against SARS-CoV-2?")
    finally:
        orch_mod.MAX_BRIEF_EDGES = original

    assert "omitted" in " ".join(brief.facts)
    # Every surviving id must appear in a surviving fact.
    joined = " ".join(brief.facts)
    assert all(eid in joined for eid in brief.citable_edge_ids)


def test_truncation_never_drops_a_warning(orch):
    """Facts are bounded; caveats are not. A caveat cut for space is the
    failure mode guards.py exists to prevent."""
    from kgav.agent import orchestrator as orch_mod

    original = orch_mod.MAX_BRIEF_EDGES
    orch_mod.MAX_BRIEF_EDGES = 1
    try:
        brief = orch.build_brief("why might chloroquine work against SARS-CoV-2?")
    finally:
        orch_mod.MAX_BRIEF_EDGES = original
    assert any("CELL CONTEXT" in w for w in brief.warnings)


# ------------------------------------------- one retry, for a dropped caveat
class _Scripted:
    """A model that returns prepared answers in order and records its prompts."""

    def __init__(self, *answers):
        self.answers, self.prompts = list(answers), []

    def __call__(self, prompt):
        self.prompts.append(prompt)
        return self.answers[min(len(self.prompts) - 1, len(self.answers) - 1)]


def _orch_with_caveat(tmp_path):
    """An orchestrator whose brief carries a required caveat."""
    from kgav.agent import guards
    o = Orchestrator.__new__(Orchestrator)

    class _Brief:
        allowed, clarification, refusal = True, None, ""
        intent, trace, facts = "explain", [], ["a fact [E1]"]
        warnings = [guards.caveat("cell_context", "SIGMAR1 route is cell-type "
                                                  "dependent.")]
        citable_edge_ids = ["E1"]
        question = "why?"

        def to_prompt(self, must_include=None):
            head = ""
            if must_include:
                head = ("REJECTED, omitted: "
                        + ";".join(str(w)[:20] for w in must_include) + "\n")
            return head + "PROMPT"

    o.build_brief = lambda q: _Brief()
    return o


def test_a_dropped_caveat_triggers_exactly_one_retry(tmp_path):
    """Withholding was right and made itself the normal outcome: a 7B model
    drops a caveat for brevity on nearly every host-directed question."""
    m = _Scripted("no caveat here", "CELL CONTEXT: SIGMAR1 is cell dependent [E1]")
    out = _orch_with_caveat(tmp_path).answer("why?", m)
    assert out["attempts"] == 2
    assert len(m.prompts) == 2
    assert "REJECTED" in m.prompts[1] and "REJECTED" not in m.prompts[0]
    assert dropped_warnings(out["answer"], _orch_with_caveat(tmp_path)
                            .build_brief("x").warnings) == []


def test_a_compliant_first_answer_is_not_retried(tmp_path):
    m = _Scripted("CELL CONTEXT: it depends on the cell type [E1]")
    out = _orch_with_caveat(tmp_path).answer("why?", m)
    assert out["attempts"] == 1 and len(m.prompts) == 1
    assert "retry_prompt" not in out


def test_the_retry_is_not_unlimited(tmp_path):
    """Re-rolling until a check passes is how a verifier becomes a sampler."""
    m = _Scripted("still no caveat", "still no caveat", "still no caveat")
    out = _orch_with_caveat(tmp_path).answer("why?", m)
    assert out["attempts"] == 2 and len(m.prompts) == 2
    # and the answer is still non-compliant, so the verifier will withhold it
    assert "cell context" not in out["answer"].lower()


def test_max_attempts_one_disables_the_retry(tmp_path):
    m = _Scripted("no caveat")
    out = _orch_with_caveat(tmp_path).answer("why?", m, max_attempts=1)
    assert out["attempts"] == 1 and len(m.prompts) == 1


def test_what_was_dropped_is_recorded(tmp_path):
    m = _Scripted("no caveat", "CELL CONTEXT: fine [E1]")
    out = _orch_with_caveat(tmp_path).answer("why?", m)
    assert out["dropped_on_first_attempt"]
    assert "CELL CONTEXT" in out["dropped_on_first_attempt"][0]


def test_the_retry_prompt_names_the_specific_label():
    """The rules already say every warning must appear, and the model dropped
    one anyway. Repeating that instruction is not a different instruction."""
    from kgav.agent import guards
    from kgav.agent.orchestrator import Brief
    b = Brief(question="why?", intent="explain", allowed=True,
              facts=["f [E1]"], citable_edge_ids=["E1"],
              warnings=[guards.caveat("cell_context", "depends on cell type")])
    first, again = b.to_prompt(), b.to_prompt(must_include=b.warnings)
    assert "CELL CONTEXT" in again.split("RULES:")[0]      # at the top
    assert "REJECTED" in again and "REJECTED" not in first


# ------------------------------------------- an assumed virus is said out loud
def test_a_question_naming_no_virus_says_which_one_it_assumed(orch):
    """build_brief reads `self._virus_in(resolved) or self.default_virus`, so
    a question naming only a compound was answered about SARS-CoV-2 with
    nothing marking the choice. The virus decides the calibration, the
    caveats and the cell-context biology, so a question meant about Zika
    answered about a coronavirus is wrong in all three at once. Latent while
    the graph held one virus; live since it holds eight across two families.
    """
    brief = orch.build_brief("why might nirmatrelvir work?")
    assert brief.allowed
    assert any(getattr(w, "key", "") == "virus_assumed" for w in brief.warnings)
    assumed = next(w for w in brief.warnings
                   if getattr(w, "key", "") == "virus_assumed")
    assert "VIRUS ASSUMED" in assumed
    # NOT "named no virus": that asserts a reason, and the reason was wrong
    # in the case that mattered -- "for dengue virus 2" resolved to nothing
    # and was answered about SARS-CoV-2 under a sentence saying the question
    # had named no virus, while listing DENV-2 two clauses later.
    assert "no virus in this question resolved" in assumed
    assert "named no virus" not in assumed
    assert any("(assumed)" in t for t in brief.trace)


def test_naming_the_virus_adds_no_assumption_caveat(orch):
    """So the caveat marks an assumption rather than appearing on every
    answer -- a warning that always fires is one a reader learns to skip."""
    brief = orch.build_brief("why might nirmatrelvir work against SARS-CoV-2?")
    assert brief.allowed
    assert not any(getattr(w, "key", "") == "virus_assumed"
                   for w in brief.warnings)


def test_the_assumption_caveat_is_enforced_like_every_other(orch):
    """It is in CAVEAT_LABELS, so verify.dropped_warnings requires its label
    in the answer. A caveat the model may drop for brevity is not a caveat."""
    from kgav.agent.verify import caveat_labels, dropped_warnings
    brief = orch.build_brief("why might nirmatrelvir work?")
    assert "VIRUS ASSUMED" in caveat_labels(brief.warnings)
    assert dropped_warnings("an answer that mentions nothing", brief.warnings)
    kept = "VIRUS ASSUMED: this is about SARS-CoV-2. " + " ".join(
        str(w) for w in brief.warnings)
    assert not any(getattr(w, "key", "") == "virus_assumed"
                   for w in dropped_warnings(kept, brief.warnings))


# ------------- routes are counted from the graph being served, not recalled
#: A second virus taxon, so a cross-virus viral pair can be built.
V_OTHER = "NCBITaxon:694009"


def _ppi_tools(tmp_path, pairs, reach=None):
    """A release holding only PHYSICALLY_INTERACTS_WITH, from `pairs`.

    Built on disk and read back through GraphTools, so routes_for is
    exercised against a real release rather than a hand-made object.
    """
    def prot(pid, viral):
        return {"id": pid, "class": "Protein",
                "properties": {"taxon_id": V if viral else "NCBITaxon:9606",
                               "is_viral": viral, "sequence_hash": "h",
                               "reviewed": True}}

    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        prot("UniProtKB:NSP5", True), prot("UniProtKB:NSP12", True),
        prot("UniProtKB:SIGMAR1", False), prot("UniProtKB:TMEM41B", False),
        # A SECOND VIRUS, because one virus cannot produce a CROSS-virus
        # viral pair and that is the shape the old check counted as a
        # host route. Every fixture here held one taxon, so the bug was
        # unreachable from this file.
        {"id": V_OTHER, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        {"id": "UniProtKB:SARS1_NSP5", "class": "Protein",
         "properties": {"taxon_id": V_OTHER, "is_viral": True,
                        "sequence_hash": "h", "reviewed": True}},
    ]
    edges = [{"subject": s, "predicate": "PHYSICALLY_INTERACTS_WITH",
              "object": o, "qualifiers": {},
              "primary_knowledge_source": "infores:virhostnet",
              "evidence_tier": 2, "first_asserted_date": "2020-06-01"}
             for s, o in pairs]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(x) for x in edges))
    return GraphTools(tmp_path, reach=reach)


def test_routes_are_counted_from_the_release(tmp_path):
    """A virus whose routes vanish in a rebuild must stop being
    answerable in the SAME run. A figure carried from a results file
    cannot do that."""
    from kgav.agent import guards
    routes = _ppi_tools(
        tmp_path, [("UniProtKB:NSP5", "UniProtKB:SIGMAR1")]
    ).routes_for(V)
    assert isinstance(routes, guards.VirusRoutes)
    assert routes.ppi == 1 and routes.routes == ("ppi",)


def test_a_host_host_edge_is_not_a_virus_route(tmp_path):
    """PHYSICALLY_INTERACTS_WITH is shared with the STRING layer, so a
    host-host pair carries the same predicate. Counting one would give a
    virus a route it does not have -- and this route is exactly what
    suppresses the no-route refusal."""
    routes = _ppi_tools(
        tmp_path, [("UniProtKB:SIGMAR1", "UniProtKB:TMEM41B")]).routes_for(V)
    assert routes.ppi == 0
    assert routes.routes == ()


def test_a_viral_viral_edge_is_not_a_virus_route(tmp_path):
    """Both endpoints are this virus, so `virus in (a, b)` would count it.
    The edge says nothing about a host protein, which is what the PPI
    channel scores, so the check is exactly-one-endpoint."""
    routes = _ppi_tools(
        tmp_path, [("UniProtKB:NSP5", "UniProtKB:NSP12")]).routes_for(V)
    assert routes.ppi == 0
    assert routes.routes == ()


def test_a_cross_virus_viral_pair_is_not_a_virus_route(tmp_path):
    """THE SHAPE THE OLD CHECK COUNTED. `(a == virus) != (b == virus)`
    over a function returning None for host proteins excludes host-host
    (None, None) and SAME-virus viral pairs (virus, virus) -- and
    counts (virus, otherVirus), because True != False.

    The comment above it claimed viral-viral was excluded "by the same
    check". It was excluded only when both proteins belonged to the
    virus being asked about, and every fixture in this file held one
    taxon. The register records SARS-CoV's 73 VirHostNet interactions as
    viral-viral, so this is live data, and the count decides whether the
    agent refuses to rank at all. Found by external audit, 2026-11-02.
    """
    routes = _ppi_tools(
        tmp_path, [("UniProtKB:NSP5", "UniProtKB:SARS1_NSP5")]).routes_for(V)
    assert routes.ppi == 0, (
        "a viral protein of another virus was counted as a host partner")
    assert routes.routes == ()
    # And symmetrically for the other virus, which must not inherit it.
    other = _ppi_tools(
        tmp_path,
        [("UniProtKB:NSP5", "UniProtKB:SARS1_NSP5")]).routes_for(V_OTHER)
    assert other.ppi == 0


def test_another_viruss_protein_meeting_a_host_protein_is_not_our_route(
        tmp_path):
    """The other half of the cross-virus case, and a surviving mutant
    until it was written. Dropping `virus in (a, b)` leaves
    "exactly one endpoint is viral", which is true of SARS-CoV's NSP5
    bound to a human protein -- a real virus-host interaction, and
    nothing at all to do with SARS-CoV-2."""
    routes = _ppi_tools(
        tmp_path,
        [("UniProtKB:SARS1_NSP5", "UniProtKB:SIGMAR1")]).routes_for(V)
    assert routes.ppi == 0, (
        "another virus's interaction was counted as a route for this one")
    assert routes.routes == ()
    # It IS a route for the virus it belongs to.
    assert _ppi_tools(
        tmp_path,
        [("UniProtKB:SARS1_NSP5", "UniProtKB:SIGMAR1")]
    ).routes_for(V_OTHER).ppi == 1


def test_only_the_host_viral_edge_of_a_mixed_release_is_counted(tmp_path):
    """All four shapes in one release: the count is 1, not 4."""
    routes = _ppi_tools(tmp_path, [
        ("UniProtKB:NSP5", "UniProtKB:SIGMAR1"),        # host-viral
        ("UniProtKB:NSP5", "UniProtKB:NSP12"),          # viral-viral, same
        ("UniProtKB:NSP5", "UniProtKB:SARS1_NSP5"),     # viral-viral, cross
        ("UniProtKB:SIGMAR1", "UniProtKB:TMEM41B"),     # host-host
    ]).routes_for(V)
    assert routes.ppi == 1


# ------------- a name inside a longer word is not a mention of that name
@pytest.mark.parametrize("text,name,want", [
    # THE LIVE FAILURE. "Depen" is a brand name for penicillamine, the
    # index holds brand names, and bare substring matching resolved it
    # out of the word "depend" -- so "which pathways does SARS-CoV-2
    # depend on?" was answered about penicillamine, fully caveated, as
    # though that were the question.
    ("which pathways does sars-cov-2 depend on?", "depen", False),
    ("which pathways does sars-cov-2 depend on?", "sars-cov-2", True),
    # The hazard the longest-first rule was written for, now enforced
    # by the matcher rather than by ordering alone.
    ("does hydroxychloroquine work", "chloroquine", False),
    ("does hydroxychloroquine work", "hydroxychloroquine", True),
    ("chloroquine and hydroxychloroquine", "chloroquine", True),
    # Punctuation is a boundary; digits and hyphens inside a name are
    # not, which is why \b is the wrong assertion for this index.
    ("rank remdesivir.", "remdesivir", True),
    ("about sars-cov-2, and", "sars-cov-2", True),
    ("the sars-cov-23 strain", "sars-cov-2", False),
    ("p53-dependent genes", "p53", True),
    ("vero e6 cells", "vero e6", True),
])
def test_a_name_must_sit_on_a_word_boundary(text, name, want):
    from kgav.agent.orchestrator import _mentions
    assert _mentions(text, name) is want


def test_a_brand_name_inside_an_english_word_is_not_resolved(tmp_path):
    """End to end, through the resolver rather than the helper: a
    question naming only a virus must resolve only the virus."""
    import json

    from kgav.agent.orchestrator import Orchestrator
    from kgav.agent.tools import GraphTools

    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        {"id": "INCHIKEY:DEPEN", "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "DEPEN",
                        "label": "depen", "chembl_id": "CHEMBL7",
                        "is_approved": True, "salt_collapsed": False,
                        "stereo_collapsed": False, "tautomer_collapsed": True}},
    ]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("")
    orch = Orchestrator(GraphTools(tmp_path))
    resolved, _problem, _omitted = orch._resolve_mentions(
        "which pathways does SARS-CoV-2 depend on?")
    assert set(resolved.values()) == {V}, (
        f"the question names one entity and resolved {resolved}")
    # And the compound IS resolved when it is actually named.
    named, _, _ = orch._resolve_mentions("what is depen")
    assert "INCHIKEY:DEPEN" in named.values()


def test_the_suppression_check_keeps_the_boundary_rule(tmp_path):
    """Both tests in the resolver have to be boundary-aware, not just
    the first.

    Masking a longer match can delete the SHORTER name's only
    boundary occurrence while leaving a non-boundary one elsewhere --
    the index holds multi-word names, so a short name can sit
    boundary-clean inside a long one. Here "ace2 receptor" is found and
    masked, which removes the clean "ace2"; what remains is "ace2x",
    and a bare `name in covered` resolves ACE2 out of it.
    """
    import json

    from kgav.agent.orchestrator import Orchestrator
    from kgav.agent.tools import GraphTools

    def prot(pid, label):
        return {"id": pid, "class": "Protein",
                "properties": {"taxon_id": "NCBITaxon:9606",
                               "is_viral": False, "sequence_hash": "h",
                               "reviewed": True, "gene_symbol": label}}

    nodes = [prot("UniProtKB:ACE2", "ace2"),
             prot("UniProtKB:ACE2R", "ace2 receptor")]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("")
    orch = Orchestrator(GraphTools(tmp_path))
    resolved, _, _ = orch._resolve_mentions("the ace2 receptor in ace2x cells")
    assert set(resolved.values()) == {"UniProtKB:ACE2R"}, (
        f"ace2 was resolved out of 'ace2x' after masking: {resolved}")


# ------------- the no-route refusal belongs to the TOOL, not to its caller
#
# guards.check() refuses a routed query only when `routes is not None`,
# and three of the four tools omitted the argument: pathways(), explain()
# and triage() passed only the calibration. So the refusal written for
# Yellow fever and DENV-4 -- labels present, no metapath terminating at
# them, nothing rankable -- was a property of the ORCHESTRATOR, which
# happens to compute routes itself, rather than of the tool. Reproduced
# by the external audit of 2026-11-02: tools.pathways on a zero-route
# virus returned ok=True.
#
# Every test of the refusal went through the orchestrator, so the gap was
# invisible from this file.
def _labelled_no_route(tmp_path) -> GraphTools:
    """A virus in the Yellow fever tier: antiviral labels, zero routes.

    No INHIBITS, no HOST_FACTOR_FOR and no virus-host PPI, so no metapath
    terminates at it. The compound has a human target and that target
    sits in a pathway, so pathways() has material to answer with and the
    refusal is the only thing stopping it.
    """
    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "Yellow fever", "family": "Flaviviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        {"id": "INCHIKEY:DRUGA", "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "DRUGA",
                        "chembl_id": "CHEMBL1"}},
        {"id": "UniProtKB:H1", "class": "Protein",
         "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False,
                        "sequence_hash": "h", "reviewed": True}},
        {"id": "REACT:P1", "class": "Pathway", "properties": {"label": "p1"}},
    ]
    edges = [
        {"subject": "INCHIKEY:DRUGA", "predicate": "TARGETS",
         "object": "UniProtKB:H1", "qualifiers": {},
         "primary_knowledge_source": "infores:chembl", "evidence_tier": 1,
         "first_asserted_date": "2020-01-01"},
        {"subject": "UniProtKB:H1", "predicate": "PARTICIPATES_IN",
         "object": "REACT:P1", "qualifiers": {},
         "primary_knowledge_source": "infores:uniprot", "evidence_tier": 1,
         "first_asserted_date": "2020-01-01"},
        # The labels the refusal discloses, and the reason the tier exists.
        {"subject": "INCHIKEY:DRUGA",
         "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST", "object": V,
         "qualifiers": {}, "primary_knowledge_source": "infores:chembl",
         "evidence_tier": 1, "first_asserted_date": "2020-01-01"},
    ]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in edges))
    return GraphTools(tmp_path)


def test_the_fixture_really_has_no_route(tmp_path):
    """GUARDS EVERY TEST BELOW. If the fixture acquires a route, the
    three refusals stop being refusals and all of them pass vacuously."""
    routes = _labelled_no_route(tmp_path).routes_for(V)
    assert routes.routes == (), routes
    assert routes.n_labels == 1, "the fixture must hold a label to disclose"


@pytest.mark.parametrize("call", [
    lambda t: t.pathways("INCHIKEY:DRUGA", V),
    lambda t: t.explain("INCHIKEY:DRUGA", V),
    lambda t: t.triage(["INCHIKEY:DRUGA"], V),
])
def test_each_tool_refuses_a_no_route_virus_on_its_own(tmp_path, call):
    """CALLED DIRECTLY, with no orchestrator in the way. A guard that
    only fires when one particular caller supplies an optional argument
    is a guard any second entry point loses silently."""
    result = call(_labelled_no_route(tmp_path))
    assert result.ok is False, (
        "the tool answered a virus no metapath terminates at; the "
        "refusal is coming from the caller, not from the tool")
    # THE REASON, NOT JUST THE FLAG. explain() returns ok=False for a
    # compound with no path too, so `ok is False` alone passed while the
    # verdict was computed and discarded -- and the user was told "most
    # compounds have an activity measurement and no target annotation",
    # a sentence about the COMPOUND, for a virus nothing can reach.
    assert "none can be reached" in (result.note or ""), (
        f"refused for the wrong reason:\n{result.note}")


def test_the_refusal_says_it_cannot_look_rather_than_that_it_found_nothing(
        tmp_path):
    """The distinction the whole module exists to keep. An empty result
    reads as 'no compound looks promising', which is the opposite of
    what is true."""
    result = _labelled_no_route(tmp_path).pathways("INCHIKEY:DRUGA", V)
    # The reason travels in `note`; verbalised is empty on a refusal,
    # which is itself the refusal-before-retrieval property.
    assert result.verbalised == [], (
        "a refusal carried an answer alongside it")
    assert "none can be reached" in result.note, result.note
    assert "not because no compound looks promising" in result.note
    # AND IT DISCLOSES THE LABELS IT IS SITTING ON. A refusal that hides
    # how much evidence exists reads as an empty graph.
    assert "antiviral measurement" in result.note, result.note


# ------------------------------------------- measured reach, plumbed in
# routes_for counts every other field FROM THE GRAPH. Reach cannot be
# counted from edges -- walking a metapath against a label set is what
# produces it -- so it arrives from labels.measured_reach. These assert
# the join: keyed by DISPLAY LABEL, and a virus the map does not cover
# resolves to None rather than to 0.
def test_reach_is_joined_by_display_label(tmp_path):
    routes = _ppi_tools(
        tmp_path, [("UniProtKB:NSP5", "UniProtKB:SIGMAR1")],
        reach={"SARS-CoV-2": 1085}).routes_for(V)
    assert routes.label == "SARS-CoV-2"
    assert routes.reached_labels == 1085


def test_a_virus_absent_from_the_reach_map_is_unmeasured(tmp_path):
    """NOT zero. The map holds the viruses metapath_comparison.json
    covers; a virus it omits has not been measured, and reporting that
    as "no route reaches a label" is the opposite claim."""
    routes = _ppi_tools(
        tmp_path, [("UniProtKB:NSP5", "UniProtKB:SIGMAR1")],
        reach={"DENV-1": 0}).routes_for(V)
    assert routes.reached_labels is None
    assert not routes.reach_measured


def test_no_reach_map_leaves_every_virus_unmeasured(tmp_path):
    """The default, and the behaviour of every caller that cannot
    compute reach: the agent says so rather than refusing or asserting."""
    tools = _ppi_tools(tmp_path, [("UniProtKB:NSP5", "UniProtKB:SIGMAR1")])
    assert tools.routes_for(V).reached_labels is None
    assert tools.routes_for(V_OTHER).reached_labels is None


def test_a_reach_of_zero_survives_the_join(tmp_path):
    """0 is a VALUE here, not an absence. A join written with `or None`
    or a truthiness guard would turn the one tier this plumbing exists
    for back into "unmeasured".

    `unreachable` is deliberately NOT asserted: this fixture writes no
    HAS_ANTIVIRAL_ACTIVITY_AGAINST edge, so n_labels is 0 and the virus
    is in the FEATURES, NOTHING TO SCORE tier, where reach of 0 means
    there was nothing to reach. The first version of this test asserted
    it anyway and failed -- a fixture that cannot produce the condition
    it claims to check, which is the defect class this suite keeps
    finding. The tier itself is checked in test_agent_guards.py, where
    the label count can be set."""
    routes = _ppi_tools(
        tmp_path, [("UniProtKB:NSP5", "UniProtKB:SIGMAR1")],
        reach={"SARS-CoV-2": 0}).routes_for(V)
    assert routes.reached_labels == 0
    assert routes.reach_measured
    assert routes.n_labels == 0 and not routes.unreachable


# ------------------------------------------- the ids are bounded too
# A13. MAX_BRIEF_EDGES capped `facts` and nothing capped
# `citable_edge_ids`. Measured on the fixture below before the fix:
# pathways() returns 29 census lines and 1,500 PARTICIPATES_IN ids, and
# Brief.to_prompt joins all 1,500 into one line -- 9,762 of the prompt's
# 14,579 characters, 67% of it, bare identifiers. That is the
# context-window problem MAX_BRIEF_EDGES exists to prevent, in the half
# of the brief it did not cover.
@pytest.fixture
def id_heavy(tmp_path):
    """A graph whose pathway census returns far more ids than facts."""
    def e(s, p, o, **q):
        return {"subject": s, "predicate": p, "object": o, "qualifiers": q,
                "primary_knowledge_source": "infores:chembl",
                "evidence_tier": 1, "first_asserted_date": "2021-01-01"}
    nodes = [{"id": V, "class": "OrganismTaxon",
              "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                             "baltimore_class": "IV", "is_enveloped": True}},
             {"id": "CHEMBL:BIG", "class": "SmallMolecule",
              "properties": {"label": "promiscuousdrug"}},
             {"id": "UniProtKB:NSP5", "class": "Protein",
              "properties": {"taxon_id": V, "is_viral": True,
                             "sequence_hash": "hv", "reviewed": True,
                             "gene_symbol": "nsp5"}}]
    edges = [e("CHEMBL:BIG", "INHIBITS", "UniProtKB:NSP5", ic50_nm=19)]
    for i in range(300):
        pid = f"UniProtKB:H{i}"
        nodes.append({"id": pid, "class": "Protein",
                      "properties": {"taxon_id": "NCBITaxon:9606",
                                     "is_viral": False,
                                     "sequence_hash": f"h{i}",
                                     "reviewed": True,
                                     "gene_symbol": f"GENE{i}"}})
        edges += [e(pid, "HOST_FACTOR_FOR", V, direction="dependency",
                    screen_type="CRISPRko"),
                  e("CHEMBL:BIG", "TARGETS", pid, direction="inhibitor")]
    for j in range(400):
        pw = f"REACT:{j}"
        nodes.append({"id": pw, "class": "Pathway",
                      "properties": {"label": f"pathway {j}"}})
        for k in range(3):
            edges.append(e(f"UniProtKB:H{(j * 3 + k) % 300}",
                           "PARTICIPATES_IN", pw))
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(x) for x in edges))
    return Orchestrator(GraphTools(tmp_path))


Q_PATHWAY = "which pathways does promiscuousdrug sit in for SARS-CoV-2"


def test_the_citable_ids_are_bounded(id_heavy):
    from kgav.agent.orchestrator import MAX_BRIEF_IDS
    b = id_heavy.build_brief(Q_PATHWAY)
    assert len(b.citable_edge_ids) <= MAX_BRIEF_IDS, len(b.citable_edge_ids)
    assert any(t.startswith("truncated_ids=") for t in b.trace), b.trace


def test_bounding_the_ids_bounds_the_prompt(id_heavy):
    """The point of the cap. Two thirds of this prompt was identifiers."""
    b = id_heavy.build_brief(Q_PATHWAY)
    prompt = b.to_prompt()
    ids_line = next(ln for ln in prompt.splitlines()
                    if ln.startswith("CITABLE EDGE IDS"))
    assert len(ids_line) < len(prompt) * 0.6, (
        f"identifiers are {len(ids_line) / len(prompt):.0%} of the prompt")


def test_the_ids_in_a_brief_are_unique(id_heavy):
    ids = id_heavy.build_brief(Q_PATHWAY).citable_edge_ids
    assert len(ids) == len(set(ids))


def test_shown_edge_ids_reads_both_citation_shapes():
    """A FACT carries "[chembl, PMID:1, 2021; E1]" and a MODEL writes
    "[E1]". Reading only the second finds nothing in any real fact,
    which is the bug the pre-existing verify tests caught when this
    change was first written."""
    from kgav.agent.orchestrator import shown_edge_ids
    assert shown_edge_ids(
        ["nirmatrelvir inhibits nsp5 [chembl, PMID:1, 2021; E335139]"]
    ) == {"E335139"}
    assert shown_edge_ids(["a claim [E7, E8]"]) == {"E7", "E8"}
    assert shown_edge_ids(["342 pathway(s) contain a factor."]) == set()


def test_an_id_whose_fact_was_shown_survives_the_trim():
    """TRIMMED BY USEFULNESS, not by iteration order -- the shown ids go
    first even when they are LAST in the input.

    TESTED DIRECTLY, not through build_brief. The first version of this
    went through the pathway fixture, where the census lines carry no
    inline ids at all: `shown` was empty, so the assertion was
    vacuously true and a mutation replacing the whole preference with
    `ids[:MAX_BRIEF_IDS]` survived it. A fixture that cannot produce the
    condition it checks, again."""
    from kgav.agent.orchestrator import MAX_BRIEF_IDS, bounded_ids
    facts = ["nirmatrelvir inhibits nsp5 [chembl, 2021; E9999]"]
    ids = [f"E{i}" for i in range(MAX_BRIEF_IDS * 2)] + ["E9999"]
    kept, dropped = bounded_ids(facts, ids)
    assert len(kept) == MAX_BRIEF_IDS
    assert dropped == len(ids) - MAX_BRIEF_IDS
    assert "E9999" in kept, "the one id with a fact behind it was trimmed"
    assert kept[0] == "E9999", "shown ids must be kept FIRST"


def test_bounded_ids_leaves_a_short_list_alone():
    from kgav.agent.orchestrator import bounded_ids
    kept, dropped = bounded_ids(["a fact [chembl, 2021; E1]"], ["E1", "E2"])
    assert kept == ["E1", "E2"] and dropped == 0


def test_bounded_ids_deduplicates():
    from kgav.agent.orchestrator import bounded_ids
    kept, _ = bounded_ids([], ["E1", "E2", "E1", "E2"])
    assert kept == ["E1", "E2"]


def test_a_brief_under_the_cap_is_left_alone(id_heavy):
    b = id_heavy.build_brief("what supports that promiscuousdrug inhibits nsp5")
    assert not any(t.startswith("truncated_ids=") for t in b.trace), b.trace


def test_no_tool_can_exceed_the_fact_bound(id_heavy):
    """THE FIRST RESULT IS ADMITTED WHOLE -- `if facts and ...` -- which
    the audit read as unbounded. It is bounded by arithmetic instead:
    the largest single result is explain()'s MAX_PATHS paths of one
    header plus M1's three edges, and _find_paths walks only M1 and M2.
    This pins that, so raising a per-tool cap past the brief bound fails
    here rather than silently doubling the prompt."""
    from kgav.agent.orchestrator import MAX_BRIEF_EDGES
    from kgav.agent import tools as T
    assert T.MAX_PATHS * 4 <= MAX_BRIEF_EDGES
    assert T.MAX_NEIGHBOURS <= MAX_BRIEF_EDGES
    assert T.MAX_PATHWAYS + T.MAX_RATIO_ROWS + 10 <= MAX_BRIEF_EDGES
    assert T.MAX_SUBGRAPH_EDGES <= MAX_BRIEF_EDGES


def test_an_over_long_brief_is_recorded_rather_than_emptied(id_heavy):
    """Dropping the first result would leave an EMPTY brief, which is
    worse than an over-long one. So it is admitted and the overflow is
    recorded where it can be seen."""
    import kgav.agent.orchestrator as O
    before = O.MAX_BRIEF_EDGES
    try:
        O.MAX_BRIEF_EDGES = 3
        b = id_heavy.build_brief(Q_PATHWAY)
        assert b.facts, "the brief was emptied rather than reported"
        assert any(t.startswith("brief_over_fact_bound=") for t in b.trace), \
            b.trace
    finally:
        O.MAX_BRIEF_EDGES = before


# ----------------------------- the cap must not decide WHICH kind survives
# A14. `if len(found) >= 6: break` sat inside a loop ordered by name
# LENGTH. Drug names are long and virus names are short, so seven named
# compounds filled the quota before "zika" (four characters) was
# reached. Reproduced end to end before the fix: the question
#
#   "rank hydroxychloroquine, dexamethasone, remdesivir, nitazoxanide,
#    camostat mesylate, ivermectin, favipiravir for Zika"
#
# resolved six compounds, dropped ivermectin AND Zika, fell back to the
# default virus and REFUSED -- "This release cannot answer that about
# NCBITaxon:2697049. It holds 0 antiviral measurement(s) for it" -- a
# refusal about a virus the user never named, quoted as a raw CURIE,
# with figures. The same question naming five compounds answers Zika
# correctly.
ZIKA = "NCBITaxon:64320"
MANY_DRUGS = ["hydroxychloroquine", "dexamethasone", "remdesivir",
              "nitazoxanide", "camostat mesylate", "ivermectin",
              "favipiravir"]


@pytest.fixture
def crowded(tmp_path):
    """One virus and more compounds than MAX_RESOLVED admits."""
    def e(s, p, o, **q):
        return {"subject": s, "predicate": p, "object": o, "qualifiers": q,
                "primary_knowledge_source": "infores:chembl",
                "evidence_tier": 1, "first_asserted_date": "2021-01-01"}
    nodes = [{"id": ZIKA, "class": "OrganismTaxon",
              "properties": {"label": "Zika", "family": "Flaviviridae",
                             "baltimore_class": "IV", "is_enveloped": True}},
             {"id": "UniProtKB:NS5", "class": "Protein",
              "properties": {"taxon_id": ZIKA, "is_viral": True,
                             "sequence_hash": "h", "reviewed": True,
                             "gene_symbol": "NS5"}},
             {"id": "KGAV:G", "class": "Gene", "properties": {"symbol": "poly"}}]
    edges = [e("UniProtKB:NS5", "ENCODED_BY", "KGAV:G"),
             e("KGAV:G", "BELONGS_TO", ZIKA)]
    for i, drug in enumerate(MANY_DRUGS):
        nodes.append({"id": f"CHEMBL:{i}", "class": "SmallMolecule",
                      "properties": {"label": drug}})
        edges.append(e(f"CHEMBL:{i}", "INHIBITS", "UniProtKB:NS5",
                       ic50_nm=100 + i))
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(x) for x in edges))
    return Orchestrator(GraphTools(tmp_path))


CROWDED_Q = "rank " + ", ".join(MANY_DRUGS) + " for Zika"


def test_the_named_virus_survives_a_crowded_question(crowded):
    """The virus decides the calibration, the routes, every caveat and
    the refusal itself. It is never the right thing to drop."""
    b = crowded.build_brief(CROWDED_Q)
    assert ZIKA in b.resolved.values(), b.resolved
    assert b.allowed, b.refusal
    assert "virus=" + ZIKA in " ".join(b.trace), b.trace
    assert "(assumed)" not in " ".join(b.trace), b.trace


def test_the_answer_is_about_the_virus_that_was_named(crowded):
    b = crowded.build_brief(CROWDED_Q)
    assert "Zika" in " ".join(str(w) for w in b.warnings)
    assert "NCBITaxon:" not in (b.refusal or ""), \
        "a raw CURIE reached the user"


def test_omitted_entities_are_named_in_a_caveat(crowded):
    """A ranking of six of seven reads as a ranking of seven: "not in
    the list" is indistinguishable from "ranked last"."""
    b = crowded.build_brief(CROWDED_Q)
    said = [str(w) for w in b.warnings if "NAMES OMITTED" in str(w)]
    assert said, [str(w)[:40] for w in b.warnings]
    assert "ivermectin" in said[0] and "favipiravir" in said[0], said[0]


def test_the_entities_kept_are_the_ones_asked_for_first(crowded):
    """FIRST-NAMED, FIRST-KEPT. Matching is ordered by descending name
    length, which is right for resolving overlapping names and arbitrary
    for deciding whose compound is dropped -- it kept "camostat
    mesylate" over "remdesivir" for no reason a user could predict."""
    b = crowded.build_brief(CROWDED_Q)
    kept = [n for n in b.resolved if n != "zika"]
    assert kept == MANY_DRUGS[:len(kept)], kept


def test_a_question_under_the_cap_omits_nothing(crowded):
    b = crowded.build_brief("rank " + ", ".join(MANY_DRUGS[:4]) + " for Zika")
    assert not any("NAMES OMITTED" in str(w) for w in b.warnings)
    assert not any(t.startswith("omitted=") for t in b.trace), b.trace


def test_the_cap_is_not_a_literal_in_the_loop():
    """The threshold's owner is the module constant, and the length
    filter's owner is virus_register -- this was the THIRD copy of
    `len(name) < 4`."""
    import inspect

    from kgav.agent.orchestrator import Orchestrator as O
    src = inspect.getsource(O._resolve_mentions)
    assert "MAX_RESOLVED" in src and "MIN_SEARCH_NAME" in src
    assert "len(found) >= 6" not in src
    assert "len(name) < 4" not in src


# ------------------------- an assumed virus is said out loud ON A REFUSAL
def test_a_refusal_says_the_virus_was_assumed(crowded):
    """THE PATH WHERE IT MATTERS MOST. The condition was `named is None
    and verdict.allowed`, so the one case where the substitution does
    the most damage -- refusing, in detail and with figures, about a
    virus the user never named -- was the case that said nothing."""
    # The default must be a virus this graph HOLDS, or the absent-default
    # guard answers first with a clarification -- which is the correct
    # behaviour there and not the path under test. The first version of
    # this test reused `crowded` as built, whose default (SARS-CoV-2) is
    # not in it, so it asserted a refusal against a fixture that cannot
    # produce one.
    orch = Orchestrator(crowded.tools, default_virus=ZIKA)
    b = orch.build_brief("which compounds should be considered")
    assert not b.allowed, (b.clarification, b.trace)
    assert "VIRUS ASSUMED" in (b.refusal or ""), b.refusal[:200]
    assert any("VIRUS ASSUMED" in str(w) for w in b.warnings)
    # and the refusal itself still follows it
    assert "prospective candidate generation" in b.refusal


def test_the_assumed_caveat_claims_nothing_about_evaluation(crowded):
    """It read "{label}, which is the only virus any protocol in this
    project evaluated" -- true of SARS-CoV-2 and interpolated with
    whatever default a caller passes. The evaluation status is already
    reported, measured, by CALIBRATION / VIRUS NOT VALIDATED in the same
    answer."""
    from kgav.agent.guards import virus_assumed_caveat
    said = virus_assumed_caveat("Zika", ["DENV-2"])
    assert "only virus" not in said, said
    assert "evaluated" not in said, said
    assert "default" in said


def test_a_default_virus_absent_from_the_graph_asks_rather_than_refuses(
        tmp_path):
    """guards.check would refuse about it and state figures -- "it holds
    0 antiviral measurement(s)" -- naming it through _label, which
    returns the raw CURIE for a node that is not there."""
    nodes = [{"id": ZIKA, "class": "OrganismTaxon",
              "properties": {"label": "Zika", "family": "Flaviviridae",
                             "baltimore_class": "IV", "is_enveloped": True}}]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("")
    # default is SARS-CoV-2, which this release does not hold
    orch = Orchestrator(GraphTools(tmp_path))
    b = orch.build_brief("rank these compounds")
    assert b.clarification, (b.allowed, b.refusal)
    assert "Zika" in b.clarification
    assert "NCBITaxon:" not in b.clarification, b.clarification
    assert "default_virus_absent" in b.trace
