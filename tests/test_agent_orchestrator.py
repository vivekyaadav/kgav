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
                        "salt_collapsed": False, "stereo_collapsed": False}},
        {"id": CYTO, "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "CYTOX", "label": "chloroquine",
                        "chembl_id": "CHEMBL2", "is_approved": True,
                        "salt_collapsed": False, "stereo_collapsed": False}},
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
def _ppi_tools(tmp_path, pairs):
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
    return GraphTools(tmp_path)


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


def test_only_the_host_viral_edge_of_a_mixed_release_is_counted(tmp_path):
    """All three shapes in one release: the count is 1, not 3."""
    routes = _ppi_tools(tmp_path, [
        ("UniProtKB:NSP5", "UniProtKB:SIGMAR1"),      # host-viral
        ("UniProtKB:NSP5", "UniProtKB:NSP12"),        # viral-viral
        ("UniProtKB:SIGMAR1", "UniProtKB:TMEM41B"),   # host-host
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
                        "stereo_collapsed": False}},
    ]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("")
    orch = Orchestrator(GraphTools(tmp_path))
    resolved, _problem = orch._resolve_mentions(
        "which pathways does SARS-CoV-2 depend on?")
    assert set(resolved.values()) == {V}, (
        f"the question names one entity and resolved {resolved}")
    # And the compound IS resolved when it is actually named.
    named, _ = orch._resolve_mentions("what is depen")
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
    resolved, _ = orch._resolve_mentions("the ace2 receptor in ace2x cells")
    assert set(resolved.values()) == {"UniProtKB:ACE2R"}, (
        f"ace2 was resolved out of 'ace2x' after masking: {resolved}")
