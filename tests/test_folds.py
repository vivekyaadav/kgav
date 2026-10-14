"""Gate for phase 2: the two-mode comparison, and the rule that a shared fold
does not license transfer across a different nucleophile.
"""

from dataclasses import dataclass
from dataclasses import field as dc_field

from kgav.emit import Emit
from kgav.folds import NUCLEOPHILES, ChainSpec, ResolvedChain, emit_classes, gate, motif_c_present


def _chain(label, family, role, cat_type, nucleophile, classes):
    spec = ChainSpec(label=label, accession="X", patterns=("x",), family=family,
                     role=role, virus_taxon="NCBITaxon:1")
    # chain_pro_id is part of being resolved: without it .ok is False and a
    # fixture would silently stop being usable by anything that filters on
    # it. A made-up but well-formed PRO_ id keeps these realistic.
    return ResolvedChain(spec=spec, chain_name=label, start=1, end=100,
                         chain_pro_id=f"PRO_{abs(hash(label)) % 10**7:07d}",
                         catalytic_type=cat_type, nucleophile=nucleophile,
                         classes=[{"accession": a, "label": n, "db": "interpro",
                                   "level": "homologous_superfamily",
                                   "start": 1, "end": 100}
                                  for a, n in classes])


# Measured in phase 0: all three families share SSF50494 and G3DSA:2.40.10.10.
SHARED_PROTEASE_FOLD = [("SSF50494", "Trypsin-like serine protease"),
                        ("G3DSA:2.40.10.10", "Trypsin-like serine protease")]
SHARED_RDRP = [("IPR043502", "DNA/RNA polymerase superfamily"),
               ("IPR007094", "RNA-directed RNA polymerase, catalytic domain")]

COV_MPRO = _chain("SARS-CoV-2 Mpro", "Coronaviridae", "protease",
                  "cysteine", "C3408", SHARED_PROTEASE_FOLD)
PICO_3C = _chain("Poliovirus 3C", "Picornaviridae", "protease",
                 "cysteine", "C1712", SHARED_PROTEASE_FOLD)
HCV_NS3 = _chain("HCV NS3", "Flaviviridae", "protease",
                 "serine", "S1165", SHARED_PROTEASE_FOLD)
COV_RDRP = _chain("SARS-CoV-2 nsp12", "Coronaviridae", "polymerase",
                  "none", "", SHARED_RDRP)
PICO_3D = _chain("Poliovirus 3D", "Picornaviridae", "polymerase",
                 "none", "", SHARED_RDRP)


# ------------------------------------------------------- nucleophile mode
def test_same_nucleophile_and_shared_fold_is_allowed():
    """Coronavirus and picornavirus proteases: cysteine both, shared fold.
    This is the transfer the conservation hypothesis predicts."""
    g = gate(COV_MPRO, PICO_3C)
    assert g.allowed and g.mode == "nucleophile" and g.strength == "strong"
    assert "SSF50494" in g.shared


def test_different_nucleophile_is_refused_despite_identical_fold():
    """THE central check. Phase 0 measured that SARS-CoV-2 Mpro and HCV NS3
    share SSF50494 and G3DSA:2.40.10.10 -- by fold they are indistinguishable.
    Without this gate M8 recommends nirmatrelvir for hepatitis C."""
    g = gate(COV_MPRO, HCV_NS3)
    assert not g.allowed
    assert g.catalytic_type_match is False
    assert "cysteine" in g.reason and "serine" in g.reason
    # the shared fold is still reported, so the refusal is explicable
    assert "SSF50494" in g.shared


def test_same_nucleophile_without_shared_fold_is_refused():
    lonely = _chain("other", "Other", "protease", "cysteine", "C10",
                    [("IPR999999", "unrelated")])
    g = gate(COV_MPRO, lonely)
    assert not g.allowed and g.strength == "none"


# ------------------------------------------------------------- motif mode
def test_polymerases_are_compared_in_motif_mode_and_marked_weak():
    """Polymerases have no catalytic nucleophile, and motif C (GDD/SDD) is
    present in all three families -- it does not discriminate the way Cys vs
    Ser does. The result must say so."""
    g = gate(COV_RDRP, PICO_3D)
    assert g.allowed and g.mode == "motif"
    assert g.strength == "weak"
    assert "WEAK" in g.reason


def test_motif_mode_reports_catalytic_domain_entries_when_present():
    g = gate(COV_RDRP, PICO_3D)
    assert "IPR007094" in g.shared


# ------------------------------------------------------------- mixed mode
def test_protease_versus_polymerase_is_refused_outright():
    """A protease and a polymerase are not comparable targets; any fold
    similarity between them is an artifact."""
    g = gate(COV_MPRO, COV_RDRP)
    assert not g.allowed and g.mode == "refused"
    assert "not comparable" in g.reason


def test_refusal_is_symmetric():
    assert gate(COV_MPRO, COV_RDRP).allowed is gate(COV_RDRP, COV_MPRO).allowed
    assert gate(COV_MPRO, HCV_NS3).allowed is gate(HCV_NS3, COV_MPRO).allowed


# ---------------------------------------------------------------- motif C
def test_motif_c_is_a_check_not_a_detector():
    """The pattern also matches CDD in the coronavirus NiRAN domain, which is
    not motif C. It may only be used on a region already known to be a
    polymerase."""
    seq = "AAAAGDDAAAA"
    assert motif_c_present(seq, 1, len(seq)) == ["GDD@5"]
    # SDD, the coronavirus variant, is matched
    assert motif_c_present("AAAASDDAAAA", 1, 11) == ["SDD@5"]
    # CDD is NOT matched by the check pattern
    assert motif_c_present("AAAACDDAAAA", 1, 11) == []


# ------------------------------------------------------------ definitions
def test_nucleophile_table_maps_to_schema_enum_values():
    from kgav.schema import load_schema
    allowed = set(load_schema().enums["catalytic_type"])
    assert set(NUCLEOPHILES.values()) <= allowed


def test_every_phase0_protease_pair_behaves_as_measured():
    """The full phase 0 matrix, encoded so a regression is visible."""
    cys = [COV_MPRO, PICO_3C]
    ser = [HCV_NS3]
    for a in cys:
        for b in cys:
            assert gate(a, b).allowed, f"{a.chain_name} vs {b.chain_name}"
    for a in cys:
        for b in ser:
            assert not gate(a, b).allowed, f"{a.chain_name} vs {b.chain_name}"


# ------------------------------- a class whose members disagree is not a gate
@dataclass
class _Chain:
    node_id: str
    catalytic_type: str
    ok: bool = True
    classes: list = dc_field(default_factory=list)


_SHARED = [{"accession": "IPR009003", "label": "Peptidase S1, PA clan",
            "level": "homologous_superfamily", "start": 1, "end": 200}]


def _emit(order):
    em = Emit()
    stats = emit_classes(em, order, "infores:interpro")
    return em.nodes["INTERPRO:IPR009003"]["properties"], stats


def test_a_mixed_nucleophile_class_is_marked_unknown_not_guessed():
    """Emit.node merges with properties.update -- last write wins, silently,
    recording only a node_merged counter. So a class holding a cysteine
    protease and a serine protease took whichever chain was written last,
    and reversing the order flipped it from serine to cysteine.

    That is the field M8's gate reads. The schema warns that without the
    gate M8 "traverses from SARS-CoV-2 Mpro to HCV NS3 -- same superfamily,
    opposite nucleophile -- and recommends nirmatrelvir for hepatitis C",
    and those two land in the same InterPro superfamily. A gate whose value
    depends on emission order is not a gate.
    """
    mpro = _Chain("UniProtKB:MPRO", "cysteine", classes=list(_SHARED))
    ns3 = _Chain("UniProtKB:NS3", "serine", classes=list(_SHARED))
    fwd, s1 = _emit([mpro, ns3])
    rev, s2 = _emit([ns3, mpro])
    assert fwd["catalytic_type"] == "unknown"
    assert rev["catalytic_type"] == "unknown"
    assert s1["classes_mixed_catalytic_type"] == 1
    assert s2["classes_mixed_catalytic_type"] == 1


def test_emission_order_cannot_change_a_class_catalytic_type():
    """The property the previous code lacked, stated directly."""
    a = _Chain("UniProtKB:A", "cysteine", classes=list(_SHARED))
    b = _Chain("UniProtKB:B", "serine", classes=list(_SHARED))
    c = _Chain("UniProtKB:C", "metallo", classes=list(_SHARED))
    assert _emit([a, b, c])[0] == _emit([c, b, a])[0]


def test_an_agreeing_class_keeps_its_catalytic_type():
    """So "unknown" marks disagreement rather than being a blanket refusal --
    a gate that never opens is as useless as one that never closes."""
    props, stats = _emit([_Chain("UniProtKB:A", "cysteine", classes=list(_SHARED)),
                          _Chain("UniProtKB:B", "cysteine", classes=list(_SHARED))])
    assert props["catalytic_type"] == "cysteine"
    assert stats["classes_mixed_catalytic_type"] == 0


def test_unknown_is_a_declared_value_of_the_enum():
    """Marking a class unknown must be sayable in the schema, not a string
    the validator will reject the moment the folds layer is assembled."""
    from pathlib import Path

    import yaml
    root = Path(__file__).resolve().parents[1]
    d = yaml.safe_load((root / "schema/kg_schema.yaml").read_text())
    assert "unknown" in d["enums"]["catalytic_type"]


# ------------------------------------------------- the id that joins the spine
def test_node_id_is_the_spine_s_pro_id_not_a_coordinate_key():
    """THE REASON THIS LAYER NEVER REACHED THE GRAPH. node_id was
    UniProtKB:{accession}_{start}_{end} while ingest_spine keys the same
    chain on its PRO_ id. Two schemes for one protein, so every
    MEMBER_OF_CLASS edge pointed at a node that does not exist -- and
    assemble_release returns 1 on a dangling edge, which is why `folds` is
    absent from LAYER_PRECEDENCE.
    """
    spec = ChainSpec(label="Mpro", accession="P0DTD1", patterns=("3C-like",),
                     family="Coronaviridae", role="protease",
                     virus_taxon="NCBITaxon:2697049")
    r = ResolvedChain(spec=spec, chain_name="3C-like proteinase",
                      start=3264, end=3569, chain_pro_id="PRO_0000449623")
    assert r.ok
    assert r.node_id == "UniProtKB:PRO_0000449623"
    assert "3264" not in r.node_id, "coordinates must not appear in the key"


def test_a_chain_without_a_pro_id_is_not_resolved():
    """Same rule the spine applies: no PRO_ id, no stable node key. Being
    part of `ok` matters because `ok` is what every caller filters on -- a
    chain that is merely 'not emitted later' would still be compared,
    counted and reported as usable."""
    spec = ChainSpec(label="X", accession="Q1", patterns=("x",),
                     family="Flaviviridae", role="protease",
                     virus_taxon="NCBITaxon:64320")
    r = ResolvedChain(spec=spec, chain_name="x", start=1, end=50)
    assert not r.ok


def test_reading_node_id_without_a_pro_id_raises_rather_than_inventing_one():
    """Returning a coordinate key here is what the previous version did, and
    it failed silently all the way through to a dangling edge."""
    import pytest as _pytest

    spec = ChainSpec(label="X", accession="Q1", patterns=("x",),
                     family="Flaviviridae", role="protease",
                     virus_taxon="NCBITaxon:64320")
    r = ResolvedChain(spec=spec, chain_name="x", start=1, end=50)
    with _pytest.raises(ValueError, match="no PRO_ id"):
        _ = r.node_id


def test_find_chain_keeps_the_feature_id_uniprot_returns():
    """Verified against rest.uniprot.org for P0DTD1 on 2026-10-15 with this
    exact field selection -- a field present in a proteome download is not
    evidence about a single-entry query, so it was checked separately."""
    from kgav import folds

    payload = {"features": [
        {"type": "Chain", "featureId": "PRO_0000449618",
         "description": "Replicase polyprotein 1ab",
         "location": {"start": {"value": 1}, "end": {"value": 7096}}},
        {"type": "Chain", "featureId": "PRO_0000449623",
         "description": "3C-like proteinase nsp5",
         "location": {"start": {"value": 3264}, "end": {"value": 3569}}},
    ]}
    real, folds._get = folds._get, lambda _u: payload
    try:
        name, s, e, note, pro = folds.find_chain("P0DTD1", ("3C-like",))
    finally:
        folds._get = real
    assert (s, e, note, pro) == (3264, 3569, "", "PRO_0000449623")
    assert "3C-like" in name


def test_a_matched_chain_with_no_pro_id_is_reported_not_keyed():
    from kgav import folds

    payload = {"features": [
        {"type": "Chain", "description": "mystery chain",
         "location": {"start": {"value": 1}, "end": {"value": 9}}}]}
    real, folds._get = folds._get, lambda _u: payload
    try:
        _name, s, e, note, pro = folds.find_chain("Q1", ("mystery",))
    finally:
        folds._get = real
    assert (s, e, pro) == (1, 9, "")
    assert "no PRO_ id" in note


def test_role_decides_mixed_mode_not_whether_a_nucleophile_resolved():
    """THE GATE HELD BY LUCK. MIXED was decided by a_nuc != b_nuc -- whether
    a nucleophile was FOUND, a property of the sequence scan. SARS-CoV-2
    PLpro resolves catalytic_type "none" because its spec matches UniProt's
    "Papain-like protease nsp3", the whole 1,945-residue multi-domain chain
    rather than the protease domain. So a PROTEASE reached MOTIF mode beside
    the polymerases, where an allowed weak transfer protease->polymerase was
    possible -- the one comparison the docstring calls an artifact.

    It never fired on v0.1 only because PLpro shares no catalytic-domain
    entry with any cross-family polymerase.
    """
    plpro = _chain("SARS-CoV-2 PLpro", "Coronaviridae", "protease",
                   "none", "", SHARED_RDRP)      # same classes as an RdRp
    g = gate(plpro, PICO_3D)
    assert not g.allowed, "a protease must not transfer to a polymerase"
    assert g.mode == "refused"
    assert "different target roles" in g.reason

    # And the old rule alone would have let it through: both lack a
    # nucleophile, and they share catalytic-domain entries.
    assert plpro.catalytic_type == PICO_3D.catalytic_type == "none"
    assert set(c["accession"] for c in plpro.classes) & \
           set(c["accession"] for c in PICO_3D.classes)


def test_matching_roles_are_still_compared():
    """The role check must not refuse everything: a gate that never opens is
    as useless as one that never closes."""
    assert gate(COV_MPRO, PICO_3C).allowed
    assert gate(COV_RDRP, PICO_3D).mode == "motif"


# ------------------------------------------- which rule read the active sites
def _machinery(sites, seq_len=6000):
    """catalytic_machinery against a stubbed UniProt response.

    `sites` is [(position, residue, description)]; the sequence is built so
    each position holds its residue.
    """
    from kgav import folds

    seq = ["A"] * seq_len
    for pos, aa, _d in sites:
        seq[pos - 1] = aa
    payload = {"sequence": {"value": "".join(seq)},
               "features": [{"location": {"start": {"value": p}},
                             "description": d} for p, _aa, d in sites]}
    real, folds._get = folds._get, lambda _u: payload
    try:
        return folds.catalytic_machinery("X", 1, seq_len)
    finally:
        folds._get = real


def test_an_annotated_nucleophile_decides_and_says_so():
    _res, nuc, typ, basis = _machinery(
        [(100, "H", "Proton acceptor; for 3CL-PRO activity"),
         (200, "C", "Nucleophile; for 3CL-PRO activity")])
    assert (nuc, typ, basis) == ("C200", "cysteine", "annotated")


def test_composition_decides_a_chymotrypsin_like_triad():
    """No annotation, but His-...-Cys with the nucleophile last is still
    readable."""
    _res, nuc, typ, basis = _machinery(
        [(100, "H", "For protease activity"),
         (200, "C", "For protease activity")])
    assert (nuc, typ, basis) == ("C200", "cysteine", "composition")


def test_a_papain_like_triad_is_read_by_composition():
    """THE PLPRO CASE. Its active sites are C1674, H1835, D1849, all
    described "For PL-PRO activity" with no nucleophile annotation. A
    papain-like protease runs Cys-His-Asp with the nucleophile FIRST, and
    the composition rule used to test the LAST residue -- so a cysteine
    protease read "none" and the gate could not evaluate it.
    """
    _res, nuc, typ, basis = _machinery(
        [(1674, "C", "For PL-PRO activity"),
         (1835, "H", "For PL-PRO activity"),
         (1849, "D", "For PL-PRO activity")])
    assert (nuc, typ, basis) == ("C1674", "cysteine", "composition")


def test_a_serine_without_a_histidine_is_still_refused():
    """THE CONDITION THAT KEEPS THE WIDENING HONEST. SARS-CoV-2 RdRp's
    sites are S5151, D5152, D5153: a serine, no histidine. Dropping the
    His test to catch PLpro would read "serine" for a polymerase -- one
    step from the "aspartic D5153" misclassification the explicit branch
    was fixed to prevent.
    """
    _res, nuc, typ, _b = _machinery(
        [(5151, "S", "For polymerase activity"),
         (5152, "D", "For polymerase activity"),
         (5153, "D", "For polymerase activity")], seq_len=6000)
    assert (nuc, typ) == ("", "none")


def test_two_nucleophile_candidates_are_refused_rather_than_picked():
    """Two C/S sites mean the set spans more than one active site and
    nothing here can say which belongs to the protease."""
    _res, nuc, typ, _b = _machinery(
        [(10, "C", "For activity A"), (20, "H", "For activity A"),
         (30, "S", "For activity B")])
    assert (nuc, typ) == ("", "none")


def test_the_widened_rule_moves_only_plpro():
    """THE BLAST RADIUS, measured before the rule changed and pinned here.

    Every chain in SPECS that has active sites, with the residues exactly
    as ingest_folds reported them on 2026-10-15. Fourteen of these were
    already decided by composition; the point of the table is that the
    position-independent rule leaves all fourteen bit-identical and moves
    exactly one chain.
    """
    cases = [
        ("SARS-CoV 3CLpro", "H3281, C3385", "cysteine"),
        ("MERS-CoV 3CLpro", "H3288, C3395", "cysteine"),
        ("HCoV-229E 3CLpro", "H3006, C3109", "cysteine"),
        ("Dengue 2 NS3pro", "H1526, D1550, S1610", "serine"),
        ("Zika NS3pro", "H1549, D1573, S1633", "serine"),
        ("West Nile NS3pro", "H1552, D1576, S1636", "serine"),
        ("HCV NS3pro", "H1083, D1107, S1165", "serine"),
        ("Yellow fever NS3pro", "H1537, D1561, S1622", "serine"),
        ("JEV NS3pro", "H1555, D1579, S1639", "serine"),
        ("EV-A71 3C", "H1588, E1619, C1695", "cysteine"),
        ("Poliovirus 3C", "H1605, E1636, C1712", "cysteine"),
        ("HAV 3C", "H1563, D1603, C1691", "cysteine"),
        ("CVB3 3C", "H1580, E1611, C1687", "cysteine"),
        ("Rhinovirus A 3C", "H1561, E1592, C1668", "cysteine"),
        # the four polymerases that must stay unreadable
        ("SARS-CoV-2 RdRp", "S5151, D5152, D5153", "none"),
        ("SARS-CoV RdRp", "S5128, D5129, D5130", "none"),
        ("Dengue 2 NS5", "K2552, D2637, K2672, E2708", "none"),
        ("Zika NS5", "K2577, D2662, K2698, E2734", "none"),
        # and the one that moves
        ("SARS-CoV-2 PLpro", "N1062, G1070, C1674, H1835, D1849", "cysteine"),
    ]
    for label, residues, expected in cases:
        sites = [(int(r[1:]), r[0], "For activity") for r in residues.split(", ")]
        _res, _nuc, typ, _b = _machinery(sites, seq_len=6000)
        assert typ == expected, f"{label}: {typ} != {expected}"


def test_no_active_sites_at_all_is_distinguishable():
    assert _machinery([])[3] == "none"
