#!/usr/bin/env python3
"""Phase 1 — schema 0.8.0: TargetClass, fold edges, and metapath M8.

WHAT PHASE 0 ESTABLISHED, and why it changes this design.

The plan assumed coronavirus and flavivirus proteases have DIFFERENT folds, so
that fold similarity alone would separate them. That is wrong. Measured at chain
level against InterPro:

    Coronaviridae 3CLpro : SSF50494, G3DSA:2.40.10.10
    Picornaviridae 3C    : SSF50494, G3DSA:2.40.10.10
    Flaviviridae NS3     : SSF50494, G3DSA:2.40.10.10

All three share the trypsin-like serine protease superfamily. What separates
them is the NUCLEOPHILE, read from the sequence at the annotated active site:

    SARS-CoV-2 Mpro   H3304, C3408          cysteine
    Poliovirus 3C     H1605, E1636, C1712   cysteine
    HCV NS3           H1083, D1107, S1165   serine
    Dengue 2 NS3      H1526, D1550, S1610   serine

Every one has histidine first and the nucleophile last, in the same spatial
arrangement. That IS the shared fold; only the final residue differs.

CONSEQUENCE FOR THE SCHEMA. The catalytic-type gate is not a safety check on top
of fold similarity -- it is the entire discriminator. Without it the breadth
score would recommend cysteine-protease inhibitors for serine proteases, since
by fold they are indistinguishable.

AND THE GATE MUST READ THE RESIDUE, NOT THE FOLD NAME. InterPro names the
coronavirus protease superfamily "trypsin-like serine protease" because that is
the architecture, regardless of which nucleophile occupies it. A classifier
reading entry names returned "ambiguous" for all six cysteine proteases. The
residue is stored explicitly for that reason.

Run from the repository root. Idempotent.
"""
from __future__ import annotations

import pathlib
import re
import sys

SCHEMA = pathlib.Path("schema/kg_schema.yaml")


def patch(text: str, anchor: str, addition: str, marker: str,
          label: str) -> tuple[str, bool]:
    if marker in text:
        print(f"  {label}: already applied")
        return text, False
    if anchor not in text:
        print(f"  {label}: ANCHOR NOT FOUND — apply by hand")
        sys.exit(1)
    print(f"  {label}: patched")
    return text.replace(anchor, addition, 1), True


def main() -> int:
    if not SCHEMA.exists():
        print(f"{SCHEMA} not found — run from the repository root")
        return 1
    t = SCHEMA.read_text()

    # ---------------------------------------------------------------- enums
    t, _ = patch(
        t,
        "  host_factor_direction: [dependency, restriction]",
        """  host_factor_direction: [dependency, restriction]

  # The nucleophile that defines a protease's catalytic class. Phase 0 measured
  # this from the sequence at UniProt's annotated active site, because fold
  # assignment cannot supply it: InterPro calls the CORONAVIRUS protease
  # superfamily "trypsin-like serine protease" -- that is the architecture, not
  # the chemistry. Reading entry names returned "ambiguous" for every cysteine
  # protease tested.
  catalytic_type: [cysteine, serine, threonine, aspartic, metallo, none, unknown]

  # How a protein was assigned to a target class.
  class_assignment_method: [interpro, scop, cath, pfam, hmm, manual]

  # How two classes were compared.
  fold_comparison_method: [interpro_shared, tmalign, dali, foldseek, manual]""",
        "catalytic_type: [cysteine",
        "enums")

    # --------------------------------------------------- Protein properties
    t, _ = patch(
        t,
        """      sequence_hash:  {type: string, required: true}
      reviewed:       {type: bool,   required: true}""",
        """      sequence_hash:  {type: string, required: true}
      reviewed:       {type: bool,   required: true}
      # Catalytic machinery, recorded per protein because it is what
      # distinguishes otherwise fold-identical proteases. catalytic_residues
      # holds positions in the PARENT accession's numbering, matching how
      # UniProt reports active sites on a polyprotein.
      catalytic_type:     {type: string, enum: catalytic_type}
      catalytic_residues: {type: "list[string]", doc: "e.g. ['H3304', 'C3408']"}
      nucleophile:        {type: string, doc: "residue and position, e.g. C3408"}""",
        "catalytic_residues:",
        "Protein properties")

    # ------------------------------------------------------ TargetClass node
    t, _ = patch(
        t,
        "  Gene:\n    id_prefixes: [NCBIGene, KGAV]",
        """  # A fold or functional family -- NOT an organism. Added in 0.8.0 because
  # every metapath previously terminated at OrganismTaxon, which made
  # "this drug inhibits 3C-like proteases" inexpressible rather than merely
  # unimplemented. Cross-family reasoning requires a node for the class itself.
  TargetClass:
    id_prefixes: [SCOP, CATH, PFAM, INTERPRO, SSF, KGAV]
    properties:
      label:          {type: string, required: true}
      source:         {type: string, required: true, enum: class_assignment_method}
      level:          {type: string, doc: "superfamily | family | domain"}
      # Identity-defining: two classes differing in catalytic type are
      # different things, never two views of one. Same role that `direction`
      # plays on HOST_FACTOR_FOR.
      catalytic_type: {type: string, enum: catalytic_type, required: true}
      scop_id:        {type: string}
      cath_id:        {type: string}
      pfam_id:        {type: string}
      interpro_id:    {type: string}

  Gene:
    id_prefixes: [NCBIGene, KGAV]""",
        "TargetClass:",
        "TargetClass node")

    # ------------------------------------------------------------- new edges
    t, _ = patch(
        t,
        "  - predicate: ENCODED_BY",
        """  - predicate: MEMBER_OF_CLASS
    subject: Protein
    object: TargetClass
    qualifiers:
      method:     {type: string, enum: class_assignment_method}
      evalue:     {type: float, min: 0}
      coverage:   {type: float, min: 0, max: 1}
      identity:   {type: float, min: 0, max: 1}
      start:      {type: int, min: 1}
      end:        {type: int, min: 1}
    required_qualifiers: [method]

  # Structural similarity between two target classes. THE GATE LIVES HERE.
  #
  # catalytic_type_match is required and is a HARD gate, not a weight. Phase 0
  # showed coronavirus, picornavirus and flavivirus proteases share SSF50494
  # and G3DSA:2.40.10.10 -- by fold they are indistinguishable. Only the
  # nucleophile separates them (Cys vs Ser). A fold similarity without this
  # check is exactly how a cysteine-protease inhibitor gets recommended for a
  # serine protease.
  - predicate: FOLD_SIMILAR_TO
    subject: TargetClass
    object: TargetClass
    qualifiers:
      tm_score:             {type: float, min: 0, max: 1}
      rmsd_angstrom:        {type: float, min: 0}
      catalytic_type_match: {type: bool}
      shared_superfamilies: {type: "list[string]"}
      method:               {type: string, enum: fold_comparison_method}
    required_qualifiers: [catalytic_type_match, method]

  - predicate: ENCODED_BY""",
        "predicate: MEMBER_OF_CLASS",
        "fold edges")

    # ------------------------------------------------------- identity quals
    t, _ = patch(
        t,
        """    identity_qualifiers: [direction]""",
        """    identity_qualifiers: [direction]
    # (see also TargetClass.catalytic_type, identity-defining for the same
    #  reason: a cysteine and a serine protease class are different facts)""",
        "identity-defining for the same",
        "identity note")

    # ------------------------------------------------------------------- M8
    t, _ = patch(
        t,
        "  M7:",
        """  # THE ROUTE v2 EXISTS FOR. A drug inhibits a protein of class X; class X is
  # structurally similar to class Y; a protein of class Y belongs to this
  # virus. This is the only metapath that reaches a virus with ZERO activity
  # data, which is the population the v1 temporal split showed to be
  # unreachable (77% of post-cutoff compounds absent from the prior graph).
  #
  # The constraint is not optional. Without it M8 traverses from SARS-CoV-2
  # Mpro to HCV NS3 -- same superfamily, opposite nucleophile -- and recommends
  # nirmatrelvir for hepatitis C.
  M8:
    doc: "cross-family direct-acting via conserved fold, catalytic type gated"
    path: [SmallMolecule, INHIBITS, Protein, MEMBER_OF_CLASS, TargetClass, FOLD_SIMILAR_TO, TargetClass, MEMBER_OF_CLASS, Protein, ENCODED_BY, Gene, BELONGS_TO, OrganismTaxon]
    constraints: {FOLD_SIMILAR_TO.catalytic_type_match: true}
    exclude_hubs: true

  M7:""",
        "M8:",
        "metapath M8")

    # --------------------------------------------------------- version bump
    m = re.search(r'schema_version:\s*"([0-9.]+)"', t)
    if not m:
        print("  version: no schema_version line found — bump by hand")
    elif m.group(1) == "0.8.0":
        print("  version: already 0.8.0")
    else:
        t = t[:m.start()] + 'schema_version: "0.8.0"' + t[m.end():]
        print(f"  version: {m.group(1)} -> 0.8.0")

    SCHEMA.write_text(t)
    print("\nschema written. Verify with:")
    print("  PYTHONPATH=src python -c \"import sys;sys.path.insert(0,'src');"
          "from kgav.schema import load_schema;s=load_schema();"
          "print('v'+s.version, len(s.metapaths),'metapaths');"
          "print('M8:',' '.join(s.metapaths['M8']['path']))\"")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
