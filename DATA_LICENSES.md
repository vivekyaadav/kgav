# Data source licences

The code here is MIT. The data is not. Each source carries its own terms, they
differ materially, and one of them prohibits commercial use outright.

No raw source data is committed. data/ is gitignored; releases are reproduced
by running the ingest against locally downloaded sources, which keeps
redistribution questions out of the repository entirely.

## The DrugBank decision

DrugBank's academic download is CC BY-NC — non-commercial only. Once its
identifiers and annotations are baked into a release, removing them later means
rebuilding, not patching.

Default for this project: build on DrugCentral and ChEMBL, leave DrugBank out.
Coverage loss is modest and the graph stays usable for any downstream purpose.
The drugbank_id property exists in the schema so the field is there if a licence
is ever obtained; it stays null until then.

Share-alike is the second thing to watch. ChEMBL (CC BY-SA 3.0) and DrugCentral
(CC BY-SA 4.0) constrain how derived releases may be redistributed, not whether
they may be used.

## Attribution

Every edge carries primary_knowledge_source as an infores CURIE, so attribution
is queryable per fact rather than only stated in aggregate. Any publication or
release derived from this graph should cite the sources that actually
contributed to it — SELECT DISTINCT primary_knowledge_source over the subgraph
used, not a blanket citation of the whole register.

## Register

Fill the verified column when you download each source, from the licence
statement on the download page itself — not from this table, and not from a
paper citing it. Terms change. Record the licence string and the access date in
that release's MANIFEST.json.

| Source | Used for | Licence (to confirm) | Commercial | verified |
|---|---|---|---|---|
| UniProt | proteomes, ID mapping | CC BY 4.0 | yes | |
| NCBI Taxonomy / Gene | taxonomy, gene ids | US public domain | yes | |
| ICTV VMR | virus reconciliation | check current release | | |
| MONDO | disease ontology | CC BY 4.0 | yes | |
| ChEMBL 37 | bioactivities, indications | CC BY-SA 3.0 | share-alike | |
| DrugBank | EXCLUDED | CC BY-NC 4.0 | NO | n/a |
| DrugCentral | indications, targets | CC BY-SA 4.0 | share-alike | |
| BindingDB | binding affinities | CC BY 3.0 | yes | |
| PubChem | bioassays | US public domain | yes | |
| STRING v12.0 | host PPI | CC BY 4.0 | yes | |
| IntAct | PPI, virus-host | CC BY 4.0 | yes | |
| BioGRID ORCS | CRISPR screens | check current terms | | |
| VirHostNet | virus-host PPI | check current terms | | |
| Reactome | pathways | CC0 | yes | |
| Gene Ontology | pathways, processes | CC BY 4.0 | yes | |
| Cellosaurus | cell lines | CC BY 4.0 | yes | |
| ClinicalTrials.gov | trials, hard negatives | US public domain | yes | |
| dbSNP / ClinVar | variants | US public domain | yes | |
