# Exact corpus PMID manifests

These files identify the exact historical PubMed records used in the study without redistributing article abstracts or raw PubMed XML.

- `retrieved_pmids_2026-08-27.txt` — **55,446 unique PMIDs** retrieved and parsed from PubMed on 2026-08-27 using the query in `../../config/pubmed_query.txt`.
- `final_corpus_pmids.txt` — **54,050 unique PMIDs** used for embedding generation and BERTopic modeling after publication-type filtering and removal of the aggregate conference record (PMID 27454254).

The final modeling corpus is a strict subset of the historical retrieval. The 1,396-record difference comprises the 1,395 publication-type exclusions plus PMID 27454254.

SHA-256 fingerprints:

- `retrieved_pmids_2026-08-27.txt`: `f2f8b539be9d12770ccb4ed165a6c51fdc01d90b1570372944c8d703594595af`
- `final_corpus_pmids.txt`: `63a5e13b776d3feec910f53ba8ada171495fc25db887087386726b28bca9648d`

Because PubMed is a live database, a later rerun of the same query may not return an identical record set. These manifests preserve the exact PMID identities used in the historical analysis.
