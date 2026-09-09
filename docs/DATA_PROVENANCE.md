# Data provenance

## Corpus

- Source: PubMed
- Retrieval date: 2026-08-27
- Query: `config/pubmed_query.txt`
- Records returned/parsed in the historical study dataset: 55,446 unique PMIDs
- After publication-type exclusions: 54,051
- Aggregate-record exclusion: PMID 27454254
- Final modeling corpus: 54,050 publications

## Exact PMID manifests

- `data/manifests/retrieved_pmids_2026-08-27.txt`: 55,446 unique PMIDs from the historical PubMed retrieval.
- `data/manifests/final_corpus_pmids.txt`: 54,050 unique PMIDs used for embedding and topic modeling.

The final manifest is a strict subset of the retrieved manifest. The difference is 1,396 PMIDs: 1,395 publication-type exclusions plus the aggregate conference record PMID 27454254. The manifests contain identifiers only and provide an exact historical corpus fingerprint despite subsequent changes to PubMed.

SHA-256:

- retrieved manifest: `f2f8b539be9d12770ccb4ed165a6c51fdc01d90b1570372944c8d703594595af`
- final-corpus manifest: `63a5e13b776d3feec910f53ba8ada171495fc25db887087386726b28bca9648d`

## Unit of analysis

One PubMed publication represented by title plus abstract.

## Expert curation

The final curated expert-annotation workbook is preserved at:

`data/validation/topic_validation_workbook_2026-09-09.xlsx`

This workbook is the source of truth for expert inclusion status, primary labels, secondary labels, coherence ratings, relevance ratings, and biological rationale.

## Final temporal analysis

The primary four-axis analysis uses expert primary labels only. Broader primary-plus-secondary topic mapping is retained as a sensitivity analysis rather than as the main biological-axis definition.

The primary annual denominator is all eligible GBM publications in each year, including HDBSCAN outliers. A clustered-only denominator sensitivity analysis is also distributed and preserves the inferential interpretation of all four axes.

## Derived outputs

`data/derived/` contains aggregate curation, primary temporal, sensitivity, sentinel-topic, and hypoxia-decomposition outputs. These allow transparent verification of the reported aggregate statistics without redistributing the abstract-bearing corpus.
