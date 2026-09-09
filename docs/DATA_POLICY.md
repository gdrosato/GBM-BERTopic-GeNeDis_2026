# Data distribution policy

## Raw PubMed data

The raw PubMed XML/XML.GZ export is **not distributed in this repository**. Parsed title–abstract corpora are also excluded.

The analysis can be reconstructed by obtaining a new PubMed export using NCBI EDirect and the exact query in `config/pubmed_query.txt`. Step-by-step instructions are provided in `docs/PUBMED_RETRIEVAL.md`. The XML dataset used for the historical analysis contained **55,446 unique PMIDs** after retrieval on **2026-08-27**. Because PubMed is updated retrospectively, a later rerun may return a slightly different count.

The reason for excluding the raw export is to avoid redistributing article abstracts for which NLM does not provide blanket third-party copyright permission.

## Not distributed

The following are excluded from Git:

- raw PubMed XML/XML.GZ;
- parsed PubMed title–abstract CSV/Parquet files;
- intermediate analytical corpora;
- full embedding matrices;
- BERTopic/HDBSCAN/UMAP serialized model objects;
- manuscript/LaTeX/BibTeX/PDF files.

## Distributed

The repository includes:

- exact retrieval query and parameters;
- analysis scripts;
- exact environment snapshot;
- curated expert-annotation workbook;
- aggregate validation and temporal outputs;
- final aggregate figures.
