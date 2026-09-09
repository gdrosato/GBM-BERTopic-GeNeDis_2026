# Data distribution policy

## Raw PubMed data

The raw PubMed XML/XML.GZ export is **not distributed in this repository**. Parsed title–abstract corpora are also excluded.

The analysis can be reconstructed by obtaining a new PubMed export using the exact query in `config/pubmed_query.txt`. The historical retrieval was performed on **2026-08-27** and returned **55,446** records.

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
