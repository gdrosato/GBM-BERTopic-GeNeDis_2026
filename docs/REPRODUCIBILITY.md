# Reproducibility notes

## Exact environment

The exact analysis environment snapshot is provided in `environment/`.

- Python 3.12.3
- Linux 6.8.0-138-generic x86_64, glibc 2.39
- PyTorch 2.13.0+cu130
- CUDA build 13.0
- NVIDIA A2 GPU

Use `environment/pip_freeze.txt` for the exact installed package versions.

## Randomness

UMAP used random state 42. The selected HDBSCAN candidate was rerun during the final BERTopic analysis and matched the tuning-run raw labels and outlier assignments for all 54,050 documents.

## Final temporal definitions

The main four-axis temporal analysis is intentionally defined from **expert primary labels only**. This avoids assigning every publication in a topic to an axis that appears only as a secondary expert association.

Two explicit robustness analyses are included:

1. `axis_definition_sensitivity.py` compares the primary-label definition with the broader primary + secondary definition.
2. `clustered_denominator_sensitivity.py` compares the primary all-eligible-GBM denominator with a clustered-publications-only denominator.

The clustered-only analysis preserved both trend direction and significance/non-significance status for all four primary axes.

## Large intermediates

Embeddings, parsed abstract-bearing corpora, and serialized model objects are intentionally not committed. They can be regenerated from a locally retrieved PubMed corpus by following `docs/PIPELINE.md`.

## PubMed mutability

PubMed is a live database. The exact bounded query and historical retrieval date are provided, but a later retrieval may differ because of retrospective indexing or record corrections. The historical analysis parsed 55,446 unique PMIDs. Exact PMID-only manifests are therefore provided in `data/manifests/` for both the historical retrieval (55,446 PMIDs) and the final modeling corpus (54,050 PMIDs). These manifests allow exact corpus identity to be checked even if a future PubMed rerun differs.
