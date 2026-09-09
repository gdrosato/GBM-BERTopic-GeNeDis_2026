# Script provenance and roles

All scripts in this repository were taken directly from the server archive `server_python_scripts.zip` supplied on 2026-09-09 and are preserved unchanged.

| Script | Role | Used in final workflow? |
|---|---|---|
| `pubmed_parser.py` | Parse PubMed XML/XML.GZ | Yes |
| `dataset_qc.py` | Initial corpus QC | Yes |
| `prepare_corpus.py` | Publication-type filtering and analytical corpus | Yes |
| `token_length_analysis.py` | BiomedBERT token-length QC | Yes |
| `generate_embeddings.py` | BiomedBERT embeddings + long-document chunking | Yes |
| `run_bertopic_baseline.py` | Baseline UMAP/HDBSCAN/BERTopic | Yes |
| `tune_hdbscan.py` | HDBSCAN candidate comparison on fixed UMAP | Yes |
| `run_bertopic_candidate_b.py` | Final selected BERTopic model | Yes |
| `build_topic_validation.py` | First validation-workbook implementation | Historical; failed at duplicate `topic` column |
| `build_topic_validation_fixed.py` | Corrected validation-workbook implementation | **Yes** |
| `make_final_figures_curated.py` | Curated validation, figures, temporal statistics | **Yes** |

The two validation scripts are both preserved for provenance. Reproduction should use `build_topic_validation_fixed.py`.
