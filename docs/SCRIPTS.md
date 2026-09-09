# Script provenance and roles

The repository contains the Python scripts used across the analysis workflow. The initial and corrected validation-workbook implementations are both retained for transparency.

| Script | Role | Used in final workflow? |
|---|---|---|
| `pubmed_parser.py` | Parse PubMed XML/XML.GZ | Yes |
| `dataset_qc.py` | Initial corpus QC | Yes |
| `prepare_corpus.py` | Publication-type filtering and analytical corpus | Yes |
| `token_length_analysis.py` | BiomedBERT token-length QC | Yes |
| `generate_embeddings.py` | Biomedical embeddings + long-document chunking | Yes |
| `run_bertopic_baseline.py` | Baseline UMAP/HDBSCAN/BERTopic | Yes |
| `tune_hdbscan.py` | HDBSCAN candidate comparison on fixed UMAP | Yes |
| `run_bertopic_candidate_b.py` | Final selected BERTopic model | Yes |
| `build_topic_validation.py` | First validation-workbook implementation | Historical; duplicate-`topic` issue |
| `build_topic_validation_fixed.py` | Corrected validation-workbook implementation | **Yes** |
| `make_final_figures_curated.py` | Curated summaries, primary-label temporal analysis, main figures, sentinel analysis | **Yes** |
| `axis_definition_sensitivity.py` | Compare primary-label axes with broader primary + secondary axes | **Yes** |
| `clustered_denominator_sensitivity.py` | Recalculate primary-axis prevalence using clustered publications only | **Yes** |

The final curation workbook is the source of truth for downstream topic labels. `build_topic_validation_fixed.py` should be used for workbook reconstruction.
