# End-to-end analysis pipeline

Several early scripts use fixed filenames and therefore should be executed from the project working directory.

## 1. Recreate the PubMed XML locally

The raw PubMed XML is not distributed. Follow [`PUBMED_RETRIEVAL.md`](PUBMED_RETRIEVAL.md), or run:

```bash
bash scripts/download_pubmed_edirect.sh
```

Historical retrieval date: **2026-08-27**  
Historical parsed XML dataset: **55,446 unique PMIDs**

Exact historical PMID set: `data/manifests/retrieved_pmids_2026-08-27.txt`

Because PubMed is a live database, a later reproduction run can differ slightly after retrospective indexing or record corrections. Record the rerun date and observed count.

## 2. Parse PubMed XML

```bash
python scripts/pubmed_parser.py \
  --input /path/to/pubmed_export.xml.gz \
  --parquet glioblastoma_pubmed.parquet \
  --csv glioblastoma_pubmed.csv
```

Historical expected output: **55,446 unique PMIDs**.

## 3. Initial QC

```bash
python scripts/dataset_qc.py
```

## 4. Prepare analytical corpus

```bash
python scripts/prepare_corpus.py
```

Expected output before the aggregate-record exclusion: **54,051 publications**.

## 5. Token-length analysis

```bash
python scripts/token_length_analysis.py
```

## 6. Generate biomedical embeddings

```bash
python scripts/generate_embeddings.py \
  --input glioblastoma_analysis_corpus_with_tokens.parquet \
  --output-dir embeddings
```

The embedding script excludes PMID **27454254** and creates the final **54,050 × 768** L2-normalized embedding matrix. The exact final PMID set is preserved at `data/manifests/final_corpus_pmids.txt`.

## 7. BERTopic baseline

```bash
python scripts/run_bertopic_baseline.py \
  --corpus glioblastoma_analysis_corpus_with_tokens.parquet \
  --embeddings embeddings/biomedbert_embeddings.npy \
  --mapping embeddings/embedding_pmids.parquet \
  --output-dir bertopic_baseline
```

Baseline HDBSCAN: `min_cluster_size=100`, `min_samples=10`.

## 8. HDBSCAN sensitivity analysis

```bash
python scripts/tune_hdbscan.py \
  --umap bertopic_baseline/umap_5d.npy \
  --mapping embeddings/embedding_pmids.parquet \
  --output-dir hdbscan_tuning
```

Selected candidate: Run B (`min_cluster_size=100`, `min_samples=2`).

## 9. Final BERTopic candidate B

```bash
python scripts/run_bertopic_candidate_b.py \
  --corpus glioblastoma_analysis_corpus_with_tokens.parquet \
  --embeddings embeddings/biomedbert_embeddings.npy \
  --mapping embeddings/embedding_pmids.parquet \
  --reference-clusters hdbscan_tuning/runs/B/assignments.parquet \
  --output-dir bertopic_candidate_b
```

Expected result:

- 110 non-outlier topics
- 18,862 outliers (34.90%)
- 100% raw-label and outlier-status agreement with tuning Run B

## 10. Topic-annotation workbook

Use the corrected workbook builder:

```bash
python scripts/build_topic_validation_fixed.py \
  --input-dir bertopic_candidate_b \
  --output-dir topic_validation \
  --overwrite
```

The final curated workbook distributed in this repository was then reviewed by the domain expert and is the source of truth for all downstream expert labels.

## 11. Curated figures, primary temporal analysis, and sentinel topics

```bash
python scripts/make_final_figures_curated.py \
  --validation-workbook data/validation/topic_validation_workbook_2026-09-09.xlsx \
  --embeddings embeddings/biomedbert_embeddings.npy \
  --assignments bertopic_candidate_b/topic_assignments.parquet \
  --prevalence bertopic_candidate_b/topic_prevalence_by_year.csv \
  --outliers bertopic_candidate_b/outliers_by_year.csv \
  --output-dir final_figures_curated
```

The main four-axis analysis uses **expert primary labels only**. The primary interval is 2011–2025. Raw annual prevalence is used for inference; a centered 3-year rolling mean is used only for visualization.

## 12. Axis-definition sensitivity

```bash
python scripts/axis_definition_sensitivity.py \
  --validation-workbook data/validation/topic_validation_workbook_2026-09-09.xlsx \
  --prevalence bertopic_candidate_b/topic_prevalence_by_year.csv \
  --output-dir axis_definition_sensitivity
```

This compares the primary-label-only axis definition with the broader primary-plus-secondary definition.

## 13. Clustered-only denominator sensitivity

```bash
python scripts/clustered_denominator_sensitivity.py \
  --membership axis_definition_sensitivity/primary_axis_membership.csv \
  --prevalence bertopic_candidate_b/topic_prevalence_by_year.csv \
  --output-dir clustered_denominator_sensitivity
```

All four primary axes preserve trend direction and significance/non-significance status under clustered-only normalization.

## 14. Repository verification

```bash
make check
```
