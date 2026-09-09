# End-to-end analysis pipeline

Several early scripts use fixed filenames and therefore should be executed from the project working directory.

## 1. Recreate the PubMed XML locally

The raw PubMed XML is not distributed. Follow the complete EDirect instructions in [`PUBMED_RETRIEVAL.md`](PUBMED_RETRIEVAL.md), or run:

```bash
bash scripts/download_pubmed_edirect.sh
```

The helper reads the exact query from `config/pubmed_query.txt`, downloads compressed XML, and reports total and unique PMID counts.

Historical study retrieval date: **2026-08-27**  
Historical parsed XML dataset: **55,446 unique PMIDs**

Because PubMed is a live database, a later reproduction run can differ slightly after retrospective indexing or record corrections. Record the rerun date and observed count.

## 2. Parse PubMed XML

```bash
python scripts/pubmed_parser.py \
  --input /path/to/pubmed_export.xml.gz \
  --parquet glioblastoma_pubmed.parquet \
  --csv glioblastoma_pubmed.csv
```

Expected parsed records: **55,446** unique PMIDs.

## 3. Initial QC

`dataset_qc.py` expects `glioblastoma_pubmed.parquet` in the current working directory:

```bash
python scripts/dataset_qc.py
```

It generates annual counts, publication-type counts, duplicate-title inspection, and short-abstract QC files.

## 4. Prepare analytical corpus

`prepare_corpus.py` also expects `glioblastoma_pubmed.parquet` in the current working directory:

```bash
python scripts/prepare_corpus.py
```

It applies the publication-type exclusions in `config/publication_type_exclusions.txt` and writes:

- `glioblastoma_analysis_corpus.parquet`
- `glioblastoma_analysis_corpus.csv`
- QC/inspection CSVs

Expected output: **54,051 publications** before the later manual aggregate-record exclusion.

## 5. Token-length analysis

```bash
python scripts/token_length_analysis.py
```

This reads `glioblastoma_analysis_corpus.parquet` and writes:

- `documents_over_512_tokens.csv`
- `glioblastoma_analysis_corpus_with_tokens.parquet`

## 6. Generate biomedical embeddings

```bash
python scripts/generate_embeddings.py \
  --input glioblastoma_analysis_corpus_with_tokens.parquet \
  --output-dir embeddings
```

The embedding script excludes PMID **27454254**, an aggregate conference abstract collection, and creates the final **54,050 × 768** L2-normalized embedding matrix.

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

## 10. Topic-validation workbook

Both the initial and corrected validation-workbook implementations are retained. The initial `build_topic_validation.py` could create duplicate `topic` columns. The successful analysis used:

```bash
python scripts/build_topic_validation_fixed.py \
  --input-dir bertopic_candidate_b \
  --output-dir topic_validation \
  --overwrite
```

The final curated expert workbook distributed in this repository was subsequently reviewed and corrected by the domain expert.

## 11. Curated final figures and temporal analysis

```bash
python scripts/make_final_figures_curated.py \
  --validation-workbook data/validation/topic_validation_workbook_2026-09-09.xlsx \
  --embeddings embeddings/biomedbert_embeddings.npy \
  --assignments bertopic_candidate_b/topic_assignments.parquet \
  --prevalence bertopic_candidate_b/topic_prevalence_by_year.csv \
  --outliers bertopic_candidate_b/outliers_by_year.csv \
  --output-dir final_figures_curated
```

Primary temporal interval: **2011–2025**. The 3-year rolling mean is used for visualization only. Spearman tests are performed on unsmoothed annual prevalence, with Benjamini–Hochberg correction within the four axis tests and separately within the four sentinel-topic tests.
