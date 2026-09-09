# GBM BERTopic — GeNeDis 2026 reproducibility repository

Code, expert curation, aggregate results, and figures supporting the study:

**AI-Driven Mapping of Hypoxia, Metabolism, Immunity, and Treatment Resistance in Glioblastoma: A BERTopic-Based Analysis of the Biomedical Literature**

Authors: **Stamatia Pouliliou** and **George Drosatos**

## Scope

This repository contains:

- the Python scripts used in the analysis workflow;
- the exact PubMed query and analysis parameters;
- PMID-only manifests identifying the exact historical retrieval and final modeling corpus;
- the exact Python/package/runtime environment export;
- the curated domain-expert topic-annotation workbook;
- aggregate derived CSV/JSON/TXT outputs used for validation and temporal analysis;
- final main and supplementary figures;
- lightweight repository-integrity tests.

It intentionally **does not contain**:

- the PubMed XML/XML.GZ export or article abstracts;
- parsed title–abstract corpora or embedding matrices;
- BERTopic/UMAP/HDBSCAN model binaries.

## Key results

- PubMed XML records used in the analysis: **55,446 unique PMIDs**
- Final modeling corpus: **54,050 publications**
- Biomedical embedding dimensionality: **768**
- BERTopic topics: **110**
- Clustered publications: **35,188**
- HDBSCAN outliers: **18,862 (34.90%)**
- Algorithmically prioritized topics for expert review: **77**
- Expert-retained topics: **41**

### Primary temporal analysis, 2011–2025

Axis membership is defined from the **expert primary label only**:

| Axis | Topics | Spearman rho | BH q | Interpretation |
|---|---:|---:|---:|---|
| Hypoxia / angiogenesis | 2 | -0.814 | 0.0003 | significant decrease |
| Metabolism / stress | 3 | 0.843 | 0.0002 | significant increase |
| Immune / tumor microenvironment | 6 | 0.868 | 0.0001 | significant increase |
| Treatment resistance | 12 | -0.443 | 0.0983 | negative, non-significant |

Two robustness analyses are included:

1. **Axis-definition sensitivity:** primary labels versus primary + secondary expert labels.
2. **Denominator sensitivity:** all eligible GBM publications versus clustered publications only.

All four primary axes preserved trend direction and significance/non-significance status under the clustered-only denominator. The immune/TME and hypoxia/angiogenesis findings were also robust to the broader axis definition; metabolism was attenuated under the broader definition, and treatment-resistance direction depended on axis composition while remaining non-significant.

## Repository structure

```text
.
├── config/                 Exact PubMed query and model/analysis parameters
├── scripts/                Retrieval, analysis, curation, and sensitivity scripts
├── environment/            Exact Python/pip/runtime snapshot
├── data/
│   ├── raw/                Documentation only; no raw PubMed data committed
│   ├── manifests/          PMID-only manifests for the exact historical corpora
│   ├── validation/         Curated expert annotation workbook
│   └── derived/            Aggregate validation/temporal/sensitivity outputs
├── results/figures/        Main and supplementary figures
├── docs/                   Pipeline, provenance, codebook, data policy
├── tests/                  Lightweight integrity checks
└── .github/workflows/      GitHub Actions integrity CI
```

## Environment

The exact analysis environment snapshot is provided under `environment/`.

Core runtime:

- Python **3.12.3**
- PyTorch **2.13.0+cu130**
- CUDA build **13.0**
- GPU **NVIDIA A2**
- BERTopic **0.17.4**
- sentence-transformers **6.0.0**
- transformers **5.16.1**
- UMAP **0.5.12**
- HDBSCAN **0.8.44**

See `environment/pip_freeze.txt` for the complete package snapshot.

## Reproducing the analysis

The raw PubMed export is not redistributed. Recreate it locally with the exact study query and NCBI EDirect. Full instructions are provided in [`docs/PUBMED_RETRIEVAL.md`](docs/PUBMED_RETRIEVAL.md).

For automated retrieval using the repository query:

```bash
bash scripts/download_pubmed_edirect.sh
```

The historical XML dataset used for analysis contained **55,446 unique PMIDs**. PubMed is a live database, so later reruns may differ slightly because of retrospective indexing or record changes. To preserve the exact historical study population, this repository includes PMID-only manifests:

- [`data/manifests/retrieved_pmids_2026-08-27.txt`](data/manifests/retrieved_pmids_2026-08-27.txt) — 55,446 PMIDs in the retrieved PubMed dataset.
- [`data/manifests/final_corpus_pmids.txt`](data/manifests/final_corpus_pmids.txt) — 54,050 PMIDs used for embedding and topic modeling.

These files contain identifiers only; they do not redistribute article titles, abstracts, or raw PubMed records. See [`data/manifests/README.md`](data/manifests/README.md) for provenance and fingerprints.

Then follow [`docs/PIPELINE.md`](docs/PIPELINE.md).

## Expert curation

The final curated workbook is:

`data/validation/topic_validation_workbook_2026-09-09.xlsx`

It is the source of truth for expert inclusion status, primary labels, secondary labels, coherence ratings, relevance ratings, and biological rationale.

## Final figures

Main figures:

- `results/figures/main/figure_topic_landscape.pdf`
- `results/figures/main/figure_temporal_axes_primary.pdf`

Selected supplementary/reproducibility figures include the expanded axis-definition sensitivity analysis, sentinel-topic trajectories, validation summaries, and outlier-rate QC.

## Data-distribution policy

The raw PubMed export and parsed abstract-bearing corpora are not included. PMID-only corpus manifests are included to identify the exact historical retrieval and final modeling corpus without redistributing abstract text. See [`docs/DATA_POLICY.md`](docs/DATA_POLICY.md).

## Integrity

```bash
make check
```

or run the tests individually:

```bash
python tests/test_repo_tree.py
python tests/test_validation_workbook.py
python tests/test_temporal_outputs.py
python tests/test_corpus_manifests.py
```

## Citation

If you use this repo in a scientific publication, we would appreciate using the following citation:

- Pouliliou, S. and Drosatos, G. (2026). AI-Driven Mapping of Hypoxia, Metabolism, Immunity, and Treatment Resistance in Glioblastoma: A BERTopic-Based Analysis of the Biomedical Literature. In GeNeDIS 2026, pages 1-15, AEMB Vol. xxxx, Springer.

and as BibTeX:

```
@InProceedings{Pouliliou_GBM_2026,
    author       = {Pouliliou, Stamatia and Drosatos, George},
    title        = {AI-Driven Mapping of Hypoxia, Metabolism, Immunity, and Treatment Resistance in Glioblastoma: A BERTopic-Based Analysis of the Biomedical Literature},
    keywords     = {Glioblastoma; BERTopic; Biomedical Literature Mining; Hypoxia; Metabolism; Tumor Microenvironment; Immunology; Treatment Resistance; BiomedBERT},
    booktitle    = {GeNeDIS 2026},
    series       = {Advances in Experimental Medicine and Biology (AEMB)},
    volume       = {xxxx},
    year         = {2026},
    pages        = {1-15},
    editor       = {Vlamos, Panagiotis},
    publisher    = {Springer Nature Switzerland},
    address      = {Cham, Switzerland},
    doi          = {},
    isbn         = {}
}
```

## Licensing

No blanket open-source/content license has been assigned. See [`docs/LICENSING.md`](docs/LICENSING.md) before public release.
