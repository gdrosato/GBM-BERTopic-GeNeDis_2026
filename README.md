# GBM BERTopic — GeNeDis 2026 reproducibility repository

Code, expert annotation, aggregate results, and figures supporting the study:

**AI-Driven Mapping of the Hypoxia–Metabolism–Immune Landscape in Glioblastoma: A BERTopic-Based Analysis of the Biomedical Literature**

Authors: **Stamatia Pouliliou** and **George Drosatos**

## Scope of this repository

This public repository intentionally contains:

- the Python scripts used in the analysis workflow;
- the exact PubMed query and analysis parameters;
- the exact Python/package/runtime environment export;
- the curated domain-expert topic-validation workbook;
- aggregate derived CSV/JSON outputs used for validation and temporal analysis;
- final main and supplementary figures.

It intentionally **does not contain**:

- the PubMed XML/XML.GZ export or article abstracts;
- parsed title–abstract corpora or embedding matrices;
- BERTopic pickle/model binaries;
- the manuscript source, LaTeX, BibTeX, or compiled manuscript PDF.

This separation keeps the repository focused on reproducible analysis code and author-generated/aggregate outputs while avoiding redistribution of abstract-bearing PubMed data and manuscript files.

## Key results

- PubMed XML records used in the analysis: **55,446 unique PMIDs**
- Final modeling corpus: **54,050 publications**
- Biomedical embedding dimensionality: **768**
- BERTopic topics: **110**
- Clustered publications: **35,188**
- HDBSCAN outliers: **18,862 (34.90%)**
- Candidate topics for expert review: **77**
- Strict expert-included topics: **41**

## Repository structure

```text
.
├── config/                 Exact PubMed query and model parameters
├── scripts/                Analysis and data-retrieval scripts
├── environment/            Exact Python/pip/runtime snapshot
├── data/
│   ├── raw/                Documentation only; no raw PubMed data committed
│   ├── validation/         Curated expert annotation workbook
│   └── derived/            Aggregate validation/temporal outputs
├── results/figures/        Main and supplementary publication figures
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

For the complete package snapshot see `environment/pip_freeze.txt`.

## Reproducing the analysis

The repository does not redistribute the raw PubMed export. Recreate it locally with the exact study query and NCBI EDirect. Full installation, API-key, download, compression, and PMID-verification instructions are provided in:

[`docs/PUBMED_RETRIEVAL.md`](docs/PUBMED_RETRIEVAL.md)

For an automated download using the repository query:

```bash
bash scripts/download_pubmed_edirect.sh
```

The historical XML dataset used for analysis contained **55,446 unique PMIDs**. PubMed is a live database, so later reruns may differ slightly because of retrospective record/indexing changes.

Then follow [`docs/PIPELINE.md`](docs/PIPELINE.md).

## Expert validation

The final curated workbook is:

`data/validation/topic_validation_workbook_2026-09-09.xlsx`

It is the final source of truth for the expert classification used by the downstream figure/temporal script.

## Data-distribution policy

The raw PubMed export and parsed abstract-bearing corpora are not included. See [`docs/DATA_POLICY.md`](docs/DATA_POLICY.md).

## Integrity

```bash
python tests/test_repo_tree.py
python tests/test_validation_workbook.py
```

## Citation

If you use this repo in a scientific publication, we would appreciate using the following citation:

- Pouliliou, S. and Drosatos, G. (2026). AI-Driven Mapping of the Hypoxia--Metabolism--Immune Landscape in Glioblastoma: A BERTopic-Based Analysis of the Biomedical Literature. In GeNeDIS 2026, pages 1-15, AEMB Vol. xxxx, Springer.

and as BibTeX:

```
@InProceedings{Pouliliou_GBM_2026,
    author       = {Pouliliou, Stamatia and Drosatos, George},
    title        = {AI-Driven Mapping of the Hypoxia--Metabolism--Immune Landscape in Glioblastoma: A BERTopic-Based Analysis of the Biomedical Literature},
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

## License

No open-source/content license has been assigned automatically. See [`docs/LICENSING.md`](docs/LICENSING.md) before public release.
