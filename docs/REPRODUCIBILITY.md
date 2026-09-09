# Reproducibility notes

## Exact environment

The exact analysis environment snapshot is provided in `environment/`.

- Python 3.12.3
- Linux 6.8.0-138-generic x86_64, glibc 2.39
- PyTorch 2.13.0+cu130
- CUDA build 13.0
- NVIDIA A2 GPU

Use `environment/pip_freeze.txt` for exact installed package versions.

## Randomness

UMAP used random state 42. The final candidate was rerun and matched the selected HDBSCAN tuning labels/outlier status exactly.

## Large intermediates

Embeddings and serialized model objects are intentionally not committed. They can be regenerated from a locally retrieved PubMed corpus by following `docs/PIPELINE.md`.
