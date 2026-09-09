# Analysis scripts

These scripts implement the complete computational workflow described in
[`../docs/PIPELINE.md`](../docs/PIPELINE.md).

For PubMed retrieval, use:

```bash
bash scripts/download_pubmed_edirect.sh
```

For topic-validation workbook generation, use `build_topic_validation_fixed.py`.
The earlier `build_topic_validation.py` is retained to document the initial
implementation, which could create a duplicate `topic` column.

Several early analysis scripts use fixed filenames and should be run from the
project working directory; see `../docs/PIPELINE.md`.
