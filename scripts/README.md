# Analysis scripts

These scripts implement the computational workflow described in [`../docs/PIPELINE.md`](../docs/PIPELINE.md).

For PubMed retrieval:

```bash
bash scripts/download_pubmed_edirect.sh
```

For topic-annotation workbook generation, use `build_topic_validation_fixed.py`. The earlier `build_topic_validation.py` is retained for provenance because it could create a duplicate `topic` column.

Final release-stage analysis uses:

- `make_final_figures_curated.py` for expert-curated summaries, primary-label temporal analysis, sentinel topics, and publication figures;
- `axis_definition_sensitivity.py` for primary-only versus primary+secondary axis definitions;
- `clustered_denominator_sensitivity.py` for the clustered-only denominator sensitivity check.

Several early scripts use fixed filenames and should be run from the project working directory; see `../docs/PIPELINE.md`.
