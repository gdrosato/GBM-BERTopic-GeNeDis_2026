# Changelog

## v1.2.0 — 2026-09-09

- Added PMID-only manifests identifying the exact 55,446-record historical PubMed retrieval and 54,050-document final modeling corpus.
- Added manifest provenance, SHA-256 fingerprints, and a corpus-manifest integrity test.
- Updated README, data policy, provenance, retrieval, pipeline, reproducibility, and release documentation to distinguish reproducible PubMed retrieval from exact historical corpus identity.

## v1.1.0 — 2026-09-09

- Updated the release workflow to use expert **primary labels only** for the main four-axis temporal analysis.
- Added axis-definition sensitivity analysis comparing primary labels with primary + secondary labels.
- Added clustered-only denominator sensitivity analysis.
- Added primary-axis and expanded-axis aggregate outputs and comparison files.
- Added post-hoc hypoxia-axis decomposition output for Topics 29 and 79.
- Replaced the main temporal figure with `figure_temporal_axes_primary` and moved the broader axis-definition figure to supplementary results.
- Updated Figure 1 terminology from "strict expert-included" to "expert-retained".
- Updated documentation, codebook, CI tests, and citation metadata to match the final analysis.
- Removed the obsolete checksum-verification step from repository maintenance.

## v1.0.0 — 2026-09-09

- Added the complete Python analysis workflow, including an EDirect retrieval helper.
- Added exact Python/package/runtime environment snapshot.
- Added EDirect retrieval documentation for recreating the PubMed XML locally.
- Added final curated expert-annotation workbook.
- Added aggregate validation and temporal outputs.
- Added main and supplementary figures.
- Explicitly excluded raw PubMed export, abstract-bearing corpora, model binaries, and manuscript files.
