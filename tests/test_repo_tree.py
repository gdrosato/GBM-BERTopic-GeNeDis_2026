from pathlib import Path

required = [
    'environment/pip_freeze.txt',
    'environment/runtime.json',
    'config/pubmed_query.txt',
    'config/analysis_parameters.yaml',
    'CITATION.cff',
    'scripts/download_pubmed_edirect.sh',
    'scripts/pubmed_parser.py',
    'scripts/dataset_qc.py',
    'scripts/prepare_corpus.py',
    'scripts/token_length_analysis.py',
    'scripts/generate_embeddings.py',
    'scripts/run_bertopic_baseline.py',
    'scripts/tune_hdbscan.py',
    'scripts/run_bertopic_candidate_b.py',
    'scripts/build_topic_validation.py',
    'scripts/build_topic_validation_fixed.py',
    'scripts/make_final_figures_curated.py',
    'scripts/axis_definition_sensitivity.py',
    'scripts/clustered_denominator_sensitivity.py',
    'data/manifests/retrieved_pmids_2026-08-27.txt',
    'data/manifests/final_corpus_pmids.txt',
    'data/manifests/README.md',
    'data/validation/topic_validation_workbook_2026-09-09.xlsx',
    'data/derived/primary_axis_membership.csv',
    'data/derived/primary_axis_prevalence_annual.csv',
    'data/derived/primary_axis_trend_stats.csv',
    'data/derived/expanded_axis_membership.csv',
    'data/derived/expanded_axis_trend_stats.csv',
    'data/derived/axis_definition_comparison.csv',
    'data/derived/denominator_sensitivity_comparison.csv',
    'data/derived/hypoxia_axis_decomposition.csv',
    'results/figures/main/figure_topic_landscape.pdf',
    'results/figures/main/figure_temporal_axes_primary.pdf',
    'results/figures/supplementary/figure_temporal_axes_expanded_sensitivity.pdf',
    'docs/PUBMED_RETRIEVAL.md',
    'docs/PIPELINE.md',
    'docs/DATA_POLICY.md',
    'docs/CODEBOOK.md',
]

missing = [x for x in required if not Path(x).exists()]
assert not missing, missing

for forbidden in ['manuscript', 'submission']:
    assert not Path(forbidden).exists(), f'Forbidden directory present: {forbidden}'

assert not list(Path('.').rglob('*.tex')), 'LaTeX source found unexpectedly'
assert not list(Path('.').rglob('*.bib')), 'BibTeX source found unexpectedly'
assert not list(Path('.').rglob('*.xml')), 'Raw XML found unexpectedly'
assert not list(Path('.').rglob('*.xml.gz')), 'Compressed raw XML found unexpectedly'

# Superseded ambiguous main-axis outputs should not be distributed.
for stale in [
    'data/derived/axis_topic_membership.csv',
    'data/derived/core_axis_prevalence_annual.csv',
    'data/derived/temporal_axis_trend_stats.csv',
    'results/figures/main/figure_temporal_core_axes.pdf',
]:
    assert not Path(stale).exists(), f'Superseded file present: {stale}'

print('Repository tree/public-safety checks: OK')
