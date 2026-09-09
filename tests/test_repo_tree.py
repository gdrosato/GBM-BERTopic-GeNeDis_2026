from pathlib import Path
required = [
 'README.md', 'CITATION.cff', 'environment/pip_freeze.txt', 'environment/runtime.json',
 'scripts/pubmed_parser.py', 'scripts/dataset_qc.py', 'scripts/prepare_corpus.py',
 'scripts/token_length_analysis.py', 'scripts/generate_embeddings.py',
 'scripts/run_bertopic_baseline.py', 'scripts/tune_hdbscan.py',
 'scripts/run_bertopic_candidate_b.py', 'scripts/build_topic_validation.py',
 'scripts/build_topic_validation_fixed.py', 'scripts/make_final_figures_curated.py',
 'data/validation/topic_validation_workbook_Tina_Currated_2026-09-09.xlsx',
 'docs/PIPELINE.md', 'docs/DATA_POLICY.md'
]
missing = [x for x in required if not Path(x).exists()]
assert not missing, missing
for forbidden in ['manuscript', 'submission']:
    assert not Path(forbidden).exists(), f'Forbidden directory present: {forbidden}'
assert not list(Path('.').rglob('*.tex')), 'LaTeX source found unexpectedly'
assert not list(Path('.').rglob('*.bib')), 'BibTeX source found unexpectedly'
print('Repository tree/public-safety checks: OK')
