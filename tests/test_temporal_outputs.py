from pathlib import Path
import math
import pandas as pd

primary = pd.read_csv(Path('data/derived/primary_axis_trend_stats.csv'))
assert len(primary) == 4

expected_topics = {
    'Hypoxia / angiogenesis': 2,
    'Metabolism / stress': 3,
    'Immune / tumor microenvironment': 6,
    'Treatment resistance': 12,
}
expected_rho = {
    'Hypoxia / angiogenesis': -0.8142857143,
    'Metabolism / stress': 0.8428571429,
    'Immune / tumor microenvironment': 0.8678571429,
    'Treatment resistance': -0.4428571429,
}

for row in primary.itertuples(index=False):
    assert expected_topics[row.axis] == int(row.n_axis_topics)
    assert math.isclose(row.spearman_rho, expected_rho[row.axis], rel_tol=0, abs_tol=1e-9)

comparison = pd.read_csv(Path('data/derived/denominator_sensitivity_comparison.csv'))
assert len(comparison) == 4
assert comparison['same_rho_sign'].astype(bool).all()
assert comparison['same_significance_status'].astype(bool).all()
assert comparison['robust_direction_and_significance'].astype(bool).all()

axis_compare = pd.read_csv(Path('data/derived/axis_definition_comparison.csv'))
assert len(axis_compare) == 4

hyp = pd.read_csv(Path('data/derived/hypoxia_axis_decomposition.csv'))
assert set(hyp['publication_year']) == set(range(2011, 2026))
row2011 = hyp[hyp['publication_year'] == 2011].iloc[0]
row2025 = hyp[hyp['publication_year'] == 2025].iloc[0]
assert math.isclose(row2011['topic29_prevalence_percent'], 1.2309920348, abs_tol=1e-9)
assert math.isclose(row2025['topic29_prevalence_percent'], 0.4382929642, abs_tol=1e-9)

print('Temporal and sensitivity outputs: OK')
