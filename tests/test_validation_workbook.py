from pathlib import Path
import pandas as pd

p = Path('data/validation/topic_validation_workbook_2026-09-09.xlsx')
df = pd.read_excel(p, sheet_name='Expert Annotation')
assert len(df) == 77
assert df['topic'].nunique() == 77
assert df['expert_include'].value_counts().to_dict() == {'Yes': 41, 'No': 24, 'Borderline': 12}
yes = df[df['expert_include'] == 'Yes']
for col in ['expert_primary_label','biological_coherence_1_to_5','relevance_to_hypoxia_metabolism_immune_resistance_1_to_5']:
    assert yes[col].notna().all(), col
print('Curated validation workbook: OK')
