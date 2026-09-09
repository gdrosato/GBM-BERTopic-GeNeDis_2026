# Derived-data codebook

## Expert curation

### `validation_all_candidates_curated.csv`
Curated rows for all 77 algorithmically prioritized candidate topics.

### `validated_topics_strict.csv`
The 41 expert-retained (`Yes`) topics used for the principal biological interpretation.

### `validation_summary.csv` / `validation_summary.json`
Aggregate curation counts and expert score summaries.

## Primary temporal analysis

### `primary_axis_membership.csv`
Mapping from expert-retained topics to the four main biological axes using **expert primary labels only**. The four axes contain 2, 3, 6, and 12 topics for hypoxia/angiogenesis, metabolism/stress, immune/tumor microenvironment, and treatment resistance, respectively.

### `primary_axis_prevalence_annual.csv`
Annual 2011–2025 primary-axis counts, total eligible GBM denominator, raw prevalence percentage, centered 3-year rolling prevalence, and topic count per axis.

### `primary_axis_trend_stats.csv`
Spearman trend statistics for the four primary-label axes with Benjamini–Hochberg-adjusted q values.

## Axis-definition sensitivity analysis

### `expanded_axis_membership.csv`
Broader mapping that allows expert primary and secondary labels to assign a topic to one or more biological axes.

### `expanded_axis_prevalence_annual.csv`
Annual prevalence for the broader primary-plus-secondary axis definition.

### `expanded_axis_trend_stats.csv`
Spearman/BH results for the broader axis definition.

### `axis_definition_comparison.csv`
Side-by-side comparison of topic counts, Spearman coefficients, q values, trend direction, and significance status for primary-only versus expanded axis definitions.

### `axis_definition_results.txt`
Human-readable summary of the axis-definition sensitivity analysis.

## Denominator sensitivity analysis

### `denominator_sensitivity_annual.csv`
Annual primary-axis prevalence calculated using both the all-eligible-GBM denominator and the clustered-publications-only denominator.

### `denominator_sensitivity_stats.csv`
Detailed Spearman/BH statistics for both denominator definitions.

### `denominator_sensitivity_comparison.csv`
Compact comparison showing whether trend direction and significance status are preserved under clustered-only normalization.

### `denominator_sensitivity_results.txt`
Human-readable denominator-sensitivity summary.

## Sentinel topics and hypoxia decomposition

### `sentinel_topic_prevalence_annual.csv`
Annual prevalence for the four preselected sentinel topics: Topic 79 (hypoxia/HIF-1), Topic 16 (metabolism/glucose), Topic 11 (macrophage/microglia), and Topic 63 (TMZ resistance).

### `sentinel_topic_trend_stats.csv`
Spearman trend statistics for the four sentinel topics, with a separate Benjamini–Hochberg correction family.

### `hypoxia_axis_decomposition.csv`
Post-hoc descriptive decomposition of the primary hypoxia/angiogenesis axis into Topic 29 (VEGF/angiogenesis) and Topic 79 (hypoxia/HIF-1) annual counts and prevalence.
