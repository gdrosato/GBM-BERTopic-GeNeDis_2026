#!/usr/bin/env python3
"""
Generate publication figures and temporal summaries for the GeNeDis 2026
glioblastoma BERTopic manuscript using the CURATED expert-validation workbook.

The validation workbook is the source-of-truth. No hard-coded list of all
included topics is required.

Main outputs
------------
1. figure_topic_landscape.pdf/png
   A separate 2D UMAP visualization for presentation only. All documents are
   shown in the background; documents in strict expert-included topics are
   highlighted by curated expert primary category.

2. figure_temporal_core_axes.pdf/png
   2011-2025 normalized prevalence of four NON-MUTUALLY-EXCLUSIVE biological
   axes derived from expert primary + secondary labels:
       Hypoxia / angiogenesis
       Metabolism / stress
       Immune / tumor microenvironment
       Treatment resistance
   A centered 3-year rolling mean is used for visualization.

3. figure_temporal_sentinel_topics.pdf/png
   Topic-specific trajectories for four mechanistically interpretable sentinel
   topics (79, 16, 11, 63), provided that they remain strict-included.

4. figure_validation_outcomes_curated.pdf/png
   Included / Borderline / Excluded validation counts.

5. figure_validation_primary_labels_curated.pdf/png
   Curated primary-label distribution among strict-included topics.

6. figure_outlier_rate_by_year.pdf/png
   HDBSCAN temporal coverage QC; recommended for Supplementary Material.

Tabular outputs
---------------
validated_topics_strict.csv
validation_summary.csv
axis_topic_membership.csv
core_axis_prevalence_annual.csv
temporal_axis_trend_stats.csv
sentinel_topic_prevalence_annual.csv

Expected analysis inputs
---------------------------
data/validation/topic_validation_workbook_2026-09-09.xlsx
embeddings/biomedbert_embeddings.npy
bertopic_candidate_b/topic_assignments.parquet
bertopic_candidate_b/topic_prevalence_by_year.csv
bertopic_candidate_b/outliers_by_year.csv

Notes
-----
- "Strict included" means expert_include == "Yes".
- Borderline topics are not used in the main biological temporal analysis.
- Core biological axes are allowed to overlap conceptually. A topic can
  contribute to more than one axis if the expert primary/secondary labels
  support multiple axes. Therefore axis lines are NOT intended to sum to 100%.
- The 2D UMAP is generated separately for visualization. The clustering model
  itself remains based on the previously fixed 5D UMAP representation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from umap import UMAP


PRIMARY_DISPLAY_MAP = {
    "Hypoxia": "Hypoxia",
    "Metabolism": "Metabolism",
    "Immune microenvironment": "Immune / tumor microenvironment",
    "Tumor microenvironment": "Immune / tumor microenvironment",
    "Treatment resistance": "Treatment resistance",
    "Therapeutic intervention": "Therapeutic intervention",
}

CORE_AXIS_ORDER = [
    "Hypoxia / angiogenesis",
    "Metabolism / stress",
    "Immune / tumor microenvironment",
    "Treatment resistance",
]

SENTINEL_TOPICS = {
    "Hypoxia / HIF-1": 79,
    "Metabolism / glucose": 16,
    "Macrophage / microglia": 11,
    "TMZ resistance": 63,
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--validation-workbook",
        type=Path,
        default=Path("data/validation/topic_validation_workbook_2026-09-09.xlsx"),
    )
    p.add_argument(
        "--embeddings",
        type=Path,
        default=Path("embeddings/biomedbert_embeddings.npy"),
    )
    p.add_argument(
        "--assignments",
        type=Path,
        default=Path("bertopic_candidate_b/topic_assignments.parquet"),
    )
    p.add_argument(
        "--prevalence",
        type=Path,
        default=Path("bertopic_candidate_b/topic_prevalence_by_year.csv"),
    )
    p.add_argument(
        "--outliers",
        type=Path,
        default=Path("bertopic_candidate_b/outliers_by_year.csv"),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("final_figures_curated"),
    )
    p.add_argument(
        "--umap2d-cache",
        type=Path,
        default=Path("bertopic_candidate_b/umap_2d_visualization.npy"),
    )
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--start-year", type=int, default=2011)
    p.add_argument("--end-year", type=int, default=2025)
    return p.parse_args()


def clean_label(value):
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def split_secondary(value):
    text = clean_label(value)
    if not text:
        return []
    return [x.strip() for x in text.split(";") if x.strip()]


def load_validation(path):
    df = pd.read_excel(path, sheet_name="Expert Annotation")

    required = {
        "topic",
        "n_documents",
        "Name",
        "expert_include",
        "expert_primary_label",
        "expert_secondary_labels",
        "biological_coherence_1_to_5",
        "relevance_to_hypoxia_metabolism_immune_resistance_1_to_5",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Validation workbook missing columns: {sorted(missing)}"
        )

    df = df.copy()
    df["topic"] = pd.to_numeric(df["topic"], errors="raise").astype(int)
    df["n_documents"] = pd.to_numeric(
        df["n_documents"], errors="raise"
    ).astype(int)
    df["expert_include"] = df["expert_include"].map(clean_label)
    df["expert_primary_label"] = df["expert_primary_label"].map(clean_label)
    df["expert_secondary_labels"] = df["expert_secondary_labels"].map(clean_label)

    allowed = {"Yes", "Borderline", "No"}
    observed = set(df["expert_include"].dropna().unique())
    unexpected = observed - allowed
    if unexpected:
        raise ValueError(
            f"Unexpected expert_include values: {sorted(unexpected)}"
        )

    strict = df[df["expert_include"] == "Yes"].copy()

    if strict[
        [
            "expert_primary_label",
            "biological_coherence_1_to_5",
            "relevance_to_hypoxia_metabolism_immune_resistance_1_to_5",
        ]
    ].isna().any().any():
        raise ValueError(
            "At least one strict-included topic lacks a primary label or score."
        )

    strict["biological_coherence_1_to_5"] = pd.to_numeric(
        strict["biological_coherence_1_to_5"], errors="raise"
    )
    strict[
        "relevance_to_hypoxia_metabolism_immune_resistance_1_to_5"
    ] = pd.to_numeric(
        strict[
            "relevance_to_hypoxia_metabolism_immune_resistance_1_to_5"
        ],
        errors="raise",
    )

    return df, strict


def topic_label_set(row):
    labels = {clean_label(row["expert_primary_label"])}
    labels.update(split_secondary(row["expert_secondary_labels"]))
    labels.discard("")
    return labels


def assign_core_axes(strict):
    rows = []

    for _, row in strict.iterrows():
        labels = topic_label_set(row)
        axes = []

        if "Hypoxia" in labels:
            axes.append("Hypoxia / angiogenesis")

        if "Metabolism" in labels:
            axes.append("Metabolism / stress")

        if (
            "Tumor microenvironment" in labels
            or "Immune microenvironment" in labels
        ):
            axes.append("Immune / tumor microenvironment")

        if "Treatment resistance" in labels:
            axes.append("Treatment resistance")

        for axis in axes:
            rows.append(
                {
                    "topic": int(row["topic"]),
                    "n_documents": int(row["n_documents"]),
                    "bertopic_name": row["Name"],
                    "expert_primary_label": row["expert_primary_label"],
                    "expert_secondary_labels": row["expert_secondary_labels"],
                    "coherence": float(row["biological_coherence_1_to_5"]),
                    "relevance": float(
                        row[
                            "relevance_to_hypoxia_metabolism_immune_resistance_1_to_5"
                        ]
                    ),
                    "axis": axis,
                }
            )

    return pd.DataFrame(rows)


def save_validation_summaries(all_validation, strict, axis_membership, outdir):
    all_validation.to_csv(
        outdir / "validation_all_candidates_curated.csv",
        index=False,
    )
    strict.to_csv(
        outdir / "validated_topics_strict.csv",
        index=False,
    )
    axis_membership.to_csv(
        outdir / "axis_topic_membership.csv",
        index=False,
    )

    status_order = ["Yes", "Borderline", "No"]
    status_counts = (
        all_validation["expert_include"]
        .value_counts()
        .reindex(status_order, fill_value=0)
    )

    primary_counts = (
        strict["expert_primary_label"]
        .value_counts()
        .sort_values(ascending=False)
    )

    summary = {
        "n_candidates": int(len(all_validation)),
        "n_strict_included": int(len(strict)),
        "n_borderline": int(status_counts["Borderline"]),
        "n_excluded": int(status_counts["No"]),
        "strict_included_percentage": float(
            100 * len(strict) / len(all_validation)
        ),
        "coherence_mean": float(
            strict["biological_coherence_1_to_5"].mean()
        ),
        "coherence_median": float(
            strict["biological_coherence_1_to_5"].median()
        ),
        "coherence_min": float(
            strict["biological_coherence_1_to_5"].min()
        ),
        "relevance_mean": float(
            strict[
                "relevance_to_hypoxia_metabolism_immune_resistance_1_to_5"
            ].mean()
        ),
        "relevance_median": float(
            strict[
                "relevance_to_hypoxia_metabolism_immune_resistance_1_to_5"
            ].median()
        ),
        "relevance_ge_4": int(
            (
                strict[
                    "relevance_to_hypoxia_metabolism_immune_resistance_1_to_5"
                ]
                >= 4
            ).sum()
        ),
        "relevance_ge_4_percentage": float(
            100
            * (
                strict[
                    "relevance_to_hypoxia_metabolism_immune_resistance_1_to_5"
                ]
                >= 4
            ).mean()
        ),
        "status_counts": {
            k: int(v) for k, v in status_counts.to_dict().items()
        },
        "primary_label_counts": {
            k: int(v) for k, v in primary_counts.to_dict().items()
        },
    }

    (outdir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    pd.DataFrame(
        [
            {"metric": "Candidate topics", "value": len(all_validation)},
            {"metric": "Strict included (Yes)", "value": len(strict)},
            {"metric": "Borderline", "value": status_counts["Borderline"]},
            {"metric": "Excluded (No)", "value": status_counts["No"]},
            {
                "metric": "Mean biological coherence",
                "value": summary["coherence_mean"],
            },
            {
                "metric": "Median biological coherence",
                "value": summary["coherence_median"],
            },
            {
                "metric": "Mean framework relevance",
                "value": summary["relevance_mean"],
            },
            {
                "metric": "Median framework relevance",
                "value": summary["relevance_median"],
            },
            {
                "metric": "Included topics with relevance >=4",
                "value": summary["relevance_ge_4"],
            },
        ]
    ).to_csv(
        outdir / "validation_summary.csv",
        index=False,
    )

    return summary


def plot_validation_figures(all_validation, strict, outdir):
    status_counts = (
        all_validation["expert_include"]
        .value_counts()
        .reindex(["Yes", "Borderline", "No"], fill_value=0)
    )

    labels = ["Included", "Borderline", "Excluded"]
    values = [
        int(status_counts["Yes"]),
        int(status_counts["Borderline"]),
        int(status_counts["No"]),
    ]

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    bars = ax.bar(labels, values)
    ax.set_ylabel("Number of topics")
    ax.set_title("Domain-expert validation outcome")
    ax.set_ylim(0, max(values) * 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    total = len(all_validation)
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.5,
            f"{value} ({100*value/total:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(
        outdir / "figure_validation_outcomes_curated.pdf",
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "figure_validation_outcomes_curated.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    primary = (
        strict["expert_primary_label"]
        .value_counts()
        .sort_values(ascending=False)
    )

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    xlabels = [x.replace(" ", "\n", 1) for x in primary.index]
    bars = ax.bar(xlabels, primary.values)
    ax.set_ylabel("Strict expert-included topics")
    ax.set_title("Primary expert labels among included topics")
    ax.set_ylim(0, max(primary.values) * 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for bar, value in zip(bars, primary.values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.25,
            str(int(value)),
            ha="center",
            va="bottom",
            fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(
        outdir / "figure_validation_primary_labels_curated.pdf",
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "figure_validation_primary_labels_curated.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def validate_assignment_alignment(assignments, embeddings):
    required = {"embedding_row", "topic"}
    missing = required - set(assignments.columns)
    if missing:
        raise ValueError(f"Assignments missing columns: {sorted(missing)}")

    assignments = assignments.sort_values(
        "embedding_row"
    ).reset_index(drop=True)

    expected = np.arange(len(assignments), dtype=np.int64)
    actual = assignments["embedding_row"].to_numpy(dtype=np.int64)

    if not np.array_equal(expected, actual):
        raise ValueError("embedding_row is not exactly 0..N-1.")

    if len(assignments) != embeddings.shape[0]:
        raise ValueError(
            f"Assignments rows ({len(assignments):,}) != "
            f"embedding rows ({embeddings.shape[0]:,})."
        )

    return assignments


def make_2d_umap(embeddings, cache_path, random_state):
    if cache_path.exists():
        coords = np.load(cache_path)
        if coords.shape == (embeddings.shape[0], 2):
            print(f"Using cached 2D UMAP: {cache_path}")
            return coords
        print("Cached 2D UMAP shape mismatch; recomputing.")

    print("Computing separate 2D UMAP for visualization...")
    model = UMAP(
        n_neighbors=15,
        n_components=2,
        min_dist=0.1,
        metric="cosine",
        random_state=random_state,
        transform_seed=random_state,
        low_memory=True,
        verbose=True,
    )

    coords = model.fit_transform(embeddings).astype(np.float32)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(cache_path, coords)
    return coords


def display_group(primary):
    return PRIMARY_DISPLAY_MAP.get(
        clean_label(primary),
        clean_label(primary) or "Other",
    )


def select_centroid_labels(strict, max_per_group=3):
    tmp = strict.copy()
    tmp["display_group"] = tmp["expert_primary_label"].map(display_group)

    selected = []

    for _, group in tmp.groupby("display_group", sort=True):
        group = group.sort_values(
            [
                "relevance_to_hypoxia_metabolism_immune_resistance_1_to_5",
                "biological_coherence_1_to_5",
                "n_documents",
            ],
            ascending=[False, False, False],
        ).head(max_per_group)
        selected.extend(group["topic"].astype(int).tolist())

    return sorted(set(selected))


def plot_topic_landscape(assignments, coords, strict, outdir):
    topic_to_group = {
        int(row.topic): display_group(row.expert_primary_label)
        for row in strict.itertuples()
    }

    topics = assignments["topic"].to_numpy(dtype=int)

    fig, ax = plt.subplots(figsize=(9.2, 7.2))

    # All documents as visual context.
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        s=2,
        alpha=0.05,
        linewidths=0,
    )

    ordered_groups = [
        "Hypoxia",
        "Metabolism",
        "Immune / tumor microenvironment",
        "Treatment resistance",
        "Therapeutic intervention",
    ]

    for group_name in ordered_groups:
        topic_ids = [
            t for t, g in topic_to_group.items()
            if g == group_name
        ]
        if not topic_ids:
            continue

        mask = np.isin(topics, topic_ids)

        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=8,
            alpha=0.68,
            linewidths=0,
            label=group_name,
        )

    for topic_id in select_centroid_labels(strict, max_per_group=3):
        mask = topics == topic_id
        if not np.any(mask):
            continue

        cx = float(np.median(coords[mask, 0]))
        cy = float(np.median(coords[mask, 1]))

        ax.text(
            cx,
            cy,
            str(topic_id),
            fontsize=7.5,
            ha="center",
            va="center",
        )

    ax.set_title(
        "Semantic landscape of the glioblastoma literature\n"
        "with strict expert-included topics highlighted"
    )
    ax.set_xlabel("UMAP dimension 1")
    ax.set_ylabel("UMAP dimension 2")
    ax.legend(frameon=False, loc="best", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    fig.savefig(
        outdir / "figure_topic_landscape.pdf",
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "figure_topic_landscape.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def validate_prevalence(prevalence):
    required = {
        "publication_year",
        "topic",
        "topic_documents",
        "all_gbm_documents",
    }
    missing = required - set(prevalence.columns)
    if missing:
        raise ValueError(
            f"Prevalence file missing columns: {sorted(missing)}"
        )

    prevalence = prevalence.copy()
    prevalence["publication_year"] = pd.to_numeric(
        prevalence["publication_year"], errors="raise"
    ).astype(int)
    prevalence["topic"] = pd.to_numeric(
        prevalence["topic"], errors="raise"
    ).astype(int)
    prevalence["topic_documents"] = pd.to_numeric(
        prevalence["topic_documents"], errors="raise"
    )
    prevalence["all_gbm_documents"] = pd.to_numeric(
        prevalence["all_gbm_documents"], errors="raise"
    )

    return prevalence


def build_axis_prevalence(prevalence, axis_membership, start_year, end_year):
    prevalence = validate_prevalence(prevalence)

    denominators = (
        prevalence[
            ["publication_year", "all_gbm_documents"]
        ]
        .drop_duplicates()
        .sort_values("publication_year")
    )

    # Safety: one denominator per year.
    if denominators["publication_year"].duplicated().any():
        raise ValueError(
            "Multiple all_gbm_documents values found for one publication year."
        )

    years = pd.DataFrame(
        {"publication_year": np.arange(start_year, end_year + 1, dtype=int)}
    )

    rows = []

    for axis in CORE_AXIS_ORDER:
        topic_ids = (
            axis_membership.loc[
                axis_membership["axis"] == axis,
                "topic",
            ]
            .astype(int)
            .unique()
        )

        tmp = prevalence[
            prevalence["topic"].isin(topic_ids)
        ].copy()

        counts = (
            tmp.groupby("publication_year", as_index=False)[
                "topic_documents"
            ]
            .sum()
            .rename(
                columns={
                    "topic_documents": "axis_topic_documents"
                }
            )
        )

        annual = (
            years
            .merge(counts, on="publication_year", how="left")
            .merge(
                denominators,
                on="publication_year",
                how="left",
                validate="one_to_one",
            )
        )

        annual["axis_topic_documents"] = (
            annual["axis_topic_documents"].fillna(0)
        )

        if annual["all_gbm_documents"].isna().any():
            missing_years = annual.loc[
                annual["all_gbm_documents"].isna(),
                "publication_year",
            ].tolist()
            raise ValueError(
                f"Missing annual GBM denominator for years: {missing_years}"
            )

        annual["axis_prevalence_percent"] = (
            100.0
            * annual["axis_topic_documents"]
            / annual["all_gbm_documents"]
        )

        annual["rolling_3y_percent"] = (
            annual["axis_prevalence_percent"]
            .rolling(
                window=3,
                center=True,
                min_periods=1,
            )
            .mean()
        )

        annual["axis"] = axis
        annual["n_axis_topics"] = len(topic_ids)

        rows.append(annual)

    return pd.concat(rows, ignore_index=True)


def bh_fdr(p_values):
    p = np.asarray(p_values, dtype=float)
    n = len(p)

    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)

    # Monotonic BH adjustment from largest rank downward.
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)

    result = np.empty(n, dtype=float)
    result[order] = adjusted
    return result


def temporal_trend_stats(axis_annual):
    rows = []

    for axis in CORE_AXIS_ORDER:
        g = (
            axis_annual[axis_annual["axis"] == axis]
            .sort_values("publication_year")
        )

        rho, p = spearmanr(
            g["publication_year"].to_numpy(),
            g["axis_prevalence_percent"].to_numpy(),
        )

        slope = np.polyfit(
            g["publication_year"].to_numpy(dtype=float),
            g["axis_prevalence_percent"].to_numpy(dtype=float),
            deg=1,
        )[0]

        rows.append(
            {
                "axis": axis,
                "n_years": len(g),
                "spearman_rho": float(rho),
                "spearman_p": float(p),
                "linear_slope_percentage_points_per_year": float(slope),
                "start_prevalence_percent": float(
                    g.iloc[0]["axis_prevalence_percent"]
                ),
                "end_prevalence_percent": float(
                    g.iloc[-1]["axis_prevalence_percent"]
                ),
            }
        )

    out = pd.DataFrame(rows)
    out["spearman_q_bh"] = bh_fdr(out["spearman_p"].to_numpy())
    return out


def sentinel_trend_stats(sentinel):
    rows = []

    for label, topic_id in SENTINEL_TOPICS.items():
        g = (
            sentinel[sentinel["sentinel_label"] == label]
            .sort_values("publication_year")
        )

        rho, p = spearmanr(
            g["publication_year"].to_numpy(),
            g["prevalence_percent"].to_numpy(),
        )

        slope = np.polyfit(
            g["publication_year"].to_numpy(dtype=float),
            g["prevalence_percent"].to_numpy(dtype=float),
            deg=1,
        )[0]

        rows.append(
            {
                "sentinel_label": label,
                "topic": topic_id,
                "n_years": len(g),
                "spearman_rho": float(rho),
                "spearman_p": float(p),
                "linear_slope_percentage_points_per_year": float(slope),
                "start_prevalence_percent": float(
                    g.iloc[0]["prevalence_percent"]
                ),
                "end_prevalence_percent": float(
                    g.iloc[-1]["prevalence_percent"]
                ),
            }
        )

    out = pd.DataFrame(rows)
    out["spearman_q_bh"] = bh_fdr(out["spearman_p"].to_numpy())
    return out


def trend_word(rho, q):
    if q < 0.05 and rho > 0:
        return "increased significantly"
    if q < 0.05 and rho < 0:
        return "decreased significantly"
    if rho > 0:
        return "showed a positive but non-significant monotonic trend"
    if rho < 0:
        return "showed a negative but non-significant monotonic trend"
    return "showed no monotonic trend"


def write_temporal_latex_summary(axis_stats, sentinel_stats, outdir):
    axis_sentences = []

    for row in axis_stats.itertuples():
        axis_sentences.append(
            f"The {row.axis} axis {trend_word(row.spearman_rho, row.spearman_q_bh)} "
            f"(Spearman $\\\\rho={row.spearman_rho:.2f}$, "
            f"$p={row.spearman_p:.3g}$, BH-adjusted $q={row.spearman_q_bh:.3g}$)."
        )

    sentinel_sentences = []

    for row in sentinel_stats.itertuples():
        sentinel_sentences.append(
            f"Topic {row.topic} ({row.sentinel_label}) "
            f"{trend_word(row.spearman_rho, row.spearman_q_bh)} "
            f"($\\\\rho={row.spearman_rho:.2f}$, "
            f"$p={row.spearman_p:.3g}$, BH-adjusted $q={row.spearman_q_bh:.3g}$)."
        )

    latex_text = (
        "During 2011--2025, "
        + " ".join(axis_sentences)
        + "\n\nSentinel-topic analysis showed that "
        + " ".join(sentinel_sentences)
        + "\n"
    )

    (outdir / "temporal_results_latex.txt").write_text(
        latex_text,
        encoding="utf-8",
    )


def plot_temporal_axes(axis_annual, outdir):
    fig, ax = plt.subplots(figsize=(9.2, 5.4))

    for axis in CORE_AXIS_ORDER:
        g = (
            axis_annual[axis_annual["axis"] == axis]
            .sort_values("publication_year")
        )

        ax.plot(
            g["publication_year"],
            g["rolling_3y_percent"],
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=axis,
        )

    ax.set_title(
        "Normalized prevalence of expert-validated biological axes, 2011-2025"
    )
    ax.set_xlabel("Publication year")
    ax.set_ylabel("Share of all eligible GBM publications (%)")
    ax.set_xticks(
        np.arange(
            axis_annual["publication_year"].min(),
            axis_annual["publication_year"].max() + 1,
            2,
        )
    )
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    fig.savefig(
        outdir / "figure_temporal_core_axes.pdf",
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "figure_temporal_core_axes.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def build_sentinel_prevalence(
    prevalence,
    strict,
    start_year,
    end_year,
):
    prevalence = validate_prevalence(prevalence)
    strict_ids = set(strict["topic"].astype(int))

    missing = [
        (label, topic)
        for label, topic in SENTINEL_TOPICS.items()
        if topic not in strict_ids
    ]
    if missing:
        raise ValueError(
            "Sentinel topic(s) are no longer strict-included: "
            + repr(missing)
        )

    denominators = (
        prevalence[
            ["publication_year", "all_gbm_documents"]
        ]
        .drop_duplicates()
        .sort_values("publication_year")
    )

    years = pd.DataFrame(
        {"publication_year": np.arange(start_year, end_year + 1, dtype=int)}
    )

    rows = []

    for label, topic_id in SENTINEL_TOPICS.items():
        tmp = prevalence[
            prevalence["topic"] == topic_id
        ][
            [
                "publication_year",
                "topic_documents",
            ]
        ].copy()

        annual = (
            years
            .merge(tmp, on="publication_year", how="left")
            .merge(
                denominators,
                on="publication_year",
                how="left",
                validate="one_to_one",
            )
        )

        annual["topic_documents"] = (
            annual["topic_documents"].fillna(0)
        )

        annual["prevalence_percent"] = (
            100.0
            * annual["topic_documents"]
            / annual["all_gbm_documents"]
        )

        annual["rolling_3y_percent"] = (
            annual["prevalence_percent"]
            .rolling(
                window=3,
                center=True,
                min_periods=1,
            )
            .mean()
        )

        annual["sentinel_label"] = label
        annual["topic"] = topic_id
        rows.append(annual)

    return pd.concat(rows, ignore_index=True)


def plot_sentinel_topics(sentinel, outdir):
    fig, ax = plt.subplots(figsize=(9.2, 5.4))

    for label in SENTINEL_TOPICS:
        g = (
            sentinel[sentinel["sentinel_label"] == label]
            .sort_values("publication_year")
        )

        ax.plot(
            g["publication_year"],
            g["rolling_3y_percent"],
            marker="o",
            markersize=3,
            linewidth=1.8,
            label=f"{label} (Topic {int(g.iloc[0]['topic'])})",
        )

    ax.set_title(
        "Temporal trajectories of sentinel expert-validated GBM topics, 2011-2025"
    )
    ax.set_xlabel("Publication year")
    ax.set_ylabel("Share of all eligible GBM publications (%)")
    ax.set_xticks(
        np.arange(
            sentinel["publication_year"].min(),
            sentinel["publication_year"].max() + 1,
            2,
        )
    )
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    fig.savefig(
        outdir / "figure_temporal_sentinel_topics.pdf",
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "figure_temporal_sentinel_topics.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_outliers(outliers, outdir):
    required = {"publication_year", "outlier_percentage"}
    missing = required - set(outliers.columns)
    if missing:
        raise ValueError(
            f"Outlier file missing columns: {sorted(missing)}"
        )

    tmp = outliers.copy()
    tmp["publication_year"] = pd.to_numeric(
        tmp["publication_year"], errors="raise"
    ).astype(int)
    tmp["outlier_percentage"] = pd.to_numeric(
        tmp["outlier_percentage"], errors="raise"
    )
    tmp = tmp.sort_values("publication_year")

    fig, ax = plt.subplots(figsize=(8.8, 4.8))

    ax.plot(
        tmp["publication_year"],
        tmp["outlier_percentage"],
        marker="o",
        markersize=3,
    )
    ax.axvline(2011, linestyle="--", linewidth=1)

    ax.set_title("HDBSCAN outlier prevalence by publication year")
    ax.set_xlabel("Publication year")
    ax.set_ylabel("Outlier publications (%)")
    ax.set_xticks(
        np.arange(
            tmp["publication_year"].min(),
            tmp["publication_year"].max() + 1,
            3,
        )
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    fig.savefig(
        outdir / "figure_outlier_rate_by_year.pdf",
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "figure_outlier_rate_by_year.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_validation, strict = load_validation(
        args.validation_workbook
    )

    axis_membership = assign_core_axes(strict)

    summary = save_validation_summaries(
        all_validation,
        strict,
        axis_membership,
        args.output_dir,
    )

    plot_validation_figures(
        all_validation,
        strict,
        args.output_dir,
    )

    print("\n=== CURATED VALIDATION QC ===")
    print(f"Candidate topics:       {summary['n_candidates']}")
    print(f"Strict included:        {summary['n_strict_included']}")
    print(f"Borderline:             {summary['n_borderline']}")
    print(f"Excluded:               {summary['n_excluded']}")
    print(
        "Coherence mean/median: "
        f"{summary['coherence_mean']:.2f} / "
        f"{summary['coherence_median']:.1f}"
    )
    print(
        "Relevance mean/median: "
        f"{summary['relevance_mean']:.2f} / "
        f"{summary['relevance_median']:.1f}"
    )
    print(
        "Relevance >=4:         "
        f"{summary['relevance_ge_4']} "
        f"({summary['relevance_ge_4_percentage']:.1f}%)"
    )

    # Semantic landscape.
    if args.embeddings.exists() and args.assignments.exists():
        embeddings = np.load(
            args.embeddings,
            mmap_mode="r",
        )
        assignments = pd.read_parquet(
            args.assignments
        )
        assignments = validate_assignment_alignment(
            assignments,
            embeddings,
        )

        coords = make_2d_umap(
            embeddings,
            args.umap2d_cache,
            args.random_state,
        )

        plot_topic_landscape(
            assignments,
            coords,
            strict,
            args.output_dir,
        )
    else:
        print(
            "Skipping semantic landscape: embeddings/assignments not found."
        )

    # Temporal analysis.
    if args.prevalence.exists():
        prevalence = pd.read_csv(
            args.prevalence
        )

        axis_annual = build_axis_prevalence(
            prevalence,
            axis_membership,
            args.start_year,
            args.end_year,
        )

        axis_annual.to_csv(
            args.output_dir / "core_axis_prevalence_annual.csv",
            index=False,
        )

        trend_stats = temporal_trend_stats(
            axis_annual
        )

        trend_stats.to_csv(
            args.output_dir / "temporal_axis_trend_stats.csv",
            index=False,
        )

        plot_temporal_axes(
            axis_annual,
            args.output_dir,
        )

        sentinel = build_sentinel_prevalence(
            prevalence,
            strict,
            args.start_year,
            args.end_year,
        )

        sentinel.to_csv(
            args.output_dir / "sentinel_topic_prevalence_annual.csv",
            index=False,
        )

        sentinel_stats = sentinel_trend_stats(sentinel)
        sentinel_stats.to_csv(
            args.output_dir / "sentinel_topic_trend_stats.csv",
            index=False,
        )

        write_temporal_latex_summary(
            trend_stats,
            sentinel_stats,
            args.output_dir,
        )

        plot_sentinel_topics(
            sentinel,
            args.output_dir,
        )

        print("\nTemporal axis trend statistics:")
        print(
            trend_stats.to_string(
                index=False,
                float_format=lambda x: f"{x:.4f}",
            )
        )

        print("\nSentinel topic trend statistics:")
        print(
            sentinel_stats.to_string(
                index=False,
                float_format=lambda x: f"{x:.4f}",
            )
        )

        print(
            "\nReady-to-paste LaTeX temporal summary: "
            f"{args.output_dir / 'temporal_results_latex.txt'}"
        )
    else:
        print(
            "Skipping temporal figures: prevalence CSV not found."
        )

    if args.outliers.exists():
        plot_outliers(
            pd.read_csv(args.outliers),
            args.output_dir,
        )

    print("\n=== FINAL FIGURES COMPLETE ===")
    print(f"Outputs: {args.output_dir}")


if __name__ == "__main__":
    main()
