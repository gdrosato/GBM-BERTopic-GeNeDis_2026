#!/usr/bin/env python3
"""
Compare strict primary-label axis definitions against the original
primary+secondary-label axis definitions for the GeNeDis 2026 GBM analysis.

No embeddings or BERTopic model are rerun. The script operates only on:
  1) the curated expert-annotation workbook; and
  2) bertopic_candidate_b/topic_prevalence_by_year.csv

Primary analysis
----------------
Axis membership is determined ONLY from expert_primary_label:
  Hypoxia / angiogenesis          <- Hypoxia
  Metabolism / stress             <- Metabolism
  Immune / tumor microenvironment<- Immune microenvironment OR Tumor microenvironment
  Treatment resistance           <- Treatment resistance

Sensitivity analysis
--------------------
Axis membership is determined from expert_primary_label PLUS
expert_secondary_labels, reproducing the previous broader definition.

Outputs
-------
primary_axis_membership.csv
expanded_axis_membership.csv
primary_axis_prevalence_annual.csv
expanded_axis_prevalence_annual.csv
primary_axis_trend_stats.csv
expanded_axis_trend_stats.csv
axis_definition_comparison.csv
figure_temporal_axes_primary.pdf/png
figure_temporal_axes_expanded_sensitivity.pdf/png
figure_temporal_axis_definition_comparison.pdf/png
axis_definition_results.txt
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


AXIS_ORDER = [
    "Hypoxia / angiogenesis",
    "Metabolism / stress",
    "Immune / tumor microenvironment",
    "Treatment resistance",
]

PRIMARY_MAP = {
    "Hypoxia": ["Hypoxia / angiogenesis"],
    "Metabolism": ["Metabolism / stress"],
    "Immune microenvironment": ["Immune / tumor microenvironment"],
    "Tumor microenvironment": ["Immune / tumor microenvironment"],
    "Treatment resistance": ["Treatment resistance"],
}

# For the sensitivity analysis, the same biological label-to-axis mapping is
# applied to both the primary and secondary expert labels.
EXPANDED_LABEL_MAP = PRIMARY_MAP


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--validation-workbook", type=Path, required=True)
    p.add_argument("--prevalence", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, default=Path("axis_definition_sensitivity"))
    p.add_argument("--start-year", type=int, default=2011)
    p.add_argument("--end-year", type=int, default=2025)
    return p.parse_args()


def clean(value) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def split_secondary(value) -> list[str]:
    text = clean(value)
    if not text:
        return []
    return [x.strip() for x in text.split(";") if x.strip()]


def load_validation(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Expert Annotation")
    required = {
        "topic",
        "n_documents",
        "Name",
        "expert_include",
        "expert_primary_label",
        "expert_secondary_labels",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing validation columns: {sorted(missing)}")

    df = df.copy()
    df["topic"] = pd.to_numeric(df["topic"], errors="raise").astype(int)
    df["expert_include"] = df["expert_include"].map(clean)
    df["expert_primary_label"] = df["expert_primary_label"].map(clean)
    df["expert_secondary_labels"] = df["expert_secondary_labels"].map(clean)

    strict = df[df["expert_include"] == "Yes"].copy()

    if len(strict) != 41:
        print(f"WARNING: expected 41 strict-included topics, found {len(strict)}")

    return strict


def build_membership(strict: pd.DataFrame, expanded: bool) -> pd.DataFrame:
    rows = []

    for row in strict.itertuples(index=False):
        labels = [clean(row.expert_primary_label)]

        if expanded:
            labels.extend(split_secondary(row.expert_secondary_labels))

        axes = set()
        for label in labels:
            axes.update(EXPANDED_LABEL_MAP.get(label, []))

        for axis in AXIS_ORDER:
            if axis in axes:
                rows.append(
                    {
                        "topic": int(row.topic),
                        "n_documents": int(row.n_documents),
                        "bertopic_name": row.Name,
                        "expert_primary_label": row.expert_primary_label,
                        "expert_secondary_labels": row.expert_secondary_labels,
                        "axis": axis,
                        "definition": (
                            "primary_plus_secondary"
                            if expanded
                            else "primary_only"
                        ),
                    }
                )

    return pd.DataFrame(rows)


def load_prevalence(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    required = {
        "publication_year",
        "topic",
        "topic_documents",
        "all_gbm_documents",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing prevalence columns: {sorted(missing)}")

    df = df.copy()
    df["publication_year"] = pd.to_numeric(
        df["publication_year"], errors="raise"
    ).astype(int)
    df["topic"] = pd.to_numeric(df["topic"], errors="raise").astype(int)
    df["topic_documents"] = pd.to_numeric(
        df["topic_documents"], errors="raise"
    )
    df["all_gbm_documents"] = pd.to_numeric(
        df["all_gbm_documents"], errors="raise"
    )

    return df


def annual_axis_prevalence(
    prevalence: pd.DataFrame,
    membership: pd.DataFrame,
    start_year: int,
    end_year: int,
) -> pd.DataFrame:
    denominators = (
        prevalence[["publication_year", "all_gbm_documents"]]
        .drop_duplicates()
        .sort_values("publication_year")
    )

    # Make sure there is exactly one denominator per year.
    if denominators["publication_year"].duplicated().any():
        raise ValueError("Multiple annual GBM denominators found for a year.")

    year_grid = pd.DataFrame(
        {"publication_year": np.arange(start_year, end_year + 1, dtype=int)}
    )

    out = []

    for axis in AXIS_ORDER:
        topic_ids = (
            membership.loc[membership["axis"] == axis, "topic"]
            .astype(int)
            .drop_duplicates()
            .tolist()
        )

        tmp = prevalence[prevalence["topic"].isin(topic_ids)].copy()

        counts = (
            tmp.groupby("publication_year", as_index=False)["topic_documents"]
            .sum()
            .rename(columns={"topic_documents": "axis_topic_documents"})
        )

        annual = (
            year_grid
            .merge(counts, on="publication_year", how="left")
            .merge(
                denominators,
                on="publication_year",
                how="left",
                validate="one_to_one",
            )
        )

        annual["axis_topic_documents"] = annual["axis_topic_documents"].fillna(0)

        if annual["all_gbm_documents"].isna().any():
            years = annual.loc[
                annual["all_gbm_documents"].isna(), "publication_year"
            ].tolist()
            raise ValueError(f"Missing GBM denominators for years: {years}")

        annual["axis_prevalence_percent"] = (
            100.0
            * annual["axis_topic_documents"]
            / annual["all_gbm_documents"]
        )

        annual["rolling_3y_percent"] = (
            annual["axis_prevalence_percent"]
            .rolling(window=3, center=True, min_periods=1)
            .mean()
        )

        annual["axis"] = axis
        annual["n_axis_topics"] = len(topic_ids)
        annual["definition"] = membership["definition"].iloc[0]

        out.append(annual)

    return pd.concat(out, ignore_index=True)


def bh_fdr(p_values):
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)

    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return out


def trend_stats(annual: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for axis in AXIS_ORDER:
        g = (
            annual[annual["axis"] == axis]
            .sort_values("publication_year")
            .copy()
        )

        rho, p = spearmanr(
            g["publication_year"].to_numpy(),
            g["axis_prevalence_percent"].to_numpy(),
        )

        rows.append(
            {
                "axis": axis,
                "n_axis_topics": int(g["n_axis_topics"].iloc[0]),
                "n_years": len(g),
                "start_prevalence_percent": float(
                    g.iloc[0]["axis_prevalence_percent"]
                ),
                "end_prevalence_percent": float(
                    g.iloc[-1]["axis_prevalence_percent"]
                ),
                "spearman_rho": float(rho),
                "spearman_p": float(p),
            }
        )

    out = pd.DataFrame(rows)
    out["spearman_q_bh"] = bh_fdr(out["spearman_p"].to_numpy())
    return out


def direction(rho, q):
    if q < 0.05 and rho > 0:
        return "significant increase"
    if q < 0.05 and rho < 0:
        return "significant decrease"
    if rho > 0:
        return "positive, non-significant"
    if rho < 0:
        return "negative, non-significant"
    return "no monotonic trend"


def compare_stats(primary: pd.DataFrame, expanded: pd.DataFrame) -> pd.DataFrame:
    a = primary.add_prefix("primary_").rename(
        columns={"primary_axis": "axis"}
    )
    b = expanded.add_prefix("expanded_").rename(
        columns={"expanded_axis": "axis"}
    )

    c = a.merge(b, on="axis", how="inner")

    c["primary_direction"] = [
        direction(r, q)
        for r, q in zip(
            c["primary_spearman_rho"],
            c["primary_spearman_q_bh"],
        )
    ]
    c["expanded_direction"] = [
        direction(r, q)
        for r, q in zip(
            c["expanded_spearman_rho"],
            c["expanded_spearman_q_bh"],
        )
    ]

    c["same_rho_sign"] = (
        np.sign(c["primary_spearman_rho"])
        == np.sign(c["expanded_spearman_rho"])
    )

    c["same_significance_status"] = (
        (c["primary_spearman_q_bh"] < 0.05)
        == (c["expanded_spearman_q_bh"] < 0.05)
    )

    c["robust_direction_and_significance"] = (
        c["same_rho_sign"] & c["same_significance_status"]
    )

    return c


def plot_single(annual: pd.DataFrame, title: str, basename: Path):
    fig, ax = plt.subplots(figsize=(9.2, 5.4))

    for axis in AXIS_ORDER:
        g = (
            annual[annual["axis"] == axis]
            .sort_values("publication_year")
        )

        # Raw annual points: data used in statistical inference.
        ax.scatter(
            g["publication_year"],
            g["axis_prevalence_percent"],
            s=18,
            alpha=0.45,
        )

        # Rolling mean: visualization only.
        ax.plot(
            g["publication_year"],
            g["rolling_3y_percent"],
            linewidth=1.8,
            label=axis,
        )

    ax.set_title(title)
    ax.set_xlabel("Publication year")
    ax.set_ylabel("Share of all eligible GBM publications (%)")
    ax.set_xticks(
        np.arange(
            annual["publication_year"].min(),
            annual["publication_year"].max() + 1,
            2,
        )
    )
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(str(basename) + ".pdf", bbox_inches="tight")
    fig.savefig(str(basename) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_comparison(
    primary: pd.DataFrame,
    expanded: pd.DataFrame,
    basename: Path,
):
    # Four separate files are generated to avoid unreadable multi-panel figures.
    # This function creates one comparison PDF per axis.
    for axis in AXIS_ORDER:
        p = (
            primary[primary["axis"] == axis]
            .sort_values("publication_year")
        )
        e = (
            expanded[expanded["axis"] == axis]
            .sort_values("publication_year")
        )

        fig, ax = plt.subplots(figsize=(7.8, 4.8))

        ax.plot(
            p["publication_year"],
            p["rolling_3y_percent"],
            marker="o",
            markersize=3,
            linewidth=2,
            label="Primary-label only",
        )
        ax.plot(
            e["publication_year"],
            e["rolling_3y_percent"],
            marker="o",
            markersize=3,
            linewidth=1.5,
            linestyle="--",
            label="Primary + secondary labels",
        )

        ax.set_title(axis)
        ax.set_xlabel("Publication year")
        ax.set_ylabel("Share of all eligible GBM publications (%)")
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()

        safe = (
            axis.lower()
            .replace(" / ", "_")
            .replace(" ", "_")
            .replace("/", "_")
        )

        fig.savefig(
            str(basename) + f"_{safe}.pdf",
            bbox_inches="tight",
        )
        fig.savefig(
            str(basename) + f"_{safe}.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)


def write_summary(
    primary_stats: pd.DataFrame,
    expanded_stats: pd.DataFrame,
    comparison: pd.DataFrame,
    outpath: Path,
):
    lines = []

    lines.append("PRIMARY ANALYSIS: EXPERT PRIMARY LABELS ONLY")
    lines.append("=" * 52)

    for row in primary_stats.itertuples(index=False):
        lines.append(
            f"{row.axis}: topics={row.n_axis_topics}; "
            f"2011={row.start_prevalence_percent:.3f}%; "
            f"2025={row.end_prevalence_percent:.3f}%; "
            f"rho={row.spearman_rho:.4f}; "
            f"p={row.spearman_p:.5g}; "
            f"BH q={row.spearman_q_bh:.5g}; "
            f"{direction(row.spearman_rho, row.spearman_q_bh)}"
        )

    lines.append("")
    lines.append("SENSITIVITY ANALYSIS: PRIMARY + SECONDARY LABELS")
    lines.append("=" * 54)

    for row in expanded_stats.itertuples(index=False):
        lines.append(
            f"{row.axis}: topics={row.n_axis_topics}; "
            f"2011={row.start_prevalence_percent:.3f}%; "
            f"2025={row.end_prevalence_percent:.3f}%; "
            f"rho={row.spearman_rho:.4f}; "
            f"p={row.spearman_p:.5g}; "
            f"BH q={row.spearman_q_bh:.5g}; "
            f"{direction(row.spearman_rho, row.spearman_q_bh)}"
        )

    lines.append("")
    lines.append("ROBUSTNESS")
    lines.append("=" * 10)

    for row in comparison.itertuples(index=False):
        lines.append(
            f"{row.axis}: "
            f"primary=[{row.primary_direction}], "
            f"expanded=[{row.expanded_direction}], "
            f"same rho sign={row.same_rho_sign}, "
            f"same significance status={row.same_significance_status}"
        )

    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    strict = load_validation(args.validation_workbook)
    prevalence = load_prevalence(args.prevalence)

    primary_membership = build_membership(strict, expanded=False)
    expanded_membership = build_membership(strict, expanded=True)

    primary_membership.to_csv(
        args.output_dir / "primary_axis_membership.csv",
        index=False,
    )
    expanded_membership.to_csv(
        args.output_dir / "expanded_axis_membership.csv",
        index=False,
    )

    primary_annual = annual_axis_prevalence(
        prevalence,
        primary_membership,
        args.start_year,
        args.end_year,
    )
    expanded_annual = annual_axis_prevalence(
        prevalence,
        expanded_membership,
        args.start_year,
        args.end_year,
    )

    primary_annual.to_csv(
        args.output_dir / "primary_axis_prevalence_annual.csv",
        index=False,
    )
    expanded_annual.to_csv(
        args.output_dir / "expanded_axis_prevalence_annual.csv",
        index=False,
    )

    primary_stats = trend_stats(primary_annual)
    expanded_stats = trend_stats(expanded_annual)

    primary_stats.to_csv(
        args.output_dir / "primary_axis_trend_stats.csv",
        index=False,
    )
    expanded_stats.to_csv(
        args.output_dir / "expanded_axis_trend_stats.csv",
        index=False,
    )

    comparison = compare_stats(primary_stats, expanded_stats)
    comparison.to_csv(
        args.output_dir / "axis_definition_comparison.csv",
        index=False,
    )

    plot_single(
        primary_annual,
        "Normalized prevalence of primary-label biological axes, 2011–2025",
        args.output_dir / "figure_temporal_axes_primary",
    )

    plot_single(
        expanded_annual,
        "Sensitivity analysis: primary + secondary expert labels, 2011–2025",
        args.output_dir / "figure_temporal_axes_expanded_sensitivity",
    )

    plot_comparison(
        primary_annual,
        expanded_annual,
        args.output_dir / "figure_axis_definition_comparison",
    )

    write_summary(
        primary_stats,
        expanded_stats,
        comparison,
        args.output_dir / "axis_definition_results.txt",
    )

    print("\n=== AXIS-DEFINITION SENSITIVITY COMPLETE ===")
    print("\nPrimary-label-only trend statistics:")
    print(
        primary_stats.to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
        )
    )

    print("\nExpanded primary+secondary trend statistics:")
    print(
        expanded_stats.to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
        )
    )

    print("\nComparison:")
    cols = [
        "axis",
        "primary_n_axis_topics",
        "expanded_n_axis_topics",
        "primary_spearman_rho",
        "primary_spearman_q_bh",
        "expanded_spearman_rho",
        "expanded_spearman_q_bh",
        "robust_direction_and_significance",
    ]
    print(
        comparison[cols].to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
        )
    )

    print(f"\nOutputs: {args.output_dir}")


if __name__ == "__main__":
    main()
