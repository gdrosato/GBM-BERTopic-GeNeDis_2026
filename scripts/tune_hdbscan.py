#!/usr/bin/env python3
"""
Fast HDBSCAN tuning on an already-computed 5D UMAP representation.

GeNeDis 2026 — Glioblastoma literature mapping

This script DOES NOT rerun:
  - BiomedBERT embeddings
  - UMAP
  - BERTopic / c-TF-IDF

It only reruns HDBSCAN on the saved UMAP coordinates.

Default configurations
----------------------
Baseline: min_cluster_size=100, min_samples=10
A:        min_cluster_size=100, min_samples=5
B:        min_cluster_size=100, min_samples=2
C:        min_cluster_size=75,  min_samples=5
D:        min_cluster_size=50,  min_samples=5

Outputs
-------
hdbscan_tuning/
├── hdbscan_tuning_summary.csv
├── hdbscan_tuning_summary.json
├── hdbscan_tuning_comparison_sorted.csv
├── config.json
├── runs/
│   ├── baseline/
│   │   ├── assignments.parquet
│   │   ├── cluster_sizes.csv
│   │   ├── cluster_persistence.csv
│   │   ├── summary.json
│   │   ├── config.json
│   │   └── hdbscan_model.joblib
│   ├── A/
│   ├── B/
│   ├── C/
│   └── D/
└── COMPLETE

Install
-------
pip install hdbscan pandas pyarrow numpy joblib

Run
---
python tune_hdbscan.py \
    --umap bertopic_baseline/umap_5d.npy \
    --mapping embeddings/embedding_pmids.parquet \
    --output-dir hdbscan_tuning

Background
----------
nohup python -u tune_hdbscan.py \
    --umap bertopic_baseline/umap_5d.npy \
    --mapping embeddings/embedding_pmids.parquet \
    --output-dir hdbscan_tuning \
    > hdbscan_tuning.log 2>&1 &
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from hdbscan import HDBSCAN


DEFAULT_CONFIGS = [
    {"name": "baseline", "min_cluster_size": 100, "min_samples": 10},
    {"name": "A", "min_cluster_size": 100, "min_samples": 5},
    {"name": "B", "min_cluster_size": 100, "min_samples": 2},
    {"name": "C", "min_cluster_size": 75, "min_samples": 5},
    {"name": "D", "min_cluster_size": 50, "min_samples": 5},
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Tune HDBSCAN on a fixed UMAP representation."
    )
    p.add_argument(
        "--umap",
        type=Path,
        default=Path("bertopic_baseline/umap_5d.npy"),
    )
    p.add_argument(
        "--mapping",
        type=Path,
        default=Path("embeddings/embedding_pmids.parquet"),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("hdbscan_tuning"),
    )
    p.add_argument("--metric", default="euclidean")
    p.add_argument(
        "--cluster-selection-method",
        choices=["eom", "leaf"],
        default="eom",
    )
    p.add_argument(
        "--near-min-factor",
        type=float,
        default=1.25,
        help=(
            "A cluster is considered near the minimum if "
            "size <= min_cluster_size * factor."
        ),
    )
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def write_json(path: Path, obj) -> None:
    path.write_text(
        json.dumps(obj, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def load_inputs(umap_path: Path, mapping_path: Path):
    if not umap_path.exists():
        raise FileNotFoundError(umap_path)
    if not mapping_path.exists():
        raise FileNotFoundError(mapping_path)

    print(f"Loading UMAP coordinates: {umap_path}")
    reduced = np.load(umap_path, mmap_mode="r")

    if reduced.ndim != 2:
        raise ValueError(f"Expected 2D UMAP matrix, got {reduced.shape}")
    if not np.isfinite(reduced).all():
        raise ValueError("UMAP coordinates contain NaN or infinity.")

    print(f"UMAP shape: {reduced.shape}, dtype={reduced.dtype}")

    print(f"Loading mapping: {mapping_path}")
    mapping = pd.read_parquet(mapping_path).copy()

    required = {"embedding_row", "pmid"}
    missing = required - set(mapping.columns)
    if missing:
        raise ValueError(f"Mapping missing: {sorted(missing)}")

    mapping["pmid"] = mapping["pmid"].astype(str).str.strip()
    mapping = mapping.sort_values("embedding_row").reset_index(drop=True)

    expected = np.arange(len(mapping), dtype=np.int64)
    actual = mapping["embedding_row"].to_numpy(dtype=np.int64)

    if not np.array_equal(expected, actual):
        raise ValueError("embedding_row is not exactly 0..N-1.")
    if mapping["pmid"].duplicated().any():
        raise ValueError("Duplicate PMIDs found in mapping.")
    if len(mapping) != reduced.shape[0]:
        raise ValueError(
            f"Mapping rows ({len(mapping):,}) != UMAP rows ({reduced.shape[0]:,})."
        )

    print(f"Alignment QC passed: {len(mapping):,} rows / PMIDs")
    return reduced, mapping


def cluster_size_table(labels: np.ndarray) -> pd.DataFrame:
    valid = labels[labels != -1]
    if len(valid) == 0:
        return pd.DataFrame(columns=["cluster", "size"])

    clusters, counts = np.unique(valid, return_counts=True)
    return (
        pd.DataFrame(
            {
                "cluster": clusters.astype(int),
                "size": counts.astype(int),
            }
        )
        .sort_values("size", ascending=False)
        .reset_index(drop=True)
    )


def persistence_table(clusterer: HDBSCAN, size_df: pd.DataFrame) -> pd.DataFrame:
    persistence = np.asarray(clusterer.cluster_persistence_, dtype=float)
    p = pd.DataFrame(
        {
            "cluster": np.arange(len(persistence), dtype=int),
            "persistence": persistence,
        }
    )
    return (
        size_df.merge(p, on="cluster", how="left", validate="one_to_one")
        .sort_values("size", ascending=False)
        .reset_index(drop=True)
    )


def safe_stats(values: np.ndarray):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"mean": None, "median": None, "min": None, "max": None}
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }


def run_configuration(
    reduced,
    mapping: pd.DataFrame,
    cfg: dict,
    output_dir: Path,
    metric: str,
    selection_method: str,
    near_min_factor: float,
):
    name = cfg["name"]
    min_cluster_size = int(cfg["min_cluster_size"])
    min_samples = int(cfg["min_samples"])

    run_dir = output_dir / "runs" / name
    run_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 72)
    print(
        f"RUN {name}: min_cluster_size={min_cluster_size}, "
        f"min_samples={min_samples}"
    )
    print("=" * 72)

    clusterer = HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        metric=metric,
        cluster_selection_method=selection_method,
        prediction_data=True,
        core_dist_n_jobs=-1,
    )

    start = time.time()
    labels = clusterer.fit_predict(reduced).astype(np.int32)
    elapsed = time.time() - start

    probabilities = np.asarray(clusterer.probabilities_, dtype=np.float32)
    outlier_scores = np.asarray(clusterer.outlier_scores_, dtype=np.float32)

    if not (
        len(labels)
        == len(probabilities)
        == len(outlier_scores)
        == len(mapping)
    ):
        raise RuntimeError(f"Length mismatch in run {name}.")

    size_df = cluster_size_table(labels)
    persistence_df = persistence_table(clusterer, size_df)

    n_docs = len(labels)
    n_outliers = int(np.sum(labels == -1))
    n_clustered = n_docs - n_outliers
    n_clusters = len(size_df)

    near_threshold = int(np.floor(min_cluster_size * near_min_factor))

    if n_clusters:
        cluster_sizes = size_df["size"].to_numpy()
        n_near_min = int(np.sum(cluster_sizes <= near_threshold))
        largest_size = int(np.max(cluster_sizes))
        smallest_size = int(np.min(cluster_sizes))
        median_size = float(np.median(cluster_sizes))
        mean_size = float(np.mean(cluster_sizes))
    else:
        n_near_min = 0
        largest_size = smallest_size = median_size = mean_size = None

    clustered_mask = labels != -1
    noise_mask = labels == -1

    membership_clustered = safe_stats(probabilities[clustered_mask])
    membership_all = safe_stats(probabilities)
    persistence_stats = safe_stats(clusterer.cluster_persistence_)
    outlier_clustered = safe_stats(outlier_scores[clustered_mask])
    outlier_noise = safe_stats(outlier_scores[noise_mask])

    assignments = mapping.copy()
    assignments["cluster"] = labels
    assignments["membership_strength"] = probabilities
    assignments["outlier_score"] = outlier_scores
    assignments.to_parquet(run_dir / "assignments.parquet", index=False)

    size_df.to_csv(run_dir / "cluster_sizes.csv", index=False)
    persistence_df.to_csv(run_dir / "cluster_persistence.csv", index=False)
    joblib.dump(clusterer, run_dir / "hdbscan_model.joblib")

    run_config = {
        "name": name,
        "min_cluster_size": min_cluster_size,
        "min_samples": min_samples,
        "metric": metric,
        "cluster_selection_method": selection_method,
        "prediction_data": True,
        "core_dist_n_jobs": -1,
        "near_min_factor": near_min_factor,
        "near_min_threshold": near_threshold,
    }
    write_json(run_dir / "config.json", run_config)

    summary = {
        "run": name,
        "min_cluster_size": min_cluster_size,
        "min_samples": min_samples,
        "n_documents": n_docs,
        "n_clusters": n_clusters,
        "n_clustered": n_clustered,
        "n_outliers": n_outliers,
        "outlier_percentage": 100.0 * n_outliers / n_docs,
        "largest_cluster_size": largest_size,
        "mean_cluster_size": mean_size,
        "median_cluster_size": median_size,
        "smallest_cluster_size": smallest_size,
        "near_min_threshold": near_threshold,
        "n_clusters_near_min_size": n_near_min,
        "pct_clusters_near_min_size": (
            100.0 * n_near_min / n_clusters if n_clusters else None
        ),
        "membership_mean_clustered": membership_clustered["mean"],
        "membership_median_clustered": membership_clustered["median"],
        "membership_min_clustered": membership_clustered["min"],
        "membership_mean_all": membership_all["mean"],
        "cluster_persistence_mean": persistence_stats["mean"],
        "cluster_persistence_median": persistence_stats["median"],
        "cluster_persistence_min": persistence_stats["min"],
        "cluster_persistence_max": persistence_stats["max"],
        "outlier_score_mean_clustered": outlier_clustered["mean"],
        "outlier_score_median_clustered": outlier_clustered["median"],
        "outlier_score_mean_noise": outlier_noise["mean"],
        "outlier_score_median_noise": outlier_noise["median"],
        "elapsed_seconds": elapsed,
    }

    write_json(run_dir / "summary.json", summary)

    print(f"Clusters:             {n_clusters:,}")
    print(
        f"Outliers:             {n_outliers:,} "
        f"({summary['outlier_percentage']:.2f}%)"
    )
    print(
        f"Cluster sizes:        min={smallest_size}, "
        f"median={median_size}, max={largest_size}"
    )
    print(
        f"Near-min clusters:    {n_near_min}/{n_clusters} "
        f"(<= {near_threshold})"
    )
    if membership_clustered["mean"] is not None:
        print(
            f"Membership strength:  mean={membership_clustered['mean']:.4f}, "
            f"median={membership_clustered['median']:.4f}"
        )
    if persistence_stats["mean"] is not None:
        print(
            f"Cluster persistence:  mean={persistence_stats['mean']:.4f}, "
            f"median={persistence_stats['median']:.4f}"
        )
    print(f"Elapsed:              {elapsed:.2f}s")

    return summary


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir

    if output_dir.exists() and any(output_dir.iterdir()) and not args.overwrite:
        raise RuntimeError(
            f"Output directory '{output_dir}' is not empty. "
            "Use a new directory or pass --overwrite."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    complete_path = output_dir / "COMPLETE"
    if complete_path.exists():
        complete_path.unlink()

    reduced, mapping = load_inputs(args.umap, args.mapping)

    master_config = {
        "umap_file": str(args.umap),
        "mapping_file": str(args.mapping),
        "n_documents": int(len(mapping)),
        "umap_shape": list(reduced.shape),
        "metric": args.metric,
        "cluster_selection_method": args.cluster_selection_method,
        "near_min_factor": args.near_min_factor,
        "configurations": DEFAULT_CONFIGS,
    }
    write_json(output_dir / "config.json", master_config)

    summaries = []
    total_start = time.time()

    for cfg in DEFAULT_CONFIGS:
        summaries.append(
            run_configuration(
                reduced=reduced,
                mapping=mapping,
                cfg=cfg,
                output_dir=output_dir,
                metric=args.metric,
                selection_method=args.cluster_selection_method,
                near_min_factor=args.near_min_factor,
            )
        )

    total_elapsed = time.time() - total_start
    summary_df = pd.DataFrame(summaries)

    summary_df.to_csv(
        output_dir / "hdbscan_tuning_summary.csv",
        index=False,
    )
    write_json(
        output_dir / "hdbscan_tuning_summary.json",
        summaries,
    )

    comparison_cols = [
        "run",
        "min_cluster_size",
        "min_samples",
        "n_clusters",
        "n_outliers",
        "outlier_percentage",
        "median_cluster_size",
        "pct_clusters_near_min_size",
        "membership_mean_clustered",
        "membership_median_clustered",
        "cluster_persistence_mean",
        "cluster_persistence_median",
        "elapsed_seconds",
    ]

    comparison = (
        summary_df[comparison_cols]
        .sort_values(
            ["outlier_percentage", "cluster_persistence_mean"],
            ascending=[True, False],
        )
        .reset_index(drop=True)
    )

    comparison.to_csv(
        output_dir / "hdbscan_tuning_comparison_sorted.csv",
        index=False,
    )

    complete_path.write_text(
        "HDBSCAN tuning completed successfully.\n",
        encoding="utf-8",
    )

    print("\n" + "=" * 72)
    print("HDBSCAN TUNING COMPLETE")
    print("=" * 72)

    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 220)

    print(
        summary_df[comparison_cols].to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
        )
    )

    print(f"\nTotal elapsed: {total_elapsed:.2f}s")
    print(f"Outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
