#!/usr/bin/env python3
"""
Full BERTopic candidate model (Run B) for the GeNeDis 2026 glioblastoma study.

Candidate clustering:
    UMAP:
        n_neighbors=15
        n_components=5
        min_dist=0.0
        metric="cosine"
        random_state=42

    HDBSCAN:
        min_cluster_size=100
        min_samples=2
        metric="euclidean"
        cluster_selection_method="eom"
        prediction_data=True

Topic representation:
    CountVectorizer(
        stop_words="english",
        ngram_range=(1, 2),
        max_features=30000
    )
    ClassTfidfTransformer(reduce_frequent_words=True)

The script uses the PRECOMPUTED BiomedBERT embeddings and:
  - aligns documents strictly through embedding_row <-> PMID mapping
  - fits the full BERTopic model
  - compares raw HDBSCAN labels to the previous tuning Run B, if available
  - saves document-topic assignments
  - saves topic info, top terms, representative papers
  - calculates outlier rate by publication year
  - calculates topic prevalence by publication year
  - saves UMAP coordinates and clustering submodels
  - saves a full BERTopic pickle for exact local reuse

IMPORTANT:
  The full pickle should only be loaded if you trust the file/source.

Example:
    python run_bertopic_candidate_b.py \
        --corpus glioblastoma_analysis_corpus_with_tokens.parquet \
        --embeddings embeddings/biomedbert_embeddings.npy \
        --mapping embeddings/embedding_pmids.parquet \
        --reference-clusters hdbscan_tuning/runs/B/assignments.parquet \
        --output-dir bertopic_candidate_b

Background:
    nohup python -u run_bertopic_candidate_b.py \
        --corpus glioblastoma_analysis_corpus_with_tokens.parquet \
        --embeddings embeddings/biomedbert_embeddings.npy \
        --mapping embeddings/embedding_pmids.parquet \
        --reference-clusters hdbscan_tuning/runs/B/assignments.parquet \
        --output-dir bertopic_candidate_b \
        > bertopic_candidate_b.log 2>&1 &
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from bertopic import BERTopic
from bertopic.vectorizers import ClassTfidfTransformer
from hdbscan import HDBSCAN
from sklearn.feature_extraction.text import CountVectorizer
from umap import UMAP


def parse_args():
    p = argparse.ArgumentParser(
        description="Run full BERTopic candidate B."
    )

    p.add_argument(
        "--corpus",
        type=Path,
        default=Path("glioblastoma_analysis_corpus_with_tokens.parquet"),
    )
    p.add_argument(
        "--embeddings",
        type=Path,
        default=Path("embeddings/biomedbert_embeddings.npy"),
    )
    p.add_argument(
        "--mapping",
        type=Path,
        default=Path("embeddings/embedding_pmids.parquet"),
    )
    p.add_argument(
        "--reference-clusters",
        type=Path,
        default=Path("hdbscan_tuning/runs/B/assignments.parquet"),
        help=(
            "Optional Run-B assignments from tune_hdbscan.py. "
            "If the file exists, raw HDBSCAN labels are compared."
        ),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("bertopic_candidate_b"),
    )

    # UMAP
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--n-components", type=int, default=5)
    p.add_argument("--min-dist", type=float, default=0.0)
    p.add_argument("--umap-metric", default="cosine")
    p.add_argument("--random-state", type=int, default=42)

    # HDBSCAN candidate B
    p.add_argument("--min-cluster-size", type=int, default=100)
    p.add_argument("--min-samples", type=int, default=2)
    p.add_argument("--hdbscan-metric", default="euclidean")
    p.add_argument(
        "--cluster-selection-method",
        choices=["eom", "leaf"],
        default="eom",
    )

    # Topic representation
    p.add_argument("--top-n-words", type=int, default=15)
    p.add_argument("--ngram-min", type=int, default=1)
    p.add_argument("--ngram-max", type=int, default=2)
    p.add_argument("--max-features", type=int, default=30000)
    p.add_argument("--top-documents-per-topic", type=int, default=10)

    p.add_argument("--overwrite", action="store_true")

    return p.parse_args()


def version_of(package):
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def write_json(path, obj):
    path.write_text(
        json.dumps(obj, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def save_environment(out):
    versions = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "bertopic": version_of("bertopic"),
        "umap-learn": version_of("umap-learn"),
        "hdbscan": version_of("hdbscan"),
        "scikit-learn": version_of("scikit-learn"),
        "numpy": version_of("numpy"),
        "pandas": version_of("pandas"),
    }

    write_json(out / "versions.json", versions)

    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            check=True,
            capture_output=True,
            text=True,
        )
        (out / "pip_freeze.txt").write_text(
            result.stdout,
            encoding="utf-8",
        )
    except Exception as exc:
        print(f"Warning: could not save pip freeze: {exc}")

    return versions


def load_aligned_data(corpus_path, embeddings_path, mapping_path):
    print("Loading embeddings...")
    embeddings = np.load(embeddings_path, mmap_mode="r")

    if embeddings.ndim != 2:
        raise ValueError(f"Expected 2D embeddings, got {embeddings.shape}")

    print(f"Embeddings shape: {embeddings.shape}")
    print(f"Embeddings dtype: {embeddings.dtype}")

    print("Loading mapping...")
    mapping = pd.read_parquet(mapping_path).copy()

    required = {"embedding_row", "pmid"}
    missing = required - set(mapping.columns)
    if missing:
        raise ValueError(f"Mapping missing columns: {sorted(missing)}")

    mapping["pmid"] = mapping["pmid"].astype(str).str.strip()
    mapping = mapping.sort_values("embedding_row").reset_index(drop=True)

    expected = np.arange(len(mapping), dtype=np.int64)
    actual = mapping["embedding_row"].to_numpy(dtype=np.int64)

    if not np.array_equal(expected, actual):
        raise ValueError("embedding_row is not exactly 0..N-1.")
    if mapping["pmid"].duplicated().any():
        raise ValueError("Duplicate PMIDs in embedding mapping.")
    if len(mapping) != embeddings.shape[0]:
        raise ValueError(
            f"Mapping rows {len(mapping):,} != embedding rows {embeddings.shape[0]:,}"
        )

    print("Loading corpus...")
    corpus = pd.read_parquet(corpus_path).copy()

    required = {"pmid", "document", "publication_year"}
    missing = required - set(corpus.columns)
    if missing:
        raise ValueError(f"Corpus missing columns: {sorted(missing)}")

    corpus["pmid"] = corpus["pmid"].astype(str).str.strip()

    if corpus["pmid"].duplicated().any():
        raise ValueError("Duplicate PMIDs in analytical corpus.")

    # Mapping defines exact embedding row order and also excludes PMID 27454254.
    aligned = mapping[["embedding_row", "pmid"]].merge(
        corpus,
        on="pmid",
        how="left",
        validate="one_to_one",
        indicator=True,
    )

    missing_docs = aligned["_merge"].ne("both")
    if missing_docs.any():
        examples = aligned.loc[missing_docs, "pmid"].tolist()[:10]
        raise ValueError(
            f"{missing_docs.sum():,} PMIDs missing from corpus. Examples: {examples}"
        )

    aligned = aligned.drop(columns="_merge")

    aligned["document"] = (
        aligned["document"].fillna("").astype(str).str.strip()
    )

    if aligned["document"].eq("").any():
        raise ValueError("Blank aligned documents found.")

    aligned["publication_year"] = pd.to_numeric(
        aligned["publication_year"],
        errors="raise",
    ).astype(int)

    norms = np.linalg.norm(embeddings, axis=1)
    print(
        "Alignment QC passed: "
        f"{len(aligned):,} documents / embeddings / PMIDs"
    )
    print(
        "Embedding norm QC: "
        f"min={norms.min():.6f}, "
        f"mean={norms.mean():.6f}, "
        f"max={norms.max():.6f}"
    )

    return aligned, embeddings


def compare_reference_labels(
    raw_labels,
    mapping,
    reference_path,
):
    result = {
        "reference_available": False,
        "n_compared": 0,
        "exact_label_matches": None,
        "exact_label_match_percentage": None,
        "same_outlier_status": None,
        "same_outlier_status_percentage": None,
    }

    if reference_path is None or not reference_path.exists():
        print(
            "Reference Run-B assignments not found; "
            "skipping clustering reproducibility comparison."
        )
        return result

    print(f"Comparing to tuning reference: {reference_path}")

    ref = pd.read_parquet(reference_path).copy()

    required = {"embedding_row", "pmid", "cluster"}
    missing = required - set(ref.columns)
    if missing:
        raise ValueError(
            f"Reference cluster file missing: {sorted(missing)}"
        )

    ref["pmid"] = ref["pmid"].astype(str).str.strip()
    ref = ref.sort_values("embedding_row").reset_index(drop=True)

    current_pmids = mapping["pmid"].astype(str).to_numpy()
    ref_pmids = ref["pmid"].to_numpy()

    if len(ref) != len(mapping):
        raise ValueError(
            f"Reference rows {len(ref):,} != current rows {len(mapping):,}"
        )

    if not np.array_equal(current_pmids, ref_pmids):
        raise ValueError(
            "PMID order differs between current model and tuning reference."
        )

    ref_labels = ref["cluster"].to_numpy(dtype=np.int32)
    raw_labels = np.asarray(raw_labels, dtype=np.int32)

    exact = raw_labels == ref_labels
    same_outlier = (raw_labels == -1) == (ref_labels == -1)

    result = {
        "reference_available": True,
        "n_compared": int(len(ref)),
        "exact_label_matches": int(exact.sum()),
        "exact_label_match_percentage": float(100 * exact.mean()),
        "same_outlier_status": int(same_outlier.sum()),
        "same_outlier_status_percentage": float(100 * same_outlier.mean()),
    }

    print(
        f"Raw HDBSCAN exact-label match: "
        f"{result['exact_label_match_percentage']:.2f}%"
    )
    print(
        f"Same outlier/non-outlier status: "
        f"{result['same_outlier_status_percentage']:.2f}%"
    )

    return result


def topic_terms_table(topic_model, topic_ids, top_n):
    rows = []

    for topic_id in topic_ids:
        terms = topic_model.get_topic(int(topic_id))
        if not terms:
            continue

        for rank, (term, score) in enumerate(terms[:top_n], start=1):
            rows.append(
                {
                    "topic": int(topic_id),
                    "rank": rank,
                    "term": term,
                    "ctfidf_score": float(score),
                }
            )

    return pd.DataFrame(rows)


def top_documents_table(assignments, n_per_topic):
    rows = []

    valid = assignments[assignments["topic"] != -1].copy()

    for topic_id, group in valid.groupby("topic", sort=True):
        group = group.sort_values(
            ["cluster_membership_strength", "publication_year"],
            ascending=[False, False],
            na_position="last",
        ).head(n_per_topic)

        for rank, (_, row) in enumerate(group.iterrows(), start=1):
            rows.append(
                {
                    "topic": int(topic_id),
                    "rank": rank,
                    "pmid": row["pmid"],
                    "publication_year": int(row["publication_year"]),
                    "cluster_membership_strength": float(
                        row["cluster_membership_strength"]
                    ),
                    "title": row.get("title", ""),
                }
            )

    return pd.DataFrame(rows)


def outliers_by_year_table(assignments):
    tmp = assignments.copy()
    tmp["is_outlier"] = tmp["topic"].eq(-1)

    result = (
        tmp.groupby("publication_year", as_index=False)
        .agg(
            n_documents=("pmid", "size"),
            n_outliers=("is_outlier", "sum"),
        )
        .sort_values("publication_year")
    )

    result["n_clustered"] = (
        result["n_documents"] - result["n_outliers"]
    )

    result["outlier_percentage"] = (
        100.0 * result["n_outliers"] / result["n_documents"]
    )

    return result


def topic_prevalence_by_year_table(assignments):
    """
    Prevalence uses ALL GBM publications in each year as denominator.
    This avoids confusing general publication growth with topic growth.
    Outlier topic -1 is excluded from the output topics but remains in
    the yearly denominator.
    """
    year_totals = (
        assignments.groupby("publication_year")
        .size()
        .rename("all_gbm_documents")
    )

    topic_year = (
        assignments[assignments["topic"] != -1]
        .groupby(["publication_year", "topic"])
        .size()
        .rename("topic_documents")
        .reset_index()
    )

    topic_year = topic_year.merge(
        year_totals.reset_index(),
        on="publication_year",
        how="left",
        validate="many_to_one",
    )

    topic_year["prevalence_all_gbm"] = (
        topic_year["topic_documents"]
        / topic_year["all_gbm_documents"]
    )

    topic_year["prevalence_percent_all_gbm"] = (
        100.0 * topic_year["prevalence_all_gbm"]
    )

    # Also provide prevalence among clustered documents as a secondary metric.
    clustered_totals = (
        assignments[assignments["topic"] != -1]
        .groupby("publication_year")
        .size()
        .rename("clustered_documents")
    )

    topic_year = topic_year.merge(
        clustered_totals.reset_index(),
        on="publication_year",
        how="left",
        validate="many_to_one",
    )

    topic_year["prevalence_clustered_only"] = (
        topic_year["topic_documents"]
        / topic_year["clustered_documents"]
    )

    topic_year["prevalence_percent_clustered_only"] = (
        100.0 * topic_year["prevalence_clustered_only"]
    )

    return topic_year.sort_values(
        ["publication_year", "topic"]
    ).reset_index(drop=True)


def main():
    args = parse_args()

    out = args.output_dir

    if out.exists() and any(out.iterdir()) and not args.overwrite:
        raise RuntimeError(
            f"Output directory '{out}' is not empty. "
            "Use a new directory or --overwrite."
        )

    out.mkdir(parents=True, exist_ok=True)

    complete_path = out / "COMPLETE"
    if complete_path.exists():
        complete_path.unlink()

    versions = save_environment(out)

    aligned, embeddings = load_aligned_data(
        args.corpus,
        args.embeddings,
        args.mapping,
    )

    docs = aligned["document"].tolist()
    n_docs = len(docs)

    config = {
        "candidate": "B",
        "n_documents": int(n_docs),
        "embedding_shape": [int(x) for x in embeddings.shape],
        "umap": {
            "n_neighbors": int(args.n_neighbors),
            "n_components": int(args.n_components),
            "min_dist": float(args.min_dist),
            "metric": args.umap_metric,
            "random_state": int(args.random_state),
        },
        "hdbscan": {
            "min_cluster_size": int(args.min_cluster_size),
            "min_samples": int(args.min_samples),
            "metric": args.hdbscan_metric,
            "cluster_selection_method": args.cluster_selection_method,
            "prediction_data": True,
        },
        "vectorizer": {
            "stop_words": "english",
            "ngram_range": [int(args.ngram_min), int(args.ngram_max)],
            "max_features": int(args.max_features),
        },
        "ctfidf": {
            "reduce_frequent_words": True,
            "bm25_weighting": False,
        },
        "bertopic": {
            "top_n_words": int(args.top_n_words),
            "calculate_probabilities": False,
            "nr_topics": None,
        },
        "versions": versions,
    }

    write_json(out / "candidate_b_config.json", config)

    umap_model = UMAP(
        n_neighbors=args.n_neighbors,
        n_components=args.n_components,
        min_dist=args.min_dist,
        metric=args.umap_metric,
        random_state=args.random_state,
        transform_seed=args.random_state,
        low_memory=True,
        verbose=True,
    )

    hdbscan_model = HDBSCAN(
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        metric=args.hdbscan_metric,
        cluster_selection_method=args.cluster_selection_method,
        prediction_data=True,
        core_dist_n_jobs=-1,
    )

    vectorizer_model = CountVectorizer(
        stop_words="english",
        ngram_range=(args.ngram_min, args.ngram_max),
        max_features=args.max_features,
        lowercase=True,
    )

    ctfidf_model = ClassTfidfTransformer(
        bm25_weighting=False,
        reduce_frequent_words=True,
    )

    topic_model = BERTopic(
        embedding_model=None,
        umap_model=umap_model,
        hdbscan_model=hdbscan_model,
        vectorizer_model=vectorizer_model,
        ctfidf_model=ctfidf_model,
        top_n_words=args.top_n_words,
        nr_topics=None,
        calculate_probabilities=False,
        verbose=True,
    )

    print("\n=== FITTING FULL BERTOPIC CANDIDATE B ===")
    start = time.time()

    topics, probabilities = topic_model.fit_transform(
        docs,
        embeddings,
    )

    elapsed = time.time() - start

    topics = np.asarray(topics, dtype=np.int32)

    raw_labels = np.asarray(
        hdbscan_model.labels_,
        dtype=np.int32,
    )

    strengths = np.asarray(
        hdbscan_model.probabilities_,
        dtype=np.float32,
    )

    outlier_scores = np.asarray(
        hdbscan_model.outlier_scores_,
        dtype=np.float32,
    )

    if not (
        len(topics)
        == len(raw_labels)
        == len(strengths)
        == len(outlier_scores)
        == n_docs
    ):
        raise RuntimeError("Output length mismatch.")

    # Save exact UMAP space from this full candidate run.
    reduced = np.asarray(
        umap_model.embedding_,
        dtype=np.float32,
    )
    np.save(out / "umap_5d.npy", reduced)

    # Compare raw HDBSCAN clustering with previous tuning Run B.
    reproducibility = compare_reference_labels(
        raw_labels=raw_labels,
        mapping=aligned[["embedding_row", "pmid"]],
        reference_path=args.reference_clusters,
    )
    write_json(
        out / "clustering_reproducibility.json",
        reproducibility,
    )

    # ------------------------------------------------------
    # Document-level assignments
    # ------------------------------------------------------
    keep_cols = [
        c for c in [
            "embedding_row",
            "pmid",
            "publication_year",
            "title",
            "doi",
            "publication_types_json",
            "mesh_descriptors_json",
            "mesh_terms_json",
            "keywords_json",
        ]
        if c in aligned.columns
    ]

    assignments = aligned[keep_cols].copy()

    # BERTopic may remap raw HDBSCAN cluster IDs by topic frequency,
    # so both are retained.
    assignments["topic"] = topics
    assignments["raw_hdbscan_cluster"] = raw_labels
    assignments["cluster_membership_strength"] = strengths
    assignments["hdbscan_outlier_score"] = outlier_scores

    assignments.to_parquet(
        out / "topic_assignments.parquet",
        index=False,
    )
    assignments.to_csv(
        out / "topic_assignments.csv",
        index=False,
    )

    # ------------------------------------------------------
    # Topic information
    # ------------------------------------------------------
    topic_info = topic_model.get_topic_info().copy()

    topic_info.to_csv(
        out / "topic_info.csv",
        index=False,
    )

    topic_ids = sorted(
        int(x) for x in np.unique(topics)
    )

    terms = topic_terms_table(
        topic_model,
        topic_ids,
        args.top_n_words,
    )
    terms.to_csv(
        out / "topic_terms.csv",
        index=False,
    )

    top_docs = top_documents_table(
        assignments,
        args.top_documents_per_topic,
    )
    top_docs.to_csv(
        out / "topic_top_documents.csv",
        index=False,
    )

    # BERTopic's own representative docs.
    representative_rows = []

    for topic_id in topic_ids:
        if topic_id == -1:
            continue

        try:
            reps = topic_model.get_representative_docs(topic_id) or []
        except Exception:
            reps = []

        for rank, text in enumerate(reps, start=1):
            representative_rows.append(
                {
                    "topic": int(topic_id),
                    "rank": int(rank),
                    "representative_document": text,
                }
            )

    pd.DataFrame(
        representative_rows
    ).to_csv(
        out / "bertopic_representative_documents.csv",
        index=False,
    )

    # ------------------------------------------------------
    # Temporal diagnostics
    # ------------------------------------------------------
    outliers_year = outliers_by_year_table(assignments)

    outliers_year.to_csv(
        out / "outliers_by_year.csv",
        index=False,
    )

    prevalence = topic_prevalence_by_year_table(assignments)

    prevalence.to_csv(
        out / "topic_prevalence_by_year.csv",
        index=False,
    )

    # ------------------------------------------------------
    # Summary
    # ------------------------------------------------------
    n_outliers = int((topics == -1).sum())
    n_topics = len([x for x in topic_ids if x != -1])

    topic_sizes = (
        assignments[assignments["topic"] != -1]
        .groupby("topic")
        .size()
        .sort_values(ascending=False)
    )

    summary = {
        "candidate": "B",
        "n_documents": int(n_docs),
        "n_topics_excluding_outlier": int(n_topics),
        "n_outliers": int(n_outliers),
        "outlier_percentage": float(100 * n_outliers / n_docs),
        "largest_topic_size": (
            int(topic_sizes.max()) if len(topic_sizes) else None
        ),
        "median_topic_size": (
            float(topic_sizes.median()) if len(topic_sizes) else None
        ),
        "smallest_topic_size": (
            int(topic_sizes.min()) if len(topic_sizes) else None
        ),
        "membership_mean_clustered": float(
            strengths[topics != -1].mean()
        ),
        "membership_median_clustered": float(
            np.median(strengths[topics != -1])
        ),
        "fit_elapsed_seconds": float(elapsed),
        "publication_year_min": int(
            assignments["publication_year"].min()
        ),
        "publication_year_max": int(
            assignments["publication_year"].max()
        ),
        "reproducibility_vs_tuning_B": reproducibility,
    }

    write_json(
        out / "candidate_b_summary.json",
        summary,
    )

    # ------------------------------------------------------
    # Save models
    # ------------------------------------------------------
    submodels = out / "submodels"
    submodels.mkdir(exist_ok=True)

    joblib.dump(
        umap_model,
        submodels / "umap_model.joblib",
    )
    joblib.dump(
        hdbscan_model,
        submodels / "hdbscan_model.joblib",
    )
    joblib.dump(
        vectorizer_model,
        submodels / "vectorizer_model.joblib",
    )
    joblib.dump(
        ctfidf_model,
        submodels / "ctfidf_model.joblib",
    )

    # Full local model serialization. Load only if you trust the file.
    print("Saving full BERTopic pickle...")
    topic_model.save(
        str(out / "bertopic_candidate_b.pkl"),
        serialization="pickle",
    )

    complete_path.write_text(
        "Full BERTopic candidate B completed successfully.\n",
        encoding="utf-8",
    )

    # ------------------------------------------------------
    # Final console report
    # ------------------------------------------------------
    print("\n=== BERTOPIC CANDIDATE B COMPLETE ===")
    print(f"Documents:              {n_docs:,}")
    print(f"Topics (excluding -1):  {n_topics:,}")
    print(f"Outliers:               {n_outliers:,}")
    print(
        f"Outlier percentage:     "
        f"{summary['outlier_percentage']:.2f}%"
    )
    print(f"Largest topic:          {summary['largest_topic_size']:,}")
    print(f"Median topic size:      {summary['median_topic_size']:.1f}")
    print(f"Smallest topic:         {summary['smallest_topic_size']:,}")
    print(
        f"Membership mean:        "
        f"{summary['membership_mean_clustered']:.4f}"
    )
    print(
        f"Membership median:      "
        f"{summary['membership_median_clustered']:.4f}"
    )

    if reproducibility["reference_available"]:
        print(
            f"Run-B label match:      "
            f"{reproducibility['exact_label_match_percentage']:.2f}%"
        )
        print(
            f"Run-B outlier match:    "
            f"{reproducibility['same_outlier_status_percentage']:.2f}%"
        )

    print("\nTop 20 topics by size:")
    display_cols = [
        c for c in ["Topic", "Count", "Name"]
        if c in topic_info.columns
    ]

    print(
        topic_info[topic_info["Topic"] != -1]
        .head(20)[display_cols]
        .to_string(index=False)
    )

    print("\nOutlier rate by year:")
    print(
        outliers_year.to_string(
            index=False,
            float_format=lambda x: f"{x:.2f}",
        )
    )

    print(f"\nOutputs saved to: {out}")


if __name__ == "__main__":
    main()
