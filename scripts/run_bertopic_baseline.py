#!/usr/bin/env python3

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
    p = argparse.ArgumentParser(description="Run BERTopic baseline on precomputed biomedical embeddings.")
    p.add_argument("--corpus", type=Path, default=Path("glioblastoma_analysis_corpus_with_tokens.parquet"))
    p.add_argument("--embeddings", type=Path, default=Path("embeddings/biomedbert_embeddings.npy"))
    p.add_argument("--mapping", type=Path, default=Path("embeddings/embedding_pmids.parquet"))
    p.add_argument("--output-dir", type=Path, default=Path("bertopic_baseline"))

    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--n-components", type=int, default=5)
    p.add_argument("--min-dist", type=float, default=0.0)
    p.add_argument("--umap-metric", default="cosine")
    p.add_argument("--random-state", type=int, default=42)

    p.add_argument("--min-cluster-size", type=int, default=100)
    p.add_argument("--min-samples", type=int, default=10)
    p.add_argument("--hdbscan-metric", default="euclidean")
    p.add_argument("--cluster-selection-method", choices=["eom", "leaf"], default="eom")

    p.add_argument("--top-n-words", type=int, default=15)
    p.add_argument("--ngram-min", type=int, default=1)
    p.add_argument("--ngram-max", type=int, default=2)
    p.add_argument("--max-features", type=int, default=30000, help="Use 0 for no vocabulary cap.")
    p.add_argument("--no-reduce-frequent-words", action="store_true")
    p.add_argument("--top-documents-per-topic", type=int, default=10)
    p.add_argument("--save-full-pickle", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def version_of(package):
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def save_environment(output_dir):
    versions = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "bertopic": version_of("bertopic"),
        "umap-learn": version_of("umap-learn"),
        "hdbscan": version_of("hdbscan"),
        "scikit-learn": version_of("scikit-learn"),
        "numpy": version_of("numpy"),
        "pandas": version_of("pandas"),
        "safetensors": version_of("safetensors"),
    }
    write_json(output_dir / "versions.json", versions)
    try:
        result = subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True, capture_output=True, text=True)
        (output_dir / "pip_freeze.txt").write_text(result.stdout, encoding="utf-8")
    except Exception as exc:
        print(f"Warning: could not save pip freeze: {exc}")
    return versions


def load_aligned_data(corpus_path, embeddings_path, mapping_path):
    print("Loading embeddings...")
    embeddings = np.load(embeddings_path, mmap_mode="r")
    if embeddings.ndim != 2:
        raise ValueError(f"Expected 2D embedding matrix, got shape {embeddings.shape}")
    print(f"Embeddings shape: {embeddings.shape}")
    print(f"Embeddings dtype: {embeddings.dtype}")

    print("Loading embedding mapping...")
    mapping = pd.read_parquet(mapping_path).copy()
    required_mapping = {"embedding_row", "pmid"}
    missing = required_mapping - set(mapping.columns)
    if missing:
        raise ValueError(f"Mapping is missing columns: {sorted(missing)}")

    mapping["pmid"] = mapping["pmid"].astype(str).str.strip()
    mapping = mapping.sort_values("embedding_row").reset_index(drop=True)
    if not np.array_equal(np.arange(len(mapping), dtype=np.int64), mapping["embedding_row"].to_numpy(dtype=np.int64)):
        raise ValueError("embedding_row is not exactly 0..N-1 after sorting.")
    if mapping["pmid"].duplicated().any():
        raise ValueError("Duplicate PMIDs found in embedding mapping.")
    if len(mapping) != embeddings.shape[0]:
        raise ValueError(f"Mapping rows ({len(mapping):,}) do not match embedding rows ({embeddings.shape[0]:,}).")

    print("Loading analytical corpus...")
    corpus = pd.read_parquet(corpus_path).copy()
    required_corpus = {"pmid", "document"}
    missing = required_corpus - set(corpus.columns)
    if missing:
        raise ValueError(f"Corpus is missing columns: {sorted(missing)}")
    corpus["pmid"] = corpus["pmid"].astype(str).str.strip()
    if corpus["pmid"].duplicated().any():
        raise ValueError("Duplicate PMIDs found in corpus.")

    aligned = mapping[["embedding_row", "pmid"]].merge(corpus, on="pmid", how="left", validate="one_to_one", indicator=True)
    missing_docs = aligned["_merge"].ne("both")
    if missing_docs.any():
        ids = aligned.loc[missing_docs, "pmid"].tolist()[:20]
        raise ValueError(f"{missing_docs.sum():,} embedding PMIDs are missing from corpus. Examples: {ids}")
    aligned = aligned.drop(columns="_merge")
    aligned["document"] = aligned["document"].fillna("").astype(str).str.strip()
    if aligned["document"].eq("").any():
        raise ValueError(f"{aligned['document'].eq('').sum():,} aligned documents are blank.")

    norms = np.linalg.norm(embeddings, axis=1)
    if not np.isfinite(norms).all():
        raise ValueError("Embedding norms contain NaN/Inf.")
    print(f"Embedding norm QC: min={norms.min():.6f}, mean={norms.mean():.6f}, max={norms.max():.6f}")
    return aligned, embeddings


def flatten_topic_terms(topic_model, topic_ids, top_n_words):
    rows = []
    for topic_id in topic_ids:
        terms = topic_model.get_topic(int(topic_id))
        if not terms:
            continue
        for rank, (term, score) in enumerate(terms[:top_n_words], start=1):
            rows.append({"topic": int(topic_id), "rank": rank, "term": term, "ctfidf_score": float(score)})
    return pd.DataFrame(rows)


def top_documents_by_topic(assignments, n_per_topic):
    rows = []
    non_outliers = assignments[assignments["topic"] != -1].copy()
    for topic_id, group in non_outliers.groupby("topic", sort=True):
        group = group.sort_values(["cluster_membership_strength", "publication_year"], ascending=[False, False], na_position="last").head(n_per_topic)
        for rank, (_, row) in enumerate(group.iterrows(), start=1):
            rows.append({
                "topic": int(topic_id),
                "rank": rank,
                "pmid": row["pmid"],
                "publication_year": row.get("publication_year"),
                "cluster_membership_strength": row["cluster_membership_strength"],
                "title": row.get("title", ""),
            })
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    out = args.output_dir
    if out.exists() and any(out.iterdir()) and not args.overwrite:
        raise RuntimeError(f"Output directory '{out}' is not empty. Use a new directory or pass --overwrite.")
    out.mkdir(parents=True, exist_ok=True)
    complete_path = out / "COMPLETE"
    if complete_path.exists():
        complete_path.unlink()

    versions = save_environment(out)
    aligned, embeddings = load_aligned_data(args.corpus, args.embeddings, args.mapping)
    docs = aligned["document"].tolist()
    n_docs = len(docs)
    max_features = None if args.max_features == 0 else args.max_features

    config = {
        "n_documents": n_docs,
        "embedding_shape": list(embeddings.shape),
        "umap": {
            "n_neighbors": args.n_neighbors,
            "n_components": args.n_components,
            "min_dist": args.min_dist,
            "metric": args.umap_metric,
            "random_state": args.random_state,
        },
        "hdbscan": {
            "min_cluster_size": args.min_cluster_size,
            "min_samples": args.min_samples,
            "metric": args.hdbscan_metric,
            "cluster_selection_method": args.cluster_selection_method,
            "prediction_data": True,
        },
        "vectorizer": {
            "stop_words": "english",
            "ngram_range": [args.ngram_min, args.ngram_max],
            "max_features": max_features,
        },
        "ctfidf": {
            "reduce_frequent_words": not args.no_reduce_frequent_words,
            "bm25_weighting": False,
        },
        "bertopic": {
            "top_n_words": args.top_n_words,
            "calculate_probabilities": False,
            "nr_topics": None,
        },
        "versions": versions,
    }
    write_json(out / "baseline_config.json", config)

    print("\n=== BASELINE CONFIGURATION ===")
    print(json.dumps(config, indent=2))

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
        max_features=max_features,
        lowercase=True,
    )

    ctfidf_model = ClassTfidfTransformer(
        bm25_weighting=False,
        reduce_frequent_words=not args.no_reduce_frequent_words,
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

    print("\n=== FITTING BERTOPIC ===")
    start = time.time()
    topics, probabilities = topic_model.fit_transform(docs, embeddings)
    elapsed = time.time() - start
    topics = np.asarray(topics, dtype=np.int32)

    strengths = getattr(hdbscan_model, "probabilities_", None)
    if strengths is None:
        if probabilities is not None and np.ndim(probabilities) == 1:
            strengths = np.asarray(probabilities, dtype=np.float32)
        else:
            strengths = np.full(n_docs, np.nan, dtype=np.float32)
    else:
        strengths = np.asarray(strengths, dtype=np.float32)

    outlier_scores = getattr(hdbscan_model, "outlier_scores_", None)
    if outlier_scores is None:
        outlier_scores = np.full(n_docs, np.nan, dtype=np.float32)
    else:
        outlier_scores = np.asarray(outlier_scores, dtype=np.float32)

    hdbscan_labels = np.asarray(hdbscan_model.labels_, dtype=np.int32)
    if len(topics) != n_docs:
        raise RuntimeError("Topic assignment count does not match corpus.")

    reduced = np.asarray(umap_model.embedding_, dtype=np.float32)
    np.save(out / "umap_5d.npy", reduced)

    assignment_columns = [c for c in [
        "embedding_row", "pmid", "publication_year", "title", "doi",
        "publication_types_json", "mesh_descriptors_json", "mesh_terms_json", "keywords_json"
    ] if c in aligned.columns]
    assignments = aligned[assignment_columns].copy()
    assignments["topic"] = topics
    assignments["hdbscan_label"] = hdbscan_labels
    assignments["cluster_membership_strength"] = strengths
    assignments["hdbscan_outlier_score"] = outlier_scores
    assignments.to_parquet(out / "topic_assignments.parquet", index=False)
    assignments.to_csv(out / "topic_assignments.csv", index=False)

    topic_info = topic_model.get_topic_info().copy()
    topic_info.to_csv(out / "topic_info.csv", index=False)

    topic_ids = sorted(int(x) for x in np.unique(topics))
    flatten_topic_terms(topic_model, topic_ids, args.top_n_words).to_csv(out / "topic_terms.csv", index=False)
    top_documents_by_topic(assignments, args.top_documents_per_topic).to_csv(out / "topic_top_documents.csv", index=False)

    rep_rows = []
    for topic_id in topic_ids:
        if topic_id == -1:
            continue
        try:
            reps = topic_model.get_representative_docs(topic_id) or []
        except Exception:
            reps = []
        for rank, text in enumerate(reps, start=1):
            rep_rows.append({"topic": topic_id, "rank": rank, "representative_document": text})
    pd.DataFrame(rep_rows).to_csv(out / "bertopic_representative_documents.csv", index=False)

    n_outliers = int(np.sum(topics == -1))
    real_topic_ids = [x for x in topic_ids if x != -1]
    n_topics = len(real_topic_ids)
    topic_sizes = assignments.loc[assignments["topic"] != -1].groupby("topic").size().sort_values(ascending=False)

    summary = {
        "n_documents": n_docs,
        "n_topics_excluding_outlier": n_topics,
        "n_outliers": n_outliers,
        "outlier_percentage": 100.0 * n_outliers / n_docs if n_docs else None,
        "largest_topic_size": int(topic_sizes.max()) if len(topic_sizes) else None,
        "smallest_topic_size": int(topic_sizes.min()) if len(topic_sizes) else None,
        "median_topic_size": float(topic_sizes.median()) if len(topic_sizes) else None,
        "fit_elapsed_seconds": elapsed,
        "configuration_file": "baseline_config.json",
    }
    write_json(out / "baseline_summary.json", summary)

    light_model_dir = out / "bertopic_model_safetensors"
    try:
        topic_model.save(str(light_model_dir), serialization="safetensors", save_ctfidf=True)
        print(f"Saved light BERTopic model: {light_model_dir}")
    except Exception as exc:
        print(f"Warning: safetensors BERTopic save failed: {type(exc).__name__}: {exc}")
        (out / "MODEL_SAVE_WARNING.txt").write_text(f"{type(exc).__name__}: {exc}\n", encoding="utf-8")

    submodels = out / "submodels"
    submodels.mkdir(exist_ok=True)
    joblib.dump(umap_model, submodels / "umap_model.joblib")
    joblib.dump(hdbscan_model, submodels / "hdbscan_model.joblib")
    joblib.dump(vectorizer_model, submodels / "vectorizer_model.joblib")
    joblib.dump(ctfidf_model, submodels / "ctfidf_model.joblib")

    if args.save_full_pickle:
        print("Saving full BERTopic pickle. Only load this file if you trust its source.")
        topic_model.save(str(out / "bertopic_full_model.pkl"), serialization="pickle")

    complete_path.write_text("BERTopic baseline completed successfully.\n", encoding="utf-8")

    print("\n=== BERTOPIC BASELINE COMPLETE ===")
    print(f"Documents:          {n_docs:,}")
    print(f"Topics (no -1):     {n_topics:,}")
    print(f"Outliers (-1):      {n_outliers:,}")
    print(f"Outlier percentage: {summary['outlier_percentage']:.2f}%")
    if len(topic_sizes):
        print(f"Largest topic:      {int(topic_sizes.max()):,}")
        print(f"Median topic size:  {topic_sizes.median():.1f}")
        print(f"Smallest topic:     {int(topic_sizes.min()):,}")

    print("\nTop topics by size:")
    display_cols = [c for c in ["Topic", "Count", "Name"] if c in topic_info.columns]
    print(topic_info[topic_info["Topic"] != -1].head(20)[display_cols].to_string(index=False))
    print(f"\nOutputs saved to: {out}")


if __name__ == "__main__":
    main()
