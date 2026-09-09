#!/usr/bin/env python3

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sentence_transformers import SentenceTransformer
from transformers import AutoTokenizer


DEFAULT_MODEL = "NeuML/biomedbert-base-embeddings"
DEFAULT_EXCLUDED_PMIDS = {
    "27454254": "Conference abstract collection / aggregate PubMed record"
}


def choose_device(requested: str) -> str:
    if requested != "auto":
        return requested
    if torch.cuda.is_available():
        return "cuda"
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def l2_normalize(x: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(x)
    if norm == 0:
        return x.astype(np.float32)
    return (x / norm).astype(np.float32)


def chunk_ids(ids, chunk_size, overlap):
    step = chunk_size - overlap
    chunks = []
    start = 0
    while start < len(ids):
        end = min(start + chunk_size, len(ids))
        chunks.append(ids[start:end])
        if end >= len(ids):
            break
        start += step
    return chunks


def encode_direct(model, texts, batch_size):
    return model.encode(
        texts,
        batch_size=batch_size,
        show_progress_bar=False,
        convert_to_numpy=True,
        normalize_embeddings=True,
    ).astype(np.float32)


def encode_long(model, tokenizer, texts, chunk_size, overlap, batch_size):
    all_chunk_texts = []
    all_weights = []
    ranges = []

    for text in texts:
        ids = tokenizer.encode(
            text,
            add_special_tokens=False,
            truncation=False,
        )
        chunks = chunk_ids(ids, chunk_size, overlap)

        start = len(all_chunk_texts)
        for chunk in chunks:
            all_chunk_texts.append(
                tokenizer.decode(
                    chunk,
                    skip_special_tokens=True,
                    clean_up_tokenization_spaces=False,
                )
            )
            all_weights.append(max(len(chunk), 1))
        end = len(all_chunk_texts)
        ranges.append((start, end))

    chunk_embeddings = model.encode(
        all_chunk_texts,
        batch_size=batch_size,
        show_progress_bar=False,
        convert_to_numpy=True,
        normalize_embeddings=False,
    ).astype(np.float32)

    weights = np.asarray(all_weights, dtype=np.float32)

    result = []
    for start, end in ranges:
        pooled = np.average(
            chunk_embeddings[start:end],
            axis=0,
            weights=weights[start:end],
        )
        result.append(l2_normalize(pooled))

    return np.vstack(result).astype(np.float32)


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input",
        type=Path,
        default=Path("glioblastoma_analysis_corpus_with_tokens.parquet"),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("embeddings"),
    )
    p.add_argument(
        "--model",
        default=DEFAULT_MODEL,
    )
    p.add_argument(
        "--chunk-size",
        type=int,
        default=480,
    )
    p.add_argument(
        "--chunk-overlap",
        type=int,
        default=50,
    )
    p.add_argument(
        "--doc-batch-size",
        type=int,
        default=256,
    )
    p.add_argument(
        "--encode-batch-size",
        type=int,
        default=64,
    )
    p.add_argument(
        "--device",
        default="auto",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
    )
    args = p.parse_args()

    if args.chunk_overlap >= args.chunk_size:
        raise ValueError("chunk-overlap must be smaller than chunk-size")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    embeddings_path = args.output_dir / "biomedbert_embeddings.npy"
    mapping_path = args.output_dir / "embedding_pmids.parquet"
    exclusions_path = args.output_dir / "manual_exclusions.csv"
    checkpoint_path = args.output_dir / "checkpoint.json"
    metadata_path = args.output_dir / "embedding_metadata.json"
    complete_path = args.output_dir / "COMPLETE"

    if complete_path.exists() and not args.overwrite:
        print("Completed output already exists.")
        print("Use --overwrite only if you intentionally want to rebuild it.")
        return

    device = choose_device(args.device)

    print(f"Loading model: {args.model}")
    print(f"Device: {device}")

    model = SentenceTransformer(args.model, device=device)
    tokenizer = AutoTokenizer.from_pretrained(args.model)

    max_len = int(model.max_seq_length)
    dim = int(model.get_sentence_embedding_dimension())
    special_tokens = tokenizer.num_special_tokens_to_add(pair=False)

    if args.chunk_size + special_tokens > max_len:
        raise ValueError(
            f"chunk-size {args.chunk_size} + {special_tokens} special tokens "
            f"exceeds model limit {max_len}"
        )

    print(f"Embedding dimension: {dim}")
    print(f"Model max sequence length: {max_len}")
    print(f"Chunk size: {args.chunk_size}")
    print(f"Chunk overlap: {args.chunk_overlap}")

    df = pd.read_parquet(args.input).copy()

    required = {"pmid", "document"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df["pmid"] = df["pmid"].astype(str).str.strip()
    df["document"] = df["document"].fillna("").astype(str).str.strip()

    excluded = df[df["pmid"].isin(DEFAULT_EXCLUDED_PMIDS)].copy()

    pd.DataFrame(
        [
            {"pmid": pmid, "reason": reason}
            for pmid, reason in DEFAULT_EXCLUDED_PMIDS.items()
            if pmid in set(df["pmid"])
        ]
    ).to_csv(exclusions_path, index=False)

    df = df[~df["pmid"].isin(DEFAULT_EXCLUDED_PMIDS)].copy()
    df = df.reset_index(drop=True)
    df["embedding_row"] = np.arange(len(df), dtype=np.int64)

    if df["pmid"].duplicated().any():
        raise RuntimeError("Duplicate PMIDs remain in corpus.")

    if (df["document"] == "").any():
        raise RuntimeError("Blank documents remain in corpus.")

    print(f"Input records: {len(df) + len(excluded):,}")
    print(f"Manual exclusions: {len(excluded):,}")
    print(f"Final embedding corpus: {len(df):,}")

    if "token_length" in df.columns:
        token_lengths = pd.to_numeric(
            df["token_length"],
            errors="coerce",
        )
    else:
        token_lengths = pd.Series(np.nan, index=df.index)

    missing_lengths = token_lengths.isna()

    if missing_lengths.any():
        print(
            f"Calculating missing token lengths for "
            f"{missing_lengths.sum():,} documents..."
        )

        texts = df.loc[missing_lengths, "document"].tolist()
        enc = tokenizer(
            texts,
            add_special_tokens=True,
            truncation=False,
            padding=False,
            return_length=True,
        )
        token_lengths.loc[missing_lengths] = enc["length"]

    token_lengths = token_lengths.astype(int)
    df["embedding_token_length"] = token_lengths

    n_total = len(df)
    n_long = int((token_lengths > max_len).sum())
    n_direct = n_total - n_long

    print(f"Direct documents: {n_direct:,}")
    print(f"Long documents:   {n_long:,}")
    print(f"Long percentage:  {100*n_long/n_total:.2f}%")

    mapping_cols = [
        c for c in [
            "embedding_row",
            "pmid",
            "publication_year",
            "title",
            "doi",
            "embedding_token_length",
        ]
        if c in df.columns
    ]
    df[mapping_cols].to_parquet(mapping_path, index=False)

    start_row = 0

    if checkpoint_path.exists() and embeddings_path.exists() and not args.overwrite:
        state = json.loads(checkpoint_path.read_text(encoding="utf-8"))

        expected = {
            "model": args.model,
            "n_documents": n_total,
            "embedding_dim": dim,
            "chunk_size": args.chunk_size,
            "chunk_overlap": args.chunk_overlap,
        }

        for key, value in expected.items():
            if state.get(key) != value:
                raise RuntimeError(
                    f"Checkpoint mismatch for {key}: "
                    f"{state.get(key)} != {value}"
                )

        start_row = int(state.get("next_row", 0))
        print(f"Resuming from row: {start_row:,}")

        emb = np.load(embeddings_path, mmap_mode="r+")
        if emb.shape != (n_total, dim):
            raise RuntimeError(
                f"Existing embedding shape {emb.shape} "
                f"does not match {(n_total, dim)}"
            )

    else:
        if embeddings_path.exists() and not args.overwrite:
            raise RuntimeError(
                f"{embeddings_path} already exists without a usable checkpoint. "
                "Use --overwrite to rebuild."
            )

        emb = np.lib.format.open_memmap(
            embeddings_path,
            mode="w+",
            dtype=np.float32,
            shape=(n_total, dim),
        )

        state = {
            "model": args.model,
            "n_documents": n_total,
            "embedding_dim": dim,
            "chunk_size": args.chunk_size,
            "chunk_overlap": args.chunk_overlap,
            "next_row": 0,
        }
        checkpoint_path.write_text(
            json.dumps(state, indent=2),
            encoding="utf-8",
        )

    t0 = time.time()

    for batch_start in range(
        start_row,
        n_total,
        args.doc_batch_size,
    ):
        batch_end = min(
            batch_start + args.doc_batch_size,
            n_total,
        )

        batch = df.iloc[batch_start:batch_end]
        lengths = token_lengths.iloc[batch_start:batch_end]

        direct_mask = (lengths <= max_len).to_numpy()
        long_mask = ~direct_mask

        direct_pos = np.flatnonzero(direct_mask)
        long_pos = np.flatnonzero(long_mask)

        if len(direct_pos):
            texts = batch.iloc[direct_pos]["document"].tolist()
            vectors = encode_direct(
                model,
                texts,
                args.encode_batch_size,
            )
            emb[batch_start + direct_pos] = vectors

        if len(long_pos):
            texts = batch.iloc[long_pos]["document"].tolist()
            vectors = encode_long(
                model,
                tokenizer,
                texts,
                args.chunk_size,
                args.chunk_overlap,
                args.encode_batch_size,
            )
            emb[batch_start + long_pos] = vectors

        emb.flush()

        state["next_row"] = batch_end
        checkpoint_path.write_text(
            json.dumps(state, indent=2),
            encoding="utf-8",
        )

        elapsed = time.time() - t0
        processed = batch_end - start_row
        rate = processed / elapsed if elapsed else 0.0

        print(
            f"Embedded {batch_end:,}/{n_total:,} "
            f"({100*batch_end/n_total:6.2f}%) | "
            f"direct={len(direct_pos):,}, "
            f"long={len(long_pos):,} | "
            f"{rate:,.1f} docs/s"
        )

    emb.flush()

    final = np.load(embeddings_path, mmap_mode="r")

    if final.shape != (n_total, dim):
        raise RuntimeError("Unexpected final embedding shape.")

    if not np.isfinite(final).all():
        raise RuntimeError("Embeddings contain NaN or infinity.")

    norms = np.linalg.norm(final, axis=1)

    metadata = {
        "model": args.model,
        "input_file": str(args.input),
        "n_documents": n_total,
        "embedding_dimension": dim,
        "dtype": "float32",
        "device": device,
        "model_max_seq_length": max_len,
        "chunk_size_content_tokens": args.chunk_size,
        "chunk_overlap_tokens": args.chunk_overlap,
        "direct_documents": n_direct,
        "long_documents": n_long,
        "manual_exclusions": DEFAULT_EXCLUDED_PMIDS,
        "embedding_norm_min": float(norms.min()),
        "embedding_norm_mean": float(norms.mean()),
        "embedding_norm_max": float(norms.max()),
    }

    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    complete_path.write_text(
        "Embedding generation completed successfully.\n",
        encoding="utf-8",
    )

    print("\n=== COMPLETE ===")
    print(f"Embeddings shape: {final.shape}")
    print(f"Embedding dtype:  {final.dtype}")
    print(
        f"L2 norms: min={norms.min():.6f}, "
        f"mean={norms.mean():.6f}, "
        f"max={norms.max():.6f}"
    )
    print(f"Embeddings:   {embeddings_path}")
    print(f"PMID mapping: {mapping_path}")
    print(f"Metadata:     {metadata_path}")
    print(f"Exclusions:   {exclusions_path}")


if __name__ == "__main__":
    main()
