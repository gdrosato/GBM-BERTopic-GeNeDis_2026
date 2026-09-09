import pandas as pd
import numpy as np
from sentence_transformers import SentenceTransformer
from transformers import AutoTokenizer


MODEL_NAME = "NeuML/biomedbert-base-embeddings"
INPUT_FILE = "glioblastoma_analysis_corpus.parquet"


# --------------------------------------------------
# Load corpus
# --------------------------------------------------

df = pd.read_parquet(INPUT_FILE)

print(f"Documents: {len(df):,}")


# --------------------------------------------------
# Load model/tokenizer
# --------------------------------------------------

model = SentenceTransformer(MODEL_NAME)

print(f"Model: {MODEL_NAME}")
print(f"SentenceTransformer max_seq_length: {model.max_seq_length}")


tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)


# --------------------------------------------------
# Count tokens WITHOUT truncation
# --------------------------------------------------

documents = df["document"].fillna("").tolist()

token_lengths = []

BATCH_SIZE = 512

for start in range(0, len(documents), BATCH_SIZE):

    batch = documents[start:start + BATCH_SIZE]

    encoded = tokenizer(
        batch,
        add_special_tokens=True,
        truncation=False,
        padding=False,
        return_length=True,
    )

    token_lengths.extend(encoded["length"])

    if start % 5000 < BATCH_SIZE:
        print(
            f"Processed "
            f"{min(start + BATCH_SIZE, len(documents)):,}"
            f"/{len(documents):,}"
        )


df["token_length"] = token_lengths


# --------------------------------------------------
# Statistics
# --------------------------------------------------

print("\n=== TOKEN LENGTH ===")

print(
    df["token_length"].describe(
        percentiles=[
            0.01,
            0.05,
            0.10,
            0.25,
            0.50,
            0.75,
            0.90,
            0.95,
            0.99,
        ]
    )
)


MAX_LENGTH = model.max_seq_length


n_over = (df["token_length"] > MAX_LENGTH).sum()

n_over_768 = (df["token_length"] > 768).sum()

n_over_1024 = (df["token_length"] > 1024).sum()


print(f"\nDocuments > {MAX_LENGTH} tokens: {n_over:,}")
print(
    f"Percentage > {MAX_LENGTH}: "
    f"{100 * n_over / len(df):.2f}%"
)

print(f"Documents > 768 tokens: {n_over_768:,}")
print(f"Documents > 1024 tokens: {n_over_1024:,}")


# --------------------------------------------------
# Save long-document inspection file
# --------------------------------------------------

long_docs = df[
    df["token_length"] > MAX_LENGTH
][
    [
        "pmid",
        "publication_year",
        "title",
        "abstract",
        "token_length",
        "publication_types_json",
    ]
].sort_values(
    "token_length",
    ascending=False
)


long_docs.to_csv(
    "documents_over_512_tokens.csv",
    index=False
)


# --------------------------------------------------
# Save updated corpus
# --------------------------------------------------

df.to_parquet(
    "glioblastoma_analysis_corpus_with_tokens.parquet",
    index=False
)


print("\nSaved:")
print("documents_over_512_tokens.csv")
print("glioblastoma_analysis_corpus_with_tokens.parquet")
