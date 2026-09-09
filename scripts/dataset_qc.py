import pandas as pd
import json

df = pd.read_parquet("glioblastoma_pubmed.parquet")

print("\n=== BASIC DATASET INFO ===")
print(f"Rows: {len(df):,}")
print(f"Unique PMIDs: {df['pmid'].nunique():,}")
print(
    f"Years: {df['publication_year'].min()} "
    f"- {df['publication_year'].max()}"
)

# --------------------------------------------------
# 1. Papers per year
# --------------------------------------------------

year_counts = (
    df.groupby("publication_year")
      .size()
      .sort_index()
)

print("\n=== PAPERS PER YEAR ===")
print(year_counts.to_string())

year_counts.to_csv(
    "papers_per_year.csv",
    header=["n_papers"]
)

# --------------------------------------------------
# 2. Abstract length
# --------------------------------------------------

df["abstract_chars"] = df["abstract"].str.len()
df["abstract_words"] = df["abstract"].str.split().str.len()

print("\n=== ABSTRACT LENGTH ===")
print(df["abstract_words"].describe(
    percentiles=[0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
))

print("\nAbstracts < 50 words:")
print((df["abstract_words"] < 50).sum())

print("\nAbstracts < 100 words:")
print((df["abstract_words"] < 100).sum())

# --------------------------------------------------
# 3. Publication types
# --------------------------------------------------

publication_type_counts = {}

for value in df["publication_types_json"]:
    types = json.loads(value)

    for publication_type in types:
        publication_type_counts[publication_type] = (
            publication_type_counts.get(publication_type, 0) + 1
        )

publication_type_df = (
    pd.DataFrame(
        publication_type_counts.items(),
        columns=["publication_type", "count"]
    )
    .sort_values("count", ascending=False)
)

print("\n=== PUBLICATION TYPES ===")
print(publication_type_df.head(30).to_string(index=False))

publication_type_df.to_csv(
    "publication_types.csv",
    index=False
)

# --------------------------------------------------
# 4. Reviews
# --------------------------------------------------

def has_pub_type(value, target):
    return target in json.loads(value)

df["is_review"] = df["publication_types_json"].apply(
    lambda x: has_pub_type(x, "Review")
)

print("\n=== REVIEWS ===")
print(f"Reviews: {df['is_review'].sum():,}")
print(
    f"Percentage: "
    f"{100 * df['is_review'].mean():.2f}%"
)

# --------------------------------------------------
# 5. Duplicate titles
# --------------------------------------------------

normalized_titles = (
    df["title"]
    .str.lower()
    .str.strip()
)

duplicate_titles = normalized_titles.duplicated(
    keep=False
)

print("\n=== DUPLICATE TITLES ===")
print(
    f"Records with duplicated normalized title: "
    f"{duplicate_titles.sum():,}"
)

if duplicate_titles.any():
    duplicates = df.loc[
        duplicate_titles,
        ["pmid", "publication_year", "title", "doi"]
    ].sort_values("title")

    duplicates.to_csv(
        "possible_duplicate_titles.csv",
        index=False
    )

# --------------------------------------------------
# 6. Very short documents
# --------------------------------------------------

short_docs = df[df["abstract_words"] < 50][
    [
        "pmid",
        "publication_year",
        "title",
        "abstract",
        "publication_types_json"
    ]
]

short_docs.to_csv(
    "short_abstracts_under_50_words.csv",
    index=False
)

print("\nQC outputs saved.")
