import json
import pandas as pd


INPUT = "glioblastoma_pubmed.parquet"

OUTPUT = "glioblastoma_analysis_corpus.parquet"
OUTPUT_CSV = "glioblastoma_analysis_corpus.csv"


df = pd.read_parquet(INPUT)

print("Initial records:", f"{len(df):,}")


# ==========================================================
# 1. Parse publication types
# ==========================================================

def parse_json_list(value):
    if pd.isna(value) or value == "":
        return []

    try:
        return json.loads(value)
    except Exception:
        return []


df["publication_types"] = (
    df["publication_types_json"]
    .apply(parse_json_list)
)


# ==========================================================
# 2. Exclusion publication types
# ==========================================================

EXCLUDED_PUBLICATION_TYPES = {
    "Retracted Publication",
    "Retraction Notice",
    "Published Erratum",
    "Comment",
    "Editorial",
    "Letter",
    "Video-Audio Media",
    "Preprint",
}


def exclusion_reason(pub_types):

    matched = sorted(
        set(pub_types)
        & EXCLUDED_PUBLICATION_TYPES
    )

    if matched:
        return "; ".join(matched)

    return ""


df["exclusion_reason"] = (
    df["publication_types"]
    .apply(exclusion_reason)
)


excluded = df[
    df["exclusion_reason"] != ""
].copy()

included = df[
    df["exclusion_reason"] == ""
].copy()


print(
    "Excluded by publication type:",
    f"{len(excluded):,}"
)

print(
    "Remaining:",
    f"{len(included):,}"
)


excluded[
    [
        "pmid",
        "publication_year",
        "title",
        "doi",
        "publication_types_json",
        "exclusion_reason",
    ]
].to_csv(
    "excluded_publication_types.csv",
    index=False,
)


# ==========================================================
# 3. Abstract/document lengths
# ==========================================================

included["abstract_words"] = (
    included["abstract"]
    .fillna("")
    .str.split()
    .str.len()
)

included["document_words"] = (
    included["document"]
    .fillna("")
    .str.split()
    .str.len()
)


print("\n=== ABSTRACT LENGTH AFTER EXCLUSIONS ===")

print(
    included["abstract_words"].describe(
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


print(
    "\nAbstracts <50 words:",
    f"{(included['abstract_words'] < 50).sum():,}"
)

print(
    "Abstracts <100 words:",
    f"{(included['abstract_words'] < 100).sum():,}"
)


# ==========================================================
# 4. Inspect very short abstracts
# ==========================================================

short_abstracts = included[
    included["abstract_words"] < 50
].copy()

short_abstracts[
    [
        "pmid",
        "publication_year",
        "title",
        "abstract",
        "abstract_words",
        "doi",
        "publication_types_json",
    ]
].to_csv(
    "short_abstracts_after_filtering.csv",
    index=False,
)


# ==========================================================
# 5. Exact duplicated documents
# ==========================================================

included["normalized_document"] = (
    included["document"]
    .fillna("")
    .str.lower()
    .str.replace(r"\s+", " ", regex=True)
    .str.strip()
)


exact_document_duplicates = (
    included["normalized_document"]
    .duplicated(keep=False)
)


print(
    "\nExact duplicate document records:",
    f"{exact_document_duplicates.sum():,}"
)


included.loc[
    exact_document_duplicates,
    [
        "pmid",
        "publication_year",
        "title",
        "doi",
        "publication_types_json",
    ]
].sort_values(
    "title"
).to_csv(
    "exact_duplicate_documents.csv",
    index=False,
)


# ==========================================================
# 6. Duplicate DOI
# ==========================================================

doi_mask = (
    included["doi"]
    .fillna("")
    .str.strip()
    .ne("")
)


duplicated_doi = (
    included.loc[doi_mask, "doi"]
    .str.lower()
    .duplicated(keep=False)
)


duplicate_doi_rows = included.loc[
    included.index[doi_mask][duplicated_doi]
].copy()


print(
    "Records sharing DOI:",
    f"{len(duplicate_doi_rows):,}"
)


duplicate_doi_rows[
    [
        "pmid",
        "publication_year",
        "title",
        "doi",
        "publication_types_json",
    ]
].sort_values(
    "doi"
).to_csv(
    "duplicate_dois.csv",
    index=False,
)


# ==========================================================
# 7. Longest abstracts
# ==========================================================

included.nlargest(
    30,
    "abstract_words"
)[
    [
        "pmid",
        "publication_year",
        "title",
        "abstract_words",
        "publication_types_json",
    ]
].to_csv(
    "longest_abstracts.csv",
    index=False,
)


# ==========================================================
# 8. Papers per year after filtering
# ==========================================================

year_counts = (
    included.groupby("publication_year")
    .size()
    .sort_index()
)


print("\n=== PAPERS PER YEAR AFTER FILTERING ===")
print(year_counts.to_string())


year_counts.to_csv(
    "papers_per_year_after_filtering.csv",
    header=["n_papers"],
)


# ==========================================================
# 9. Save analytical corpus
# ==========================================================

# Drop helper fields we don't need to persist twice
included = included.drop(
    columns=[
        "normalized_document",
    ]
)


included.to_parquet(
    OUTPUT,
    index=False,
)


included.to_csv(
    OUTPUT_CSV,
    index=False,
)


print("\nFinal analytical corpus:")
print(f"{len(included):,} records")

print(
    "Unique PMIDs:",
    f"{included['pmid'].nunique():,}"
)

print(f"\nSaved: {OUTPUT}")
print(f"Saved: {OUTPUT_CSV}")
