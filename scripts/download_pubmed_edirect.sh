#!/usr/bin/env bash
set -euo pipefail

# Recreate the local PubMed XML input used by the analysis.
# Raw PubMed XML is intentionally ignored by Git.

QUERY_FILE="${QUERY_FILE:-config/pubmed_query.txt}"
OUTPUT="${OUTPUT:-glioblastoma_pubmed_2000_2026.xml.gz}"
HISTORICAL_UNIQUE_PMIDS=55446

for cmd in esearch efetch xtract gzip sort wc tr; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: required command '$cmd' was not found." >&2
    echo "Install NCBI EDirect first; see docs/PUBMED_RETRIEVAL.md." >&2
    exit 1
  fi
done

if [[ ! -f "$QUERY_FILE" ]]; then
  echo "ERROR: query file not found: $QUERY_FILE" >&2
  exit 1
fi

if [[ -z "${NCBI_API_KEY:-}" ]]; then
  echo "NOTE: NCBI_API_KEY is not set. EDirect will use the lower default request rate." >&2
fi

QUERY="$(tr '\n' ' ' < "$QUERY_FILE")"

echo "Query file: $QUERY_FILE"
echo "Output:     $OUTPUT"
echo

echo "Checking current PubMed count..."
CURRENT_COUNT="$({ esearch -db pubmed -query "$QUERY" \
  | xtract -pattern ENTREZ_DIRECT -element Count; } | tr -d '[:space:]')"
echo "Current PubMed search count: ${CURRENT_COUNT:-unknown}"

echo
echo "Downloading PubMed XML..."
esearch -db pubmed -query "$QUERY" \
  | efetch -format xml \
  | gzip -c \
  > "$OUTPUT"

echo "Download complete."

# Linux uses zcat; macOS may need gzcat.
if command -v gzcat >/dev/null 2>&1; then
  ZCAT="gzcat"
elif command -v zcat >/dev/null 2>&1; then
  ZCAT="zcat"
else
  echo "ERROR: neither zcat nor gzcat is available." >&2
  exit 1
fi

TOTAL_PMIDS="$($ZCAT "$OUTPUT" \
  | xtract -pattern PubmedArticle -element MedlineCitation/PMID \
  | wc -l | tr -d '[:space:]')"

UNIQUE_PMIDS="$($ZCAT "$OUTPUT" \
  | xtract -pattern PubmedArticle -element MedlineCitation/PMID \
  | sort -u \
  | wc -l | tr -d '[:space:]')"

echo "Total PubmedArticle PMIDs: $TOTAL_PMIDS"
echo "Unique PMIDs:              $UNIQUE_PMIDS"

if [[ "$UNIQUE_PMIDS" == "$HISTORICAL_UNIQUE_PMIDS" ]]; then
  echo "Historical unique PMID count matched: $HISTORICAL_UNIQUE_PMIDS"
else
  echo "WARNING: historical analysis used $HISTORICAL_UNIQUE_PMIDS unique PMIDs." >&2
  echo "A different count can occur because PubMed records/indexing can change retrospectively." >&2
  echo "Record the retrieval date and observed count for this rerun." >&2
fi
