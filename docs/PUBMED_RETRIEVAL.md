# Recreating the PubMed XML input with NCBI EDirect

The raw PubMed XML is intentionally **not distributed** in this repository.
Users can recreate the study input locally with NCBI Entrez Direct (EDirect)
and the exact query stored in [`config/pubmed_query.txt`](../config/pubmed_query.txt).

> **Important:** PubMed is a live database. Retrospective indexing corrections,
> retractions, and record updates can cause a later rerun of the same bounded
> query to return a slightly different count. Record the retrieval date and the
> observed counts for every reproduction run.

## 1. Install NCBI EDirect

On Linux or macOS, use the official NCBI installer:

```bash
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
```

For the current shell session:

```bash
export PATH=${HOME}/edirect:${PATH}
```

To make EDirect available in future shell sessions, add the corresponding PATH
line to your shell profile (for example `~/.bashrc`, `~/.bash_profile`, or
`~/.zshrc`). On Windows, use a Unix-compatible environment such as WSL.

Official documentation:

- https://www.ncbi.nlm.nih.gov/books/NBK179288/

## 2. Optional: configure an NCBI API key

NCBI supports up to **3 E-utility requests/second without an API key** and up to
**10 requests/second with an API key** by default. Create a key under **My NCBI
→ Account Settings → API Key Management**.

Set it only in your local shell environment:

```bash
export NCBI_API_KEY="YOUR_NCBI_API_KEY"
```

**Never commit an API key to GitHub, shell scripts, notebooks, logs, or
configuration files tracked by Git.**

Official API-key documentation:

- https://www.ncbi.nlm.nih.gov/books/NBK25497/
- https://www.ncbi.nlm.nih.gov/books/NBK53593/

## 3. Use exactly the study query

The query used in the study is:

```text
("glioblastoma"[MeSH Terms] OR "glioblastoma"[Title/Abstract])
AND ("2000/01/01"[Date - Publication] : "2026/08/27"[Date - Publication])
AND hasabstract
AND english[Language]
```

Do not maintain a second independently edited copy of the query in code. For
command-line reproduction, load it directly from the repository:

```bash
QUERY="$(tr '\n' ' ' < config/pubmed_query.txt)"
```

## 4. Check the current PubMed count before downloading

```bash
esearch -db pubmed -query "$QUERY" \
| xtract -pattern ENTREZ_DIRECT -element Count
```

During the original study workflow, the search interface initially showed
approximately **55,506** hits, while the XML dataset ultimately parsed for the
analysis contained **55,446 unique PubMed records**. Because PubMed can change
retrospectively, a later rerun should not be expected to reproduce the live
search count exactly.

## 5. Download the complete PubMed records as XML

Uncompressed XML:

```bash
esearch -db pubmed -query "$QUERY" \
| efetch -format xml \
> glioblastoma_pubmed_2000_2026.xml
```

Then compress it:

```bash
gzip glioblastoma_pubmed_2000_2026.xml
```

Alternatively, write the compressed file directly:

```bash
esearch -db pubmed -query "$QUERY" \
| efetch -format xml \
| gzip -c \
> glioblastoma_pubmed_2000_2026.xml.gz
```

The raw XML/XML.GZ file is a **local analysis input** and is ignored by Git in
this repository.

## 6. Verify downloaded PubMed records

Count PubMed articles:

```bash
zcat glioblastoma_pubmed_2000_2026.xml.gz \
| xtract -pattern PubmedArticle -element MedlineCitation/PMID \
| wc -l
```

Count unique PMIDs:

```bash
zcat glioblastoma_pubmed_2000_2026.xml.gz \
| xtract -pattern PubmedArticle -element MedlineCitation/PMID \
| sort -u \
| wc -l
```

For the dataset used in this study, both counts were:

```text
55446
```

If your rerun differs, record the new count and retrieval date rather than
silently forcing the historical number.

On macOS systems where `zcat` behaves differently, `gzcat` can be used instead.

## 7. Parse the XML

Once the XML has been recreated locally, continue with:

```bash
python scripts/pubmed_parser.py \
  --input glioblastoma_pubmed_2000_2026.xml.gz \
  --parquet glioblastoma_pubmed.parquet \
  --csv glioblastoma_pubmed.csv
```

The historical study dataset produced **55,446 parsed records with 55,446
unique PMIDs**.

## Automated helper

The repository also includes:

```bash
bash scripts/download_pubmed_edirect.sh
```

This helper reads the exact query from `config/pubmed_query.txt`, reports the
current PubMed hit count, downloads a compressed XML file, and reports total and
unique PMID counts. It warns if the reproduced unique count differs from the
historical study dataset, but does not fail solely because PubMed has changed.
