# Raw PubMed input

Raw PubMed XML/XML.GZ is intentionally not distributed in this repository and
is ignored by Git.

Recreate the input locally using NCBI EDirect:

- [`../../docs/PUBMED_RETRIEVAL.md`](../../docs/PUBMED_RETRIEVAL.md)

or run:

```bash
bash scripts/download_pubmed_edirect.sh
```

The historical XML dataset used in the study contained **55,446 unique PMIDs**.
A later PubMed rerun can differ because the database is updated retrospectively.
