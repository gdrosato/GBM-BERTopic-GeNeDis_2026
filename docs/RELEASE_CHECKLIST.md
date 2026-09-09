# Release checklist

## Repository contents

- [x] Complete Python analysis workflow included
- [x] PubMed EDirect retrieval instructions included
- [x] Exact environment export included
- [x] PMID-only manifests for retrieved and final modeling corpora included
- [x] Curated expert-annotation workbook included
- [x] Primary-label temporal outputs included
- [x] Axis-definition sensitivity outputs included
- [x] Denominator sensitivity outputs included
- [x] Main and supplementary figures included
- [x] Raw PubMed export excluded
- [x] Parsed abstract-bearing corpus excluded
- [x] Embeddings/model binaries excluded
- [x] Manuscript/LaTeX/BibTeX/PDF excluded
- [x] Lightweight GitHub Actions CI included
- [x] Machine-readable `CITATION.cff` included

## Before public release after acceptance

- [ ] Confirm final paper title in `README.md`, `CITATION.cff`, and `config/analysis_parameters.yaml`
- [ ] Add final journal/book-series citation, DOI, volume, and page range when available
- [ ] Confirm repository visibility is changed from private to public
- [ ] Decide whether to assign an explicit code/content license; if not, retain `docs/LICENSING.md`
- [ ] Run `make check` and confirm all tests pass
- [ ] Create the tagged public release
