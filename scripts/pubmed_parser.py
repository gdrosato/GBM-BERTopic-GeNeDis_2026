#!/usr/bin/env python3
"""
Streaming PubMed XML parser for BERTopic / bibliometric analysis.

Input:
    PubMed XML exported with NCBI EDirect / EFetch.
    Supports .xml and .xml.gz.

Output:
    Parquet (recommended) and/or CSV with one row per PMID.

Example:
    python pubmed_parser.py \
        --input glioblastoma_pubmed_2000_2026.xml.gz \
        --parquet glioblastoma_pubmed.parquet \
        --csv glioblastoma_pubmed.csv

Dependencies:
    pip install lxml pandas pyarrow
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path
from typing import Dict, Iterator, List, Optional, TextIO, BinaryIO

import pandas as pd
from lxml import etree


YEAR_RE = re.compile(r"\b(18|19|20|21)\d{2}\b")
WS_RE = re.compile(r"\s+")


def clean_text(value: Optional[str]) -> str:
    """Normalize whitespace and return an empty string for None."""
    if not value:
        return ""
    return WS_RE.sub(" ", value).strip()


def element_text(elem: Optional[etree._Element]) -> str:
    """Return all text inside an XML element, including nested inline tags."""
    if elem is None:
        return ""
    return clean_text("".join(elem.itertext()))


def first_text(node: etree._Element, xpath: str) -> str:
    elem = node.find(xpath)
    return element_text(elem)


def extract_year(text: str) -> Optional[int]:
    """Extract the first plausible 4-digit year from strings such as '2024 Sep-Oct'."""
    if not text:
        return None
    match = YEAR_RE.search(text)
    return int(match.group(0)) if match else None


def parse_pubdate(article: etree._Element) -> Dict[str, object]:
    """
    Canonical publication date.

    Priority:
      1. JournalIssue/PubDate (the main PubMed bibliographic publication date)
      2. ArticleDate
      3. PubMed history date with PubStatus='pubmed'
      4. DateCompleted

    Keeps the raw bibliographic PubDate as well.
    """
    pubdate = article.find("./MedlineCitation/Article/Journal/JournalIssue/PubDate")

    raw_pub_date = ""
    year = None
    month = ""
    day = ""
    date_source = ""

    if pubdate is not None:
        y = first_text(pubdate, "./Year")
        m = first_text(pubdate, "./Month")
        d = first_text(pubdate, "./Day")
        medline_date = first_text(pubdate, "./MedlineDate")

        if y:
            year = extract_year(y)
            raw_pub_date = " ".join(x for x in [y, m, d] if x)
            month = m
            day = d
            date_source = "JournalIssue/PubDate"
        elif medline_date:
            year = extract_year(medline_date)
            raw_pub_date = medline_date
            date_source = "JournalIssue/PubDate/MedlineDate"

    # Fallback 1: ArticleDate
    if year is None:
        article_date = article.find("./MedlineCitation/Article/ArticleDate")
        if article_date is not None:
            y = first_text(article_date, "./Year")
            m = first_text(article_date, "./Month")
            d = first_text(article_date, "./Day")
            if y:
                year = extract_year(y)
                raw_pub_date = " ".join(x for x in [y, m, d] if x)
                month = m
                day = d
                date_source = "ArticleDate"

    # Fallback 2: PubMed history date
    if year is None:
        for hist_date in article.findall(
            "./PubmedData/History/PubMedPubDate"
        ):
            if hist_date.get("PubStatus") == "pubmed":
                y = first_text(hist_date, "./Year")
                m = first_text(hist_date, "./Month")
                d = first_text(hist_date, "./Day")
                if y:
                    year = extract_year(y)
                    raw_pub_date = " ".join(x for x in [y, m, d] if x)
                    month = m
                    day = d
                    date_source = "PubMedPubDate:pubmed"
                    break

    # Fallback 3: DateCompleted
    if year is None:
        completed = article.find("./MedlineCitation/DateCompleted")
        if completed is not None:
            y = first_text(completed, "./Year")
            m = first_text(completed, "./Month")
            d = first_text(completed, "./Day")
            if y:
                year = extract_year(y)
                raw_pub_date = " ".join(x for x in [y, m, d] if x)
                month = m
                day = d
                date_source = "DateCompleted"

    return {
        "publication_year": year,
        "publication_month": month,
        "publication_day": day,
        "publication_date_raw": raw_pub_date,
        "publication_date_source": date_source,
    }


def parse_abstract(article: etree._Element) -> str:
    """
    Concatenate all AbstractText sections.
    Preserves labels such as BACKGROUND, METHODS, RESULTS when present.
    """
    parts: List[str] = []

    for abstract_text in article.findall(
        "./MedlineCitation/Article/Abstract/AbstractText"
    ):
        text = element_text(abstract_text)
        if not text:
            continue

        label = clean_text(abstract_text.get("Label", ""))
        nlm_category = clean_text(abstract_text.get("NlmCategory", ""))

        prefix = label or nlm_category
        if prefix and not text.upper().startswith(prefix.upper() + ":"):
            parts.append(f"{prefix}: {text}")
        else:
            parts.append(text)

    # Some records may use OtherAbstract
    if not parts:
        for abstract_text in article.findall(
            "./MedlineCitation/OtherAbstract/AbstractText"
        ):
            text = element_text(abstract_text)
            if text:
                parts.append(text)

    return clean_text(" ".join(parts))


def parse_authors(article: etree._Element) -> List[str]:
    authors: List[str] = []

    for author in article.findall(
        "./MedlineCitation/Article/AuthorList/Author"
    ):
        collective = first_text(author, "./CollectiveName")
        if collective:
            authors.append(collective)
            continue

        last = first_text(author, "./LastName")
        fore = first_text(author, "./ForeName")
        initials = first_text(author, "./Initials")

        if last and fore:
            authors.append(f"{last} {fore}")
        elif last and initials:
            authors.append(f"{last} {initials}")
        elif last:
            authors.append(last)

    return authors


def parse_affiliations(article: etree._Element) -> List[str]:
    seen = set()
    affiliations: List[str] = []

    for aff in article.findall(
        "./MedlineCitation/Article/AuthorList/Author/AffiliationInfo/Affiliation"
    ):
        text = element_text(aff)
        if text and text not in seen:
            seen.add(text)
            affiliations.append(text)

    return affiliations


def parse_mesh(article: etree._Element) -> Dict[str, List[str]]:
    descriptors: List[str] = []
    qualified_terms: List[str] = []

    for heading in article.findall("./MedlineCitation/MeshHeadingList/MeshHeading"):
        descriptor = heading.find("./DescriptorName")
        descriptor_text = element_text(descriptor)

        if descriptor_text:
            descriptors.append(descriptor_text)

        qualifiers = [
            element_text(q)
            for q in heading.findall("./QualifierName")
            if element_text(q)
        ]

        if descriptor_text and qualifiers:
            for qualifier in qualifiers:
                qualified_terms.append(f"{descriptor_text}/{qualifier}")
        elif descriptor_text:
            qualified_terms.append(descriptor_text)

    # preserve order, remove duplicates
    descriptors = list(dict.fromkeys(descriptors))
    qualified_terms = list(dict.fromkeys(qualified_terms))

    return {
        "mesh_descriptors": descriptors,
        "mesh_terms": qualified_terms,
    }


def parse_keywords(article: etree._Element) -> List[str]:
    keywords = []
    for kw in article.findall("./MedlineCitation/KeywordList/Keyword"):
        text = element_text(kw)
        if text:
            keywords.append(text)
    return list(dict.fromkeys(keywords))


def parse_publication_types(article: etree._Element) -> List[str]:
    values = []
    for elem in article.findall(
        "./MedlineCitation/Article/PublicationTypeList/PublicationType"
    ):
        text = element_text(elem)
        if text:
            values.append(text)
    return list(dict.fromkeys(values))


def parse_article_ids(article: etree._Element) -> Dict[str, str]:
    ids: Dict[str, str] = {}

    for elem in article.findall("./PubmedData/ArticleIdList/ArticleId"):
        id_type = clean_text(elem.get("IdType", "")).lower()
        value = element_text(elem)
        if id_type and value and id_type not in ids:
            ids[id_type] = value

    # DOI fallback from ELocationID
    if not ids.get("doi"):
        for elem in article.findall("./MedlineCitation/Article/ELocationID"):
            eid_type = clean_text(elem.get("EIdType", "")).lower()
            if eid_type == "doi":
                value = element_text(elem)
                if value:
                    ids["doi"] = value
                    break

    return ids


def parse_journal(article: etree._Element) -> Dict[str, str]:
    journal = article.find("./MedlineCitation/Article/Journal")
    if journal is None:
        return {
            "journal_title": "",
            "journal_iso": "",
            "issn": "",
            "volume": "",
            "issue": "",
        }

    return {
        "journal_title": first_text(journal, "./Title"),
        "journal_iso": first_text(journal, "./ISOAbbreviation"),
        "issn": first_text(journal, "./ISSN"),
        "volume": first_text(journal, "./JournalIssue/Volume"),
        "issue": first_text(journal, "./JournalIssue/Issue"),
    }


def parse_record(article: etree._Element) -> Dict[str, object]:
    citation = article.find("./MedlineCitation")
    article_node = article.find("./MedlineCitation/Article")

    pmid = first_text(article, "./MedlineCitation/PMID")
    title = first_text(article, "./MedlineCitation/Article/ArticleTitle")
    abstract = parse_abstract(article)

    authors = parse_authors(article)
    affiliations = parse_affiliations(article)
    mesh = parse_mesh(article)
    keywords = parse_keywords(article)
    publication_types = parse_publication_types(article)
    article_ids = parse_article_ids(article)
    date_fields = parse_pubdate(article)
    journal_fields = parse_journal(article)

    language_values = [
        element_text(x)
        for x in article.findall("./MedlineCitation/Article/Language")
        if element_text(x)
    ]

    # The BERTopic input document.
    if title and abstract:
        document = f"{title}. {abstract}"
    else:
        document = title or abstract

    record = {
        "pmid": pmid,
        **date_fields,
        "title": title,
        "abstract": abstract,
        "document": clean_text(document),
        **journal_fields,
        "doi": article_ids.get("doi", ""),
        "pmc": article_ids.get("pmc", ""),
        "pii": article_ids.get("pii", ""),
        "authors_json": json.dumps(authors, ensure_ascii=False),
        "first_author": authors[0] if authors else "",
        "last_author": authors[-1] if authors else "",
        "affiliations_json": json.dumps(affiliations, ensure_ascii=False),
        "publication_types_json": json.dumps(publication_types, ensure_ascii=False),
        "mesh_descriptors_json": json.dumps(
            mesh["mesh_descriptors"], ensure_ascii=False
        ),
        "mesh_terms_json": json.dumps(mesh["mesh_terms"], ensure_ascii=False),
        "keywords_json": json.dumps(keywords, ensure_ascii=False),
        "languages_json": json.dumps(
            list(dict.fromkeys(language_values)), ensure_ascii=False
        ),
    }

    return record


def open_xml(path: Path) -> BinaryIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def iter_pubmed_records(path: Path) -> Iterator[Dict[str, object]]:
    """
    Stream PubmedArticle records with low memory usage.
    """
    with open_xml(path) as handle:
        context = etree.iterparse(
            handle,
            events=("end",),
            tag="PubmedArticle",
            recover=True,
            huge_tree=True,
        )

        for _, elem in context:
            yield parse_record(elem)

            # Aggressively release parsed XML from memory.
            elem.clear()
            parent = elem.getparent()
            if parent is not None:
                while elem.getprevious() is not None:
                    del parent[0]

        del context


def parse_all(path: Path, progress_every: int = 5000) -> pd.DataFrame:
    rows = []
    for i, record in enumerate(iter_pubmed_records(path), start=1):
        rows.append(record)
        if progress_every and i % progress_every == 0:
            print(f"Parsed {i:,} PubMed records...")

    df = pd.DataFrame(rows)
    print(f"Finished: {len(df):,} records.")
    return df


def quality_control(df: pd.DataFrame) -> pd.DataFrame:
    """
    Basic QC appropriate for the current glioblastoma corpus:
      - remove blank PMIDs
      - remove duplicate PMIDs
      - require non-empty title and abstract
      - report publication-year issues
    """
    before = len(df)

    df = df[df["pmid"].astype(str).str.strip().ne("")].copy()
    df = df.drop_duplicates(subset=["pmid"], keep="first").copy()

    df["title"] = df["title"].fillna("").astype(str).str.strip()
    df["abstract"] = df["abstract"].fillna("").astype(str).str.strip()
    df["document"] = df["document"].fillna("").astype(str).str.strip()

    df = df[(df["title"] != "") & (df["abstract"] != "")].copy()

    # Nullable integer keeps missing years as <NA>.
    df["publication_year"] = pd.to_numeric(
        df["publication_year"], errors="coerce"
    ).astype("Int64")

    print("\nQC summary")
    print("----------")
    print(f"Input rows:              {before:,}")
    print(f"Rows after QC:           {len(df):,}")
    print(f"Unique PMIDs:            {df['pmid'].nunique():,}")
    print(f"Missing publication year:{df['publication_year'].isna().sum():,}")
    print(f"Missing DOI:             {(df['doi'].fillna('') == '').sum():,}")

    if df["publication_year"].notna().any():
        print(
            "Publication-year range:  "
            f"{df['publication_year'].min()}–{df['publication_year'].max()}"
        )

    return df


def save_outputs(
    df: pd.DataFrame,
    parquet_path: Optional[Path],
    csv_path: Optional[Path],
) -> None:
    if parquet_path:
        parquet_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(parquet_path, index=False)
        print(f"Saved Parquet: {parquet_path}")

    if csv_path:
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(csv_path, index=False, encoding="utf-8")
        print(f"Saved CSV:     {csv_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse PubMed XML/XML.GZ into a BERTopic-ready dataset."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--parquet", type=Path)
    parser.add_argument("--csv", type=Path)
    parser.add_argument(
        "--progress-every",
        type=int,
        default=5000,
        help="Print progress every N records (default: 5000).",
    )
    args = parser.parse_args()

    if not args.parquet and not args.csv:
        parser.error("Specify at least one of --parquet or --csv.")

    df = parse_all(args.input, progress_every=args.progress_every)
    df = quality_control(df)

    # Sort only after parsing so the raw XML order does not matter.
    if "publication_year" in df.columns:
        df = df.sort_values(
            ["publication_year", "pmid"],
            na_position="last",
        ).reset_index(drop=True)

    save_outputs(df, args.parquet, args.csv)


if __name__ == "__main__":
    main()
