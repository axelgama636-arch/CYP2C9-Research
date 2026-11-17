#!/usr/bin/env python3
"""Retrieve CYP2C9 sequencing/genotyping literature with resilient API fallback.

Primary source: Europe PMC. Fallback: NCBI PubMed E-utilities.
Search records are discovery candidates; method labels are text-derived and
must be checked against each full text before evidence synthesis.
"""
import csv
import json
import re
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

EPMC_API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
QUERY = (
    'CYP2C9 AND (sequencing OR "next generation sequencing" OR NGS '
    'OR Sanger OR pyrosequencing OR genotyping OR haplotype OR allele)'
)
PUBMED_QUERY = (
    'CYP2C9 AND (sequencing OR "next generation sequencing" OR NGS '
    'OR Sanger OR pyrosequencing OR genotyping OR haplotype OR allele)'
)
RETRYABLE = {429, 500, 502, 503, 504}
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "CYP2C9-Research/1.0 (literature discovery)"})

def clean(value):
    if value is None:
        return ""
    if isinstance(value, list):
        return "; ".join(str(x) for x in value)
    return re.sub(r"\\s+", " ", str(value)).strip()

def get_json(url, params, attempts=5):
    last_error = None
    for attempt in range(attempts):
        try:
            response = SESSION.get(url, params=params, timeout=(20, 90))
            if response.status_code in RETRYABLE:
                raise requests.HTTPError(
                    f"Temporary HTTP {response.status_code}", response=response
                )
            response.raise_for_status()
            return response.json()
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            status = getattr(getattr(exc, "response", None), "status_code", None)
            if status is not None and status not in RETRYABLE:
                raise
            if attempt == attempts - 1:
                break
            wait = min(60, 2 ** attempt * 3)
            print(f"API temporalmente no disponible ({exc}); reintento en {wait}s.")
            time.sleep(wait)
    raise RuntimeError(f"Fallaron {attempts} intentos de consulta: {last_error}")

def search_europe_pmc():
    params = {
        "query": QUERY, "format": "json", "pageSize": 1000,
        "resultType": "core", "sort": "RELEVANCE"
    }
    payload = get_json(EPMC_API, params)
    return payload.get("resultList", {}).get("result", []), payload

def search_pubmed():
    search = get_json(PUBMED_ESEARCH, {
        "db": "pubmed", "term": PUBMED_QUERY, "retmode": "json",
        "retmax": 1000, "sort": "relevance"
    })
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return [], {"source": "PubMed", "query": PUBMED_QUERY, "idlist": []}
    summary = get_json(PUBMED_ESUMMARY, {
        "db": "pubmed", "id": ",".join(ids), "retmode": "json"
    })
    result = summary.get("result", {})
    records = []
    for pmid in ids:
        item = result.get(pmid, {})
        authors = item.get("authors", [])
        records.append({
            "title": clean(item.get("title")),
            "authorString": "; ".join(clean(a.get("name")) for a in authors),
            "pubYear": clean(item.get("pubdate", "")[:4]),
            "journalTitle": clean(item.get("fulljournalname") or item.get("source")),
            "pmid": pmid, "pmcid": "", "doi": "",
            "pubTypeList": {"pubType": item.get("pubtype", [])},
            "abstractText": ""
        })
    return records, {
        "source": "PubMed E-utilities", "query": PUBMED_QUERY,
        "idlist": ids, "result": result
    }

def main():
    out = Path("data/processed/literature")
    raw_dir = Path("data/raw/europe_pmc")
    out.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    source = "Europe PMC"
    try:
        results, payload = search_europe_pmc()
    except Exception as epmc_error:
        print(f"Europe PMC no respondió correctamente: {epmc_error}")
        print("Activando búsqueda alternativa en PubMed.")
        results, payload = search_pubmed()
        source = "PubMed E-utilities (fallback)"
        (raw_dir / "cyp2c9_sequencing_search_pubmed.json").write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )
    else:
        (raw_dir / "cyp2c9_sequencing_search.json").write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )

    records = []
    seen = set()
    for item in results:
        pmid = clean(item.get("pmid"))
        doi = clean(item.get("doi")).lower()
        title = clean(item.get("title"))
        key = ("doi:" + doi) if doi else (("pmid:" + pmid) if pmid else
               "title:" + title.lower())
        if key in seen:
            continue
        seen.add(key)

        abstract = clean(item.get("abstractText"))
        text = (title + " " + abstract).lower()
        methods = []
        for term in ["next-generation sequencing", "next generation sequencing",
                     "ngs", "sanger sequencing", "pyrosequencing",
                     "whole-genome sequencing", "whole exome sequencing",
                     "targeted sequencing", "genotyping", "haplotype"]:
            if term in text:
                methods.append(term)

        records.append({
            "title": title,
            "authors": clean(item.get("authorString")),
            "year": clean(item.get("pubYear")),
            "journal": clean(item.get("journalTitle")),
            "pmid": pmid,
            "pmcid": clean(item.get("pmcid")),
            "doi": clean(item.get("doi")),
            "publication_type": clean(item.get("pubTypeList", {}).get("pubType", [])),
            "methods_terms_detected": "; ".join(methods),
            "abstract": abstract,
            "source": source,
            "record_url": (
                "https://pubmed.ncbi.nlm.nih.gov/" + pmid + "/" if pmid else
                ("https://doi.org/" + doi if doi else "")
            )
        })

    records.sort(key=lambda x: (x["year"], x["title"]), reverse=True)
    fields = [
        "title","authors","year","journal","pmid","pmcid","doi",
        "publication_type","methods_terms_detected","abstract","source","record_url"
    ]
    with (out / "cyp2c9_sequencing_literature.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)

    summary = {
        "database_used": source,
        "primary_database": "Europe PMC",
        "fallback_database": "PubMed E-utilities",
        "query": QUERY,
        "retrieved_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "records_returned": len(results),
        "records_after_deduplication": len(records),
        "deduplication_priority": "DOI, then PMID, then normalized title",
        "method_terms_are": "automatically detected text terms, not manually verified method assignments",
        "limitations": [
            "Search sensitivity depends on title/abstract/indexing terms.",
            "PubMed fallback metadata may not include abstracts or DOI values.",
            "A detected method term is not proof that the study used that method experimentally.",
            "Full-text review is required to extract sample sizes, platforms, coverage and validation metrics.",
            "The dataset is a literature discovery corpus, not a systematic review or meta-analysis."
        ]
    }
    (out / "cyp2c9_sequencing_literature_metadata.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    if not records:
        raise RuntimeError("La búsqueda terminó correctamente, pero no devolvió registros.")

if __name__ == "__main__":
    main()
