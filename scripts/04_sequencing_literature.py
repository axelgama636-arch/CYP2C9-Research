#!/usr/bin/env python3
"""Retrieve and curate CYP2C9 sequencing/genotyping literature from Europe PMC.

The script saves raw API records, a deduplicated CSV, and JSON metadata.
It does not infer study quality or clinical validity automatically.
"""
import csv
import json
import re
import time
from pathlib import Path

import requests

API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
QUERY = (
    'CYP2C9 AND (sequencing OR "next generation sequencing" OR NGS '
    'OR Sanger OR pyrosequencing OR genotyping OR haplotype OR allele)'
)

def clean(value):
    if value is None:
        return ""
    if isinstance(value, list):
        return "; ".join(str(x) for x in value)
    return re.sub(r"\s+", " ", str(value)).strip()

def main():
    out = Path("data/processed/literature")
    raw_dir = Path("data/raw/europe_pmc")
    out.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    params = {
        "query": QUERY,
        "format": "json",
        "pageSize": 1000,
        "resultType": "core",
        "sort": "RELEVANCE"
    }
    response = requests.get(API, params=params, timeout=120)
    response.raise_for_status()
    payload = response.json()
    results = payload.get("resultList", {}).get("result", [])

    (raw_dir / "cyp2c9_sequencing_search.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )

    records = []
    seen = set()
    for item in results:
        pmid = clean(item.get("pmid"))
        doi = clean(item.get("doi")).lower()
        key = ("doi:" + doi) if doi else (("pmid:" + pmid) if pmid else
               "title:" + clean(item.get("title")).lower())
        if key in seen:
            continue
        seen.add(key)

        title = clean(item.get("title"))
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
            "source": "Europe PMC",
            "record_url": (
                "https://europepmc.org/article/MED/" + pmid if pmid else
                ("https://doi.org/" + doi if doi else "")
            )
        })

    records.sort(key=lambda x: (x["year"], x["title"]), reverse=True)
    fields = list(records[0].keys()) if records else [
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
        "database": "Europe PMC",
        "query": QUERY,
        "retrieved_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "api_total_hits": payload.get("hitCount"),
        "records_returned": len(results),
        "records_after_deduplication": len(records),
        "deduplication_priority": "DOI, then PMID, then normalized title",
        "method_terms_are": "automatically detected text terms, not manually verified method assignments",
        "limitations": [
            "Search sensitivity depends on title/abstract/indexing terms.",
            "A detected method term is not proof that the study used that method experimentally.",
            "Full-text review is required to extract sample sizes, platforms, coverage and validation metrics.",
            "The dataset is a literature search corpus, not a systematic review or meta-analysis."
        ]
    }
    (out / "cyp2c9_sequencing_literature_metadata.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
