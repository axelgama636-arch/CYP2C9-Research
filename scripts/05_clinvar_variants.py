#!/usr/bin/env python3
"""Build a reproducible CYP2C9 variant catalogue from NCBI ClinVar.

ClinVar assertions are reported as submitted/reviewed by the database.
This script does not infer star alleles, pathogenicity, or population frequencies.
"""
import csv
import json
import re
import time
from pathlib import Path

import requests

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
QUERY = "CYP2C9[gene]"
BATCH_SIZE = 100
RETRYABLE = {429, 500, 502, 503, 504}
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "CYP2C9-Research/1.0 (variant catalogue)"})
STAR_ALLELES = {
    "rs1799853": "CYP2C9*2",
    "rs1057910": "CYP2C9*3",
}

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
            delay = min(60, 3 * (2 ** attempt))
            print(f"NCBI temporalmente no disponible ({exc}); reintento en {delay}s.")
            time.sleep(delay)
    raise RuntimeError(f"NCBI no respondió tras {attempts} intentos: {last_error}")

def fetch_uids():
    first = get_json(ESEARCH, {
        "db": "clinvar", "term": QUERY, "retmode": "json",
        "retmax": 0, "usehistory": "n"
    })
    count = int(first["esearchresult"]["count"])
    ids = []
    for start in range(0, count, 500):
        page = get_json(ESEARCH, {
            "db": "clinvar", "term": QUERY, "retmode": "json",
            "retstart": start, "retmax": min(500, count - start)
        })
        ids.extend(page["esearchresult"].get("idlist", []))
        time.sleep(0.34)
    return count, ids

def fetch_summaries(uids):
    summaries = []
    for start in range(0, len(uids), BATCH_SIZE):
        batch = uids[start:start + BATCH_SIZE]
        payload = get_json(ESUMMARY, {
            "db": "clinvar", "id": ",".join(batch), "retmode": "json"
        })
        result = payload.get("result", {})
        for uid in result.get("uids", []):
            if uid in result:
                summaries.append(result[uid])
        time.sleep(0.34)
    return summaries

def clean(value):
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return re.sub(r"\\s+", " ", str(value)).strip()

def extract_rsids(item):
    serialized = json.dumps(item, ensure_ascii=False)
    return sorted(set(re.findall(r"\\brs[0-9]+\\b", serialized)))

def main():
    raw_dir = Path("data/raw/clinvar")
    out_dir = Path("data/processed/variants")
    raw_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    total, uids = fetch_uids()
    summaries = fetch_summaries(uids)
    raw = {
        "database": "NCBI ClinVar",
        "query": QUERY,
        "retrieved_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "search_count": total,
        "uids_requested": len(uids),
        "records_retrieved": len(summaries),
        "records": summaries,
    }
    (raw_dir / "cyp2c9_clinvar_esummary.json").write_text(
        json.dumps(raw, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    rows = []
    for item in summaries:
        rsids = extract_rsids(item)
        star = "; ".join(STAR_ALLELES[x] for x in rsids if x in STAR_ALLELES)
        genes = item.get("genes", [])
        gene_names = sorted({
            g.get("symbol", "") for g in genes
            if isinstance(g, dict) and g.get("symbol")
        })
        clinical = item.get("clinical_significance", {})
        rows.append({
            "uid": clean(item.get("uid")),
            "accession": clean(item.get("accession")),
            "title": clean(item.get("title")),
            "object_type": clean(item.get("obj_type")),
            "gene_symbols": "; ".join(gene_names),
            "rsids_detected_in_record": "; ".join(rsids),
            "known_star_allele_mapping": star,
            "clinical_significance_description": clean(
                clinical.get("description") if isinstance(clinical, dict) else ""
            ),
            "review_status": clean(
                clinical.get("review_status") if isinstance(clinical, dict) else ""
            ),
            "last_evaluated": clean(
                clinical.get("last_evaluated") if isinstance(clinical, dict) else ""
            ),
            "variation_set_json": clean(item.get("variation_set", [])),
            "record_url": (
                "https://www.ncbi.nlm.nih.gov/clinvar/" +
                clean(item.get("accession")) if item.get("accession") else
                "https://www.ncbi.nlm.nih.gov/clinvar/?term=" + clean(item.get("uid"))
            ),
        })

    rows.sort(key=lambda x: (x["known_star_allele_mapping"] == "", x["title"]))
    fields = [
        "uid", "accession", "title", "object_type", "gene_symbols",
        "rsids_detected_in_record", "known_star_allele_mapping",
        "clinical_significance_description", "review_status", "last_evaluated",
        "variation_set_json", "record_url"
    ]
    csv_path = out_dir / "cyp2c9_clinvar_variants.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    known = {
        "CYP2C9*2": {"rsid": "rs1799853", "HGVS_coding": "NM_000771.4:c.430C>T"},
        "CYP2C9*3": {"rsid": "rs1057910", "HGVS_coding": "NM_000771.4:c.1075A>C"},
    }
    metadata = {
        "database": "NCBI ClinVar",
        "query": QUERY,
        "retrieved_at_utc": raw["retrieved_at_utc"],
        "search_count": total,
        "records_retrieved": len(summaries),
        "csv_records": len(rows),
        "known_variant_reference_points": known,
        "interpretation_rules": [
            "ClinVar clinical significance is preserved as reported; it is not reclassified by this script.",
            "Review status and last evaluated date must be considered when interpreting each assertion.",
            "Star-allele labels are assigned only for the two explicitly curated rsID mappings above.",
            "A ClinVar record is not necessarily an independent study or a population-frequency observation.",
            "No allele frequency is calculated because this pipeline does not retrieve population genotype denominators."
        ],
        "limitations": [
            "The search is gene-index based and may include records with differing molecular scopes.",
            "rsIDs are detected from the serialized ClinVar summary; inspect the raw JSON for complex records.",
            "Full HGVS normalization and PharmVar haplotype assignment require additional transcript and haplotype review.",
            "This catalogue is descriptive and is not a clinical interpretation."
        ]
    }
    (out_dir / "cyp2c9_clinvar_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2, ensure_ascii=False))
    if not summaries:
        raise RuntimeError("ClinVar no devolvió registros para CYP2C9.")

if __name__ == "__main__":
    main()
