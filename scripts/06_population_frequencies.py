#!/usr/bin/env python3
"""Retrieve population-level allele frequency observations for selected CYP2C9 variants.

Source: Ensembl REST variation endpoint with population data enabled.
Frequencies are preserved per population and allele; they are not pooled or averaged.
"""
import csv
import json
import time
from pathlib import Path

import requests

BASE = "https://rest.ensembl.org/variation/human"
VARIANTS = {
    "rs1799853": {"star_allele": "*2", "known_variant_allele": "T"},
    "rs1057910": {"star_allele": "*3", "known_variant_allele": "C"},
    "rs28371686": {"star_allele": "*5", "known_variant_allele": ""},
    "rs9332131": {"star_allele": "*6", "known_variant_allele": ""},
    "rs7900194": {"star_allele": "*8", "known_variant_allele": ""},
    "rs28371685": {"star_allele": "*11", "known_variant_allele": ""},
}
SESSION = requests.Session()
SESSION.headers.update({
    "Accept": "application/json",
    "User-Agent": "CYP2C9-Research/1.0 (population-frequency analysis)"
})

def fetch_variant(rsid):
    url = f"{BASE}/{rsid}"
    response = SESSION.get(
        url, params={"pops": 1}, timeout=(20, 60)
    )
    response.raise_for_status()
    return response.json()

def main():
    raw_dir = Path("data/raw/ensembl")
    out_dir = Path("data/processed/population_frequencies")
    raw_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    observations = []
    metadata = {
        "source": "Ensembl REST API",
        "endpoint": "/variation/human/{rsID}?pops=1",
        "retrieved_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "variants_requested": len(VARIANTS),
        "variants_retrieved": 0,
        "frequency_unit": "proportion (0-1)",
        "interpretation": (
            "Each row is one allele-frequency observation for one variant "
            "in one named population/source. Do not pool or average across "
            "populations or datasets. Population labels and sample denominators "
            "are retained as returned by Ensembl when available."
        ),
        "limitations": [
            "Ensembl aggregates frequency observations from multiple source datasets.",
            "Not every variant has observations for every population.",
            "A missing row means no observation was returned, not a frequency of zero.",
            "The rsID-to-star-allele mapping is a project annotation and should be "
            "checked against current PharmVar/ClinVar records before publication.",
            "These data are not a representative estimate for Mexico or any local population."
        ],
        "variants": {}
    }

    for rsid, info in VARIANTS.items():
        try:
            payload = fetch_variant(rsid)
            (raw_dir / f"{rsid}_ensembl_variation.json").write_text(
                json.dumps(payload, indent=2), encoding="utf-8"
            )
            metadata["variants_retrieved"] += 1
            metadata["variants"][rsid] = {
                "star_allele": info["star_allele"],
                "known_variant_allele": info["known_variant_allele"],
                "minor_allele": payload.get("minor_allele"),
                "MAF": payload.get("MAF"),
                "evidence": payload.get("evidence", []),
                "most_severe_consequence": payload.get("most_severe_consequence"),
                "mappings": payload.get("mappings", [])
            }
            for pop in payload.get("populations", []):
                observations.append({
                    "star_allele": info["star_allele"],
                    "rsid": rsid,
                    "population": pop.get("population", ""),
                    "allele": pop.get("allele", ""),
                    "frequency": pop.get("frequency", ""),
                    "frequency_percent": (
                        float(pop["frequency"]) * 100
                        if isinstance(pop.get("frequency"), (int, float)) else ""
                    ),
                    "count": pop.get("count", ""),
                    "dataset_source": pop.get("population", "").split(":")[0],
                    "known_variant_allele": info["known_variant_allele"],
                    "is_known_variant_allele": (
                        str(pop.get("allele", "")).upper() == info["known_variant_allele"]
                        if info["known_variant_allele"] else ""
                    )
                })
            print(f"{rsid}: {len(payload.get('populations', []))} population observations")
        except requests.RequestException as exc:
            metadata["variants"][rsid] = {
                "star_allele": info["star_allele"],
                "error": str(exc)
            }
            print(f"WARNING: {rsid} could not be retrieved: {exc}")
        time.sleep(0.35)

    fields = [
        "star_allele", "rsid", "population", "allele", "frequency",
        "frequency_percent", "count", "dataset_source",
        "known_variant_allele", "is_known_variant_allele"
    ]
    csv_path = out_dir / "cyp2c9_population_allele_frequencies.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(observations)

    metadata["population_observations"] = len(observations)
    metadata["variants_with_frequency_observations"] = len({
        row["rsid"] for row in observations
    })
    (out_dir / "cyp2c9_population_frequency_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))
    if not observations:
        raise RuntimeError("No population frequency observations were returned.")

if __name__ == "__main__":
    main()
