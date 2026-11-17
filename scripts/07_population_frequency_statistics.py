#!/usr/bin/env python3
"""Descriptive population-frequency analysis for CYP2C9 defining variants.

Uses Ensembl population observations already retrieved in Analysis 05.
No cross-dataset pooling and no inferential tests are performed because
allele counts/denominators are not consistently supplied in the input.
"""
import csv
import math
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path("data/processed/population_frequencies")
INPUT = BASE / "cyp2c9_population_allele_frequencies.csv"
ALLELES = {
    "rs1799853": ("*2", "T", "SNV"),
    "rs1057910": ("*3", "C", "SNV"),
    "rs28371686": ("*5", "G", "SNV"),
    "rs9332131": ("*6", "A", "deletion; Ensembl representation A vs AA"),
    "rs7900194": ("*8", "A", "SNV"),
    "rs28371685": ("*11", "T", "SNV"),
}
# 1000 Genomes phase 3 population codes (individual populations only).
POP_LABELS = {
    "ACB":"African Caribbean in Barbados", "ASW":"African ancestry in SW USA",
    "BEB":"Bengali in Bangladesh", "CDX":"Chinese Dai", "CEU":"Utah residents with NW European ancestry",
    "CHB":"Han Chinese in Beijing", "CHS":"Southern Han Chinese", "CLM":"Colombian in Medellin",
    "ESN":"Esan in Nigeria", "FIN":"Finnish in Finland", "GBR":"British in England/Scotland",
    "GIH":"Gujarati Indian in Houston", "GWD":"Gambian in Western Division", "IBS":"Iberian in Spain",
    "ITU":"Indian Telugu in UK", "JPT":"Japanese in Tokyo", "KHV":"Kinh in Ho Chi Minh City",
    "LWK":"Luhya in Webuye, Kenya", "MSL":"Mende in Sierra Leone", "MXL":"Mexican ancestry in Los Angeles",
    "PEL":"Peruvian in Lima", "PJL":"Punjabi in Lahore", "PUR":"Puerto Rican in Puerto Rico",
    "STU":"Sri Lankan Tamil in UK", "TSI":"Toscani in Italy", "YRI":"Yoruba in Ibadan"
}
AGGREGATES = {"ALL","AFR","AMR","EAS","EUR","SAS"}

def read_csv():
    with INPUT.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def write_csv(path, rows, fields):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def main():
    rows = read_csv()
    harmonized = []
    for r in rows:
        rsid = r["rsid"]
        allele = ALLELES.get(rsid, ("", "", ""))[1]
        r["known_variant_allele_verified"] = allele
        r["variant_allele_frequency"] = (
            r["frequency"] if allele and r["allele"].upper() == allele else ""
        )
        r["variant_allele_match"] = (
            str(r["allele"]).upper() == allele if allele else ""
        )
        r["variant_definition"] = ALLELES.get(rsid, ("","",""))[2]
        harmonized.append(r)

    write_csv(BASE / "cyp2c9_population_frequencies_harmonized.csv", harmonized,
              list(harmonized[0].keys()))

    # Only individual 1000 Genomes phase-3 populations; aggregate super-populations
    # and ALL are excluded to avoid treating nested groups as independent samples.
    pop_rows = []
    for r in harmonized:
        pop = r["population"]
        if not pop.startswith("1000GENOMES:phase_3:"):
            continue
        code = pop.rsplit(":", 1)[-1]
        if code in AGGREGATES or code not in POP_LABELS:
            continue
        if r["variant_allele_frequency"] == "":
            continue
        pop_rows.append({
            "star_allele": r["star_allele"], "rsid": r["rsid"],
            "variant_allele": r["known_variant_allele_verified"],
            "population_code": code, "population_description": POP_LABELS[code],
            "frequency": float(r["variant_allele_frequency"]),
            "frequency_percent": float(r["variant_allele_frequency"]) * 100,
            "dataset": "1000 Genomes phase 3"
        })
    write_csv(BASE / "cyp2c9_1000genomes_variant_allele_frequencies.csv",
              pop_rows, list(pop_rows[0].keys()))

    summary = []
    for rsid, (star, allele, definition) in ALLELES.items():
        vals = [x["frequency_percent"] for x in pop_rows if x["rsid"] == rsid]
        if not vals:
            continue
        summary.append({
            "star_allele": star, "rsid": rsid, "defining_allele": allele,
            "variant_type": definition, "n_1000G_individual_populations": len(vals),
            "min_frequency_percent": min(vals), "max_frequency_percent": max(vals),
            "unweighted_mean_across_populations_percent": statistics.mean(vals),
            "median_across_populations_percent": statistics.median(vals),
            "interpretation": "Descriptive only; populations are not sample-size weighted"
        })
    write_csv(BASE / "cyp2c9_1000genomes_descriptive_summary.csv",
              summary, list(summary[0].keys()))

    # Targeted AMR superpopulation and MXL individual population, retained separately.
    targeted = []
    for r in harmonized:
        if not r["population"].startswith("1000GENOMES:phase_3:"):
            continue
        code = r["population"].rsplit(":", 1)[-1]
        if code not in {"AMR", "MXL", "CLM", "PEL", "PUR"}:
            continue
        if r["variant_allele_frequency"] == "":
            continue
        targeted.append({
            "star_allele": r["star_allele"], "rsid": r["rsid"],
            "variant_allele": r["known_variant_allele_verified"],
            "population_code": code, "population_description": POP_LABELS.get(code, "1000G AMR aggregate"),
            "frequency_percent": float(r["variant_allele_frequency"]) * 100,
            "note": "AMR is a superpopulation aggregate; MXL/CLM/PEL/PUR are distinct cohorts"
        })
    write_csv(BASE / "cyp2c9_1000genomes_americas_focus.csv", targeted,
              list(targeted[0].keys()))

    # Heatmap across individual 1000G populations. Missing observations remain blank.
    codes = [c for c in POP_LABELS if any(x["population_code"] == c for x in pop_rows)]
    stars = [v[0] for v in ALLELES.values()]
    matrix = np.full((len(stars), len(codes)), np.nan)
    for i, star in enumerate(stars):
        for j, code in enumerate(codes):
            vals = [x["frequency_percent"] for x in pop_rows
                    if x["star_allele"] == star and x["population_code"] == code]
            if vals:
                matrix[i, j] = vals[0]
    fig_w = max(12, len(codes) * 0.48)
    fig, ax = plt.subplots(figsize=(fig_w, 5.6))
    masked = np.ma.masked_invalid(matrix)
    im = ax.imshow(masked, aspect="auto", cmap="magma", vmin=0,
                   vmax=max(5, float(np.nanmax(matrix))))
    ax.set_yticks(range(len(stars)), stars)
    ax.set_xticks(range(len(codes)), codes, rotation=55, ha="right", fontsize=8)
    ax.set_xlabel("Población individual (1000 Genomes, fase 3)")
    ax.set_ylabel("Alelo estrella / variante definitoria")
    ax.set_title("Frecuencia del alelo definitorio de CYP2C9 por población")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Frecuencia alélica (%)")
    fig.tight_layout()
    fig.savefig(BASE / "CYP2C9_1000Genomes_frequency_heatmap.png", dpi=300)
    plt.close(fig)

    print(f"Observaciones originales: {len(rows)}")
    print(f"Filas armonizadas: {len(harmonized)}")
    print(f"Observaciones 1000G individuales: {len(pop_rows)}")
    print(f"Variantes con resumen: {len(summary)}")
    print("Resultados descriptivos guardados; no se ejecutaron pruebas inferenciales.")

if __name__ == "__main__":
    main()
