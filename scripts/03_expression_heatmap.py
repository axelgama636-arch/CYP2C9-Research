#!/usr/bin/env python3
"""Build a CYP2C9-centered 50-gene x 6-tissue HPA RNA expression heatmap.

Uses the official Human Protein Atlas gene-level RNA tissue matrix.
Missing measurements remain missing (never imputed). Z-scores are calculated
across the six tissues within each gene after log2(nTPM + 1).
"""
import argparse
import io
import json
import time
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

HPA_URL = "https://www.proteinatlas.org/download/tsv/rna_tissue_hpa.tsv.zip"
HPA_VERSION = "26.0"
TISSUES = ["liver", "kidney", "small intestine", "lung", "cerebral cortex", "heart muscle"]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", default="data/processed/interactome/nodes.csv")
    parser.add_argument("--out", default="data/processed/expression")
    parser.add_argument("--hpa-url", default=HPA_URL)
    args = parser.parse_args()

    nodes_path = Path(args.nodes)
    if not nodes_path.exists():
        raise SystemExit(f"Missing {nodes_path}. Run scripts/02_interactome_string.py first.")
    node_table = pd.read_csv(nodes_path)
    gene_col = "gene" if "gene" in node_table.columns else node_table.columns[0]
    genes = node_table[gene_col].dropna().astype(str).drop_duplicates().tolist()
    if "CYP2C9" not in genes:
        genes.insert(0, "CYP2C9")
    genes = genes[:50]

    response = requests.get(args.hpa_url, timeout=180)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        tsv_files = [n for n in archive.namelist() if n.lower().endswith((".tsv", ".txt"))]
        if not tsv_files:
            raise RuntimeError("HPA ZIP contains no TSV/TXT data file.")
        with archive.open(tsv_files[0]) as handle:
            hpa = pd.read_csv(handle, sep="\t", low_memory=False)

    required = {"Gene name", "Tissue", "nTPM"}
    if not required.issubset(hpa.columns):
        raise RuntimeError(f"Unexpected HPA columns. Required {sorted(required)}; got {list(hpa.columns)}")

    selected = hpa[
        hpa["Gene name"].astype(str).isin(genes)
        & hpa["Tissue"].astype(str).isin(TISSUES)
    ][["Gene name", "Tissue", "nTPM"]].copy()
    selected["nTPM"] = pd.to_numeric(selected["nTPM"], errors="coerce")
    pivot = selected.pivot_table(
        index="Gene name", columns="Tissue", values="nTPM", aggfunc="mean"
    ).reindex(index=genes, columns=TISSUES)
    pivot.index.name = "Gene"
    pivot.columns.name = "Tissue"

    log_matrix = np.log2(pivot + 1)
    means = log_matrix.mean(axis=1, skipna=True)
    stds = log_matrix.std(axis=1, skipna=True, ddof=0)
    z_matrix = log_matrix.sub(means, axis=0).div(stds.replace(0, np.nan), axis=0)
    z_matrix = z_matrix.mask(stds.eq(0), 0)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    pivot.to_csv(out / "cyp2c9_network_50x6_ntpm.csv", float_format="%.6g")
    log_matrix.to_csv(out / "cyp2c9_network_50x6_log2_ntpm.csv", float_format="%.6g")
    z_matrix.to_csv(out / "cyp2c9_network_50x6_zscore.csv", float_format="%.6g")

    present = pivot.notna().any(axis=1)
    coverage = {
        "requested_gene_count": len(genes),
        "genes_with_any_tissue_value": int(present.sum()),
        "genes_without_values": present.index[~present].tolist(),
        "missing_cells": int(pivot.isna().sum().sum()),
        "total_cells": int(pivot.shape[0] * pivot.shape[1]),
        "tissues": TISSUES,
    }
    (out / "coverage.json").write_text(json.dumps(coverage, indent=2), encoding="utf-8")

    plt.rcParams["font.family"] = "Cambria"
    plt.rcParams["font.size"] = 9
    fig, ax = plt.subplots(figsize=(10, max(10, len(genes) * 0.23)))
    image = ax.imshow(z_matrix.to_numpy(dtype=float), aspect="auto",
                      interpolation="nearest", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(TISSUES)))
    ax.set_xticklabels(TISSUES, rotation=35, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=6.5)
    ax.set_title(f"CYP2C9-centered expression profile | HPA v{HPA_VERSION}\nlog2(nTPM + 1), row-wise Z-score")
    ax.set_xlabel("Tissue")
    ax.set_ylabel("Genes in CYP2C9 STRING network")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.03)
    colorbar.set_label("Z-score per gene")
    fig.tight_layout()
    fig.savefig(out / "CYP2C9_expression_heatmap_50x6.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    metadata = {
        "analysis": "CYP2C9-centered tissue expression heatmap",
        "source": "Human Protein Atlas",
        "hpa_version": HPA_VERSION,
        "hpa_release_url": args.hpa_url,
        "organism": "Homo sapiens",
        "gene_set_source": "STRING network nodes CSV",
        "requested_gene_count": len(genes),
        "tissues": TISSUES,
        "source_measure": "nTPM",
        "transform": "log2(nTPM + 1)",
        "scaling": "row-wise Z-score across selected tissues; ddof=0",
        "duplicate_aggregation": "mean by gene and tissue",
        "missing_data": "retained as NA; no imputation",
        "retrieved_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "coverage": coverage,
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))

if __name__ == "__main__":
    main()
