#!/usr/bin/env python3
"""Retrieve a reproducible human CYP2C9 functional-association network from STRING.

STRING edges are functional associations and are not necessarily direct
physical protein-protein binding events.
"""
import argparse
import csv
import json
import time
from pathlib import Path

import requests

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=int, default=50)
    parser.add_argument("--score", type=int, default=700)
    parser.add_argument("--out", default="data/processed/interactome")
    args = parser.parse_args()

    if not 2 <= args.nodes <= 500:
        raise SystemExit("--nodes must be between 2 and 500")
    if not 150 <= args.score <= 1000:
        raise SystemExit("--score must be between 150 and 1000")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    params = {
        "identifiers": "CYP2C9",
        "species": 9606,
        "required_score": args.score,
        "add_nodes": args.nodes - 1,
        "network_type": "functional",
        "caller_identity": "CYP2C9_Research"
    }
    response = requests.get(
        "https://string-db.org/api/json/network",
        params=params,
        timeout=120
    )
    response.raise_for_status()
    data = response.json()
    if not data:
        raise RuntimeError("STRING returned no network edges.")

    nodes = {}
    edges = []
    for item in data:
        a, b = item["preferredName_A"], item["preferredName_B"]
        nodes[a] = {"gene": a, "string_id": item.get("stringId_A")}
        nodes[b] = {"gene": b, "string_id": item.get("stringId_B")}
        edges.append({
            "source": a,
            "target": b,
            "combined_score": item.get("score"),
            "database": "STRING",
            "evidence_type": "functional_association",
            "physical_binding_confirmed": False
        })

    node_rows = sorted(nodes.values(), key=lambda x: x["gene"])
    edge_rows = sorted(edges, key=lambda x: (x["source"], x["target"]))
    for filename, content in (("nodes.json", node_rows), ("edges.json", edge_rows)):
        (out / filename).write_text(json.dumps(content, indent=2), encoding="utf-8")

    with (out / "nodes.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "string_id"])
        writer.writeheader()
        writer.writerows(node_rows)

    edge_fields = ["source", "target", "combined_score", "database",
                   "evidence_type", "physical_binding_confirmed"]
    with (out / "edges.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=edge_fields)
        writer.writeheader()
        writer.writerows(edge_rows)

    metadata = {
        "gene": "CYP2C9",
        "taxon": 9606,
        "database": "STRING",
        "network_type": "functional",
        "minimum_combined_score": args.score,
        "requested_nodes_including_query": args.nodes,
        "observed_nodes": len(node_rows),
        "observed_edges": len(edge_rows),
        "retrieved_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "query_parameters": params,
        "interpretation": "Functional association does not necessarily mean direct physical binding."
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    if len(node_rows) != args.nodes:
        print("NOTE: STRING returned a different number of nodes; do not add invented nodes.")

if __name__ == "__main__":
    main()
