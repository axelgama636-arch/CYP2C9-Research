
"""
Proyecto: CYP2C9-Research
Módulo: Recuperación bibliográfica desde PubMed

Recupera artículos científicos y genera:
1. CSV con información bibliográfica.
2. JSON con metadatos de búsqueda.

Fuente: NCBI E-utilities.
"""

import argparse
import csv
import json
import time
from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET

import requests

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

DEFAULT_QUERY = (
    '"CYP2C9"[Title/Abstract] OR '
    '"cytochrome P450 2C9"[Title/Abstract] OR '
    "rs1799853 OR rs1057910"
)


def obtener_texto(elemento, ruta):
    encontrado = elemento.find(ruta)
    if encontrado is None:
        return ""
    return "".join(encontrado.itertext()).strip()


def buscar_pmids(query, limite, email):
    parametros = {
        "db": "pubmed",
        "term": query,
        "retmax": limite,
        "retmode": "json",
        "sort": "relevance",
        "tool": "CYP2C9_Research"
    }

    if email:
        parametros["email"] = email

    respuesta = requests.get(
        BASE_URL + "esearch.fcgi",
        params=parametros,
        timeout=60
    )
    respuesta.raise_for_status()

    datos = respuesta.json()["esearchresult"]
    return datos["idlist"], int(datos["count"])


def recuperar_articulos(pmids, email):
    if not pmids:
        return []

    parametros = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "tool": "CYP2C9_Research"
    }

    if email:
        parametros["email"] = email

    respuesta = requests.get(
        BASE_URL + "efetch.fcgi",
        params=parametros,
        timeout=120
    )
    respuesta.raise_for_status()

    raiz = ET.fromstring(respuesta.content)
    articulos = []

    for articulo in raiz.findall(".//PubmedArticle"):
        medline = articulo.find("MedlineCitation")
        if medline is None:
            continue

        pmid = obtener_texto(medline, ".//PMID")
        datos = medline.find("Article")

        if datos is None:
            continue

        titulo = obtener_texto(datos, "ArticleTitle")
        revista = obtener_texto(datos, ".//Journal/Title")
        anio = obtener_texto(datos, ".//Journal/JournalIssue/PubDate/Year")

        if not anio:
            anio = obtener_texto(
                datos, ".//Journal/JournalIssue/PubDate/MedlineDate"
            )[:4]

        autores = []
        for autor in datos.findall(".//AuthorList/Author"):
            apellido = obtener_texto(autor, "LastName")
            nombre = obtener_texto(autor, "ForeName")
            colectivo = obtener_texto(autor, "CollectiveName")

            nombre_completo = (
                colectivo or " ".join(x for x in [nombre, apellido] if x)
            )
            if nombre_completo:
                autores.append(nombre_completo)

        doi = ""
        for identificador in articulo.findall(
            ".//PubmedData/ArticleIdList/ArticleId"
        ):
            if identificador.attrib.get("IdType") == "doi":
                doi = (identificador.text or "").strip()

        resumen = " ".join(
            "".join(parrafo.itertext()).strip()
            for parrafo in datos.findall(".//Abstract/AbstractText")
        )

        articulos.append({
            "PMID": pmid,
            "DOI": doi,
            "Titulo": titulo,
            "Autores": "; ".join(autores),
            "Anio": anio,
            "Revista": revista,
            "Resumen": resumen,
            "URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        })

    return articulos


def main():
    parser = argparse.ArgumentParser(
        description="Recuperación bibliográfica de CYP2C9 desde PubMed"
    )

    parser.add_argument("--query", default=DEFAULT_QUERY)
    parser.add_argument("--limit", type=int, default=500)
    parser.add_argument("--email", default="")
    parser.add_argument("--output", default="data/raw")

    args = parser.parse_args()

    if args.limit < 1 or args.limit > 10000:
        raise ValueError("El límite debe estar entre 1 y 10000.")

    carpeta = Path(args.output)
    carpeta.mkdir(parents=True, exist_ok=True)

    print("Buscando publicaciones en PubMed...")

    pmids, total = buscar_pmids(
        args.query, args.limit, args.email
    )

    print(f"Resultados encontrados: {total}")
    print(f"PMID recuperados: {len(pmids)}")

    articulos = recuperar_articulos(pmids, args.email)

    fecha = datetime.now(timezone.utc).isoformat()
    archivo_csv = carpeta / "cyp2c9_pubmed_articles.csv"
    archivo_json = carpeta / "pubmed_search_metadata.json"

    campos = [
        "PMID", "DOI", "Titulo", "Autores",
        "Anio", "Revista", "Resumen", "URL"
    ]

    with archivo_csv.open(
        "w", newline="", encoding="utf-8-sig"
    ) as archivo:
        escritor = csv.DictWriter(
            archivo, fieldnames=campos
        )
        escritor.writeheader()
        escritor.writerows(articulos)

    metadatos = {
        "gene": "CYP2C9",
        "database": "PubMed",
        "source": "NCBI E-utilities",
        "query": args.query,
        "retrieved_at_utc": fecha,
        "total_matching_records": total,
        "pmids_requested": len(pmids),
        "articles_recovered": len(articulos),
        "limit": args.limit
    }

    with archivo_json.open("w", encoding="utf-8") as archivo:
        json.dump(
            metadatos, archivo,
            ensure_ascii=False, indent=4
        )

    print("\nProceso terminado.")
    print(f"Artículos recuperados: {len(articulos)}")
    print(f"CSV: {archivo_csv}")
    print(f"Metadatos: {archivo_json}")


if __name__ == "__main__":
    main()
