# Análisis 02. Perfil de expresión tisular de CYP2C9 y genes asociados

## Objetivo
Cuantificar y visualizar la expresión de CYP2C9 y de los genes incluidos en su red STRING en seis tejidos humanos, usando datos RNA-seq del Human Protein Atlas (HPA).

## Fuente de datos
Human Protein Atlas, archivo oficial de expresión génica RNA por tejido:
- Archivo: rna_tissue_hpa.tsv.zip
- Versión documentada: HPA 26.0; Ensembl 109.
- Medida primaria: normalized transcript expression (nTPM).
- Descarga oficial: https://www.proteinatlas.org/download/tsv/rna_tissue_hpa.tsv.zip
- Descripción y cobertura: https://www.proteinatlas.org/humanproteome/tissue/data

El archivo contiene identificador Ensembl, símbolo génico, tejido, TPM, pTPM y nTPM. La expresión génica resumida se presenta para 40 tejidos HPA.

## Selección de tejidos
Se seleccionan seis tejidos relevantes para contextualizar la expresión de genes relacionados con metabolismo y disposición de compuestos:
1. Hígado (liver)
2. Riñón (kidney)
3. Intestino delgado (small intestine)
4. Pulmón (lung)
5. Corteza cerebral (cerebral cortex)
6. Corazón (heart muscle)

Esta selección es un subconjunto analítico deliberado de los tejidos disponibles, no una afirmación de que CYP2C9 tenga la misma relevancia funcional en todos ellos.

## Conjunto de genes
Los símbolos se obtienen del archivo de nodos producido por el análisis interactómico de STRING. Se analizan CYP2C9 y los vecinos recuperados, hasta un máximo de 50 genes en total. No se completará el conjunto con genes inventados. Si STRING devuelve menos de 50 nodos, el heatmap tendrá ese número real de filas.

## Procesamiento
1. Descargar el archivo oficial HPA RNA tissue gene data.
2. Filtrar por los símbolos de genes de la red y los seis tejidos definidos.
3. Si existieran filas repetidas para una combinación gen-tejido, resumirlas mediante la media aritmética.
4. Construir una matriz nTPM genes × tejidos.
5. Transformar cada celda mediante log2(nTPM + 1).
6. Calcular Z-score por fila (cada gen, comparando sus seis tejidos), con desviación estándar poblacional (ddof=0).
7. Conservar los valores ausentes como NA; no imputar expresión cero.
8. Generar matrices CSV, cobertura y heatmap PNG a 300 dpi.

## Interpretación y límites
- nTPM es la medida normalizada de expresión proporcionada por HPA.
- log2(nTPM+1) reduce asimetría y mantiene definidos los valores cero.
- El Z-score por fila muestra en qué tejidos un gen está relativamente más o menos expresado respecto de su propio perfil; no compara abundancia absoluta entre genes.
- El heatmap de RNA no demuestra abundancia proteica, actividad enzimática ni metabolismo de un fármaco.
- La expresión tisular no sustituye estudios funcionales ni datos farmacocinéticos.

## Archivos producidos
- cyp2c9_network_50x6_ntpm.csv: matriz original nTPM.
- cyp2c9_network_50x6_log2_ntpm.csv: matriz transformada.
- cyp2c9_network_50x6_zscore.csv: matriz estandarizada por gen.
- CYP2C9_expression_heatmap_50x6.png: figura final.
- coverage.json: cobertura y genes sin dato.
- metadata.json: versión, fuente, parámetros y fecha de extracción.

## Ejecución local
    python -m pip install pandas numpy matplotlib requests
    python scripts/02_interactome_string.py --nodes 50 --score 700
    python scripts/03_expression_heatmap.py

## Referencia
The Human Protein Atlas. RNA gene expression data, tissue (HPA version 26.0; Ensembl version 109). https://www.proteinatlas.org/humanproteome/tissue/data
