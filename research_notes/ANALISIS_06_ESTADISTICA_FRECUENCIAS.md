# ANÁLISIS 06. Análisis estadístico descriptivo de frecuencias poblacionales de CYP2C9

## Objetivo
Armonizar la identificación del alelo definitorio de seis variantes CYP2C9 y describir su frecuencia en las poblaciones individuales del panel 1000 Genomes Phase 3.

## Definición de alelos
La asignación se verificó con la tabla de alelos Tier 1 de AMP/CAP y la nomenclatura CYP2C9:
- *2, rs1799853: T (c.430C>T)
- *3, rs1057910: C (c.1075A>C)
- *5, rs28371686: G (c.1080C>G)
- *6, rs9332131: deleción de A (c.818del; en la representación de Ensembl aparece como A frente a AA)
- *8, rs7900194: A (c.449G>A)
- *11, rs28371685: T (c.1003C>T)

La orientación se refiere a la nomenclatura del transcrito CYP2C9 y a la representación alélica de Ensembl; *6 se trata como indel y no como SNV.

## Diseño estadístico
Se realizó análisis descriptivo por variante y población. Para 1000 Genomes Phase 3 se conservaron las 26 poblaciones individuales y se excluyeron de la matriz de comparación los agregados superpoblacionales (AFR, AMR, EAS, EUR, SAS) y ALL para evitar tratar grupos anidados como observaciones independientes. Se calcularon mínimo, máximo, media aritmética no ponderada y mediana entre poblaciones con dato disponible. También se exportó un subconjunto de América (AMR agregado y cohortes MXL, CLM, PEL y PUR), manteniéndolas separadas.

## Resultados y productos
Los archivos generados por el script son:
- cyp2c9_population_frequencies_harmonized.csv: observaciones originales con alelo definitorio armonizado.
- cyp2c9_1000genomes_variant_allele_frequencies.csv: frecuencias del alelo definitorio por población individual.
- cyp2c9_1000genomes_descriptive_summary.csv: resumen descriptivo por variante.
- cyp2c9_1000genomes_americas_focus.csv: AMR y cohortes americanas, sin fusionarlas.
- CYP2C9_1000Genomes_frequency_heatmap.png: mapa de calor de frecuencias.

## Interpretación y límites
La media y mediana entre poblaciones son descripciones no ponderadas de frecuencias poblacionales, no estimaciones combinadas de individuos. No se calcularon intervalos de confianza, pruebas chi-cuadrada ni metaanálisis porque el CSV recuperado no proporciona denominadores/conteos de alelos consistentes para cada observación. No se debe inferir significancia estadística a partir de diferencias descriptivas. Las cohortes MXL (mexicanos de Los Ángeles), CLM (colombianos de Medellín), PEL (peruanos de Lima) y PUR (puertorriqueños) no representan la población nacional de sus países.

## Referencias
Pratt VM et al. (2019). Recommendations for Clinical CYP2C9 Genotyping Allele Selection: A Joint Recommendation of the Association for Molecular Pathology and College of American Pathologists. *The Journal of Molecular Diagnostics*. Tabla de alelos Tier 1 y variantes definitorias. https://pmc.ncbi.nlm.nih.gov/articles/PMC7057225/

Ensembl. REST API, endpoint GET /variation/human/{rsID}?pops=1. https://rest.ensembl.org/documentation/info/variation_id
