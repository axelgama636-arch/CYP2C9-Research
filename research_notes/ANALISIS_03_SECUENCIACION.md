# Análisis 03. Literatura sobre secuenciación y genotipificación de CYP2C9

## Objetivo
Identificar y organizar publicaciones que hayan utilizado secuenciación, genotipificación, análisis de haplotipos o caracterización de variantes de CYP2C9.

## Estrategia de búsqueda
Base de datos: Europe PMC.
Consulta:
\`CYP2C9 AND (sequencing OR "next generation sequencing" OR NGS OR Sanger OR pyrosequencing OR genotyping OR haplotype OR allele)\`

La consulta se ejecuta contra el índice bibliográfico de Europe PMC. Se conservan los metadatos, el resumen, DOI, PMID y PMCID cuando están disponibles. La deduplicación prioriza DOI, después PMID y finalmente título normalizado.

## Estudios prioritarios identificados en la búsqueda bibliográfica inicial

| Estudio | Año | Diseño / método | Identificadores |
|---|---:|---|---|
| Devarajan et al. Pharmacogenomic Next-Generation DNA Sequencing: Lessons from the Identification and Functional Characterization of Variants of Unknown Significance in CYP2C9 and CYP2C19 | 2019 | NGS; 1,013 individuos; identificación y caracterización funcional de variantes | PMID: 30745309; PMCID: PMC6423619; DOI: 10.1124/dmd.118.084269 |
| CYP2C8, CYP2C9, and CYP2C19 Characterization Using Next-Generation Sequencing and Haplotype Analysis: A GeT-RM Collaborative Project | 2022 | NGS y análisis de haplotipos; 137 muestras previamente caracterizadas | PMID: 35134542; PMCID: PMC9069873; DOI: 10.1016/j.jmoldx.2021.12.011 |
| ClinPharmSeq: A targeted sequencing panel for clinical pharmacogenetics implementation | 2022 | Panel de secuenciación dirigida farmacogenética | PMID: 35901010; PMCID: PMC9333201; DOI: 10.1371/journal.pone.0272129 |
| Recommendations for Clinical CYP2C9 Genotyping Allele Selection: A Joint Recommendation of AMP and CAP | 2019 | Recomendación de alelos y diseño de paneles de genotipificación clínica | PMID: 31075510 |
| Current status of the analytical validation of next generation sequencing applications for pharmacogenetic profiling | 2023 | Revisión sistemática de validación analítica NGS farmacogenética | PMID: 37787843 |

## Síntesis metodológica inicial

### Secuenciación NGS
Devarajan et al. analizaron datos de secuenciación de 1,013 individuos e identificaron siete variantes no caracterizadas en CYP2C9 y siete en CYP2C19, con caracterización funcional de las proteínas codificadas. El estudio ilustra que la secuenciación puede detectar variantes raras o poco estudiadas que no se incluyen en paneles dirigidos.

El proyecto GeT-RM reanalizó 137 muestras con secuenciación dirigida y de genoma completo. Los autores describieron alelos/subalelos novedosos y diferencias entre algoritmos de asignación de diplotipos. Esto destaca la importancia de la fase haplotípica, los materiales de referencia y la concordancia entre herramientas.

### Métodos dirigidos frente a secuenciación
La recomendación AMP/CAP propone conjuntos mínimos y ampliados de alelos para pruebas clínicas de CYP2C9. Los paneles dirigidos son útiles para variantes preseleccionadas, pero no identifican variantes no incluidas en el diseño. La secuenciación puede interrogar todas las variantes de la región cubierta, aunque exige validación analítica e interpretación cuidadosa.

### Validación analítica
La revisión de 2023 concluye que la información publicada sobre validación analítica de NGS farmacogenético sigue siendo limitada y que se requiere informar al menos exactitud, precisión, límite de detección y especificidad.

## Extracción de variables para la siguiente fase
Para cada artículo se deberán revisar manualmente:
- población, tamaño muestral y materiales de referencia;
- región genómica y transcrito de referencia;
- plataforma y química de secuenciación;
- profundidad y cobertura;
- pipeline bioinformático y versión;
- variantes y alelos detectados;
- método de asignación de diplotipos/fase;
- validación ortogonal (por ejemplo, Sanger);
- exactitud, concordancia, sensibilidad y especificidad cuando se reporten;
- DOI, PMID, PMCID y disponibilidad de texto completo.

## Archivos generados por el pipeline
- \`data/raw/europe_pmc/cyp2c9_sequencing_search.json\`: respuesta original de Europe PMC.
- \`data/processed/literature/cyp2c9_sequencing_literature.csv\`: corpus deduplicado con metadatos y resúmenes.
- \`data/processed/literature/cyp2c9_sequencing_literature_metadata.json\`: consulta, fecha, conteos y limitaciones.

## Límites
La detección automática de términos en títulos y resúmenes no equivale a extracción manual de métodos. El corpus no es por sí solo una revisión sistemática: se requiere un protocolo de inclusión/exclusión, cribado de texto completo y evaluación de calidad. No se calculará un metaanálisis hasta contar con resultados comparables y denominadores adecuados.

## Referencias principales
1. Devarajan, S. et al. (2019). Pharmacogenomic Next-Generation DNA Sequencing: Lessons from the Identification and Functional Characterization of Variants of Unknown Significance in CYP2C9 and CYP2C19. *Drug Metabolism and Disposition*, 47(4), 425–435. https://doi.org/10.1124/dmd.118.084269
2. (2022). CYP2C8, CYP2C9, and CYP2C19 Characterization Using Next-Generation Sequencing and Haplotype Analysis: A GeT-RM Collaborative Project. *The Journal of Molecular Diagnostics*, 24(4), 337–350. https://doi.org/10.1016/j.jmoldx.2021.12.011
3. Lee, S.-B. et al. (2022). ClinPharmSeq: A targeted sequencing panel for clinical pharmacogenetics implementation. *PLOS ONE*, 17(7), e0272129. https://doi.org/10.1371/journal.pone.0272129
4. Pratt, V.M. et al. (2019). Recommendations for Clinical CYP2C9 Genotyping Allele Selection: A Joint Recommendation of AMP and CAP. https://pubmed.ncbi.nlm.nih.gov/31075510/
5. (2023). Current status of the analytical validation of next generation sequencing applications for pharmacogenetic profiling. https://pubmed.ncbi.nlm.nih.gov/37787843/

## Referencias de búsqueda
- Europe PMC REST API: https://www.ebi.ac.uk/europepmc/webservices/rest/
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=CYP2C9+sequencing
