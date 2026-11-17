# Análisis 04. Catálogo de variantes CYP2C9 en ClinVar

## Objetivo
Construir un inventario reproducible de registros de variantes asociados con CYP2C9 en NCBI ClinVar, conservando la anotación clínica y el estado de revisión tal como aparecen en la fuente.

## Variantes de referencia prioritarias
- CYP2C9*2: rs1799853; NM_000771.4:c.430C>T.
- CYP2C9*3: rs1057910; NM_000771.4:c.1075A>C.

Estas dos correspondencias se usan como referencias conocidas. No se asignan alelos estrella a otras variantes por inferencia automática.

## Fuente y estrategia
- Base principal: NCBI ClinVar.
- Consulta E-utilities: `CYP2C9[gene]`.
- ESearch recupera los UID de ClinVar; ESummary recupera resúmenes por lotes.
- Se conserva el JSON original y se genera una tabla CSV normalizada.
- Se aplican reintentos con espera progresiva ante errores HTTP temporales.

## Variables extraídas
UID, accession, título, tipo de objeto, símbolo génico, rsID detectados, mapeo conocido *2/*3, significancia clínica reportada, estado de revisión, fecha de evaluación y resumen de variación.

## Interpretación
La significancia clínica de ClinVar se presenta como dato reportado por los remitentes/revisores de ClinVar; no se reclasifica. El estado de revisión (review status) es esencial para contextualizar la solidez de la afirmación. Los registros de ClinVar no equivalen a estudios independientes ni a frecuencias poblacionales.

## Archivos generados
- `data/raw/clinvar/cyp2c9_clinvar_esummary.json`
- `data/processed/variants/cyp2c9_clinvar_variants.csv`
- `data/processed/variants/cyp2c9_clinvar_metadata.json`

## Limitaciones
1. La búsqueda depende de la indexación de ClinVar por gen y puede recuperar registros con distintos alcances moleculares.
2. La detección de rsID se realiza sobre el resumen serializado; los registros complejos deben revisarse en el JSON original.
3. La normalización HGVS completa y la asignación de haplotipos requieren revisión adicional de transcrito y PharmVar.
4. No se calculan frecuencias alélicas: este módulo no obtiene denominadores genotípicos poblacionales.
5. Este catálogo es descriptivo, no constituye una interpretación clínica.

## Próximo paso
Contrastar las variantes de ClinVar con PharmVar y dbSNP, revisar manualmente las anotaciones de *2 y *3 y, solo después, preparar un análisis de frecuencias poblacionales con datos y denominadores verificables.
