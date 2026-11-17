# ANÁLISIS 05. Frecuencias alélicas por población para variantes de CYP2C9

## Objetivo
Recuperar y organizar observaciones de frecuencia alélica para seis variantes asociadas con alelos estrella de CYP2C9 (*2, *3, *5, *6, *8 y *11), conservando la población y la fuente reportadas.

## Fuente y método
Se consulta el endpoint de variación de Ensembl REST para cada rsID con el parámetro `pops=1`, que solicita datos de frecuencia por población. La documentación oficial indica que este endpoint devuelve información de variación y, opcionalmente, frecuencias poblacionales. Las observaciones se exportan sin agregarlas entre poblaciones o bases de datos.

## Variables exportadas
- Alelo estrella anotado en el proyecto y rsID.
- Identificador de población tal como lo devuelve Ensembl.
- Alelo y frecuencia (proporción y porcentaje).
- Conteo, cuando la fuente lo proporciona.
- Fuente/dataset inferida del prefijo del identificador de población.
- Indicador de coincidencia con el alelo variante conocido, solo para *2 (T) y *3 (C).

## Criterios de interpretación
1. Una frecuencia ausente no equivale a frecuencia cero.
2. No se promedian frecuencias de poblaciones o datasets distintos.
3. Las etiquetas de población deben conservarse literalmente; no se extrapolan a México.
4. La frecuencia global de un rsID no sustituye la frecuencia del alelo funcional específico si la orientación del alelo no está confirmada.
5. La asignación rsID–alelo estrella es una anotación de trabajo y debe verificarse con fuentes de nomenclatura farmacogenética antes de una publicación.

## Archivos de salida
- `data/processed/population_frequencies/cyp2c9_population_allele_frequencies.csv`
- `data/processed/population_frequencies/cyp2c9_population_frequency_metadata.json`
- Respuestas originales en `data/raw/ensembl/`

## Limitaciones
Ensembl agrega observaciones provenientes de diferentes conjuntos de datos. La disponibilidad, el tamaño muestral y la cobertura pueden variar entre variantes y poblaciones. Estos datos no constituyen una estimación representativa de la población mexicana ni un metaanálisis. Se requiere revisar denominadores y armonizar alelos antes de cualquier comparación cuantitativa.

## Referencia de la API
Ensembl. REST API, endpoint `GET /variation/:species/:id`, parámetro opcional `pops=1`. Documentación: https://rest.ensembl.org/documentation/info/variation_id
