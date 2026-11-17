# Análisis 01. Mapa interactómico de CYP2C9

**Organismo:** Homo sapiens (taxón NCBI 9606)  
**Gen central:** CYP2C9 (NCBI Gene 1559; Ensembl ENSG00000138109)  
**Objetivo:** construir una red reproducible de hasta 50 nodos, conservando el tipo y procedencia de cada evidencia.

## Alcance y distinciones de evidencia

- **STRING:** asociaciones funcionales integradas. Una arista no demuestra necesariamente unión física directa.
- **IntAct:** interacciones moleculares curadas; se deben conservar método experimental, identificadores y publicación de respaldo.
- **miRTarBase:** relaciones miRNA–ARNm validadas; constituyen una capa regulatoria distinta de las interacciones proteína-proteína.
- **HIPPIE/HPA:** recursos integradores; sus puntuaciones y criterios no son intercambiables con STRING o IntAct.

No deben combinarse las puntuaciones de STRING, MIscore e HIPPIE en una misma escala.

## Socios candidatos para revisión

Una consulta secundaria de Sugi Atlas resume las siguientes asociaciones STRING para CYP2C9. Se incluyen como **semillas de revisión**, no como sustituto de una descarga y verificación directa de STRING:

| Socio | Score STRING reportado |
|---|---:|
| VKORC1 | 0.982 |
| PPIG | 0.961 |
| SLCO1B1 | 0.907 |
| GGCX | 0.899 |
| UGT1A6 | 0.877 |
| NR1I2 | 0.823 |
| UGT1A4 | 0.794 |
| UGT1A10 | 0.774 |
| UGT1A8 | 0.774 |
| UGT1A1 | 0.773 |
| UGT1A7 | 0.773 |
| UGT2B7 | 0.770 |
| SLC35A2 | 0.751 |
| UGT1A9 | 0.722 |
| ABCB1 | 0.719 |
| TPMT | 0.719 |

La misma consulta secundaria informa seis socios IntAct (CDC37, MAST1, RAB3D, RNF32, STAMBP y STAT5A) con MIscore 0.370 y un registro cada uno. Los registros primarios y sus métodos deben verificarse antes de interpretar estas aristas.

HIPPIE resume POR (0.88), CYP2C18 (0.83), TOMM40 (0.82), ADGRE5 (0.82), CYP2C19 (0.82), HNRNPA1L2 (0.72) y SSR3 (0.72). Estas puntuaciones son propias de HIPPIE.

## Interpretación inicial

POR es un socio biológicamente plausible por su función en la transferencia de electrones a enzimas CYP. CYP2C18 y CYP2C19 aparecen en más de un recurso. Los socios UGT y transportadores incluidos por STRING deben interpretarse como asociaciones funcionales, no como demostración automática de contacto físico.

## Procedimiento reproducible

El script `scripts/02_interactome_string.py` consulta la API de STRING para humano, solicita CYP2C9 más hasta 49 vecinos y usa un umbral inicial de score combinado de 0.700. El número final depende de la respuesta real de la API; no se completará artificialmente hasta 50.

Instalación y ejecución:

```bash
python -m pip install requests
python scripts/02_interactome_string.py --nodes 50 --score 700
```

Genera nodos y aristas en CSV/JSON y metadatos con parámetros, fecha y tamaño real de la red.

## Validaciones pendientes para cerrar el módulo

1. Ejecutar la consulta viva de STRING y archivar la respuesta.
2. Verificar que cada extremo de arista exista en la tabla de nodos y que no haya duplicados.
3. Recuperar los registros primarios de IntAct y sus métodos/PMID.
4. Consultar miRTarBase para CYP2C9 como gen diana y conservar el tipo de ensayo.
5. Construir la figura final de 50 nodos únicamente después de la recuperación real.

## Fuentes

- STRING CYP2C9: https://string-db.org/network/9606.ENSP00000260682
- STRING API: https://string-db.org/help/api/
- STRING score interpretation: https://string-db.org/help/scores/
- Human Protein Atlas, CYP2C9 interactions: https://www.proteinatlas.org/ENSG00000138109-CYP2C9/interaction
- Sugi Atlas, CYP2C9: https://www.sugi.bio/atlas/gene/CYP2C9/
- HIPPIE: https://cbdm-01.zdv.uni-mainz.de/~mschaefer/tipi/public/query.php?s=CP2C9_HUMAN
- Cui et al. (2025). miRTarBase 2025. Nucleic Acids Research, 53(D1), D147–D156. https://doi.org/10.1093/nar/gkae1072
