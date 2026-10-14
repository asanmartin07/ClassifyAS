# ClassifyAS — anotación estática de sitios de splicing

Clasifica cada sitio 5′ss/3′ss por **fuerza**, **uso** (constitutivo / alternativo / críptico)
y **tipo de evento** que resulta del no-uso (IR / MXE / alt5ss / alt3ss / ES), en UNA
condición. Corre por encima de la salida de SSE_local.

## Method (as used in the associated manuscript)

For each donor and acceptor site, we classified the dominant type of alternative
splicing event using ClassifyAS, a custom static annotation tool that assigns event
type as a property of the site's local junction topology within a single condition.
For each site with sufficient total coverage (≥20 supporting reads), usage was first
summarized as Splice Site Efficiency (SSE = α/(α+β1+β2), α = canonical junction
reads, β1 = intron-retention reads, β2 = reads supporting a competing, alternative
junction) and classified as constitutive (SSE ≥ 0.95), alternative (0.05 ≤ SSE < 0.95)
or cryptic (SSE < 0.05). For alternative and cryptic sites, the non-canonical fraction
of reads (1 − SSE) was further partitioned among its underlying causes: β1 was
attributed to intron retention (IR), while each individual β2 competing junction
(retained if it reached ≥5 supporting reads and ≥5% of total site coverage) was
assigned to mutually exclusive exon usage (MXE), alternative 5' or 3' splice-site
usage (alt5'ss/alt3'ss), or exon skipping (ES) based on the relative genomic position
of its competing splice site. A site's dominant event type was reported together with
a purity score (fraction of non-canonical reads explained by that single event type);
sites with purity below 0.70 were flagged as mixed events.

## Citation

If you use ClassifyAS, please cite: [San Martín et al., in preparation — Plant
Physiology]. This README will be updated with the final citation upon publication.

## Archivos a pasar

- `classifyas_portable.py` — el script (versión con bloque de config editable).
- este README.

## Requisitos en la máquina destino

1. **SSE_local ya corriendo** (lo tienen). ClassifyAS importa `derive_SSE_from_bams.py`
   de ahí — no lo reimplementa.
2. **Python con `numpy`, `pandas` y `openpyxl`** (el mismo env que usa SSE_local sirve;
   `openpyxl` solo si van a leer/escribir xlsx — el script en sí escribe CSV).
3. Las **entradas que produce SSE_local**: la carpeta `junctions/` con los
   `<muestra>.junctions.bed`, el `targets.csv` (con columna `condition`), y la tabla
   `*.SSE.tsv` (con columnas `<muestra>_alpha` y `<muestra>_beta`).

## Configuración (editar SOLO el bloque CONFIG arriba del script)

Todo lo específico de máquina está en un único bloque marcado
`CONFIGURACIÓN — EDITAR SOLO ESTE BLOQUE`. Hay que setear:

| variable | qué es |
|---|---|
| `SSE_LOCAL_SCRIPTS` | carpeta donde está `derive_SSE_from_bams.py` |
| `CFG_JUNCTIONS` | carpeta con los `.junctions.bed` |
| `CFG_TARGETS` | `targets.csv` (samplesheet con columna `condition`) |
| `CFG_SSE_TSV` | la tabla `*.SSE.tsv` de SSE_local |
| `CFG_STRENGTH_DONOR` / `CFG_STRENGTH_ACCEPTOR` | tablas de fuerza (ver abajo) |
| `CFG_ALL_SAMPLES` | **lista completa** de nombres de muestra, idénticos a los del `.SSE.tsv` y los `.bed` |
| `CFG_OUTDIR` | carpeta de salida |

## Las tablas de fuerza (lo más importante para que sirva en otro dataset)

`CFG_STRENGTH_DONOR` y `CFG_STRENGTH_ACCEPTOR` son **específicas del organismo y de la
anotación** (las provistas son Arabidopsis TAIR10). Están indexadas por `chrom / posición /
strand`, así que:

- **Mismo organismo + misma anotación (TAIR10):** pasarles los dos xlsx tal cual
  (`arabidopsis_5ss_strength_vs_usage_22C_Beckel.xlsx`,
  `arabidopsis_3ss_strength_vs_usage_22C.xlsx`).
- **Otro organismo / otra anotación:** hay que regenerar esas tablas para su genoma. El
  script **corre igual sin ellas** (detecta que faltan y deja las columnas de fuerza
  vacías); uso y tipo de evento se calculan igual.

Formato esperado de las columnas: donores → `chrom, donor_pos_+1, strand, fuerza_E,
percentil_fuerza, is_canonical_GT`; aceptores → `chrom, accpos, strand, fuerza_E_acc,
percentil_fuerza, is_AG`.

## Uso

```bash
python classifyas_portable.py --samples m1,m2,m3 --etiqueta MiCond
python classifyas_portable.py --gene AT1G77080 --samples m1,m2,m3   # un solo locus (validación)
```

- `--samples`: las réplicas de la condición a anotar (subconjunto de `CFG_ALL_SAMPLES`).
- `--etiqueta`: nombre para el archivo de salida.
- `--gene`: opcional; sin él corre sobre todo el genoma.

## Puntos que suelen romper (chequear primero)

1. **Nombres de muestra**: `CFG_ALL_SAMPLES`, las columnas `_alpha`/`_beta` del `.SSE.tsv`
   y los nombres de los `.bed` tienen que coincidir EXACTAMENTE.
2. **Coordenadas / strand**: la regla de MXE y el manejo de hebra fueron validados; el
   sitio se ubica derivándolo de sus propias junturas, no del tipo — no tocar `lado_sitio`.
3. **La columna de fuerza vacía** no es un error: es que faltan las tablas de fuerza para
   ese genoma (ver arriba).

## Cómo leer la salida (dos advertencias con las que ya nos tropezamos)

- **Nunca leer `evento_dominante` solo.** Mirar `pureza` y la columna `eventos`: ~19 % de
  los sitios no constitutivos son mezclas (p. ej. `MXE 52%; IR 47%`).
- **`b1_es_IR = NO`** marca sitios donde el β₁ puede ser cobertura exónica de otra
  isoforma y no retención real; el % de IR reportado es un techo, no una estimación.
