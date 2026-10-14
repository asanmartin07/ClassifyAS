# -*- coding: utf-8 -*-
"""
ClassifyAS v1.0 - anotacion ESTATICA de sitios de splicing.

Describe la arquitectura de cada sitio 5'ss/3'ss en UNA condicion, sin comparar
condiciones. Es la pieza estatica; el analisis diferencial (respuesta a genotipo o a
ambiente) va aparte, y la combinacion de ambos es un tercer paso.

Motivacion de la separacion: el tipo de evento es una propiedad de la TOPOLOGIA del
locus, no del contraste. Derivarlo de una comparacion hacia que el mismo sitio saliera
MXE en un contraste e IR en otro (caso FLM/AT1G77080), lo cual es incoherente.

Por cada sitio informa:
  - posicion, tipo (donor/acceptor), dinucleotido canonico
  - FUERZA: dG de hibridacion (ViennaRNA), cruda y en percentil DENTRO DE SU CLASE
      donores  : 9-mer (-3..+6) contra U1
      aceptores: 23-mer (tracto polipirimidinico + AG)
    Las escalas de donor y aceptor NO son comparables entre si: comparar por percentil.
  - alpha, beta1, beta2 y SSE = alpha/(alpha+beta1+beta2)
  - uso: constitutivo (SSE>=0.95) / alternativo (0.05<=SSE<0.95) / criptico (SSE<0.05)
  - sitios alternativos: posicion, ancla compartida y lecturas de cada competidor de beta2
  - reparto del NO-USO (1-SSE) en tipos de evento: IR desde beta1; MXE / alt5ss /
    alt3ss / ES(n) desde cada competidor de beta2, clasificados por geometria

LIMITACIONES CONOCIDAS (no borrar; son reales y estan medidas):
  1. No calcula fuerza de sitios NO ANOTADOS. Las tablas de fuerza derivan de TAIR10, y
     los sitios cripticos suelen no estar en ellas -> quedan sin dato de fuerza.
  2. beta1 no siempre es retencion. Si una isoforma alternativa deja al sitio EXONICO,
     su cobertura se cuenta como beta1 y simula IR. Se senaliza de dos formas
     (columnas 'exonizantes' y 'asim_b1' -> 'b1_es_IR'), pero AMBAS tienen un punto
     ciego: cuando el exon alternativo cubre los DOS extremos del intron, beta1 sube
     simetricamente y ninguna de las dos lo detecta (caso AT1G09140:2943430). El % de
     IR reportado es por lo tanto un TECHO, no una estimacion.
  3. La regla MXE se valido sobre FLM (hebra +). En el censo genomico detecta en ambas
     hebras (513 + / 358 -) pero con un sesgo ~5 sigma hacia hebra + que no esta
     explicado y no aparece en las otras categorias.
  4. Depende de SSE-local aguas arriba (alpha/beta2 y la tabla de beta1). No es autonomo.

DECISIONES DE UMBRAL (todas distribucionales, no biofisicas):
  - fuerza: percentil <25 debil, >75 fuerte. Elegido por continuidad con los analisis
    previos de PRMT5 del proyecto.
  - constitutivo 0.95 / criptico 0.05.
  - evidencia: T>=20; via presente si >=5 lecturas Y >=5% de T.

Se probaron medidas de fuerza secundarias y se descartaron (ver spec):
  donores  : recalibracion de Beckel 2023 -> cambia 11.6% de las clases, NUNCA invierte
             debil<->fuerte, y predice el uso identico (rho +0.132 vs +0.130). No aporta.
  aceptores: % de pirimidinas del tracto -> discrepa mucho (51% cambia de clase) pero es
             PEOR predictor del uso (+0.080 vs +0.234 en sitios alternativos). No es una
             medida alternativa legitima, es una peor.
"""
import sys, argparse, bisect, collections, numpy as np, pandas as pd

# ============================================================================
#  CONFIGURACION — EDITAR SOLO ESTE BLOQUE PARA CORRER EN OTRA MAQUINA
# ============================================================================
# (1) Carpeta con los scripts de SSE_local (donde esta derive_SSE_from_bams.py):
SSE_LOCAL_SCRIPTS = "/ruta/a/SSE_local/scripts"

# (2) Entradas producidas por tu corrida de SSE_local:
CFG_JUNCTIONS = "/ruta/a/sse_run/junctions"      # carpeta con <muestra>.junctions.bed
CFG_TARGETS   = "/ruta/a/sse_run/targets.csv"    # samplesheet con columna 'condition'
CFG_SSE_TSV   = "/ruta/a/sse_run/salida.SSE.tsv" # tabla SSE (columnas <muestra>_alpha, <muestra>_beta)

# (3) Tablas de FUERZA del sitio (especificas del organismo/anotacion; para Arabidopsis
#     TAIR10 son las provistas). Si no las tenes, el script corre igual pero deja la
#     columna de fuerza vacia. Formato: donor -> cols chrom, donor_pos_+1, strand,
#     fuerza_E, percentil_fuerza, is_canonical_GT ; acceptor -> chrom, accpos, strand,
#     fuerza_E_acc, percentil_fuerza, is_AG.
CFG_STRENGTH_DONOR    = "/ruta/a/arabidopsis_5ss_strength_vs_usage_22C_Beckel.xlsx"
CFG_STRENGTH_ACCEPTOR = "/ruta/a/arabidopsis_3ss_strength_vs_usage_22C.xlsx"

# (4) TODAS las muestras del experimento (mismos nombres que en el SSE.tsv y en los .bed).
#     Definen el filtro de soporte por replica. Editar la lista COMPLETA:
CFG_ALL_SAMPLES = ["muestra1", "muestra2", "muestra3", "..."]

# (5) Carpeta de salida:
CFG_OUTDIR = "/ruta/a/salida"
# ============================================================================

sys.path.insert(0, SSE_LOCAL_SCRIPTS)
import derive_SSE_from_bams as D

# rutas y muestras: TODAS vienen del bloque CONFIG de arriba
JUNCT = CFG_JUNCTIONS; TARGETS = CFG_TARGETS; SSE = CFG_SSE_TSV; RES = CFG_OUTDIR
ALL42 = CFG_ALL_SAMPLES

MIN_T        = 20     # evidencia minima para anotar el sitio
MIN_PRESENT  = 5      # lecturas minimas para considerar presente una via
MIN_FRAC     = 0.05   # ... y ademas fraccion minima de T
CONSTITUTIVO = 0.95
CRIPTICO     = 0.05
P_DEBIL, P_FUERTE = 25, 75
MAX_EXON     = 3000
ASIM_SOSP    = 5      # asimetria de beta1 entre extremos que vuelve sospechosa la IR
PUREZA_MIX   = 0.70   # por debajo de esto el sitio se marca 'mixto'

ap = argparse.ArgumentParser(description="ClassifyAS v1.0")
ap.add_argument("--gene", default=None, help="sin este flag corre sobre TODOS los sitios")
ap.add_argument("--samples", default="Col0_1,Col0_2,Col0_3")
ap.add_argument("--etiqueta", default="Col22")
a = ap.parse_args()
SMP = a.samples.split(",")

# ---------------------------------------------------------------- fuerza (solo dG)
print("1) tablas de fuerza...")
import os
_hay_fuerza = os.path.exists(CFG_STRENGTH_DONOR) and os.path.exists(CFG_STRENGTH_ACCEPTOR)
if not _hay_fuerza:
    print("   [aviso] tablas de fuerza no encontradas -> columnas dG/percentil/clase quedaran vacias")
    FD = {}; FA = {}
else:
  dn = pd.read_excel(CFG_STRENGTH_DONOR)
  dn = dn.rename(columns={"donor_pos_+1": "donor_pos"})
  FD = {(str(r.chrom), int(r.donor_pos), r.strand): (r.fuerza_E, r.percentil_fuerza, r.is_canonical_GT)
        for r in dn.itertuples()}
  ac = pd.read_excel(CFG_STRENGTH_ACCEPTOR)
  FA = {(str(r.chrom), int(r.accpos), r.strand): (r.fuerza_E_acc, r.percentil_fuerza, r.is_AG)
        for r in ac.itertuples()}
print(f"   donores {len(FD)}  aceptores {len(FA)}")

# ---------------------------------------------------------------- junturas
print("2) junturas...")
jc = collections.defaultdict(lambda: collections.defaultdict(int))
for s in ALL42:
    for j in D.parse_junctions_bed(f"{JUNCT}/{s}.junctions.bed", min_reads=1):
        jc[(j["chrom"], j["intron_start"], j["intron_end"], j["strand"])][s] += j["count"]
jg = D.filter_by_replicate_support(dict(jc), D.load_condition_map(TARGETS, "condition"), 2)
by_ep, _ = D.build_junction_indices(jg)
jm = {k: float(np.mean([v.get(s, 0) for s in SMP])) for k, v in jg.items()}

by_start = collections.defaultdict(set); by_end = collections.defaultdict(set); jset = set()
srt_S = collections.defaultdict(list); srt_E = collections.defaultdict(list)
for (ch, S, E, sd) in jg:
    by_start[(ch, sd, S)].add(E); by_end[(ch, sd, E)].add(S); jset.add((ch, S, E, sd))
    srt_S[(ch, sd)].append(S); srt_E[(ch, sd)].append(E)
for d in (srt_S, srt_E):
    for k in d: d[k] = sorted(set(d[k]))


def _v(l, lo, hi): return l[bisect.bisect_left(l, lo):bisect.bisect_right(l, hi)]


def lado_sitio(ajs, p0):
    """Lado GENOMICO del sitio dentro de su intron, deducido de los datos y no del tipo:
    'izq' si p0 == intron_start, 'der' si p0 == intron_end - 1.
    NO deducirlo de SiteType: en hebra + el donor es el extremo izquierdo y en hebra - es
    el derecho. Deducirlo del tipo devolvia None en hebra - y deshabilitaba EN SILENCIO
    la regla MXE en la mitad del genoma.
    Display, tambien independiente de hebra: izq -> Site = S+1, der -> Site = E."""
    if all(S == p0 for (c, S, E, d) in ajs): return "izq", p0
    if all(E == p0 + 1 for (c, S, E, d) in ajs): return "der", p0 + 1
    return None, None


def mxe_partner(ch, sd, pa, lado, jc_):
    """El competidor define un exon mutuamente excluyente con el exon que contiene al sitio:
    anclas izquierda y derecha compartidas, exones disjuntos, y NINGUNA juntura que los una."""
    _, s2, e2, _ = jc_
    if lado == "izq":
        if e2 not in by_start[(ch, sd, pa)]: return False
        lim = min(pa, s2)
        for S in _v(srt_S[(ch, sd)], lim - 2 * MAX_EXON, lim):
            ends = by_start[(ch, sd, S)]
            for x1 in ends:
                if not (S < x1 < pa and pa - x1 <= MAX_EXON): continue
                for y1 in ends:
                    if y1 == x1 or not (S < y1 < s2 and s2 - y1 <= MAX_EXON): continue
                    if not (pa <= y1 or s2 <= x1): continue
                    if (ch, pa, y1, sd) in jset or (ch, s2, x1, sd) in jset: continue
                    return True
        return False
    if s2 not in by_end[(ch, sd, pa)]: return False
    lim = max(pa, e2)
    for E in _v(srt_E[(ch, sd)], lim, lim + 2 * MAX_EXON):
        starts = by_end[(ch, sd, E)]
        for x2 in starts:
            if not (pa < x2 < E and x2 - pa <= MAX_EXON): continue
            for y2 in starts:
                if y2 == x2 or not (e2 < y2 < E and y2 - e2 <= MAX_EXON): continue
                if not (x2 <= e2 or y2 <= pa): continue
                if (ch, x2, e2, sd) in jset or (ch, y2, pa, sd) in jset: continue
                return True
    return False


def n_exones(ch, sd, S, E):
    """Exones completos dentro del intron [S,E). ES exige >=1: sin eso, 'S < p0 < E' se
    cumple en cualquier sitio alternativo corriente y se etiquetaban ES los alt-3'ss rio
    abajo / alt-5'ss rio arriba, de forma asimetrica."""
    ends_in = _v(srt_E[(ch, sd)], S + 1, E - 1); starts_in = _v(srt_S[(ch, sd)], S + 1, E - 1)
    n = 0; ult = S
    for x1 in ends_in:
        if x1 < ult: continue
        x2 = next((s for s in starts_in if x1 < s <= E and s - x1 <= MAX_EXON), None)
        if x2 is None: continue
        n += 1; ult = x2
    return n


print("3) alpha/beta2...")
sites, idx = D.build_sites_from_junctions(jg)
jg_sub = {jk: {s: c for s, c in v.items() if s in set(SMP)} for jk, v in jg.items()}
alpha, beta2 = D.compute_alpha_beta2(sites, idx, jg_sub, by_ep, None, SMP, "shared_anchor")
alpha = np.array(alpha, float).mean(1); beta2 = np.array(beta2, float).mean(1)
kidx = {f"{s['chrom']}:{s['pos']+1}:{s['strand']}": i for i, s in enumerate(sites)}

print("4) beta1...")
T = pd.read_csv(SSE, sep="\t", usecols=["Region", "Site", "Strand", "Gene", "SiteType"]
                + [f"{s}_beta" for s in SMP], low_memory=False)
T["key"] = T.Region.astype(str) + ":" + T.Site.astype(str) + ":" + T.Strand
T = T.set_index("key"); T = T[~T.index.duplicated(keep="first")]
bt = T[[f"{s}_beta" for s in SMP]].mean(axis=1).to_dict()
keys = (list(T[T.Gene.astype(str).str.contains(a.gene, na=False)].index) if a.gene else list(T.index))
print(f"5) anotando {len(keys)} sitios ({a.gene or 'GENOMA'}) en {a.etiqueta}...")

rows = []
for n_i, k in enumerate(keys):
    if not a.gene and n_i and n_i % 20000 == 0:
        print(f"   ... {n_i}/{len(keys)}", flush=True)
    if k not in kidx: continue
    i = kidx[k]; ch, pos, strand = k.split(":"); pos = int(pos); p0 = pos - 1
    tp = T.at[k, "SiteType"]; gene = T.at[k, "Gene"]
    al = float(alpha[i]); b2 = float(beta2[i]); b1 = max(0.0, float(bt.get(k, 0.0)) - b2)
    tot = al + b1 + b2
    fz = FD.get((ch, pos, strand)) if tp == "donor" else FA.get((ch, pos, strand))
    dG, pct, canon = fz if fz else (np.nan, np.nan, np.nan)
    base = dict(Gene=gene, chrom=ch, Site=pos, Strand=strand, tipo=tp, canonico=canon,
                dG=dG, pct_dG=pct, alpha=round(al, 1), beta1=round(b1, 1), beta2=round(b2, 1))

    if tot < MIN_T:
        rows.append(dict(base, clase_fuerza="", SSE=np.nan, uso="sin evidencia (T<20)",
                         no_uso=np.nan, pareja="", b1_pareja=np.nan, asim_b1=np.nan,
                         exonizantes="", b1_es_IR="-", sitios_alternativos="",
                         eventos="", evento_dominante="sin evidencia",
                         pureza=np.nan, mixto=np.nan))
        continue

    sse = al / tot
    uso = ("constitutivo" if sse >= CONSTITUTIVO else
           "criptico" if sse < CRIPTICO else "alternativo")
    cf = ("" if pct != pct else
          "debil" if pct < P_DEBIL else "fuerte" if pct > P_FUERTE else "intermedia")

    ajs = by_ep.get((ch, p0, strand), set()); lado, pa = lado_sitio(ajs, p0)
    other = set()
    for jk in ajs:
        cc, S, E, sd = jk
        (dp, ds), (apo, as_) = D.site_keys_for_junction(S, E, sd)
        if (dp, ds) == (p0, strand): other.add((cc, apo, as_))
        elif (apo, as_) == (p0, strand): other.add((cc, dp, ds))
    comp = {jk for ep in other for jk in by_ep.get(ep, set()) if jk not in ajs}

    # beta1 en el OTRO extremo del mismo intron: la retencion real es simetrica,
    # el artefacto exonizante de un solo extremo no lo es
    b1p = np.nan; pareja = ""
    if ajs:
        cc, S, E, sd = max(ajs, key=lambda j: jm.get(j, 0.0))
        (dp, ds), (apo, as_) = D.site_keys_for_junction(S, E, sd)
        for pp in (dp + 1, apo + 1):
            if pp != pos:
                pk = f"{ch}:{pp}:{strand}"
                if pk in kidx and pk in bt:
                    b1p = max(0.0, float(bt[pk]) - float(beta2[kidx[pk]])); pareja = str(pp)
    asim = (b1 + 1) / (b1p + 1) if b1p == b1p else np.nan

    # competidores exonizantes: intron contenido en el alfa y que NO cubre el sitio
    exz = []
    for (cc, S, E, sd) in ajs:
        for jk2 in comp:
            c2, s2, e2, d2 = jk2
            if c2 != cc or d2 != sd: continue
            if s2 >= S and e2 <= E and not (s2 <= p0 < e2) and jm[jk2] >= MIN_PRESENT:
                exz.append(f"{(s2+1) if lado=='izq' else e2}({jm[jk2]:.0f}r)")
    exz = sorted(set(exz))

    # ---- reparto del NO-USO ----
    parts = collections.defaultdict(float); alts = []
    if b1 >= MIN_PRESENT and b1 / tot >= MIN_FRAC:
        parts["IR"] += b1 / tot
    for jk in sorted(comp):
        cc, S, E, sd = jk
        r = jm[jk]
        if r < MIN_PRESENT or r / tot < MIN_FRAC: continue
        if lado == "izq": sitio_alt, ancla = S + 1, E     # comparte el ancla derecha
        else:             sitio_alt, ancla = E, S + 1     # ... o la izquierda
        if pa is not None and mxe_partner(cc, sd, pa, lado, jk):
            cat = "MXE"
        elif S < p0 < E and n_exones(cc, sd, S, E) >= 1:
            cat = f"ES(n={n_exones(cc, sd, S, E)})"
        else:
            cat = "alt5ss" if tp == "donor" else "alt3ss"
        parts[cat] += r / tot
        alts.append(f"{sitio_alt}[{cat},{r:.0f}r,ancla {ancla}]")

    nu = 1 - sse
    orden = sorted(parts.items(), key=lambda x: -x[1])
    ev = "; ".join(f"{c} {v/nu:.0%}" for c, v in orden) if (nu > 0.01 and parts) else ""
    # etiqueta explicita en vez de celda vacia: distinto es 'constitutivo' de
    # 'cambia pero nada supero los filtros'
    if uso == "constitutivo":      dom = "constitutivo"
    elif orden:                    dom = orden[0][0].split("(")[0]
    else:                          dom = "sin atribuir"
    # PUREZA: 'evento_dominante' solo es engañoso en las mezclas. En el censo, el 18.7% de
    # los sitios no constitutivos tienen pureza <70%, y hay casos de empate (FLM 28958491
    # sale 'IR' con IR 36% / MXE 36% / alt5ss 27%). Reportar siempre pureza junto al dominante.
    pur = (orden[0][1] / nu) if (orden and nu > 0.01) else np.nan
    mix = (pur < PUREZA_MIX) if pur == pur else np.nan
    b1ok = "-" if b1 < MIN_PRESENT else ("NO" if (exz or (asim == asim and asim >= ASIM_SOSP)) else "si")
    rows.append(dict(base, clase_fuerza=cf, SSE=round(sse, 3), uso=uso, no_uso=round(nu, 3),
                     pareja=pareja, b1_pareja=round(b1p, 1) if b1p == b1p else np.nan,
                     asim_b1=round(asim, 1) if asim == asim else np.nan,
                     exonizantes=" ".join(exz), b1_es_IR=b1ok,
                     sitios_alternativos=" ".join(alts), eventos=ev, evento_dominante=dom,
                     pureza=round(pur, 3) if pur == pur else np.nan, mixto=mix))

R = pd.DataFrame(rows).sort_values(["chrom", "Site"])
out = f"{RES}/classifyas_{a.gene or 'GENOMA'}_{a.etiqueta}.csv"
R.to_csv(out, index=False)

if a.gene:
    cols = ["Site", "tipo", "dG", "pct_dG", "clase_fuerza", "alpha", "beta1", "beta2",
            "SSE", "uso", "b1_es_IR", "evento_dominante", "sitios_alternativos"]
    with pd.option_context("display.width", 300, "display.max_columns", 40,
                           "display.max_rows", 300, "display.max_colwidth", 55):
        print(R[cols].to_string(index=False))
else:
    ev_ = R[~R.uso.isin(["constitutivo", "sin evidencia (T<20)"])]
    print(f"\n=== CENSO {a.etiqueta}: {len(R)} sitios ===")
    print(R.uso.value_counts().to_string())
    print(f"\n-- atribuciones entre los {len(ev_)} no-constitutivos --")
    cnt = collections.Counter()
    for e in ev_.eventos.dropna():
        for tok in str(e).split(";"):
            t = tok.strip().split(" ")[0].split("(")[0]
            if t: cnt[t] += 1
    tt = sum(cnt.values())
    for c, n in cnt.most_common():
        print(f"   {c:<8} {n:>7}  ({100*n/max(1,tt):.1f}%)")
    print(f"\n-- evento dominante --"); print(ev_.evento_dominante.value_counts().to_string())
    print(f"\n-- uso por tipo de sitio --"); print(pd.crosstab(R.tipo, R.uso).to_string())
    j = R[R.b1_es_IR != "-"]
    print(f"\n-- beta1 juzgable (n={len(j)}): 'NO' = {int((j.b1_es_IR=='NO').sum())} "
          f"({100*(j.b1_es_IR=='NO').mean():.1f}%) -> el % de IR es un TECHO --")
print(f"\n-- clase de fuerza --"); print(R[R.clase_fuerza != ""].clase_fuerza.value_counts().to_string())
print(f"\ntabla -> {out}")
