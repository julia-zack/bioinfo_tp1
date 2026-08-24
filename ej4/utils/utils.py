# aux.py -- funciones auxiliares del ejercicio 4a:
# estadistica de aminoacidos, muestreo de residuos, convergencia y carga de
# secuencias reales. La "historia" (main + graficos) vive en 4-a.py.

import os
import json
import time
import random
from collections import Counter

import numpy as np
from Bio import Entrez, SeqIO
from Bio.Data.IUPACData import protein_letters

Entrez.email = "julia.zack@rootstocklabs.com"

AA = list(protein_letters)   # 20 aminoacidos estandar


# ---------------------------------------------------------------------------
# Estadistica de aminoacidos
# ---------------------------------------------------------------------------

def natural_statistics(seqs):
    """Devuelve (frecuencias_aa, lista_de_largos) a partir de una lista de secuencias."""
    lengths = [len(s) for s in seqs]
    counts = Counter()
    for s in seqs:
        counts.update(s)
    total = sum(counts.values())
    frecuencias = {aa: counts[aa] / total for aa in AA}
    return frecuencias, lengths


# Distancia de variacion total entre dos distribuciones de aa (0 a 1)
def distributions_distance(dist1, dist2):
    diffs_sum = 0.0
    for aa in AA:
        freq_dist1 = dist1.get(aa, 0)   # 0 si ese aa no esta en dist1
        freq_dist2 = dist2.get(aa, 0)   # 0 si ese aa no esta en dist2
        diff = abs(freq_dist1 - freq_dist2)
        diffs_sum += diff
    return 0.5 * diffs_sum


#  Frecuencias de aa de una secuencia (str) o lista de secuencias
def freqs(seq_or_list):
    if isinstance(seq_or_list, str):
        return freqs_seq(seq_or_list)
    return freqs_list(list(seq_or_list))

def freqs_seq(seq):
    frecuencias, length = natural_statistics([seq])
    return frecuencias

def freqs_list(seqs_list):
    frecuencias, length = natural_statistics(seqs_list)
    return frecuencias


# ---------------------------------------------------------------------------
# Generacion / muestreo de residuos
# ---------------------------------------------------------------------------

def random_residues(n):
    """n residuos de una proteina uniforme al azar (cada aa con prob 1/20)."""
    residuos_al_azar = random.choices(AA, k=n)
    return "".join(residuos_al_azar)


def concatenate_til_n(seqs, n):
    """Concatena proteinas hasta juntar ~n residuos (recorta a n)."""
    res = ""
    for seq in seqs:
        res += seq
        if len(res) >= n:
            break
    return res[:n]


NUM_SAMPLE_SIZES = 10   # cuantos tamaños de muestra probar

def sample_sizes(res_length):
    """~NUM_SAMPLE_SIZES tamaños de muestra espaciados en escala logaritmica,
    desde 50 hasta res_length. DINAMICO: mas datos -> tamaños mas grandes. El
    ultimo valor es res_length, asi la curva llega a usar TODOS los residuos."""
    minimo = 50
    puntos = np.geomspace(minimo, res_length, NUM_SAMPLE_SIZES)  # log-espaciados
    tamanos = []
    for p in puntos:
        t = int(round(p))
        if t not in tamanos:   # evitamos duplicados que puedan surgir del redondeo
            tamanos.append(t)
    return tamanos


# ---------------------------------------------------------------------------
# Convergencia de la distribucion
# ---------------------------------------------------------------------------

REPS = 5    # repeticiones para suavizar el azar

def smooth_curve(res_list, grid, final_distribution):
    """Repite REPS veces (mezclando el orden) y promedia las curvas punto a
    punto, para obtener una curva suave que no dependa del orden particular."""
    curvas_por_repeticion = []
    for rep in range(REPS):
        random.shuffle(res_list)
        shuffle_seq = "".join(res_list)

        curva = []
        for R in grid:
            prefijo = shuffle_seq[:R]              # primeros R residuos
            distribucion_prefijo = freqs(prefijo)
            distancia = distributions_distance(distribucion_prefijo, final_distribution)
            curva.append(distancia)
        curvas_por_repeticion.append(curva)

    curva_promedio = list(np.mean(curvas_por_repeticion, axis=0))
    return curva_promedio


def curva_convergencia(res):
    """Grilla de tamaños de muestra y su curva (distancia a la distribucion
    final) para una fuente de residuos."""
    res_list = list(res)              # lista de residuos, para poder mezclarla
    final_distribution = freqs(res)   # distribucion usando TODOS los residuos
    res_length = len(res_list)
    grid = sample_sizes(res_length)   # tamaños de muestra dinamicos (hasta res_length)

    curva = smooth_curve(res_list, grid, final_distribution)
    return grid, curva


UMBRAL = 0.01   # "estabilizada": distancia por debajo de este valor

def encontrar_N(grid, ds):
    """Primer tamano de la grilla con distancia por debajo del umbral."""
    for R, d in zip(grid, ds):
        if d < UMBRAL:
            return R
    return grid[-1]   # no se estabilizo dentro del rango probado


# ---------------------------------------------------------------------------
# Descarga / carga de secuencias reales
# ---------------------------------------------------------------------------

def fetch_proteins(organism, n=300, reviewed_only=True):
    """Hasta `n` proteinas (SeqRecord) del organismo, via esearch + efetch."""
    term = f'"{organism}"[Organism]'
    if reviewed_only:
        term += " AND refseq[filter]"   # RefSeq: secuencias curadas

    handle = Entrez.esearch(db="protein", term=term, retmax=n, usehistory="y")
    search = Entrez.read(handle)
    handle.close()

    webenv = search["WebEnv"]
    query_key = search["QueryKey"]
    total = min(int(search["Count"]), n)
    print(f"{organism}: {search['Count']} proteinas encontradas, bajando {total}")

    records = []
    batch = 200
    for start in range(0, total, batch):
        handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text",
                               retstart=start, retmax=batch,
                               webenv=webenv, query_key=query_key)
        records.extend(SeqIO.parse(handle, "fasta"))
        handle.close()
        time.sleep(0.4)   # respeta el limite de NCBI (3 req/s sin API key)
    return records


def clean_sequences(records):
    """Solo cadenas de aa validas: descarta caracteres no estandar (X, B, Z...)
    y saca el '*' de stop del final."""
    seqs = []
    for r in records:
        s = str(r.seq).upper().rstrip("*")
        if s and set(s) <= set(AA):
            seqs.append(s)
    return seqs


def get_real_sequences(organism, cache):
    """Secuencias de proteinas reales. Usa el archivo cache si existe; si no,
    las descarga de NCBI y las guarda para no volver a bajarlas."""
    if os.path.exists(cache):
        print(f"Usando secuencias cacheadas de {cache}")
        with open(cache) as f:
            return json.load(f)
    print(f"Descargando proteinas de {organism}...")
    seqs = clean_sequences(fetch_proteins(organism, n=300))
    os.makedirs(os.path.dirname(cache), exist_ok=True)   # crea data/ si no existe
    with open(cache, "w") as f:
        json.dump(seqs, f)
    return seqs
