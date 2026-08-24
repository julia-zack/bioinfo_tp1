# utils.py -- helper functions for exercise 4a: amino acid statistics, residue
# sampling, convergence, and loading of real sequences. The "story" (main +
# plots) lives in main.py.

import os
import json
import time
import random
from collections import Counter

import numpy as np
from Bio import Entrez, SeqIO
from Bio.Data.IUPACData import protein_letters

Entrez.email = "julia.zack@rootstocklabs.com"

AA = list(protein_letters)   # the 20 standard amino acids


# ---------------------------------------------------------------------------
# Amino acid statistics
# ---------------------------------------------------------------------------

def natural_statistics(seqs):
    """Return (aa_frequencies, list_of_lengths) for a list of sequences."""
    lengths = [len(s) for s in seqs]
    counts = Counter()
    for s in seqs:
        counts.update(s)
    total = sum(counts.values())
    frequencies = {aa: counts[aa] / total for aa in AA}
    return frequencies, lengths


# Total variation distance between two aa distributions (0 to 1)
def distributions_distance(dist1, dist2):
    diffs_sum = 0.0
    for aa in AA:
        freq_dist1 = dist1.get(aa, 0)   # 0 if that aa is not in dist1
        freq_dist2 = dist2.get(aa, 0)   # 0 if that aa is not in dist2
        diff = abs(freq_dist1 - freq_dist2)
        diffs_sum += diff
    return 0.5 * diffs_sum


#  Amino acid frequencies of a sequence (str) or a list of sequences
def freqs(seq_or_list):
    if isinstance(seq_or_list, str):
        return freqs_seq(seq_or_list)
    return freqs_list(list(seq_or_list))

def freqs_seq(seq):
    frequencies, length = natural_statistics([seq])
    return frequencies

def freqs_list(seqs_list):
    frequencies, length = natural_statistics(seqs_list)
    return frequencies


# ---------------------------------------------------------------------------
# Residue generation / sampling
# ---------------------------------------------------------------------------

def random_residues(n):
    """n residues of a uniformly random protein (each aa with probability 1/20)."""
    random_aa = random.choices(AA, k=n)
    return "".join(random_aa)


def concatenate_til_n(seqs, n):
    """Concatenate proteins until reaching ~n residues (truncated to n)."""
    res = ""
    for seq in seqs:
        res += seq
        if len(res) >= n:
            break
    return res[:n]


NUM_SAMPLE_SIZES = 10   # how many sample sizes to try

def sample_sizes(res_length):
    """~NUM_SAMPLE_SIZES sample sizes spaced on a log scale, from 50 up to
    res_length. DYNAMIC: more data -> larger sizes. The last value is
    res_length, so the curve ends up using ALL the residues."""
    minimum = 50
    points = np.geomspace(minimum, res_length, NUM_SAMPLE_SIZES)  # log-spaced
    sizes = []
    for p in points:
        size = int(round(p))
        if size not in sizes:   # skip duplicates introduced by rounding
            sizes.append(size)
    return sizes


# ---------------------------------------------------------------------------
# Convergence of the distribution
# ---------------------------------------------------------------------------

REPS = 5    # repetitions, to smooth out the randomness

def smooth_curve(res_list, grid, final_distribution):
    """Repeat REPS times (shuffling the order) and average the curves point by
    point, to get a smooth curve that does not depend on any particular
    ordering."""
    curves_per_rep = []
    for rep in range(REPS):
        random.shuffle(res_list)
        shuffle_seq = "".join(res_list)

        curve = []
        for R in grid:
            prefix = shuffle_seq[:R]              # first R residues
            prefix_distribution = freqs(prefix)
            distance = distributions_distance(prefix_distribution, final_distribution)
            curve.append(distance)
        curves_per_rep.append(curve)

    average_curve = list(np.mean(curves_per_rep, axis=0))
    return average_curve


def convergence_curve(res):
    """Grid of sample sizes and its curve (distance to the final distribution)
    for one residue source."""
    res_list = list(res)              # list of residues, so it can be shuffled
    final_distribution = freqs(res)   # distribution using ALL the residues
    res_length = len(res_list)
    grid = sample_sizes(res_length)   # dynamic sample sizes (up to res_length)

    curve = smooth_curve(res_list, grid, final_distribution)
    return grid, curve


THRESHOLD = 0.01   # "stabilised": distance below this value

def find_N(grid, ds):
    """First size in the grid whose distance falls below the threshold."""
    for R, d in zip(grid, ds):
        if d < THRESHOLD:
            return R
    return grid[-1]   # did not stabilise within the range tested


# ---------------------------------------------------------------------------
# Downloading / loading real sequences
# ---------------------------------------------------------------------------

def fetch_proteins(organism, n=300, reviewed_only=True):
    """Up to `n` proteins (SeqRecord) for the organism, via esearch + efetch."""
    term = f'"{organism}"[Organism]'
    if reviewed_only:
        term += " AND refseq[filter]"   # RefSeq: curated sequences

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
        time.sleep(0.4)   # respects the NCBI limit (3 req/s without an API key)
    return records


def clean_sequences(records):
    """Keep only valid aa strings: drop non-standard characters (X, B, Z...)
    and strip the trailing '*' stop."""
    seqs = []
    for r in records:
        s = str(r.seq).upper().rstrip("*")
        if s and set(s) <= set(AA):
            seqs.append(s)
    return seqs


def get_real_sequences(organism, cache):
    """Real protein sequences. Uses the cache file if it exists; otherwise
    downloads them from NCBI and saves them so they are not fetched again."""
    if os.path.exists(cache):
        print(f"Usando secuencias cacheadas de {cache}")
        with open(cache) as f:
            return json.load(f)
    print(f"Descargando proteinas de {organism}...")
    seqs = clean_sequences(fetch_proteins(organism, n=300))
    os.makedirs(os.path.dirname(cache), exist_ok=True)   # create data/ if missing
    with open(cache, "w") as f:
        json.dump(seqs, f)
    return seqs
