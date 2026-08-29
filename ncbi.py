"""Everything that talks to NCBI: protein and nucleotide fetching, plus the
local cache so the same records are not downloaded twice."""

import os
import json
import time
import random

from Bio import Entrez, SeqIO

from sequences import AA

# NCBI requires a contact address on every request.
Entrez.email = 'arielzingman@gmail.com'

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')


def cache_path(filename):
    return os.path.join(DATA_DIR, filename)


# ---------------------------------------------------------------------------
# Proteins (used by exercise 4)
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
    print(f"{organism}: {search['Count']} proteins found, downloading {total}")

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
        print(f"Using cached sequences from {cache}")
        with open(cache) as f:
            return json.load(f)
    print(f"Downloading proteins for {organism}...")
    seqs = clean_sequences(fetch_proteins(organism, n=300))
    os.makedirs(os.path.dirname(cache), exist_ok=True)   # create data/ if missing
    with open(cache, "w") as f:
        json.dump(seqs, f)
    return seqs


# ---------------------------------------------------------------------------
# Nucleotide sequences (used by exercise 3)
# ---------------------------------------------------------------------------

NUCLEOTIDE_CACHE = cache_path("ncbi_sample.fasta")

# A random sample of this query, not of all GenBank: the query is part of the method.
NUCLEOTIDE_QUERY = '"Homo sapiens"[Organism] AND biomol_mrna[PROP] AND 500:5000[SLEN]'


def fetch_random_records(n=100, query=NUCLEOTIDE_QUERY,
                         cache=NUCLEOTIDE_CACHE, seed=None):
    """Download n sequences sampled at random from the results of a query.

    Writes a multi-FASTA to `cache` and reuses it on later runs, so the
    requests to NCBI are not repeated.
    """
    if os.path.exists(cache):
        return list(SeqIO.parse(cache, "fasta"))

    # esearch returns at most 10000 ids per request; that is the pool the
    # sample is drawn from.
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=10000)
    id_pool = Entrez.read(handle)["IdList"]
    handle.close()
    if len(id_pool) < n:
        raise ValueError(f"The query returned {len(id_pool)} ids, fewer than the {n} requested.")

    sampled_ids = random.Random(seed).sample(id_pool, n)

    # In batches, respecting the 3 requests per second limit without an API key.
    records = []
    batch_size = 20
    for start in range(0, len(sampled_ids), batch_size):
        batch = sampled_ids[start:start + batch_size]
        handle = Entrez.efetch(db="nucleotide", id=",".join(batch),
                               rettype="fasta", retmode="text")
        records.extend(SeqIO.parse(handle, "fasta"))
        handle.close()
        print(f"  downloaded {len(records)}/{len(sampled_ids)}")
        time.sleep(0.4)

    os.makedirs(os.path.dirname(cache), exist_ok=True)
    SeqIO.write(records, cache, "fasta")
    return records
