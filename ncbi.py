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

def fetch_proteins(organism, n=None, reviewed_only=True, seed=0):
    """Every protein (SeqRecord) the organism has in Swiss-Prot.

    Pass `n` to take a random sample of that size instead; `seed` fixes which
    ones, so a deleted cache re-downloads the same sample.

    The default is everything because for these organisms the reviewed set is
    not a sample of the proteome, it IS the proteome (6.074 for E. coli K-12,
    7.923 for yeast, 20.616 for human). Holding the whole population removes
    the sampling question from the reference distribution entirely, and it is
    what makes exercise 4b answerable: measuring "how many proteins are
    enough?" against a reference built from the same 400 proteins mostly
    measures the 400.
    """
    term = f'"{organism}"[Organism]'
    if reviewed_only:
        # Swiss-Prot (UniProtKB reviewed), mirrored by NCBI so the same Entrez
        # client works. It is the manually curated half of UniProt: a human
        # reads the literature and writes the entry.
        #
        # The alternative, refseq[filter], is what this used to be, and it is
        # genome-centric: an annotation pipeline emits one protein record per
        # predicted gene per sequenced genome. That gives three problems for
        # a study of amino acid composition:
        #   - redundancy: one RefSeq entry per transcript variant, so a single
        #     human gene contributes five near-identical "proteins";
        #   - unverified entries: every called ORF becomes a record, hence the
        #     "hypothetical protein" and "partial" hits;
        #   - strain duplication: '"Escherichia coli"[Organism]' matches 6.6M
        #     RefSeq proteins, one per strain annotated, against 23k in
        #     Swiss-Prot.
        # Swiss-Prot keeps one entry per gene per organism, with isoforms as
        # annotations inside the entry rather than as separate records.
        term += " AND swissprot[filter]"

    pool = fetch_id_pool(term)
    if n is None or n >= len(pool):
        chosen = pool
    else:
        # Sampled at random, NOT the first n. esearch returns its ids in
        # NCBI's default order, which is by deposit date, so the head of the
        # list is whatever was added most recently. Measured on E. coli K-12:
        # the first 400 entries average 256 residues against 304 for 400 drawn
        # at random, i.e. taking the head biases towards proteins 19% shorter
        # than the organism's.
        chosen = random.Random(seed).sample(pool, n)
    print(f"{organism}: {len(pool)} proteins found, downloading {len(chosen)}")

    records = []
    batch = 200
    for start in range(0, len(chosen), batch):
        handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text",
                               id=",".join(chosen[start:start + batch]))
        records.extend(SeqIO.parse(handle, "fasta"))
        handle.close()
        time.sleep(0.4)   # respects the NCBI limit (3 req/s without an API key)
    return records


ESEARCH_PAGE = 10000   # esearch returns at most 10.000 ids per request


def fetch_id_pool(term):
    """Every id matching `term`, paging esearch until the list is complete.

    The paging is not an optimisation, it is what makes the sampling honest:
    the human reviewed proteome has 20.616 entries and a single esearch call
    would return only the first 10.000 of them. Sampling from that truncated
    pool would still be a draw from the most recently deposited half.
    """
    ids = []
    while True:
        handle = Entrez.esearch(db="protein", term=term,
                                retstart=len(ids), retmax=ESEARCH_PAGE)
        page = Entrez.read(handle)
        handle.close()

        ids.extend(page["IdList"])
        if not page["IdList"] or len(ids) >= int(page["Count"]):
            return ids
        time.sleep(0.4)


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
    downloads them from Swiss-Prot and saves them so they are not fetched again.

    The cache file name carries the source (see a.py), so switching databases
    cannot silently reuse sequences downloaded from the previous one.
    """
    if os.path.exists(cache):
        print(f"Using cached sequences from {cache}")
        with open(cache) as f:
            return json.load(f)
    print(f"Downloading proteins for {organism}...")
    seqs = clean_sequences(fetch_proteins(organism))
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
