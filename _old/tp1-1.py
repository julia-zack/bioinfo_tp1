# Valina (Val, V)
# Leucina (Leu, L)
# Treonina (Thr, T)
# Lisina (Lys, K)
# Triptófano (Trp, W)
# Histidina (His, H)
# Fenilalanina (Phe, F)
# Isoleucina (Ile, I)
# Arginina (Arg, R)
# Metionina (Met, M)
# Alanina (Ala, A)
# Prolina (Pro, P)
# Glicina (Gly, G)
# Serina (Ser, S)
# Cisteína (Cys, C)
# Asparagina (Asn, N)
# Glutamina (Gln, Q)
# Tirosina (Tyr, Y)
# Ácido aspártico (Asp, D)
# Ácido glutámico (Glu, E)

import os
import random
import time
import numpy as np
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio import Entrez, SeqIO

Entrez.email = 'arielzingman@gmail.com'

###############################################################################
# Constants
###############################################################################

aa = ['V','L','T','K','W','H','F','I','R','M','A','P','G','S','C','N','Q','Y','D','E']

# Standard genetic code, DNA alphabet. '*' = stop codon.
codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',

    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

###############################################################################
# Helppers
###############################################################################

def generate_weighted_sequence(length, aa_weights):
    return ''.join(random.choices(aa, weights=aa_weights, k=length))

def generate_random_aa_sequence(length):
    return generate_weighted_sequence(length, [1 / len(aa)] * len(aa))

def gc_fraction(nt_seq):
    """Fraction of G and C, ignoring ambiguity codes (N, R, Y, ...)."""
    counts = {base: str(nt_seq).upper().count(base) for base in 'ACGT'}
    total = sum(counts.values())
    return (counts['G'] + counts['C']) / total if total else 0.0

def generate_random_nt_sequence(length, gc=0.5):
    """Random sequence with the requested GC content (gc=0.5 => uniform)."""
    nucleotides = ['A', 'T', 'C', 'G']
    weights = [(1 - gc) / 2, (1 - gc) / 2, gc / 2, gc / 2]
    return ''.join(random.choices(nucleotides, weights=weights, k=length))

def calculate_aa_frequencies(sequence):
    counts = {}
    for residue in sequence:
        if residue in counts:
            counts[residue] += 1
        else:
            counts[residue] = 1
    return {residue: count / len(sequence) for residue, count in counts.items()}

def get_complementary_sequence(nt_seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in nt_seq)

def get_reverse_complementary_sequence(nt_seq):
    # Reversed so the result reads 5' -> 3', as translation requires.
    return get_complementary_sequence(nt_seq)[::-1]

def get_aa_from_codon(codon):
    codon = codon.upper()
    if codon not in codon_table:
        raise ValueError(f"Codon {codon} is not valid.")
    return codon_table[codon]

def get_aa_seq_from_nt_seq(nt_seq):
    aa_seq = ''
    # Trailing 1-2 nt do not form a codon, so they are discarded.
    last_full_codon_end = len(nt_seq) - len(nt_seq) % 3
    for i in range(0, last_full_codon_end, 3):
        codon = nt_seq[i:i+3]
        aa_seq += get_aa_from_codon(codon)

    return aa_seq

def plot_frequencies(frequencies, title, categories=None, expected=None):
    categories = categories or aa
    values = [frequencies.get(c, 0) for c in categories]

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(categories, values, color='#4a6fa5', width=0.7)

    if expected is not None:
        ax.axhline(expected, color='#b0b0b0', linestyle='--', linewidth=1, zorder=0)
        ax.text(len(categories) - 0.4, expected, f'  expected = {expected:.3f}',
                va='center', fontsize=8, color='#707070')

    ax.set_title(title)
    ax.set_xlabel('Amino acid')
    ax.set_ylabel('Relative frequency')
    ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)

    fig.tight_layout()
    plt.show()   # or fig.savefig('freqs.png', dpi=150) for the report

def plot_frequencies_comparison(series, title, categories=None):
    """Grouped bar chart comparing several frequency dicts.

    series: dict of {label: frequencies}, one group of bars per category.
    """
    categories = categories or aa
    colors = ['#4a6fa5', '#d1852f']
    bar_width = 0.8 / len(series)

    fig, ax = plt.subplots(figsize=(11, 4.5))
    for i, (label, frequencies) in enumerate(series.items()):
        values = [frequencies.get(c, 0) for c in categories]
        offsets = [pos - 0.4 + bar_width * (i + 0.5) for pos in range(len(categories))]
        ax.bar(offsets, values, width=bar_width * 0.9,
               color=colors[i % len(colors)], label=label)

    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories)
    ax.set_title(title)
    ax.set_xlabel('Amino acid')
    ax.set_ylabel('Relative frequency')
    ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False)

    fig.tight_layout()
    plt.show()


def get_six_frames_from_nt_seq(nt_seq):
    """Translate the sequence in all six reading frames.

    Keys are '+1'..'+3' for the forward strand and '-1'..'-3' for the reverse
    complement; both strands are read 5' -> 3'. Each frame drops the trailing
    1-2 nucleotides that do not complete a codon.
    """
    sequences = {}
    forward = Seq(str(nt_seq))
    strands = [('+', forward), ('-', forward.reverse_complement())]

    for sign, strand in strands:
        for frame in range(3):
            shifted = strand[frame:]
            # Trim explicitly, otherwise translate() warns about a partial codon.
            trimmed = shifted[:len(shifted) - len(shifted) % 3]
            sequences[f'{sign}{frame + 1}'] = str(trimmed.translate())

    return sequences


def get_orf_sizes(seq):
    """Return a list of the sizes of all ORFs in the sequence.

    An ORF is defined as a sequence that starts with a start codon (M) and ends
    with a stop codon (*). The size is the number of amino acids, including the
    start and stop codons.
    """
    orf_sizes = []
    current_orf_size = 0
    in_orf = False

    for aa in seq:
        if aa == 'M' and not in_orf:
            in_orf = True
            current_orf_size = 1  # Start counting from the start codon
        elif aa == '*' and in_orf:
            current_orf_size += 1  # Count the stop codon
            orf_sizes.append(current_orf_size)
            in_orf = False
            current_orf_size = 0
        elif in_orf:
            current_orf_size += 1

    return orf_sizes

def plot_orf_sizes(orf_sizes, title):
    """Strip plot of the ORF sizes (in amino acids) per reading frame.

    orf_sizes: dict {frame_label: [sizes]}, like the one ex_2b() builds by
    applying get_orf_sizes() to each frame of get_six_frames_from_nt_seq().
    """
    frame_labels = list(orf_sizes.keys())
    # A separate RNG for the jitter, so it does not disturb the random module
    # state used to generate the sequences.
    jitter_rng = random.Random()

    fig, ax = plt.subplots(figsize=(9, 4.5))

    for i, label in enumerate(frame_labels):
        sizes = orf_sizes[label]
        if not sizes:
            continue
        xs = [i + jitter_rng.uniform(-0.15, 0.15) for _ in sizes]
        ax.scatter(xs, sizes, color='#4a6fa5', alpha=0.7,
                   edgecolors='white', linewidths=0.5, zorder=3)

    ax.set_xticks(range(len(frame_labels)))
    ax.set_xticklabels(frame_labels)
    ax.set_title(title)
    ax.set_xlabel('Reading frame')
    ax.set_ylabel('ORF size (amino acids, incl. start/stop)')
    ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)

    fig.tight_layout()
    plt.show()

###############################################################################
# NCBI: sampling real sequences
###############################################################################

NCBI_CACHE = "ncbi_sample.fasta"

# A random sample of this query, not of all GenBank: the query is part of the method.
NCBI_QUERY = '"Homo sapiens"[Organism] AND biomol_mrna[PROP] AND 500:5000[SLEN]'


def fetch_random_records(n=100, query=NCBI_QUERY, cache_path=NCBI_CACHE, seed=None):
    """Download n sequences sampled at random from the results of a query.

    Writes a multi-FASTA to cache_path and reuses it on later runs, so the
    requests to NCBI are not repeated.
    """
    if os.path.exists(cache_path):
        return list(SeqIO.parse(cache_path, "fasta"))

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

    SeqIO.write(records, cache_path, "fasta")
    return records


def collect_orf_sizes(nt_seq):
    """Every ORF from the six reading frames, in a single list."""
    frames = get_six_frames_from_nt_seq(nt_seq)
    return [size for aa_seq in frames.values() for size in get_orf_sizes(aa_seq)]


def plot_orf_size_comparison(real_sizes, control_sizes, title):
    """Compare two ORF size distributions.

    Left: density, to show the shape. Right: P(size >= k) on a log scale,
    where a geometric distribution appears as a straight line.
    """
    colors = {'real': '#4a6fa5', 'control': '#d1852f'}
    fig, (ax_hist, ax_tail) = plt.subplots(1, 2, figsize=(12, 4.5))

    bins = range(0, max(max(real_sizes), max(control_sizes)) + 10, 5)
    for label, sizes in (('real', real_sizes), ('control', control_sizes)):
        ax_hist.hist(sizes, bins=bins, density=True, histtype='step',
                     linewidth=2, color=colors[label], label=label)

    ax_hist.set_title('Density', fontsize=10)
    ax_hist.set_xlabel('ORF size (amino acids)')
    ax_hist.set_ylabel('Density')
    ax_hist.legend(frameon=False)

    for label, sizes in (('real', real_sizes), ('control', control_sizes)):
        ordered = sorted(sizes)
        survival = [1 - i / len(ordered) for i in range(len(ordered))]
        ax_tail.plot(ordered, survival, linewidth=2, color=colors[label], label=label)

    ax_tail.set_yscale('log')
    ax_tail.set_title('Tail: P(size >= k)', fontsize=10)
    ax_tail.set_xlabel('ORF size (amino acids)')
    ax_tail.set_ylabel('P(size >= k)')
    ax_tail.legend(frameon=False)

    for ax in (ax_hist, ax_tail):
        ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)

    fig.suptitle(title)
    fig.tight_layout()
    plt.show()


###############################################################################
# Exercises. Each one runs on its own; comment out the calls in main() that you
# do not want to run.
###############################################################################

def ex_1a_i():
    """Random amino acid sequence, every residue equally likely."""
    input_length = int(input("Enter the length of the amino acid sequence: "))
    random_sequence = generate_random_aa_sequence(input_length)
    print(random_sequence)


def ex_1a_ii():
    """Random amino acid sequence with a custom weight per residue."""
    input_length = int(input("Enter the length of the amino acid sequence: "))
    # aa_weights = [0.05] * len(aa)
    aa_weights = [1] + [0.00] * (len(aa) - 1)
    weighted_sequence = generate_weighted_sequence(input_length, aa_weights)
    print(weighted_sequence)


def ex_1b_i():
    """Amino acid frequencies, from a random protein and from random DNA.

    Returns the two frequency dicts and their lengths so 1B-II can plot them.
    """
    # from aminoacids
    aa_input_length = int(input("Enter the length of the amino acid sequence: "))
    random_aa_sequence = generate_random_aa_sequence(aa_input_length)
    print(random_aa_sequence)
    from_aa_frequencies = calculate_aa_frequencies(random_aa_sequence)
    print(from_aa_frequencies)

    # from nucleotides
    nt_input_length = int(input("Enter the length of the nucleotide sequence: "))
    random_nt_sequence = generate_random_nt_sequence(nt_input_length)
    print(random_nt_sequence)
    aa_seq = get_aa_seq_from_nt_seq(random_nt_sequence)
    print(aa_seq)
    from_nt_seq_frequencies = calculate_aa_frequencies(aa_seq)
    print(from_nt_seq_frequencies)

    return (from_aa_frequencies, aa_input_length,
            from_nt_seq_frequencies, nt_input_length)


def ex_1b_ii():
    """Bar charts for both frequency distributions of 1B-I."""
    (from_aa_frequencies, aa_input_length,
     from_nt_seq_frequencies, nt_input_length) = ex_1b_i()

    plot_frequencies(from_aa_frequencies,
                     f'Amino acid frequencies from random amino acid sequence (n = {aa_input_length}, uniform)',
                     expected=1 / len(aa))

    plot_frequencies(from_nt_seq_frequencies,
                     f'Amino acid frequencies from random nucleotide sequence (n = {nt_input_length}, uniform)',
                     categories=aa + ['*'])


def ex_1b_iii():
    """Written analysis of a run with n = 1000 for both sequences.

    The printed text stays in Spanish: it is an answer for the report.
    """
    print("""
Comparación entre las dos secuencias aleatorias, basada en una corrida del
programa con el mismo tamaño para ambas secuencias (1000 aminoácidos y 3000
nucleótidos, es decir 1000 codones).

Secuencia aleatoria de aminoácidos (uniforme), n = 1000:
{'V': 0.043, 'N': 0.043, 'I': 0.052, 'W': 0.053, 'E': 0.061, 'F': 0.054,
 'R': 0.044, 'P': 0.05,  'S': 0.051, 'L': 0.055, 'D': 0.038, 'A': 0.055,
 'H': 0.037, 'C': 0.056, 'Q': 0.047, 'K': 0.062, 'M': 0.051, 'T': 0.055,
 'Y': 0.04,  'G': 0.053}

Todos los valores estan entre 0.037 y 0.062, dispersos alrededor de 1/20 = 0.05.

Secuencia aleatoria de nucleótidos traducida, n = 1000 codones:
{'A': 0.063, 'G': 0.067, 'L': 0.102, 'F': 0.027, 'R': 0.091, 'I': 0.045,
 'C': 0.035, 'E': 0.034, 'V': 0.054, 'T': 0.065, 'M': 0.012, '*': 0.043,
 'H': 0.036, 'W': 0.02,  'D': 0.035, 'S': 0.091, 'N': 0.025, 'Q': 0.035,
 'P': 0.06,  'K': 0.028, 'Y': 0.032}

Claramente no uniforme, con un rango mucho más amplio (0.012 a 0.102). La
frecuencia esperada de cada aminoácido es (número de codones que lo
representa)/64.

Tabla comparativa de frecuencias esperadas y observadas hecha por Claude:
┌──────────┬───────────────────────────┬──────────┬─────────────────────┐
│ Codones  │        Aminoácidos        │ Esperado │      Observado      │
├──────────┼───────────────────────────┼──────────┼─────────────────────┤
│ 6        │ L, R, S                   │ 0.094    │ 0.102, 0.091, 0.091 │
├──────────┼───────────────────────────┼──────────┼─────────────────────┤
│ 4        │ V, P, T, A, G             │ 0.063    │ 0.054–0.067         │
├──────────┼───────────────────────────┼──────────┼─────────────────────┤
│ 2        │ C, D, E, F, H, K, N, Q, Y │ 0.031    │ 0.025–0.036         │
├──────────┼───────────────────────────┼──────────┼─────────────────────┤
│ 1        │ M, W                      │ 0.016    │ 0.012, 0.020        │
├──────────┼───────────────────────────┼──────────┼─────────────────────┤
│ 3 (stop) │ *                         │ 0.047    │ 0.043               │
└──────────┴───────────────────────────┴──────────┴─────────────────────┘

Conclusión: aunque los nucleótidos se eligen aleatoria y uniformemente, la
traducción a aminoácidos no es uniforme, porque la cantidad de codones que
representa a cada aminoácido es diferente.
""")


def ex_2a():
    """Exercise 2a."""

    nt_input_length = int(input("Enter the length of the nucleotide sequence: "))
    random_nt_sequence = generate_random_nt_sequence(nt_input_length)

    nt_seq = Seq(random_nt_sequence)
    aa_seq = nt_seq.translate()
    print(aa_seq)

    aa_seq_frequencies = calculate_aa_frequencies(str(aa_seq))
    aa_seq_frequencies = dict(sorted(aa_seq_frequencies.items()))
    print(aa_seq_frequencies)

    rev_comp = nt_seq.reverse_complement()
    aa_seq_rev_comp = rev_comp.translate()
    print(aa_seq_rev_comp)
    aa_seq_rev_comp_frequencies = calculate_aa_frequencies(str(aa_seq_rev_comp))
    aa_seq_rev_comp_frequencies = dict(sorted(aa_seq_rev_comp_frequencies.items()))
    print(aa_seq_rev_comp_frequencies)

    plot_frequencies_comparison(
        {'Forward strand': aa_seq_frequencies,
         'Reverse complement': aa_seq_rev_comp_frequencies},
        f'Amino acid frequencies: forward strand vs reverse complement (n = {nt_input_length} nt)',
        categories=aa + ['*'])


def ex_2b():
    """Exercise 2b."""
    nt_input_length = int(input("Enter the length of the nucleotide sequence: "))
    random_nt_sequence = generate_random_nt_sequence(nt_input_length)

    nt_seq = Seq(random_nt_sequence)
    frames = get_six_frames_from_nt_seq(nt_seq)

    orf_sizes = {}
    for label, aa_seq in frames.items():
        orf_sizes[label] = get_orf_sizes(aa_seq)

    for label, sizes in orf_sizes.items():
        print(f'{label}: {sizes}')

    plot_orf_sizes(orf_sizes,
                f'ORF sizes by reading frame (n = {nt_input_length} nt)')

def ex_3b():

    handle = Entrez.efetch(db="nucleotide", id="AY851612", rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    print(record.id, len(record.seq))
    print(record.description)
    print(record.seq)

def ex_3c():
    """100 real NCBI sequences against an equivalent random control."""
    records = fetch_random_records(n=100, seed=0)
    print(f"{len(records)} sequences, mean length "
          f"{sum(len(r.seq) for r in records) / len(records):.0f} nt")

    real_sizes = []
    control_sizes = []
    for record in records:
        real_sizes.extend(collect_orf_sizes(record.seq))
        # Control: same length and same %GC as the real sequence, but random.
        control = generate_random_nt_sequence(len(record.seq), gc=gc_fraction(record.seq))
        control_sizes.extend(collect_orf_sizes(control))

    for label, sizes in (("real", real_sizes), ("control", control_sizes)):
        ordered = sorted(sizes)
        print(f"{label:>8}: {len(sizes):5d} ORFs  "
              f"mean {sum(sizes) / len(sizes):6.1f}  "
              f"median {ordered[len(ordered) // 2]:4d}  "
              f"max {max(sizes):5d}")

    plot_orf_size_comparison(
        real_sizes, control_sizes,
        f'ORF sizes: {len(records)} NCBI sequences vs random control (same length and %GC)')


def main():
    # ex_1a_i()
    # ex_1a_ii()
    # ex_1b_i()
    # ex_1b_ii()
    # ex_1b_iii()

    # ex_2a()
    # ex_2b()

    # ex_3b()
    ex_3c()

if __name__ == "__main__":
    main()
