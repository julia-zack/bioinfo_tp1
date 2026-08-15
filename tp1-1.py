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

import random
import numpy as np
import matplotlib.pyplot as plt

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

def generate_random_nt_sequence(length):
    nucleotides = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(length))

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

# def get_six_aa_seq_from_nt_seq(nt_seq):
#     sequences = []
#     frames =
#     for frame in range(3):
#         aa_seq = get_aa_seq_from_nt_seq(nt_seq[frame:])
#         sequences.append(aa_seq)
#     return sequences


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
    """Análisis escrito, sobre una corrida con n = 1000 para ambas secuencias."""
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


def main():
    # ex_1a_i()
    # ex_1a_ii()
    # ex_1b_i()
    # ex_1b_ii()
    ex_1b_iii()


if __name__ == "__main__":
    main()
