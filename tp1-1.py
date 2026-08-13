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

aa = ['V','L','T','K','W','H','F','I','R','M','A','P','G','S','C','N','Q','Y','D','E']
# aa_weights = [0.00, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06, 0.04]
aa_weights = [1, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]

weighted_aa = np.random.choice(aa, size=len(aa), p=aa_weights)

def generate_weighted_sequence(length):
    return ''.join(random.choice(weighted_aa) for _ in range(length))


def generate_sequence(length):
    return ''.join(random.choice(aa) for _ in range(length))

def calculate_aa_frequencies(sequence):
    counts = {}
    for residue in sequence:
        if residue in counts:
            counts[residue] += 1
        else:
            counts[residue] = 1
    return {residue: count / len(sequence) for residue, count in counts.items()}


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

def get_complementary_sequence(nt_seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in nt_seq)

def get_aa_from_codon(codon):
    codon = codon.upper()
    if codon not in codon_table:
        raise ValueError(f"Codon {codon} is not valid.")
    return codon_table[codon]


def get_aa_seq_from_nt_seq(nt_seq):
    aa_seq = ''
    for i in range(0, len(nt_seq), 3):
        codon = nt_seq[i:i+3]
        aa_seq += get_aa_from_codon(codon)

    return aa_seq

# def get_six_aa_seq_from_nt_seq(nt_seq):
#     sequences = []
#     frames =
#     for frame in range(3):
#         aa_seq = get_aa_seq_from_nt_seq(nt_seq[frame:])
#         sequences.append(aa_seq)
#     return sequences



input_length = int(input("Enter the length of the amino acid sequence: "))
random_sequence = generate_sequence(input_length)
print (random_sequence)
print (calculate_aa_frequencies(random_sequence))

weighted_sequence = generate_weighted_sequence(input_length)
print(weighted_sequence)
print(calculate_aa_frequencies(weighted_sequence))