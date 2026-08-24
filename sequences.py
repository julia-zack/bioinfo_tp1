"""Sequence handling shared by every exercise: the genetic code, random
sequence generation, translation, reading frames and ORFs."""

import random

from Bio.Seq import Seq
from Bio.Data.IUPACData import protein_letters

# The 20 standard amino acids, in alphabetical order (taken from Biopython so
# the list is not hand-typed).
AA = list(protein_letters)

NUCLEOTIDES = ['A', 'T', 'C', 'G']

# Standard genetic code, DNA alphabet. '*' = stop codon.
CODON_TABLE = {
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


# ---------------------------------------------------------------------------
# Random sequence generation
# ---------------------------------------------------------------------------

def generate_random_aa_sequence(length):
    """Random protein sequence, every residue equally likely."""
    return ''.join(random.choices(AA, k=length))


def generate_weighted_aa_sequence(length, aa_weights):
    """Random protein sequence with one weight per residue (same order as AA).

    The weights do not need to add up to 1.
    """
    return ''.join(random.choices(AA, weights=aa_weights, k=length))


def generate_random_nt_sequence(length, gc=0.5):
    """Random DNA sequence with the requested GC content (gc=0.5 => uniform)."""
    weights = [(1 - gc) / 2, (1 - gc) / 2, gc / 2, gc / 2]
    return ''.join(random.choices(NUCLEOTIDES, weights=weights, k=length))


def gc_fraction(nt_seq):
    """Fraction of G and C, ignoring ambiguity codes (N, R, Y, ...)."""
    counts = {base: str(nt_seq).upper().count(base) for base in 'ACGT'}
    total = sum(counts.values())
    return (counts['G'] + counts['C']) / total if total else 0.0


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------

def calculate_aa_frequencies(sequence):
    """Relative frequency of each symbol present in the sequence."""
    counts = {}
    for residue in sequence:
        if residue in counts:
            counts[residue] += 1
        else:
            counts[residue] = 1
    return {residue: count / len(sequence) for residue, count in counts.items()}


# ---------------------------------------------------------------------------
# Complement and translation
# ---------------------------------------------------------------------------

def get_complementary_sequence(nt_seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in nt_seq)


def get_reverse_complementary_sequence(nt_seq):
    # Reversed so the result reads 5' -> 3', as translation requires.
    return get_complementary_sequence(nt_seq)[::-1]


def get_aa_from_codon(codon):
    codon = codon.upper()
    if codon not in CODON_TABLE:
        raise ValueError(f"Codon {codon} is not valid.")
    return CODON_TABLE[codon]


def get_aa_seq_from_nt_seq(nt_seq):
    """Translate one reading frame with our own codon table.

    Kept alongside the Biopython version on purpose: exercise 1 asks for the
    translation to be written by hand.
    """
    aa_seq = ''
    # Trailing 1-2 nt do not form a codon, so they are discarded.
    last_full_codon_end = len(nt_seq) - len(nt_seq) % 3
    for i in range(0, last_full_codon_end, 3):
        codon = nt_seq[i:i+3]
        aa_seq += get_aa_from_codon(codon)

    return aa_seq


# ---------------------------------------------------------------------------
# Reading frames and ORFs
# ---------------------------------------------------------------------------

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
    start and stop codons. ORFs do not overlap: once one is open, later M's are
    absorbed into it, and a run that never reaches a stop is discarded.
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


def collect_orf_sizes(nt_seq):
    """Every ORF from the six reading frames, in a single list."""
    frames = get_six_frames_from_nt_seq(nt_seq)
    return [size for aa_seq in frames.values() for size in get_orf_sizes(aa_seq)]
