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


def generate_random_nt_sequence(length):
    """Random DNA sequence, each base equally likely."""
    return ''.join(random.choices(NUCLEOTIDES, k=length))


def generate_random_aa_sequence_from_dna(length):
    """Random protein obtained by translating random DNA.

    The residues carry no information, but they follow the genetic code, so
    amino acids with more codons (L, R, S) come out more often. Stop codons
    are dropped, because a protein has no internal stops.
    """
    residues = ''
    while len(residues) < length:
        # 3 nt per codon, plus a margin for the ~4.7% that translate to a stop
        missing = length - len(residues)
        nt_seq = generate_random_nt_sequence(int(missing * 3 * 1.2) + 3)
        residues += get_aa_seq_from_nt_seq(nt_seq).replace('*', '')
    return residues[:length]


def expected_aa_frequencies():
    """Exact amino acid frequencies of the random-DNA null model.

    With every base equally likely each codon is equally likely, so this is
    (number of codons)/61 per amino acid.
    """
    base_prob = {base: 0.25 for base in NUCLEOTIDES}

    frequencies = {}
    for codon, residue in CODON_TABLE.items():
        if residue == '*':
            continue
        probability = base_prob[codon[0]] * base_prob[codon[1]] * base_prob[codon[2]]
        frequencies[residue] = frequencies.get(residue, 0) + probability

    total = sum(frequencies.values())
    return {residue: p / total for residue, p in frequencies.items()}


def stop_codon_probability():
    """Probability that a random codon is a stop: 3 of the 64 codons."""
    stops = sum(1 for residue in CODON_TABLE.values() if residue == '*')
    return stops / len(CODON_TABLE)


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

    Exercise 1 asks for the translation to be written by hand;
    get_six_frames_from_nt_seq() uses Biopython instead.
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


def find_orfs(nt_seq, min_length=1):
    """Every ORF in the six frames, with frame, position, length and protein.

    Positions are on the forward strand, so a reverse-strand ORF still points
    at the region of `nt_seq` it came from.
    """
    orfs = []
    total = len(nt_seq)
    frames = get_six_frames_from_nt_seq(nt_seq)

    for label, aa_seq in frames.items():
        offset = int(label[1]) - 1          # 0, 1 or 2 nucleotides into the strand
        reverse = label.startswith('-')

        residue_start = None
        for i, residue in enumerate(aa_seq):
            if residue == 'M' and residue_start is None:
                residue_start = i
            elif residue == '*' and residue_start is not None:
                peptide = aa_seq[residue_start:i]      # without the stop
                if len(peptide) >= min_length:
                    # back to nucleotide coordinates on the strand that was read
                    strand_start = offset + residue_start * 3
                    strand_end = offset + (i + 1) * 3   # includes the stop codon
                    if reverse:
                        start, end = total - strand_end, total - strand_start
                    else:
                        start, end = strand_start, strand_end
                    orfs.append({
                        'frame': label,
                        'start': start,
                        'end': end,
                        'length': len(peptide),
                        'protein': peptide,
                    })
                residue_start = None

    return orfs


def collect_orf_sizes(nt_seq):
    """Every ORF from the six reading frames, in a single list."""
    frames = get_six_frames_from_nt_seq(nt_seq)
    return [size for aa_seq in frames.values() for size in get_orf_sizes(aa_seq)]
