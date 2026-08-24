"""Objetivo 2) Marcos de lectura y ORFs.

2a) Traducir una secuencia de ADN al azar y su reversa complementaria, y
    comparar las distribuciones de aminoacidos de ambas.
2b) Obtener los seis marcos de lectura y determinar los ORFs de cada uno.
"""

from Bio.Seq import Seq

from sequences import (
    AA,
    calculate_aa_frequencies,
    generate_random_nt_sequence,
    get_orf_sizes,
    get_six_frames_from_nt_seq,
)
from plots import plot_frequencies_comparison, plot_orf_sizes


def a():
    """Forward strand vs reverse complement, translated and compared."""
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
        categories=AA + ['*'])


def b():
    """ORF sizes in each of the six reading frames."""
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
