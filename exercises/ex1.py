"""Objetivo 1) Generar y analizar secuencias de Proteina al azar.

1a) i.  Generar una secuencia de aminoacidos al azar.
    ii. Generar una secuencia al azar con distinta probabilidad por aminoacido.
1b) i.  Calcular la distribucion de aminoacidos, generada desde aminoacidos y
        desde nucleotidos.
    ii. Graficar ambas distribuciones.
    iii. Analizar y discutir los resultados obtenidos.
"""

from sequences import (
    AA,
    calculate_aa_frequencies,
    generate_random_aa_sequence,
    generate_random_nt_sequence,
    generate_weighted_aa_sequence,
    get_aa_seq_from_nt_seq,
)
from plots import plot_frequencies


def a_i():
    """Random amino acid sequence, every residue equally likely."""
    input_length = int(input("Enter the length of the amino acid sequence: "))
    random_sequence = generate_random_aa_sequence(input_length)
    print(random_sequence)


def a_ii():
    """Random amino acid sequence with a custom weight per residue."""
    input_length = int(input("Enter the length of the amino acid sequence: "))
    # aa_weights = [0.05] * len(AA)
    aa_weights = [1] + [0.00] * (len(AA) - 1)
    weighted_sequence = generate_weighted_aa_sequence(input_length, aa_weights)
    print(weighted_sequence)


def b_i():
    """Amino acid frequencies, from a random protein and from random DNA.

    Returns the two frequency dicts and their lengths so 1b-ii can plot them.
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


def b_ii():
    """Bar charts for both frequency distributions of 1b-i."""
    (from_aa_frequencies, aa_input_length,
     from_nt_seq_frequencies, nt_input_length) = b_i()

    plot_frequencies(from_aa_frequencies,
                     f'Amino acid frequencies from random amino acid sequence (n = {aa_input_length}, uniform)',
                     expected=1 / len(AA))

    plot_frequencies(from_nt_seq_frequencies,
                     f'Amino acid frequencies from random nucleotide sequence (n = {nt_input_length}, uniform)',
                     categories=AA + ['*'])


def b_iii():
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
