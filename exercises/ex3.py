"""Objetivo 3) Secuencias reales bajadas de una base de datos.

3b) Bajar una secuencia de GenBank y mostrar sus datos.
3c) Bajar 100 secuencias al azar y comparar la distribucion de tamanos de ORF
    contra un control aleatorio del mismo largo.
"""

from Bio import Entrez, SeqIO

from sequences import collect_orf_sizes, generate_random_nt_sequence
from ncbi import fetch_random_records
from plots import plot_orf_size_comparison


def b():
    """One GenBank record: id, description and sequence."""
    handle = Entrez.efetch(db="nucleotide", id="AY851612",
                           rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    print(record.id, len(record.seq))
    print(record.description)
    print(record.seq)


def c():
    """100 real NCBI sequences against an equivalent random control."""
    records = fetch_random_records(n=100, seed=0)
    print(f"{len(records)} sequences, mean length "
          f"{sum(len(r.seq) for r in records) / len(records):.0f} nt")

    real_sizes = []
    control_sizes = []
    for record in records:
        real_sizes.extend(collect_orf_sizes(record.seq))
        # Control: same length as the real sequence, bases drawn at random.
        control = generate_random_nt_sequence(len(record.seq))
        control_sizes.extend(collect_orf_sizes(control))

    for label, sizes in (("real", real_sizes), ("control", control_sizes)):
        ordered = sorted(sizes)
        print(f"{label:>8}: {len(sizes):5d} ORFs  "
              f"mean {sum(sizes) / len(sizes):6.1f}  "
              f"median {ordered[len(ordered) // 2]:4d}  "
              f"max {max(sizes):5d}")

    plot_orf_size_comparison(
        real_sizes, control_sizes,
        f'ORF sizes: {len(records)} NCBI sequences vs random control of the same length')
