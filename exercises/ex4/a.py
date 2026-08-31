"""4a) Compare Distribucion de aminoacidos de secuencia de proteinas al azar
vs. secuencias de proteinas reales.

Sources compared, all with exactly the same number of residues:
  - Random DNA: information-free residues, obtained by translating random DNA,
    so they still obey the genetic code.
  - E. coli, yeast, human: real proteins downloaded from NCBI.
  - Natural: the three organisms pooled in equal parts. 4e needs a single
    reference distribution to score an ORF against, and this is it.

Produces:  aa_distribution.png  and  aa_natural_vs_random.png
"""

import os
import random

import numpy as np
import matplotlib.pyplot as plt

from sequences import AA, expected_aa_frequencies, generate_random_aa_sequence_from_dna
from stats import freqs
from ncbi import cache_path, get_real_sequences

HERE = os.path.dirname(os.path.abspath(__file__))

# real organisms to compare (name -> (organism name, cache file))
#
# The cache files are named after the source database. The previous
# sequences_*.json held RefSeq downloads, and reusing those under the same name
# would have silently mixed the two sources.
# E. coli is pinned to K-12: '"Escherichia coli"[Organism]' matches every
# sequenced strain (23.300 Swiss-Prot entries), which reintroduces the strain
# redundancy the switch to Swiss-Prot was meant to remove. K-12 is one
# organism's proteome (6.074). Yeast and human are already single-organism.
# Swiss-Prot assigns entries to the strain level, not the substrain, so
# "...str. K-12 substr. MG1655" matches nothing.
REAL_ORGANISMS = {
    "E. coli": ("Escherichia coli K-12",     cache_path("swissprot_ecoli_k12.json")),
    "Yeast":   ("Saccharomyces cerevisiae",  cache_path("swissprot_yeast.json")),
    "Human":   ("Homo sapiens",              cache_path("swissprot_human.json")),
}

# Null model: residues obtained by translating random DNA. Information-free,
# but subject to the genetic code, so a stricter control than uniform
# residues. The uniform case is the 1/20 line drawn on each figure.
RANDOM_SOURCES = ["Random DNA"]

COLOURS = {
    "Random DNA": "#d9613d",   # warm: the control
    "Natural":    "#2e8b6f",   # the pooled real profile (figure 3 only)
    "E. coli":    "#97baed",   # cool: the three real organisms
    "Yeast":      "#7de2ef",
    "Human":      "#98efaf",
}

# Which proteins end up in a truncated sample must not depend on the order
# NCBI happened to return them in, so the lists are shuffled once, with a
# fixed seed to keep runs reproducible.
SAMPLE_SEED = 0


def run():
    real_seqs = download_seqs()
    # The smallest proteome sizes the random source and each organism's
    # share of the pooled profile.
    data_length = calculate_data_length(real_seqs)
    sources = build_sources(real_seqs, data_length)
    natural = build_natural_source(real_seqs, data_length)

    print("\nSources compared:")
    for name, seq in sources.items():
        print(f"  {name:<11} {len(seq):>7} residues")
    print(f"  {'Natural':<11} {len(natural):>7} residues  (pooled reference)")

    plot_distributions(sources)
    plot_natural_vs_random(natural, real_seqs, data_length)

    plt.show()
    print("\nSaved: aa_distribution.png and aa_natural_vs_random.png")


def download_seqs():
    """Load (or download and cache) the real sequences for each organism.

    The protein list of each organism is shuffled, because later steps keep
    only the first `data_length` residues: shuffling is what makes that a
    random sample of the organism and not a slice of NCBI's ordering.
    """
    shuffler = random.Random(SAMPLE_SEED)
    real_seqs = {}
    for name, (org, cache) in REAL_ORGANISMS.items():
        seqs = list(get_real_sequences(org, cache))
        shuffler.shuffle(seqs)
        real_seqs[name] = seqs
    return real_seqs


def calculate_data_length(real_seqs):
    """Residue count of the SMALLEST real source. Every source is levelled to
    this total (TOT), so the comparison is fair: all of them are measured with
    the same amount of data."""
    real_totals = []
    for seqs in real_seqs.values():
        total_residues = 0
        for protein in seqs:
            total_residues += len(protein)
        real_totals.append(total_residues)
    return min(real_totals)


def build_sources(real_seqs, random_length):
    """Every organism's whole proteome, plus a random source of `random_length`.

    """
    sources = {"Random DNA": generate_random_aa_sequence_from_dna(random_length)}
    for name, seqs in real_seqs.items():
        sources[name] = "".join(seqs)
    return sources


def build_natural_source(real_seqs, per_organism):
    """One pooled "natural" string, `per_organism` residues from each.

    4e scores an ORF against a single reference distribution, and this is it.
    Equal parts keep the largest proteome from dominating the profile.
    """
    return "".join("".join(seqs)[:per_organism] for seqs in real_seqs.values())


def plot_distributions(sources):
    """Figure 1: amino acid distribution of every source, using all the data.

    Sorted by mean frequency across the organisms, so the random source can be
    read against them residue by residue.
    """
    names = list(sources)
    distributions = {name: freqs(seq) for name, seq in sources.items()}
    organisms = [n for n in names if n not in RANDOM_SOURCES]
    order = sorted(
        AA,
        key=lambda a: sum(distributions[n].get(a, 0) for n in organisms) / len(organisms),
        reverse=True)

    positions = np.arange(len(order))
    bar_width = 0.8 / len(names)

    fig, ax = plt.subplots(figsize=(13, 5))
    for j, name in enumerate(names):
        heights = [distributions[name].get(a, 0) for a in order]
        ax.bar(positions - 0.4 + bar_width * (j + 0.5), heights,
               bar_width * 0.9, color=COLOURS[name],
               label=f'{name}  (n = {len(sources[name]):,})')

    ax.axhline(1 / len(AA), color='#b0b0b0', linestyle='--', linewidth=1, zorder=0)
    ax.set_xticks(positions)
    ax.set_xticklabels(order)
    ax.set_xlabel('AA')
    ax.set_ylabel('Relative frequency of residues')
    ax.set_title('Amino acid distribution: random control vs real proteins')
    ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False, ncol=len(names))

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'aa_distribution.png'), dpi=150)


def plot_natural_vs_random(natural, real_seqs, data_length):
    """Figure 3: the pooled natural profile against the random control.

    The vertical bars show the min-max range across the three organisms, so
    the spread of "natural" is visible next to the gap to random. If that
    range overlaps the random bar for an amino acid, that residue carries no
    discriminating power on its own.
    """
    natural_freqs = freqs(natural)
    # The control here is a reference, not a sample: use its exact frequencies
    # so no sampling noise of its own enters the comparison.
    random_freqs = expected_aa_frequencies()
    per_organism = {name: freqs("".join(seqs)[:data_length])
                    for name, seqs in real_seqs.items()}

    order = sorted(AA, key=lambda a: natural_freqs.get(a, 0), reverse=True)
    positions = np.arange(len(order))
    bar_width = 0.38

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(positions - bar_width / 2,
           [natural_freqs.get(a, 0) for a in order],
           bar_width, color=COLOURS['Natural'],
           label=f'Natural (pooled)  (n = {len(natural):,})')
    ax.bar(positions + bar_width / 2,
           [random_freqs.get(a, 0) for a in order],
           bar_width, color=COLOURS['Random DNA'],
           label='Random DNA (exact)')

    # min-max across organisms, drawn over the natural bar
    for x, a in zip(positions, order):
        values = [d.get(a, 0) for d in per_organism.values()]
        ax.plot([x - bar_width / 2] * 2, [min(values), max(values)],
                color='#333333', linewidth=1.2, zorder=4)

    ax.axhline(1 / len(AA), color='#b0b0b0', linestyle='--', linewidth=1, zorder=0)
    ax.set_xticks(positions)
    ax.set_xticklabels(order)
    ax.set_xlabel('AA')
    ax.set_ylabel('Relative frequency of residues')
    ax.set_title('Natural profile vs random control '
                 '(vertical bars: range across E. coli, yeast and human)')
    ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'aa_natural_vs_random.png'), dpi=150)
