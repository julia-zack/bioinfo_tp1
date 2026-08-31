"""4d) Escriba un codigo que dadas dos distribuciones las compare y obtenga la
metrica correspondiente. Utilicela para:

  i)  ver como evoluciona las distribuciones obtenidas en 4a al aumentar el
      tamano de la muestra (evaluar cuando la distribucion deja de cambiar).
  ii) estudiar de manera sistematica las diferencias entre las distribuciones
      obtenidas en 4a para diferentes muestras (comparar distribuciones de
      secuencias naturales y al azar).

Part (i) is exercise 4b: the convergence curve applies the metric to a source
against its own final distribution.

This module is part (ii). It compares the sources against each other, and then
asks the question 4e depends on: how well does the metric separate a single
real protein from a single random one?

Produces:  orf_discrimination.png
"""

import itertools
import os
import random

import matplotlib.pyplot as plt

from sequences import AA, expected_aa_frequencies, generate_random_aa_sequence_from_dna
from stats import freqs, distributions_distance
from exercises.ex4.a import (
    build_natural_source,
    build_sources,
    calculate_data_length,
    download_seqs,
)

HERE = os.path.dirname(os.path.abspath(__file__))

# Amino acid subsets to compare as discriminators. Restricting the metric to a
# few residues is not obviously better than using all 20, so it is measured.
SUBSETS = {
    "all 20": AA,
    "E R K N C": list("ERKNC"),
    "R C": list("RC"),
}

FRAGMENT_LENGTHS = [16, 30, 60, 120, 250, 400]
SAMPLE_SEED = 0


def distance_matrix(distributions):
    """Every pairwise distance between the given distributions."""
    names = list(distributions)
    return {(a, b): distributions_distance(distributions[a], distributions[b])
            for a in names for b in names}


def print_distance_matrix(distributions):
    names = list(distributions)
    matrix = distance_matrix(distributions)
    print("\nTotal variation distance between sources:\n")
    print(" " * 12 + "".join(f"{n:>12}" for n in names))
    for a in names:
        print(f"{a:>12}" + "".join(f"{matrix[(a, b)]:12.3f}" for b in names))


def subset_distance(seq, reference, subset):
    """Distance restricted to a subset of amino acids."""
    frequencies = freqs(seq)
    return 0.5 * sum(abs(frequencies.get(a, 0) - reference.get(a, 0)) for a in subset)


def hit_rate(real_scores, control_scores):
    """Probability that a real fragment scores closer to natural than a random one.

    Takes one real fragment and one random one, and counts how often the real
    one wins, with ties worth half. 0.5 is chance, 1.0 is perfect separation.
    """
    wins = 0.0
    for r in real_scores:
        wins += sum(1 for c in control_scores if r < c)
        wins += 0.5 * sum(1 for c in control_scores if r == c)
    return wins / (len(real_scores) * len(control_scores))


def discrimination_by_length(proteins, natural_freqs, rng, samples=400):
    """Hit rate of each subset, as a function of the fragment length scored.

    An ORF detector scores one ORF at a time, so what matters is not how well
    the metric separates two large corpora but how well it separates two short
    fragments.
    """
    results = {name: [] for name in SUBSETS}
    for length in FRAGMENT_LENGTHS:
        usable = [p for p in proteins if len(p) > length]
        chosen = usable[:samples]
        fragments = []
        for protein in chosen:
            start = rng.randrange(0, len(protein) - length)
            fragments.append(protein[start:start + length])
        controls = [generate_random_aa_sequence_from_dna(length) for _ in chosen]

        for name, subset in SUBSETS.items():
            real_scores = [subset_distance(f, natural_freqs, subset) for f in fragments]
            control_scores = [subset_distance(c, natural_freqs, subset) for c in controls]
            results[name].append(hit_rate(real_scores, control_scores))
    return results


def print_discrimination(results):
    print("\nHit rate separating a real fragment from a random one, by fragment length:\n")
    print(f"{'length':>8}" + "".join(f"{name:>12}" for name in results))
    for i, length in enumerate(FRAGMENT_LENGTHS):
        print(f"{length:>8}" + "".join(f"{results[n][i]:12.3f}" for n in results))


def print_single_residues(proteins, natural_freqs, rng, length=250, samples=300):
    """How much each amino acid discriminates on its own."""
    usable = [p for p in proteins if len(p) > length][:samples]
    fragments = [p[rng.randrange(0, len(p) - length):][:length] for p in usable]
    controls = [generate_random_aa_sequence_from_dna(length) for _ in usable]

    scores = []
    for a in AA:
        real_scores = [subset_distance(f, natural_freqs, [a]) for f in fragments]
        control_scores = [subset_distance(c, natural_freqs, [a]) for c in controls]
        scores.append((hit_rate(real_scores, control_scores), a))

    scores.sort(reverse=True)
    print(f"\nHit rate of each amino acid on its own (fragments of {length}):\n")
    print("  " + "   ".join(f"{a}: {v:.3f}" for v, a in scores[:8]))


# ---------------------------------------------------------------------------
# Leave-one-organism-out: does the choice of amino acids generalise?
# ---------------------------------------------------------------------------

HOLDOUT_LENGTH = 250   # fragment length used for the held-out test


def best_pair(train_seqs, reference, rng, samples=200):
    """The pair of amino acids that separates best on the training organisms.

    Searches all 190 pairs; which one wins is part of what the test checks.
    """
    proteins = [p for seqs in train_seqs.values() for p in seqs
                if len(p) > HOLDOUT_LENGTH]
    rng.shuffle(proteins)
    chosen = proteins[:samples]

    fragments = [p[rng.randrange(0, len(p) - HOLDOUT_LENGTH):][:HOLDOUT_LENGTH]
                 for p in chosen]
    controls = [generate_random_aa_sequence_from_dna(HOLDOUT_LENGTH) for _ in chosen]

    best = None
    for a, b in itertools.combinations(AA, 2):
        pair = [a, b]
        rate = hit_rate([subset_distance(f, reference, pair) for f in fragments],
                        [subset_distance(c, reference, pair) for c in controls])
        if best is None or rate > best[0]:
            best = (rate, pair)
    return best[1], best[0]


def holdout_by_organism(real_seqs, rng, samples=200):
    """For each organism: choose the amino acids on the other two, score on it.

    The reference profile is rebuilt from the two training organisms, so the
    held-out one stays outside both the choice and the profile.
    """
    results = []
    for held_out in real_seqs:
        train_seqs = {name: seqs for name, seqs in real_seqs.items()
                      if name != held_out}
        reference = freqs("".join("".join(seqs) for seqs in train_seqs.values()))

        pair, train_rate = best_pair(train_seqs, reference, rng, samples)

        proteins = [p for p in real_seqs[held_out] if len(p) > HOLDOUT_LENGTH]
        rng.shuffle(proteins)
        chosen = proteins[:samples]
        fragments = [p[rng.randrange(0, len(p) - HOLDOUT_LENGTH):][:HOLDOUT_LENGTH]
                     for p in chosen]
        controls = [generate_random_aa_sequence_from_dna(HOLDOUT_LENGTH)
                    for _ in chosen]

        test_rate = hit_rate(
            [subset_distance(f, reference, pair) for f in fragments],
            [subset_distance(c, reference, pair) for c in controls])

        results.append((held_out, pair, train_rate, test_rate))
    return results


def print_holdout(results):
    print(f"\nLeave-one-out: pair chosen on two organisms, scored on the third "
          f"(fragments of {HOLDOUT_LENGTH}):\n")
    print(f"{'held out':>10} {'pair':>8} {'on training':>13} {'on held out':>13}")
    for held_out, pair, train_rate, test_rate in results:
        print(f"{held_out:>10} {''.join(pair):>8} {train_rate:>13.3f} {test_rate:>13.3f}")


def plot_discrimination(results):
    """Hit rate against fragment length, one line per subset."""
    fig, ax = plt.subplots(figsize=(9, 5))
    colours = ['#d9613d', '#4a6fa5', '#2e8b6f']
    for (name, values), colour in zip(results.items(), colours):
        ax.plot(FRAGMENT_LENGTHS, values, marker='o', color=colour, label=name)

    ax.axhline(0.5, color='#b0b0b0', linestyle='--', linewidth=1, zorder=0)
    ax.text(FRAGMENT_LENGTHS[-1], 0.5, '  chance', va='center', fontsize=8, color='#707070')
    ax.set_xscale('log')
    ax.set_xlabel('Fragment length (residues)')
    ax.set_ylabel('Hit rate: real vs random')
    ax.set_title('How well amino acid composition separates real from random')
    ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'orf_discrimination.png'), dpi=150)


def run():
    rng = random.Random(SAMPLE_SEED)
    # The random controls are generated through the random module, so it needs
    # seeding too, otherwise the reported numbers move on every run.
    random.seed(SAMPLE_SEED)

    real_seqs = download_seqs()
    data_length = calculate_data_length(real_seqs)
    sources = build_sources(real_seqs, data_length)
    natural = build_natural_source(real_seqs, data_length)

    distributions = {name: freqs(seq) for name, seq in sources.items()}
    distributions["Natural"] = freqs(natural)
    distributions["Random DNA"] = expected_aa_frequencies()
    print_distance_matrix(distributions)

    proteins = [p for seqs in real_seqs.values() for p in seqs]
    rng.shuffle(proteins)

    results = discrimination_by_length(proteins, distributions["Natural"], rng)
    print_discrimination(results)
    print_single_residues(proteins, distributions["Natural"], rng)

    print_holdout(holdout_by_organism(real_seqs, rng))

    plot_discrimination(results)
    plt.show()
    print("\nSaved: orf_discrimination.png")
