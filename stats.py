"""Distribution statistics: amino acid frequencies over sets of sequences,
a distance between distributions, and the convergence analysis of exercise 4a."""

import random
from collections import Counter

import numpy as np

from sequences import AA


# ---------------------------------------------------------------------------
# Amino acid statistics
# ---------------------------------------------------------------------------

def natural_statistics(seqs):
    """Return (aa_frequencies, list_of_lengths) for a list of sequences."""
    lengths = [len(s) for s in seqs]
    counts = Counter()
    for s in seqs:
        counts.update(s)
    total = sum(counts.values())
    frequencies = {aa: counts[aa] / total for aa in AA}
    return frequencies, lengths


def freqs(seq_or_list):
    """Amino acid frequencies of a sequence (str) or a list of sequences."""
    if isinstance(seq_or_list, str):
        return freqs_seq(seq_or_list)
    return freqs_list(list(seq_or_list))


def freqs_seq(seq):
    frequencies, length = natural_statistics([seq])
    return frequencies


def freqs_list(seqs_list):
    frequencies, length = natural_statistics(seqs_list)
    return frequencies


# ---------------------------------------------------------------------------
# Metric for comparing two distributions (exercise 4c)
# ---------------------------------------------------------------------------

def distributions_distance(dist1, dist2):
    """Total variation distance between two aa distributions, from 0 to 1.

    Chosen as the metric for 4c: it is half the sum of the absolute
    differences, so it reads directly as "fraction of probability mass that
    would have to move to turn one distribution into the other".
    """
    diffs_sum = 0.0
    for aa in AA:
        freq_dist1 = dist1.get(aa, 0)   # 0 if that aa is not in dist1
        freq_dist2 = dist2.get(aa, 0)   # 0 if that aa is not in dist2
        diff = abs(freq_dist1 - freq_dist2)
        diffs_sum += diff
    return 0.5 * diffs_sum


# ---------------------------------------------------------------------------
# Residue sampling
# ---------------------------------------------------------------------------

def concatenate_til_n(seqs, n):
    """Concatenate proteins until reaching ~n residues (truncated to n)."""
    res = ""
    for seq in seqs:
        res += seq
        if len(res) >= n:
            break
    return res[:n]


NUM_SAMPLE_SIZES = 10   # how many sample sizes to try


def sample_sizes(res_length):
    """~NUM_SAMPLE_SIZES sample sizes spaced on a log scale, from 50 up to
    res_length. DYNAMIC: more data -> larger sizes. The last value is
    res_length, so the curve ends up using ALL the residues."""
    minimum = 50
    points = np.geomspace(minimum, res_length, NUM_SAMPLE_SIZES)  # log-spaced
    sizes = []
    for p in points:
        size = int(round(p))
        if size not in sizes:   # skip duplicates introduced by rounding
            sizes.append(size)
    return sizes


# ---------------------------------------------------------------------------
# Convergence of the distribution
# ---------------------------------------------------------------------------

REPS = 5    # repetitions, to smooth out the randomness


def smooth_curve(res_list, grid, final_distribution):
    """Repeat REPS times (shuffling the order) and average the curves point by
    point, to get a smooth curve that does not depend on any particular
    ordering."""
    curves_per_rep = []
    for rep in range(REPS):
        random.shuffle(res_list)
        shuffle_seq = "".join(res_list)

        curve = []
        for R in grid:
            prefix = shuffle_seq[:R]              # first R residues
            prefix_distribution = freqs(prefix)
            distance = distributions_distance(prefix_distribution, final_distribution)
            curve.append(distance)
        curves_per_rep.append(curve)

    average_curve = list(np.mean(curves_per_rep, axis=0))
    return average_curve


def convergence_curve(res):
    """Grid of sample sizes and its curve (distance to the final distribution)
    for one residue source."""
    res_list = list(res)              # list of residues, so it can be shuffled
    final_distribution = freqs(res)   # distribution using ALL the residues
    res_length = len(res_list)
    grid = sample_sizes(res_length)   # dynamic sample sizes (up to res_length)

    curve = smooth_curve(res_list, grid, final_distribution)
    return grid, curve


THRESHOLD = 0.01   # "stabilised": distance below this value


def find_N(grid, ds):
    """First size in the grid whose distance falls below the threshold."""
    for R, d in zip(grid, ds):
        if d < THRESHOLD:
            return R
    return grid[-1]   # did not stabilise within the range tested
