
# 4a) Compare distribucion de aminoacidos de secuencia de proteinas al azar vs. secuencias de proteinas reales
#
# Design, in two steps:
# Step 1
# Pick the real proteins to use: E. coli and human.
# When does the distribution stabilise?
# -> N is the number of residues beyond which the frequency estimate stops changing.
#
# Step 2
# Compare the three distributions (random, E. coli, human) at sample sizes
# << N, ~ N and >> N, to show graphically that below N the estimate is
# unreliable and above N it is stable.

# Run from this folder:  python3 main.py
# Produces:  4a_convergencia.png  and  4a_regimenes.png
# The real sequences ship pre-downloaded in data/. If you delete those files,
# the first run fetches them from NCBI again and re-caches them there.
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # to import the utils package (inside ej4/)

from utils.utils import (
    AA,
    THRESHOLD,
    freqs,
    random_residues,
    concatenate_til_n,
    convergence_curve,
    find_N,
    get_real_sequences,
)

# real organisms to compare (name -> (NCBI name, cache file))
REAL_ORGANISMS = {
    "E. coli": ("Escherichia coli", os.path.join(HERE, "utils", "data", "sequences_ecoli.json")),
    "Human":   ("Homo sapiens",     os.path.join(HERE, "utils", "data", "sequences_human.json")),
}
COLOURS = {"Random": "#c1666b", "E. coli": "#4a6fa5", "Human": "#6a9955"}


def download_seqs():
    """Load (or download and cache) the real sequences for each organism."""
    real_seqs = {}
    for name, (org, cache) in REAL_ORGANISMS.items():
        real_seqs[name] = get_real_sequences(org, cache)
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


def build_sources(real_seqs, length):
    """Build, for each source, a string of EXACTLY `length` residues: one
    random sequence, and the real ones truncated to `length`."""
    sources = {"Random": random_residues(length)}
    for name, seqs in real_seqs.items():
        full_sequence = "".join(seqs)
        sources[name] = full_sequence[:length]
    return sources


def analyze_convergence(sources):
    """For each source compute: its convergence curve, its N, and its final
    distribution (using every residue)."""
    curves = {}
    N_by_source = {}
    final_distributions = {}
    for name, res in sources.items():
        grid, ds = convergence_curve(res)
        curves[name] = (grid, ds)
        N_by_source[name] = find_N(grid, ds)
        final_distributions[name] = freqs(res)
    return curves, N_by_source, final_distributions


def print_N_by_source(N_by_source):
    print("\nTamano de estabilizacion N (residuos) por fuente:")
    for name, N in N_by_source.items():
        print(f"  {name}: N = {N}")


def representative_N(curves):
    """N where the AVERAGE of every source's curve stabilises. The three
    converge together, so a single N summarises "when it is enough for all of
    them" (and avoids the noise of a per-source N jumping between runs)."""
    all_curves = []
    for grid, ds in curves.values():
        all_curves.append(ds)
    mean_curve = list(np.mean(all_curves, axis=0))
    common_grid = curves["Random"][0]   # every source shares the same grid
    return find_N(common_grid, mean_curve)


def mark_N_on_xaxis(ax, N):
    """Add N as one more tick on the X axis, formatted like the 10^k ticks but
    highlighted (bold + a fourth colour) and only once."""
    # 1) build the "decade" ticks (10^2, 10^3, ...) that fall inside the range
    lo, hi = ax.get_xlim()
    min_exponent = int(np.ceil(np.log10(lo)))
    max_exponent = int(np.floor(np.log10(hi)))
    decades = []
    for k in range(min_exponent, max_exponent + 1):
        decades.append(10 ** k)
    # 2) add N to the tick list (without duplicating it if it lands on a decade)
    ticks = sorted(set(decades + [N]))
    ax.set_xticks(ticks)
    # 3) label: N as a plain number; the decades as 10^k
    labels = []
    for t in ticks:
        if t == N:
            labels.append(str(int(N)))
        else:
            exponent = int(round(np.log10(t)))
            labels.append(r'$10^{%d}$' % exponent)
    ax.set_xticklabels(labels)
    # 4) highlight N's label (bold + colour)
    for lbl, t in zip(ax.get_xticklabels(), ticks):
        if t == N:
            lbl.set_color('#d1462f')       # a fourth colour, so it stands out
            lbl.set_fontweight('bold')


def plot_convergence(curves, N):
    """Figure 1: one convergence curve per source, with N marked."""
    fig, ax = plt.subplots(figsize=(9, 5))
    for name, (grid, ds) in curves.items():
        ax.plot(grid, ds, marker='o', color=COLOURS[name], label=name)
    ax.axhline(THRESHOLD, color='#b0b0b0', linestyle='--', linewidth=1, zorder=0,
               label=f'threshold = {THRESHOLD}')
    ax.axvline(N, color='#d1462f', linestyle=':', linewidth=1.3)
    ax.set_xscale('log')
    ax.set_xlabel('Residues analysed')
    ax.set_ylabel('Distance to the final distribution')
    ax.set_title('4a - Step 1: when does the aa distribution stabilise?')
    ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False)

    mark_N_on_xaxis(ax, N)

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, '4a_convergencia.png'), dpi=150)


def choose_regimes(N, data_length):
    """Three sample sizes to compare: far below N, around N, and all the data
    available (far above N)."""
    R_small = 1000                # far below N: still noisy (try 500 for more noise)
    R_med = min(N, data_length)   # around N
    R_large = data_length         # all the data available (>> N)
    return [
        (f'<< N  (R = {R_small})', R_small),
        (f'~ N   (R = {R_med})',   R_med),
        (f'>> N  (R = {R_large})', R_large),
    ]


def plot_regimes(real_seqs, final_distributions, N, data_length):
    """Figure 2: aa distribution of the three sources at three sample sizes
    (<< N, ~ N, >> N)."""
    regimes = choose_regimes(N, data_length)

    # sort the aa by their frequency in E. coli, so the plot reads better
    order = sorted(AA, key=lambda a: final_distributions["E. coli"].get(a, 0), reverse=True)
    base_positions = np.arange(len(order))   # one x position per amino acid
    bar_width = 0.27

    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)
    for ax, (regime_label, R) in zip(axes, regimes):
        # aa distribution of each source using R residues
        samples = {
            "Random":  freqs(random_residues(R)),
            "E. coli": freqs(concatenate_til_n(real_seqs["E. coli"], R)),
            "Human":   freqs(concatenate_til_n(real_seqs["Human"], R)),
        }
        # one batch of bars per source, shifted left/centre/right
        for j, name in enumerate(["Random", "E. coli", "Human"]):
            distribution = samples[name]
            heights = []
            for a in order:
                heights.append(distribution.get(a, 0))
            positions = base_positions + (j - 1) * bar_width
            ax.bar(positions, heights, bar_width, label=name, color=COLOURS[name])

        ax.axhline(1 / len(AA), color='#b0b0b0', linestyle='--', linewidth=1, zorder=0)
        ax.set_title(f'Sample {regime_label} residues')
        ax.set_ylabel('Frequency')
        ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)

    axes[0].legend(frameon=False, ncol=3)
    axes[-1].set_xticks(base_positions)
    axes[-1].set_xticklabels(order)
    axes[-1].set_xlabel('Amino acid (ordered by frequency in E. coli)')
    fig.suptitle('4a - Step 2: aa distribution at different sample sizes')
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, '4a_regimenes.png'), dpi=150)


def main():
    real_seqs = download_seqs()
    data_length = calculate_data_length(real_seqs)
    sources = build_sources(real_seqs, data_length)

    # Step 1: when does the aa distribution stabilise
    curves, N_by_source, final_distributions = analyze_convergence(sources)
    print_N_by_source(N_by_source)
    N = representative_N(curves)
    print(f"\nN representativo (promedio de las tres fuentes): {N}")
    plot_convergence(curves, N)

    # Step 2: compare the distributions at << N, ~ N, >> N
    plot_regimes(real_seqs, final_distributions, N, data_length)

    plt.show()
    print("\nGuardados: 4a_convergencia.png y 4a_regimenes.png")


if __name__ == "__main__":
    main()
