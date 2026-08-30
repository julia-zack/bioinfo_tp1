"""4b) Analice como cambian las distribuciones al aumentar el tamano de la
secuencia "al azar" analizada, y al incrementar el numero de secuencias reales
analizadas. Cuando es suficiente?

Ayuda: seaborn.histplot() o matplotlib.pyplot.hist() pueden ayudar con las
comparaciones.

La consigna nombra DOS ejes, y no son el mismo: se puede crecer en cantidad de
RESIDUOS o en cantidad de PROTEINAS. Este archivo mide los dos.

Step 1, axis 1 (residues)
When does the distribution stabilise? N is the number of residues beyond which
the frequency estimate stops changing: the first sample size whose distance to
the source's own final distribution falls below THRESHOLD. Applies to every
source, the random control included.

Step 1, axis 2 (number of real sequences)
K is the number of whole proteins beyond which an organism's distribution stops
changing. Only the real sources have this axis: the random control has no
"sequences", it is one arbitrarily long string.

Step 2
Show the distributions at sample sizes << N, ~ N and >> N, so it is visible
that below N the estimate is noise and above N it is stable.

Produces:  aa_stabilized.png  and  aa_sample_sizes.png
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from sequences import AA
from stats import (
    THRESHOLD,
    convergence_curve,
    count_convergence_curve,
    find_N,
    freqs,
)
from exercises.ex4.a import (
    COLOURS,
    RANDOM_SOURCES,
    build_sources,
    calculate_data_length,
    download_seqs,
)

HERE = os.path.dirname(os.path.abspath(__file__))

def run():
    real_seqs = download_seqs()
    data_length = calculate_data_length(real_seqs)
    sources = build_sources(real_seqs, data_length)

    # Step 1, axis 1: how many RESIDUES until the aa distribution stabilises
    curves, N_by_source, final_distributions = analyze_convergence(sources)
    print_N_by_source(N_by_source)
    N = representative_N(N_by_source)
    print(f"\nRepresentative N (enough for every source): {N}")

    # Step 1, axis 2: how many PROTEINS until an organism's distribution
    # stabilises. Uses the full protein lists, not the truncated sources.
    count_curves, K_by_organism = analyze_protein_count(real_seqs)
    print_K_by_organism(K_by_organism, real_seqs)
    K = representative_N(K_by_organism)
    print(f"\nRepresentative K (enough for every organism): {K}")

    plot_convergence(curves, N, count_curves, K_by_organism, K)

    # Step 2: compare the distributions at << N, ~ N, >> N
    plot_regimes(sources, final_distributions, N, data_length)

    plt.show()
    print("\nSaved: aa_stabilized.png and aa_sample_sizes.png")


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
    print("\nStabilisation size N (residues) per source:")
    for name, N in N_by_source.items():
        print(f"  {name}: N = {N}")


def representative_N(N_by_source):
    """The size where ALL sources have already stabilised: the largest
    per-source value, i.e. the slowest one.

    Used on both axes: residues (N) on axis 1, whole proteins (K) on axis 2.
    """
    return max(N_by_source.values())


def analyze_protein_count(real_seqs):
    """Second axis: for each organism, its curve against the NUMBER of whole
    proteins analysed, and the K at which it crosses the threshold.

    The random control is left out on purpose: it is generated as one string of
    residues, so "number of sequences" is not defined for it.
    """
    count_curves = {}
    K_by_organism = {}
    for name, proteins in real_seqs.items():
        counts, curve = count_convergence_curve(proteins)
        count_curves[name] = (counts, curve)
        K_by_organism[name] = find_N(counts, curve)
    return count_curves, K_by_organism


def print_K_by_organism(K_by_organism, real_seqs):
    """K per organism, and what K is worth in residues.

    The residue equivalent is the point of the comparison: K proteins carry far
    more residues than the N of axis 1, because residues inside one protein are
    correlated and a whole protein is therefore worth less than the same number
    of residues drawn independently.
    """
    print("\nStabilisation size K (whole proteins) per organism:")
    for name, K in K_by_organism.items():
        proteins = real_seqs[name]
        mean_length = sum(len(p) for p in proteins) / len(proteins)
        print(f"  {name}: K = {K} of {len(proteins)} proteins available "
              f"(~{K * mean_length:,.0f} residues, mean length {mean_length:.0f})")


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


def tidy(ax):
    """Recessive grid and no top/right spines, like the rest of the 4a/4b plots."""
    ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def plot_convergence(curves, N, count_curves, K_by_organism, K):
    """Figure 1: the two axes of the consigna, side by side.

    Left: distance to the final distribution as the number of RESIDUES grows,
    one curve per source, with N marked. Right: the same, as the number of
    whole PROTEINS grows, one curve per organism, with K marked. Both panels
    highlight their stabilisation point the same way, so the two axes can be
    read off the figure directly.
    """
    fig, (ax_residues, ax_proteins) = plt.subplots(1, 2, figsize=(13, 5))
    plot_residue_axis(ax_residues, curves, N)
    plot_protein_axis(ax_proteins, count_curves, K_by_organism, K)

    fig.suptitle('AA distribution stabilization: residues vs number of proteins')
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'aa_stabilized.png'), dpi=150)


def plot_residue_axis(ax, curves, N):
    """Left panel: one convergence curve per source, with N marked."""
    for name, (grid, ds) in curves.items():
        # Dashed for the two controls, solid for the real organisms, so the
        # two groups stay apart even in black and white.
        style = '--' if name in RANDOM_SOURCES else '-'
        ax.plot(grid, ds, marker='o', linestyle=style,
                color=COLOURS[name], label=name)
    ax.axhline(THRESHOLD, color='#b0b0b0', linestyle='--', linewidth=1, zorder=0,
               label=f'threshold = {THRESHOLD}')
    ax.axvline(N, color='#d1462f', linestyle=':', linewidth=1.3)
    ax.set_xscale('log')
    ax.set_xlabel('Residues analysed')
    ax.set_ylabel('Distance to the final distribution')
    ax.set_title(f'Axis 1: residues (enough at ~{N})', fontsize=10)
    tidy(ax)
    ax.legend(frameon=False)

    mark_N_on_xaxis(ax, N)


def plot_protein_axis(ax, count_curves, K_by_organism, K):
    """Right panel: one curve per organism against the number of whole
    proteins, with the representative K marked.

    Each organism's own K goes in the legend; only the representative one is
    drawn on the axis, the same way the left panel marks N.
    """
    for name, (counts, ds) in count_curves.items():
        ax.plot(counts, ds, marker='o', linestyle='-',
                color=COLOURS[name], label=f'{name}  (K = {K_by_organism[name]})')
    ax.axhline(THRESHOLD, color='#b0b0b0', linestyle='--', linewidth=1, zorder=0,
               label=f'threshold = {THRESHOLD}')
    ax.axvline(K, color='#d1462f', linestyle=':', linewidth=1.3)
    ax.set_xscale('log')
    ax.set_xlabel('Number of whole proteins analysed')
    ax.set_ylabel('Distance to the final distribution')
    ax.set_title(f'Axis 2: number of real sequences (enough at ~{K})', fontsize=10)
    tidy(ax)
    ax.legend(frameon=False)

    mark_N_on_xaxis(ax, K)


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


def plot_regimes(sources, final_distributions, N, data_length):
    """Figure 2: aa distribution of every source at three sample sizes
    (<< N, ~ N, >> N).

    Each sample is a prefix of the string built in build_sources(), so it uses
    exactly the same data as the convergence curves.
    """
    regimes = choose_regimes(N, data_length)
    names = list(sources)

    # sort the aa by their frequency in E. coli, so the plot reads better
    order = sorted(AA, key=lambda a: final_distributions["E. coli"].get(a, 0), reverse=True)
    base_positions = np.arange(len(order))   # one x position per amino acid
    bar_width = 0.8 / len(names)

    fig, axes = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
    for ax, (regime_label, R) in zip(axes, regimes):
        # one batch of bars per source, spread around each amino acid position
        for j, name in enumerate(names):
            distribution = freqs(sources[name][:R])
            heights = []
            for a in order:
                heights.append(distribution.get(a, 0))
            positions = base_positions - 0.4 + bar_width * (j + 0.5)
            ax.bar(positions, heights, bar_width * 0.9,
                   label=name, color=COLOURS[name])

        ax.axhline(1 / len(AA), color='#b0b0b0', linestyle='--', linewidth=1, zorder=0)
        ax.set_title(f'Sample {regime_label} residues')
        ax.set_ylabel('Relative frequency of residues')
        # sharex hides the tick labels on the upper panels; put them back, so
        # each panel can be read on its own.
        ax.set_xticks(base_positions)
        ax.set_xticklabels(order)
        ax.tick_params(labelbottom=True)
        ax.yaxis.grid(True, color='#e5e5e5', linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)

    axes[0].legend(frameon=False, ncol=len(names))
    axes[-1].set_xlabel('AA')
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'aa_sample_sizes.png'), dpi=150)
