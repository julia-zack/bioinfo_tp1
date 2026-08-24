"""Plot helpers reused across exercises.

Plots that only one exercise ever draws live in that exercise's own module.
"""

import random

import matplotlib.pyplot as plt

from sequences import AA

BLUE = '#4a6fa5'
ORANGE = '#d1852f'
GREY = '#b0b0b0'
GRID = '#e5e5e5'


def _tidy(ax):
    """Recessive grid and no top/right spines, so the data stands out."""
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def plot_frequencies(frequencies, title, categories=None, expected=None):
    """Bar chart of one frequency distribution.

    `expected` draws a dashed reference line, for the cases where a flat
    expectation applies (a uniform draw over the 20 amino acids).
    """
    categories = categories or AA
    values = [frequencies.get(c, 0) for c in categories]

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(categories, values, color=BLUE, width=0.7)

    if expected is not None:
        ax.axhline(expected, color=GREY, linestyle='--', linewidth=1, zorder=0)
        ax.text(len(categories) - 0.4, expected, f'  expected = {expected:.3f}',
                va='center', fontsize=8, color='#707070')

    ax.set_title(title)
    ax.set_xlabel('Amino acid')
    ax.set_ylabel('Relative frequency')
    _tidy(ax)

    fig.tight_layout()
    plt.show()


def plot_frequencies_comparison(series, title, categories=None):
    """Grouped bar chart comparing several frequency dicts.

    series: dict of {label: frequencies}, one group of bars per category.
    """
    categories = categories or AA
    colors = [BLUE, ORANGE]
    bar_width = 0.8 / len(series)

    fig, ax = plt.subplots(figsize=(11, 4.5))
    for i, (label, frequencies) in enumerate(series.items()):
        values = [frequencies.get(c, 0) for c in categories]
        offsets = [pos - 0.4 + bar_width * (i + 0.5) for pos in range(len(categories))]
        ax.bar(offsets, values, width=bar_width * 0.9,
               color=colors[i % len(colors)], label=label)

    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories)
    ax.set_title(title)
    ax.set_xlabel('Amino acid')
    ax.set_ylabel('Relative frequency')
    _tidy(ax)
    ax.legend(frameon=False)

    fig.tight_layout()
    plt.show()


def plot_orf_sizes(orf_sizes, title):
    """Strip plot of the ORF sizes (in amino acids) per reading frame.

    orf_sizes: dict {frame_label: [sizes]}, like the one exercise 2b builds by
    applying get_orf_sizes() to each frame of get_six_frames_from_nt_seq().
    """
    frame_labels = list(orf_sizes.keys())
    # A separate RNG for the jitter, so it does not disturb the random module
    # state used to generate the sequences.
    jitter_rng = random.Random()

    fig, ax = plt.subplots(figsize=(9, 4.5))

    for i, label in enumerate(frame_labels):
        sizes = orf_sizes[label]
        if not sizes:
            continue
        xs = [i + jitter_rng.uniform(-0.15, 0.15) for _ in sizes]
        ax.scatter(xs, sizes, color=BLUE, alpha=0.7,
                   edgecolors='white', linewidths=0.5, zorder=3)

    ax.set_xticks(range(len(frame_labels)))
    ax.set_xticklabels(frame_labels)
    ax.set_title(title)
    ax.set_xlabel('Reading frame')
    ax.set_ylabel('ORF size (amino acids, incl. start/stop)')
    _tidy(ax)

    fig.tight_layout()
    plt.show()


def plot_orf_size_comparison(real_sizes, control_sizes, title):
    """Compare two ORF size distributions.

    Left: density, to show the shape. Right: P(size >= k) on a log scale,
    where a geometric distribution appears as a straight line.
    """
    colors = {'real': BLUE, 'control': ORANGE}
    fig, (ax_hist, ax_tail) = plt.subplots(1, 2, figsize=(12, 4.5))

    bins = range(0, max(max(real_sizes), max(control_sizes)) + 10, 5)
    for label, sizes in (('real', real_sizes), ('control', control_sizes)):
        ax_hist.hist(sizes, bins=bins, density=True, histtype='step',
                     linewidth=2, color=colors[label], label=label)

    ax_hist.set_title('Density', fontsize=10)
    ax_hist.set_xlabel('ORF size (amino acids)')
    ax_hist.set_ylabel('Density')
    ax_hist.legend(frameon=False)

    for label, sizes in (('real', real_sizes), ('control', control_sizes)):
        ordered = sorted(sizes)
        survival = [1 - i / len(ordered) for i in range(len(ordered))]
        ax_tail.plot(ordered, survival, linewidth=2, color=colors[label], label=label)

    ax_tail.set_yscale('log')
    ax_tail.set_title('Tail: P(size >= k)', fontsize=10)
    ax_tail.set_xlabel('ORF size (amino acids)')
    ax_tail.set_ylabel('P(size >= k)')
    ax_tail.legend(frameon=False)

    for ax in (ax_hist, ax_tail):
        _tidy(ax)

    fig.suptitle(title)
    fig.tight_layout()
    plt.show()
