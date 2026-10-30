"""Plots figure S1c: Mutation Waterfall"""

from os.path import abspath, dirname

import numpy as np
import seaborn as sns

from pilot.data_import import import_meta
from pilot.figures.figure_setup import get_setup

FILE_DIR = dirname(abspath(__file__))
RACE_COLORS = {"Black": "tab:red", "White": "tab:purple"}


def make_figure():
    # Import meta data
    meta = import_meta()
    meta = meta.drop("FLT3", axis=1)
    meta = meta.loc[
        ~meta.index.str.endswith("Bridge"),
        :
    ]

    # Split Black and White patients
    mutations = (meta.loc[
        :,
        "ASXL1":"ZRSR2"
    ] == "Mutant").astype(int)

    # Trim to mutations present in at least 10 patients in either Race
    mutations = mutations.loc[
        :,
        mutations.sum(axis=0) > 10
    ]

    # Denote WES vs. clinical calls
    mutations.loc[meta.loc[mutations.index, "WES_call"], :] *= 2
    
    # Setup figure
    fig, axes = get_setup(
        1,
        2,
        fig_params={
            "figsize": (8, 3),
            "width_ratios": (4, 1)
        }
    )

    # Get waterfall and bar plot axes
    waterfall_ax = axes[0]
    bar_ax = axes[1]

    # Get sum of mutations
    mutation_sums = (mutations > 0).sum(axis=0).sort_values(ascending=False)

    # Sort mutations by frequency
    mutations = mutations.loc[
        :,
        mutation_sums.index
    ]
    mutations = mutations.sort_values(
        by=list(mutations.columns),
        ascending=False,
        axis=0
    ).T

    # Plot waterfall
    sns.heatmap(
        mutations,
        ax=waterfall_ax,
        cmap="Greys",
        linewidths=0.1,
        linecolor="tab:grey",
        cbar=False
    )

    # Format waterfall plot, label ticks and axes
    waterfall_ax.set_xticks([])
    waterfall_ax.set(
        ylabel="Mutation",
        yticks=np.arange(0.5, mutations.shape[0]),
        yticklabels=mutations.index
    )

    # Plot number of patients with each mutation
    bar_ax.barh(
        np.arange(0.5, len(mutation_sums), 1),
        mutation_sums / mutations.shape[1],
        zorder=3
    )

    # Format bar plot ticks and limits
    bar_ax.set(
        xticks=np.arange(0, 0.5, 0.1),
        xlim=(0, 0.5),
        ylim=(0, mutations.shape[0]),
        yticks=[]
    )
    bar_ax.grid(True, zorder=0)

    # Turn off frame, invert y-axis to match heatmap
    bar_ax.set_frame_on(False)
    bar_ax.yaxis.set_inverted(True)

    # Denote mutation percentage
    lim = bar_ax.get_xlim()
    lim = lim[1] - lim[0]
    for offset, gene in enumerate(mutation_sums.index):
        bar_ax.text(
            (mutation_sums.loc[gene] / mutations.shape[1]) + (lim * 0.01),
            0.5 + offset,
            ha="left",
            ma="left",
            va="center",
            s=int(mutation_sums.loc[gene])
        )

    return fig


if __name__ == "__main__":
    make_figure()
