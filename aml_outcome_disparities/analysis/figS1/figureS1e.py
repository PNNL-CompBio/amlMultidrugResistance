"""Plots figure S1e: Mutation Waterfall"""

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

    # Split Black and White patients
    black_mutations = (meta.loc[
        meta.loc[:, "Race"] == "Black",
        "ASXL1":"ZRSR2"
    ] == "Mutant").astype(int)
    white_mutations = (meta.loc[
        meta.loc[:, "Race"] == "White",
        "ASXL1":"ZRSR2"
    ] == "Mutant").astype(int)

    # Trim to mutations present in at least 5 patients in either Race
    black_mutations = black_mutations.loc[
        :,
        np.logical_or(
            black_mutations.sum(axis=0) > 5,
            white_mutations.sum(axis=0) > 5
        )
    ]
    white_mutations = white_mutations.loc[:, black_mutations.columns]

    # Setup figure
    fig, axes = get_setup(
        2,
        2,
        fig_params={
            "figsize": (8, 6),
            "width_ratios": (4, 1)
        }
    )

    # Iterate through Black, White patients
    for row_index, (race, dataset) in enumerate(
        zip(
            ["Black", "White"],
            [black_mutations, white_mutations]
        )
    ):
        # Get waterfall and bar plot axes
        waterfall_ax = axes[row_index, 0]
        bar_ax = axes[row_index, 1]

        # Sort mutations by frequency
        dataset = dataset.loc[
            :,
            dataset.sum(axis=0).sort_values(ascending=False).index
        ]
        dataset = dataset.sort_values(
            by=list(dataset.columns),
            ascending=False
        ).T

        # Plot waterfall
        sns.heatmap(
            dataset,
            ax=waterfall_ax,
            cmap="Greys",
            linewidths=0.1,
            linecolor="tab:grey",
            cbar=False
        )

        # Format waterfall plot, label ticks and axes
        waterfall_ax.set_xticks([])
        waterfall_ax.set(
            title=f"{race} Patients",
            yticks=np.arange(0.5, dataset.shape[0]),
            yticklabels=dataset.index
        )

        # Plot number of patients with each mutation
        mutation_sums = dataset.sum(axis=1)
        bar_ax.barh(
            np.arange(0.5, dataset.shape[0], 1),
            mutation_sums
        )

        # Format bar plot ticks and limits
        bar_ax.set(
            ylim=(0, dataset.shape[0]),
            yticks=[]
        )

        # Label bottom barplot
        if row_index == 1:
            bar_ax.set(
                xlabel="Number of patients\nwith mutation"
            )
        else:
            bar_ax.set_xticks([])

        # Turn off frame, invert y-axis to match heatmap
        bar_ax.set_frame_on(False)
        bar_ax.yaxis.set_inverted(True)

        # Denote mutation percentage
        lim = bar_ax.get_xlim()
        lim = lim[1] - lim[0]
        for offset, gene in enumerate(mutation_sums.index):
            bar_ax.text(
                dataset.loc[gene].sum() + lim * 0.01,
                0.5 + offset,
                ha="left",
                ma="left",
                va="center",
                s=f"{
                    round(
                        mutation_sums.loc[gene] / dataset.shape[1] * 100,
                        1
                    )
                }%"
            )

    return fig


if __name__ == "__main__":
    make_figure()
