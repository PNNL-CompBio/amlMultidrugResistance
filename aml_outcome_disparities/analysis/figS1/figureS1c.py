"""Plots figure S1c: Mutation Proportions"""

from os.path import abspath, dirname

import numpy as np

from pilot.data_import import import_meta
from pilot.figures.figure_setup import get_setup

FILE_DIR = dirname(abspath(__file__))
RACE_COLORS = {"Black": "tab:red", "White": "tab:purple"}


def make_figure():
    # Import meta data
    meta = import_meta()
    meta = meta.drop("FLT3", axis=1)

    # Split Black and White patients
    black_mutations = meta.loc[
        meta.loc[:, "Race"] == "Black",
        "ASXL1":"ZRSR2"
    ] == "Mutant"
    white_mutations = meta.loc[
        meta.loc[:, "Race"] == "White",
        "ASXL1":"ZRSR2"
    ] == "Mutant"

    # Get sum of patients with each mutation
    black_mutations = black_mutations.sum(axis=0)
    white_mutations = white_mutations.sum(axis=0)

    # Reorder mutations, most to least present in Black patients
    black_mutations = black_mutations.sort_values(ascending=False)
    white_mutations = white_mutations.loc[black_mutations.index]

    # Trim to mutations present in at least 10 patients in either Race
    black_mutations = black_mutations.loc[
        np.logical_or(
            black_mutations > 10,
            white_mutations > 10
        )
    ]
    white_mutations = white_mutations.loc[black_mutations.index]

    # Setup figure
    fig, ax = get_setup(
        1,
        1,
        fig_params={
            "figsize": (6, 3)
        }
    )

    # Plot mutation counts
    ax.bar(
        np.arange(0, len(black_mutations) * 3, 3),
        black_mutations,
        color=RACE_COLORS["Black"],
        width=0.9,
        label="Black Patients",
        zorder=3
    )
    ax.bar(
        np.arange(1, len(white_mutations) * 3, 3),
        white_mutations,
        color=RACE_COLORS["White"],
        width=0.9,
        label="White Patients",
        zorder=3
    )
    ax.legend()

    for index, mutation in enumerate(black_mutations.index):
        ax.text(
            3 * index,
            black_mutations.loc[mutation] + 0.1,
            ha="center",
            ma="center",
            va="bottom",
            s=black_mutations.loc[mutation]
        )
        ax.text(
            3 * index + 1,
            white_mutations.loc[mutation] + 0.1,
            ha="center",
            ma="center",
            va="bottom",
            s=white_mutations.loc[mutation]
        )


    ax.set(
        xlim=(-1, 3 * len(black_mutations) - 1),
        ylim=(0, 125),
        xticks=np.arange(
            0.5, len(black_mutations) * 3, 3
        ),
        ylabel="Number of Patients with Mutation"
    )
    ax.set_xticklabels(
        black_mutations.index,
        rotation=45,
        ha="right",
        ma="right",
        va="top"
    )
    ax.grid(True, zorder=0)

    return fig


if __name__ == "__main__":
    make_figure()
