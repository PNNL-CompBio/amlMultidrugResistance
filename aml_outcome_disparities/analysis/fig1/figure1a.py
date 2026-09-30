"""Plots figure 1a: Mutation and Sex across Race."""

import numpy as np

from pilot.data_import import import_meta
from pilot.figures.figure_setup import get_setup


def make_figure():
    # Import metadata
    meta = import_meta()

    # Setup figure
    fig, axes = get_setup(
        2,
        2,
        fig_params={
            "figsize": (4, 4)
        }
    )
    for ax_row, race in enumerate(["Black", "White"]):
        # Trim to race
        _meta = meta.loc[meta.loc[:, "Race"] == race, :]

        # Get counts of NPM1, FLT3-ITD, double mutants
        mutation_counts = _meta.loc[:, "NPM1"].replace(
            {
                "Mutant": "NPM1",
                "WT": "Wild Type"
            }
        )
        mutation_counts.loc[
            np.logical_and(
                _meta.loc[:, "NPM1"] == "Mutant",
                _meta.loc[:, "FLT3_ITD"] == "Mutant"
            )
        ] = "Both"
        mutation_counts.loc[
            np.logical_and(
                _meta.loc[:, "NPM1"] != "Mutant",
                _meta.loc[:, "FLT3_ITD"] == "Mutant"
            )
        ] = "FLT3-ITD"
        mutation_counts = mutation_counts.value_counts()

        # Pie chart for mutations
        ax = axes[ax_row, 0]
        ax.pie(
            mutation_counts,
            labels=mutation_counts.index,
            autopct="%1.1f%%"
        )
        ax.set(
            title=f"{race}: Mutation"
        )

        # Pie chart for sex
        ax = axes[ax_row, 1]
        sex_counts = _meta.loc[:, "Sex"].fillna("Unknown").value_counts()
        ax.pie(
            sex_counts,
            labels=sex_counts.index,
            autopct="%1.1f%%"
        )
        ax.set(
            title=f"{race}: Sex"
        )

    return fig


if __name__ == "__main__":
    make_figure()
