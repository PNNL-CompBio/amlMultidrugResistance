"""Plots figure 1c: Age Comparisons."""

import numpy as np

from pilot.data_import import import_meta
from pilot.figures.figure_setup import get_setup


def make_figure():
    # Import metadata
    meta = import_meta()

    # Setup figure
    fig, axes = get_setup(
        2,
        1,
        fig_params={
            "figsize": (4, 4),
            "sharex": True
        }
    )

    # Iterate through race, axes
    for ax, race in zip(axes, ["Black", "White"]):
        # Trim to patients with matching race
        _meta = meta.loc[meta.loc[:, "Race"] == race, :]

        # Plot all patients with matching race
        ax.barh(
            0,
            _meta.loc[:, "Age"].mean(),
            xerr=_meta.loc[:, "Age"].std(),
            capsize=5,
            height=0.9
        )

        # Plot non-mutant patients
        mutation_meta = _meta.loc[
            np.logical_and(
                _meta.loc[:, "NPM1"] != "Mutant",
                _meta.loc[:, "FLT3_ITD"] != "Mutant"
            ),
            "Age"
        ]
        ax.barh(
            1,
            mutation_meta.mean(),
            xerr=mutation_meta.std(),
            capsize=5,
            height=0.9
        )

        # Plot NPM1 mutant patients
        mutation_meta = _meta.loc[
            np.logical_and(
                _meta.loc[:, "NPM1"] == "Mutant",
                _meta.loc[:, "FLT3_ITD"] != "Mutant"
            ),
            "Age"
        ]
        ax.barh(
            2,
            mutation_meta.mean(),
            xerr=mutation_meta.std(),
            capsize=5,
            height=0.9
        )

        # Plot FLT3-ITD patients
        mutation_meta = _meta.loc[
            np.logical_and(
                _meta.loc[:, "NPM1"] != "Mutant",
                _meta.loc[:, "FLT3_ITD"] == "Mutant"
            ),
            "Age"
        ]
        ax.barh(
            3,
            mutation_meta.mean(),
            xerr=mutation_meta.std(),
            capsize=5,
            height=0.9
        )

        # Plot patients with both NPM1 and FLT3-ITD mutations
        mutation_meta = _meta.loc[
            np.logical_and(
                _meta.loc[:, "NPM1"] == "Mutant",
                _meta.loc[:, "FLT3_ITD"] == "Mutant"
            ),
            "Age"
        ]
        ax.barh(
            4,
            mutation_meta.mean(),
            xerr=mutation_meta.std(),
            capsize=5,
            height=0.9
        )

        # Invert axes, set labels
        ax.invert_yaxis()
        ax.set(
            xticks=np.arange(0, 90, 10),
            yticks=np.arange(5),
            yticklabels=[
                f"{race} Patients",
                f"{race}: Wild Type",
                f"{race}: NPM1",
                f"{race}: FLT3-ITD",
                f"{race}: Both"
            ]
        )

    # Label bottom x-axis with age
    axes[1].set_xlabel("Age")

    return fig


if __name__ == "__main__":
    make_figure()
