"""Plots figure 1c: Age Comparisons."""

from collections import OrderedDict

import numpy as np
from matplotlib.patches import Patch

from pilot.data_import import import_meta
from pilot.figures.figure_setup import get_setup


def make_figure():
    # Import metadata
    meta = import_meta()

    # Setup figure
    fig, ax = get_setup(
        1,
        1,
        fig_params={
            "figsize": (4, 4)
        }
    )

    # Setup groups for comparison
    mutation_groups = OrderedDict([
        (
            "All",
            meta.loc[meta.loc[:, "Race"].isin(["Black", "White"]), :].index
        ),
        (
            "WT",
            meta.loc[
                np.logical_and(
                    meta.loc[:, "NPM1"] != "Mutant",
                    meta.loc[:, "FLT3_ITD"] != "Mutant"
                ),
                :
            ].index
        ),
        (
            "NPM1",
            meta.loc[
                np.logical_and(
                    meta.loc[:, "NPM1"] == "Mutant",
                    meta.loc[:, "FLT3_ITD"] != "Mutant"
                ),
            :
            ].index
        ),
        (
            "FLT3-ITD",
            meta.loc[
                np.logical_and(
                    meta.loc[:, "NPM1"] != "Mutant",
                    meta.loc[:, "FLT3_ITD"] == "Mutant"
                ),
                :
            ].index
        ),
        (
            "NPM1/FLT3-ITD",
            meta.loc[
                np.logical_and(
                    meta.loc[:, "NPM1"] == "Mutant",
                    meta.loc[:, "FLT3_ITD"] == "Mutant"
                ),
                :
            ].index
        )
    ])

    # Iterate through mutation groups
    for offset, (group_name, patients) in enumerate(mutation_groups.items()):
        # Trim to patients matching mutation group
        _meta = meta.loc[patients, :]

        # Plot Black patients
        ax.barh(
            offset * 3,
            _meta.loc[_meta.loc[:, "Race"] == "Black", "Age"].mean(),
            xerr=_meta.loc[
                _meta.loc[:, "Race"] == "Black",
                "Age"
            ].std() * 1.96 / np.sqrt(sum(_meta.loc[:, "Race"] == "Black")),
            capsize=5,
            height=0.9,
            color="tab:red",
            zorder=3
        )
        # Plot White patients
        ax.barh(
            offset * 3 + 1,
            _meta.loc[_meta.loc[:, "Race"] == "White", "Age"].mean(),
            xerr=_meta.loc[
                _meta.loc[:, "Race"] == "White",
                "Age"
            ].std() * 1.96 / np.sqrt(sum(_meta.loc[:, "Race"] == "White")),
            capsize=5,
            height=0.9,
            color="tab:purple",
            zorder=3
        )

    # Add custom legend
    black_patients = Patch(color="tab:red", label="Black Patients")
    white_patients = Patch(color="tab:purple", label="White Patients")
    ax.legend(handles=[black_patients, white_patients])

    # Setup and label axes, add grid
    ax.set(
        xlim=(0, 80),
        xlabel="Age",
        yticks=np.arange(0.5, len(mutation_groups) * 3, 3),
        yticklabels=mutation_groups.keys()
    )
    ax.invert_yaxis()
    ax.grid(True, zorder=0)

    return fig


if __name__ == "__main__":
    make_figure()
