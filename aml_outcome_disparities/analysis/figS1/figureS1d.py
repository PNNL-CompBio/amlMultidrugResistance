"""Plots figure S1d: NPM1/FLT3-ITD Justification"""

from os.path import abspath, dirname

import numpy as np
import pandas as pd

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
        ["NPM1", "FLT3_ITD", "DNMT3A", "NRAS"]
    ]
    white_mutations = meta.loc[
        meta.loc[:, "Race"] == "White",
        black_mutations.columns
    ]

    # Setup figure
    fig, axes = get_setup(
        1,
        3,
        fig_params={
            "figsize": (9, 3)
        }
    )

    # Plot counts in NPM1/FLT3-ITD
    ax = axes[0]

    # Iterate through Black, White patients
    for offset, (race, dataset) in enumerate(zip(
        ["Black", "White"],
        [black_mutations, white_mutations]
    )):
        # Get mutation counts
        # Only consider patients with the one mutation!
        _mutation_status = pd.Series(
            dtype=str,
            index=dataset.index
        )
        _mutation_status.loc[
            np.logical_and(
                dataset.loc[:, "FLT3_ITD"] == "WT",
                dataset.loc[:, "NPM1"] == "WT"
            )
        ] = "WT"
        _mutation_status.loc[
            np.logical_and(
                dataset.loc[:, "FLT3_ITD"] == "Mutant",
                dataset.loc[:, "NPM1"] == "WT"
            )
        ] = "FLT3-ITD"
        _mutation_status.loc[
            np.logical_and(
                dataset.loc[:, "FLT3_ITD"] == "WT",
                dataset.loc[:, "NPM1"] == "Mutant"
            )
        ] = "NPM1"
        _mutation_status.loc[
            np.logical_and(
                dataset.loc[:, "FLT3_ITD"] == "Mutant",
                dataset.loc[:, "NPM1"] == "Mutant"
            )
        ] = "NPM1/FLT3-ITD"

        # Drop patients with multiple mutations, reorder for consistency
        _mutation_status.dropna(inplace=True)
        _mutation_status = _mutation_status.value_counts()
        _mutation_status = _mutation_status.loc[
            [
                "WT",
                "FLT3-ITD",
                "NPM1",
                "NPM1/FLT3-ITD"
            ]
        ]

        # Plot bars, label quantities
        ax.bar(
            np.arange(offset, len(_mutation_status) * 3, 3),
            _mutation_status,
            color=RACE_COLORS[race],
            width=0.9,
            label=f"{race} Patients",
            zorder=3
        )
        for index, mutation_type in enumerate(_mutation_status.index):
            ax.text(
                3 * index + offset,
                _mutation_status.loc[mutation_type] + 0.1,
                ha="center",
                ma="center",
                va="bottom",
                s=_mutation_status.loc[mutation_type]
            )

    # Label axes and ticks, set legend and grid
    ax.set(
        xlim=(-1, 3 * 4 - 1),
        ylim=(0, 100),
        xticks=np.arange(
            0.5, 3 * 4, 3
        ),
        ylabel="Number of Patients\nwith Mutation"
    )
    ax.set_xticklabels(
        [
            "WT",
            "FLT3-ITD",
            "NPM1",
            "NPM1/\nFLT3-ITD"
        ],
        rotation=45,
        ha="right",
        ma="right",
        va="top"
    )
    ax.grid(True, zorder=0)
    ax.legend()

    # Compare single mutations when including two mutations with reasonable
    # quantities
    for ax, mutation in zip(axes[1:], ["DNMT3A", "NRAS"]):
        for offset, (race, dataset) in enumerate(zip(
            ["Black", "White"],
            [black_mutations, white_mutations]
        )):
            # Get counts of patients with only one mutation or NPM1/FLT3-ITD
            _mutation_status = pd.Series(
                dtype=str,
                index=dataset.index
            )
            _mutation_status.loc[
                np.logical_and(
                    dataset.loc[:, mutation] == "WT",
                    np.logical_and(
                        dataset.loc[:, "FLT3_ITD"] == "WT",
                        dataset.loc[:, "NPM1"] == "WT"
                    )
                )
            ] = "WT"
            _mutation_status.loc[
                np.logical_and(
                    dataset.loc[:, mutation] == "Mutant",
                    np.logical_and(
                        dataset.loc[:, "FLT3_ITD"] == "WT",
                        dataset.loc[:, "NPM1"] == "WT"
                    )
                )
            ] = mutation
            _mutation_status.loc[
                np.logical_and(
                    dataset.loc[:, mutation] == "WT",
                    np.logical_and(
                        dataset.loc[:, "FLT3_ITD"] == "Mutant",
                        dataset.loc[:, "NPM1"] == "WT"
                    )
                )
            ] = "FLT3-ITD"
            _mutation_status.loc[
                np.logical_and(
                    dataset.loc[:, mutation] == "WT",
                    np.logical_and(
                        dataset.loc[:, "FLT3_ITD"] == "WT",
                        dataset.loc[:, "NPM1"] == "Mutant"
                    )
                )
            ] = "NPM1"
            _mutation_status.loc[
                np.logical_and(
                    dataset.loc[:, mutation] == "WT",
                    np.logical_and(
                        dataset.loc[:, "FLT3_ITD"] == "Mutant",
                        dataset.loc[:, "NPM1"] == "Mutant"
                    )
                )
            ] = "NPM1/FLT3-ITD"

            # Drop patients with multiple mutations
            _mutation_status.dropna(inplace=True)
            _mutation_status = _mutation_status.value_counts()
            _mutation_status = _mutation_status.loc[
                [
                    "WT",
                    mutation,
                    "FLT3-ITD",
                    "NPM1",
                    "NPM1/FLT3-ITD"
                ]
            ]

            # Plot bars, label quantities
            ax.bar(
                np.arange(offset, len(_mutation_status) * 3, 3),
                _mutation_status,
                color=RACE_COLORS[race],
                width=0.9,
                label=f"{race} Patients",
                zorder=3
            )
            for index, mutation_type in enumerate(_mutation_status.index):
                ax.text(
                    3 * index + offset,
                    _mutation_status.loc[mutation_type] + 0.1,
                    ha="center",
                    ma="center",
                    va="bottom",
                    s=_mutation_status.loc[mutation_type]
                )

        # Label axes and ticks, set legend and grid
        ax.set(
            xlim=(-1, 3 * 5 - 1),
            ylim=(0, 100),
            xticks=np.arange(
                0.5, 5 * 3, 3
            ),
            ylabel="Number of Patients\nwith Mutation"
        )
        ax.set_xticklabels(
            [
                "WT",
                mutation,
                "FLT3-ITD",
                "NPM1",
                "NPM1/\nFLT3-ITD"
            ],
            rotation=45,
            ha="right",
            ma="right",
            va="top"
        )
        ax.grid(True, zorder=0)
        ax.legend()

    return fig


if __name__ == "__main__":
    make_figure()
