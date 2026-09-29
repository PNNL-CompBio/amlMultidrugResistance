"""Plots figure S1f: NPM1/FLT3 and Mutation Enrichment"""

from decimal import Decimal
from os.path import abspath, dirname

import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact

from pilot.data_import import import_meta
from pilot.figures.figure_setup import get_setup

FILE_DIR = dirname(abspath(__file__))
RACE_COLORS = {"Black": "tab:red", "White": "tab:purple"}
MUTATIONS = [
    "DNMT3A",
    "NRAS"
]


def make_figure():
    # Import meta data
    meta = import_meta()
    meta = meta.loc[meta.loc[:, "Race"].isin(["Black", "White"])]
    meta = meta.loc[
        ~meta.index.str.endswith("Bridge"),
        :
    ]

    # Convert to bool where True if mutation is present
    all_patients = meta.loc[
        :,
        "ASXL1":"ZRSR2"
    ] == "Mutant"

    # Split Black and White patients
    black_mutations = meta.loc[
        meta.loc[:, "Race"] == "Black",
        :
    ] == "Mutant"
    white_mutations = meta.loc[
        meta.loc[:, "Race"] == "White",
        :
    ] == "Mutant"

    # Trim to mutations present in at least 5 patients in either Race
    black_mutations = black_mutations.loc[
        :,
        np.logical_or(
            black_mutations.sum(axis=0) > 5,
            white_mutations.sum(axis=0) > 5
        )
    ]
    white_mutations = white_mutations.loc[:, black_mutations.columns]
    all_patients = all_patients.loc[:, black_mutations.columns]

    # Put mutation sets into an iterable list
    datasets = [black_mutations, white_mutations, all_patients]

    # Setup figure
    fig, axes = get_setup(
        len(MUTATIONS),
        2 * len(datasets),
        fig_params={
            "figsize": (4 * len(datasets), len(MUTATIONS) * 2),
        }
    )

    # Iterate through Black, White patients
    for race_index, (race, dataset) in enumerate(
        zip(
            ["Black", "White", "All"],
            datasets
        )
    ):
        # Iterate through FLT3, NPM1 mutations
        for mutation_index, base_gene in enumerate(["FLT3_ITD", "NPM1"]):
            # Iterate through comparison mutations
            for row_index, comp_gene in enumerate(MUTATIONS):
                # Get ax
                ax = axes[
                    row_index,
                    2 * race_index + mutation_index
                ]

                # Setup contingency table, fill values
                table = pd.DataFrame(
                    0,
                    dtype=int,
                    index=["WT", "Mutant"],  # Comparison gene
                    columns=["WT", "Mutant"]  # NPM1/FLT3-ITD
                )
                table.loc[
                    "WT",
                    "WT"
                ] = sum(
                    np.logical_and(
                        ~dataset.loc[:, base_gene],
                        ~dataset.loc[:, comp_gene]
                    )
                )
                table.loc[
                    "Mutant",
                    "WT"
                ] = sum(
                    np.logical_and(
                        ~dataset.loc[:, base_gene],
                        dataset.loc[:, comp_gene]
                    )
                )
                table.loc[
                    "WT",
                    "Mutant"
                ] = sum(
                    np.logical_and(
                        dataset.loc[:, base_gene],
                        ~dataset.loc[:, comp_gene]
                    )
                )
                table.loc[
                    "Mutant",
                    "Mutant"
                ] = sum(
                    np.logical_and(
                        dataset.loc[:, base_gene],
                        dataset.loc[:, comp_gene]
                    )
                )

                fet = fisher_exact(table)

                # Plot contingency table
                sns.heatmap(
                    table,
                    annot=True,
                    cmap="Reds",
                    fmt="d",
                    cbar=False,
                    annot_kws={"size": 30},
                    ax=ax
                )

                # Label axes and ticks
                ax.set(
                    xticks=np.arange(0.5, table.shape[0], 1),
                    yticks=np.arange(0.5, table.shape[1], 1),
                    xticklabels=table.columns,
                    yticklabels=table.index,
                    xlabel=base_gene,
                    ylabel=comp_gene
                )

                if row_index == 0:
                    ax.set_title(f"{race}: {base_gene}")

                # Include Fisher's Exact Test result
                ax.text(
                    0.99,
                    0.01,
                    s=f"Fisher's Exact: {round(fet.statistic, 2)}\n"
                      f"p-value: {'{:.2E}'.format(Decimal(fet.pvalue))}",
                    transform=ax.transAxes,
                    ha="right",
                    ma="right",
                    va="bottom",
                    fontsize=6,
                    color="black",
                )

    return fig


if __name__ == "__main__":
    make_figure()
