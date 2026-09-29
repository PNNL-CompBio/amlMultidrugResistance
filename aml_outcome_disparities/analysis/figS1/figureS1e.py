"""Plots figure S1e: Mutation Comparison Across Race"""

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
    "NRAS",
    "NPM1",
    "FLT3_ITD"
]


def make_figure():
    # Import meta data
    meta = import_meta()
    meta = meta.loc[meta.loc[:, "Race"].isin(["Black", "White"])]
    meta = meta.loc[
        ~meta.index.str.endswith("Bridge"),
        :
    ]

    # Setup figure
    fig, axes = get_setup(
        1,
        len(MUTATIONS),
        fig_params={
            "figsize": (2 * len(MUTATIONS), 2),
        }
    )

    # Iterate through mutations
    for ax, gene in zip(axes, MUTATIONS):
        # Setup contingency table, fill values
        table = pd.DataFrame(
            0,
            dtype=int,
            index=["Black", "White"],
            columns=["WT", "Mutant"]
        )
        table.loc[
            "Black",
            "WT"
        ] = sum(
            np.logical_and(
                meta.loc[:, "Race"] == "Black",
                meta.loc[:, gene] != "Mutant"
            )
        )
        table.loc[
            "White",
            "WT"
        ] = sum(
            np.logical_and(
                meta.loc[:, "Race"] == "White",
                meta.loc[:, gene] != "Mutant"
            )
        )
        table.loc[
            "Black",
            "Mutant"
        ] = sum(
            np.logical_and(
                meta.loc[:, "Race"] == "Black",
                meta.loc[:, gene] == "Mutant"
            )
        )
        table.loc[
            "White",
            "Mutant"
        ] = sum(
            np.logical_and(
                meta.loc[:, "Race"] == "White",
                meta.loc[:, gene] == "Mutant"
            )
        )

        # Run Fisher's Exact
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
            ylabel="Race",
            title=gene
        )

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
