"""Plots figure S1c: Myelodysplasia investigation."""

import sys
from decimal import Decimal
from os.path import abspath, dirname, join

sys.path.append(
    join(dirname(dirname(dirname(abspath(__file__)))), "src", "python")
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact

from pilot.data_import import import_meta
from pilot.figures.figure_setup import get_setup

AGE_THRESHOLD = 40
REPO_PATH = abspath(dirname(dirname(dirname(__file__))))
MYELODYSPLASIA_GENES = [
    "ASXL1", "BCOR", "EZH2", "RUNX1", "SF3B1", "SRSF2", "STAG2", "U2AF1",
    "ZRSR2"
]



def make_figure():
    # Import meta-data
    meta = import_meta()

    # Trim to Black/White patients
    meta = meta.loc[
        meta.loc[:, "Race"].isin(["Black", "White"]),
        :
    ]

    # Trim to young AML patients
    meta = meta.loc[
        meta.loc[:, "Age"].between(15, AGE_THRESHOLD, inclusive="left"),
        :
    ]
    meta = meta.loc[~meta.index.str.endswith("Bridge"), :]

    # Identify MDS-like AML cases
    # Defined as presence of any MDS-like mutation WITHOUT mutations in NPM1,
    # FLT3-ITD, or TP53
    meta.loc[:, "MDS-like"] = "WT"
    meta.loc[np.logical_and(
        ~np.any(
            meta.loc[
                :,
                [
                    "NPM1",
                    "FLT3_ITD",
                    "TP53"
                ]
            ] == "Mutant",
            axis=1
        ),
        np.any(
            meta.loc[:, MYELODYSPLASIA_GENES] == "Mutant",
            axis=1
        )
    ), "MDS-like"] = "Mutant"

    # Remove double NPM1/FLT3-ITD cases from NPM1/FLT3-ITD calls
    meta.loc[
        np.logical_and(
            meta.loc[:, "FLT3_ITD"] == "Mutant",
            meta.loc[:, "NPM1"] == "Mutant"
        ),
        ["NPM1", "FLT3_ITD"]
    ] = "WT"

    # Setup plots
    fig, axes = get_setup(1, 3, {"figsize": (9, 3)})

    for mutation, ax in zip(["NPM1", "FLT3_ITD", "MDS-like"], axes):
        # Build table for FET
        table = pd.DataFrame(
            0,
            dtype=int,
            index=[f"Non-{mutation}", mutation],
            columns=["White", "Black"]
        )
        table.loc[f"Non-{mutation}", "White"] = sum(
            np.logical_and(
                meta.loc[:, "Race"] == "White",
                meta.loc[:, mutation] != "Mutant"
            )
        )
        table.loc[mutation, "White"] = sum(
            np.logical_and(
                meta.loc[:, "Race"] == "White",
                meta.loc[:, mutation] == "Mutant"
            )
        )
        table.loc[f"Non-{mutation}", "Black"] = sum(
            np.logical_and(
                meta.loc[:, "Race"] == "Black",
                meta.loc[:, mutation] != "Mutant"
            )
        )
        table.loc[mutation, "Black"] = sum(
            np.logical_and(
                meta.loc[:, "Race"] == "Black",
                meta.loc[:, mutation] == "Mutant"
            )
        )
    
        # Run Fisher's Exact
        fisher_res = fisher_exact(table)
        
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
        )

        # Include Fisher's Exact Test result
        ax.text(
            0.99,
            0.01,
            s=f"Fisher's Exact: {round(fisher_res.statistic, 2)}\n"
            f"p-value: {'{:.2E}'.format(Decimal(fisher_res.pvalue))}",
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
