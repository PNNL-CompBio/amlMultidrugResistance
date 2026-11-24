"""Plots figure 7c: Metabolic Dysregulation."""

from decimal import Decimal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from pilot.data_import import (import_acetyl, import_global, import_meta,
                               import_metabolites, syn_login)
from pilot.figures.figure_setup import get_setup, run_ols

COLORS = {
    "NDUFB6-K66k": "tab:green",
    "NDUFB6": "tab:blue"
}


def make_figure():
    # Import datasets
    syn = syn_login()
    acetyl = pd.concat(import_acetyl(syn))
    prot = pd.concat(import_global(syn))
    metab = pd.concat(import_metabolites(syn))

    # Trim to NDUFB6, patients with all meaasurements
    acetyl = acetyl.loc[
        acetyl.index.intersection(prot.index),
        "NDUFB6-K66k"
    ].dropna()
    acetyl = acetyl.loc[
        acetyl.index.intersection(metab.index)
    ]
    prot = prot.loc[
        acetyl.index,
        "NDUFB6"
    ]
    metab = metab.loc[
        acetyl.index,
        :
    ]

    # Derive correlations
    a_res = pearsonr(
        metab,
        pd.concat([acetyl] * metab.shape[1], axis=1)
    )
    p_res = pearsonr(
        metab,
        pd.concat([prot] * metab.shape[1], axis=1)
    )
    correlations = pd.DataFrame(
        [
            a_res.statistic,
            p_res.statistic
        ],
        index=["NDUFB6-K66k", "NDUFB6"],
        columns=metab.columns
    )
    p_values = pd.DataFrame(
        [
            a_res.pvalue,
            p_res.pvalue
        ],
        index=["NDUFB6-K66k", "NDUFB6"],
        columns=metab.columns
    )

    # Concatenate results and sort by correlation
    correlations = correlations.loc[
        :,
        abs(correlations).max(axis=0).sort_values().index[-10:]
    ]
    correlations = correlations.sort_values(
        by="NDUFB6-K66k",
        axis=1
    )

    # Setup figure
    fig, ax = get_setup(
        1,
        1,
        fig_params={
            "figsize": (3, 3)
        }
    )

    # Iterate through acetylation, proteomic correlations
    for x_offset, d_type in enumerate(correlations.index):
        # Plot correlations between metabolites, NDUFB6
        ax.bar(
            np.arange(x_offset, correlations.shape[1] * 3, 3),
            correlations.loc[d_type, :],
            color=COLORS[d_type],
            label=d_type
        )

    # Add dashed line at 0
    x_lims = ax.get_xlim()
    ax.plot(
        [-1, 100],
        [0, 0],
        linestyle="--",
        color="black"
    )

    # Format axes, add legend
    ax.set(
        xticks=np.arange(0.5, correlations.shape[1] * 3, 3),
        xlim=x_lims,
        ylabel="Pearson Correlation",
    )
    ax.set_xticklabels(
        correlations.columns.str[6:].str.replace(".", "-"),
        rotation=45,
        ha="right",
        ma="right",
        va="top"
    )
    ax.legend()

    return fig


if __name__ == "__main__":
    make_figure()
