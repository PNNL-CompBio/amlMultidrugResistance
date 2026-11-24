"""Plots figure 3b: S100A8 Acetylation & ILKAP."""

from decimal import Decimal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pilot.data_import import (import_acetyl, import_global, import_rna,
                               syn_login)
from pilot.figures.figure_setup import get_setup, run_ols

GENE = "ILKAP"


def make_figure():
    syn = syn_login()
    acetyl = pd.concat(import_acetyl(syn))
    prot = pd.concat(import_global(syn))
    rna = np.log2(pd.concat(import_rna(syn)) + 1)

    prot = prot.loc[
        prot.index.intersection(acetyl.index),
        :
    ]
    prot = prot.loc[
        prot.index.intersection(rna.index),
        :
    ]
    rna = rna.loc[prot.index, :]
    acetyl = acetyl.loc[
        prot.index,
        acetyl.columns.str.startswith("S100A8")
    ]

    fig, axes = get_setup(
        1,
        2,
        fig_params={
            "figsize": (3.5, 3),
            "width_ratios": [3.375, 0.125]
        }
    )

    ax = axes[0]
    colorbar_ax = axes[1]

    scatter = ax.scatter(
        rna.loc[:, GENE],
        prot.loc[:, GENE],
        s=3,
        c=acetyl.loc[
            :,
            acetyl.columns.str.startswith("S100A8")
        ].mean(axis=1),
        cmap="coolwarm"
    )
    ols = run_ols(rna.loc[:, GENE], prot.loc[:, GENE], ax)

    ax.text(
        0.99,
        0.01,
        ha="right",
        ma="right",
        va="bottom",
        s=f"R-squared: {round(ols.rsquared_adj, 3)}\nCoefficient p-value: "
          f"{'{:.2E}'.format(Decimal(ols.pvalues.iloc[0]))}",
        transform=ax.transAxes
    )

    plt.colorbar(
        scatter,
        label="S100A8 Acetylation",
        cax=colorbar_ax
    )
    ax.set(
        xlabel=f"RNA: {GENE}",
        ylabel=f"Protein: {GENE}"
    )

    return fig


if __name__ == "__main__":
    make_figure()
