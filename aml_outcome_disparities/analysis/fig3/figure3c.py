"""Plots figure 3c: S100A8 Acetylation vs. ILKAP Gap."""

from decimal import Decimal

import numpy as np
import pandas as pd

from pilot.data_import import (import_acetyl, import_global, import_rna,
                               syn_login)
from pilot.figures.figure_setup import get_setup, run_ols

GENE = "ILKAP"


def make_figure():
    # Import meta-data, data, log-transform TPM
    syn = syn_login()
    acetyl = pd.concat(import_acetyl(syn, global_correct=False))
    prot = pd.concat(import_global(syn))
    rna = np.log2(pd.concat(import_rna(syn)) + 1)

    # Trim to patients with all measurements
    prot = prot.loc[
        prot.index.intersection(acetyl.index),
        :
    ]
    prot = prot.loc[prot.index.intersection(rna.index), :]
    rna = rna.loc[prot.index, GENE]
    diff = prot.loc[:, GENE] - rna

    # Trim to S100A8 sites
    acetyl = acetyl.loc[
        prot.index,
        acetyl.columns.str.startswith("S100A8")
    ].mean(axis=1).dropna()

    # Setup figure
    fig, ax = get_setup(
        1,
        1,
        fig_params={
            "figsize": (3, 3)
        }
    )

    # Plot acetylation vs. global/RNA differences
    ax.scatter(
        acetyl,
        diff.loc[acetyl.index],
        s=3
    )

    # Add best fit line
    ols = run_ols(acetyl, diff.loc[acetyl.index], ax)

    # Write OLS results
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

    # Label axes
    ax.set(
        xlabel=f"S100A8 Acetylation",
        ylabel=f"{GENE}: Global - RNA"
    )

    return fig


if __name__ == "__main__":
    make_figure()
