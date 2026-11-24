"""Plots figure 7d: H3 and Mitochondrial Metabolism."""

from decimal import Decimal

import pandas as pd

from pilot.data_import import import_acetyl, import_global, syn_login
from pilot.figures.figure_setup import get_setup, run_ols

MITO_GENE = "NDUFB6"
HIST_GENE = "HIST1H3A"


def make_figure():
    # Import meta-data, data
    syn = syn_login()
    prot = pd.concat(import_global(syn))
    acetyl = pd.concat(import_acetyl(syn))

    # Trim to NDUFB6
    mito_acetyl = acetyl.loc[
        :,
        acetyl.columns.str.startswith(MITO_GENE)
    ].dropna()
    mito_prot = prot.loc[
        :,
        MITO_GENE
    ].dropna()
    hist_acetyl = acetyl.loc[
        :,
        acetyl.columns.str.startswith(HIST_GENE)
    ].mean(axis=1).dropna()

    # Trim to patients with no missingness in NDUFB6 across measurements
    mito_acetyl = mito_acetyl.loc[
        mito_acetyl.index.intersection(hist_acetyl.index),
    ]
    mito_acetyl = mito_acetyl.loc[
        mito_acetyl.index.intersection(mito_prot.index)
    ]
    mito_prot = mito_prot.loc[mito_acetyl.index]
    hist_acetyl = hist_acetyl.loc[mito_acetyl.index]

    # Setup figure
    fig, axes = get_setup(
        2,
        1,
        fig_params={
            "figsize": (3, 6)
        }
    )

    for ax, name, mito_exp in zip(
        axes,
        ["NDUFB6-K66k", "NDUFB6"],
        [mito_acetyl, mito_prot]
    ):
        # Plot NDUFB6 vs. H3 acetylation
        ax.scatter(
            mito_exp,
            hist_acetyl,
            s=3
        )

        # Add best fit line
        ols = run_ols(
            mito_exp,
            hist_acetyl,
            ax
        )

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
            xlabel=name,
            ylabel=f"{HIST_GENE} Acetylation"
        )

    return fig


if __name__ == "__main__":
    make_figure()
