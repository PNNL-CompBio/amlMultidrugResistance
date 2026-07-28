"""Plots figure 3d: S100A8 Acetylation and Inflammation, nF-KB Signaling."""

from decimal import Decimal

import pandas as pd

from pilot.data_import import import_acetyl, import_global, syn_login
from pilot.figures.figure_setup import get_setup, run_ols

ACETYL_GENE = "S100A8"
COMPARISON_GENES = [
    "GRN",
    "PYCARD",
    "GSK3B",
    "BCL10",
    "ILKAP"
]


def make_figure():
    # Import data
    syn = syn_login()
    acetyl = pd.concat(import_acetyl(syn))
    prot = pd.concat(import_global(syn))

    # Trim to patients with both acetylation and proteomic measurements, take
    # mean of acetylation across S100A8 sites
    prot = prot.loc[
        prot.index.intersection(acetyl.index),
        :
    ]
    acetyl = acetyl.loc[
        prot.index,
        acetyl.columns.str.startswith(ACETYL_GENE)
    ].mean(axis=1)
    acetyl = acetyl.dropna()
    prot = prot.loc[acetyl.index, :]

    # Setup figure
    fig, axes = get_setup(
        2,
        len(COMPARISON_GENES),
        fig_params={
            "figsize": (3 * len(COMPARISON_GENES), 6)
        }
    )

    # Iterate through genes to compare against acetylation, proteomic
    for ax_col, comparison_gene in enumerate(COMPARISON_GENES):
        # Compare S100A8 with proteomics
        ax = axes[0, ax_col]
        ax.scatter(
            prot.loc[:, ACETYL_GENE],
            prot.loc[:, comparison_gene],
            s=3
        )
        ols = run_ols(
            prot.loc[:, ACETYL_GENE],
            prot.loc[:, comparison_gene],
            ax=ax
        )
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
        ax.set(
            xlabel=f"Global: {ACETYL_GENE}",
            ylabel=comparison_gene
        )

        # Compare S100A8 acetylation with proteomics
        ax = axes[1, ax_col]
        ax.scatter(
            acetyl,
            prot.loc[:, comparison_gene],
            s=3
        )
        ols = run_ols(
            acetyl,
            prot.loc[:, comparison_gene],
            ax=ax
        )
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
        ax.set(
            xlabel=f"{ACETYL_GENE} Acetylation",
            ylabel=comparison_gene
        )

    return fig


if __name__ == "__main__":
    make_figure()
