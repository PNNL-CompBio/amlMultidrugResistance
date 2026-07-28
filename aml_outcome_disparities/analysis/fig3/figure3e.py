"""Plots figure 3e: S100A8 Acetylation & MYC."""

from decimal import Decimal

import gseapy as gp
import pandas as pd

from pilot.data_import import import_acetyl, import_global, syn_login
from pilot.figures.figure_setup import get_setup, run_ols

ACETYL_GENE = "S100A8"


def make_figure():
    # Import data
    syn = syn_login()
    acetyl = pd.concat(import_acetyl(syn))
    prot = pd.concat(import_global(syn))

    # Load MYC targets
    myc_targets = gp.get_library("MSigDB_Hallmark_2020")["Myc Targets V2"]
    myc_targets = prot.loc[
        :,
        prot.columns.intersection(myc_targets)
    ].columns
    myc_mean = prot.loc[:, myc_targets].mean(axis=1)

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
    myc_mean = myc_mean.loc[acetyl.index]

    # Setup figure
    fig, axes = get_setup(
        1,
        2,
        fig_params={
            "figsize": (6, 3)
        }
    )

    # Compare S100A8 with proteomics
    ax = axes[0]
    ax.scatter(
        prot.loc[:, ACETYL_GENE],
        myc_mean,
        s=3
    )
    ols = run_ols(
        prot.loc[:, ACETYL_GENE],
        myc_mean,
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
        ylabel="MYC Targets"
    )

    # Compare S100A8 acetylation with proteomics
    ax = axes[1]
    ax.scatter(
        acetyl,
        myc_mean,
        s=3
    )
    ols = run_ols(
        acetyl,
        myc_mean,
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
        ylabel="MYC Targets"
    )

    return fig


if __name__ == "__main__":
    make_figure()
