"""Plots figure 5d: AURKB and CDK1"""

import sys
from decimal import Decimal
from os.path import abspath, dirname, join

import numpy as np
from scipy.stats import pearsonr

from pilot.data_import import import_meta
from pilot.figures.figure_setup import get_setup, run_ols
from pilot.gene_analysis import get_kinase_score

sys.path.append(
    join(dirname(dirname(dirname(abspath(__file__)))), "src", "python")
)


def make_figure():
    # Import meta-data, reformat mutations
    meta = import_meta()
    meta.loc[:, "ALT":] = meta.loc[:, "ALT":].replace(
        {
            np.nan: "WT",
            "Not measured": "WT"
        }
    )
    meta = meta.loc[:, ["Race", "NPM1"]]

    # Score kinase activity
    aurkb = get_kinase_score("AURKB")
    cdk1 = get_kinase_score("CDK1", force_mean=True)

    # Trim to patients with both kinase measurements
    cdk1 = cdk1.loc[cdk1.index.intersection(aurkb.index)]
    cdk1 = cdk1.loc[cdk1.index.intersection(meta.index)]
    aurkb = aurkb.loc[cdk1.index]
    meta = meta.loc[cdk1.index, :]

    # Setup plot
    fig, axes = get_setup(2, 2, fig_params={"figsize": (4, 4)})

    # Iterate through race, NPM1 mutations
    for ax_row, race in enumerate(["Black", "White"]):
        for ax_col, mutation in enumerate(["WT", "Mutant"]):
            ax = axes[ax_row, ax_col]

            # Subset by race, mutation
            _aurkb = aurkb.loc[
                np.logical_and(
                    meta.loc[:, "Race"] == race,
                    meta.loc[:, "NPM1"] == mutation
                )
            ]
            _cdk1 = cdk1.loc[
                np.logical_and(
                    meta.loc[:, "Race"] == race,
                    meta.loc[:, "NPM1"] == mutation
                )
            ]

            # Get correlations, compare kinase scores
            res = pearsonr(_aurkb, _cdk1)
            ax.scatter(
                _aurkb,
                _cdk1,
                s=3
            )
            run_ols(_aurkb, _cdk1, ax)

            # Denote results in plot
            ax.text(
                0.99,
                0.01,
                s=f"p-value: {'{:.2E}'.format(Decimal(res.pvalue))}\n"
                  f"Corr: {round(res.statistic, 3)}",
                ha="right",
                ma="right",
                va="bottom",
                fontsize=8,
                transform=ax.transAxes,
            )
            ax.set(
                xlabel="AURKB",
                ylabel="CDK1",
                title=f"{race}: NPM1 {mutation}"
            )

    return fig


if __name__ == "__main__":
    make_figure()
