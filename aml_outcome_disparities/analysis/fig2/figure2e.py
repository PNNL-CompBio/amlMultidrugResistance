import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import false_discovery_control, pearsonr
import pandas as pd

from pilot.data_import import syn_login, import_lipids, import_metabolites
from pilot.figures.figure_setup import get_setup
from pilot.gene_analysis import get_prkaca_score


def make_figure():
    # Get PRKACA scores
    prkaca_scores = get_prkaca_score()

    # Load metabolites, lipids
    syn = syn_login()
    metabolites = import_metabolites(syn)
    lipids = import_lipids(syn)
    metabolites = pd.concat(metabolites, axis=0)
    lipids = pd.concat(lipids, axis=0)

    # Trim to patients with PRKACA scores
    metabolites = metabolites.loc[
        metabolites.index.intersection(prkaca_scores.index),
        :
    ]
    lipids = lipids.loc[
        lipids.index.intersection(prkaca_scores.index),
        :
    ]

    # Setup dataframes to store results
    lipid_stats = pd.DataFrame(
        0,
        dtype=float,
        index=lipids.columns,
        columns=["Rho", "FDR q-value"]
    )
    meta_stats = pd.DataFrame(
        0,
        dtype=float,
        index=metabolites.columns,
        columns=["Rho", "FDR q-value"]
    )

    # Compare PRKACA and metabolites, correct p-values via FDR
    meta_res = pearsonr(
        pd.concat(
            [prkaca_scores.loc[metabolites.index]] * metabolites.shape[1],
            axis=1
        ),
        np.log2(metabolites)
    )
    meta_stats.loc[:, "Rho"] = meta_res.statistic
    meta_stats.loc[:, "FDR q-value"] = false_discovery_control(
        meta_res.pvalue
    )

    # Compare PRKACA and lipids
    # Needs to be done in for loop to account for missing lipid measurements
    for lipid in lipids.columns:
        measurement = lipids.loc[:, lipid].dropna()
        lipid_res = pearsonr(
            measurement,
            prkaca_scores.loc[measurement.index]
        )
        lipid_stats.loc[lipid, "Rho"] = lipid_res.statistic
        lipid_stats.loc[lipid, "FDR q-value"] = lipid_res.pvalue

    # Correct p-values via FDR
    lipid_stats.loc[:, "FDR q-value"] = false_discovery_control(
        lipid_stats.loc[:, "FDR q-value"]
    )

    # Setup plot
    fig, axes = get_setup(
        1,
        2,
        fig_params={"figsize": (6, 3)},
    )

    # Iterate through 'omes
    for name, stat, ax in zip(
        ["Lipids", "Metabolites"],
        [lipid_stats, meta_stats],
        axes
    ):
        # Color measurements with significant differences in expression
        colors = pd.Series(
            ["grey"] * stat.shape[0],
            index=stat.index
        )
        colors.loc[
            np.logical_and(
                stat.loc[:, "Rho"] > 0,
                stat.loc[:, "FDR q-value"] < 0.05,
            )
        ] = "tab:red"
        colors.loc[
            np.logical_and(
                stat.loc[:, "Rho"] < 0,
                stat.loc[:, "FDR q-value"] < 0.05,
            )
        ] = "tab:blue"

        # Plot volcano
        ax.scatter(
            stat.loc[:, "Rho"],
            -np.log10(stat.loc[:, "FDR q-value"]),
            c=colors,
            alpha=0.5
        )
        x_lims, y_lims = ax.get_xlim(), ax.get_ylim()

        # Add line for q-value threshold
        ax.plot(
            [-1, 1],
            [-np.log10(0.05), -np.log10(0.05)],
            linestyle="--",
            color="black",
            zorder=1,
        )
        ax.plot(
            [0, 0],
            [-1, 10],
            linestyle="--",
            color="black",
            zorder=1,
        )

        # Format plot, enable grid, label axes
        ax.grid(True)
        ax.set_axisbelow(True)
        ax.set(
            xlim=x_lims,
            ylim=y_lims,
            xlabel="Pearson Rho",
            ylabel="-log10(FDR q-value)",
            title=name
        )

    return fig


if __name__ == "__main__":
    make_figure()
