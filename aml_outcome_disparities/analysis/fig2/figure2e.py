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
        metabolites.index.intersection(prkaca_scores.index), :
    ]
    lipids = lipids.loc[lipids.index.intersection(prkaca_scores.index), :]

    # Setup dataframes to store results
    lipid_stats = pd.DataFrame(
        0, dtype=float, index=lipids.columns, columns=["Rho", "FDR q-value"]
    )
    meta_stats = pd.DataFrame(
        0,
        dtype=float,
        index=metabolites.columns,
        columns=["Rho", "FDR q-value"],
    )

    # Compare PRKACA and metabolites, correct p-values via FDR
    meta_res = pearsonr(
        pd.concat(
            [prkaca_scores.loc[metabolites.index]] * metabolites.shape[1],
            axis=1,
        ),
        np.log2(metabolites),
    )
    meta_stats.loc[:, "Rho"] = meta_res.statistic
    meta_stats.loc[:, "FDR q-value"] = false_discovery_control(meta_res.pvalue)

    # Compare PRKACA and lipids
    # Needs to be done in for loop to account for missing lipid measurements
    for lipid in lipids.columns:
        measurement = lipids.loc[:, lipid].dropna()
        lipid_res = pearsonr(measurement, prkaca_scores.loc[measurement.index])
        lipid_stats.loc[lipid, "Rho"] = lipid_res.statistic
        lipid_stats.loc[lipid, "FDR q-value"] = lipid_res.pvalue

    # Correct p-values via FDR
    lipid_stats.loc[:, "FDR q-value"] = false_discovery_control(
        lipid_stats.loc[:, "FDR q-value"]
    )

    # Reformat lipid, metabolite names
    lipid_stats.index = lipid_stats.index.str[4:]
    lipid_stats.index = [
        lipid_name.split("__")[0] for lipid_name in lipid_stats.index
    ]
    lipid_stats.index = [
        lipid_name.split("|")[1] if ("|" in lipid_name) else lipid_name for
        lipid_name in lipid_stats.index
    ]
    meta_stats.index = meta_stats.index.str[6:]

    # Setup plot
    fig, axes = get_setup(
        1,
        2,
        fig_params={"figsize": (8, 3)},
    )

    # Iterate through 'omes
    for name, stat, ax in zip(
        ["Lipids", "Metabolites"], [lipid_stats, meta_stats], axes
    ):
        # Trim to top 30
        stat = stat.sort_values(by="FDR q-value", ascending=True)
        stat = stat.iloc[:25, :].sort_values(by="Rho", ascending=True)

        # Color measurements with significant differences in expression
        colors = pd.Series(["grey"] * stat.shape[0], index=stat.index)
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
        ax.bar(
            np.arange(stat.shape[0]),
            stat.loc[:, "Rho"],
            color=colors
        )

        # Plot dashed line at 0
        x_lims, y_lims = ax.get_xlim(), ax.get_ylim()
        ax.plot(
            [-2, 100],
            [0, 0],
            color="black",
            linestyle="--",
        )

        # Format plot, enable grid, label axes
        ax.grid(True)
        ax.set_axisbelow(True)
        ax.set(
            xlim=x_lims,
            ylim=y_lims,
            xticks=np.arange(stat.shape[0]),
            ylabel="Pearson Correlation"
        )
        ax.set_xticklabels(
            stat.index,
            rotation=45,
            ha="right",
            ma="right",
            va="top",
            fontsize=5
        )

    return fig


if __name__ == "__main__":
    make_figure()
