from decimal import Decimal

from scipy.stats import false_discovery_control, pearsonr, fisher_exact
import numpy as np
import pandas as pd
import seaborn as sns

from pilot.data_import import syn_login, import_metabolites
from pilot.figures.figure_setup import get_setup
from pilot.gene_analysis import get_prkaca_score


def make_figure():
    # Get PRKACA scores
    prkaca_scores = get_prkaca_score()

    # Load metabolites, lipids
    syn = syn_login()
    metabolites = import_metabolites(syn)
    metabolites = pd.concat(metabolites, axis=0)

    # Trim to patients with PRKACA scores
    metabolites = metabolites.loc[
        metabolites.index.intersection(prkaca_scores.index), :
    ]
    meta_stats = pd.DataFrame(
        0,
        dtype=float,
        index=metabolites.columns,
        columns=["Rho", "FDR q-value"],
    )

    # Pearson correlations between PRKACA and metabolites
    meta_res = pearsonr(
        pd.concat(
            [prkaca_scores.loc[metabolites.index]] * metabolites.shape[1],
            axis=1,
        ),
        np.log2(metabolites),
    )
    meta_stats.loc[:, "Rho"] = meta_res.statistic
    meta_stats.loc[:, "FDR q-value"] = false_discovery_control(meta_res.pvalue)
    meta_stats.index = meta_stats.index.str.lower()

    # Compare if carnitine is specifically enriched
    # Setup fisher's exact test table
    carnitine_table = pd.DataFrame(
        0,
        dtype=int,
        index=["Other", "Carnitine"],
        columns=["Non-Significant", "Significant"],
    )
    sig_meta = meta_stats.loc[
        np.logical_and(
            meta_stats.loc[:, "FDR q-value"] < 0.05,
            meta_stats.loc[:, "Rho"] > 0,
        ),
        :,
    ]
    ns_meta = meta_stats.drop(sig_meta.index)

    # Fill in contingency table
    carnitine_table.loc["Carnitine", "Significant"] = sum(
        sig_meta.index.str.lower().str.contains("carnitine")
    )
    carnitine_table.loc["Carnitine", "Non-Significant"] = sum(
        ns_meta.index.str.lower().str.contains("carnitine")
    )
    carnitine_table.loc["Other", "Significant"] = sum(
        ~sig_meta.index.str.lower().str.contains("carnitine")
    )
    carnitine_table.loc["Other", "Non-Significant"] = sum(
        ~ns_meta.index.str.lower().str.contains("carnitine")
    )

    # Rename columns for plotting
    carnitine_table.columns = [
        "Not Significantly\nPositively Correlated",
        "Significantly\nPositively Correlated",
    ]

    # Run Fisher's Exact
    fisher_res = fisher_exact(carnitine_table)

    # Plot contingency table
    fig, ax = get_setup(1, 1, {"figsize": (3, 3)})
    sns.heatmap(
        carnitine_table,
        annot=True,
        cmap="Reds",
        cbar=False,
        annot_kws={"size": 30},
    )

    # Label axes and ticks
    ax.set(
        xticks=np.arange(0.5, carnitine_table.shape[0], 1),
        yticks=np.arange(0.5, carnitine_table.shape[1], 1),
        xticklabels=carnitine_table.columns,
        yticklabels=carnitine_table.index,
    )

    # Include Fisher's Exact Test result
    ax.text(
        0.99,
        0.99,
        s=f"Fisher's Exact: {round(fisher_res.statistic, 2)}\n"
        f"p-value: {'{:.2E}'.format(Decimal(fisher_res.pvalue))}",
        transform=ax.transAxes,
        ha="right",
        ma="right",
        va="top",
        fontsize=6,
        color="black",
    )

    return fig


if __name__ == "__main__":
    make_figure()
