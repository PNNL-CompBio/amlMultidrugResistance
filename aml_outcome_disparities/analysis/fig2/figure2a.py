"""Plots figure 2a: KSEA kinase selection."""

from decimal import Decimal
import sys

from matplotlib.patches import Patch
from os.path import abspath, dirname, join
from sklearn.preprocessing import LabelEncoder
from statsmodels.multivariate.manova import MANOVA
import numpy as np
import pandas as pd

from pilot.data_import import import_meta, import_phospho, syn_login
from pilot.figures.figure_setup import get_setup, compare_means

sys.path.append(
    join(dirname(dirname(dirname(abspath(__file__)))), "src", "python")
)

ANALYSIS_DIR = dirname(abspath(__file__))


def make_figure():
    # Load KSEA results
    ksea_table = pd.read_csv(join(ANALYSIS_DIR, "ksea_results.csv"))
    substrates = pd.read_csv(join(ANALYSIS_DIR, "Kinase-Substrate Links.csv"))

    # Trim to significantly enriched Kinases with >10 substrates
    ksea_table = ksea_table.loc[
        np.logical_and(
            ksea_table.loc[:, "m"] > 10,
            ksea_table.loc[:, "FDR"] < 0.05
        )
    ]
    substrates = substrates.loc[
        substrates.loc[:, "Kinase.Gene"].isin(ksea_table["Kinase.Gene"]),
        :
    ]

    # Reformat names to match Synapse data, get counts of each Kinase's
    # substrates
    substrates.index = (
        substrates.loc[:, "Substrate.Gene"]
        + "-"
        + substrates.loc[:, "Substrate.Mod"]
        + substrates.loc[:, "Substrate.Mod"].str[0].str.lower()
    )
    kinase_counts = substrates.loc[:, "Kinase.Gene"].value_counts()

    # Load data
    syn = syn_login()
    phospho = import_phospho(syn)
    phospho = pd.concat(phospho, axis=0)

    # Trim KSEA table to phospho substrates that also appear in phospho
    # measurements, trim to substrates with >90% coverage across patients
    shared = substrates.index.intersection(phospho.columns)
    phospho = phospho.loc[:, shared]
    phospho = phospho.loc[
        :,
        phospho.isna().mean(axis=0) < 0.01
    ]

    # Get counts of high-coverage substrates
    hc_counts = substrates.loc[
        phospho.columns,
        "Kinase.Gene"
    ].value_counts()
    hc_counts = hc_counts.loc[kinase_counts.index]

    fig, ax = get_setup(
        1,
        1,
        fig_params={"figsize": (6, 4)},
    )

    ax.barh(
        np.arange(0, 3 * len(kinase_counts), 3),
        kinase_counts,
        height=1,
        color="tab:blue",
        label="Substrates"
    )
    ax.barh(
        np.arange(1, 3 * len(kinase_counts), 3),
        hc_counts,
        height=1,
        color="tab:orange",
        label="High-Coverage Substrates"
    )

    ax.set(
        xlabel="Substrates",
        ylabel="Kinase",
        yticks=np.arange(0.5, 3 * len(kinase_counts), 3)
    )
    ax.set_yticklabels(
        kinase_counts.index,
        ha="right",
        ma="right",
        va="center"
    )
    ax.legend()
    ax.set()

    return fig


if __name__ == "__main__":
    make_figure()
