"""Plots figure 1c: Cell Type Score Comparisons"""

import sys
from os.path import abspath, dirname, join

sys.path.append(
    join(dirname(dirname(dirname(abspath(__file__)))), "src", "python")
)

import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import false_discovery_control, ttest_ind
from sklearn.preprocessing import scale

from pilot.data_import import import_meta, import_rna, syn_login
from pilot.figures.figure_setup import get_setup
from pilot.gene_analysis import cell_type_scores

FILE_DIR = dirname(abspath(__file__))
RACE_COLORS = {"Black": "tab:red", "White": "tab:purple"}


def make_figure():
    # Import rna
    syn = syn_login()
    rna = import_rna(syn)
    rna = pd.concat(rna, axis=0)
    meta = import_meta(syn)

    rna = rna.loc[rna.index.intersection(meta.index), :]
    rna.loc[:] = scale(rna, axis=1)
    meta = meta.loc[rna.index, :]
    races = meta.loc[:, "Race"]

    # Get cell type scores
    ct_scores = cell_type_scores()
    ct_scores = ct_scores.loc[races.index, :]

    # Setup figure
    fig, axes = get_setup(
        2,
        3,
        fig_params={
            "figsize": (9, 6),
        },
    )
    axes = axes.flatten()

    # Run t-test to compare cell type scores between Black and White patients
    t_res = ttest_ind(
        ct_scores.loc[races == "Black", :],
        ct_scores.loc[races == "White", :],
        axis=0,
    )
    t_stats = pd.Series(t_res.statistic, dtype=float, index=ct_scores.columns)
    p_vals = pd.Series(
        false_discovery_control(t_res.pvalue),
        dtype=float,
        index=ct_scores.columns,
    )

    # Add Race column to data, filter to only black and white patients
    ct_scores.loc[:, "Race"] = races
    ct_scores = ct_scores.loc[
        ct_scores.loc[:, "Race"].isin(["Black", "White"]), :
    ]
    for ax, cell_type in zip(axes, ct_scores.columns[:6]):
        sns.stripplot(
            ct_scores,
            x="Race",
            y=cell_type,
            hue="Race",
            palette=RACE_COLORS,
            ax=ax,
        )

        white_scores = ct_scores.loc[
            ct_scores.loc[:, "Race"] == "White", cell_type
        ]
        black_scores = ct_scores.loc[
            ct_scores.loc[:, "Race"] == "Black", cell_type
        ]
        ax.errorbar(
            0,
            white_scores.mean(),
            yerr=1.96 * (white_scores.std() / np.sqrt(len(white_scores))),
            capsize=10,
            markersize=20,
            linewidth=3,
            markeredgewidth=3,
            marker="_",
            color="k",
            zorder=10,
        )
        ax.errorbar(
            1,
            white_scores.mean(),
            yerr=1.96 * (black_scores.std() / np.sqrt(len(black_scores))),
            capsize=10,
            markersize=20,
            linewidth=3,
            markeredgewidth=3,
            marker="_",
            color="k",
            zorder=10,
        )

        y_lim = ax.get_ylim()
        ax.set(
            ylabel=f"{cell_type} score",
            xlim=(-0.5, 1.5),
            ylim=(y_lim[0], y_lim[1] + 0.1),
        )
        ax.text(
            0.99,
            0.99,
            s=f"t-stat: {round(t_stats.loc[cell_type], 3)}\n"
            f"q-value: {round(p_vals.loc[cell_type], 3)}\n",
            ha="right",
            ma="right",
            va="top",
            transform=ax.transAxes,
        )

    return fig


if __name__ == "__main__":
    make_figure()
