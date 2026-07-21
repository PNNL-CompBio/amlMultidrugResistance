"""Plots figure 2d: PRKACA scoring validation."""

from scipy.stats import ttest_ind
import numpy as np
import seaborn as sns

from pilot.data_import import import_meta
from pilot.figures.figure_setup import get_setup
from pilot.gene_analysis import get_prkaca_score


def make_figure():
    # Load PRKACA scores, trim meta-data to Black/White patients with scores
    prkaca_scores = get_prkaca_score()
    meta = import_meta()
    race = meta.loc[:, "Race"]
    race = race.loc[race.isin(["Black", "White"])]
    prkaca_scores = prkaca_scores.loc[
        prkaca_scores.index.intersection(race.index),
    ]
    race = race.loc[prkaca_scores.index]

    # Setup figure
    fig, ax = get_setup(1, 1, fig_params={"figsize": (3, 3)})

    # T-test on PCA scores between Black and White patients
    res = ttest_ind(
        prkaca_scores.loc[race == "Black"], prkaca_scores.loc[race == "White"]
    )

    # Plot strip plots of PRKACA scores, subset by race
    sns.stripplot(
        x=0,
        y=prkaca_scores.loc[race == "White"],
        alpha=0.5,
        size=3,
        ax=ax,
        color="tab:purple",
    )
    sns.stripplot(
        x=1,
        y=prkaca_scores.loc[race == "Black"],
        alpha=0.5,
        size=3,
        ax=ax,
        color="tab:red",
    )

    # Plot 95% confidence intervals for PRKACA scores
    ax.errorbar(
        0,
        prkaca_scores.loc[race == "White"].mean(),
        yerr=1.96
        * prkaca_scores.loc[race == "White"].std()
        / np.sqrt(sum(race == "White")),
        capsize=2,
        color="black",
        zorder=3,
        linewidth=1,
        markersize=3,
        marker="_",
    )
    ax.errorbar(
        1,
        prkaca_scores.loc[race == "Black"].mean(),
        yerr=1.96
        * prkaca_scores.loc[race == "Black"].std()
        / np.sqrt(sum(race == "Black")),
        capsize=2,
        color="black",
        zorder=3,
        linewidth=1,
        markersize=3,
        marker="_",
    )

    # Label axes, annotate with t-test
    ax.set(
        xlim=(-0.5, 1.5),
        xticks=(0, 1),
        xticklabels=("White", "Black"),
        ylabel="PRKACA Score",
    )
    ax.text(
        0.01,
        0.01,
        f"t-stat: {round(res.statistic, 3)}\n"
        f"p-value: {round(res.pvalue, 3)}",
        ha="left",
        ma="left",
        va="bottom",
        transform=ax.transAxes,
        fontsize=8,
    )

    return fig


if __name__ == "__main__":
    make_figure()
