"""
Plots heatmap of loadings for DIABLO components.
"""
from os.path import abspath, dirname, join

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from pilot.data_import import import_meta
from pilot.figure_setup import get_setup
from pilot.utils import reorder_table

REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def main():
    meta = import_meta()
    phospho_scores = pd.read_csv(
        join(REPO_PATH, "data", "diablo_phospho_scores.csv"), index_col=0
    )
    global_scores = pd.read_csv(
        join(REPO_PATH, "data", "diablo_global_scores.csv"), index_col=0
    )
    phospho_loadings = pd.read_csv(
        join(REPO_PATH, "data", "diablo_phospho_loadings.csv"), index_col=0
    )
    global_loadings = pd.read_csv(
        join(REPO_PATH, "data", "diablo_global_loadings.csv"), index_col=0
    )

    bridge_index = phospho_scores.index.str.contains("Bridge")
    meta = meta.loc[phospho_scores.index.str.replace("-Bridge", ""), :]
    meta.loc[~meta.loc[:, "Race"].isin(["White", "Black"]), "Race"] = "White"
    race = meta.loc[:, "Race"].squeeze()
    race.loc[bridge_index] = "Bridge"

    bridge_samples = race.loc[race == "Bridge"].index
    bridge_samples = pd.Series(bridge_samples + "-Bridge", index=bridge_samples)
    race.rename(index=bridge_samples, inplace=True)
    race = race.sort_values(ascending=True)

    result_sets = [(phospho_scores, phospho_loadings), (global_scores, global_loadings)]
    fig, axes = get_setup(1, 4, {"figsize": (8, 4)})

    for col_index, (scores, loadings) in enumerate(result_sets):
        scores = scores.loc[race.index, :]
        scores /= abs(scores).max(axis=0)
        loadings = reorder_table(loadings)
        sns.heatmap(
            scores,
            vmin=-1,
            vmax=1,
            cmap="coolwarm",
            ax=axes[col_index],
            cbar=col_index == len(result_sets) - 1,
        )
        sns.heatmap(loadings, center=0, cmap="coolwarm", ax=axes[col_index + 2])
        axes[col_index].set_xticklabels(np.arange(scores.shape[1]) + 1, rotation=0)
        axes[col_index + 2].set_xticklabels(
            np.arange(loadings.shape[1]) + 1, rotation=0
        )

    for ax in axes[:2]:
        ax.set_yticks(
            [
                0,
                sum(race == "Black"),
                sum(race != "White"),
                len(race) - sum(race == "White"),
                len(race),
            ]
        )

    plt.show()


if __name__ == "__main__":
    main()
