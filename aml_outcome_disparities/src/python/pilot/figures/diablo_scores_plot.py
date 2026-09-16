"""
Plots swarmplots of DIABLO scores, colored by patient race.
"""

from os.path import abspath, dirname, join

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from pilot.data_import import import_meta
from pilot.figure_setup import get_setup

REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def main():
    meta = import_meta()
    phospho_scores = pd.read_csv(
        join(REPO_PATH, "data", "diablo_phospho_scores.csv"), index_col=0
    )
    global_scores = pd.read_csv(
        join(REPO_PATH, "data", "diablo_global_scores.csv"), index_col=0
    )

    meta = meta.loc[phospho_scores.index, :]
    race = meta.loc[meta.loc[:, "Race"].isin(["White", "Black"]), "Race"].squeeze()

    fig, axes = get_setup(
        1, phospho_scores.shape[1], {"figsize": (phospho_scores.shape[1] * 2, 2)}
    )
    for ax, comp in zip(axes, phospho_scores.columns):
        sns.swarmplot(
            x=0,
            y=phospho_scores.loc[:, comp],
            hue=race,
            hue_order=["Black", "White"],
            size=3,
            palette={"Black": "tab:blue", "White": "tab:orange"},
            alpha=0.5,
            ax=ax,
        )
        sns.swarmplot(
            x=1,
            y=global_scores.loc[:, comp],
            hue_order=["Black", "White"],
            hue=race,
            legend=False,
            size=3,
            palette={"Black": "tab:blue", "White": "tab:orange"},
            edgecolor="black",
            alpha=0.5,
            ax=ax,
        )
        ax.set_xticklabels(
            ["Phosphoproteomic", "Global"],
            rotation=0,
            ha="center",
            ma="center",
            va="top",
        )
        ax.set_ylabel("DIABLO Score")
        ax.set_title(f"Component {comp[-1]}")

    plt.show()


if __name__ == "__main__":
    main()
