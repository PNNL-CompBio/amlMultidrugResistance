"""
Plots race-associated differences in expression of CDK4 substrates.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.weightstats import ztest

from pilot.data_import import import_meta, import_phospho
from pilot.figure_setup import get_setup


def main():
    meta = import_meta()
    ptrc, pilot = import_phospho(corrected=True)
    phospho = pd.concat([ptrc, pilot])

    cdk4_phospho = ["RB1-S795s", "RB1-T356t", "RB1-T373t", "RBL1-T369t"]
    fig, axes = get_setup(1, len(cdk4_phospho), {"figsize": (2 * len(cdk4_phospho), 2)})
    axes = axes.flatten()
    for ax, target in zip(axes, cdk4_phospho):
        black_target = phospho.loc[meta.loc[:, "Race"] == "Black", target]
        white_target = phospho.loc[meta.loc[:, "Race"] == "White", target]
        stat, p_val = ztest(black_target.dropna(), white_target.dropna())

        sns.swarmplot(x=0, y=black_target, alpha=0.5, size=3, ax=ax)
        sns.swarmplot(x=1, y=white_target, alpha=0.5, size=3, ax=ax)
        ax.errorbar(
            0,
            black_target.mean(),
            yerr=1.96 * black_target.std() / np.sqrt(len(black_target)),
            capsize=2,
            color="black",
            zorder=3,
            linewidth=1,
            markersize=3,
            marker="_",
        )
        ax.errorbar(
            1,
            white_target.mean(),
            yerr=1.96 * white_target.std() / np.sqrt(len(white_target)),
            capsize=2,
            color="black",
            zorder=3,
            linewidth=1,
            markersize=3,
            marker="_",
        )

        p_val = f"{p_val: .2E}"
        ax.text(
            0.99,
            0.01,
            s=f"p-value: {p_val}",
            ha="right",
            ma="right",
            va="bottom",
            transform=ax.transAxes,
        )
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Black\nPatients", "White\nPatients"])
        ax.set_xlim([-0.5, 1.5])

    fig.suptitle("CDK4 Substrates")
    plt.show()


if __name__ == "__main__":
    main()
