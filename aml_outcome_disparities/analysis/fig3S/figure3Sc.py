"""Plots figure 3Sc: PLSR, Race, & NPM1."""

from decimal import Decimal

from scipy.stats import false_discovery_control, ttest_ind
import numpy as np
import seaborn as sns

from pilot.data_import import import_meta
from pilot.decomposition import run_acetyl_plsr
from pilot.figures.figure_setup import get_setup

PALETTE = {
    "Black: WT": "tab:red",
    "Black: Mutant": "#EB9394",
    "White: WT": "tab:purple",
    "White: Mutant": "#CAB3DE",
}


def make_figure():
    # Import meta-data, run acetyl PLSR
    meta = import_meta()
    components, _, r2xs = run_acetyl_plsr()

    # Extract scores, reformat component names
    a_scores = components["Acetyl Scores"]
    a_scores.columns = [
        f"PLSR {i}" for i in np.arange(1, 5)
    ]

    # Drop FLT3-ITD mutants
    meta = meta.loc[
        np.logical_and(
            meta.loc[:, "Race"].isin(["Black", "White"]),
            meta.loc[:, "FLT3_ITD"] != "Mutant"
        ),
        :
    ]
    meta = meta.loc[
        a_scores.index.intersection(meta.index),
        :
    ]
    a_scores = a_scores.loc[
        meta.index,
        :
    ]

    # Compare NPM1 across within Black and White populations
    black_res = ttest_ind(
        a_scores.loc[
            np.logical_and(
                meta.loc[:, "NPM1"] == "Mutant",
                meta.loc[:, "Race"] == "Black"
            )
        ],
        a_scores.loc[
            np.logical_and(
                meta.loc[:, "NPM1"] != "Mutant",
                meta.loc[:, "Race"] == "Black"
            )
        ],
        nan_policy="omit"
    )
    white_res = ttest_ind(
        a_scores.loc[
            np.logical_and(
                meta.loc[:, "NPM1"] == "Mutant",
                meta.loc[:, "Race"] == "White"
            )
        ],
        a_scores.loc[
            np.logical_and(
                meta.loc[:, "NPM1"] != "Mutant",
                meta.loc[:, "Race"] == "White"
            )
        ],
        nan_policy="omit"
    )

    # Correct p-values
    p_values = np.vstack(
        [white_res.pvalue, black_res.pvalue]
    )
    p_values = false_discovery_control(p_values)

    # Setup figure
    fig, ax = get_setup(
        1,
        1,
        fig_params={"figsize": (8, 4)},
    )

    # Reformat for plotting via seaborn
    a_scores.loc[:, "Group"] = meta.loc[:, "Race"] + ": " + meta.loc[:, "NPM1"]
    a_scores = a_scores.melt(
        id_vars="Group",
        var_name="Component"
    )

    # Plot boxplots
    sns.violinplot(
        a_scores,
        x="Component",
        y="value",
        hue="Group",
        hue_order=[
            "White: WT",
            "White: Mutant",
            "Black: WT",
            "Black: Mutant"
        ],
        palette=PALETTE,
        ax=ax
    )

    # Reformat tick labels, label axes
    ticks = [
        f"PLSR {comp + 1}\n"
        f"White p-value: {'{:.2E}'.format(Decimal(p_values[0, comp]))}\n"
        f"Black p-value: {'{:.2E}'.format(Decimal(p_values[1, comp]))}"
        for comp in np.arange(4)
    ]
    ax.set_xticklabels(
        ticks,
        ha="center",
        ma="center",
        va="top"
    )
    ax.set(
        xlabel="",
        ylabel="Acetyl Score"
    )

    return fig


if __name__ == "__main__":
    make_figure()
