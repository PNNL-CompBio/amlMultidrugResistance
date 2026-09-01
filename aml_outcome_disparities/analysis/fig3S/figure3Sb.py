"""Plots figure 3Sb: PLSR & Race Correlations."""

from decimal import Decimal

from scipy.stats import false_discovery_control, ttest_ind
import numpy as np
import seaborn as sns

from pilot.data_import import import_meta
from pilot.decomposition import run_acetyl_plsr
from pilot.figures.figure_setup import get_setup


def make_figure():
    # Run acetyl PLSR, load meta-data, and get acetyl scores
    components, _, r2x = run_acetyl_plsr()
    a_scores = components["Acetyl Scores"]
    meta = import_meta()
    meta = meta.loc[a_scores.index, :]

    # Compare components across race, correct for FDR
    race_res = ttest_ind(
        a_scores.loc[
            meta.loc[:, "Race"] == "Black"
        ],
        a_scores.loc[
            meta.loc[:, "Race"] == "White"
        ],
        nan_policy="omit"
    )
    p_values = false_discovery_control(race_res.pvalue)

    # Trim to Black/White patients for plotting
    meta = meta.loc[meta.loc[:, "Race"].isin(["Black", "White"]), :]
    a_scores = a_scores.loc[meta.index.intersection(a_scores.index), :]
    meta = meta.loc[a_scores.index, :]

    # Reformat component names
    a_scores.columns = [f"PLSR {i}" for i in np.arange(a_scores.shape[1]) + 1]
    a_scores.loc[:, "Race"] = meta.loc[:, "Race"]
    a_scores = a_scores.melt(id_vars="Race", var_name="Component")

    # Setup figure
    fig, ax = get_setup(1, 1, fig_params={"figsize": (6, 3)})
    sns.violinplot(
        a_scores,
        x="Component",
        y="value",
        hue="Race",
        palette=["tab:purple", "tab:red"]
    )

    # Reformat tick labels, label axes
    xtick_labels = [
        f"{tick_label.get_text()}\np-val: {'{:.2E}'.format(Decimal(p_val))}" for
        tick_label, p_val in
        zip(ax.get_xticklabels(), p_values)
    ]
    ax.set(
        ylabel="PLSR Score",
        xticklabels=xtick_labels
    )

    return fig


if __name__ == "__main__":
    make_figure()
