"""Plots figure 7a: MKI67/H3 Acetylation and Race/NPM1."""

from decimal import Decimal

from matplotlib.patches import Patch
from scipy.stats import ttest_ind
import numpy as np
import pandas as pd

from pilot.data_import import import_acetyl, import_meta, syn_login
from pilot.figures.figure_setup import get_setup

COLORS = {
    "White": [
        "#9467bd",
        "#b495d1"
    ],
    "Black": [
        "#d62728",
        "#e36768"
    ]
}


def make_figure():
    # Import meta-data, acetyl
    syn = syn_login()
    meta = import_meta(syn)
    acetyl = pd.concat(import_acetyl(syn))

    # Trim to S100A8 sites
    acetyl = acetyl.loc[
        :,
        np.logical_or(
            acetyl.columns.str.startswith("MKI67"),
            acetyl.columns.str.startswith("HIST1H3")
        )
    ]

    # Trim to Black/White AML patients, non-FLT3-ITD mutants
    meta = meta.loc[
        np.logical_and(
            meta.loc[:, "Race"].isin(["White", "Black"]),
            meta.loc[:, "FLT3_ITD"] != "Mutant"
        ),
        :
    ]

    # Trim to only Black and White patients
    acetyl = acetyl.loc[meta.index.intersection(acetyl.index), :].mean(axis=1)
    meta = meta.loc[acetyl.index, :]

    # Setup figure
    fig, ax = get_setup(1, 1, fig_params={"figsize": (2, 3)})

    # Setup storage for p-values
    p_values = pd.Series(
        np.nan,
        index=["Black", "White"]
    )

    # Iterate through each acetyl site, compare Black and White patients
    for x_index, race in enumerate(["Black", "White"]):
        _meta = meta.loc[meta.loc[:, "Race"] == race, :]
        _acetyl = acetyl.loc[
            _meta.index
        ].dropna()
        _meta = _meta.loc[_acetyl.index, :]

        wt_bp = ax.boxplot(
            _acetyl.loc[_meta.loc[:, "NPM1"] != "Mutant"],
            positions=[3 * x_index],
            widths=0.8,
            patch_artist=True,
            notch=True,
            medianprops={"color": "black"},
            flierprops={"color": COLORS[race][0]},
        )
        mutant_bp = ax.boxplot(
            _acetyl.loc[_meta.loc[:, "NPM1"] == "Mutant"],
            positions=[3 * x_index + 1],
            widths=0.8,
            patch_artist=True,
            notch=True,
            medianprops={"color": "black"},
            flierprops={"color": COLORS[race][1]},
        )
        wt_bp["boxes"][0].set_facecolor(COLORS[race][0])
        mutant_bp["boxes"][0].set_facecolor(COLORS[race][1])

        res = ttest_ind(
            _acetyl.loc[_meta.loc[:, "NPM1"] == "Mutant"],
            _acetyl.loc[_meta.loc[:, "NPM1"] != "Mutant"],
            nan_policy="omit"
        )
        p_values.loc[race] = res.pvalue

    # Reformat ticks and labels
    ax.set(
        xticks=[0.5, 3.5],
        xticklabels=[
            f"Black: {'{:.2E}'.format(Decimal(p_values.loc['Black']))}",
            f"White: {'{:.2E}'.format(Decimal(p_values.loc['White']))}"
        ],
        ylabel="MKI67 + H3 Acetylation",
        ylim=(-4, 4)
    )

    # Add legend for Black & White patients
    legend_elements = [
        Patch(
            facecolor=COLORS["White"][0],
            edgecolor=COLORS["White"][0],
            label="White: WT",
        ),
        Patch(
            facecolor=COLORS["White"][1],
            edgecolor=COLORS["White"][1],
            label="White: NPM1+",
        ),
        Patch(
            facecolor=COLORS["Black"][0],
            edgecolor=COLORS["Black"][0],
            label="Black: WT",
        ),
        Patch(
            facecolor=COLORS["Black"][1],
            edgecolor=COLORS["Black"][1],
            label="Black: NPM1+",
        ),
    ]
    ax.legend(handles=legend_elements, loc="lower right")

    return fig


if __name__ == "__main__":
    make_figure()
