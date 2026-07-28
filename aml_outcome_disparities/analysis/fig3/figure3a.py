"""Plots figure 3a: S100A8 Acetylation vs. Race."""

from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from pilot.data_import import import_acetyl, import_meta, syn_login
from pilot.figures.figure_setup import compare_means, get_setup


def make_figure():
    # Import meta-data, acetyl
    syn = syn_login()
    meta = import_meta(syn)
    acetyl = pd.concat(import_acetyl(syn))

    # Trim to S100A8 sites
    acetyl = acetyl.loc[
        :,
        acetyl.columns.str.startswith("S100A8")
    ]

    # Trim to only Black and White patients
    race = meta.loc[acetyl.index, "Race"]
    acetyl = acetyl.loc[race.isin(["White", "Black"]), :]
    race = race.loc[acetyl.index]

    # Setup figure
    fig, ax = get_setup(1, 1, fig_params={"figsize": (2, 3)})

    # Iterate through each acetyl site, compare Black and White patients
    for x_index, acetyl_site in enumerate(acetyl.columns):
        acetyl_col = acetyl.loc[:, acetyl_site].dropna()
        _race = race.loc[acetyl_col.index]
        white_bp = ax.boxplot(
            acetyl_col.loc[_race == "White"],
            positions=[3 * x_index],
            widths=0.8,
            patch_artist=True,
            notch=True,
            medianprops={"color": "black"},
            flierprops={"color": "tab:purple"},
        )
        black_bp = ax.boxplot(
            acetyl_col.loc[_race == "Black"],
            positions=[3 * x_index + 1],
            widths=0.8,
            patch_artist=True,
            notch=True,
            medianprops={"color": "black"},
            flierprops={"color": "tab:red"},
        )
        white_bp["boxes"][0].set_facecolor("tab:purple")
        black_bp["boxes"][0].set_facecolor("tab:red")

    # Compare means between acetyl sites
    ax = compare_means(
        acetyl,
        race,
        1,
        3,
        ax,
        star_only=True,
        alternative="greater"
    )

    # Reformat ticks and labels
    ax.set(
        xticks=np.arange(0.5, 3 * acetyl.shape[1], 3),
        ylabel="Log2 Expression",
        ylim=(-8, 8)
    )
    ax.set_xticklabels(
        acetyl.columns, rotation=45, ha="right", ma="right", va="top"
    )

    # Add legend for Black & White patients
    legend_elements = [
        Patch(
            facecolor="tab:purple",
            edgecolor="tab:purple",
            label="White Patients",
        ),
        Patch(facecolor="tab:red", edgecolor="tab:red", label="Black Patients")
    ]
    ax.legend(handles=legend_elements, loc="lower left")

    return fig


if __name__ == "__main__":
    make_figure()
