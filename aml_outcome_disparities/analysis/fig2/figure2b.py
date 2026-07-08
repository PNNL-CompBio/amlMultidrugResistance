"""Plots figure 2b: KSEA substrate comparison."""

import sys

from matplotlib.patches import Patch
from os.path import abspath, dirname, join
import numpy as np
import pandas as pd

from pilot.data_import import import_meta, import_phospho, syn_login
from pilot.figures.figure_setup import get_setup, compare_means

sys.path.append(
    join(dirname(dirname(dirname(abspath(__file__)))), "src", "python")
)

ANALYSIS_DIR = dirname(abspath(__file__))


def make_figure():
    # Load KSEA substrates, trim to PRKACA substrates
    substrates = pd.read_csv(join(ANALYSIS_DIR, "Kinase-Substrate Links.csv"))
    substrates = substrates.loc[
        substrates.loc[:, "Kinase.Gene"] == "PRKACA", :
    ].sort_values(by="log2FC", ascending=True)

    # Reformat names to match Synapse data
    substrates.index = (
        substrates.loc[:, "Substrate.Gene"]
        + "-"
        + substrates.loc[:, "Substrate.Mod"]
        + substrates.loc[:, "Substrate.Mod"].str[0].str.lower()
    )

    # Load data
    syn = syn_login()
    phospho = import_phospho(syn)
    phospho = pd.concat(phospho, axis=0)
    meta = import_meta(syn)
    meta = meta.loc[phospho.index, :]

    # Trim KSEA table to phospho substrates that also appear in phospho
    # measurements. Not all are included since the KSEA splits phosphosites with
    # multiple phosphorylations into individual columns!
    shared = substrates.index.intersection(phospho.columns)
    phospho = phospho.loc[:, shared]

    # Drop phospho measurements with high missingness, isolate to Black and
    # White patients
    # MANOVA requires no missing values, so those with a missingness above
    # 1% are dropped to avoid removing too many samples
    phospho = phospho.loc[:, phospho.isna().mean(axis=0) < 0.01]
    race = meta.loc[:, "Race"]
    race = race.loc[race.isin(["Black", "White"])]
    phospho = phospho.loc[race.index, :]

    # Setup plot
    fig, ax = get_setup(1, 1, fig_params={"figsize": (8, 4)})

    # Plot each PRKACA substrate as boxplots comparing expression between Black
    # and White patients
    for x_index, phosphosite in enumerate(phospho.columns):
        phospho_col = phospho.loc[:, phosphosite].dropna()
        _race = race.loc[phospho_col.index]
        white_bp = ax.boxplot(
            phospho_col.loc[_race == "White"],
            positions=[3 * x_index],
            widths=0.5,
            patch_artist=True,
            notch=True,
            medianprops={"color": "black"},
            flierprops={"color": "tab:purple"},
        )
        black_bp = ax.boxplot(
            phospho_col.loc[_race == "Black"],
            positions=[3 * x_index + 1],
            widths=0.5,
            patch_artist=True,
            notch=True,
            medianprops={"color": "black"},
            flierprops={"color": "tab:red"},
        )
        white_bp["boxes"][0].set_facecolor("tab:purple")
        black_bp["boxes"][0].set_facecolor("tab:red")

    ax = compare_means(
        phospho,
        race,
        1,
        3,
        ax,
        star_only=True,
        alternative="less"
    )

    # Reformat plot, add labels, note MANOVA result
    ax.set(
        xticks=np.arange(0.5, 3 * phospho.shape[1], 3),
        ylabel="Log2 Expression",
        ylim=(-8, 8),
    )
    ax.set_xticklabels(
        phospho.columns, rotation=45, ha="right", ma="right", va="top"
    )

    # Add legend for Black & White patients
    legend_elements = [
        Patch(
            facecolor="tab:purple",
            edgecolor="tab:purple",
            label="White Patients",
        ),
        Patch(facecolor="tab:red", edgecolor="tab:red", label="Black Patients"),
    ]
    ax.legend(handles=legend_elements, loc="upper left")

    return fig


if __name__ == "__main__":
    make_figure()
