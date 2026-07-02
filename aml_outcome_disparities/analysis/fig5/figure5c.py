"""Plots figure 5c: KSEA substrate comparison across mutations."""

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

KINASES = [
    "AURKB",
    "CDK1",
    "ATR",
    "PRKCI"
]
COLORS = {
    "White": [
        "#9467bd",
        "#b495d1",
        "#d4c2e5"
    ],
    "Black": [
        "#d62728",
        "#e36768",
        "#efa8a8"
    ]
}


def make_figure():
    # Load Phosphosite+ database, phospho, meta
    syn = syn_login()
    phospho = pd.concat(import_phospho(syn))
    meta = import_meta(syn)
    substrates = pd.read_csv(syn.get("syn73849653").path)

    # Reformat names to match Synapse data, trim to kinase names
    substrates.index = (
        substrates.loc[:, "SUB_GENE"]
        + "-"
        + substrates.loc[:, "SUB_MOD_RSD"]
        + substrates.loc[:, "SUB_MOD_RSD"].str[0].str.lower()
    )
    substrates = substrates.loc[:, "GENE"]

    # Trim to phosphosites in both datasets
    substrates = substrates.loc[substrates.index.intersection(phospho.columns)]
    phospho = phospho.loc[
        :,
        substrates.index[~substrates.index.duplicated()]
    ]

    # Drop duplicates
    # Index is reset to check index-value combo
    substrates = substrates.loc[
        ~substrates.reset_index().duplicated().values
    ]

    # Format mutations, trim to WT, NPM1, and Double mutants
    meta = meta.loc[meta.loc[:, "Race"].isin(["Black", "White"])]
    meta = meta.loc[
        ~np.logical_and(
            meta.loc[:, "FLT3_ITD"] == "Mutant",
            meta.loc[:, "NPM1"] != "Mutant"
        ),
        :
    ]
    meta.loc[:, "Mutation"] = "WT"
    meta.loc[
        meta.loc[:, "NPM1"] == "Mutant",
        "Mutation"
    ] = "NPM1"
    meta.loc[
        np.logical_and(
            meta.loc[:, "NPM1"] == "Mutant",
            meta.loc[:, "FLT3_ITD"] == "Mutant"
        ),
        "Mutation"
    ] = "Both"

    # Trim meta to patients in phospho
    meta = meta.loc[
        phospho.index.intersection(meta.index),
        :
    ]
    phospho = phospho.loc[
        meta.index,
        :
    ]

    # Setup plot
    fig, axes = get_setup(
        len(KINASES),
        2,
        fig_params={
            "figsize": (16, 4 * len(KINASES)),
            "sharey": True
        }
    )

    # Iterate through kinases
    for row_index, kinase in enumerate(KINASES):
        # Trim to kinase substrates
        kinase_phospho = phospho.loc[
            :,
            substrates.loc[substrates == kinase].index
        ].sort_index(axis=1)

        # If more than 10 phosphosites, keep only 10 with least missingness
        if kinase_phospho.shape[1] > 10:
            kinase_phospho = kinase_phospho.loc[
                :,
                kinase_phospho.notna().mean(axis=0).sort_values(
                    ascending=False
                ).index[:10]
            ]

        # Iterate through Black, White patients
        for col_index, race in enumerate(["Black", "White"]):
            ax = axes[row_index, col_index]
            race_phospho = kinase_phospho.loc[
                meta.loc[:, "Race"] == race
            ]
            race_mutations = meta.loc[
                race_phospho.index,
                "Mutation"
            ]

            # Plot each substrate, comparing across mutations
            for x_index, phosphosite in enumerate(race_phospho.columns):
                phospho_col = race_phospho.loc[:, phosphosite].dropna()
                _mutations = race_mutations.loc[phospho_col.index]
                wt_bp = ax.boxplot(
                    phospho_col.loc[_mutations == "WT"],
                    positions=[4 * x_index],
                    widths=0.5,
                    patch_artist=True,
                    notch=True,
                    medianprops={"color": "black"},
                    flierprops={
                        "color": f"{COLORS[race][0]}"
                    },
                )
                npm1_bp = ax.boxplot(
                    phospho_col.loc[_mutations == "NPM1"],
                    positions=[4 * x_index + 1],
                    widths=0.5,
                    patch_artist=True,
                    notch=True,
                    medianprops={"color": "black"},
                    flierprops={
                        "color": f"{COLORS[race][1]}"
                    }
                )
                both_bp = ax.boxplot(
                    phospho_col.loc[_mutations == "Both"],
                    positions=[4 * x_index + 2],
                    widths=0.5,
                    patch_artist=True,
                    notch=True,
                    medianprops={"color": "black"},
                    flierprops={
                        "color": f"{COLORS[race][2]}"
                    }
                )

                wt_bp["boxes"][0].set_facecolor(
                    f"{COLORS[race][0]}"
                )
                npm1_bp["boxes"][0].set_facecolor(
                    f"{COLORS[race][1]}"
                )
                both_bp["boxes"][0].set_facecolor(
                    f"{COLORS[race][2]}"
                )

            for spacer, comparison in enumerate([
                ["WT", "NPM1"],
                ["WT", "Both"]
            ]):
                ax = compare_means(
                    race_phospho.loc[race_mutations.isin(comparison)],
                    race_mutations.loc[race_mutations.isin(comparison)],
                    1 + spacer,
                    4,
                    ax,
                    star_only=True,
                    bracket_space=spacer / 10 + 0.05
                )

            # Reformat plot, add labels
            ax.set(
                xticks=np.arange(0.5, 4 * race_phospho.shape[1], 4),
                ylabel="Log2 Expression",
                title=f"{race}: {kinase}"
            )
            ax.set_xticklabels(
                race_phospho.columns,
                rotation=45,
                ha="right",
                ma="right",
                va="top"
            )

            # Add dashed line at 0
            x_lims = ax.get_xlim()
            ax.plot(
                [-1000, 1000],
                [0, 0],
                linestyle="--",
                color="k",
                zorder=-3
            )
            ax.set_xlim(x_lims)

            # Add legend
            legend_elements = [
                Patch(
                    facecolor=f"{COLORS[race][index]}",
                    edgecolor=f"{COLORS[race][index]}",
                    label=mutation,
                ) for index, mutation in enumerate(
                    ["WT", "NPM1", "Both"]
                )
            ]
            ax.legend(handles=legend_elements)

    return fig


if __name__ == "__main__":
    make_figure()
