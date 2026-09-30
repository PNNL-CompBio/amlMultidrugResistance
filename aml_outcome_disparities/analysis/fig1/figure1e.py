"""Plots figure 1e: Dataset Missingness Across Variables."""

from collections import OrderedDict
from os.path import abspath, dirname

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

from pilot.data_import import (
    import_acetyl,
    import_global,
    import_lipids,
    import_meta,
    import_metabolites,
    import_phospho,
    import_rna,
    syn_login,
)
from pilot.figures.figure_setup import get_setup

FILE_DIR = dirname(abspath(__file__))
META_VARS = ["Study", "Source", "NPM1", "FLT3_ITD", "Sex", "Race", "Age"]
COLORS = {
    "Race": OrderedDict(
        [
            ("Black", "tab:red"),
            ("White", "tab:purple"),
            ("HispNative","#393b79"),
            ("Asian", "#ce6dbd"),
            ("Multiracial", "#ad494a"),
            ("Pacific Islander", "tab:cyan"),
            ("Unknown", "lightgrey")
        ]
    ),
    "Sex": OrderedDict(
        [
            ("Male", "#8c6d31"),
            ("Female", "#bd9e39")
        ]
    ),
    "Source": OrderedDict(
        [
            ("BeatAML", "#2ca02c"),
            ("newSamples", "#98df8a")
        ]
    ),
    "Study": OrderedDict(
        [
            ("BeatAML", "#1f77b4"),
            ("pilotStudy", "#ff7f0e")
        ]
    ),
    "NPM1": OrderedDict(
        [
            ("Mutant", "#6b6ecf"),
            ("WT", "#9c9ede")
        ]
    ),
    "FLT3_ITD": OrderedDict(
        [
            ("Mutant", "#637939"),
            ("WT", "#8ca252")
        ]
    )
}


def make_figure():
    # Get synapse client, import metadata
    syn = syn_login()
    meta = import_meta(syn)

    # Sort patients by Study, Source, Race, then Age
    meta.sort_values(
        META_VARS,
        inplace=True
    )
    meta = meta.loc[:, META_VARS]

    # Get datasets
    datasets = {
        "RNA-seq": import_rna,
        "Global Proteomics": import_global,
        "Acetylomics": import_acetyl,
        "Phosphoproteomics": import_phospho,
        "Metabolomics": import_metabolites,
        "Lipidomics": import_lipids
    }

    # Track missingness
    present_data = pd.DataFrame(
        np.nan,
        index=list(datasets.keys()),
        columns=meta.index,
        dtype=float
    )
    for dataset_name, dataset_func in datasets.items():
        data = pd.concat(dataset_func(syn))
        present_data.loc[
            dataset_name,
            data.index.intersection(present_data.columns)
        ] = 1

    # Setup figure
    fig, axes = get_setup(
        len(META_VARS) + 1,
        1,
        fig_params={
            "dpi": 500,
            "figsize": (7, 2.5),
            "height_ratios": [1] * len(META_VARS) + [6]
        }
    )

    # Plot categorical variables
    for ax, variable in zip(axes, META_VARS[:-1]):
        # Get variable, recolor
        var_col = meta.loc[:, variable]
        colors = COLORS[variable]
        sub = {key: index for index, key in enumerate(colors.keys())}
        var_col = var_col.replace(sub).to_frame().astype(float).T

        # Convert variable coloring to colormap for heatmap
        cmap = LinearSegmentedColormap.from_list(
            "categories",
            list(colors.values())
        )

        # Plot variable
        heatmap = sns.heatmap(
            var_col,
            ax=ax,
            cmap=cmap,
            cbar=False,
            linewidths=0.05
        )
        heatmap.set_facecolor("lightgrey")

        # Remove ticks
        ax.set(
            xticks=[],
            yticks=[0.5]
        )
        ax.set_yticklabels([variable], rotation=0)

    # Plot age
    ax = axes[-2]
    heatmap = sns.heatmap(
        meta.loc[:, "Age"].astype(float).to_frame().T,
        ax=ax,
        cmap="Oranges",
        linewidths=0.05
    )
    heatmap.set_facecolor("lightgrey")
    ax.set(
        xticks=[],
        yticks=[0.5],
    )
    ax.set_yticklabels(["Age"], rotation=0)

    # Plot data missingness
    ax = axes[-1]
    heatmap = sns.heatmap(
        present_data,
        ax=ax,
        cmap="gray",
        cbar=False,
        linewidths=0.05,
        vmin=0,
        vmax=100
    )
    heatmap.set_facecolor("lightgrey")
    ax.set(
        xticks=[],
        yticks=np.arange(0.5, present_data.shape[0], 1),
        yticklabels=present_data.index
    )

    return fig


if __name__ == "__main__":
    make_figure()
