"""Plots figure S1a: Age-Race interaction comparison."""
import os
import shutil
import sys
from os.path import abspath, dirname, join

sys.path.append(
    join(dirname(dirname(dirname(abspath(__file__)))), "src", "python")
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from pilot.data_import import (import_acetyl, import_global, import_lipids,
                               import_meta, import_metabolites, import_phospho,
                               import_rna, syn_login)
from pilot.figures.figure_setup import get_setup

REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def make_figure():
    # Make directory for R volcanoes
    os.makedirs(
        join(
            REPO_PATH,
            "analysis",
            "figS1",
            "age_volcano",
        ),
        exist_ok=True
    )

    # Load datasets
    syn = syn_login()
    datasets = {
        "Acetylomics": import_acetyl(syn),
        "Phosphoproteomics": import_phospho(syn),
        "Lipidomics": import_lipids(syn),
        "Global Proteomics": import_global(syn),
        "Metabolites": import_metabolites(syn),
        "Transcriptomics": import_rna(syn)
    }
    meta = import_meta(syn)

    # Convert meta-data to R
    with localconverter(ro.default_converter + pandas2ri.converter):
        meta_r = ro.conversion.py2rpy(meta)

    # Source R function
    r_source = ro.r["source"]
    r_source(
        join(
            REPO_PATH,
            "src",
            "r",
            "proteomics",
            "limma_age.R"
        )
    )
    compare_age = ro.globalenv["compare_age"]

    # Setup figure for plotting
    fig, axes = get_setup(2, 3, fig_params={"figsize": (9, 6)})
    axes = axes.flatten()

    # Iterate through datasets, plot R volcanoes
    for ax, (name, dataset) in zip(axes, datasets.items()):
        # Concatenate BeatAML, Pilot cohorts
        dataset = pd.concat(dataset)
        if name == "Transcriptomics":
            dataset = np.log2(dataset + 1)

        # Convert dataset to R
        with localconverter(ro.default_converter + pandas2ri.converter):
            dataset = ro.conversion.py2rpy(dataset)

        # Run R function
        res = compare_age(
            dataset,
            meta_r,
            0.05,
            name
        )
        name = name.lower().replace(" ", "_")

        # Move R image to local analysis directory
        shutil.move(
            join(
                REPO_PATH,
                "src",
                "r",
                "proteomics",
                f"{name}_volcano.png"
            ),
            join(
                REPO_PATH,
                "analysis",
                "figS1",
                "age_volcano",
                f"{name}_volcano.png"
            )
        )

        # Read R image file, show
        volcano_image = plt.imread(
            join(
                REPO_PATH,
                "analysis",
                "figS1",
                "age_volcano",
                f"{name}_volcano.png"
            )
        )
        ax.imshow(volcano_image)
        ax.set_axis_off()

    return fig


if __name__ == "__main__":
    make_figure()
