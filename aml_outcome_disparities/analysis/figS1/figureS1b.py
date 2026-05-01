"""Plots figure 3c: Compares young AML cases."""
import os
import sys
from os.path import abspath, dirname, join

sys.path.append(
    join(dirname(dirname(dirname(abspath(__file__)))), "src", "python")
)

import numpy as np
import pandas as pd
import rpy2.robjects as ro
import seaborn as sns
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from pilot.data_import import (import_acetyl,import_global, import_lipids,
                               import_meta, import_metabolites, import_phospho,
                               import_rna, syn_login)
from pilot.gene_analysis import preranked_enrichment
from pilot.figures.figure_setup import get_setup

AGE_THRESHOLD = 40
REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def make_figure():
    # Make directories for enrichment results
    os.makedirs(
        join(
            REPO_PATH,
            "analysis",
            "figS1",
            "young_aml",
            "rna"
        ),
        exist_ok=True
    )
    os.makedirs(
        join(
            REPO_PATH,
            "analysis",
            "figS1",
            "young_aml",
            "global"
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
        "Transcriptomics": import_rna(syn),
    }
    meta = import_meta(syn)

    # Trim to patients in the age range of 15-39, inclusive
    meta = meta.loc[
        meta.loc[:, "Age"].between(15, AGE_THRESHOLD, inclusive="left"),
        :
    ]
    meta = meta.loc[
        meta.loc[:, "Race"].isin(["Black", "White"]),
        :
    ]

    # Source R LIMMA function
    r_source = ro.r["source"]
    r_source(
        join(
            REPO_PATH,
            "src",
            "r",
            "proteomics",
            "limma_correlations.R"
        )
    )
    run_limma = ro.globalenv["run_limma"]

    # Store differentially expressed measurements across 'omes
    diffs = pd.DataFrame(
        0,
        dtype=int,
        index=list(datasets.keys()),
        columns=["Black", "White"]
    )

    # Iterate through datasets
    for name, dataset in datasets.items():
        # Concatenate BeatAML, Pilot samples, trim meta-data to samples in
        # present dataset
        dataset = pd.concat(dataset)
        dataset = dataset.loc[
            dataset.index.intersection(meta.index),
            :
        ]
        _meta = meta.loc[dataset.index, :]

        # Log scale RNA-Seq
        if name == "Transcriptomics":
            dataset = np.log2(dataset + 1)

        # Convert to R
        with localconverter(ro.default_converter + pandas2ri.converter):
            dataset = ro.conversion.py2rpy(dataset)
            _meta = ro.conversion.py2rpy(_meta)

        # Run R LIMMA function
        res = run_limma(
            dataset,
            _meta,
            0.05
        )

        # Convert result back to Python
        with localconverter(ro.default_converter + pandas2ri.converter):
            res, _ = ro.conversion.rpy2py(res)

        # Update differentially expressed counts for 'ome
        diffs.loc[name, "Black"] = sum(
            np.logical_and(
                res.loc[:, "logFC"] > 0,
                res.loc[:, "adj.P.Val"] < 0.05
            )
        )
        diffs.loc[name, "White"] = sum(
            np.logical_and(
                res.loc[:, "logFC"] < 0,
                res.loc[:, "adj.P.Val"] < 0.05
            )
        )

        # If RNA-seq or Global, store genes for enrichment analysis
        if name == "Transcriptomics":
            rna_genes = res.loc[:, "logFC"].squeeze()
            rna_genes = rna_genes * -np.log10(
                res.loc[:, "P.Value"]
            ).sort_values()
        elif name == "Global Proteomics":
            global_genes = res.loc[:, "logFC"].squeeze()
            global_genes = global_genes * -np.log10(
                res.loc[:, "P.Value"]
            )

    # Setup plot
    fig, ax = get_setup(1, 1, fig_params={"figsize": (3, 2)})

    # Plot heatmap comparing counts of differentially-expressed measurements
    # across Black and White patients
    sns.heatmap(
        diffs,
        annot=True,
        ax=ax,
        vmin=0,
        cmap="Reds",
        cbar_kws={
            "label": "Differentially Expressed\nMeasurements"
        }
    )

    # Reformat axes
    ax.set(
        xticks=np.arange(0.5, diffs.shape[1], 1),
        yticks=np.arange(0.5, diffs.shape[0], 1),
        xticklabels=diffs.columns,
        yticklabels=diffs.index
    )

    # Perform enrichment on RNA-seq, Global genes
    preranked_enrichment(
        rna_genes,
        join(
            REPO_PATH,
            "analysis",
            "figS1",
            "young_aml",
            "rna"
        ),
        [
            "GO_Biological_Process_2025",
            "Metabolomics_Workbench_Metabolites_2022",
            "GO_Cellular_Component_2025",
            "GO_Molecular_Function_2025",
            "MSigDB_Hallmark_2020",
            "MSigDB_Oncogenic_Signatures"
        ]
    )
    preranked_enrichment(
        global_genes,
        join(
            REPO_PATH,
            "analysis",
            "figS1",
            "young_aml",
            "global"
        ),
        [
            "GO_Biological_Process_2025",
            "Metabolomics_Workbench_Metabolites_2022",
            "GO_Cellular_Component_2025",
            "GO_Molecular_Function_2025",
            "MSigDB_Hallmark_2020",
            "MSigDB_Oncogenic_Signatures"
        ]
    )

    return fig


if __name__ == "__main__":
    make_figure()
