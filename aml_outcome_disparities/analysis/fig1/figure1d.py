"""Plots figure 1d: Race/LIMMA Comparisons."""

from os.path import abspath, dirname, join
import os

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from pilot.data_import import (import_acetyl, import_global, import_meta,
                               import_phospho, import_rna, import_metabolites,
                               import_lipids)

DATASET_FUNCS = {
    "Transcriptomics": import_rna,
    "Global Proteomics": import_global,
    "Lipidomics": import_lipids,
    "Metabolomics": import_metabolites,
    "Acetylomics": import_acetyl,
    "Phosphoproteomics": import_phospho
}
REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def make_figure():
    # Import metadata
    meta = import_meta()

    # Setup directory for svg files
    os.makedirs(
        "volcano_svgs",
        exist_ok=True
    )

    # Source LIMMA function
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

    # Iterate through datasets
    for dataset_name, data_func in DATASET_FUNCS.items():
        # Use DESeq2 for RNA-seq
        if dataset_name == "Transcriptomics":
            # Import RNA-seq, remove genes without symbols
            data = pd.concat(data_func(tpm=False))
            data = data.loc[:, ~data.columns.str.startswith("ENSG")].astype(int)

            # Trim to patients with metadata
            data_meta = meta.loc[
                data.index.intersection(meta.index),
                :
            ]
            data = data.loc[data_meta.index, :]

            # Convert Python objects to R
            with localconverter(ro.default_converter + pandas2ri.converter):
                data_r = ro.conversion.py2rpy(data)
                data_meta_r = ro.conversion.py2rpy(data_meta)

            # Source DESeq2 function
            r_source = ro.r["source"]
            r_source(
                join(
                    REPO_PATH,
                    "src",
                    "r",
                    "rna",
                    "differential_expression.R"
                )
            )
            run_deseq = ro.globalenv["run_deseq2"]

            # Run DESeq2, save volcano plot
            run_deseq(
                data_r,
                data_meta_r,
                p_cutoff=0.05,
                save_plot=f"volcano_svgs/{dataset_name}.svg"
            )
        # Use LIMMA for other 'omes
        else:
            # Import dataset
            data = pd.concat(data_func())
            data = data.loc[
                data.index.intersection(meta.index),
                :
            ]
            data_meta = meta.loc[data.index, :]

            # Convert Python objects to R
            with localconverter(ro.default_converter + pandas2ri.converter):
                data_r = ro.conversion.py2rpy(data)
                data_meta_r = ro.conversion.py2rpy(data_meta)

            # Run LIMMA, save volcano plot
            run_limma(
                data_r,
                data_meta_r,
                p_cutoff=0.05,
                save_plot=f"volcano_svgs/{dataset_name}.svg"
            )

    return None


if __name__ == "__main__":
    make_figure()
