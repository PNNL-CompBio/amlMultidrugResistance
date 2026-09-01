"""Plots figure 3Se: Acetylation Enrichment."""

from os.path import abspath, dirname

from scipy.stats import pearsonr
import numpy as np
import pandas as pd

from pilot.data_import import import_global
from pilot.decomposition import run_acetyl_plsr
from pilot.gene_analysis import preranked_enrichment

FIG_PATH = abspath(dirname(__file__))


def make_figure():
    # Run Acetyl PLSR
    components, _, r2x = run_acetyl_plsr()
    a_scores = components["Acetyl Scores"]
    p_loadings = components["Proteomic Loadings"]

    # Import data
    prot = pd.concat(import_global())
    prot = prot.loc[a_scores.index, :]

    # Get correlations
    corr = pd.DataFrame(
        np.nan,
        index=a_scores.columns,
        columns=prot.columns,
        dtype=float
    )
    for gene in prot.columns:
        gene_col = prot.loc[:, gene].dropna()
        res = pearsonr(
            pd.concat([gene_col] * a_scores.shape[1], axis=1),
            a_scores.loc[gene_col.index, :]
        )
        corr.loc[:, gene] =  res.statistic

    # Iterate through components, perform enrichment
    # We do two sets of enrichment here--one that runs against the loadings,
    # and one that runs against the correlation of the scores with genes.
    # The former gives a sense of what shifts the acetylation causes, while the
    # latter provides more insight into how the changes manifest in patients.
    for comp in p_loadings.columns:
        # Run enrichment using loadings
        preranked_enrichment(
            p_loadings.loc[:, comp],
            f"{FIG_PATH}/{comp}_deflated_loading_enrichment",
            gene_libraries=[
                "GO_Cellular_Component_2025",
                "GO_Molecular_Function_2025",
                "KEGG_2026",
                "GO_Biological_Process_2025",
                "Reactome_Pathways_2024",
                "MSigDB_Hallmark_2020"
            ]
        )
        # Run enrichment using correlations
        preranked_enrichment(
            corr.loc[comp, :],
            f"{FIG_PATH}/{comp}_global_corr_enrichment",
            gene_libraries=[
                "GO_Cellular_Component_2025",
                "GO_Molecular_Function_2025",
                "KEGG_2026",
                "GO_Biological_Process_2025",
                "Reactome_Pathways_2024",
                "MSigDB_Hallmark_2020"
            ]
        )

    return None


if __name__ == "__main__":
    make_figure()
