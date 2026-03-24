"""Plots figure 2g: PRKACA enrichment analyses."""

from scipy.stats import false_discovery_control, pearsonr
import numpy as np
import pandas as pd

from pilot.data_import import import_global, import_rna
from pilot.gene_analysis import get_prkaca_score, preranked_enrichment


def make_figure():
    # Get PRKACA scores
    prkaca_scores = get_prkaca_score()

    # Load data
    rna = pd.concat(import_rna())
    rna = rna.loc[rna.index.intersection(prkaca_scores.index), :]
    rna = rna.loc[:, ~rna.columns.duplicated()]
    prot = pd.concat(import_global())
    prot = prot.loc[prot.index.intersection(prkaca_scores.index), :]

    for name, dataset in zip(
        ["global", "rna"],
        [prot, rna]
    ):
        ranks = pd.Series(
            0,
            dtype=float,
            index=dataset.columns
        )
        p_vals = ranks.copy()

        # Run Pearson correlation for pre-ranking
        for gene in dataset.columns:
            gene_col = dataset.loc[:, gene].dropna()
            res = pearsonr(
                gene_col,
                prkaca_scores.loc[gene_col.index]
            )
            ranks.loc[gene] = res.statistic
            p_vals.loc[gene] = res.pvalue

        p_vals.loc[:] = false_discovery_control(p_vals)

        # Run enrichment, save dotplots and results to /enrichment
        preranked_enrichment(
            ranks,
            output_path=f"enrichment/{name}/",
            gene_libraries=[
                "GO_Biological_Process_2025",
                "Metabolomics_Workbench_Metabolites_2022",
                "GO_Cellular_Component_2025",
                "GO_Molecular_Function_2025",
                "MSigDB_Hallmark_2020",
                "MSigDB_Oncogenic_Signatures"
            ],
        )


if __name__ == "__main__":
    make_figure()
