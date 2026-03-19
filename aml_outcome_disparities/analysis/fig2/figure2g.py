"""Plots figure 2g: PRKACA enrichment analyses."""

from scipy.stats import pearsonr
import pandas as pd

from pilot.data_import import import_global, import_rna
from pilot.gene_analysis import get_prkaca_score, preranked_enrichment


def make_figure():
    # Get PRKACA scores
    prkaca_scores = get_prkaca_score()

    # Load data
    rna = pd.concat(import_rna())
    rna = rna.loc[rna.index.intersection(prkaca_scores.index), :]
    prot = pd.concat(import_global())
    prot = prot.loc[prot.index.intersection(prkaca_scores.index), :]

    for name, dataset in zip(
        ["global", "rna"],
        [prot, rna]
    ):
        # Run Pearson correlation for pre-ranking
        res = pearsonr(
            dataset,
            pd.concat(
                [prkaca_scores.loc[dataset.index]] * dataset.shape[1],
                axis=1
            )
        )
        ranks = pd.Series(res.statistic, index=dataset.columns)

        # Run enrichment, save dotplots and results to /enrichment
        preranked_enrichment(
            ranks,
            output_path=f"enrichment/{name}/",
            gene_libraries=[
                "GO_Biological_Process_2025",
                "Metabolomics_Workbench_Metabolites_2022",
            ],
        )


if __name__ == "__main__":
    make_figure()
