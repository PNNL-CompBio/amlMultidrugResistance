"""
Runs proteomic enrichment, using DIABLO scores as pre-ranking.
"""

from os.path import abspath, dirname, join

import pandas as pd

from pilot.gene_analysis import preranked_enrichment

LIBRARIES = [
    "CellMarker_2024",
    "KEGG_2021_Human",
    "GO_Biological_Process_2025",
    "GO_Molecular_Function_2025",
]
REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def main():
    global_loadings = pd.read_csv(
        join(REPO_PATH, "data", "diablo_global_loadings.csv"), index_col=0
    )
    for comp in global_loadings.columns:
        genes = global_loadings.loc[:, comp]
        genes = genes.sort_values(ascending=True)
        preranked_enrichment(
            genes,
            join(REPO_PATH, "output", "diablo_enrichment", f"diablo_{comp[-1]}"),
            LIBRARIES,
        )


if __name__ == "__main__":
    main()
