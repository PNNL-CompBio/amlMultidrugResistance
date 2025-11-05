"""
Plots genes targeted by phosphorylation.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pilot.data_import import import_meta, import_phospho
from pilot.figure_setup import get_setup
from pilot.gene_analysis import calculate_fc


def compare_distributions(genes: pd.Series, base_genes: pd.Series, name: str):
    """
    Compare composition of target genes to base_genes.

    Args:
        genes (pd.Series): Gene counts for target group.
        base_genes (pd.Series): Gene counts in base dataset.
        name (str): Group name for save files.

    Returns:
         None.
    """
    fig, axes = get_setup(
        1, 2, fig_params=dict(figsize=(8, 4), constrained_layout=True)
    )
    genes = genes.sort_values(ascending=False)

    axes[0].bar(np.arange(0, 30, 3), genes.iloc[:10], width=1, color="tab:blue")
    axes[0].bar(
        np.arange(1, 30, 3), base_genes.loc[genes.index[:10]], width=1, color="tab:grey"
    )

    axes[0].set_xticks(np.arange(0.5, 30, 3))
    axes[0].set_xticklabels(
        genes.index[:10], rotation=45, ha="right", ma="right", va="top"
    )
    axes[0].set_ylabel("Number of Phosphosites")
    axes[0].legend(
        [
            f"Significantly Elevated in {name} Patients",
            "Total Phosphosites in Target Gene",
        ]
    )

    gene_prop = genes / base_genes.loc[genes.index]
    gene_prop = gene_prop.loc[base_genes.loc[gene_prop.index] > 10]
    gene_prop = gene_prop.sort_values(ascending=False)

    axes[1].bar(np.arange(10), gene_prop.iloc[:10], color="tab:blue")

    axes[1].set_xticks(np.arange(10))
    axes[1].set_xticklabels(
        gene_prop.index[:10], rotation=45, ha="right", ma="right", va="top"
    )
    axes[1].set_ylabel(
        "Proportion of Gene Phosphosites with\n"
        f"Significant Elevation in {name} Patients"
    )


def main():
    meta = import_meta()

    ptrc, pilot = import_phospho(corrected=True)
    data = pd.concat([ptrc, pilot])

    data = data.loc[data.index.intersection(meta.index), :]
    meta = meta.loc[data.index, :]

    black_patients = data.loc[meta.loc[:, "Race"] == "Black", :]
    white_patients = data.loc[meta.loc[:, "Race"] == "White", :]

    lfc, p_values = calculate_fc(black_patients, white_patients)
    black_over_exp = lfc.loc[np.logical_and(lfc > 0.5, p_values < 0.05)]
    white_over_exp = lfc.loc[np.logical_and(lfc < 0.5, p_values < 0.05)]

    black_over_exp = black_over_exp.index.str.split("-", expand=True).get_level_values(
        0
    )
    white_over_exp = white_over_exp.index.str.split("-", expand=True).get_level_values(
        0
    )
    base_genes = pilot.columns.str.split("-", expand=True).get_level_values(0)

    black_genes = black_over_exp.value_counts()
    white_genes = white_over_exp.value_counts()
    base_genes = base_genes.value_counts()

    compare_distributions(black_genes, base_genes, "Black")
    plt.show()

    compare_distributions(white_genes, base_genes, "White")
    plt.show()


if __name__ == "__main__":
    main()
