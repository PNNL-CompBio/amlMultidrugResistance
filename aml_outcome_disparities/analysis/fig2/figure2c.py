from os.path import abspath, dirname, join
from scipy.stats import false_discovery_control, pearsonr
import numpy as np
import pandas as pd
import seaborn as sns

from pilot.data_import import import_meta, import_phospho, import_rna, syn_login
from pilot.figures.figure_setup import get_setup
from pilot.gene_analysis import cell_type_scores, get_has_scores

ANALYSIS_DIR = dirname(abspath(__file__))


def make_figure():
    # Load KSEA substrate links, trim to PRKACA substrates
    substrates = pd.read_csv(
        join(
            ANALYSIS_DIR,
            "Kinase-Substrate Links.csv"
        )
    )
    substrates = substrates.loc[
        substrates.loc[:, "Kinase.Gene"] == "PRKACA",
        :
    ].sort_values(
        by="log2FC",
        ascending=True
    )

    # Reformat names to match Synapse data
    substrates.index = (substrates.loc[:, "Substrate.Gene"] + "-" +
        substrates.loc[:, "Substrate.Mod"] +
        substrates.loc[:, "Substrate.Mod"].str[0].str.lower()
    )

    # Load data
    syn = syn_login()
    phospho = import_phospho(syn, pre_corrected=False)
    phospho = pd.concat(phospho, axis=0)
    rna = import_rna(syn)
    rna = pd.concat(rna, axis=0)

    meta = import_meta(syn)
    meta = meta.loc[phospho.index, :]
    rna = rna.loc[rna.index.intersection(meta.index), :]

    # Trim KSEA table to phospho substrates that also appear in phospho
    # measurements. Not all are included since the KSEA splits phosphosites with
    # multiple phosphorylations into individual columns!
    shared = substrates.index.intersection(phospho.columns)
    phospho = phospho.loc[:, shared]

    # Collect cell-type scores
    ct_scores = cell_type_scores()

    # Drop phospho measurements with high missingness to match figure 2b
    phospho = phospho.loc[:, phospho.isna().mean(axis=0) < 0.01]

    ############################################################################
    ### Bottomly et al Cell Type Comparisons
    ############################################################################

    # Setup figure
    fig, axes = get_setup(
        2,
        1,
        fig_params={
            "figsize": (8, 5),
            "height_ratios": [3, 1]
        }
    )

    # Setup storage for correlation results
    rhos = pd.DataFrame(
        index=ct_scores.columns,
        columns=phospho.columns,
        dtype=float
    )
    p_values = rhos.copy(deep=True)

    # Evaluate correlations between each phosphosite and each cell type score
    for phosphosite in phospho.columns:
        # Repeats phosphosite as columns for pearson correlation
        phospho_col = phospho.loc[
            phospho.index.intersection(ct_scores.index),
            phosphosite
        ].dropna()
        res = pearsonr(
            pd.concat(
                [phospho_col] * ct_scores.shape[1],
                axis=1
            ),
            ct_scores.loc[phospho_col.index, :]
        )
        rhos.loc[:, phosphosite] = res.statistic
        p_values.loc[:, phosphosite] = res.pvalue

    # False discovery control for p-values
    # Might be a little stringent, as this controls for all p-values across
    # phosphosites and cell types, but correlation is very strong
    p_values.loc[:] = false_discovery_control(p_values)

    # Setup markers for significance
    annot_df = np.empty(rhos.shape, dtype=np.dtype('U100'))
    annot_df[p_values.values < 0.05] = "*"
    annot_df[p_values.values < 0.01] = "**"
    annot_df[p_values.values < 0.001] = "***"

    # Plot heatmap of correlations between phosphosites and cell type scores
    ax = axes[0]
    sns.heatmap(
        rhos,
        center=0,
        cmap="coolwarm",
        annot=annot_df,
        fmt="s",
        cbar_kws={"label": "Pearson Correlation"},
        ax=ax
    )

    # Label axes, add tick labels
    ax.set(
        xticks=np.arange(0.5, rhos.shape[1], 1),
        yticks=np.arange(0.5, rhos.shape[0], 1),
        xlabel="Phosphosite (PRKACA Substrates)",
        ylabel="Cell Types Score\n(Bottomly et al)"
    )
    ax.set_xticklabels(
        rhos.columns,
        rotation=45,
        ha="right",
        ma="right",
        va="top"
    )
    ax.set_yticklabels(
        rhos.index,
        ha="right",
        ma="right",
        va="center",
        rotation=0
    )

    ############################################################################
    ### Cheng et al Aging Score Comparisons
    ############################################################################

    # Load HAS scores
    has_scores = get_has_scores()

    # Setup correlation storage for HAS comparisons
    has_rhos = pd.DataFrame(
        0,
        dtype=np.float64,
        index=phospho.columns,
        columns=["Weighted HAS", "Unweighted HAS"]
    )
    has_p = has_rhos.copy(deep=True)

    # Correlate HAS scores with each phosphosite
    for phosphosite in phospho.columns:
        phospho_col = phospho.loc[
            phospho.index.intersection(ct_scores.index),
            phosphosite
        ].dropna()
        res = pearsonr(
            pd.concat(
                [phospho_col] * 2,
                axis=1
            ),
            has_scores.loc[phospho_col.index, :]
        )
        has_rhos.loc[phosphosite, :] = res.statistic
        has_p.loc[phosphosite, :] = res.pvalue

    # False discovery control for p-values
    has_p.loc[:] = false_discovery_control(has_p)
    has_rhos = has_rhos.T
    has_p = has_p.T

    # Setup markers for significance
    has_annot = np.empty(has_rhos.shape, dtype=np.dtype('U100'))
    has_annot[has_p.values < 0.05] = "*"
    has_annot[has_p.values < 0.01] = "**"
    has_annot[has_p.values < 0.001] = "***"

    # Plot heatmap of correlations between phosphosites and HAS scores
    ax = axes[1]
    sns.heatmap(
        has_rhos,
        center=0,
        cmap="coolwarm",
        annot=has_annot,
        fmt="s",
        cbar_kws={"label": "Pearson Correlation"},
        ax=ax
    )
    ax.set(
        xticks=np.arange(0.5, has_rhos.shape[1], 1),
        yticks=np.arange(0.5, has_rhos.shape[0], 1),
        xlabel="Phosphosite (PRKACA Substrates)",
        ylabel="HAS Score\n(Cheng et al)"
    )
    ax.set_xticklabels(
        has_rhos.columns,
        rotation=45,
        ha="right",
        ma="right",
        va="top"
    )
    ax.set_yticklabels(
        has_rhos.index,
        ha="right",
        ma="right",
        va="center",
        rotation=0
    )

    return fig


if __name__ == '__main__':
    make_figure()
