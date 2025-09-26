from os import PathLike
from os.path import abspath, dirname, join, splitext
from typing import Iterable

import datashader as ds
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

REPO_PATH = dirname(dirname(abspath(__file__)))


def preranked_enrichment(
    genes: pd.Series,
    output_path: str | PathLike | None = None,
    gene_libraries: Iterable[str] | None = None,
):
    """
    Performs enrichment analysis via Enrichr.

    Args:
        genes (pd.Series): Gene rankings; index is gene names, values are
            ranking.
        output_path (PathLike, optional): Path to output file. Defaults to
            'data/enrichr_results'. Saves both a .png and a .csv.
        gene_libraries (list[str]): Libraries to compare against. Library
            names are found here: https://maayanlab.cloud/Enrichr/#libraries

    Returns:
        None.

    Notes:
        Stores dotplot, enrichment results to output_path.
    """
    if output_path is None:
        output_path = join(REPO_PATH, "data", "enrichr_results")
    if "." in output_path:
        output_path = splitext(output_path)[:-1][0]
    if gene_libraries is None:
        gene_libraries = ["GO_Biological_Process_2025"]

    for library in gene_libraries:
        result = gp.prerank(genes, gene_sets=[library], outdir=None).res2d

        over = result.loc[result.loc[:, "NES"] > 0, :]
        under = result.loc[result.loc[:, "NES"] < 0, :]

        for exp, df in zip(["over", "under"], [over, under]):
            df = df.loc[df.loc[:, "FDR q-val"] < 0.05, :]
            if df.shape[0] == 0:
                continue

            out_file = output_path + "_" + library + f"_{exp}expressed"
            df.to_csv(out_file + ".csv")

            if library == "GO_Biological_Process_2025":
                df.loc[:, "Term"] = df.loc[:, "Term"].str[:-13]

            df.loc[:, "Term"] = (
                df.loc[:, "Term"].str.split("__", expand=True).iloc[:, 1]
            )
            for row in df.index:
                term = df.loc[row, "Term"]
                char_index = 0
                current_length = 0
                while char_index < len(term):
                    if current_length >= 20 and term[char_index] == " ":
                        term = term[:char_index] + "\n" + term[char_index + 1 :]
                        current_length = 0

                    char_index += 1
                    current_length += 1

                df.loc[row, "Term"] = term

            gp.dotplot(
                df,
                column="FDR q-val",
                title=library.replace("_", " "),
                cutoff=0.25,
                top_term=8,
                cmap=plt.cm.viridis,
                ofname=out_file + ".png",
            )


def enrichment_analysis(genes: list[str], output_path: PathLike | str | None = None):
    """
    Performs enrichment analysis via Enrichr.

    Args:
        genes (list[str]): List of gene names.
        output_path (PathLike, optional): Path to output file. Defaults to
            'data/enrichr_results'. Saves both a .png and a .csv.

    Returns:
        None.

    Notes:
        Stores dotplot, enrichment results to output_path.
    """
    if output_path is None:
        output_path = join(REPO_PATH, "data", "enrichr_results")
    if "." in output_path:
        output_path = splitext(output_path)[:-1][0]

    result = gp.enrichr(
        genes, gene_sets=["GO_Biological_Process_2025"], outdir=None
    ).res2d

    result.loc[:, "-log(p)"] = -np.log10(result.loc[:, "Adjusted P-value"])
    result.to_csv(output_path + ".csv")

    gp.dotplot(
        result,
        title="GO_Biological_Process_2025",
        cmap=plt.cm.viridis,
        ofname=output_path + "_dotplot.png",
    )


def calculate_fc(dataset_1: pd.DataFrame, dataset_2: pd.DataFrame, log: bool = True):
    """
    Runs t-test and calculates fold-change in measurements between datasets.

    Args:
        dataset_1 (pd.DataFrame): First measurement set; over-expressed
            measurements will be positive.
        dataset_2 (pd.DataFrame): Second measurement set; over-expressed
            measurements will be negative.
        log (bool, default:True): Defines if measurements are log-transformed.

    Returns:
        pd.Series: Fold-change in measurement differences between datasets.
        pd.Series: Corrected p-values (FDR) for measurements between datasets.
    """
    t_result = ttest_ind(dataset_1, dataset_2, axis=0)
    corrected_p = pd.Series(
        multipletests(np.nan_to_num(t_result.pvalue), method="fdr_bh")[1],
        index=dataset_1.columns,
    )

    if log:
        fc = dataset_1.mean(axis=0) - dataset_2.mean(axis=0)
    else:
        fc = dataset_1.mean(axis=0) / dataset_2.mean(axis=0)

    return fc, corrected_p


def volcano_plot(
    fc: pd.Series,
    p_values: pd.Series,
    ax: plt.Axes | None = None,
    fc_min: float = 0.5,
    p_max: float = 0.05,
    x_max: int = 2,
    y_max: int = 20,
):
    """
    Creates volcano plot comparing datasets.

    Args:
        fc (pd.Series): Fold-change in measurements.
        p_values (pd.Series): Corrected p-values.
        ax (plt.Axes): Axes to plot on; providing None skips plotting.
        fc_min (float): Minimum fold-change for significance.
        p_max (float): Maximum p-value for significance.
        x_max (int): Maximum FC for plot.
        y_max (int): Maximum -log10(p) for plot.

    Returns:
        list[str]: Over-expressed genes for dataset_1.
        list[str]: Over-expressed genes for dataset_2.
        plt.Axes: Axes containing volcano plot.
    """
    colors = pd.Series("lightgrey", index=p_values.index)
    colors.loc[np.logical_and(fc.values < -fc_min, p_values.values < p_max)] = "blue"
    colors.loc[np.logical_and(fc.values > fc_min, p_values.values < p_max)] = "red"

    matrix = pd.concat([fc, -np.log10(p_values)], axis=1)
    matrix.columns = matrix.columns.astype(str)
    matrix.loc[:, "label"] = pd.Categorical(colors)
    scatter_colors = {color: color for color in colors.unique()}

    cvs = ds.Canvas(
        plot_width=200, plot_height=200, x_range=(-x_max, x_max), y_range=(0, y_max)
    )
    agg = cvs.points(
        matrix, matrix.columns[0], matrix.columns[1], agg=ds.count_cat("label")
    )
    result = ds.tf.shade(agg, color_key=scatter_colors, how="eq_hist", min_alpha=255)
    # result = tf.dynspread(result, threshold=0.95, max_px=10)
    result = ds.tf.set_background(result, "white")
    img_rev = result.data[::-1]
    mpl_img = np.dstack(
        [img_rev & 0x0000FF, (img_rev & 0x00FF00) >> 8, (img_rev & 0xFF0000) >> 16]
    )

    ax.imshow(mpl_img)
    ax.set_xlabel("Log-fold Change")
    ax.set_ylabel("-log(p-value)")

    ax.set_xticks(np.arange(0, 250, 50))
    ax.set_yticks(np.arange(0, 250, 50))
    ax.set_xticklabels(np.linspace(-x_max, x_max, 5))
    ax.set_yticklabels(np.linspace(y_max, 0, 5))

    over_exp_1 = list(colors.loc[colors == "red"].index)
    over_exp_2 = list(colors.loc[colors == "blue"].index)

    return over_exp_1, over_exp_2, ax
