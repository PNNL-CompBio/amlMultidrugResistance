from os import PathLike
from os.path import abspath, dirname, join, splitext
from typing import Iterable

import datashader as ds
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ttest_ind
from sklearn.preprocessing import scale
from statsmodels.multivariate.pca import PCA

REPO_PATH = dirname(dirname(abspath(__file__)))


def preranked_enrichment(
    genes: pd.Series,
    output_path: str | PathLike | None = None,
    gene_libraries: Iterable[str] | None = None,
) -> None:
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
    # Setup output directory for results
    if output_path is None:
        output_path = join(REPO_PATH, "../data", "enrichr_results")
    if "." in output_path:
        output_path = splitext(output_path)[:-1][0]
    if gene_libraries is None:
        gene_libraries = ["GO_Biological_Process_2025"]

    # Iterate through queued libraries
    for library in gene_libraries:
        # Runs pre-ranked enrichment
        result = gp.prerank(genes, gene_sets=[library], outdir=None).res2d

        # Separate results for over-/under-expressed gene sets
        over = result.loc[result.loc[:, "NES"] > 0, :]
        under = result.loc[result.loc[:, "NES"] < 0, :]
        for exp, df in zip(["over", "under"], [over, under]):
            # Trim to only those with FDR < 0.05
            df = df.loc[df.loc[:, "FDR q-val"] < 0.05, :]
            if df.shape[0] == 0:
                continue

            # Export results to .csv
            out_file = output_path + "_" + library + f"_{exp}expressed"
            df.to_csv(out_file + ".csv")

            # Trims term tags if using GO Biological Process
            if library == "GO_Biological_Process_2025":
                df.loc[:, "Term"] = df.loc[:, "Term"].str[:-13]

            # Drops prefix (name of library)
            df.loc[:, "Term"] = (
                df.loc[:, "Term"].str.split("__", expand=True).iloc[:, 1]
            )
            # Insert breaks if term names are too long
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

            # Plot dotplot
            gp.dotplot(
                df,
                column="FDR q-val",
                title=library.replace("_", " "),
                cutoff=0.25,
                top_term=8,
                cmap=plt.cm.viridis,
                ofname=out_file + ".png",
            )


def enrichment_analysis(
    genes: list[str], output_path: PathLike | str | None = None
) -> None:
    """
    Performs enrichment analysis via Enrichr.

    Args:
        genes (list[str]): List of gene names.
        output_path (PathLike | str | None, default: None): Path to output file.
        Defaults to 'data/enrichr_results' if None. Saves both a .png and a
        .csv.

    Returns:
        None.

    Notes:
        Stores dotplot, enrichment results to output_path.
    """
    # Defines output directory for results
    if output_path is None:
        output_path = join(REPO_PATH, "../data", "enrichr_results")
    if "." in output_path:
        output_path = splitext(output_path)[:-1][0]

    # Runs enrichment via Enrichr
    result = gp.enrichr(
        genes, gene_sets=["GO_Biological_Process_2025"], outdir=None
    ).res2d

    # Adds -log(p) column to results
    result.loc[:, "-log(p)"] = -np.log10(result.loc[:, "Adjusted P-value"])
    result.to_csv(output_path + ".csv")

    # Plots dotplot
    gp.dotplot(
        result,
        title="GO_Biological_Process_2025",
        cmap=plt.cm.viridis,
        ofname=output_path + "_dotplot.png",
    )


def calculate_fc(
    dataset_1: pd.DataFrame, dataset_2: pd.DataFrame, log: bool = True
) -> tuple[pd.Series, pd.Series]:
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
    # Runs t-test
    t_result = ttest_ind(dataset_1, dataset_2, axis=0, nan_policy="omit")
    corrected_p = pd.Series(
        t_result.pvalue, index=dataset_1.columns, dtype=np.float64
    )

    # Drops results where t-test couldn't be performed (too few samples) and
    # corrects via Benjamini-Hochberg for FDR
    missing_p = corrected_p.isna()
    corrected_p = corrected_p.loc[~missing_p]
    corrected_p.loc[:] = false_discovery_control(
        corrected_p.values, method="bh"
    )

    # Derive fold-change
    if log:
        fc = dataset_1.mean(axis=0) - dataset_2.mean(axis=0)
    else:
        fc = dataset_1.mean(axis=0) / dataset_2.mean(axis=0)

    fc = fc.loc[~missing_p]
    return fc, corrected_p


def cell_type_scores(data: pd.DataFrame) -> pd.DataFrame:
    """
    Identifies AML cell types in samples using cell types defined in van Galen et al.
    (10.1016/j.cell.2019.01.031)

    Args:
        data (pd.DataFrame): Expression data.

    Returns:
        pd.DataFrame: AML cell type scores for each sample.
    """
    # Loads van Galen genes
    vg_genes = pd.read_csv(join(REPO_PATH, "data", "van_Galen_genes.csv"))
    vg_scores = pd.DataFrame(
        0, dtype=float, index=data.index, columns=vg_genes.columns
    )

    # Calculates PC1 scores for each set of cell-type genes
    for cell_type in vg_genes.columns:
        # Trims to genes included in both van Galen list and in dataset
        cell_type_data = data.loc[
            :,
            data.columns.intersection(vg_genes.loc[:, cell_type]),
        ]
        cell_type_data.loc[:] = scale(cell_type_data, axis=0)
        pca = PCA(
            cell_type_data,
            ncomp=1,
            missing="fill-em",
            demean=False,
            standardize=False,
            method="nipals",
        )
        vg_scores.loc[:, cell_type] = pca.scores.iloc[:, 0].values

    return vg_scores


def volcano_plot(
    fc: pd.Series,
    p_values: pd.Series,
    ax: plt.Axes | None = None,
    fc_min: float = 1,
    p_max: float = 0.05,
    x_max: int = 2,
    y_max: int = 20,
    left_color: str = "red",
    right_color: str = "purple",
) -> tuple[list, list, plt.Axes]:
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
        left_color (str): Color for significantly low FCs.
        right_color (str): Color for significantly high FCs.

    Returns:
        list[str]: Over-expressed genes for dataset_1.
        list[str]: Over-expressed genes for dataset_2.
        plt.Axes: Axes containing volcano plot.
    """
    # Sets points with high FC and significant p-values to be colored
    colors = pd.Series("lightgrey", index=p_values.index)
    colors.loc[np.logical_and(fc.values < -fc_min, p_values.values < p_max)] = (
        left_color
    )
    colors.loc[np.logical_and(fc.values > fc_min, p_values.values < p_max)] = (
        right_color
    )

    # Concatenates FC, p-values for datashader
    matrix = pd.concat([fc, -np.log10(p_values)], axis=1)
    matrix.columns = matrix.columns.astype(str)
    matrix.loc[:, "label"] = pd.Categorical(colors)
    scatter_colors = {color: color for color in colors.unique()}

    # Plots via datashader
    # datashader is useful for datasets with lots of variables (like
    # transcriptomics) where points are likely to overlap
    cvs = ds.Canvas(
        plot_width=200,
        plot_height=200,
        x_range=(-x_max, x_max),
        y_range=(0, y_max),
    )
    agg = cvs.points(
        matrix, matrix.columns[0], matrix.columns[1], agg=ds.count_cat("label")
    )
    result = ds.tf.shade(
        agg, color_key=scatter_colors, how="eq_hist", min_alpha=255
    )
    result = ds.tf.set_background(result, "white")
    img_rev = result.data[::-1]
    mpl_img = np.dstack(
        [
            img_rev & 0x0000FF,
            (img_rev & 0x00FF00) >> 8,
            (img_rev & 0xFF0000) >> 16,
        ]
    )

    ax.imshow(mpl_img)
    ax.set_xlabel("Log-fold Change")
    ax.set_ylabel("-log(p-value)")

    # Set up ticks--needs some transformation given how datashader handles
    # axes
    ax.set_xticks(np.arange(0, 250, 50))
    ax.set_yticks(np.arange(0, 250, 50))
    ax.set_xticklabels(np.linspace(-x_max, x_max, 5))
    ax.set_yticklabels(np.linspace(y_max, 0, 5))

    over_exp_1 = list(colors.loc[colors == right_color].index)
    over_exp_2 = list(colors.loc[colors == left_color].index)

    return over_exp_1, over_exp_2, ax
