import requests
from os import makedirs, PathLike
from os.path import abspath, dirname, exists, join, splitext
from pathlib import Path
from typing import Iterable

import datashader as ds
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ttest_ind
from sklearn.preprocessing import scale
from statsmodels.multivariate.pca import PCA

from pilot.data_import import import_global, import_rna, import_phospho

REPO_PATH = dirname(dirname(abspath(__file__)))
HAS_SCORES = pd.Series(
    {
        "FGG": 0.173,
        "FGA": 0.144,
        "FGB": 0.141,
        # "CALCRL": 0.091,  # Not measured in Global
        "S100A9": 0.090,
        "THBS1": 0.083,
        "S100A8": 0.077,
        "VNN2": 0.072,
        "FPR1": 0.070,
        "MMRN1": 0.062,
        "CAVIN2": 0.061,
        "RBPMS": 0.060,
        # "PF4V1": 0.059,   # Very low coverage in Global
        "MS4A3": 0.056,
        "PF4": 0.054,
        "VCAN": 0.047,
        "S100A12": 0.041,
        "MEIS1": 0.041,
        # "PSAP": -0.19,
    }
)
LSC_17 = pd.Series(
    {
        "DNMT3B": 0.0874,
        "ZBTB46": -0.0347,
        "NYNRIN": 0.00865,
        "ARHGAP22": -0.0138,
        "LAPTM4B": 0.00582,
        "MMRN1": 0.0258,
        "DPYSL3": 0.0284,
        "ENSG00000226777": 0.0196,  # ENSEMBL for KIAA0125
        "CDK6": -0.0704,
        "CPXM1": -0.0258,
        "SOCS2": 0.0271,
        "SMIM24": -0.0226,
        "EMP1": 0.0146,
        "BEX3": 0.0465,  # Updated NGFRAP1
        "CD34": 0.0338,
        "AKR1C3": -0.0402,
        "ADGRG1": 0.0501  # Updated GPR56
        # "PSAP": -0.19,
    }
)


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

    # Uses default white background
    plt.style.use("default")

    # Sets to output_path to PathLike, if not already
    if not isinstance(output_path, PathLike):
        output_path = Path(output_path)

    makedirs(output_path, exist_ok=True)

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
            out_file = output_path / f"{library}_{exp}expressed"
            df.to_csv(out_file.with_suffix(".csv"))

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
                    if current_length >= 30 and term[char_index] == " ":
                        term = term[:char_index] + "\n" + term[char_index + 1 :]
                        current_length = 0

                    char_index += 1
                    current_length += 1

                df.loc[row, "Term"] = term

            # Catch in case all q-values are very near 0
            if (df.loc[:, "FDR q-val"] == 0).all():
                df.loc[:, "FDR q-val"] += 1E-6

            # Plot dotplot
            gp.dotplot(
                df,
                column="FDR q-val",
                title=library.replace("_", " "),
                cutoff=0.25,
                top_term=8,
                cmap=plt.cm.viridis,
                ofname=out_file.with_suffix(".png"),
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

    # Uses default white background
    plt.style.use("default")

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


def cell_type_scores() -> pd.DataFrame:
    """
    Identifies AML cell types in samples using cell types defined in van
    Galen et al. (doi.org/10.1016/j.cell.2019.01.031) Uses SVD-based scoring
    approach defined in Bottomly et al. (doi.org/10.1016/j.ccell.2022.07.002)

    Args:
        None.

    Returns:
        pd.DataFrame: AML cell type scores for each sample.
    """
    # Load RNA-seq data
    data = pd.concat(import_rna(), axis=0)

    # Loads van Galen genes
    vg_genes = pd.read_excel(
        join(REPO_PATH, "data", "mmc3.xlsx"), header=1, index_col=0
    )

    gmp_like = vg_genes.loc[:, "GMP-like"].iloc[:-2]
    vg_genes = vg_genes.iloc[:-2, -7:-1]
    vg_genes.rename(columns={"GMP-like.1": "GMP-like"}, inplace=True)
    vg_genes.loc[:, "GMP-like"] = gmp_like
    vg_genes.index = vg_genes.index.astype(int)

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


def cell_type_scores_vg(
    tpm: bool = True,
    puram: bool = False
) -> pd.DataFrame:
    """
    Calculates AML subtype scores using approaches outlined in van Galen et al.
    (doi.org/10.1016/j.cell.2019.01.031) and Puram et al.
    (https://doi.org/10.1016/j.cell.2017.10.044).

    Args:
        tpm (bool): TPM measurements prior to score derivation.
        puram (bool): Use binning-based approach in Puram et al. False uses
            approach outlined in van Galen et al. instead.

    Returns:
        pd.DataFrame: AML cell type scores for each RNA-seq sample.
    """
    # Load RNA-seq data
    data = np.log2(
        pd.concat(
            import_rna(tpm=tpm),
            axis=0
        ) + 1
    )

    # Loads van Galen genes
    vg_genes = pd.read_excel(
        join(REPO_PATH, "data", "mmc3.xlsx"), header=1, index_col=0
    )

    # Reformat GMP-like column
    gmp_like = vg_genes.loc[:, "GMP-like"].iloc[:-2]
    vg_genes = vg_genes.iloc[:-2, -7:-1]
    vg_genes.rename(columns={"GMP-like.1": "GMP-like"}, inplace=True)
    vg_genes.loc[:, "GMP-like"] = gmp_like
    vg_genes.index = vg_genes.index.astype(int)

    # Initialize DataFrame for storage
    ct_scores = pd.DataFrame(
        np.nan,
        dtype=float,
        index=data.index,
        columns=vg_genes.columns
    )

    if puram:
        # Setup RNG
        rng = np.random.default_rng(20170215)

        # Define bins for expression
        bins = pd.qcut(data.mean(axis=0), 25, labels=np.arange(25))

        for cell_type in vg_genes.columns:
            # Get signature genes
            ct_data = data.loc[
                :,
                data.columns.intersection(
                    vg_genes.loc[
                        :,
                        cell_type
                    ].values
                )
            ]

            # Adjust each gene with genes in same expression bin
            for gene in ct_data.columns:
                ct_data.loc[:, gene] -= data.loc[
                    :,
                    rng.choice(bins.loc[bins == bins.loc[gene]].index, 100)
                ].mean(axis=1)

            # Store scores
            ct_scores.loc[:, cell_type] = ct_data.mean(axis=1)

    else:
        # Pre-compute means for each gene
        gene_means = data.mean(axis=0)

        for cell_type in vg_genes.columns:
            # Get signature genes
            ct_data = data.loc[
                :,
                data.columns.intersection(
                    vg_genes.loc[
                        :,
                        cell_type
                    ].values
                )
            ]

            # Adjust each gene with genes with similar mean expression
            for gene in ct_data.columns:
                similar = abs(
                    gene_means - gene_means.loc[gene]
                ).nsmallest(n=101).drop(gene).index
                ct_data.loc[:, gene] -= gene_means.loc[similar].mean()

            # Store scores
            ct_scores.loc[:, cell_type] = ct_data.mean(axis=1)

    return ct_scores


def get_has_scores():
    """
    Calculates Hematopoietic Aging Signature (HAS) scores via Cheng et al.
    (doi.org/10.1182/blood.2024027692)

    Args:
        None.

    Returns:
        pd.DataFrame: Weighted & unweighted HAS scores for each sample.
            Unweighted considers arithmetic mean of measurements; weighted
            includes coefficients from Cheng et al.
    """
    data = pd.concat(import_global(), axis=0)
    data = data.loc[:, HAS_SCORES.index]
    data = 2 ** data

    has_scores = pd.DataFrame(
        0,
        dtype=float,
        index=data.index,
        columns=["Weighted HAS", "Unweighted HAS"],
    )

    # Includes HAS coefficients
    has_scores.loc[:, "Weighted HAS"] = (
        data * HAS_SCORES.loc[data.columns]
    ).mean(axis=1)

    # HAS score without coefficients
    # Not in-place, so the Weighted HAS score stays in data
    # PSAP is negatively associated with HAS, so removed
    # has_scores.loc[:, "Unweighted HAS"] = data.drop("PSAP", axis=1).mean(axis=1)
    has_scores.loc[:, "Unweighted HAS"] = data.mean(axis=1)

    # Try a PCA-based scoring
    data.loc[:] = scale(np.log2(data), axis=0)
    pca = PCA(
        data,
        ncomp=1,
        missing="fill-em",
        demean=False,
        standardize=False,
        method="nipals",
    )
    has_scores.loc[:, "PCA HAS"] = pca.scores.iloc[:, 0].values

    return has_scores


def get_prkaca_score() -> pd.Series:
    """
    Calculates PRKACA score using approach described in Bottomly et al.
    (doi.org/10.1016/j.ccell.2022.07.002)

    Args:
        None.

    Returns:
        pd.Series: PRKACA scores for each sample.
    """
    # Loads substrate idea
    substrates = pd.read_csv(
        join(REPO_PATH, "data", "Kinase-Substrate Links.csv")
    )
    substrates = substrates.loc[
        substrates.loc[:, "Kinase.Gene"] == "PRKACA", :
    ].sort_values(by="log2FC", ascending=True)

    # Reformat names to match Synapse data
    substrates.index = (
        substrates.loc[:, "Substrate.Gene"]
        + "-"
        + substrates.loc[:, "Substrate.Mod"]
        + substrates.loc[:, "Substrate.Mod"].str[0].str.lower()
    )

    # Load phosphoproteomic measurements
    phospho = import_phospho()
    phospho = pd.concat(phospho, axis=0)

    # Trim KSEA table to phospho substrates that also appear in phospho
    # measurements. Not all are included since the KSEA splits phosphosites with
    # multiple phosphorylations into individual columns!
    shared = substrates.index.intersection(phospho.columns)
    phospho = phospho.loc[:, shared]

    # Drop phospho measurements with high missingness, trim to patients with
    # all measurements
    phospho = phospho.loc[:, phospho.isna().mean(axis=0) < 0.01]
    phospho.dropna(axis=0, inplace=True)

    # Fit PCA
    phospho.loc[:] = scale(phospho)
    pca = PCA(
        phospho,
        ncomp=1,
        missing="fill-em",
        demean=False,
        standardize=False,
        method="nipals",
    )

    # Store scores for patients
    scores = pd.Series(pca.factors.squeeze(), index=phospho.index)

    return scores


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


def get_iscore():
    """
    Calculates iScores via Lasry et al. (doi.org/10.1038/s43018-022-00480-0)
    from gene expression.

    Args:
        None.

    Returns:
        pd.DataFrame: iScores for each patient.

    Notes:
        Loads iScore coefficients from Lasry et al if not found in python/data.
    """
    if not exists(join(REPO_PATH, "data", "lasry_supplement.xlsx")):
        # Load response from URL with .xlsx
        response = requests.get(
            "https://static-content.springer.com/esm/"
            "art%3A10.1038%2Fs43018-022-00480-0/MediaObjects/"
            "43018_2022_480_MOESM2_ESM.xlsx"
        )

        # Save to data directory
        with open(
            join(REPO_PATH, "data", "lasry_supplement.xlsx"), "wb"
        ) as xlsx_file:
            xlsx_file.write(response.content)

        xlsx_file.close()

    # Load RNA-seq data
    data = pd.concat(
        import_rna(return_symbols=True, batch_correct=True, tpm=True)
    )

    # Load iScore genes, trim to genes in both datasets
    iscore_genes = pd.read_excel(
        join(REPO_PATH, "data", "lasry_supplement.xlsx"),
        sheet_name="Supp11. Adult iScore genes",
        index_col=0,
    )
    iscore_genes = iscore_genes.loc[
        iscore_genes.index.intersection(data.columns), "beta mean"
    ]
    data = np.log2(data.loc[:, iscore_genes.index])

    iscores = data * iscore_genes
    iscores.replace({-np.inf: 0}, inplace=True)
    iscores = iscores.sum(axis=1)

    return iscores
