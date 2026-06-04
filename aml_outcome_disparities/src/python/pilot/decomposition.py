from os.path import abspath, dirname, join
from typing import Iterable

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from sklearn.cross_decomposition import PLSRegression
from sklearn.impute import KNNImputer
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from pilot.data_import import (import_acetyl, import_global, import_rna,
                               syn_login)

REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def run_plsr(
    data: pd.DataFrame, labels: pd.Series, n_comp: int | Iterable[int] = 2
) -> tuple[pd.Series, list[PLSRegression]]:
    """
    Runs PLSR.

    Args:
        data (pd.DataFrame): Dataframe containing measurements.
        labels (Iterable): Regression target.
        n_comp (int, default:2): Number of components to use. If an array is
             provided, each number of components will be tested.

    Returns:
        pd.Series: Cross-validation accuracy across ranks.
        sklearn.PLSR: Fitted PLSR model.
    """
    if isinstance(n_comp, int):
        n_comp = [n_comp]

    scaler = StandardScaler()
    skf = StratifiedKFold(n_splits=10)
    predictions = pd.DataFrame(0, dtype=float, index=data.index, columns=n_comp)
    if data.isna().any().any():
        knn = KNNImputer()

    for train_index, test_index in skf.split(data, labels):
        train_index = data.index[train_index]
        test_index = data.index[test_index]

        train_data = data.loc[train_index, :]
        test_data = data.loc[test_index, :]
        train_labels = labels.loc[train_index]

        if data.isna().any().any():
            train_data[:] = knn.fit_transform(train_data)
            test_data[:] = knn.transform(test_data)

        train_data = scaler.fit_transform(train_data)
        test_data = scaler.transform(test_data)

        for rank in n_comp:
            plsr = PLSRegression(n_components=rank)
            plsr.fit(train_data, train_labels)
            predictions.loc[test_index, rank] = plsr.predict(test_data)

    accuracies = pd.Series(0, dtype=float, index=n_comp)
    plsr_models = []
    if data.isna().any().any():
        data[:] = knn.fit_transform(data)

    for rank in n_comp:
        accuracies.loc[rank] = accuracy_score(
            labels, predictions.loc[:, rank].round()
        )
        plsr = PLSRegression(n_components=rank)
        plsr.fit(data, labels)
        plsr_models.append(plsr)

    return accuracies, plsr_models


def run_acetyl_plsr() -> tuple[
    dict[str, pd.DataFrame],
    dict[str, pd.DataFrame],
    dict[str, pd.DataFrame]
]:
    """
    Runs stacked PLSR approach to interpret acetyl.

    Args:
        None.

    Returns:
        dict[str, pd.DataFrame]: Acetyl/Global components, following deflation.
        dict[str, pd.DataFrame]: Global reconstructed from RNA, RNA + acetyl.
        dict[str, pd.DataFrame]: Variance explained in deflation, residual
            steps.
    """
    # Import data
    syn = syn_login()
    acetyl = pd.concat(import_acetyl(syn))
    rna = pd.concat(import_rna(syn))
    prot = pd.concat(import_global(syn))

    # Log-scale RNA
    rna.loc[:] = np.log10(rna + 1)
    rna = rna.loc[
        :,
        ~rna.columns.duplicated()
    ]

    # Trim to genes in both global and RNA measurements
    prot = prot.loc[:, rna.columns.intersection(prot.columns)]
    rna = rna.loc[:, prot.columns.intersection(rna.columns)]

    # Trim to shared RNA/proteomic samples
    rna = rna.loc[prot.index.intersection(rna.index), :]
    external_prot = prot.drop(rna.index)
    train_prot = prot.loc[rna.index, :]

    # Scale measurements
    # We use StandardScaler to restore to native measurement scale post-PLSR
    rna_scaler = StandardScaler()
    prot_scaler = StandardScaler()
    acetyl_scaler = StandardScaler()
    rna.loc[:] = rna_scaler.fit_transform(rna)
    prot.loc[:] = prot_scaler.fit_transform(prot)
    acetyl.loc[:] = acetyl_scaler.fit_transform(acetyl)

    # Convert to R
    with localconverter(ro.default_converter + pandas2ri.converter):
        rna_r = ro.conversion.py2rpy(rna)
        prot_r = ro.conversion.py2rpy(train_prot)
        extra_prot_r = ro.conversion.py2rpy(external_prot)

    # Source R PLSR function
    r_source = ro.r["source"]
    r_source(
        join(
            REPO_PATH,
            "r",
            "analysis_tools.R"
        )
    )
    run_plsr = ro.globalenv["run_plsr"]

    # Run PLSR between acetyl and RNA
    res = run_plsr(
        prot_r,
        rna_r,
        4,
        extra_prot_r
    )

    # Get PLSR components
    with localconverter(ro.default_converter + pandas2ri.converter):
        rna_p_scores, _, rna_p_loadings, _, deflate_r2x = ro.conversion.rpy2py(
            res
        )

    # Deflate variance covered by transcriptomics
    rna_reconstructed = rna_p_scores @ rna_p_loadings.T
    deflated = prot - rna_reconstructed

    # Scale deflated measurements
    deflated_scaler = StandardScaler()
    deflated.loc[:] = deflated_scaler.fit_transform(deflated)

    # Convert deflated, acetyl measurements to R
    with localconverter(ro.default_converter + pandas2ri.converter):
        deflated_prot = ro.conversion.py2rpy(deflated)
        acetyl_r = ro.conversion.py2rpy(
            acetyl.loc[deflated.index, :]
        )

    # Run PLSR on deflated measurements
    deflated_res = run_plsr(
        deflated_prot,
        acetyl_r,
        4
    )

    # Convert PLSR results to Python
    with localconverter(ro.default_converter + pandas2ri.converter):
        p_scores, a_scores, p_loadings, a_loadings, r2x = ro.conversion.rpy2py(
            deflated_res
        )

    # Reconstruct global from RNA, RNA + acetyl
    acetyl_reconstructed = p_scores @ p_loadings.T
    acetyl_reconstructed.loc[:] = prot_scaler.inverse_transform(
        deflated_scaler.inverse_transform(
            acetyl_reconstructed
        ) + rna_reconstructed
    )
    rna_reconstructed.loc[:] = prot_scaler.inverse_transform(rna_reconstructed)

    # Store results
    components = {
        "Acetyl Scores": a_scores,
        "Proteomic Scores": p_scores,
        "Acetyl Loadings": a_loadings,
        "Proteomic Loadings": p_loadings
    }
    reconstructed = {
        "RNA": rna_reconstructed,
        "RNA + Acetyl": acetyl_reconstructed
    }
    r2x = {
        "Deflation": deflate_r2x,
        "Residual": r2x
    }

    return components, reconstructed, r2x
