from typing import Iterable

import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.impute import KNNImputer
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler


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
        accuracies.loc[rank] = accuracy_score(labels, predictions.loc[:, rank].round())
        plsr = PLSRegression(n_components=rank)
        plsr.fit(data, labels)
        plsr_models.append(plsr)

    return accuracies, plsr_models
