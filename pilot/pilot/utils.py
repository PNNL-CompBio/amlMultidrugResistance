import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch


def nan_normalize(data: pd.DataFrame, axis: int = 1) -> pd.DataFrame:
    """
    Normalizes samples, ignoring NaN values.

    Args:
        data (pd.DataFrame): Data to normalize.
        axis (int, default: 1): Axis to normalize.

    Returns:
        pd.DataFrame: Normalized data.
    """
    data = data.T / np.sqrt(np.nansum(data**2, axis=axis))
    return data.T


def reorder_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Reorder a DataFrame's rows using hierarchical clustering.

    Args:
        df (pandas.DataFrame): data to be clustered; rows are treated as samples
            to be clustered.

    Returns:
        df (pandas.DataFrame): data with rows reordered via hierarchical
            clustering.
    """
    y = sch.linkage(df.to_numpy(), method="centroid")
    index = sch.dendrogram(y, orientation="right", no_plot=True)["leaves"]
    return df.iloc[index]
