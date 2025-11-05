"""
Runs PCA to compare variation w/r to batch and other covariates.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer
from sklearn.preprocessing import LabelEncoder, scale

from pilot.data_import import (batch_correct, import_global, import_meta,
                               import_phospho, syn_login)

TAB10 = plt.get_cmap("tab10")
COLORS = [TAB10(i) for i in range(10)]


def main():
    # Import measurements
    syn = syn_login()
    ptrc_phospho, pilot_phospho = import_phospho(syn, corrected=True)
    ptrc_global, pilot_global = import_global(syn)
    meta = import_meta(syn)
    phospho = pd.concat([ptrc_phospho, pilot_phospho])
    global_prot = pd.concat([ptrc_global, pilot_global])

    # Trim meta-data, ensure patient order
    meta = meta.dropna(axis=0, how="all")
    meta = meta.loc[:, ["Race", "study", "Age"]]
    phospho = phospho.loc[meta.index, :]
    global_prot = global_prot.loc[meta.index, :]
    batches = meta.loc[:, "study"]
    batches = list(batches)

    # Impute missing values for PCA
    knn = KNNImputer()
    phospho.loc[:] = knn.fit_transform(phospho)
    global_prot.loc[:] = knn.fit_transform(global_prot)

    # Batch correct via ComBAT-seq
    bridge_corrected = batch_correct(2**phospho.T, batches)
    bridge_corrected[:] = np.log2(bridge_corrected)

    bridge_global = batch_correct(2**global_prot.T, batches)
    bridge_global[:] = np.log2(bridge_global)

    names = ["Plex Phospho", "Bridge Phospho", "Plex Global", "Bridge Global"]
    datasets = [phospho, bridge_corrected, global_prot, bridge_global]

    fig, axes = plt.subplots(
        len(datasets),
        meta.shape[1],
        figsize=(3 * meta.shape[1], 2 * len(datasets)),
        constrained_layout=True,
    )

    for row, (name, data) in enumerate(zip(names, datasets)):
        _data = scale(data, axis=1)
        _data = scale(_data, axis=0)
        pca = PCA(n_components=2)
        factors = pca.fit_transform(_data)
        encoder = LabelEncoder()
        for index, column in enumerate(meta.columns):
            if column != "Age":
                encoder.fit(meta.loc[:, column])
                groups = encoder.classes_
                for group, color in zip(groups, COLORS):
                    axes[row, index].scatter(
                        factors[meta.loc[:, column] == group, 0],
                        factors[meta.loc[:, column] == group, 1],
                        c=color,
                        s=5,
                    )
                    if row == len(datasets) - 1:
                        axes[row, index].set_xlabel("PC 1")
                    if index == 0:
                        axes[row, index].set_ylabel("PC 2")

            else:
                colors = meta.loc[:, column] / meta.loc[:, column].max()
                axes[row, index].scatter(
                    factors[:, 0], factors[:, 1], c=colors, cmap="coolwarm", s=5
                )
                axes[-1, index].set_xlabel("PC 1")

            axes[0, index].set_title(column)

    plt.show()


if __name__ == "__main__":
    main()
