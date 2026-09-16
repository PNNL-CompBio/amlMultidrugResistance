"""
Plots PCA comparison of BeatAML/Pilot study cohort phospho and global samples.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer
from sklearn.preprocessing import scale

from pilot.data_import import (import_global, import_meta, import_phospho,
                               syn_login)
from pilot.figure_setup import confidence_ellipse, get_setup
from pilot.utils import nan_normalize

RACE_COLORS = pd.Series(
    ["tab:blue", "tab:orange", "grey"], index=["Black", "White", "Not Black or White"]
)
BATCH_COLORS = pd.Series(
    ["tab:red", "tab:purple", "black"], index=["BeatAML", "Pilot", "Bridging"]
)


def main():
    syn = syn_login()
    ptrc_phospho, pilot_phospho = import_phospho(syn, corrected=True)
    ptrc_global, pilot_global = import_global(syn)
    meta = import_meta(syn)

    batches = pd.Series(
        np.repeat(
            ["BeatAML", "Pilot"], [ptrc_phospho.shape[0], pilot_phospho.shape[0]]
        ),
        index=ptrc_phospho.index.append(pilot_phospho.index),
    )
    batches.loc[meta.index.str.contains("Bridge")] = "Bridging"

    merged_phospho = pd.concat([ptrc_phospho, pilot_phospho])
    merged_global = pd.concat([ptrc_global, pilot_global])
    knn = KNNImputer()
    merged_phospho.loc[:] = knn.fit_transform(merged_phospho)
    merged_global.loc[:] = knn.fit_transform(merged_global)

    meta = meta.reindex(merged_phospho.index)
    meta.loc[:, "Batch"] = batches
    meta = meta.loc[:, ["Sex", "Race", "Age", "Batch"]]
    meta.loc[~meta.loc[:, "Race"].isin(["White", "Black"]), "Race"] = (
        "Not Black or White"
    )

    names = ["Phosphoproteomics", "Global Proteomics"]
    datasets = [merged_phospho, merged_global]
    factor_set = []
    for data in datasets:
        _data = nan_normalize(data, axis=1)
        _data = scale(_data, axis=0)
        pca = PCA(n_components=2)
        factor_set.append(pca.fit_transform(_data))

    for column in meta.columns:
        fig, axes = get_setup(
            1,
            len(factor_set),
            {"figsize": (4 * len(factor_set), 4)},
        )
        for index, (name, factors) in enumerate(zip(names, factor_set)):
            if column != "Age":
                if column == "Batch":
                    colors = BATCH_COLORS
                else:
                    colors = RACE_COLORS

                for group in colors.index:
                    axes[index].scatter(
                        factors[meta.loc[:, column] == group, 0],
                        factors[meta.loc[:, column] == group, 1],
                        c=colors.loc[group],
                        edgecolors=None,
                        linewidths=0,
                        label=f"{group} " f"(n={sum(meta.loc[:, column] == group)})",
                        s=15,
                    )
                    confidence_ellipse(
                        factors[meta.loc[:, column] == group, 0],
                        factors[meta.loc[:, column] == group, 1],
                        axes[index],
                        edgecolor=colors.loc[group],
                        linestyle="--",
                        alpha=0.9,
                        label="_" + group,
                    )

                axes[index].legend(fontsize=8)

            else:
                colors = meta.loc[:, column]
                scatter = axes[index].scatter(
                    factors[:, 0], factors[:, 1], c=colors, cmap="coolwarm", s=15
                )
                if index == len(axes) - 1:
                    plt.colorbar(scatter, ax=axes[index])

            axes[index].set_title(name)
            axes[index].set_xlabel("PC 1")
            axes[index].set_ylabel("PC 2")

        plt.show()


if __name__ == "__main__":
    main()
