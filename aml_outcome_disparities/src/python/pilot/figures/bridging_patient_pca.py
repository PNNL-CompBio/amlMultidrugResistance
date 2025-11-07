"""
Overlays PCA with changes in bridging patient samples.
"""

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer
from sklearn.preprocessing import scale

from pilot.data_import import (import_global, import_meta, import_phospho,
                               syn_login)
from pilot.figure_setup import confidence_ellipse, get_setup
from pilot.utils import nan_normalize

BATCH_COLORS = pd.Series(
    ["tab:red", "tab:purple", "black"], index=["BeatAML", "Pilot", "Bridging"]
)


def main():
    syn = syn_login()
    ptrc_phospho, pilot_phospho = import_phospho(syn, corrected=True)
    ptrc_global, pilot_global = import_global(syn)
    meta = import_meta(syn)

    bridge_patients = meta.loc[meta.index.str.contains("Bridge"), :].index
    beataml_ids = bridge_patients.str[:-7]
    meta.loc[bridge_patients, "study"] = "Bridging"

    merged_phospho = pd.concat([ptrc_phospho, pilot_phospho])
    merged_global = pd.concat([ptrc_global, pilot_global])
    knn = KNNImputer()
    merged_phospho.loc[:] = knn.fit_transform(merged_phospho)
    merged_global.loc[:] = knn.fit_transform(merged_global)

    names = ["Phosphoproteomics", "Global Proteomics"]
    datasets = [merged_phospho, merged_global]
    factor_set = []
    for data in datasets:
        _data = nan_normalize(data, axis=1)
        _data = scale(_data, axis=0)
        pca = PCA(n_components=2)
        factor_set.append(
            pd.DataFrame(pca.fit_transform(_data), index=meta.index, columns=[0, 1])
        )

    fig, axes = get_setup(
        1,
        len(factor_set),
        {"figsize": (4 * len(factor_set), 4)},
    )
    for index, (name, factors) in enumerate(zip(names, factor_set)):
        for group in BATCH_COLORS.index:
            axes[index].scatter(
                factors.loc[meta.loc[:, "study"] == group, 0],
                factors.loc[meta.loc[:, "study"] == group, 1],
                c=BATCH_COLORS.loc[group],
                edgecolors=None,
                linewidths=0,
                label=f"{group} " f"(n={sum(meta.loc[:, "study"] == group)})",
                s=15,
                alpha=0.25,
            )
            confidence_ellipse(
                factors.loc[meta.loc[:, "study"] == group, 0],
                factors.loc[meta.loc[:, "study"] == group, 1],
                axes[index],
                edgecolor=BATCH_COLORS.loc[group],
                linestyle="--",
                alpha=0.9,
                label="_" + group,
            )

        for bridge_patient, beataml_id in zip(bridge_patients, beataml_ids):
            axes[index].scatter(
                factors.loc[beataml_id, 0],
                factors.loc[beataml_id, 1],
                c=BATCH_COLORS.loc["BeatAML"],
                edgecolors="black",
                linewidths=0.5,
                s=15,
            )
            axes[index].scatter(
                factors.loc[bridge_patient, 0],
                factors.loc[bridge_patient, 1],
                c=BATCH_COLORS.loc["Bridging"],
                edgecolors="black",
                linewidths=0.5,
                s=15,
            )
            axes[index].annotate(
                "",
                xytext=(factors.loc[beataml_id, 0], factors.loc[beataml_id, 1]),
                xy=(factors.loc[bridge_patient, 0], factors.loc[bridge_patient, 1]),
                arrowprops=dict(arrowstyle="->"),
            )

        axes[index].legend(fontsize=8)
        axes[index].set_title(name)
        axes[index].set_xlabel("PC 1")
        axes[index].set_ylabel("PC 2")

    plt.show()


if __name__ == "__main__":
    main()
