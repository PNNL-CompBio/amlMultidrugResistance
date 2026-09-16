"""
Clusters cell-type deconvolution and plots cell composition for Pilot study.
"""

import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from leidenalg import ModularityVertexPartition, find_partition
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph

from pilot.data_import import import_deconvolution
from pilot.figure_setup import get_setup
from pilot.utils import reorder_table

COLORS = pd.Series(
    {
        "LSPC-Primed": "salmon",
        "LSPC-Quiescent": "tomato",
        "LSPC-Cycle": "red",
        "ProMono-like": "cornflowerblue",
        "Mono-like": "royalblue",
        "Monocyte": "blue",
        "cDC": "darkblue",
        "cDC-like": "navy",
        "GMP-like": "purple",
        "NK": "whitesmoke",
        "Plasma": "lightgrey",
        "T": "silver",
        "B": "grey",
        "CTL": "black",
    }
)


def main():
    deconv = import_deconvolution()
    deconv = deconv.loc[:, COLORS.index]

    pca = PCA(n_components=2)

    components = pd.DataFrame(
        pca.fit_transform(deconv), index=deconv.index, columns=[1, 2]
    )

    mono_proportion = deconv.loc[
        :, ["Mono-like", "Monocyte", "ProMono-like", "cDC", "cDC-like"]
    ].sum(axis=1)
    gmp_proportion = deconv.loc[:, "GMP-like"]
    lspc_proportion = deconv.loc[:, deconv.columns.str.contains("LSPC")].sum(axis=1)

    names = ["Monocyte", "GMP", "LSPC"]
    proportions = [mono_proportion, gmp_proportion, lspc_proportion]
    fig, axes = get_setup(1, 3, {"figsize": (9, 3)})

    for ax, proportion, name in zip(axes, proportions, names):
        ax.scatter(
            components.loc[:, 1],
            components.loc[:, 2],
            s=5,
            c=proportion,
            cmap="coolwarm",
            vmin=0,
            vmax=1,
        )
        ax.set_title(name)

    plt.show()

    knn_graph = kneighbors_graph(
        components, n_neighbors=20, mode="distance", include_self=False
    )
    graph = ig.Graph.Weighted_Adjacency(knn_graph)
    partition = find_partition(graph, ModularityVertexPartition)
    cluster_labels = np.array(partition.membership) + 1

    fig, ax = get_setup(1, 1, {"figsize": (3, 3)})
    for cluster in np.unique(cluster_labels):
        ax.scatter(
            components.loc[cluster_labels == cluster, 1],
            components.loc[cluster_labels == cluster, 2],
            s=5,
            label=f"Cluster {cluster}",
        )

    ax.set_xlabel("PCA 1")
    ax.set_ylabel("PCA 2")
    ax.legend()
    plt.show()

    fig, axes = get_setup(2, 3, {"figsize": (5, 4), "width_ratios": [2, 2, 1]})

    gs = axes[0, 2].get_gridspec()
    for ax in axes[:, 2]:
        ax.remove()

    legend_ax = fig.add_subplot(gs[:, 2])
    legend_ax.set_frame_on(False)
    legend_ax.set_xticks([])
    legend_ax.set_yticks([])

    axes = axes[:, :2].flatten()

    for ax, cluster in zip(axes, np.unique(cluster_labels)):
        cluster_cells = reorder_table(deconv.loc[cluster_labels == cluster, :])
        for cell_type in COLORS.index:
            ax.bar(
                np.arange(cluster_cells.shape[0]),
                cluster_cells.loc[:, cell_type],
                bottom=cluster_cells.loc[:, :cell_type].iloc[:, :-1].sum(axis=1),
                color=COLORS.loc[cell_type],
                width=1,
                label=cell_type,
            )

        ax.set_xlim([-0.5, cluster_cells.shape[0] - 0.5])
        ax.set_ylim([0, 1])
        ax.set_xticks([])
        ax.set_xlabel(f"Patients (n={cluster_cells.shape[0]})")
        ax.set_ylabel("Cell Type\nProportion")
        ax.set_title(f"Cluster {cluster}")

        if ax == axes[-1]:
            handles, labels = ax.get_legend_handles_labels()

    legend_ax.legend(handles, labels, loc="center left")

    plt.show()


if __name__ == "__main__":
    main()
