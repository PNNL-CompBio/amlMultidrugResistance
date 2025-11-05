"""
Plots correlations between LCP1 and global proteins.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pilot.data_import import import_global, import_meta, import_phospho
from pilot.figure_setup import get_setup, run_ols

SUBSTRATE = "LCP1-S5s"
PROTEINS = ["FCGR1A", "FCGR2A"]


def main():
    meta = import_meta()
    ptrc, pilot = import_phospho(corrected=True)
    phospho = pd.concat([ptrc, pilot])

    ptrc, pilot = import_global()
    global_prot = pd.concat([ptrc, pilot])
    global_prot = global_prot.loc[~global_prot.index.duplicated(), :]

    patients = phospho.index.intersection(meta.index)
    phospho = phospho.loc[patients, SUBSTRATE].squeeze()
    global_prot = global_prot.loc[patients, PROTEINS]

    fig, axes = get_setup(
        1, global_prot.shape[1], {"figsize": (global_prot.shape[1] * 2, 2)}
    )
    for ax, gene in zip(axes, global_prot.columns):
        non_missing = np.logical_and(~phospho.isna(), ~global_prot.loc[:, gene].isna())
        ax.scatter(phospho, global_prot.loc[:, gene], s=5)
        run_ols(phospho.loc[non_missing], global_prot.loc[non_missing, gene], ax)
        ax.set_xlabel(SUBSTRATE)
        ax.set_ylabel(gene)

    plt.show()


if __name__ == "__main__":
    main()
