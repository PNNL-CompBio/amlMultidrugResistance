"""
Plots comparison of LCP1s and inflammation score.
"""

from os.path import abspath, dirname, join

import matplotlib.pyplot as plt
import pandas as pd

from pilot.data_import import import_global, import_meta, import_phospho
from pilot.figure_setup import get_setup, run_ols

REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def main():
    meta = import_meta()
    ptrc, pilot = import_phospho(corrected=True)
    phospho = pd.concat([ptrc, pilot])

    ptrc, pilot = import_global()
    global_prot = pd.concat([ptrc, pilot])
    global_prot = global_prot.loc[~global_prot.index.duplicated(), :]
    iscore_coef = pd.read_csv(join(REPO_PATH, "data", "iscore_coef.csv"), index_col=0)
    shared_proteins = global_prot.columns.intersection(iscore_coef.index)
    iscore_coef = iscore_coef.loc[shared_proteins, "beta mean"]

    patients = phospho.index.intersection(meta.index)
    phospho = phospho.loc[patients, :]
    global_prot = global_prot.loc[patients, :]
    i_scores = global_prot.loc[:, iscore_coef.index].fillna(0) @ iscore_coef

    fig, ax = get_setup(1, 1, {"figsize": (4, 4)})

    ax.scatter(phospho.loc[:, "LCP1-S5s"], i_scores, s=15)
    ols = run_ols(
        phospho.loc[:, "LCP1-S5s"].dropna(),
        i_scores.loc[~phospho.loc[:, "LCP1-S5s"].isna()],
        ax,
    )

    ax.text(
        0.99,
        0.01,
        s=f"p-value: {round(ols.pvalues.loc['LCP1-S5s'], 5)}\n"
        f"R2: {round(ols.rsquared_adj, 3)}",
        ha="right",
        ma="right",
        va="bottom",
        fontsize=8,
        transform=ax.transAxes,
    )
    ax.set_xlabel("LCP1-S5s")
    ax.set_ylabel("iScore")

    plt.show()


if __name__ == "__main__":
    main()
