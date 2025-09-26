import numpy as np
import pandas as pd

from pilot.data_import import get_ksea_table, import_phospho


def main():
    ksea_table = get_ksea_table()
    phospho_loadings = pd.read_csv("data/diablo_phospho_loadings.csv", index_col=0)
    phospho_scores = pd.read_csv("data/diablo_phospho_scores.csv", index_col=0)
    phospho_loadings.index = phospho_loadings.index.str.replace(".", "-")
    ksea_table = ksea_table.loc[
        ksea_table.index.intersection(phospho_loadings.index), :
    ]
    ptrc, pilot = import_phospho()
    phospho = pd.concat([ptrc, pilot])
    phospho = phospho.loc[phospho_scores.index, :]

    for comp in phospho_loadings.columns:
        phospho_scores = phospho_scores.sort_values(comp, ascending=False)
        fc = pd.Series(
            (
                np.nanmean(phospho.loc[phospho_scores.loc[:, comp] < 0, :], axis=0) ** 2
                - np.nanmean(phospho.loc[phospho_scores.loc[:, comp] > 0, :], axis=0)
                ** 2
            ),
            index=phospho.columns,
        )
        ksea_table.loc[:, "FC"] = fc.loc[ksea_table.index]
        ksea_table.to_csv(f"data/diablo_{comp}.csv", index_label=False)


if __name__ == "__main__":
    main()
