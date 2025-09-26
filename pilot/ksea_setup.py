import gseapy as gp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pilot.data_import import (get_ksea_table, import_phospho,
                               import_pilot_meta, import_ptrc_meta)
from pilot.gene_analysis import enrichment_analysis, volcano_plot


def main():
    ptrc = import_ptrc_meta()
    pilot = import_pilot_meta()
    meta = pd.concat([ptrc, pilot])
    meta.loc[:, "batch"] = np.repeat(["PTRC", "Pilot"], [ptrc.shape[0], pilot.shape[0]])

    ptrc, pilot = import_phospho(corrected=True)
    data = pd.concat([ptrc, pilot])

    data = data.loc[data.index.intersection(meta.index), :]
    meta = meta.loc[data.index, :]

    black_patients = data.loc[meta.loc[:, "Race_label"] == "Black", :]
    white_patients = data.loc[meta.loc[:, "Race_label"] == "White", :]

    black_ksea = get_ksea_table()
    black_patients = 2**black_patients
    white_patients = 2**white_patients
    black_ksea.loc[:, "FC"] = black_patients.mean(axis=0) / white_patients.mean(axis=0)

    white_ksea = black_ksea.copy(deep=True)
    white_ksea.loc[:, "FC"] = 1 / white_ksea.loc[:, "FC"]

    black_ksea.to_csv("black_ksea.csv", index=False)
    white_ksea.to_csv("white_ksea.csv", index=False)


if __name__ == "__main__":
    main()
