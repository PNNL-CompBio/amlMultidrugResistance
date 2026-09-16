"""
Creates volcano plots comparing phosphoproteomic measurements across race.
"""

from os.path import abspath, dirname, join

import gseapy as gp
import matplotlib.pyplot as plt
import pandas as pd

from pilot.data_import import import_meta, import_phospho
from pilot.figure_setup import get_setup
from pilot.gene_analysis import calculate_fc, enrichment_analysis, volcano_plot

GO_LIBRARY = gp.parser.get_library("GO_Biological_Process_2025")
GO_TERMS = [
    "Negative Regulation of Immune Response (GO:0050777)",
    "Negative Regulation of Innate Immune Response (GO:0045824)",
    "Negative Regulation of Defense Response (GO:0031348)",
    "Negative Regulation of Cytokine Production (GO:0001818)",
    "Negative Regulation of Lymphocyte Activation (GO:0051250)",
    "Inflammatory Response (GO:0006954)",
    "Positive Regulation of Inflammatory Response (GO:0050729)",
]
REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def main():
    meta = import_meta()
    ptrc, pilot = import_phospho(corrected=True)
    data = pd.concat([ptrc, pilot])

    data = data.loc[data.index.intersection(meta.index), :]
    meta = meta.loc[data.index, :]

    black_patients = data.loc[meta.loc[:, "Race"] == "Black", :]
    white_patients = data.loc[meta.loc[:, "Race"] == "White", :]
    fig, ax = get_setup(1, 1, {"figsize": (4, 4)})

    fc, p_values = calculate_fc(black_patients, white_patients)
    black_over_exp, white_over_exp, ax = volcano_plot(fc, p_values, ax=ax)

    black_over_exp = [gene.split("-")[0] for gene in black_over_exp]
    white_over_exp = [gene.split("-")[0] for gene in white_over_exp]

    enrichment_analysis(black_over_exp, join(REPO_PATH, "output", "black_enrichr"))
    enrichment_analysis(white_over_exp, join(REPO_PATH, "output", "white_enrichr"))

    plt.show()


if __name__ == "__main__":
    main()
