"""
Plots heatmap of enriched metabolites.
"""

from os.path import abspath, dirname, join

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from pilot.figure_setup import get_setup
from pilot.utils import reorder_table

REPO_PATH = abspath(dirname(dirname(dirname(__file__))))


def main():
    rp = pd.read_excel(
        join(REPO_PATH, "data", "beataml_pilot_bvsw_rp_metabolites.xlsx")
    )
    hilic = pd.read_excel(
        join(REPO_PATH, "data", "beataml_pilot_bvsw_hilic_metabolites.xlsx")
    )
    metabolites = pd.concat([rp, hilic])
    metabolites = metabolites.loc[
        metabolites.loc[:, "P_value_A_Black_vs_White"] < 0.05, :
    ]
    metabolites = metabolites.sort_values("Fold_change_Black_vs_White", ascending=True)
    metabolites = metabolites.loc[
        ~metabolites.loc[:, "unformatted_name"].duplicated(keep="last"), :
    ]
    known_metabolites = metabolites.loc[
        ~metabolites["unformatted_name"].str.contains("Unknown"), :
    ]
    known_metabolites.set_index("unformatted_name", inplace=True)
    known_metabolites = known_metabolites.loc[:, ["Mean_Black", "Mean_White"]]
    known_metabolites = reorder_table(known_metabolites).T

    fig, ax = get_setup(1, 1, {"figsize": (4, 1.5)})
    sns.heatmap(
        known_metabolites, ax=ax, cmap="Greys", cbar_kws={"label": "log2\nexpression"}
    )

    ax.set_xticklabels(
        known_metabolites.columns, rotation=45, ha="right", ma="right", va="top"
    )
    ax.set_yticklabels(["Black Patients", "White Patients"], rotation=0)
    ax.set_xlabel("")

    plt.show()


if __name__ == "__main__":
    main()
