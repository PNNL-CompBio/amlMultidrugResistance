"""Plots figure 1a: PCA for batch comparison."""

import sys
from os.path import abspath, dirname, join

sys.path.append(
    join(dirname(dirname(dirname(abspath(__file__)))), "src", "python")
)

import pandas as pd
from sklearn.preprocessing import scale
from statsmodels.multivariate.pca import PCA

from pilot.data_import import (
    import_acetyl,
    import_global,
    import_lipids,
    import_meta,
    import_metabolites,
    import_phospho,
    import_rna,
    syn_login,
)
from pilot.figures.figure_setup import confidence_ellipse, get_setup

FILE_DIR = dirname(abspath(__file__))


def make_figure():
    # Import data
    syn = syn_login()
    meta = import_meta(syn)
    metabolites = import_metabolites(syn)
    lipids = import_lipids(syn)
    phospho = import_phospho(syn)
    global_prot = import_global(syn)
    acetyl = import_acetyl(syn)
    rna = import_rna(syn, return_symbols=True, batch_correct=True, tpm=True)

    cohort_colors = pd.Series(
        {
            "BeatAML 210": "tab:blue",
            "Pilot": "tab:orange",
            "Bridge": "tab:green",
        }
    )
    race_colors = pd.Series(
        {
            "Black": "tab:red",
            "White": "tab:purple",
        }
    )

    # Reformat acetyl
    bridge_ids = (
        phospho[1].loc[phospho[1].index.str.contains("Bridge"), :].index
    )
    bridge_ids = {bridge_id[:-7]: bridge_id for bridge_id in bridge_ids}
    acetyl.rename(index=bridge_ids, inplace=True)
    acetyl = (None, acetyl)

    # Format matplotlib axes
    datasets = {
        "Metabolites": metabolites,
        "Lipidomics": lipids,
        "Transcriptomics": rna,
        "Acetylomics": acetyl,
        "Phosphoproteomics": phospho,
        "Global Proteomics": global_prot,
    }

    fig, axes = get_setup(
        2,
        len(datasets),
        fig_params={
            "figsize": (len(datasets) * 3, 6),
        },
    )

    for col_index, (dataset_name, dataset) in enumerate(datasets.items()):
        # Merge 210 and Pilot cohorts, then scale features
        ptrc, pilot = dataset

        if ptrc is None:
            data = pilot
        elif dataset_name == "Transcriptomics":
            data = pd.concat(dataset, axis=0)
            data = data.loc[data.index.intersection(meta.index), :]
            data.loc[:] = scale(data, axis=1)
        else:
            data = pd.concat(dataset, axis=0)
            assert len(data) == len(pilot) + len(ptrc)

        data.loc[:] = scale(data, axis=0)

        # Label sample cohorts, race
        cohorts = pd.Series("BeatAML 210", index=data.index)
        cohorts.iloc[-pilot.shape[0] :] = "Pilot"
        cohorts.loc[cohorts.index.str.contains("Bridge")] = "Bridge"
        races = meta.loc[data.index, "Race"]

        # Run PCA
        pca = PCA(
            data,
            ncomp=2,
            missing="fill-em",
            demean=False,
            standardize=False,
            method="nipals",
        )

        # Plot scores by cohort
        ax = axes[0, col_index]
        for cohort in cohorts.unique():
            cohort_scores = pca.scores.loc[cohorts == cohort, :]
            confidence_ellipse(
                cohort_scores.iloc[:, 0],
                cohort_scores.iloc[:, 1],
                ax,
                facecolor="None",
                edgecolor=cohort_colors.loc[cohort],
                linestyle="--",
                alpha=0.9,
            )
            ax.scatter(
                cohort_scores.iloc[:, 0],
                cohort_scores.iloc[:, 1],
                c=cohort_colors.loc[cohort],
                s=12,
                edgecolor="black",
                linewidths=1,
                label=cohort,
            )

        ax.set(xlabel="PC 1", ylabel="PC 2", title=f"{dataset_name}: Cohort")
        ax.legend()

        # Plot scores by race
        ax = axes[1, col_index]
        for race in race_colors.index:
            race_scores = pca.scores.loc[races == race, :]
            confidence_ellipse(
                race_scores.iloc[:, 0],
                race_scores.iloc[:, 1],
                ax,
                facecolor="None",
                edgecolor=race_colors.loc[race],
                linestyle="--",
                alpha=0.9,
            )
            ax.scatter(
                race_scores.iloc[:, 0],
                race_scores.iloc[:, 1],
                c=race_colors.loc[race],
                s=12,
                edgecolor="black",
                linewidths=1,
                label=race,
            )

        ax.set(xlabel="PC 1", ylabel="PC 2", title=f"{dataset_name}: Race")
        ax.legend()

    return fig


if __name__ == "__main__":
    make_figure()
