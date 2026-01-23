"""Plots figure 1b: Preliminary Volcano Plots"""

import sys
from os.path import abspath, dirname, join

sys.path.append(join(dirname(dirname(dirname(abspath(__file__)))), "src", "python"))

import numpy as np
import pandas as pd
from sklearn.preprocessing import scale

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
from pilot.figures.figure_setup import get_setup
from pilot.gene_analysis import calculate_fc, volcano_plot

FILE_DIR = dirname(abspath(__file__))


def make_figure():
    # Import data
    syn = syn_login()
    meta = import_meta(syn)
    metabolites = import_metabolites(syn)
    lipids = import_lipids(syn)
    phospho = import_phospho(syn, pre_corrected=False)
    global_prot = import_global(syn, pre_corrected=True)
    acetyl = import_acetyl(syn)
    rna = import_rna(syn, return_symbols=True)
    races = meta.loc[:, "Race"]

    # Reformat acetyl
    bridge_ids = phospho[1].loc[phospho[1].index.str.contains("Bridge"), :].index
    bridge_ids = {bridge_id[:-7]: bridge_id for bridge_id in bridge_ids}
    acetyl.rename(index=bridge_ids, inplace=True)
    acetyl = (None, acetyl)

    datasets = {
        "Metabolites": metabolites,
        "Lipids": lipids,
        "Transcriptomics": rna,
        "Acetylomics": acetyl,
        "Phosphoproteomics": phospho,
        "Global Proteomics": global_prot,
    }

    # Setup figure
    fig, axes = get_setup(
        2,
        3,
        fig_params={
            "figsize": (9, 6),
        },
    )
    axes = axes.flatten()

    # Make volcano plots for each dataset
    for row_index, (dataset_name, dataset) in enumerate(datasets.items()):
        # Merge 210 and Pilot cohorts
        (ptrc, pilot) = dataset

        if ptrc is None:
            data = pilot
        elif dataset_name == "Transcriptomics":
            data = pd.concat(dataset, axis=0)
            data = data.loc[data.index.intersection(meta.index), :]
            data.loc[:] = scale(data, axis=1)
            data.loc[:] = np.log2(data)
            data.replace({np.inf: np.nan, -np.inf: np.nan}, inplace=True)
        else:
            data = pd.concat(dataset, axis=0)

        # Separate Black and White patients, then derive log-fold change
        # and t-test (adjusted) p-values
        black_patients = data.loc[races.loc[data.index] == "Black", :]
        white_patients = data.loc[races.loc[data.index] == "White", :]
        fc, p_values = calculate_fc(
            black_patients,
            white_patients,
        )
        p_values.loc[p_values == 0] = p_values.loc[p_values != 0].min()

        # Plot the volcano plot
        ax = axes[row_index]
        black_over_exp, white_over_exp, _ = volcano_plot(
            fc,
            p_values,
            ax=ax,
            x_max=int(np.ceil(abs(fc).max())),
            y_max=int(np.ceil(-np.log10(p_values.min()))),
        )

        # Include counts of over-expressed counts
        ax.text(
            0.99,
            0.99,
            s=f"Black over-expressed: {len(black_over_exp)}\n"
            f"White over-expressed: {len(white_over_exp)}",
            transform=ax.transAxes,
            ha="right",
            ma="right",
            va="top",
        )
        ax.set(xlabel="log2(Fold Change)", ylabel="-log10(p)", title=dataset_name)

    return fig


if __name__ == "__main__":
    make_figure()
