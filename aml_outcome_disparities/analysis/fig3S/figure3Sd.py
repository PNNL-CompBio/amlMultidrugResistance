"""Plots figure 3Sd: PLSR loadings."""

import numpy as np

from pilot.decomposition import run_acetyl_plsr
from pilot.figures.figure_setup import get_setup


def make_figure():
    # Run acetyl PLSR, extract acetylation loadings
    components, _, _ = run_acetyl_plsr()

    # Setup figure
    fig, axes = get_setup(
        4,
        4,
        fig_params={
            "figsize": (12, 8),
        }
    )

    # Iterate through acetyl and proteomic loadings
    for col_offset, df in enumerate(
        [
            components["Acetyl Loadings"],
            components["Proteomic Loadings"]
        ]
    ):
        # Iterate through components, plotting most negatively and positively
        # weighted acetylations
        for row_index, comp in enumerate(df.columns):
            comp_loadings = df.loc[:, comp].sort_values(
                ascending=True
            )
            for col_index, loadings in enumerate([
                comp_loadings.iloc[:10],
                comp_loadings.iloc[-10:]
            ]):
                ax = axes[row_index, col_offset * 2 + col_index]
                ax.bar(
                    np.arange(len(loadings)),
                    loadings,
                )
                ax.set_xticks(np.arange(len(loadings)))
                if col_index == 0:
                    ax.set_ylabel(f"PLSR {row_index + 1} Loading")
                ax.set_xticklabels(
                    loadings.index,
                    rotation=45,
                    ha="right",
                    ma="right",
                    va="top"
                )

    return fig


if __name__ == "__main__":
    make_figure()
