"""Plots figure 3Sa: Acetylation Variance Explained."""

import numpy as np

from pilot.decomposition import run_acetyl_plsr
from pilot.figures.figure_setup import get_setup

COLORS = {
    "RNA-Seq": "tab:orange",
    "Acetyl": "tab:green"
}


def make_figure():
    # Run Acetyl PLSR
    components, _, r2x = run_acetyl_plsr()
    deflate_r2x = r2x["Deflation"]
    a_r2x = r2x["Residual"]

    # Setup figure
    fig, axes = get_setup(
        1,
        2,
        fig_params={
            "figsize": (6, 3)
        }
    )
    
    # Plot each R2X across components
    # There's two variances to consider--that explained by the deflation step, 
    # and that explained by the acetylation on the residual
    for ax, r2x, stage in zip(
        axes,
        [deflate_r2x, a_r2x],
        ["RNA-Seq", "Acetyl"]
    ):
        # Plot global variance explained
        ax.bar(
            np.arange(0, 3 * r2x.shape[0], 3),
            r2x.loc[:, "X"],
            width=1,
            label="Global",
            color="tab:blue"
        )
        # Plot RNA-seq/acetyl variance explained
        ax.bar(
            np.arange(1, 3 * r2x.shape[0], 3),
            r2x.loc[:, "Y"],
            width=1,
            label=stage,
            color=COLORS[stage]
        )

        # Format title and label axes
        title = f"Global x. {stage}"
        if ax == axes[-1]:
            title = "Deflated " + title

        ax.set(
            xlabel="PLSR Component",
            ylabel="Variance Explained",
            title=title,
            xticks=np.arange(0.5, 3 * r2x.shape[0], 3),
            xticklabels=np.arange(r2x.shape[0]) + 1
        )
        ax.legend()

    return fig


if __name__ == "__main__":
    make_figure()
