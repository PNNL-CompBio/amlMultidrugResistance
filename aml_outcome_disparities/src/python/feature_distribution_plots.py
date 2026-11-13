from enum import StrEnum
from typing import Optional

import polars as pl
from matplotlib import (
    pyplot as plt, 
    rcParams, 
    patches as mpatches, 
    colors as mcolors, 
    cm
)
from matplotlib.axes import Axes
import numpy as np
from scipy.stats import gaussian_kde


rcParams["font.family"] = "Arial"
rcParams["font.size"] = 7


def _compute_kde(data, bw_adjust=1.0, pad=0.1, n_points=50):
    """
    Computes (x, y) coordinates for a 1D KDE curve.
    
    Parameters:
    -----------
    data : array-like
        1D list or array of numeric values.
    x_grid : array-like
        Predefined evaluation points
    bw_adjust : float
        Multiplicative factor to adjust the default bandwidth (Scott's rule).
        Values < 1.0 make the curve tighter/wigglier; > 1.0 make it smoother.
    """
    data = np.asarray(data)
    data = data[~np.isnan(data)]  # Drop NaNs
    
    # fit KDE
    kde = gaussian_kde(data)
    
    # adjust bandwidth if desired
    if bw_adjust != 1.0:
        kde.set_bandwidth(kde.factor * bw_adjust)  # type: ignore

    # define range to evaluate over
    data_min, data_max = data.min(), data.max()
    range_span = data_max - data_min
    x_min = data_min - pad * range_span
    x_max = data_max + pad * range_span
    x_grid = np.linspace(x_min, x_max, n_points)
    
    # evaluate densities
    y_dens = kde(x_grid)
    
    return x_grid, y_dens


def two_group_ridge(
        df: pl.DataFrame,
        color_a: str,
        color_b: str,
        figname: Optional[str] = None,
        x_min: float = -1.,
        x_max: float = 1.,
        kde_pad: float = 0.1,
        kde_points: int = 50,
        kde_height: float = 1.1,
        fig_height_mult: float = 0.3,
        fig_width: float = 3.,
        xlabel: Optional[str] = None,
        bw_adjust: float = 1.,
        white_outline: bool = False,
        transparent: bool = True,
): 
    """
    dataframe must have columns: 
    - "index" (type: int)
    - "feature" (type: str)
    - "A" (type: list[float])
    - "B" (type: list[float])
    """
    n = len(df)

    fig, ax = plt.subplots(
        figsize=(fig_width, fig_height_mult * n),
    )

    # populate the bubble panel
    ax.grid()
    ax.set_axisbelow(True)
    ax.axvline(0, c="k", lw=0.75)

    for i, a, b in df.select("index", "A", "B").iter_rows():
        # compute and scale KDEs from LFC data
        xa, ya = _compute_kde(a, bw_adjust=bw_adjust, pad=kde_pad, n_points=kde_points)
        ya = (ya / ya.max()) * kde_height
        xb, yb = _compute_kde(b, bw_adjust=bw_adjust, pad=kde_pad, n_points=kde_points)
        yb = (yb / yb.max()) * kde_height

        # fill and outline
        ax.fill_between(xa, i, i + ya, color=color_a, alpha=0.5, zorder=1 + n - i)
        ax.plot(xa, i + ya, color="w" if white_outline else color_a, lw=0.75, zorder=1.5 + n - i)
        ax.fill_between(xb, i, i + yb, color=color_b, alpha=0.5, zorder=1 + n - i)
        ax.plot(xb, i + yb, color="w" if white_outline else color_b, lw=0.75, zorder=1.5 + n - i)
        
    for d in ["top", "left", "right"]: 
        ax.spines[d].set_visible(False)
    #ax.xaxis.tick_top()
    #ax.xaxis.set_label_position('top')
    if xlabel is not None: 
        ax.set_xlabel(xlabel)
    ax.set_ylim((0, n -0.9 + kde_height))
    ax.set_yticks(df["index"])
    ax.set_yticklabels(df["feature"])
    ax.set_xlim((x_min, x_max))

    if figname is not None:
        plt.savefig(
            figname,
            dpi=400,
            bbox_inches="tight",
            transparent=transparent
        )
    plt.show()
    plt.close()

