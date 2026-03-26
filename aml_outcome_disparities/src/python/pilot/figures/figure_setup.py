"""Functions for figure setup and consistent style."""

from typing import Any, Dict, Iterable

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
import statsmodels.api as sm
from matplotlib.patches import Ellipse
from numpy.typing import ArrayLike
from scipy.stats import false_discovery_control, ttest_ind
from sklearn.preprocessing import LabelEncoder
from statsmodels.regression.linear_model import RegressionResults

matplotlib.rcParams["axes.labelsize"] = 8
matplotlib.rcParams["axes.linewidth"] = 0.6
matplotlib.rcParams["axes.titlesize"] = 10
matplotlib.rcParams["font.family"] = ["sans-serif"]
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams["font.size"] = 8
matplotlib.rcParams["grid.linestyle"] = "dotted"
matplotlib.rcParams["legend.borderpad"] = 0.35
matplotlib.rcParams["legend.fontsize"] = 8
matplotlib.rcParams["legend.framealpha"] = 0.5
matplotlib.rcParams["legend.handlelength"] = 0.5
matplotlib.rcParams["legend.handletextpad"] = 0.5
matplotlib.rcParams["legend.labelspacing"] = 0.2
matplotlib.rcParams["legend.markerscale"] = 0.7
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["xtick.labelsize"] = 8
matplotlib.rcParams["xtick.major.pad"] = 1.0
matplotlib.rcParams["xtick.minor.pad"] = 0.9
matplotlib.rcParams["ytick.labelsize"] = 8
matplotlib.rcParams["ytick.major.pad"] = 1.0
matplotlib.rcParams["ytick.minor.pad"] = 0.9


def confidence_ellipse(
    x: ArrayLike,
    y: ArrayLike,
    ax: plt.Axes,
    n_std: float = 3.0,
    facecolor: str | None = None,
    **kwargs,
):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.
    Taken from: https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The Axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    facecolor : str or None, default None
        Fill color for ellipses. Defaults to None, which is transparent.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse(
        (0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs,
    )

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = (
        transforms.Affine2D()
        .rotate_deg(45)
        .scale(scale_x, scale_y)
        .translate(mean_x, mean_y)
    )

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def run_ols(
    data: Iterable[float], labels: Iterable[float], ax: plt.Axes | None = None
) -> RegressionResults:
    """
    Runs OLS to regress data against labels.

    Args:
        data (Iterable[float]): Data to regress.
        labels (Iterable[float]): Labels to regress data against.
        ax (plt.Axes): Axes to add regression line to.

    Returns:
        OLS model fit results.
    """
    data = sm.add_constant(data, prepend=False)
    model = sm.OLS(labels, data)
    result = model.fit()

    if ax is not None:
        lims = ax.get_xlim()
        ax.plot(
            lims,
            [
                result.params.loc[data.columns[0]] * lims[0]
                + result.params.loc["const"],
                result.params.loc[data.columns[0]] * lims[1]
                + result.params.loc["const"],
            ],
            linestyle="--",
            color="black",
        )
        ax.set_xlim(lims)

    return result


def get_setup(
    n_rows: int, n_cols: int, fig_params: Dict[str, Any] | None = None
) -> tuple[plt.Figure, np.ndarray | plt.Axes]:
    """
    Builds subplot figure and axes.

    Args:
        n_rows (int): Number of rows.
        n_cols (int): Number of columns.
        fig_params (Dict[str: Any]): Matplotlib subplot params.

    Returns
        plt.Figure: Matplotlib figure.
        plt.Axes: Matplotlib axes.
    """
    if fig_params is None:
        fig_params = {"constrained_layout": True, "dpi": 300}
    if "dpi" not in fig_params.keys():
        fig_params["dpi"] = 300
    if "constrained_layout" not in fig_params.keys():
        fig_params["constrained_layout"] = True

    fig, axes = plt.subplots(n_rows, n_cols, **fig_params)

    return fig, axes


def compare_means(
    data: pd.DataFrame,
    groups: Iterable,
    x_within: float,
    x_between: float,
    ax: plt.Axes,
    alternative: str = "two-sided",
    star_only: bool = False,
    bracket_height: float = 0.01,
    bracket_space: float = 0.05
) -> plt.Axes:
    """
    Compares means of two groups via t-test and denotes p-value in plot.

    Args:
        data (pd.DataFrame): Data to compare between groups.
        groups (Iterable): Group membership for each row in data.
        x_within (float): Space between groups within each compared measurement.
        x_between (float): Space between each measurement.
        ax (plt.Axes): Axes to add significance test to.
        alternative (str): Alternative hypothesis to test. Accepts same values
            as scipy.stats.ttest_ind.
        star_only (bool, default: False): Use * in place of exact p-values.
        bracket_height (float, default: 0.05): Height of bracket.
        bracket_space (float, default: 0.05): Space between bracket and max
            measurement value.

    Returns
        plt.Axes: Matplotlib axes.
    """
    if alternative not in ["two-sided", "less", "greater"]:
        raise ValueError("alternative must be one of 'two-sided', 'less', "
                         "or 'greater'.")

    # Encode groups for standardized use
    le = LabelEncoder()
    groups = le.fit_transform(groups)

    # Perform t-test, correct for FDR
    res = ttest_ind(
        data.loc[groups == 0, :],
        data.loc[groups == 1, :],
        nan_policy="omit",
        alternative=alternative
    )
    p_values = false_discovery_control(res.pvalue)

    # Convert p-values to *
    if star_only:
        stars = np.repeat("n.s.", len(p_values))
        stars[p_values < 0.05] = "*"
        stars[p_values < 0.01] = "**"
        stars[p_values < 0.001] = "***"

    # Add comparison brackets
    vert_scale = ax.get_ylim()
    vert_scale = vert_scale[1] - vert_scale[0]

    for index, (p_value, column) in enumerate(zip(p_values, data.columns)):
        col_max = data.loc[:, column].max()
        x_pos = [index * x_between] * 2 + [index * x_between + x_within] * 2
        y_pos = [
            col_max + vert_scale * bracket_space,
            col_max + vert_scale * (bracket_space + bracket_height),
            col_max + vert_scale * (bracket_space + bracket_height),
            col_max + vert_scale * bracket_space
        ]
        ax.plot(x_pos, y_pos, color="black")
        if star_only:
            ax.text(
                index * x_between + x_within / 2,
                col_max + vert_scale * (bracket_space + bracket_height),
                s=stars[index],
                ha="center",
                va="bottom"
            )
        else:
            ax.text(
                index * x_between + x_within / 2,
                col_max + vert_scale * (bracket_space + bracket_height),
                s=round(p_values[index], 4),
                ha="center",
                va="bottom"
            )

    return ax
