from enum import StrEnum
from typing import Optional, Literal

import polars as pl
from matplotlib import (
    pyplot as plt, 
    rcParams, 
    patches as mpatches, 
    colors as mcolors, 
    cm
)
from matplotlib.axes import Axes

from colors import RACE_MUT_TO_COLOR


rcParams["font.family"] = "Arial"
rcParams["font.size"] = 7


class _COLORS(StrEnum):
    LABEL_TEXT = "#FFF"
    EDGES = "#000"
    CMAP_MIN = "#44C"
    CMAP_MID = "#CCC"
    CMAP_MAX = "#C44"
    POP_BLACK = "#000"
    POP_WHITE = "#FFF"


# Custom colormap and normalizer centered around 0 in range [-N, +N]
N = 100
CMAP = mcolors.LinearSegmentedColormap.from_list(
    "custom_map", 
    [_COLORS.CMAP_MIN, _COLORS.CMAP_MID, _COLORS.CMAP_MAX]
)
NORM = mcolors.Normalize(vmin=-N, vmax=N, clip=True)


TRI_INNER_PAD = 0.025
LW = 0.75
LABEL_SIZE = 7.5
TEXT_INSET = 0.3
MUT_RADIUS = 0.125
MUT_LABEL_SIZE = 5.5


def _get_mapped_color(count: int) -> str:
    return mcolors.to_hex(CMAP(NORM(count)))


def _plot_counts_triangles(
        ax: Axes,
        n_up: Optional[int],
        n_down: Optional[int],
):
    # setup
    ax.set_axis_off()
    ax.set_aspect("equal")
    ax.set_xlim((0 - TRI_INNER_PAD, 1 + TRI_INNER_PAD))
    ax.set_ylim((0 - TRI_INNER_PAD, 1 + TRI_INNER_PAD))

    if n_up is not None:
        # top triangle
        ax.add_patch(mpatches.Polygon(
            [(0, 0), (0, 1), (1, 1)],
            ec=_COLORS.EDGES,
            fc=_get_mapped_color(n_up),
            lw=LW,
            joinstyle="round",
        ))
        # top label
        ax.text(
            TEXT_INSET, 1. - TEXT_INSET,
            f"{n_up:,}",
            fontsize=LABEL_SIZE,
            color=_COLORS.LABEL_TEXT,
            fontweight="bold",
            ha="center", va="center",
            rotation=45
        )
    else:
        # top triangle
        ax.add_patch(mpatches.Polygon(
            [(0, 0), (0, 1), (1, 1)],
            ec=_COLORS.CMAP_MID,
            fc="none",
            lw=LW,
            joinstyle="round",
        ))
        # top label
        ax.text(
            TEXT_INSET, 1. - TEXT_INSET,
            "0",
            fontsize=LABEL_SIZE,
            color=_COLORS.CMAP_MID,
            fontweight="bold",
            ha="center", va="center",
            rotation=45
        )

    if n_down is not None:
        # bottom triangle
        ax.add_patch(mpatches.Polygon(
            [(0, 0), (1, 1), (1, 0)],
            ec=_COLORS.EDGES,
            fc=_get_mapped_color(-n_down),
            lw=LW,
            joinstyle="round",
        ))
        # bottom label
        ax.text(
            1 - TEXT_INSET, TEXT_INSET,
            f"{n_down:,}",
            fontsize=LABEL_SIZE,
            color=_COLORS.LABEL_TEXT,
            fontweight="bold",
            ha="center", va="center",
            rotation=45
        )
    else:
        # bottom triangle
        ax.add_patch(mpatches.Polygon(
            [(0, 0), (1, 1), (1, 0)],
            ec=_COLORS.CMAP_MID,
            fc="none",
            lw=LW,
            joinstyle="round",
        ))
        # bottom label
        ax.text(
            1 - TEXT_INSET, TEXT_INSET,
            "0",
            fontsize=LABEL_SIZE,
            color=_COLORS.CMAP_MID,
            fontweight="bold",
            ha="center", va="center",
            rotation=45
        )


def _plot_comparison_metadata(
        ax: Axes,
        c_left: str, c_right: str,
        top_npm1: bool,
        top_flt3: bool,
        bottom_npm1: bool,
        bottom_flt3: bool,

):
    # setup
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal")
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 1))

    # box indicating population
    box_inset = 0.08
    box_height = 0.2
    w = (1 - 4 * box_inset) / 2
    ax.add_patch(mpatches.Rectangle(
        (box_inset, 1 - box_inset / 2 - box_height),
        w, 
        box_height,
        ec=_COLORS.EDGES,
        fc=c_left,
        lw=LW,
    ))
    ax.add_patch(mpatches.Rectangle(
        (w + 3 * box_inset, 1 - box_inset / 2 - box_height),
        w, 
        box_height,
        ec=_COLORS.EDGES,
        fc=c_right,
        lw=LW,
    ))

    # mutations
    metadata_divider_x_inset = 0.1
    xmin = metadata_divider_x_inset
    xmax = 1 - metadata_divider_x_inset
    ymin = metadata_divider_x_inset
    divider_height = 0.5
    ymax = ymin + divider_height
    ax.plot(
        [xmin, xmax],
        [ymin, ymax],
        ls="--",
        lw=LW,
        c="k",
    )

    dx = 0.1 
    dy = 0.2
    for point, flag, mut in zip(
        [
            (xmin + dx / 2, ymin + 1.5 * dy), 
            (xmin + 3 * dx, ymin + 2 * dy), 
            (xmax - dx / 2, ymax - 1.5 * dy), 
            (xmax - 3 * dx, ymax - 2 * dy)
        ],
        [top_npm1, top_flt3, bottom_flt3, bottom_npm1],
        "NFFN"
    ):
        ax.add_patch(mpatches.Circle(
            point,
            MUT_RADIUS,
            ec=_COLORS.EDGES,
            fc=_COLORS.EDGES if flag else _COLORS.LABEL_TEXT,
            lw=LW
        ))
        ax.text(
            point[0], point[1] - 0.01,
            mut,
            fontsize=MUT_LABEL_SIZE,
            color=_COLORS.LABEL_TEXT if flag else _COLORS.EDGES,
            fontweight="bold",
            ha="center", va="center"
        )
    
    ax.set_xticks([])
    ax.set_yticks([])


def plot_sig_up_down_feature_counts(
    results: pl.DataFrame, 
    lfc: float, 
    p_value: float, 
    p_adj: bool,
    figname: str
): 
    # key: (contrast, omics type, direction)
    # value: feature count
    counts = {
        tuple(k): v
        for *k, v
        in (
            results
            .filter(
                (pl.col("logFC").abs() > lfc) 
                & (pl.col("adj.P.Val" if p_adj else "P.Value") <= p_value)
            )
            .with_columns(
                pl.when(pl.col("logFC") > 0)
                .then(pl.lit("up"))
                .otherwise(pl.lit("down"))
                .alias("Direction")
            )
            .group_by("OmicsType", "Contrast", "Direction")
            .len()
            .sort("Contrast", "OmicsType", "Direction")
            .select("Contrast", "OmicsType", "Direction", "len")
            .iter_rows()
        )
    }

    SUBPLOT_PAD = 0.1
    SCALE = 0.55

    fig, axs = plt.subplots(
        figsize=(12 * SCALE, 7 * SCALE),
        ncols=12, nrows=7,
    )


    # first columns: comparison of white vs. black (within specific mutations)
    x = 0
    for mutations in [
        (False, False),
        (True, False),
        (False, True),
        (True, True),
    ]:
        muts = (
            f"NPM1{'+' if mutations[0] else '-'}/"
            f"FLT3-ITD{'+' if mutations[1] else '-'}"
        )
        comparison = (
            "Black_"
            + muts
            + "_vs_White_"
            + muts
        )
        _plot_comparison_metadata(
            axs[0][x], 
            RACE_MUT_TO_COLOR[("Black", muts)], 
            RACE_MUT_TO_COLOR[("White", muts)], 
            *mutations, 
            *mutations
        )
        for y, omics_data_type in enumerate(
                [
                    "Transcriptomics",
                    "Proteomics",
                    "Phosphoproteomics",
                    "Acetylomics",
                    "Metabolomics",
                    "Lipidomics"
                ], 
                start=1
            ):
                _plot_counts_triangles(
                    axs[y][x],
                    counts.get((comparison, omics_data_type, "up")),
                    counts.get((comparison, omics_data_type, "down")),
                    
                )
        # next column
        x += 1

    for mutations in [
        # NPM1 (top), FLT3-ITD (top), NPM1 (bottom), FLT3-ITD (bottom)
        (True, False, False, False),
        (False, True, False, False),
        (True, True, True, False),
        (True, True, False, True),
    ]:
        for pop in ["Black", "White"]:
            muts_top = (
                f"NPM1{'+' if mutations[0] else '-'}/"
                f"FLT3-ITD{'+' if mutations[1] else '-'}"
            )
            muts_bot = (
                f"NPM1{'+' if mutations[2] else '-'}/"
                f"FLT3-ITD{'+' if mutations[3] else '-'}"
            )
            # metadata description
            _plot_comparison_metadata(
                axs[0][x], 
                RACE_MUT_TO_COLOR[(pop, muts_top)],
                RACE_MUT_TO_COLOR[(pop, muts_bot)],
                *mutations
            )
            comparison = (
                f"{pop}_"
                + muts_top
                + f"_vs_{pop}_"
                + muts_bot
            )
            # triangles
            for y, omics_data_type in enumerate(
                [
                    "Transcriptomics",
                    "Proteomics",
                    "Phosphoproteomics",
                    "Acetylomics",
                    "Lipidomics",
                    "Metabolomics",
                ], 
                start=1
            ):
                _plot_counts_triangles(
                    axs[y][x],
                    counts.get((comparison, omics_data_type, "up")),
                    counts.get((comparison, omics_data_type, "down")),
                    
                )
            # next column
            x += 1

    fig.subplots_adjust(wspace=SUBPLOT_PAD, hspace=SUBPLOT_PAD)
    plt.savefig(
        figname,
        dpi=400,
        bbox_inches="tight",
        transparent=True
    )
    plt.show()
    plt.close()


def plot_sig_up_down_pathway_counts(
        results: pl.DataFrame, 
        figname: str
): 
    # key: (contrast, omics type, direction)
    # value: feature count
    counts = {
        tuple(k): v
        for *k, v
        in (
            results
            .with_columns(
                pl.when(pl.col("NES") > 0)
                .then(pl.lit("up"))
                .otherwise(pl.lit("down"))
                .alias("Direction")
            )
            .group_by("OmicsType", "Contrast", "Direction")
            .len()
            .sort("Contrast", "OmicsType", "Direction")
            .select("Contrast", "OmicsType", "Direction", "len")
            .iter_rows()
        )
    }

    SUBPLOT_PAD = 0.1
    SCALE = 0.55

    fig, axs = plt.subplots(
        figsize=(12 * SCALE, 3 * SCALE),
        ncols=12, nrows=3,
    )


    # first columns: comparison of white vs. black (within specific mutations)
    x = 0
    for mutations in [
        (False, False),
        (True, False),
        (False, True),
        (True, True),
    ]:
        muts = (
            f"NPM1{'+' if mutations[0] else '-'}/"
            f"FLT3-ITD{'+' if mutations[1] else '-'}"
        )
        comparison = (
            "Black_"
            + muts
            + "_vs_White_"
            + muts
        )
        _plot_comparison_metadata(
            axs[0][x], 
            RACE_MUT_TO_COLOR[("Black", muts)], 
            RACE_MUT_TO_COLOR[("White", muts)], 
            *mutations, 
            *mutations
        )
        for y, omics_data_type in enumerate(
                ["Transcriptomics", "Proteomics"], 
                start=1
            ):
                _plot_counts_triangles(
                    axs[y][x],
                    counts.get((comparison, omics_data_type, "up")),
                    counts.get((comparison, omics_data_type, "down")),
                    
                )
        # next column
        x += 1

    for mutations in [
        # NPM1 (top), FLT3-ITD (top), NPM1 (bottom), FLT3-ITD (bottom)
        (True, False, False, False),
        (False, True, False, False),
        (True, True, True, False),
        (True, True, False, True),
    ]:
        for pop in ["Black", "White"]:
            muts_top = (
                f"NPM1{'+' if mutations[0] else '-'}/"
                f"FLT3-ITD{'+' if mutations[1] else '-'}"
            )
            muts_bot = (
                f"NPM1{'+' if mutations[2] else '-'}/"
                f"FLT3-ITD{'+' if mutations[3] else '-'}"
            )
            # metadata description
            _plot_comparison_metadata(
                axs[0][x], 
                RACE_MUT_TO_COLOR[(pop, muts_top)],
                RACE_MUT_TO_COLOR[(pop, muts_bot)],
                *mutations
            )
            comparison = (
                f"{pop}_"
                + muts_top
                + f"_vs_{pop}_"
                + muts_bot
            )
            # triangles
            for y, omics_data_type in enumerate(
                ["Transcriptomics", "Proteomics"], 
                start=1
            ):
                _plot_counts_triangles(
                    axs[y][x],
                    counts.get((comparison, omics_data_type, "up")),
                    counts.get((comparison, omics_data_type, "down")),
                    
                )
            # next column
            x += 1

    fig.subplots_adjust(wspace=SUBPLOT_PAD, hspace=SUBPLOT_PAD)
    plt.savefig(
        figname,
        dpi=400,
        bbox_inches="tight",
        transparent=True
    )
    plt.show()
    plt.close()


def pw_bubble_plot(
        results_df: pl.DataFrame, 
        contrast: str, 
        gene_set_collection: str, 
        figname: str,
        limit_top_n: int = 20
):
    subset = (
        results_df
        .filter(
            (pl.col("Contrast") == contrast)
            & (pl.col("GeneSet") == gene_set_collection)
        )
        .select(
            # shorten OmicsType to P or T
            pl.col("OmicsType").str.slice(0, 1), 
            "Term", 
            pl.col("gene_set_size").alias("Size"),
            "NES", 
            pl.col("Tag %").alias("Coverage"),
            # convert FDR q-values to -log10 scale
            # clip 0s to minimum 1e-5 for sensible colorbar scale (and no infs)
            -(pl.col("FDR q-val").clip(1e-5).log10()).alias("log10FDR"),
        )
        # combine by term, keep separate NES/log10FDR
        .pivot(
            index="Term",
            on="OmicsType",
            values=["NES", "log10FDR", "Coverage", "Size"],
        )
        # combine size columns
        .with_columns(
            (
                pl.when(pl.col("Size_T").is_null())
                .then(pl.col("Size_P").cast(pl.String))
                .when(pl.col("Size_P").is_null())
                .then(pl.col("Size_T").cast(pl.String))
                .when(pl.col("Size_P") == pl.col("Size_T"))
                .then(pl.col("Size_T").cast(pl.String))
                .otherwise(
                    pl.concat_str("Size_T", "Size_P", separator="|")
                )
            )
            .alias("Size")
        )
        # sort by ascending |NES| 
        # (averaged between P and T if both present)
        .sort(pl.mean_horizontal("NES_T", "NES_P").abs(), descending=True)
        .head(limit_top_n)
        .with_row_index()
        .with_columns(
            (limit_top_n - pl.col("index")).alias("index")
        )
        .sort("index")
        # reorder columns
        .select(
            "index", "Term", "Size", 
            "NES_T", "log10FDR_T", "Coverage_T", 
            "NES_P", "log10FDR_P", "Coverage_P"
        )
    )

    height_mult = 0.2

    fig, axs = plt.subplot_mosaic(
        [
            ["bubble", "cbar"],
            ["bubble", "sizes"]
        ],
        width_ratios=[15, 1],
        height_ratios=[1, 3],
        figsize=(4, height_mult * limit_top_n),
        gridspec_kw={"wspace": 0.5}
    )

    def _marker_size(coverage: float, tri: bool = False) -> float : 
        max_sz = 200.
        tri_sz_scale = 1.1 if tri else 1.
        return (coverage / 100.)**2 * max_sz * tri_sz_scale

    # populate the bubble panel
    axs["bubble"].grid()
    axs["bubble"].set_axisbelow(True)
    axs["bubble"].axvline(0, c="k", lw=0.75)
    for idx, *_, nes_t, fdr_t, cov_t, nes_p, fdr_p, cov_p in subset.iter_rows():
        if nes_t is not None:
            axs["bubble"].scatter(
                nes_t, idx,
                # square coverage so marker area scales properly
                s=_marker_size(cov_t, tri=True),  
                c=fdr_t,
                cmap="viridis", vmin=0, vmax=5,
                marker="^",
                linewidths=0.5,
                edgecolors="k",
                # adjust zorder to avoid large circle covering small triangle
                zorder=2 if (cov_p is not None) and (cov_p > cov_t) else 1
            )
        if nes_p is not None:
            axs["bubble"].scatter(
                nes_p, idx,
                # square coverage so marker area scales properly
                s=_marker_size(cov_p),  # type: ignore
                c=fdr_p,
                cmap="viridis", vmin=0, vmax=5,
                marker="o",
                linewidths=0.5,
                edgecolors="k"
            )
        
    for d in ["bottom", "left", "right"]: 
        axs["bubble"].spines[d].set_visible(False)
    axs["bubble"].xaxis.tick_top()
    axs["bubble"].xaxis.set_label_position('top')
    axs["bubble"].set_xlabel("NES")
    axs["bubble"].set_xlim((-4.5, 4.5))
    second_y = axs["bubble"].secondary_yaxis("right")
    for ax in [axs["bubble"], second_y]:
        ax.set_yticks(subset["index"])
        ax.set_ylim((-1, limit_top_n + 1))  # type: ignore
    axs["bubble"].set_yticklabels(subset["Term"])
    second_y.spines["right"].set_visible(False)
    second_y.set_yticklabels(subset["Size"])

    # FDR q-value colorbar
    cbar = fig.colorbar(
        cm.ScalarMappable(norm=mcolors.Normalize(vmin=0, vmax=5), cmap="viridis"), 
        cax=axs["cbar"], 
        orientation="vertical", 
        shrink=False
    )
    cbar.set_label("-log10(FDR q-value)")

    # sizes panel
    axs["sizes"].text(
        0.5, 0.9,
        f"coverage",
        ha="center",
        va="top"
    )
    _y2, _d = 0.5, 0.1
    for i, cov in enumerate([25, 50, 75]):
        _dy = i**1.25 * _d
        _ty = 0.07 * cov / 75
        axs["sizes"].scatter(
            0.5, _y2 + _dy,
            s=_marker_size(cov, tri=True),
            marker="^",
            c="w",
            linewidths=0.5,
            edgecolors="k"
        )
        axs["sizes"].text(
            0.5, _y2 + _dy - _ty,
            f"{cov}%",
            ha="center",
            va="top"
        )
        axs["sizes"].scatter(
            0.5, _y2 + _dy,
            s=_marker_size(cov),
            marker="o",
            c="#00000000",
            linewidths=0.5,
            edgecolors="k",
        )
    axs["sizes"].set_xlim((-0.2, 1.2))
    axs["sizes"].set_ylim((0, 1))
    axs["sizes"].set_axis_off()

    plt.savefig(
        figname,
        dpi=400,
        bbox_inches="tight",
        transparent=True
    )
    plt.show()
    plt.close()
