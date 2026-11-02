#!/usr/bin/env python3
"""Render a grouped-bar alternative to the component-ablation figure."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import NullFormatter
import numpy as np


HERE = Path(__file__).resolve()
REPO = HERE.parents[3]
RESULTS = REPO / "invariant_region_mlp/results/final_component_ablation"
FIGURES = RESULTS / "figures"

PDES = ("hjb_log", "joint_concave", "joint_convex", "nonconvex_coupled")
PDE_LABELS = {
    "hjb_log": "Radial Log-Terminal\nHJB",
    "joint_concave": "Concave Cubic–\nQuadratic",
    "joint_convex": "Convex Cubic–\nQuadratic",
    "nonconvex_coupled": "Oscillatory Value–\nGradient",
}
ABS_ORDER = ("a0_full", "a2_no_terminal_anchor", "a1_no_mcv", "a3_no_updates")
ABS_LABELS = ("Full", "w/o\nanchor", "w/o\nMCV", "w/o\nupdates")

# Component-consistent color families.  The exact endpoint colors are reused
# between the absolute bars (left) and their matching diagnostics (right).
MCV_COLORS = ("#F0B2A3", "#DD8C86", "#C66B73")
ANCHOR_COLOR = "#AA9BC4"
ACTOR_UPDATE_COLORS = ("#BBDCE8", "#98C5D9", "#74ACC8", "#4F91B2", "#2F7095")
ABS_COLORS = (
    ACTOR_UPDATE_COLORS[-1],  # Full <-> selected/final Actor A4.
    ANCHOR_COLOR,             # w/o anchor <-> w/o anchor.
    MCV_COLORS[-1],           # w/o MCV <-> w/o all MCV.
    ACTOR_UPDATE_COLORS[0],   # w/o updates <-> initial Actor A0.
)
ABS_HATCHES = (None, None, None, None)

RATIO_X = np.asarray((0.0, 1.0, 2.0, 4.0, 6.0, 7.0, 8.0, 9.0, 10.0))
RATIO_LABELS = (
    "w/o\ninfer",
    "w/o\ntrain",
    "w/o\nall",
    "w/o\nanchor",
    "$A_0$",
    "$A_1$",
    "$A_2$",
    "$A_3$",
    "$A_4$",
)
RATIO_COLORS = (
    *MCV_COLORS,
    ANCHOR_COLOR,
    *ACTOR_UPDATE_COLORS,
)
RATIO_HATCHES = (None, None, None, None, None, None, None, None, None)
RATIO_BAR_FLOOR = 0.70


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def soft_background(ax: plt.Axes) -> None:
    gradient = np.linspace(0.0, 1.0, 256).reshape(-1, 1)
    cmap = LinearSegmentedColormap.from_list(
        "ablation_soft_blue_rose",
        (
            (0.00, "#E2EEF5"),
            (0.25, "#E8F0F3"),
            (0.50, "#F1F1F0"),
            (0.75, "#F3EAEB"),
            (1.00, "#F4E6E8"),
        ),
    )
    ax.imshow(
        gradient,
        cmap=cmap,
        origin="lower",
        aspect="auto",
        extent=(0, 1, 0, 1),
        transform=ax.transAxes,
        interpolation="bicubic",
        zorder=-30,
    )


def sci(value: float) -> str:
    return f"{value:.1e}".replace("e-0", "e-").replace("e+0", "e+")


def main() -> None:
    summary_rows = read_csv(RESULTS / "ABLATION_SUMMARY.csv")
    seed_rows = read_csv(RESULTS / "ABLATION_SEED_METRICS.csv")
    mcv_rows = read_csv(RESULTS / "MCV_INFERENCE_ONLY_SUMMARY.csv")
    iteration_rows = read_csv(RESULTS / "ACTOR_ITERATION_SUMMARY.csv")

    summary = {(row["variant"], row["pde"]): row for row in summary_rows}
    mcv = {row["pde"]: row for row in mcv_rows}

    primary = [
        row
        for row in seed_rows
        if row["primary_replicate"].lower() == "true"
    ]

    def primary_by_seed(pde: str, variant: str) -> dict[int, float]:
        return {
            int(row["seed"]): float(row["value_rmse"])
            for row in primary
            if row["pde"] == pde and row["variant"] == variant
        }

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 0.58,
            "xtick.major.width": 0.55,
            "ytick.major.width": 0.55,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(7.35, 5.05), facecolor="white")
    grid = fig.add_gridspec(
        4,
        2,
        left=0.17,
        right=0.963,
        bottom=0.135,
        top=0.925,
        hspace=0.28,
        wspace=0.18,
        width_ratios=(1.22, 3.05),
    )
    abs_axes: list[plt.Axes] = []
    ratio_axes: list[plt.Axes] = []

    for pde_index, pde in enumerate(PDES):
        bottom = pde_index == len(PDES) - 1
        abs_ax = fig.add_subplot(grid[pde_index, 0])
        ratio_ax = fig.add_subplot(grid[pde_index, 1])
        abs_axes.append(abs_ax)
        ratio_axes.append(ratio_ax)
        soft_background(abs_ax)
        soft_background(ratio_ax)

        # Absolute RMSE bars.
        means = np.asarray(
            [float(summary[(variant, pde)]["rmse_mean"]) for variant in ABS_ORDER]
        )
        spreads = np.asarray(
            [float(summary[(variant, pde)]["rmse_std"]) for variant in ABS_ORDER]
        )
        abs_x = np.arange(4, dtype=float)
        absolute_bars = abs_ax.bar(
            abs_x,
            means,
            width=0.67,
            color=ABS_COLORS,
            edgecolor="#465763",
            linewidth=0.55,
            yerr=spreads,
            capsize=2.0,
            error_kw={
                "ecolor": "#26343E",
                "elinewidth": 0.72,
                "capthick": 0.72,
                "zorder": 7,
            },
            zorder=4,
        )
        for bar, hatch in zip(absolute_bars, ABS_HATCHES, strict=True):
            bar.set_hatch(hatch)
        abs_upper = float(np.max(means + spreads) * 1.30)
        abs_ax.set_xlim(-0.58, 3.58)
        abs_ax.set_ylim(0.0, abs_upper)
        abs_ax.set_xticks(abs_x)
        abs_ax.set_xticklabels(ABS_LABELS if bottom else (), fontsize=5.2)
        abs_ax.tick_params(axis="x", length=0, pad=3)
        abs_ax.set_yticks((0.0, abs_upper))
        abs_ax.set_yticklabels(("0", f"{abs_upper:.1e}"), fontsize=5.6)
        abs_ax.tick_params(axis="y", length=2.0, pad=1.4)
        abs_ax.grid(axis="y", which="major", color="white", linewidth=0.75,
                    alpha=0.95, zorder=1)
        for x_value, mean, spread in zip(abs_x, means, spreads, strict=True):
            abs_ax.annotate(
                sci(mean),
                xy=(x_value, mean + spread),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=4.25,
                color="#273641",
                zorder=9,
            )

        # Relative effects, now encoded entirely with filled bars.
        full = primary_by_seed(pde, "a0_full")
        anchor = primary_by_seed(pde, "a2_no_terminal_anchor")
        anchor_samples = np.asarray(
            [anchor[seed] / full[seed] for seed in sorted(full)], dtype=float
        )
        ratio_values = [
            float(mcv[pde]["plain_over_mcv_rmse_ratio_mean"]),
            float(mcv[pde]["without_training_over_full_ratio_mean"]),
            float(mcv[pde]["without_all_over_full_ratio_mean"]),
            float(np.mean(anchor_samples)),
        ]
        ratio_spreads = [
            float(mcv[pde]["plain_over_mcv_rmse_ratio_std"]),
            float(mcv[pde]["without_training_over_full_ratio_std"]),
            float(mcv[pde]["without_all_over_full_ratio_std"]),
            float(np.std(anchor_samples, ddof=1)),
        ]

        local_iterations = [row for row in iteration_rows if row["pde"] == pde]
        selected_by_seed = {
            int(row["seed"]): float(row["value_rmse"])
            for row in local_iterations
            if int(row["iteration"]) == int(row["selected_iteration"])
        }
        for iteration in range(5):
            current = {
                int(row["seed"]): float(row["value_rmse"])
                for row in local_iterations
                if int(row["iteration"]) == iteration
            }
            samples = np.asarray(
                [current[seed] / selected_by_seed[seed] for seed in sorted(selected_by_seed)],
                dtype=float,
            )
            ratio_values.append(float(np.mean(samples)))
            ratio_spreads.append(float(np.std(samples, ddof=1)))

        values = np.asarray(ratio_values)
        errors = np.asarray(ratio_spreads)
        for index, (x_value, value, spread, color, hatch) in enumerate(
            zip(RATIO_X, values, errors, RATIO_COLORS, RATIO_HATCHES, strict=True)
        ):
            ratio_ax.bar(
                x_value,
                max(value - RATIO_BAR_FLOOR, 0.0),
                bottom=RATIO_BAR_FLOOR,
                width=0.68 if index < 4 else 0.62,
                color=color,
                edgecolor="#465763",
                linewidth=0.58,
                hatch=hatch,
                yerr=spread,
                capsize=1.8,
                error_kw={
                    "ecolor": "#26343E",
                    "elinewidth": 0.68,
                    "capthick": 0.68,
                    "zorder": 8,
                },
                zorder=4,
            )
            ratio_ax.annotate(
                f"{value:.2f}×",
                xy=(x_value, value + spread),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=4.45,
                color="#273641",
                zorder=9,
            )

        ratio_ax.set_yscale("log")
        ratio_ax.set_ylim(0.68, 30.0)
        ratio_ax.set_xlim(-0.62, 10.62)
        ratio_ax.axhline(
            1.0,
            color="#527FA0",
            linewidth=0.75,
            linestyle=(0, (2, 2)),
            alpha=0.9,
            zorder=3,
        )
        ratio_ax.axvline(3.05, color="#98A5AE", linewidth=0.55, alpha=0.7)
        ratio_ax.axvline(5.05, color="#98A5AE", linewidth=0.55, alpha=0.7)
        ratio_ax.set_yticks((1.0, 3.0, 10.0, 30.0))
        ratio_ax.set_yticklabels(("1×", "3×", "10×", "30×"), fontsize=5.7)
        ratio_ax.yaxis.set_minor_formatter(NullFormatter())
        ratio_ax.tick_params(axis="y", which="major", length=2.0, pad=1.4)
        ratio_ax.tick_params(axis="y", which="minor", length=1.0)
        ratio_ax.set_xticks(RATIO_X)
        ratio_ax.set_xticklabels(RATIO_LABELS if bottom else (), fontsize=5.2)
        ratio_ax.tick_params(axis="x", length=0, pad=3)
        ratio_ax.grid(axis="y", which="major", color="#A7B2BA", linewidth=0.46,
                      alpha=0.64, zorder=1)
        ratio_ax.grid(axis="y", which="minor", color="#B9C1C7", linewidth=0.30,
                      alpha=0.42, linestyle=":", zorder=1)

        for axis in (abs_ax, ratio_ax):
            axis.set_axisbelow(True)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
            axis.spines["left"].set_color("#65737E")
            axis.spines["bottom"].set_color("#65737E")

        box = abs_ax.get_position()
        fig.text(
            0.012,
            (box.y0 + box.y1) / 2.0,
            PDE_LABELS[pde],
            ha="left",
            va="center",
            fontsize=6.25,
            fontweight="bold",
            linespacing=1.0,
            color="#223746",
        )

    abs_axes[0].set_title("(a)  Absolute RMSE", fontsize=8.2, pad=7,
                          color="#233B4C", fontweight="bold")
    ratio_axes[0].set_title("(b)  Relative component effects", fontsize=8.2, pad=7,
                            color="#233B4C", fontweight="bold")
    abs_axes[0].text(0.03, 0.94, "RMSE", transform=abs_axes[0].transAxes,
                     ha="left", va="top", fontsize=5.0, color="#65717A")
    ratio_axes[0].text(
        0.02,
        0.94,
        r"$\mathrm{RMSE}/\mathrm{RMSE}_{\mathrm{Full}}$",
        transform=ratio_axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=5.0,
        color="#65717A",
    )

    # One quiet bracket identifies the shared monotone Actor trajectory.
    actor_ax = ratio_axes[0]
    actor_ax.plot(
        (5.68, 5.68, 10.32, 10.32),
        (20.5, 24.0, 24.0, 20.5),
        color="#52636F",
        linewidth=0.72,
        clip_on=False,
        zorder=10,
    )
    actor_ax.text(
        8.0,
        23.75,
        "strictly decreasing",
        ha="center",
        va="bottom",
        fontsize=6.8,
        color="#43545F",
    )

    # The right-side bracket summarizes the common direction of every true
    # component-removal diagnostic; A4 remains the 1x complete-model anchor.
    top_box = ratio_axes[0].get_position()
    bottom_box = ratio_axes[-1].get_position()
    bracket_x = 0.974
    bracket_tick = 0.006
    bracket_top = top_box.y1
    bracket_bottom = bottom_box.y0
    fig.add_artist(
        Line2D(
            (
                bracket_x - bracket_tick,
                bracket_x,
                bracket_x,
                bracket_x - bracket_tick,
            ),
            (bracket_top, bracket_top, bracket_bottom, bracket_bottom),
            transform=fig.transFigure,
            color="#52636F",
            linewidth=0.78,
            clip_on=False,
        )
    )
    fig.text(
        0.986,
        0.5 * (bracket_top + bracket_bottom),
        "all component removals > 1×",
        ha="center",
        va="center",
        rotation=270,
        fontsize=7.0,
        color="#43545F",
    )

    legend_handles = (
        Patch(facecolor=MCV_COLORS[-1], edgecolor="#465763", label="MCV removal (light→dark)"),
        Patch(facecolor=ANCHOR_COLOR, edgecolor="#465763", label="terminal anchoring"),
        Patch(facecolor=ACTOR_UPDATE_COLORS[-1], edgecolor="#465763", label="actor updates A₀→A₄ (light→dark)"),
    )
    legend = fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.61, 0.025),
        ncol=3,
        frameon=True,
        fancybox=False,
        edgecolor="#C4CBD1",
        facecolor="white",
        fontsize=5.8,
        handlelength=1.5,
        columnspacing=1.2,
    )
    legend.get_frame().set_linewidth(0.5)

    FIGURES.mkdir(parents=True, exist_ok=True)
    pdf = FIGURES / "component_ablation_grouped_bars_preview.pdf"
    png = FIGURES / "component_ablation_grouped_bars_preview.png"
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {pdf}")
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
