#!/usr/bin/env python3
"""Render a vector/raster preview of the three-PDE single-point table."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


OUT = Path(__file__).resolve().parent
DIMS = (20, 50, 100, 200, 400)
PDES = (
    "Full-Rank Log-Sum-Exp HJB",
    "Fourier–Tanh Coupled PDE",
    "Joint Nonseparable PDE",
)
METHODS = ("Tangent-Actor", "FH-MLP", "Deep BSDE")

# Each method stores RMSE followed by RelMAE.
DATA = {
    "Full-Rank Log-Sum-Exp HJB": {
        "Tangent-Actor": (
            ("7.13e-4", "1.42e-3", "1.59e-3", "4.17e-3", "3.63e-3"),
            ("7.84e-4", "8.76e-4", "9.04e-4", "1.33e-3", "1.25e-3"),
        ),
        "FH-MLP": (
            ("8.56e-2", "2.71e-1", "5.89e-1", "1.93e0", "3.19e0"),
            ("1.16e-1", "2.03e-1", "4.32e-1", "7.49e-1", "1.24e0"),
        ),
        "Deep BSDE": (
            ("1.13e-1", "2.53e-1", "2.82e-1", "3.72e-1", "3.24e-1"),
            ("1.58e-1", "1.94e-1", "2.11e-1", "1.49e-1", "1.30e-1"),
        ),
    },
    "Fourier–Tanh Coupled PDE": {
        "Tangent-Actor": (
            ("4.25e-3", "4.32e-3", "3.72e-3", "8.69e-3", "1.40e-2"),
            ("6.06e-3", "7.25e-3", "2.33e-3", "1.98e-2", "3.15e-2"),
        ),
        "FH-MLP": (
            ("3.60e-2", "4.23e-2", "1.33e-1", "1.43e-1", "2.93e-1"),
            ("4.98e-2", "7.21e-2", "1.07e-1", "3.56e-1", "6.37e-1"),
        ),
        "Deep BSDE": (
            ("4.38e-2", "7.79e-2", "1.13e-1", "7.52e-2", "7.91e-2"),
            ("6.95e-2", "1.45e-1", "9.37e-2", "1.94e-1", "1.73e-1"),
        ),
    },
    "Joint Nonseparable PDE": {
        "Tangent-Actor": (
            ("1.18e-3", "1.06e-3", "1.04e-3", "1.22e-3", "1.24e-3"),
            ("1.68e-3", "1.54e-3", "1.49e-3", "1.70e-3", "1.77e-3"),
        ),
        "FH-MLP": (
            ("1.14e-2", "2.56e-2", "5.15e-2", "1.02e-1", "2.13e-1"),
            ("1.46e-2", "3.44e-2", "7.39e-2", "1.41e-1", "3.00e-1"),
        ),
        "Deep BSDE": (
            ("2.59e-3", "3.76e-3", "3.48e-2", "NF", "NF"),
            ("3.35e-3", "4.46e-3", "1.77e-2", "NF", "NF"),
        ),
    },
}


def render() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(14.2, 4.35))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    left, right = 0.022, 0.992
    method_w = 0.105
    data_left = left + method_w
    cell_w = (right - data_left) / 15
    top = 0.862
    group_h, dim_h, row_h = 0.090, 0.076, 0.092
    rows_top = top - group_h - dim_h
    bottom = rows_top - 6 * row_h

    navy = "#17324D"
    ta_fill = "#E4F0F7"
    header_fill = "#F4F6F8"
    grid = "#D9DEE3"
    group_rule = "#9EA8B2"
    text = "#1D2630"

    ax.text(
        left,
        0.965,
        "Single-point accuracy across three PDEs",
        ha="left",
        va="center",
        fontsize=15.0,
        fontweight="bold",
        color=navy,
    )
    ax.text(
        left,
        0.913,
        "RMSE above RelMAE  ·  32 frozen query states  ·  bold = best valid value",
        ha="left",
        va="center",
        fontsize=8.4,
        color="#65717D",
    )

    ax.add_patch(Rectangle((left, top - group_h - dim_h), right - left, group_h + dim_h,
                           facecolor=header_fill, edgecolor="none"))
    ax.add_patch(Rectangle((left, rows_top - 2 * row_h), right - left, 2 * row_h,
                           facecolor=ta_fill, edgecolor="none"))

    ax.text(left + method_w / 2, top - (group_h + dim_h) / 2,
            "Method", ha="center", va="center", fontsize=8.5,
            fontweight="bold", color=text)

    for pde_i, pde in enumerate(PDES):
        x0 = data_left + pde_i * 5 * cell_w
        ax.text(x0 + 2.5 * cell_w, top - group_h / 2, pde,
                ha="center", va="center", fontsize=9.0,
                fontweight="bold", color=navy)
        for dim_i, dim in enumerate(DIMS):
            ax.text(x0 + (dim_i + 0.5) * cell_w, top - group_h - dim_h / 2,
                    str(dim), ha="center", va="center", fontsize=7.7,
                    color=text)

    for method_i, method in enumerate(METHODS):
        y_pair_top = rows_top - method_i * 2 * row_h
        ax.text(left + method_w / 2, y_pair_top - row_h,
                method, ha="center", va="center", fontsize=8.0,
                fontweight="bold" if method == "Tangent-Actor" else "normal",
                color=navy if method == "Tangent-Actor" else text)
        for metric_i in range(2):
            y = y_pair_top - (metric_i + 0.5) * row_h
            for pde_i, pde in enumerate(PDES):
                vals = DATA[pde][method][metric_i]
                for dim_i, val in enumerate(vals):
                    x = data_left + (pde_i * 5 + dim_i + 0.5) * cell_w
                    is_ta = method == "Tangent-Actor"
                    is_nf = val == "NF"
                    ax.text(
                        x,
                        y,
                        val,
                        ha="center",
                        va="center",
                        fontsize=7.25,
                        fontweight="bold" if is_ta else "normal",
                        fontstyle="italic" if is_nf else "normal",
                        color=navy if is_ta else ("#7E8790" if is_nf else text),
                    )

    # Rules: restrained, with stronger separators between PDE groups.
    ax.plot([left, right], [top, top], color=navy, lw=1.15)
    ax.plot([left, right], [top - group_h, top - group_h], color=grid, lw=0.55)
    ax.plot([left, right], [rows_top, rows_top], color=navy, lw=0.85)
    for row_idx in range(1, 7):
        y = rows_top - row_idx * row_h
        lw = 0.75 if row_idx in (2, 4, 6) else 0.35
        color = group_rule if row_idx in (2, 4, 6) else grid
        ax.plot([left, right], [y, y], color=color, lw=lw)
    ax.plot([left + method_w, left + method_w], [bottom, top], color=group_rule, lw=0.65)
    for col_idx in range(1, 15):
        x = data_left + col_idx * cell_w
        strong = col_idx in (5, 10)
        ax.plot([x, x], [bottom, top], color=group_rule if strong else grid,
                lw=1.05 if strong else 0.35)

    ax.text(
        right,
        0.047,
        "NF: numerical failure in the frozen Deep BSDE configuration.",
        ha="right",
        va="center",
        fontsize=7.2,
        color="#65717D",
    )

    fig.savefig(OUT / "single_point_three_pdes.pdf", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUT / "single_point_three_pdes.png", dpi=300,
                bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


if __name__ == "__main__":
    render()
