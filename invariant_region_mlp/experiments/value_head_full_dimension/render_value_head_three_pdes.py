#!/usr/bin/env python3
"""Render the Allen-/Oscillatory-free direct-value-head comparison."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


REPO = Path(__file__).resolve().parents[3]
RESULTS = REPO / "invariant_region_mlp/results/value_head_full_dimension"
CSV_PATH = RESULTS / "VALUE_HEAD_FULL_DIMENSION.csv"
TABLE_DIR = RESULTS / "tables"

PDES = (
    ("joint_concave", "Concave Cubic–Quadratic PDE"),
    ("joint_convex", "Convex Cubic–Quadratic PDE"),
    ("dense_fourier", "Fourier–Tanh Coupled PDE"),
)
DIMS = (20, 50, 100, 200, 400)
METHODS = (
    ("Direct value head", "direct_value_head_rmse"),
    ("DPI", "dpi_rmse"),
    ("PINN-HTE", "pinn_hte_rmse"),
    ("SDGD-PINN", "sdgd_pinn_rmse"),
)


def sci(value: float) -> str:
    return f"{value:.2e}".replace("e-0", "e-").replace("e+0", "e+")


def render() -> None:
    with CSV_PATH.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    lookup = {(row["pde"], int(row["dimension"])): row for row in rows}

    values = {}
    winner = {}
    runner_up = {}
    for pde, _ in PDES:
        for dim in DIMS:
            row = lookup[(pde, dim)]
            cell = [float(row[key]) for _, key in METHODS]
            values[(pde, dim)] = cell
            ranking = sorted(range(len(cell)), key=cell.__getitem__)
            winner[(pde, dim)] = ranking[0]
            runner_up[(pde, dim)] = ranking[1]

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(14.2, 3.25))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    left, right = 0.022, 0.992
    method_w = 0.125
    data_left = left + method_w
    cell_w = (right - data_left) / 15
    top = 0.80
    group_h, dim_h, row_h = 0.105, 0.088, 0.118
    rows_top = top - group_h - dim_h
    bottom = rows_top - len(METHODS) * row_h

    navy = "#17324D"
    text = "#1D2630"
    muted = "#65717D"
    header_fill = "#F4F6F8"
    head_fill = "#E4F0F7"
    grid = "#D9DEE3"
    group_rule = "#9EA8B2"

    ax.text(left, 0.952, "Direct value head across three PDEs",
            ha="left", va="center", fontsize=15.0,
            fontweight="bold", color=navy)
    ax.text(left, 0.895,
            "Whole-function value RMSE on matched frozen state banks  ·  lower is better",
            ha="left", va="center", fontsize=8.4, color=muted)

    ax.add_patch(Rectangle((left, top - group_h - dim_h),
                           right - left, group_h + dim_h,
                           facecolor=header_fill, edgecolor="none"))
    ax.add_patch(Rectangle((left, rows_top - row_h),
                           right - left, row_h,
                           facecolor=head_fill, edgecolor="none"))

    ax.text(left + method_w / 2, top - (group_h + dim_h) / 2,
            "Method", ha="center", va="center", fontsize=8.5,
            fontweight="bold", color=text)
    for pde_i, (_, label) in enumerate(PDES):
        x0 = data_left + pde_i * 5 * cell_w
        ax.text(x0 + 2.5 * cell_w, top - group_h / 2, label,
                ha="center", va="center", fontsize=9.0,
                fontweight="bold", color=navy)
        for dim_i, dim in enumerate(DIMS):
            ax.text(x0 + (dim_i + 0.5) * cell_w,
                    top - group_h - dim_h / 2,
                    str(dim), ha="center", va="center",
                    fontsize=7.8, color=text)

    for method_i, (method, _) in enumerate(METHODS):
        y = rows_top - (method_i + 0.5) * row_h
        ax.text(left + method_w / 2, y, method,
                ha="center", va="center", fontsize=8.0,
                fontweight="bold" if method_i == 0 else "normal",
                color=navy if method_i == 0 else text)
        for pde_i, (pde, _) in enumerate(PDES):
            for dim_i, dim in enumerate(DIMS):
                val = values[(pde, dim)][method_i]
                x = data_left + (pde_i * 5 + dim_i + 0.5) * cell_w
                is_best = winner[(pde, dim)] == method_i
                is_second = runner_up[(pde, dim)] == method_i
                ax.text(x, y, sci(val), ha="center", va="center",
                        fontsize=7.35,
                        fontweight="bold" if is_best else "normal",
                        color=navy if is_best or method_i == 0 else text)
                if is_second:
                    ax.plot(
                        [x - 0.34 * cell_w, x + 0.34 * cell_w],
                        [y - 0.020, y - 0.020],
                        color=navy if method_i == 0 else text,
                        lw=0.65,
                        solid_capstyle="butt",
                    )

    ax.plot([left, right], [top, top], color=navy, lw=1.15)
    ax.plot([left, right], [top - group_h, top - group_h], color=grid, lw=0.55)
    ax.plot([left, right], [rows_top, rows_top], color=navy, lw=0.85)
    for row_idx in range(1, len(METHODS) + 1):
        y = rows_top - row_idx * row_h
        ax.plot([left, right], [y, y],
                color=group_rule if row_idx == len(METHODS) else grid,
                lw=0.75 if row_idx == len(METHODS) else 0.45)
    ax.plot([data_left, data_left], [bottom, top], color=group_rule, lw=0.65)
    for col_idx in range(1, 15):
        x = data_left + col_idx * cell_w
        strong = col_idx in (5, 10)
        ax.plot([x, x], [bottom, top],
                color=group_rule if strong else grid,
                lw=1.05 if strong else 0.35)

    ax.text(left, 0.070,
            "Direct-head pairwise wins: 35/45 (77.8%)  ·  vs DPI 7/15  ·  vs PINN-HTE 13/15  ·  vs SDGD-PINN 15/15",
            ha="left", va="center", fontsize=7.8,
            fontweight="bold", color="#334858")
    ax.text(right, 0.070, "Bold = best  ·  underline = second best",
            ha="right", va="center", fontsize=7.2, color=muted)

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(TABLE_DIR / "value_head_three_pdes.pdf",
                bbox_inches="tight", pad_inches=0.04)
    fig.savefig(TABLE_DIR / "value_head_three_pdes.png", dpi=300,
                bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


if __name__ == "__main__":
    render()
