#!/usr/bin/env python3
"""Render the 25-cell direct-value-head comparison as PDF and PNG."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


REPO = Path(__file__).resolve().parents[3]
RESULTS = REPO / "invariant_region_mlp/results/value_head_full_dimension"
CSV_PATH = RESULTS / "VALUE_HEAD_FULL_DIMENSION.csv"
TABLE_DIR = RESULTS / "tables"

PDE_ORDER = (
    "allen_cahn",
    "joint_concave",
    "joint_convex",
    "nonconvex_coupled",
    "dense_fourier",
)
PDE_LABELS = {
    "allen_cahn": "Allen–Cahn",
    "joint_concave": "Concave Cubic–Quadratic",
    "joint_convex": "Convex Cubic–Quadratic",
    "nonconvex_coupled": "Oscillatory Value–Gradient",
    "dense_fourier": "Fourier–Tanh Coupled",
}
DIMENSIONS = (20, 50, 100, 200, 400)


def sci(value):
    return "{:.2e}".format(float(value)).replace("e-0", "e-").replace("e+0", "e+")


def main():
    with open(CSV_PATH, newline="") as handle:
        rows = list(csv.DictReader(handle))
    lookup = {(row["pde"], int(row["dimension"])): row for row in rows}

    columns = [
        "PDE",
        "d",
        "Direct value head",
        "DPI",
        "PINN-HTE",
        "SDGD-PINN",
        "Head wins",
    ]
    body = []
    winners = []
    for pde in PDE_ORDER:
        for index, dimension in enumerate(DIMENSIONS):
            row = lookup[(pde, dimension)]
            values = [
                float(row["direct_value_head_rmse"]),
                float(row["dpi_rmse"]),
                float(row["pinn_hte_rmse"]),
                float(row["sdgd_pinn_rmse"]),
            ]
            winners.append(values.index(min(values)) + 2)
            body.append(
                [
                    PDE_LABELS[pde] if index == 0 else "",
                    str(dimension),
                    sci(values[0]),
                    sci(values[1]),
                    sci(values[2]),
                    sci(values[3]),
                    row["direct_head_pairwise_wins"] + "/3",
                ]
            )

    fig = plt.figure(figsize=(12.4, 10.0), facecolor="white")
    ax = fig.add_axes([0.035, 0.075, 0.93, 0.865])
    ax.axis("off")
    ax.text(
        0.0,
        1.035,
        "Direct Tangent-Actor value head across dimensions",
        transform=ax.transAxes,
        fontsize=15.5,
        fontweight="semibold",
        va="bottom",
        color="#17212b",
    )
    ax.text(
        0.0,
        1.003,
        "Whole-function value RMSE on matched frozen state banks · lower is better",
        transform=ax.transAxes,
        fontsize=9.2,
        va="bottom",
        color="#5d6874",
    )

    table = ax.table(
        cellText=body,
        colLabels=columns,
        colWidths=[0.245, 0.055, 0.16, 0.125, 0.14, 0.145, 0.09],
        cellLoc="center",
        loc="upper left",
        bbox=[0.0, 0.045, 1.0, 0.945],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)

    header_bg = "#23364d"
    direct_bg = "#edf5fa"
    block_bg = ("#ffffff", "#f7f8fa")
    edge = "#d8dde3"
    for col in range(len(columns)):
        cell = table[(0, col)]
        cell.set_facecolor(header_bg)
        cell.set_text_props(color="white", fontweight="semibold")
        cell.set_edgecolor(header_bg)
        cell.set_linewidth(0.8)

    for row_index in range(1, len(body) + 1):
        pde_block = (row_index - 1) // len(DIMENSIONS)
        for col in range(len(columns)):
            cell = table[(row_index, col)]
            cell.set_facecolor(block_bg[pde_block % 2])
            cell.set_edgecolor(edge)
            cell.set_linewidth(0.45)
            if col == 0:
                cell.set_text_props(ha="left", color="#27313d", fontweight="semibold")
            else:
                cell.set_text_props(color="#27313d")
        table[(row_index, 2)].set_facecolor(direct_bg)
        table[(row_index, winners[row_index - 1])].set_text_props(
            fontweight="bold", color="#0b253a"
        )
        if (row_index - 1) % len(DIMENSIONS) == 0:
            for col in range(len(columns)):
                table[(row_index, col)].set_linewidth(0.9)
                table[(row_index, col)].set_edgecolor("#aeb7c1")

    ax.text(
        0.0,
        0.012,
        "Pairwise wins: 47/75 (62.7%)  ·  vs DPI 8/25  ·  vs PINN-HTE 15/25  ·  vs SDGD-PINN 24/25",
        transform=ax.transAxes,
        fontsize=9.2,
        fontweight="semibold",
        color="#263746",
        va="bottom",
    )
    ax.text(
        1.0,
        0.012,
        "Bold = lowest RMSE in row",
        transform=ax.transAxes,
        fontsize=8.2,
        color="#68737f",
        va="bottom",
        ha="right",
    )

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    pdf = TABLE_DIR / "value_head_full_dimension_preview.pdf"
    png = TABLE_DIR / "value_head_full_dimension_preview.png"
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", pdf)
    print("wrote", png)


if __name__ == "__main__":
    main()
