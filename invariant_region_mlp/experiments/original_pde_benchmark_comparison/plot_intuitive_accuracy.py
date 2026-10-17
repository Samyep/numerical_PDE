"""Draw a compact, task-matched accuracy comparison for the final report."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, stdev

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np


EXPERIMENT_DIR = Path(__file__).resolve().parent
INPUT_CSV = EXPERIMENT_DIR / "function_approximation_comparison.csv"
FIGURE_DIR = EXPERIMENT_DIR / "figures"


@dataclass(frozen=True)
class Point:
    label: str
    mean: float
    sd: float | None
    kind: str


def _rows() -> list[dict[str, str]]:
    with INPUT_CSV.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _values(
    rows: list[dict[str, str]],
    *,
    benchmark: str,
    method: str,
    role: str,
    variant: str | None = None,
    time_steps: int | None = None,
) -> list[float]:
    selected: list[float] = []
    for row in rows:
        if row["benchmark"] != benchmark or row["method"] != method:
            continue
        if row["evaluation_role"] != role:
            continue
        if variant is not None and row["variant"] != variant:
            continue
        if time_steps is not None and int(row["time_steps"] or -1) != time_steps:
            continue
        selected.append(float(row["value_rmse"]))
    if not selected:
        raise RuntimeError(
            f"missing accuracy rows for {benchmark=} {method=} {role=} {variant=}"
        )
    return selected


def _point(label: str, values: list[float], kind: str) -> Point:
    return Point(
        label=label,
        mean=mean(values),
        sd=stdev(values) if len(values) > 1 else None,
        kind=kind,
    )


def _panels(rows: list[dict[str, str]]) -> list[tuple[str, list[Point], str]]:
    actor = "terminal_aware_tangent_actor_mcv"
    full_history = "full_history_mlp_raw"
    return [
        (
            "HJB-Log (100D)",
            [
                _point(
                    "Ours: Tangent-Actor + MCV",
                    _values(
                        rows,
                        benchmark="hjb_log_d100",
                        method=actor,
                        role="frozen_headline",
                        variant="focused/a4",
                        time_steps=64,
                    ),
                    "ours",
                ),
                _point(
                    "Full-History MLP (3,24)",
                    _values(
                        rows,
                        benchmark="hjb_log_d100",
                        method=full_history,
                        role="frozen_control_grid",
                        variant="n3/M24",
                    ),
                    "full_history",
                ),
                _point(
                    "SDGD-PINN",
                    _values(
                        rows,
                        benchmark="hjb_log_d100",
                        method="sdgd_pinn_official",
                        role="official_fixed",
                    ),
                    "official",
                ),
            ],
            "20.5% lower RMSE than Full-History\n35.1× lower than SDGD-PINN",
        ),
        (
            "OU / Gaussian-mixture HJB (100D)",
            [
                _point(
                    "Ours: Tangent-Actor + MCV",
                    _values(
                        rows,
                        benchmark="dpi_ou_gmm_d100",
                        method=actor,
                        role="frozen_headline",
                        variant="focused/a4",
                        time_steps=64,
                    ),
                    "ours",
                ),
                _point(
                    "Deep Picard Iteration",
                    _values(
                        rows,
                        benchmark="dpi_ou_gmm_d100",
                        method="dpi",
                        role="official_fixed",
                    ),
                    "official",
                ),
                _point(
                    "D-DBSDE",
                    _values(
                        rows,
                        benchmark="dpi_ou_gmm_d100",
                        method="d_dbsde",
                        role="official_fixed",
                    ),
                    "official",
                ),
                _point(
                    "PINN-HTE",
                    _values(
                        rows,
                        benchmark="dpi_ou_gmm_d100",
                        method="pinn_hte",
                        role="official_fixed",
                    ),
                    "official",
                ),
            ],
            "79.9× lower RMSE than Deep Picard\n(best official neural baseline)",
        ),
        (
            "Allen–Cahn (100D)",
            [
                _point(
                    "Ours: Tangent-Actor",
                    _values(
                        rows,
                        benchmark="allen_cahn_d100",
                        method=actor,
                        role="frozen_headline",
                        variant="focused/a4",
                        time_steps=64,
                    ),
                    "ours",
                ),
                _point(
                    "Full-History MLP (3,24)",
                    _values(
                        rows,
                        benchmark="allen_cahn_d100",
                        method=full_history,
                        role="frozen_control_grid",
                        variant="n3/M24",
                    ),
                    "full_history",
                ),
            ],
            "18.8% lower RMSE than Full-History\n(reference-limited regime)",
        ),
    ]


def main() -> None:
    panels = _panels(_rows())
    colors = {
        "ours": "#D1495B",
        "full_history": "#3A6EA5",
        "official": "#707782",
    }
    markers = {"ours": "*", "full_history": "o", "official": "s"}
    sizes = {"ours": 18, "full_history": 10, "official": 9}

    fig, axes = plt.subplots(1, 3, figsize=(16.8, 6.6), constrained_layout=False)
    fig.subplots_adjust(left=0.16, right=0.985, top=0.78, bottom=0.31, wspace=0.66)
    fig.suptitle(
        "Value accuracy on original 100D PDE benchmarks",
        fontsize=19,
        fontweight="bold",
        y=0.96,
    )
    fig.text(
        0.5,
        0.895,
        "Held-out Value RMSE on matched states within each panel  •  farther left is more accurate",
        ha="center",
        fontsize=12.5,
        color="#374151",
    )

    for axis, (title, points, callout) in zip(axes, panels, strict=True):
        y = np.arange(len(points))[::-1]
        values = np.asarray([point.mean for point in points])
        low = float(values.min() / 1.85)
        high = float(values.max() * 3.0)

        axis.axhspan(y[0] - 0.43, y[0] + 0.43, color="#FDECEF", zorder=0)
        for yi, point in zip(y, points, strict=True):
            xerr = point.sd if point.sd is not None else None
            axis.errorbar(
                point.mean,
                yi,
                xerr=xerr,
                fmt=markers[point.kind],
                markersize=sizes[point.kind],
                markerfacecolor=colors[point.kind],
                markeredgecolor="white" if point.kind == "ours" else colors[point.kind],
                markeredgewidth=1.1,
                ecolor=colors[point.kind],
                elinewidth=1.8,
                capsize=4 if xerr is not None else 0,
                zorder=3,
            )
            value_label = f"{point.mean:.3g}"
            if point.sd is not None:
                value_label += f" ± {point.sd:.2g}"
            axis.annotate(
                value_label,
                xy=(point.mean, yi),
                xytext=(9, 0),
                textcoords="offset points",
                va="center",
                fontsize=10.5,
                fontweight="bold" if point.kind == "ours" else "normal",
                color=colors[point.kind],
            )

        axis.set_xscale("log")
        axis.set_xlim(low, high)
        axis.set_ylim(-0.65, len(points) - 0.35)
        axis.set_yticks(y, [point.label for point in points])
        axis.set_title(title, fontsize=14, fontweight="bold", pad=18)
        axis.set_xlabel("Value RMSE  (log scale, lower is better)", fontsize=10.5)
        axis.grid(axis="x", which="both", alpha=0.22, linewidth=0.8)
        axis.xaxis.set_minor_formatter(mticker.NullFormatter())
        if title.startswith("Allen"):
            axis.xaxis.set_major_locator(mticker.FixedLocator([1e-4, 2e-4, 5e-4]))
            axis.xaxis.set_major_formatter(mticker.LogFormatterSciNotation())
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(axis="y", length=0, labelsize=10.5)
        axis.tick_params(axis="x", labelsize=9.5)
        axis.text(
            0.5,
            -0.24,
            callout,
            transform=axis.transAxes,
            ha="center",
            va="top",
            fontsize=10.5,
            color="#1F2937",
            bbox={
                "boxstyle": "round,pad=0.38",
                "facecolor": "#F3F4F6",
                "edgecolor": "#D1D5DB",
            },
        )

    fig.text(
        0.5,
        0.035,
        "Error bars: sample SD over 3 seeds; official one-run neural baselines have none.  "
        "Full-History: lowest-RMSE member of the frozen grid.\n"
        "Deep BSDE is omitted because its published task is single-point estimation rather than function approximation.",
        ha="center",
        va="bottom",
        fontsize=9.2,
        color="#4B5563",
    )

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    png = FIGURE_DIR / "intuitive_accuracy_comparison.png"
    pdf = FIGURE_DIR / "intuitive_accuracy_comparison.pdf"
    fig.savefig(png, dpi=240, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
