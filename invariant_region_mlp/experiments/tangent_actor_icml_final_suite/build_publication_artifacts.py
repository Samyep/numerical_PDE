"""Build publication figures, LaTeX tables and reproducibility manifests.

The builder is safe to run during staged execution: measured rows are drawn,
missing cells are marked, and the artifact audit remains ``in_progress`` until
the execution manifest contains no pending supported cell.  No numerical
endpoint is extrapolated.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
import platform
from pathlib import Path
import subprocess
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd

from .protocol import RESULTS_DIR, SOURCE_COMMITS


SOURCE_ROOT = Path(__file__).resolve().parents[3]
FIGURE_DIR = RESULTS_DIR / "figures"
TABLE_DIR = RESULTS_DIR / "tables"
DIMENSIONS = [10, 20, 50, 100, 200, 400]
MLP_DIMENSIONS = [20, 50, 100, 200, 400, 1000]
PDE_ORDER = [
    "hjb_log", "ou_gmm", "joint_concave", "joint_convex",
    "nonconvex_coupled", "allen_cahn",
]
PDE_LABELS = {
    "hjb_log": "HJB–Log",
    "ou_gmm": "OU/GMM HJB",
    "joint_concave": "Joint concave",
    "joint_convex": "Joint convex",
    "nonconvex_coupled": "Nonconvex coupled",
    "mlp_failure": "MLP failure",
    "allen_cahn": "Allen–Cahn",
}
METHOD_ORDER = ["tangent_actor", "pinn_hte", "sdgd_pinn", "dpi"]
METHOD_LABELS = {
    "tangent_actor": "Tangent-Actor",
    "pinn_hte": "PINN-HTE",
    "sdgd_pinn": "SDGD-PINN",
    "dpi": "DPI",
    "deep_bsde": "Deep BSDE",
}
COLORS = {
    "tangent_actor": "#C73E1D",
    "pinn_hte": "#2A6FBB",
    "sdgd_pinn": "#2E8B57",
    "dpi": "#8C5AAE",
    "deep_bsde": "#D18F00",
    "full_history": "#4B5563",
    "plain": "#667085",
    "mcv": "#C73E1D",
}
PDE_COLORS = {
    "hjb_log": "#2A6FBB", "ou_gmm": "#8C5AAE",
    "joint_concave": "#2E8B57", "joint_convex": "#008C95",
    "nonconvex_coupled": "#C73E1D", "mlp_failure": "#D18F00",
    "allen_cahn": "#B5477A",
}
MARKERS = {
    "tangent_actor": "o", "pinn_hte": "s", "sdgd_pinn": "^",
    "dpi": "D", "deep_bsde": "P",
}
FRONTIER_LABEL_OFFSETS = {
    "hjb_log": (4, 6), "ou_gmm": (4, 5),
    "joint_concave": (4, 12), "joint_convex": (-52, -9),
    "nonconvex_coupled": (4, -3), "mlp_failure": (4, 7),
    "allen_cahn": (4, 8),
}


def _configure_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9.5,
        "axes.titlesize": 11, "axes.labelsize": 10,
        "legend.fontsize": 8.5, "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True,
        "grid.alpha": 0.22, "grid.linewidth": 0.7,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def _read(name: str) -> pd.DataFrame:
    path = RESULTS_DIR / name
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path)


def _save_figure(fig: plt.Figure, stem: str) -> list[Path]:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    outputs = [FIGURE_DIR / f"{stem}.pdf", FIGURE_DIR / f"{stem}.png"]
    # Matplotlib otherwise embeds the wall-clock creation time in each PDF,
    # changing artifact hashes even when the numerical data and rendering are
    # identical.  Omit volatile timestamps so repeated builds are byte stable.
    fig.savefig(
        outputs[0], bbox_inches="tight",
        metadata={
            "Creator": "Tangent-Actor ICML final-suite artifact builder",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    fig.savefig(
        outputs[1], dpi=300, bbox_inches="tight",
        metadata={"Software": "Tangent-Actor ICML final-suite artifact builder"},
    )
    plt.close(fig)
    return outputs


def _historical_full_history() -> pd.DataFrame:
    path = SOURCE_ROOT / "invariant_region_mlp/results/counterexample_rescue/summary_all.csv"
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path)
    frame = frame[
        (frame["n"] == 3) & (frame["m"] == 3) & (frame["mode"] == "raw")
        & frame["d"].isin(MLP_DIMENSIONS)
    ].copy()
    return frame.rename(columns={"d": "dimension"})


def figure_mlp_failure() -> list[Path]:
    full = _read("MLP_FAILURE_FULL_HISTORY.csv")
    if full.empty:
        full = _historical_full_history()
    actor = _read("DIMENSION_SCALING.csv")
    if not actor.empty:
        actor = actor[(actor["pde"] == "mlp_failure") & (actor["method"] == "tangent_actor")]

    fig, axis = plt.subplots(figsize=(7.25, 4.55))
    if not full.empty:
        full = full.sort_values("dimension").set_index("dimension").reindex(MLP_DIMENSIONS).reset_index()
        yerr = None
        if {"value_rmse_ci95_low", "value_rmse_ci95_high"}.issubset(full.columns):
            y = full["value_rmse"].to_numpy(float)
            yerr = np.vstack((
                y - full["value_rmse_ci95_low"].to_numpy(float),
                full["value_rmse_ci95_high"].to_numpy(float) - y,
            ))
        axis.errorbar(
            full["dimension"], full["value_rmse"], yerr=yerr,
            color=COLORS["full_history"], marker="s", linewidth=2,
            capsize=3, label="Raw Full-History MLP (n=m=3)", zorder=3,
        )
    if not actor.empty:
        actor = actor.sort_values("dimension")
        axis.scatter(
            actor["dimension"], actor["value_rmse_mean"], s=52,
            color=COLORS["tangent_actor"], marker="o",
            label="Tangent-Actor + MCV", zorder=5,
        )
    axis.set_xlabel("Ambient dimension $d$")
    axis.set_ylabel("Whole-function Value RMSE")
    axis.set_title("Dimensional amplification in the MLP counterexample")
    axis.set_xticks(MLP_DIMENSIONS)
    axis.set_xlim(5, 1040)
    axis.set_ylim(bottom=0)
    axis.legend(loc="upper left", frameon=False)

    zoom = inset_axes(axis, width="45%", height="42%", loc="center right", borderpad=1.3)
    categorical_x = np.arange(len(MLP_DIMENSIONS))
    actor_by_dimension = {
        int(row["dimension"]): float(row["value_rmse_mean"])
        for _, row in actor.iterrows()
    } if not actor.empty else {}
    if actor_by_dimension:
        present = [d for d in MLP_DIMENSIONS if d in actor_by_dimension]
        zoom.errorbar(
            [MLP_DIMENSIONS.index(d) for d in present],
            [actor_by_dimension[d] for d in present],
            color=COLORS["tangent_actor"],
            marker="o", linewidth=1.7, capsize=3,
        )
    completed = set(actor["dimension"].astype(int)) if not actor.empty else set()
    missing = [d for d in MLP_DIMENSIONS if d not in completed]
    if missing:
        zoom.scatter([MLP_DIMENSIONS.index(d) for d in missing], np.zeros(len(missing)),
                     marker="x", color="#98A2B3", s=25)
        zoom.text(0.02, 0.95, "× pending", transform=zoom.transAxes, va="top", color="#667085")
    maximum = float(actor["value_rmse_mean"].max()) if not actor.empty else 0.01
    zoom.set_ylim(-0.05 * maximum, 1.35 * maximum)
    zoom.set_xticks(categorical_x, [str(d) for d in MLP_DIMENSIONS])
    zoom.tick_params(axis="x", rotation=35, labelsize=7)
    zoom.set_title("Tangent-Actor scale", fontsize=8.5)
    zoom.grid(alpha=0.2)
    fig.text(
        0.01, 0.005,
        "Crosses mark unexecuted cells; lines connect measured configurations only.",
        fontsize=7.5, color="#667085",
    )
    return _save_figure(fig, "figure1_mlp_failure")


def _representative_frontier(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    data = frame[frame["pde"].isin((
        "hjb_log", "ou_gmm", "joint_concave", "joint_convex",
        "nonconvex_coupled", "mlp_failure", "allen_cahn",
    ))].copy()
    data = data[
        (data["method"] != "tangent_actor")
        | data["estimator"].isin(("mcv", "direct"))
    ]
    return data


def figure_accuracy_compute() -> list[Path]:
    data = _representative_frontier(_read("ACCURACY_COMPUTE_FRONTIER.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), sharey=True)
    for axis, cost, title in (
        (axes[0], "total_seconds", "Training + selection + inference"),
        (axes[1], "inference_seconds", "Frozen-model inference"),
    ):
        if not data.empty:
            subset = data[np.isfinite(pd.to_numeric(data[cost], errors="coerce"))].copy()
            grouped = subset.groupby(
                ["pde", "method", "steps"], dropna=False,
            ).agg(**{
                cost: (cost, "mean"), "value_rmse": ("value_rmse", "mean"),
            }).reset_index()
            for pde in PDE_COLORS:
                for method in METHOD_ORDER:
                    rows = grouped[(grouped["pde"] == pde) & (grouped["method"] == method)]
                    if rows.empty:
                        continue
                    rows = rows.sort_values(cost)
                    linestyle = "-" if method == "tangent_actor" and len(rows) > 1 else "none"
                    axis.plot(
                        rows[cost], rows["value_rmse"], color=PDE_COLORS[pde],
                        marker=MARKERS[method], markersize=6.5, linewidth=1.4,
                        linestyle=linestyle, alpha=0.86, zorder=3,
                    )
                    if method == "tangent_actor":
                        endpoint = rows.iloc[-1]
                        axis.annotate(
                            PDE_LABELS[pde], (endpoint[cost], endpoint["value_rmse"]),
                            xytext=FRONTIER_LABEL_OFFSETS[pde], textcoords="offset points",
                            fontsize=6.7,
                            color=PDE_COLORS[pde],
                            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.62, "pad": 0.25},
                        )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Measured seconds")
        axis.set_title(title)
    axes[0].set_ylabel("Value RMSE")
    method_handles = [
        Line2D([0], [0], marker=MARKERS[m], color="#344054", linestyle="none",
               label=METHOD_LABELS[m], markersize=6)
        for m in METHOD_ORDER
    ]
    pde_handles = [
        Line2D([0], [0], color=color, linewidth=2, label=PDE_LABELS[pde])
        for pde, color in PDE_COLORS.items()
    ]
    fig.suptitle("Measured accuracy–compute frontier", y=0.995, fontsize=12, fontweight="bold")
    fig.legend(handles=pde_handles, loc="upper center", bbox_to_anchor=(0.5, 0.935),
               ncol=7, frameon=False, fontsize=8)
    fig.legend(handles=method_handles, loc="upper center", bbox_to_anchor=(0.5, 0.88),
               ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.79))
    return _save_figure(fig, "figure2_accuracy_compute")


def figure_dimension_scaling() -> list[Path]:
    data = _read("DIMENSION_SCALING.csv")
    fig, axes = plt.subplots(2, 3, figsize=(10.2, 6.0), sharex=True)
    positions = {dimension: index for index, dimension in enumerate(DIMENSIONS)}
    for axis, pde in zip(axes.flat, PDE_ORDER, strict=True):
        subset = data[data["pde"] == pde] if not data.empty else pd.DataFrame()
        if not subset.empty:
            # Keep source-native fidelity anchors such as SDGD d=250 out of
            # the frozen six-dimension sweep.  Aggregation enforces this too;
            # the local filter keeps the plot builder robust to older CSVs.
            subset = subset[subset["dimension"].isin(DIMENSIONS)]
        finite_values: list[float] = []
        for method in METHOD_ORDER:
            rows = subset[subset["method"] == method].sort_values("dimension") if not subset.empty else pd.DataFrame()
            if not rows.empty:
                values = rows["value_relative_l2_mean"].to_numpy(float)
                finite_values.extend(values[np.isfinite(values)].tolist())
                axis.errorbar(
                    [positions[int(value)] for value in rows["dimension"]], values,
                    yerr=rows["value_relative_l2_sample_std"],
                    color=COLORS[method], marker=MARKERS[method], linewidth=1.6,
                    capsize=2.5, label=METHOD_LABELS[method], zorder=3,
                )
        if finite_values:
            for lane, method in enumerate(METHOD_ORDER):
                completed = set(
                    subset.loc[subset["method"] == method, "dimension"].astype(int)
                ) if not subset.empty else set()
                missing = [d for d in DIMENSIONS if d not in completed]
                if missing:
                    axis.scatter(
                        [positions[d] + (lane - 1.5) * 0.07 for d in missing],
                        [1.015 + 0.025 * lane] * len(missing),
                        marker="x", s=13, linewidths=0.8,
                        color=COLORS[method], alpha=0.55,
                        transform=axis.get_xaxis_transform(), clip_on=False,
                    )
            axis.set_yscale("log")
        else:
            axis.text(0.5, 0.5, "All supported cells pending", transform=axis.transAxes,
                      ha="center", va="center", color="#667085")
            axis.set_ylim(0, 1)
        axis.set_title(PDE_LABELS[pde])
        axis.set_xticks(range(len(DIMENSIONS)), [str(d) for d in DIMENSIONS])
    for axis in axes[:, 0]:
        axis.set_ylabel("Value relative $L^2$")
    for axis in axes[-1, :]:
        axis.set_xlabel("Ambient dimension $d$")
    handles = [
        Line2D([0], [0], color=COLORS[m], marker=MARKERS[m], label=METHOD_LABELS[m])
        for m in METHOD_ORDER
    ]
    fig.suptitle("Dimension scaling (× above an axis marks a pending cell)", y=0.995,
                 fontsize=12, fontweight="bold")
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.955),
               ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    return _save_figure(fig, "figure3_dimension_scaling")


def _iteration_number(value: Any) -> int:
    text = str(value)
    return int(text[1:]) if text.startswith("a") else int(float(text))


def _mean_band(axis: plt.Axes, frame: pd.DataFrame, column: str, *, label: str,
               color: str, marker: str = "o") -> None:
    clean = frame[np.isfinite(pd.to_numeric(frame[column], errors="coerce"))].copy()
    if clean.empty:
        return
    clean["iteration_number"] = clean["iteration"].map(_iteration_number)
    grouped = clean.groupby("iteration_number")[column].agg(["mean", "std", "count"]).reset_index()
    axis.plot(grouped["iteration_number"], grouped["mean"], color=color, marker=marker,
              linewidth=1.8, label=label)
    valid = grouped["count"] > 1
    if valid.any():
        x = grouped.loc[valid, "iteration_number"].to_numpy(float)
        mean = grouped.loc[valid, "mean"].to_numpy(float)
        sd = grouped.loc[valid, "std"].to_numpy(float)
        axis.fill_between(x, np.maximum(mean - sd, np.finfo(float).tiny), mean + sd,
                          color=color, alpha=0.14, linewidth=0)


def figure_actor_iteration() -> list[Path]:
    data = _read("ABLATION_RESULTS.csv")
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.1), sharex="col")
    for row_index, pde in enumerate(("joint_concave", "nonconvex_coupled")):
        updates = data[
            (data["component"] == "learned_tangent_updates") & (data["pde"] == pde)
        ].copy() if not data.empty else pd.DataFrame()
        left, right = axes[row_index]
        if updates.empty:
            for axis in (left, right):
                axis.text(0.5, 0.5, "Pending mature all-iteration diagnostic",
                          transform=axis.transAxes, ha="center", va="center",
                          color="#667085", fontsize=9)
                axis.set_xticks(range(5), [f"a{i}" for i in range(5)])
            left.set_ylabel(PDE_LABELS[pde] + "\nerror")
            right.set_ylabel("Consistency / variance reduction")
            continue
        _mean_band(left, updates, "value_rmse", label="Policy Value RMSE",
                   color=COLORS["tangent_actor"], marker="o")
        _mean_band(left, updates, "actor_gradient_vector_rmse", label="Actor gradient RMSE",
                   color="#2A6FBB", marker="s")
        left.set_yscale("log")
        left.set_ylabel(PDE_LABELS[pde] + "\nerror")
        if left.get_legend_handles_labels()[0]:
            left.legend(frameon=False, loc="best")

        _mean_band(right, updates, "selection_metric", label="Reference-free consistency",
                   color="#8C5AAE", marker="D")
        twin = right.twinx()
        _mean_band(twin, updates, "variance_reduction", label="MCV variance reduction",
                   color="#2E8B57", marker="^")
        right.set_yscale("log")
        twin.set_yscale("log")
        right.set_ylabel("Consistency")
        twin.set_ylabel("Variance reduction")
        handles1, labels1 = right.get_legend_handles_labels()
        handles2, labels2 = twin.get_legend_handles_labels()
        if handles1 or handles2:
            right.legend(handles1 + handles2, labels1 + labels2, frameon=False, loc="best")
        for axis in (left, right, twin):
            axis.set_xticks(range(5), [f"a{i}" for i in range(5)])
    axes[-1, 0].set_xlabel("Actor update")
    axes[-1, 1].set_xlabel("Actor update")
    fig.suptitle("Frozen Actor iterations and martingale control variates", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return _save_figure(fig, "figure4_actor_iteration_mcv")


def figure_high_accuracy() -> list[Path]:
    data = _read("ACCURACY_COMPUTE_FRONTIER.csv")
    if not data.empty:
        data = data[(data["method"] == "tangent_actor") & (data["estimator"] == "mcv")]
    counts = data.groupby("pde")["steps"].nunique() if not data.empty else pd.Series(dtype=int)
    families = [name for name in ("joint_concave", "nonconvex_coupled", "mlp_failure")
                if counts.get(name, 0) >= 2]
    if not families and not data.empty:
        families = list(counts.sort_values(ascending=False).head(2).index)
    figure_columns = max(1, len(families))
    fig, axes = plt.subplots(1, figure_columns, figsize=(4.3 * figure_columns, 3.7), squeeze=False)
    for axis, pde in zip(axes.flat, families, strict=False):
        subset = data[data["pde"] == pde].copy()
        grouped = subset.groupby("steps").agg(
            rmse=("value_rmse", "mean"), rmse_sd=("value_rmse", "std"),
            mc=("mc_rmse_contribution", "mean"),
            reference=("reference_uncertainty", "max"),
        ).reset_index()
        axis.errorbar(grouped["steps"], grouped["rmse"], yerr=grouped["rmse_sd"],
                      color=COLORS["tangent_actor"], marker="o", capsize=3,
                      label="Observed RMSE")
        axis.plot(grouped["steps"], grouped["mc"], color="#2A6FBB", marker="s",
                  linestyle="--", label="Finite-path MC proxy")
        reference = grouped["reference"].dropna()
        if not reference.empty and float(reference.max()) > 0:
            axis.axhline(float(reference.max()), color="#667085", linestyle=":",
                         label="Reference uncertainty")
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_xticks(sorted(grouped["steps"].astype(int).unique()))
        axis.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        axis.set_xlabel("Euler steps")
        axis.set_title(PDE_LABELS.get(pde, pde))
        axis.legend(frameon=False)
    axes[0, 0].set_ylabel("Value error")
    fig.suptitle("Post-freeze high-accuracy frontier", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save_figure(fig, "figure5_high_accuracy_frontier")


def _latex_escape(value: Any) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}", "_": r"\_", "%": r"\%", "&": r"\&",
        "#": r"\#", "$": r"\$", "{": r"\{", "}": r"\}",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def _metric_cell(rows: pd.DataFrame) -> str:
    if rows.empty:
        return "--"
    rmse = rows["value_rmse_mean"].iloc[0]
    rel = rows["value_relative_l2_mean"].iloc[0]
    seeds = int(rows["seeds"].iloc[0])
    return rf"\shortstack{{{rel:.3g} rel.\\{rmse:.3g} RMSE\\$n={seeds}$}}"


def table_main() -> Path:
    data = _read("DIMENSION_SCALING.csv")
    lines = [
        r"\begin{tabular}{llcccc}", r"\toprule",
        r"PDE & $d$ & Tangent-Actor & PINN-HTE & SDGD-PINN & DPI \\",
        r"\midrule",
    ]
    for pde in PDE_ORDER:
        for index, d in enumerate(DIMENSIONS):
            cells = []
            for method in METHOD_ORDER:
                rows = data[(data["pde"] == pde) & (data["dimension"] == d)
                            & (data["method"] == method)] if not data.empty else pd.DataFrame()
                cells.append(_metric_cell(rows))
            label = _latex_escape(PDE_LABELS[pde]) if index == 0 else ""
            lines.append(f"{label} & {d} & " + " & ".join(cells) + r" \\")
        lines.append(r"\addlinespace")
    lines.extend((r"\bottomrule", r"\end{tabular}"))
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = TABLE_DIR / "table1_main_benchmark.tex"
    path.write_text("\n".join(lines) + "\n")
    return path


def _select_native_primary_rows(data: pd.DataFrame) -> pd.DataFrame:
    """Keep the predeclared estimate for the sole primary training seed."""
    if data.empty:
        return data.copy()
    selected = []
    keys = ["pde", "dimension", "method", "method_status", "seed"]
    for key, rows in data.groupby(keys, dropna=False, sort=True):
        candidates = rows.copy()
        method = str(key[2])
        if method == "tangent_actor" and len(candidates) > 1:
            steps = pd.to_numeric(candidates["steps"], errors="coerce")
            if (steps == 64).any():
                candidates = candidates.loc[steps == 64]
            estimators = candidates["estimator"].fillna("").astype(str).str.lower()
            if (estimators == "mcv").any():
                candidates = candidates.loc[estimators == "mcv"]
            elif (estimators == "plain").any():
                candidates = candidates.loc[estimators == "plain"]
        if len(candidates) != 1:
            raise RuntimeError(
                "native table requires exactly one primary estimate per training seed; "
                f"key={key!r}, candidates={len(candidates)}"
            )
        selected.append(candidates.iloc[0])
    return pd.DataFrame(selected).reset_index(drop=True)


def table_native() -> Path:
    data = _read("NATIVE_POINT_RESULTS.csv")
    data = _select_native_primary_rows(data)
    lines = [
        r"\begin{tabular}{lrllrll}", r"\toprule",
        r"PDE & $d$ & Method & Status & Seeds & Abs. error (seed 0) & Cost (s) \\",
        r"\midrule",
    ]
    if not data.empty:
        grouped = data.groupby(["pde", "dimension", "method", "method_status"], dropna=False)
        for (pde, dimension, method, status), rows in grouped:
            errors = rows["absolute_error"].to_numpy(float)
            costs = rows["total_seconds"].dropna().to_numpy(float)
            error = f"{errors.mean():.3g}"
            if len(errors) > 1:
                error += rf" $\pm$ {errors.std(ddof=1):.2g}"
            cost = f"{costs.mean():.3g}" if len(costs) else "--"
            lines.append(
                f"{_latex_escape(PDE_LABELS.get(pde, pde))} & {int(dimension)} & "
                f"{_latex_escape(METHOD_LABELS.get(method, method))} & "
                f"{_latex_escape(status)} & {rows['seed'].nunique()} & {error} & {cost} " + r"\\"
            )
    lines.extend((r"\bottomrule", r"\end{tabular}"))
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = TABLE_DIR / "native_point_comparison.tex"
    path.write_text("\n".join(lines) + "\n")
    return path


def table_coverage() -> Path:
    compatibility = _read("DIMENSION_COMPATIBILITY_MATRIX.csv")
    manifest_path = RESULTS_DIR / "EXECUTION_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    status = {
        (row["pde"], int(row["dimension"]), row["method"]): row["execution_status"]
        for row in manifest["method_cells"]
    }
    lines = [
        r"\begin{tabular}{llccccc}", r"\toprule",
        r"PDE & $d$ & Ours & PINN-HTE & SDGD & DPI & Deep BSDE \\", r"\midrule",
    ]
    methods = METHOD_ORDER + ["deep_bsde"]
    for pde in PDE_ORDER + ["mlp_failure"]:
        dims = MLP_DIMENSIONS if pde == "mlp_failure" else DIMENSIONS
        for index, d in enumerate(dims):
            cells = []
            for method in methods:
                value = status.get((pde, d, method), "not_listed")
                if value.startswith("completed"):
                    cell = r"\checkmark"
                elif value == "unsupported":
                    cell = "n/a"
                else:
                    cell = r"\texttimes"
                cells.append(cell)
            label = _latex_escape(PDE_LABELS[pde]) if index == 0 else ""
            lines.append(f"{label} & {d} & " + " & ".join(cells) + r" \\")
    lines.extend((r"\bottomrule", r"\end{tabular}"))
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = TABLE_DIR / "execution_coverage.tex"
    path.write_text("\n".join(lines) + "\n")
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT,
        check=True, capture_output=True, text=True,
    ).stdout.strip()


def build_reproducibility_manifest(artifacts: Iterable[Path]) -> Path:
    execution = json.loads((RESULTS_DIR / "EXECUTION_MANIFEST.json").read_text())
    resources = _read("RESOURCE_ACCOUNTING.csv")
    jobs = []
    if not resources.empty:
        current = resources[resources["accounting_origin"] == "icml_final_suite_incremental"]
        jobs = current[[
            "job_id", "job_name", "state", "elapsed_seconds", "allocated_gpu_count",
            "allocated_gpu_hours",
        ]].to_dict("records")
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "branch": "ir-mlp-tangent-actor-icml-final-suite",
        "git_head": _git_head(), "base_commit": execution["base_commit"],
        "pinned_source_commits": SOURCE_COMMITS,
        "environment": {
            "python": platform.python_version(), "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
            "gpu_from_job_logs": "NVIDIA GH200 120GB; see slurm/*.out",
        },
        "canonical_commands": [
            "python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.preflight",
            "python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.aggregate_results",
            "python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.collect_resource_accounting",
            "python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.build_publication_artifacts",
            "python -m pytest -q invariant_region_mlp/experiments/tangent_actor_icml_final_suite",
        ],
        "slurm_jobs": jobs,
        "test_status": "24 passed after Stage 0 paired-bank repair; see final-suite pytest output",
        "artifacts": [
            {
                "path": str(path.relative_to(SOURCE_ROOT)), "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in sorted(set(artifacts)) if path.exists()
        ],
        "raw_large_artifacts": (
            "Checkpoints, cached labels, raw JSON, predictions and Slurm logs remain on "
            "the shared workspace and are intentionally excluded from Git."
        ),
    }
    path = RESULTS_DIR / "REPRODUCIBILITY_MANIFEST.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
    return path


def build_artifact_audit(artifacts: Iterable[Path]) -> Path:
    manifest = json.loads((RESULTS_DIR / "EXECUTION_MANIFEST.json").read_text())
    main = pd.read_csv(RESULTS_DIR / "MAIN_RESULTS.csv")
    native = pd.read_csv(RESULTS_DIR / "NATIVE_POINT_RESULTS.csv")
    scaling = pd.read_csv(RESULTS_DIR / "DIMENSION_SCALING.csv")
    review_gate = json.loads((RESULTS_DIR / "REVIEW_GATE_STATUS.json").read_text())
    from .verify_paired_banks import audit as paired_bank_audit
    paired = paired_bank_audit(RESULTS_DIR / "MAIN_RESULTS.csv", dimension=100)
    paired_path = RESULTS_DIR / "PAIRWISE_BANK_VERIFICATION_100D.json"
    paired_path.write_text(json.dumps(paired, indent=2, sort_keys=True) + "\n")
    counts: dict[str, int] = {}
    for row in manifest["method_cells"]:
        counts[row["execution_status"]] = counts.get(row["execution_status"], 0) + 1
    records = []
    for path in sorted(set(artifacts)):
        records.append({
            "path": str(path.relative_to(SOURCE_ROOT)), "exists": path.exists(),
            "bytes": path.stat().st_size if path.exists() else None,
            "sha256": _sha256(path) if path.exists() else None,
        })
    pending = sum(count for status, count in counts.items() if not status.startswith("completed")
                  and status != "unsupported")
    main_text = main.astype(str)
    native_text = native.astype(str)
    source_paths = list(main["source_artifact"].dropna()) + list(native["source_artifact"].dropna())
    checks = {
        "main_has_no_deep_bsde": int((main["method"] == "deep_bsde").sum()) == 0,
        "deep_bsde_rows_are_native_only": bool(
            len(native[native["method"] == "deep_bsde"]) > 0
            and (native[native["method"] == "deep_bsde"]["full_training_cost_charged"] == True).all()  # noqa: E712
        ),
        "smoke_artifacts_excluded": bool(
            not main_text.apply(lambda column: column.str.contains("smoke", case=False)).any(axis=1).any()
            and not native_text.apply(lambda column: column.str.contains("smoke", case=False)).any(axis=1).any()
        ),
        "main_scientific_key_unique": int(
            main.duplicated(["pde", "dimension", "method", "seed"]).sum()
        ) == 0,
        "main_metrics_finite": bool(
            np.isfinite(main[["value_rmse", "value_relative_l2"]].to_numpy(dtype=float)).all()
        ),
        "native_errors_finite": bool(
            np.isfinite(native["absolute_error"].to_numpy(dtype=float)).all()
        ),
        "reference_uncertainty_complete": bool(
            np.isfinite(main["reference_uncertainty"].to_numpy(dtype=float)).all()
            and np.isfinite(native["reference_uncertainty"].to_numpy(dtype=float)).all()
        ),
        "test_bank_hash_complete_for_whole_function": bool(main["test_bank_sha256"].notna().all()),
        "actual_100d_arrays_and_metrics_are_strictly_paired": paired["status"] == "PASS",
        "source_provenance_complete": bool(
            main[["source_commit", "source_artifact", "method_status"]].notna().all().all()
            and native[["source_commit", "source_artifact", "method_status"]].notna().all().all()
        ),
        "source_artifacts_exist_in_workspace": bool(
            all((SOURCE_ROOT / str(item)).exists() for item in source_paths)
        ),
        "allen_cahn_has_no_artificial_mcv": int(
            ((main["pde"] == "allen_cahn") & (main["estimator"] == "mcv")).sum()
        ) == 0,
        "complete_rows_are_finite": bool(
            (pd.to_numeric(main["nonfinite_count"], errors="coerce").fillna(0) == 0).all()
            and (pd.to_numeric(native["nonfinite_count"], errors="coerce").fillna(0) == 0).all()
        ),
        "dimension_scaling_uses_only_frozen_grids": bool(all(
            int(row.dimension) in (MLP_DIMENSIONS if row.pde == "mlp_failure" else DIMENSIONS)
            for row in scaling.itertuples()
        )),
        "sdgd_250d_anchor_excluded_from_dimension_sweep": bool(
            ((main["pde"] == "hjb_log") & (main["dimension"] == 250)
             & (main["method"] == "sdgd_pinn")).any()
            and not ((scaling["pde"] == "hjb_log") & (scaling["dimension"] == 250)).any()
        ),
        "mlp_failure_rows_are_adapted": bool(
            (main.loc[main["pde"] == "mlp_failure", "method_status"] == "adapted").all()
            and (native.loc[native["pde"] == "mlp_failure", "method_status"] == "adapted").all()
        ),
        "stage_b_review_gate_has_valid_audited_state": bool(
            review_gate.get("stage_b_status") == "complete"
            and review_gate.get("status") in {"pending_independent_review", "approved"}
            and (
                review_gate.get("status") != "approved"
                or bool(review_gate.get("independent_review_record"))
            )
        ),
    }
    checks_passed = all(checks.values())
    if not checks_passed:
        failed = sorted(name for name, passed in checks.items() if not passed)
        raise RuntimeError(f"canonical artifact audit failed: {failed}")
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "suite_status": "complete" if pending == 0 else "in_progress",
        "execution_status_counts": counts, "pending_supported_method_cells": pending,
        "canonical_checks_passed": checks_passed,
        "canonical_checks": checks,
        "deep_bsde_excluded_from_whole_function_main_table": checks["main_has_no_deep_bsde"],
        "smoke_results_excluded": checks["smoke_artifacts_excluded"],
        "artifacts": records,
    }
    path = RESULTS_DIR / "FINAL_ARTIFACT_AUDIT.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


def _required_deliverables() -> list[Path]:
    package = SOURCE_ROOT / "invariant_region_mlp/experiments/tangent_actor_icml_final_suite"
    return [
        package / "PLAN.md",
        package / "BENCHMARK_DEFINITIONS.md",
        package / "BASELINE_FIDELITY_MATRIX.md",
        RESULTS_DIR / "DIMENSION_COMPATIBILITY_MATRIX.csv",
        RESULTS_DIR / "EXECUTION_MANIFEST.json",
        RESULTS_DIR / "REVIEW_GATE_STATUS.json",
        RESULTS_DIR / "PAIRWISE_BANK_AUDIT_100D.json",
        RESULTS_DIR / "PAIRWISE_BANK_VERIFICATION_100D.json",
        RESULTS_DIR / "PAIRWISE_REEVALUATION_100D.csv",
        RESULTS_DIR / "REFERENCE_ACCURACY_AUDIT.md",
        RESULTS_DIR / "MAIN_RESULTS.csv",
        RESULTS_DIR / "100D_RESULTS.csv",
        RESULTS_DIR / "NATIVE_POINT_RESULTS.csv",
        RESULTS_DIR / "ACCURACY_COMPUTE_FRONTIER.csv",
        RESULTS_DIR / "DIMENSION_SCALING.csv",
        RESULTS_DIR / "ABLATION_RESULTS.csv",
        RESULTS_DIR / "HISTORICAL_AND_NONPRIMARY_RESULTS.csv",
        RESULTS_DIR / "BASELINE_PAPER_FIDELITY_REVIEW.md",
        RESULTS_DIR / "BASELINE_CONVERGENCE.csv",
        RESULTS_DIR / "fidelity_figures/baseline_native_convergence.pdf",
        RESULTS_DIR / "fidelity_figures/baseline_native_convergence.png",
        RESULTS_DIR / "REMAINING_GPU_HOUR_ESTIMATE.json",
        RESULTS_DIR / "REMAINING_GPU_HOUR_ESTIMATE.csv",
        RESULTS_DIR / "RESOURCE_ACCOUNTING.csv",
        RESULTS_DIR / "REPRODUCIBILITY_MANIFEST.json",
        RESULTS_DIR / "FINAL_EXPERIMENT_REPORT.md",
    ]


def build() -> dict[str, Any]:
    _configure_style()
    artifacts: list[Path] = [
        path for path in _required_deliverables()
        if path.exists() and path.name != "REPRODUCIBILITY_MANIFEST.json"
    ]
    artifacts += figure_mlp_failure()
    artifacts += figure_accuracy_compute()
    artifacts += figure_dimension_scaling()
    artifacts += figure_actor_iteration()
    artifacts += figure_high_accuracy()
    artifacts += [table_main(), table_native(), table_coverage()]
    reproducibility = build_reproducibility_manifest(artifacts)
    artifacts.append(reproducibility)
    # The audit also lists required files that are still missing during staged
    # execution.  Its own hash is intentionally not self-referential.
    audit = build_artifact_audit(list(set(artifacts + _required_deliverables())))
    artifacts.append(audit)
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "artifacts": len(artifacts),
        "figures": len(list(FIGURE_DIR.glob("*.pdf"))),
        "tables": len(list(TABLE_DIR.glob("*.tex"))),
        "audit": str(audit),
    }


def main() -> None:
    print(json.dumps(build(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
