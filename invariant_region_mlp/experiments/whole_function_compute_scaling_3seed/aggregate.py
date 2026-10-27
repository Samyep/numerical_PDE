"""Aggregate three independent training seeds and draw the scaling figure."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np

from invariant_region_mlp.experiments.main_highdim_expansion.benchmarks import (
    array_sha256,
)
from invariant_region_mlp.experiments.whole_function_compute_frontier.common import (
    file_sha256,
    write_json_atomic,
)

from .engine import trajectory_path
from .protocol import (
    CALIBRATION_MANIFEST,
    EXPECTED_INPUTS,
    FIGURE_DIR,
    REPORT_PATH,
    RESULTS_DIR,
    RUN_MANIFEST,
    RUNS_DIR,
    STATUS_DIR,
    TA_INFERENCE_REPLICATES,
    TRAINING_SEEDS,
    formal_specs,
)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _raw_metric_rows() -> list[dict[str, Any]]:
    rows = []
    for path in sorted(RUNS_DIR.glob("*/predictions/**/*.json")):
        row = json.loads(path.read_text(encoding="utf-8"))
        row["metric_file"] = str(path)
        row["metric_file_sha256"] = file_sha256(path)
        rows.append(row)
    return rows


def seed_metrics() -> list[dict[str, Any]]:
    grouped: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for row in _raw_metric_rows():
        key = (
            row["pde"], row["method"], row["sweep"],
            int(row["parameter_value"]), int(row["training_seed"]),
        )
        grouped[key].append(row)
    output: list[dict[str, Any]] = []
    for key, rows in sorted(grouped.items()):
        pde, method, sweep, value, seed = key
        expected = TA_INFERENCE_REPLICATES if method == "tangent_actor" else 1
        if len(rows) != expected:
            raise RuntimeError(f"{key}: expected {expected} evaluation rows, got {len(rows)}")
        rmse = np.asarray([float(row["rmse"]) for row in rows])
        output.append({
            "pde": pde, "method": method, "sweep": sweep,
            "parameter_name": rows[0]["parameter_name"],
            "parameter_value": value, "training_seed": seed,
            "training_seconds": float(rows[0]["training_seconds"]),
            "selection_seconds": float(rows[0].get("selection_seconds", 0.0)),
            "training_and_selection_seconds": float(
                rows[0].get("training_and_selection_seconds", rows[0]["training_seconds"])
            ),
            "rmse": float(rmse.mean()),
            "inference_replicate_rmse_std": (
                float(rmse.std(ddof=1)) if len(rmse) > 1 else 0.0
            ),
            "relative_l2": float(np.mean([float(row["relative_l2"]) for row in rows])),
            "relmae": float(np.mean([float(row["relmae"]) for row in rows])),
            "mae": float(np.mean([float(row["mae"]) for row in rows])),
            "bias": float(np.mean([float(row["bias"]) for row in rows])),
            "nonfinite_count": int(sum(int(row.get("nonfinite_count", 0)) for row in rows)),
            "inference_seconds": float(np.sum([float(row["inference_seconds"]) for row in rows])),
            "inference_replicates": len(rows),
            "training_peak_gpu_bytes": rows[0].get("training_peak_gpu_bytes"),
            "evaluation_peak_gpu_bytes": max(
                int(row.get("evaluation_peak_gpu_bytes") or 0) for row in rows
            ),
            "training_scheduler": json.dumps(
                rows[0].get("training_scheduler", {}), sort_keys=True
            ),
            "evaluation_schedulers": json.dumps(
                [row.get("evaluation_scheduler", {}) for row in rows], sort_keys=True
            ),
            "prediction_files": json.dumps([row["prediction_file"] for row in rows]),
            "prediction_sha256s": json.dumps([row["prediction_sha256"] for row in rows]),
            "checkpoint": rows[0]["checkpoint"],
            "checkpoint_sha256": rows[0]["checkpoint_sha256"],
            "run_id": rows[0]["run_id"],
        })
    return output


def scaling_summary(seed_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for row in seed_rows:
        key = (row["pde"], row["method"], row["sweep"], int(row["parameter_value"]))
        grouped[key].append(row)
    output = []
    for key, rows in sorted(grouped.items()):
        pde, method, sweep, value = key
        seeds = sorted(int(row["training_seed"]) for row in rows)
        if seeds != list(TRAINING_SEEDS):
            raise RuntimeError(f"{key}: expected seeds {TRAINING_SEEDS}, got {seeds}")
        times = np.asarray([float(row["training_seconds"]) for row in rows])
        rmses = np.asarray([float(row["rmse"]) for row in rows])
        inference_times = np.asarray([float(row["inference_seconds"]) for row in rows])
        selection_times = np.asarray([float(row["selection_seconds"]) for row in rows])
        total_method_times = np.asarray([
            float(row["training_and_selection_seconds"]) for row in rows
        ])
        output.append({
            "pde": pde, "method": method, "sweep": sweep,
            "parameter_name": rows[0]["parameter_name"],
            "parameter_value": value, "training_seed_count": len(rows),
            "mean_training_seconds": float(times.mean()),
            "std_training_seconds": float(times.std(ddof=1)),
            "mean_selection_seconds": float(selection_times.mean()),
            "mean_training_and_selection_seconds": float(total_method_times.mean()),
            "mean_inference_seconds": float(inference_times.mean()),
            "std_inference_seconds": float(inference_times.std(ddof=1)),
            "mean_training_peak_gpu_bytes": float(np.mean([
                float(row.get("training_peak_gpu_bytes") or 0) for row in rows
            ])),
            "max_training_peak_gpu_bytes": int(max(
                int(row.get("training_peak_gpu_bytes") or 0) for row in rows
            )),
            "max_evaluation_peak_gpu_bytes": int(max(
                int(row.get("evaluation_peak_gpu_bytes") or 0) for row in rows
            )),
            "mean_rmse": float(rmses.mean()),
            "std_rmse": float(rmses.std(ddof=1)),
            "mean_relative_l2": float(np.mean([float(row["relative_l2"]) for row in rows])),
            "mean_relmae": float(np.mean([float(row["relmae"]) for row in rows])),
            "mean_inference_replicate_rmse_std": float(np.mean([
                float(row["inference_replicate_rmse_std"]) for row in rows
            ])),
            "total_nonfinite_count": int(sum(
                int(row.get("nonfinite_count", 0)) for row in rows
            )),
            "seed_rmses": json.dumps(rmses.tolist()),
            "seed_training_seconds": json.dumps(times.tolist()),
        })
    return output


def _plot(rows: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    colors = {
        "tangent_actor": "#1f77b4", "dpi": "#d55e00", "pinn_hte": "#2ca25f",
    }
    names = {"tangent_actor": "Tangent-Actor", "dpi": "DPI", "pinn_hte": "PINN-HTE"}
    styles = {"C0": "-", "C1": "--", "C2": ":"}
    titles = {"b6_fullrank": "B6 Full-Rank HJB (100D)", "dense_fourier": "Dense Fourier (100D)"}
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), sharey=False)
    for axis, pde in zip(axes, ("b6_fullrank", "dense_fourier")):
        for method in colors:
            for sweep in styles:
                line = sorted(
                    [r for r in rows if r["pde"] == pde and r["method"] == method and r["sweep"] == sweep],
                    key=lambda row: int(row["parameter_value"]),
                )
                if not line:
                    continue
                x = np.asarray([float(r["mean_training_seconds"]) / 60.0 for r in line])
                y = np.asarray([float(r["mean_rmse"]) for r in line])
                yerr = np.asarray([float(r["std_rmse"]) for r in line])
                axis.errorbar(
                    x, y, yerr=yerr, color=colors[method], linestyle=styles[sweep],
                    marker="o", markersize=4.5, linewidth=2.0, capsize=2.5,
                    alpha=0.96,
                )
        axis.set_title(titles[pde])
        axis.set_xlabel("Mean training time per seed (min, log scale)")
        axis.set_ylabel("Whole-function RMSE")
        axis.set_xscale("log")
        axis.set_yscale("log")
        pde_times = [
            float(row["mean_training_seconds"]) / 60.0
            for row in rows if row["pde"] == pde
        ]
        low, high = min(pde_times), max(pde_times)
        candidate_ticks = (3, 5, 10, 20, 30, 50, 100)
        ticks = [
            tick for tick in candidate_ticks
            if low * 0.92 <= tick <= high * 1.08
        ]
        axis.set_xlim(low * 0.88, high * 1.10)
        axis.xaxis.set_major_locator(FixedLocator(ticks))
        axis.xaxis.set_major_formatter(
            FuncFormatter(lambda value, _position: f"{value:g}")
        )
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.grid(True, which="both", alpha=0.22)
    method_handles = [
        Line2D([0], [0], color=colors[m], marker="o", lw=2.5, label=names[m])
        for m in colors
    ]
    sweep_handles = [
        Line2D([0], [0], color="#444444", linestyle=styles[s], lw=2.2,
               label={"C0": "C0: iterations", "C1": "C1: estimator", "C2": "C2: optimization/sample"}[s])
        for s in styles
    ]
    fig.legend(handles=method_handles + sweep_handles, loc="lower center", ncol=6, frameon=False)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(FIGURE_DIR / "two_pde_three_seed_scaling.png", dpi=220)
    fig.savefig(FIGURE_DIR / "two_pde_three_seed_scaling.pdf")
    plt.close(fig)


def _plot_faceted(rows: list[dict[str, Any]]) -> None:
    """A paper-facing small-multiple view of the same frozen results.

    The original two-panel plot is intentionally retained above.  Here the
    sweep identity is encoded by the column rather than by another set of
    line styles, leaving only the three methods to distinguish within each
    panel.  Open markers identify observations beyond the intended 35-minute
    calibration window; those observations remain at their measured times.
    """
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    methods = ("tangent_actor", "dpi", "pinn_hte")
    colors = {
        "tangent_actor": "#0072B2",
        "dpi": "#D55E00",
        "pinn_hte": "#009E73",
    }
    markers = {"tangent_actor": "o", "dpi": "s", "pinn_hte": "^"}
    names = {
        "tangent_actor": "Tangent-Actor",
        "dpi": "DPI",
        "pinn_hte": "PINN-HTE",
    }
    pdes = ("b6_fullrank", "dense_fourier")
    pde_names = {
        "b6_fullrank": "B6 Full-Rank HJB",
        "dense_fourier": "Dense Fourier",
    }
    sweeps = ("C0", "C1", "C2")
    sweep_titles = {
        "C0": "C0  ·  Iterations",
        "C1": "C1  ·  Estimator accuracy",
        "C2": "C2  ·  Optimization / batch",
    }

    fig, axes = plt.subplots(
        2, 3, figsize=(11.7, 6.75), sharey=True,
        gridspec_kw={"hspace": 0.24, "wspace": 0.12},
    )
    for row_index, pde in enumerate(pdes):
        for column_index, sweep in enumerate(sweeps):
            axis = axes[row_index, column_index]
            panel_rows = [
                row for row in rows
                if row["pde"] == pde and row["sweep"] == sweep
            ]
            for method in methods:
                line = sorted(
                    [row for row in panel_rows if row["method"] == method],
                    key=lambda row: int(row["parameter_value"]),
                )
                if not line:
                    continue
                x = np.asarray([
                    float(row["mean_training_seconds"]) / 60.0 for row in line
                ])
                y = np.asarray([float(row["mean_rmse"]) for row in line])
                yerr = np.asarray([float(row["std_rmse"]) for row in line])
                axis.errorbar(
                    x,
                    y,
                    yerr=yerr,
                    color=colors[method],
                    linestyle="-",
                    marker=markers[method],
                    markerfacecolor=colors[method],
                    markeredgecolor="white",
                    markeredgewidth=0.65,
                    markersize=5.2,
                    linewidth=1.8,
                    elinewidth=0.9,
                    capsize=2.0,
                    alpha=0.93,
                    zorder=3,
                )
                over_budget = x > 35.0
                if np.any(over_budget):
                    axis.scatter(
                        x[over_budget],
                        y[over_budget],
                        marker=markers[method],
                        s=38,
                        facecolor="white",
                        edgecolor=colors[method],
                        linewidth=1.25,
                        zorder=5,
                    )

            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set_ylim(4.2e-4, 3.3e-1)
            panel_times = [
                float(row["mean_training_seconds"]) / 60.0
                for row in panel_rows
            ]
            low, high = min(panel_times), max(panel_times)
            axis.set_xlim(max(2.6, low * 0.82), high * 1.14)
            candidate_ticks = (3, 5, 10, 20, 30, 50, 100)
            ticks = [
                tick for tick in candidate_ticks
                if axis.get_xlim()[0] <= tick <= axis.get_xlim()[1]
            ]
            axis.xaxis.set_major_locator(FixedLocator(ticks))
            axis.xaxis.set_major_formatter(
                FuncFormatter(lambda value, _position: f"{value:g}")
            )
            axis.xaxis.set_minor_formatter(NullFormatter())
            axis.set_axisbelow(True)
            axis.grid(True, which="major", color="#DDDDDD", linewidth=0.70)
            axis.grid(False, which="minor")
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
            axis.tick_params(axis="both", labelsize=8.8)
            if row_index == 0:
                axis.set_title(sweep_titles[sweep], fontsize=10.5, pad=9)
            if column_index == 0:
                axis.set_ylabel(
                    f"{pde_names[pde]}\nWhole-function RMSE",
                    fontsize=9.5,
                    labelpad=8,
                )
            if any(time > 35.0 for time in panel_times):
                axis.axvline(
                    35.0, color="#8A8A8A", linestyle=(0, (3, 3)),
                    linewidth=0.9, alpha=0.8, zorder=1,
                )

    method_handles = [
        Line2D(
            [0], [0], color=colors[method], marker=markers[method],
            markerfacecolor=colors[method], markeredgecolor="white",
            markeredgewidth=0.65, linewidth=2.0, markersize=6.2,
            label=names[method],
        )
        for method in methods
    ]
    fig.legend(
        handles=method_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.005),
        ncol=3,
        frameon=False,
        columnspacing=2.2,
        handlelength=2.4,
        fontsize=10.0,
    )
    fig.text(
        0.995,
        0.012,
        "Mean ± s.d. across three training seeds; dashed line: 35-min boundary; open markers: beyond boundary",
        ha="right",
        va="bottom",
        fontsize=8.0,
        color="#666666",
    )
    fig.supxlabel("Mean training time per seed (min, log scale)", y=0.047, fontsize=10.0)
    fig.subplots_adjust(left=0.095, right=0.985, bottom=0.125, top=0.90)
    fig.savefig(
        FIGURE_DIR / "two_pde_three_seed_scaling_faceted.png",
        dpi=240,
        bbox_inches="tight",
    )
    fig.savefig(
        FIGURE_DIR / "two_pde_three_seed_scaling_faceted.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)


def _audit(seed_rows: list[dict[str, Any]], summary: list[dict[str, Any]]) -> dict[str, Any]:
    specs = formal_specs()
    expected_runs = {spec.run_id for spec in specs}
    trajectory_failures = []
    for spec in specs:
        path = trajectory_path(spec)
        if not path.is_file():
            trajectory_failures.append({"run_id": spec.run_id, "failure": "missing trajectory"})
            continue
        record = json.loads(path.read_text(encoding="utf-8"))
        if record.get("status") != "COMPLETE":
            trajectory_failures.append({"run_id": spec.run_id, "failure": record.get("status")})
        if record.get("ground_truth_used_for_training_selection_or_stopping") is not False:
            trajectory_failures.append({"run_id": spec.run_id, "failure": "truth-policy flag"})
        randomness = record.get("randomness", {})
        if int(randomness.get("training_seed", -1)) != int(spec.seed):
            trajectory_failures.append({"run_id": spec.run_id, "failure": "training-seed metadata"})
        if len(str(randomness.get("initialization_key_sha256", ""))) != 64:
            trajectory_failures.append({"run_id": spec.run_id, "failure": "initialization fingerprint"})
        if spec.method in {"dpi", "pinn_hte"}:
            expanded = record["configuration"]["expanded_config"]
            if "FREQ: null" not in expanded:
                trajectory_failures.append({"run_id": spec.run_id, "failure": "exact callback not disabled"})
    artifact_failures = []
    for row in seed_rows:
        checkpoint = Path(row["checkpoint"])
        if not checkpoint.is_file() or file_sha256(checkpoint) != row["checkpoint_sha256"]:
            artifact_failures.append({"run_id": row["run_id"], "failure": "checkpoint hash"})
        for name, digest in zip(json.loads(row["prediction_files"]), json.loads(row["prediction_sha256s"])):
            path = Path(name)
            if not path.is_file() or file_sha256(path) != digest:
                artifact_failures.append({"run_id": row["run_id"], "failure": "prediction hash"})
                continue
            with np.load(path, allow_pickle=False) as payload:
                if array_sha256(np.asarray(payload["t"]), np.asarray(payload["x"])) != EXPECTED_INPUTS[row["pde"]]["test_input_sha256"]:
                    artifact_failures.append({"run_id": row["run_id"], "failure": "test-input hash"})
                if array_sha256(np.asarray(payload["truth"])) != EXPECTED_INPUTS[row["pde"]]["truth_value_sha256"]:
                    artifact_failures.append({"run_id": row["run_id"], "failure": "truth hash"})
    observed_runs = {row["run_id"] for row in seed_rows}
    if observed_runs != expected_runs:
        artifact_failures.append({
            "failure": "formal run-set mismatch",
            "missing": sorted(expected_runs - observed_runs),
            "unexpected": sorted(observed_runs - expected_runs),
        })
    expected_points = {
        (spec.pde, spec.method, spec.sweep, int(value))
        for spec in specs
        for value in (
            spec.checkpoint_values if spec.sweep == "C0"
            else (spec.parameter_value,)
        )
    }
    observed_points = {
        (row["pde"], row["method"], row["sweep"], int(row["parameter_value"]))
        for row in summary
    }
    if observed_points != expected_points:
        artifact_failures.append({
            "failure": "summary point-set mismatch",
            "missing": sorted(expected_points - observed_points),
            "unexpected": sorted(observed_points - expected_points),
        })

    provenance_failures = []
    gate_path = STATUS_DIR / "formal_training_complete.json"
    if not RUN_MANIFEST.is_file():
        provenance_failures.append({"failure": "missing run_manifest.json"})
    if not gate_path.is_file():
        provenance_failures.append({"failure": "missing global truth-gate record"})
    if RUN_MANIFEST.is_file() and gate_path.is_file():
        gate = json.loads(gate_path.read_text(encoding="utf-8"))
        if gate.get("run_manifest_sha256") != file_sha256(RUN_MANIFEST):
            provenance_failures.append({"failure": "run-manifest gate hash mismatch"})
    resource_path = RESULTS_DIR / "resource_accounting.json"
    resource = (
        json.loads(resource_path.read_text(encoding="utf-8"))
        if resource_path.is_file() else None
    )
    if resource is None:
        provenance_failures.append({"failure": "missing resource accounting"})
    elif resource.get("within_90_h100_gpu_hour_cap") is not True:
        provenance_failures.append({"failure": "90 H100 GPU-hour cap exceeded"})
    chunk_validation_path = RESULTS_DIR / "dpi_chunking_validation.json"
    if not chunk_validation_path.is_file():
        provenance_failures.append({"failure": "missing DPI chunking validation"})
    else:
        chunk_validation = json.loads(
            chunk_validation_path.read_text(encoding="utf-8")
        )
        if chunk_validation.get("status") != "PASS":
            provenance_failures.append({"failure": "DPI chunking validation failed"})

    numerical_failures = [
        {
            "pde": row["pde"], "method": row["method"],
            "sweep": row["sweep"], "parameter_value": row["parameter_value"],
            "total_nonfinite_count": row["total_nonfinite_count"],
        }
        for row in summary if int(row.get("total_nonfinite_count", 0)) > 0
    ]
    payload = {
        "overall_pass": (
            not trajectory_failures
            and not artifact_failures
            and not provenance_failures
            and not numerical_failures
        ),
        "expected_formal_runs": len(expected_runs),
        "observed_formal_runs": len(observed_runs),
        "expected_summary_points": len(expected_points),
        "observed_summary_points": len(observed_points),
        "trajectory_failures": trajectory_failures,
        "artifact_failures": artifact_failures,
        "provenance_failures": provenance_failures,
        "numerical_failures": numerical_failures,
        "summary_points": len(summary),
        "aggregation_is_mean_of_per_seed_rmse_not_prediction_ensemble": True,
        "running_minimum_or_test_selected_envelope_used": False,
        "calibration_manifest_sha256": file_sha256(CALIBRATION_MANIFEST),
    }
    write_json_atomic(RESULTS_DIR / "artifact_audit.json", payload)
    return payload


def _report(summary: list[dict[str, Any]], audit: dict[str, Any]) -> None:
    calibration = json.loads(CALIBRATION_MANIFEST.read_text(encoding="utf-8"))
    resource_path = RESULTS_DIR / "resource_accounting.json"
    resource = (
        json.loads(resource_path.read_text(encoding="utf-8"))
        if resource_path.is_file() else None
    )
    run_manifest = (
        json.loads(RUN_MANIFEST.read_text(encoding="utf-8"))
        if RUN_MANIFEST.is_file() else {}
    )
    chunk_validation_path = RESULTS_DIR / "dpi_chunking_validation.json"
    chunk_validation = (
        json.loads(chunk_validation_path.read_text(encoding="utf-8"))
        if chunk_validation_path.is_file() else {}
    )
    lines = [
        "# Three-seed hyperparameter-specific accuracy--compute scaling",
        "",
        "## Completeness and audit",
        "",
        "| Check | Result |",
        "|---|---:|",
        f"| Artifact audit | **{'PASS' if audit['overall_pass'] else 'FAIL'}** |",
        f"| Formal runs | {audit['observed_formal_runs']} / {audit['expected_formal_runs']} |",
        f"| Aggregated curve points | {audit['observed_summary_points']} / {audit['expected_summary_points']} |",
        "| Training seeds per point | 3 (0, 1, 2) |",
        f"| Frozen calibration grid | {calibration.get('status', 'MISSING')} |",
        f"| Resource cap | {('PASS' if resource and resource.get('within_90_h100_gpu_hour_cap') else 'FAIL / NOT COLLECTED')} |",
        "",
        "Each point is the mean of three separately computed per-seed RMSE values; error bars are their sample standard deviation. Predictions are not ensembled.",
        "",
        "## Descriptive best observed points",
        "",
        "These post-evaluation minima summarize the measured grid only; they were not used for training, stopping, checkpoint selection, or a test-selected envelope.",
        "",
        "| PDE | Method | Sweep | Parameter | Mean time (min) | Mean RMSE | Ratio to TA best |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for pde in ("b6_fullrank", "dense_fourier"):
        pde_rows = [row for row in summary if row["pde"] == pde]
        ta_best = min(
            float(row["mean_rmse"])
            for row in pde_rows if row["method"] == "tangent_actor"
        )
        for method in ("tangent_actor", "dpi", "pinn_hte"):
            row = min(
                (item for item in pde_rows if item["method"] == method),
                key=lambda item: float(item["mean_rmse"]),
            )
            lines.append(
                f"| {pde} | {method} | {row['sweep']} | {row['parameter_value']} | "
                f"{float(row['mean_training_seconds']) / 60:.3f} | "
                f"{float(row['mean_rmse']):.8g} | "
                f"{float(row['mean_rmse']) / ta_best:.2f}x |"
            )
    lines += [
        "",
        "## Main results",
        "",
        "| PDE | Method | Sweep | Parameter | Mean time (min) | Mean RMSE | Seed SD |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in summary:
        lines.append(
            f"| {row['pde']} | {row['method']} | {row['sweep']} | {row['parameter_value']} | "
            f"{float(row['mean_training_seconds']) / 60:.3f} | {float(row['mean_rmse']):.8g} | {float(row['std_rmse']):.3g} |"
        )
    lines += [
        "", "## Interpretation boundary", "",
        "These are one-parameter scaling curves, not globally optimized Pareto frontiers. C0 points are checkpoints from one trajectory per seed; C1/C2 points are independent fresh runs. No running minimum, smoothing, interpolation, or test-selected envelope is used.",
        "The plotted x-axis is synchronized active training time, as frozen. Reference-free TA selector time, inference time, scheduler wall time, and allocated GPU time are retained separately rather than silently folded into that axis.",
        "", "## Timing calibration and frozen grids", "",
        f"Calibration used no accuracy and projected {float(calibration['projected_formal_gpu_hours']):.2f} formal H100 GPU-hours before formal training.",
        "",
        "| PDE | Method | Sweep | Frozen parameter values |",
        "|---|---|---|---|",
    ]
    for row in calibration.get("selected_lines", []):
        lines.append(
            f"| {row['pde']} | {row['method']} | {row['sweep']} | "
            f"{', '.join(str(v) for v in row['parameter_values'])} |"
        )
    lines += [
        "", "## Resource accounting and failures", "",
    ]
    if resource:
        lines += [
            f"Total allocated H100 time (including calibration failures and retries): **{float(resource['allocated_gpu_hours']):.3f} GPU-hours**.",
            "",
            f"Scheduler rows not marked COMPLETED: {len(resource.get('failed_or_cancelled_rows', []))}.",
        ]
    else:
        lines.append("Resource accounting was not collected; this is an audit failure.")
    failures = (
        audit.get("trajectory_failures", [])
        + audit.get("artifact_failures", [])
        + audit.get("provenance_failures", [])
        + audit.get("numerical_failures", [])
    )
    lines += ["", "Audit failures and recorded deviations:", ""]
    if failures:
        lines.extend(f"- `{json.dumps(item, sort_keys=True)}`" for item in failures)
    else:
        lines.append("- Final artifact audit: no failures.")
    if resource:
        noncompleted = resource.get("failed_or_cancelled_rows", [])
        state_counts = Counter(
            str(row.get("state", "UNKNOWN")).split()[0]
            for row in noncompleted
        )
        lines.append(
            "- Scheduler history retained and counted: "
            + ", ".join(
                f"{state}={count}" for state, count in sorted(state_counts.items())
            )
            + "."
        )
        if state_counts.get("FAILED", 0) or state_counts.get("CANCELLED", 0):
            lines.append(
                "- The FAILED and CANCELLED scheduler rows are calibration/"
                "preflight attempts (mostly zero- or few-second startup "
                "failures on subsequently excluded nodes, plus curtailed "
                "calibration arrays); they do not represent missing formal "
                "runs."
            )
        formal_timeouts = [
            row for row in noncompleted
            if row.get("state") == "TIMEOUT"
        ]
        if formal_timeouts:
            lines.append(
                "- Formal deviation: the three B6 PINN-HTE C2 runs at "
                "batch size 16384 timed out at 70:12, then completed with "
                "the identical numerical configuration after only the Slurm "
                "wall-time limit was raised; both attempts are included in "
                "resource accounting."
            )
    over_window = [
        row for row in summary
        if float(row["mean_training_seconds"]) / 60.0 > 35.0
    ]
    if over_window:
        details = ", ".join(
            f"{row['pde']}/{row['method']}/{row['sweep']}={float(row['mean_training_seconds']) / 60:.1f} min"
            for row in over_window
        )
        lines.append(
            "- Timing-calibration deviation: points beyond the intended "
            f"25--35 minute endpoint are retained ({details})."
        )
    pre_fix_audit = RESULTS_DIR / "artifact_audit_pre_dpi_c0_evaluation_fix.json"
    if pre_fix_audit.is_file():
        lines.append(
            "- A pre-fix artifact audit failed because the evaluator skipped "
            "pre-registered DPI C0 checkpoints at 80/160/320. The FAIL record "
            "is preserved in `artifact_audit_pre_dpi_c0_evaluation_fix.json`; "
            "the evaluator-only fix did not retrain models, and the six C0 "
            "trajectories were re-evaluated in Slurm job 6619464."
        )
    archived_chunk_dir = CALIBRATION_MANIFEST.parent / "calibration" / "archived_pre_final_chunk"
    if archived_chunk_dir.is_dir():
        archived_count = sum(path.is_dir() for path in archived_chunk_dir.iterdir())
        lines.append(
            f"- {archived_count} pre-final DPI chunk-timing directories were "
            "moved to `calibration/archived_pre_final_chunk/` and preserved; "
            "they were not used to freeze the final grid."
        )
    lines += [
        "", "## DPI large-sample chunking validation", "",
        f"Status: **{chunk_validation.get('status', 'MISSING')}**; all nominal samples used: "
        f"`{chunk_validation.get('all_nominal_samples_used', False)}`; one-chunk maximum absolute difference "
        f"`{chunk_validation.get('single_chunk_max_abs', 'NA')}`; true multi-chunk maximum absolute difference "
        f"`{chunk_validation.get('multi_chunk_max_abs', 'NA')}`.",
        "", "## Dense Fourier DPI provenance", "",
        "The preceding study's historical C0 consistency check remains FAIL. Its RNG-only diagnostic explained 96.72% of the discrepancy but was not bitwise. This study disables the historical exact-solution callback (`EVAL.FREQ=null`) and averages new truth-gated runs over seeds 0, 1, and 2.",
        "", "## Provenance", "",
        f"Calibration manifest SHA-256: `{audit['calibration_manifest_sha256']}`.",
        f"Run manifest SHA-256: `{file_sha256(RUN_MANIFEST) if RUN_MANIFEST.is_file() else 'MISSING'}`.",
        f"Checked-out training revision: `{run_manifest.get('checked_out_revision', 'MISSING')}`.",
        f"Frozen input hashes: `{json.dumps(calibration.get('frozen_input_hashes', {}), sort_keys=True)}`.",
        f"Seed metrics SHA-256: `{file_sha256(RESULTS_DIR / 'seed_metrics.csv')}`.",
        f"Scaling summary SHA-256: `{file_sha256(RESULTS_DIR / 'scaling_summary.csv')}`.",
        f"Figure PDF SHA-256: `{file_sha256(FIGURE_DIR / 'two_pde_three_seed_scaling.pdf')}`.",
        f"Faceted figure PDF SHA-256: `{file_sha256(FIGURE_DIR / 'two_pde_three_seed_scaling_faceted.pdf')}`.",
        f"DPI chunking validation SHA-256: `{file_sha256(RESULTS_DIR / 'dpi_chunking_validation.json') if (RESULTS_DIR / 'dpi_chunking_validation.json').is_file() else 'MISSING'}`.",
    ]
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def aggregate() -> dict[str, Any]:
    seed_rows = seed_metrics()
    summary = scaling_summary(seed_rows)
    _write_csv(RESULTS_DIR / "seed_metrics.csv", seed_rows)
    _write_csv(RESULTS_DIR / "scaling_summary.csv", summary)
    _plot(summary)
    _plot_faceted(summary)
    audit = _audit(seed_rows, summary)
    _report(summary, audit)
    payload = {
        "seed_metric_rows": len(seed_rows),
        "summary_points": len(summary),
        "audit_pass": audit["overall_pass"],
        "figure": str(FIGURE_DIR / "two_pde_three_seed_scaling.pdf"),
        "faceted_figure": str(FIGURE_DIR / "two_pde_three_seed_scaling_faceted.pdf"),
        "report": str(REPORT_PATH),
    }
    write_json_atomic(RESULTS_DIR / "analysis_summary.json", payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=2, sort_keys=True))
