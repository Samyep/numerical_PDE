"""Aggregate three independent training seeds and draw the scaling figure."""

from __future__ import annotations

from collections import defaultdict
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
    RUNS_DIR,
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
            "rmse": float(rmse.mean()),
            "inference_replicate_rmse_std": (
                float(rmse.std(ddof=1)) if len(rmse) > 1 else 0.0
            ),
            "relative_l2": float(np.mean([float(row["relative_l2"]) for row in rows])),
            "relmae": float(np.mean([float(row["relmae"]) for row in rows])),
            "mae": float(np.mean([float(row["mae"]) for row in rows])),
            "bias": float(np.mean([float(row["bias"]) for row in rows])),
            "inference_seconds": float(np.sum([float(row["inference_seconds"]) for row in rows])),
            "inference_replicates": len(rows),
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
        output.append({
            "pde": pde, "method": method, "sweep": sweep,
            "parameter_name": rows[0]["parameter_name"],
            "parameter_value": value, "training_seed_count": len(rows),
            "mean_training_seconds": float(times.mean()),
            "std_training_seconds": float(times.std(ddof=1)),
            "mean_rmse": float(rmses.mean()),
            "std_rmse": float(rmses.std(ddof=1)),
            "mean_relative_l2": float(np.mean([float(row["relative_l2"]) for row in rows])),
            "mean_relmae": float(np.mean([float(row["relmae"]) for row in rows])),
            "mean_inference_replicate_rmse_std": float(np.mean([
                float(row["inference_replicate_rmse_std"]) for row in rows
            ])),
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
        axis.set_xlabel("Mean training time per seed (min)")
        axis.set_ylabel("Whole-function RMSE")
        axis.set_yscale("log")
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
    payload = {
        "overall_pass": not trajectory_failures and not artifact_failures,
        "expected_formal_runs": len(expected_runs),
        "observed_formal_runs": len({row["run_id"] for row in seed_rows}),
        "trajectory_failures": trajectory_failures,
        "artifact_failures": artifact_failures,
        "summary_points": len(summary),
        "aggregation_is_mean_of_per_seed_rmse_not_prediction_ensemble": True,
        "running_minimum_or_test_selected_envelope_used": False,
        "calibration_manifest_sha256": file_sha256(CALIBRATION_MANIFEST),
    }
    write_json_atomic(RESULTS_DIR / "artifact_audit.json", payload)
    return payload


def _report(summary: list[dict[str, Any]], audit: dict[str, Any]) -> None:
    calibration = json.loads(CALIBRATION_MANIFEST.read_text(encoding="utf-8"))
    lines = [
        "# Three-seed hyperparameter-specific accuracy--compute scaling",
        "",
        f"Artifact audit: **{'PASS' if audit['overall_pass'] else 'FAIL'}**. ",
        "Each point is the mean of three separately computed per-seed RMSE values; error bars are their sample standard deviation. Predictions are not ensembled.",
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
        "", "## Timing calibration", "",
        f"Calibration used no accuracy and projected {float(calibration['projected_formal_gpu_hours']):.2f} formal H100 GPU-hours before formal training.",
        "", "## Dense Fourier DPI provenance", "",
        "The preceding study's historical C0 consistency check remains FAIL. Its RNG-only diagnostic explained 96.72% of the discrepancy but was not bitwise. This study disables the historical exact-solution callback (`EVAL.FREQ=null`) and averages new truth-gated runs over seeds 0, 1, and 2.",
        "", "## Provenance", "",
        f"Calibration manifest SHA-256: `{audit['calibration_manifest_sha256']}`.",
    ]
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def aggregate() -> dict[str, Any]:
    seed_rows = seed_metrics()
    summary = scaling_summary(seed_rows)
    _write_csv(RESULTS_DIR / "seed_metrics.csv", seed_rows)
    _write_csv(RESULTS_DIR / "scaling_summary.csv", summary)
    _plot(summary)
    audit = _audit(seed_rows, summary)
    _report(summary, audit)
    payload = {
        "seed_metric_rows": len(seed_rows),
        "summary_points": len(summary),
        "audit_pass": audit["overall_pass"],
        "figure": str(FIGURE_DIR / "two_pde_three_seed_scaling.pdf"),
        "report": str(REPORT_PATH),
    }
    write_json_atomic(RESULTS_DIR / "analysis_summary.json", payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=2, sort_keys=True))
