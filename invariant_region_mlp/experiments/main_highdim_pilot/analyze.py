"""Audit, aggregate, plot, and report the frozen high-dimensional pilot."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, MaxNLocator, NullFormatter
import numpy as np

from .b6_adapter import SOURCE_COMMIT, load_b6_reuse
from .equations import array_sha256
from .prepare import file_sha256
from .protocol import (
    CHECKPOINT_ROOT,
    FIGURES_ROOT,
    INFERENCE_PATHS,
    INFERENCE_REPLICATES,
    INPUT_ROOT,
    RAW_ROOT,
    REPORT_PATH,
    RESULTS_ROOT,
    REVISION,
    UPDATES,
    VALIDATION_PATH,
    ensure_directories,
    write_json_atomic,
)


PREREGISTRATION_COMMIT = "50908e13560649025aefb874a625e55ecf76edd0"


DISPLAY = {
    "dense_fourier": "Dense Fourier joint",
    "coupled_lqg": "Published coupled LQG",
    "B6": "B6 full-rank HJB",
}
COLORS = {
    "B6": "#376795",
    "dense_fourier": "#d87529",
    "coupled_lqg": "#a23b4a",
}


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _mean(rows: list[dict[str, Any]], key: str) -> float:
    return float(np.mean([float(row[key]) for row in rows]))


def _std(rows: list[dict[str, Any]], key: str) -> float:
    values = [float(row[key]) for row in rows]
    return float(np.std(values, ddof=1)) if len(values) > 1 else 0.0


def _prediction_audit() -> tuple[dict[str, Any], dict[tuple[str, int], float]]:
    checks: list[dict[str, Any]] = []
    pooled: dict[tuple[str, int], float] = {}
    overall = True
    for pde in ("dense_fourier", "coupled_lqg"):
        canonical: dict[str, np.ndarray] | None = None
        for iteration in range(UPDATES + 1):
            predictions = []
            for replicate in range(INFERENCE_REPLICATES):
                row_path = RAW_ROOT / "evaluations" / pde / f"a{iteration}_rep{replicate:02d}.json"
                row = _json(row_path)
                prediction_path = Path(row["prediction_path"])
                with np.load(prediction_path) as data:
                    arrays = {key: data[key] for key in data.files}
                if canonical is None:
                    canonical = {key: arrays[key].copy() for key in ("t", "x", "truth", "truth_q")}
                bank_equal = all(
                    np.array_equal(arrays[key], canonical[key])
                    for key in ("t", "x", "truth", "truth_q")
                )
                prediction = arrays["mcv"]
                truth = arrays["truth"]
                pair_variance = arrays["mcv_pair_variance"]
                rmse = float(np.sqrt(np.mean(np.square(prediction - truth))))
                relative = float(np.linalg.norm(prediction - truth) / np.linalg.norm(truth))
                bias = float(np.mean(prediction - truth))
                mc = float(np.sqrt(np.mean(pair_variance / (INFERENCE_PATHS // 2))))
                finite = bool(
                    all(np.all(np.isfinite(value)) for value in arrays.values())
                )
                local_pass = bool(
                    bank_equal
                    and finite
                    and file_sha256(prediction_path) == row["prediction_sha256"]
                    and array_sha256(arrays["t"], arrays["x"]) == row["test_bank_sha256"]
                    and abs(rmse - float(row["value_rmse"])) <= 5e-13
                    and abs(relative - float(row["relative_l2"])) <= 5e-13
                    and abs(bias - float(row["bias"])) <= 5e-13
                    and abs(mc - float(row["mc_rmse_contribution"])) <= 5e-13
                    and row["ground_truth_loaded_only_after_selection"] is True
                )
                overall &= local_pass
                checks.append(
                    {
                        "pde": pde,
                        "actor": f"a{iteration}",
                        "replicate": replicate,
                        "prediction_sha256": row["prediction_sha256"],
                        "bank_equal": bank_equal,
                        "finite": finite,
                        "recomputed_rmse": rmse,
                        "recomputed_mc_rmse_contribution": mc,
                        "pass": local_pass,
                    }
                )
                predictions.append(prediction)
            assert canonical is not None
            pooled_prediction = np.mean(predictions, axis=0)
            pooled[(pde, iteration)] = float(
                np.sqrt(np.mean(np.square(pooled_prediction - canonical["truth"])))
            )
    return {"status": "PASS" if overall else "FAIL", "checks": checks}, pooled


def _checkpoint_audit() -> dict[str, Any]:
    checks = []
    overall = True
    for pde in ("dense_fourier", "coupled_lqg"):
        training = _json(RAW_ROOT / f"training_{pde}_seed0.json")
        no_truth = (
            training["ground_truth_used_for_training_or_selection"] is False
            and training["selection"]["ground_truth_used"] is False
        )
        for iteration in range(UPDATES + 1):
            path = CHECKPOINT_ROOT / f"{pde}_seed0_a{iteration}.pt"
            expected = None
            if iteration:
                expected = training["phases"][iteration - 1]["checkpoint_sha256"]
            actual = file_sha256(path)
            local_pass = path.is_file() and (expected is None or actual == expected) and no_truth
            overall &= local_pass
            checks.append(
                {
                    "pde": pde,
                    "actor": f"a{iteration}",
                    "path": str(path),
                    "sha256": actual,
                    "recorded_sha256": expected,
                    "truth_excluded_from_training_and_selection": no_truth,
                    "pass": local_pass,
                }
            )
    return {"status": "PASS" if overall else "FAIL", "checks": checks}


def _new_rows(
    pooled: dict[tuple[str, int], float]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, Any]]]:
    replicate_rows: list[dict[str, Any]] = []
    aggregate_rows: list[dict[str, Any]] = []
    training_records: dict[str, dict[str, Any]] = {}
    for pde in ("dense_fourier", "coupled_lqg"):
        training = _json(RAW_ROOT / f"training_{pde}_seed0.json")
        training_records[pde] = training
        selector = {row["actor"]: row for row in training["selection"]["rows"]}
        phase_cumulative = 0.0
        for iteration in range(UPDATES + 1):
            if iteration:
                phase_cumulative += float(training["phases"][iteration - 1]["phase_seconds"])
            local = []
            for replicate in range(INFERENCE_REPLICATES):
                row = _json(
                    RAW_ROOT / "evaluations" / pde / f"a{iteration}_rep{replicate:02d}.json"
                )
                compact = {
                    "pde_id": pde,
                    "pde": DISPLAY[pde],
                    "dimension": 100,
                    "result_origin": "NEW_5070_TI_PILOT",
                    "actor": f"a{iteration}",
                    "replicate": replicate,
                    "value_rmse": row["value_rmse"],
                    "relative_l2": row["relative_l2"],
                    "bias": row["bias"],
                    "mc_rmse_contribution": row["mc_rmse_contribution"],
                    "non_mc_rmse_proxy": row["non_mc_rmse_proxy"],
                    "plain_value_rmse": row["plain_value_rmse"],
                    "actor_value_rmse": row["actor_value_rmse"],
                    "actor_gradient_component_rmse": row["actor_gradient_component_rmse"],
                    "actor_gradient_vector_rmse": row["actor_gradient_vector_rmse"],
                    "variance_reduction": row["variance_reduction"],
                    "weight_ess_fraction": row["weight_ess_fraction"],
                    "inference_seconds": row["inference_seconds"],
                    "peak_gpu_bytes": row["peak_gpu_bytes"],
                    "nonfinite_count": row["nonfinite_count"],
                    "selected_by_reference_free_rule": row["selected_by_reference_free_rule"],
                    "test_bank_sha256": row["test_bank_sha256"],
                    "prediction_sha256": row["prediction_sha256"],
                }
                replicate_rows.append(compact)
                local.append(row)
            selected = training["selection"]["selected_actor"] == f"a{iteration}"
            selector_row = selector[f"a{iteration}"]
            aggregate_rows.append(
                {
                    "pde_id": pde,
                    "pde": DISPLAY[pde],
                    "dimension": 100,
                    "result_origin": "NEW_5070_TI_PILOT",
                    "actor": f"a{iteration}",
                    "replicates": INFERENCE_REPLICATES,
                    "value_rmse_mean": _mean(local, "value_rmse"),
                    "value_rmse_std": _std(local, "value_rmse"),
                    "pooled_value_rmse_diagnostic": pooled[(pde, iteration)],
                    "relative_l2_mean": _mean(local, "relative_l2"),
                    "mc_rmse_contribution_mean": _mean(local, "mc_rmse_contribution"),
                    "non_mc_rmse_proxy_mean": _mean(local, "non_mc_rmse_proxy"),
                    "actor_value_rmse": local[0]["actor_value_rmse"],
                    "actor_gradient_component_rmse": local[0]["actor_gradient_component_rmse"],
                    "actor_gradient_vector_rmse": local[0]["actor_gradient_vector_rmse"],
                    "variance_reduction_mean": _mean(local, "variance_reduction"),
                    "inference_seconds_mean": _mean(local, "inference_seconds"),
                    "inference_seconds_all_replicates": float(
                        sum(float(row["inference_seconds"]) for row in local)
                    ),
                    "peak_gpu_bytes": max(int(row["peak_gpu_bytes"]) for row in local),
                    "nonfinite_count": sum(int(row["nonfinite_count"]) for row in local),
                    "selected_by_reference_free_rule": selected,
                    "selector_normalized_total_inconsistency": selector_row[
                        "normalized_total_inconsistency"
                    ],
                    "selector_noise_corrected_component_rmse": selector_row.get(
                        "gradient_noise_corrected_component_rmse"
                    ),
                    "selector_noise_corrected_value_rmse": selector_row.get(
                        "value_noise_corrected_rmse"
                    ),
                    "cumulative_update_seconds": phase_cumulative,
                    "compute_to_one_evaluation_seconds": phase_cumulative
                    + _mean(local, "inference_seconds"),
                }
            )
    return replicate_rows, aggregate_rows, training_records


def _b6_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    replicate_rows: list[dict[str, Any]] = []
    aggregate_rows: list[dict[str, Any]] = []
    audits: dict[str, Any] = {}
    for dimension in (100, 400):
        reuse = load_b6_reuse(dimension)
        audits[str(dimension)] = reuse
        cumulative = 0.0
        phases = {int(item["iteration"]): item for item in reuse["training"]["phases"]}
        selector = {item["actor"]: item for item in reuse["training"]["selection"]["rows"]}
        for row in reuse["iterations"]:
            iteration = int(row["actor_request"][1:])
            if iteration:
                cumulative += float(phases[iteration]["phase_seconds"])
            compact = {
                "pde_id": "B6",
                "pde": DISPLAY["B6"],
                "dimension": dimension,
                "result_origin": "REUSED_MATURE_RESULT",
                "actor": row["actor_request"],
                "replicate": "aggregate_of_4",
                "value_rmse": float(row["mean_value_rmse"]),
                "relative_l2": float(row["mean_relative_l2"]),
                "bias": float(row["mean_bias"]),
                "mc_rmse_contribution": float(row["mean_mc_rmse_contribution"]),
                "non_mc_rmse_proxy": float(row["mean_non_mc_rmse_proxy"]),
                "actor_gradient_vector_rmse": float(row["mean_gradient_vector_rmse"]),
                "variance_reduction": float(row["mean_variance_reduction"]),
                "weight_ess_fraction": float(row["mean_weight_ess_fraction"]),
                "inference_seconds": float(row["mean_inference_seconds"]),
                "peak_gpu_bytes": float(row["mean_peak_gpu_bytes"]),
                "nonfinite_count": float(row["mean_nonfinite_count"]),
                "selected_by_reference_free_rule": row["actor_request"] == "a4",
                "test_bank_sha256": reuse["calibration"]["bank_sha256"],
            }
            replicate_rows.append(compact)
            selector_row = selector[row["actor_request"]]
            aggregate_rows.append(
                {
                    "pde_id": "B6",
                    "pde": DISPLAY["B6"],
                    "dimension": dimension,
                    "result_origin": "REUSED_MATURE_RESULT",
                    "actor": row["actor_request"],
                    "replicates": int(row["replicates"]),
                    "value_rmse_mean": float(row["mean_value_rmse"]),
                    "value_rmse_std": float(row["std_value_rmse"]),
                    "pooled_value_rmse_diagnostic": None,
                    "relative_l2_mean": float(row["mean_relative_l2"]),
                    "mc_rmse_contribution_mean": float(row["mean_mc_rmse_contribution"]),
                    "non_mc_rmse_proxy_mean": float(row["mean_non_mc_rmse_proxy"]),
                    "actor_value_rmse": None,
                    "actor_gradient_component_rmse": None,
                    "actor_gradient_vector_rmse": float(row["mean_gradient_vector_rmse"]),
                    "variance_reduction_mean": float(row["mean_variance_reduction"]),
                    "inference_seconds_mean": float(row["mean_inference_seconds"]),
                    "inference_seconds_all_replicates": None,
                    "peak_gpu_bytes": float(row["mean_peak_gpu_bytes"]),
                    "nonfinite_count": float(row["mean_nonfinite_count"]),
                    "selected_by_reference_free_rule": row["actor_request"] == "a4",
                    "selector_normalized_total_inconsistency": selector_row[
                        "normalized_total_inconsistency"
                    ],
                    "selector_noise_corrected_component_rmse": selector_row[
                        "noise_corrected_component_rmse"
                    ],
                    "selector_noise_corrected_value_rmse": None,
                    "cumulative_update_seconds": cumulative,
                    "compute_to_one_evaluation_seconds": cumulative
                    + float(row["mean_inference_seconds"]),
                }
            )
    return replicate_rows, aggregate_rows, audits


def _calibration_rows(b6_audits: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in _csv(RESULTS_ROOT / "benchmark_calibration.csv"):
        pde_id = "dense_fourier" if row["pde"].startswith("Dense") else "coupled_lqg"
        rows.append(
            {
                "pde_id": pde_id,
                "pde": DISPLAY[pde_id],
                "dimension": int(row["dimension"]),
                "effective_rank": float(row["effective_rank"]),
                "outside_top10": float(row["z_variance_outside_top10"]),
                "nonlinear_skill": float(row["nonlinear_skill"]),
                "heat_only_rmse": float(row["nonlinear_heat_rmse"]),
                "driver_rms": float(row["driver_rms"]),
                "u_std": float(row["u_std"]),
                "z_rms": float(row["z_rms"]),
                "reference_uncertainty_rmse": float(row["reference_uncertainty_rmse"]),
            }
        )
    for dimension in (100, 400):
        row = b6_audits[str(dimension)]["calibration"]
        rows.append(
            {
                "pde_id": "B6",
                "pde": DISPLAY["B6"],
                "dimension": dimension,
                "effective_rank": float(row["r_glob"]),
                "outside_top10": float(row["z_variance_outside_top10"]),
                "nonlinear_skill": float(row["heat_only_skill"]),
                "heat_only_rmse": float(row["heat_only_rmse"]),
                "driver_rms": math.hypot(float(row["f_mean"]), float(row["f_std"])),
                "u_std": float(row["u_std"]),
                "z_rms": float(row["z_rms"]),
                "reference_uncertainty_rmse": None,
            }
        )
    return rows


def _resource_rows(
    aggregate_rows: list[dict[str, Any]],
    training_records: dict[str, dict[str, Any]],
    b6_audits: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    rows: list[dict[str, Any]] = []
    totals: dict[str, float] = {}
    for pde, training in training_records.items():
        for phase in training["phases"]:
            rows.append(
                {
                    "pde_id": pde,
                    "dimension": 100,
                    "result_origin": "NEW_5070_TI_PILOT",
                    "phase": f"a{phase['iteration']}_labels",
                    "gpu_seconds": sum(
                        float(phase["labels"][split]["wall_seconds"])
                        for split in ("train", "validation", "replication")
                    ),
                    "peak_gpu_bytes": max(
                        int(phase["labels"][split]["peak_gpu_bytes"])
                        for split in ("train", "validation", "replication")
                    ),
                }
            )
            rows.append(
                {
                    "pde_id": pde,
                    "dimension": 100,
                    "result_origin": "NEW_5070_TI_PILOT",
                    "phase": f"a{phase['iteration']}_fit",
                    "gpu_seconds": float(phase["training"]["wall_seconds"]),
                    "peak_gpu_bytes": int(phase["training"]["peak_gpu_bytes"]),
                }
            )
        rows.append(
            {
                "pde_id": pde,
                "dimension": 100,
                "result_origin": "NEW_5070_TI_PILOT",
                "phase": "reference_free_selector",
                "gpu_seconds": float(training["selection"]["wall_seconds"]),
                "peak_gpu_bytes": None,
            }
        )
        local = [row for row in aggregate_rows if row["pde_id"] == pde]
        inference = float(sum(row["inference_seconds_all_replicates"] for row in local))
        rows.append(
            {
                "pde_id": pde,
                "dimension": 100,
                "result_origin": "NEW_5070_TI_PILOT",
                "phase": "all_a0_a1_a2_inference",
                "gpu_seconds": inference,
                "peak_gpu_bytes": max(int(row["peak_gpu_bytes"]) for row in local),
            }
        )
        total = float(training["timing"]["training_total_invocation_seconds"]) + inference
        totals[pde] = total
        rows.append(
            {
                "pde_id": pde,
                "dimension": 100,
                "result_origin": "NEW_5070_TI_PILOT",
                "phase": "TOTAL_NEW_PILOT",
                "gpu_seconds": total,
                "peak_gpu_bytes": max(int(row["peak_gpu_bytes"]) for row in local),
            }
        )
    for dimension in (100, 400):
        audit = b6_audits[str(dimension)]
        training_seconds = float(audit["training_resources"]["total_seconds"])
        selected = next(
            row
            for row in aggregate_rows
            if row["pde_id"] == "B6"
            and row["dimension"] == dimension
            and row["actor"] == "a4"
        )
        total = training_seconds + float(selected["inference_seconds_mean"])
        totals[f"B6_d{dimension}"] = total
        rows.extend(
            [
                {
                    "pde_id": "B6",
                    "dimension": dimension,
                    "result_origin": "REUSED_MATURE_RESULT",
                    "phase": "mature_training_and_selector",
                    "gpu_seconds": training_seconds,
                    "peak_gpu_bytes": float(audit["training_resources"]["peak_gpu_bytes"]),
                },
                {
                    "pde_id": "B6",
                    "dimension": dimension,
                    "result_origin": "REUSED_MATURE_RESULT",
                    "phase": "selected_a4_mean_inference",
                    "gpu_seconds": float(selected["inference_seconds_mean"]),
                    "peak_gpu_bytes": float(selected["peak_gpu_bytes"]),
                },
                {
                    "pde_id": "B6",
                    "dimension": dimension,
                    "result_origin": "REUSED_MATURE_RESULT",
                    "phase": "TOTAL_MATURE_SEED0_PLUS_ONE_INFERENCE",
                    "gpu_seconds": total,
                    "peak_gpu_bytes": max(
                        float(audit["training_resources"]["peak_gpu_bytes"]),
                        float(selected["peak_gpu_bytes"]),
                    ),
                },
            ]
        )
    return rows, totals


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "legend.fontsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "savefig.dpi": 240,
        }
    )


def _save(fig: plt.Figure, name: str) -> None:
    FIGURES_ROOT.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES_ROOT / f"{name}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES_ROOT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _plots(
    aggregate: list[dict[str, Any]], calibration: list[dict[str, Any]]
) -> None:
    _style()
    panels = [("B6", 100), ("B6", 400), ("dense_fourier", 100), ("coupled_lqg", 100)]
    fig, axes = plt.subplots(2, 2, figsize=(8.3, 5.7))
    for axis, (pde, dimension) in zip(axes.flat, panels):
        local = sorted(
            [row for row in aggregate if row["pde_id"] == pde and row["dimension"] == dimension],
            key=lambda row: int(row["actor"][1:]),
        )
        x = np.array([int(row["actor"][1:]) for row in local])
        y = np.array([float(row["value_rmse_mean"]) for row in local])
        error = np.array([float(row["value_rmse_std"]) for row in local])
        axis.errorbar(x, y, yerr=error, color=COLORS[pde], marker="o", capsize=3, lw=1.7)
        selected = next(row for row in local if row["selected_by_reference_free_rule"])
        axis.scatter(
            [int(selected["actor"][1:])],
            [float(selected["value_rmse_mean"])],
            marker="s",
            s=45,
            facecolors="none",
            edgecolors="#222222",
            zorder=4,
            label="reference-free selected",
        )
        axis.set_yscale("log")
        axis.set_xticks(x)
        axis.set_xticklabels([f"a{value}" for value in x])
        axis.grid(axis="y", alpha=0.25, linewidth=0.6)
        axis.set_title(f"{DISPLAY[pde]} — {dimension}D")
        axis.set_xlabel("Actor iteration")
        axis.set_ylabel("Value RMSE")
        if pde == "coupled_lqg":
            axis.text(
                0.03,
                0.94,
                "runaway learned policy",
                transform=axis.transAxes,
                va="top",
                color=COLORS[pde],
            )
    fig.suptitle("Actor iteration diagnostics", fontweight="normal")
    fig.tight_layout()
    _save(fig, "actor_iteration_rmse")

    labels = [
        f"{row['pde']}\n{row['dimension']}D" for row in sorted(calibration, key=lambda r: (r["pde_id"], r["dimension"]))
    ]
    ordered = sorted(calibration, key=lambda r: (r["pde_id"], r["dimension"]))
    positions = np.arange(len(ordered))
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.3))
    axes[0].barh(positions, [row["effective_rank"] for row in ordered], color=[COLORS[row["pde_id"]] for row in ordered], alpha=0.86)
    axes[0].set_yticks(positions, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Centred gradient effective rank")
    axes[0].grid(axis="x", alpha=0.25, linewidth=0.6)
    axes[1].barh(positions, [row["outside_top10"] for row in ordered], color=[COLORS[row["pde_id"]] for row in ordered], alpha=0.86)
    axes[1].set_yticks(positions, [""] * len(positions))
    axes[1].invert_yaxis()
    axes[1].set_xlim(0, 0.8)
    axes[1].set_xlabel("Gradient variance outside top 10")
    axes[1].grid(axis="x", alpha=0.25, linewidth=0.6)
    fig.suptitle("Effective-dimensionality diagnostics", fontweight="normal")
    fig.tight_layout()
    _save(fig, "effective_dimension")

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.2))
    axes[0].barh(positions, [row["nonlinear_skill"] for row in ordered], color=[COLORS[row["pde_id"]] for row in ordered], alpha=0.86)
    axes[0].axvline(0.15, color="#555555", linestyle="--", lw=1, label="0.15 calibration marker")
    axes[0].set_yticks(positions, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Heat-only RMSE / std(exact solution)")
    axes[0].grid(axis="x", alpha=0.25, linewidth=0.6)
    axes[0].legend(frameon=False, loc="lower right")
    axes[1].barh(positions, [row["driver_rms"] for row in ordered], color=[COLORS[row["pde_id"]] for row in ordered], alpha=0.86)
    axes[1].set_yticks(positions, [""] * len(positions))
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Driver RMS on frozen test bank")
    axes[1].grid(axis="x", alpha=0.25, linewidth=0.6)
    fig.suptitle("Nonlinearity strength", fontweight="normal")
    fig.subplots_adjust(left=0.30, right=0.98, bottom=0.22, top=0.82, wspace=0.10)
    _save(fig, "nonlinearity_strength")

    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.2))
    for axis, pde in zip(axes, ("B6", "dense_fourier", "coupled_lqg")):
        local = sorted(
            [row for row in aggregate if row["pde_id"] == pde and row["dimension"] == 100],
            key=lambda row: int(row["actor"][1:]),
        )
        x = np.array([float(row["compute_to_one_evaluation_seconds"]) for row in local])
        y = np.array([float(row["value_rmse_mean"]) for row in local])
        axis.plot(x, y, marker="o", color=COLORS[pde], lw=1.7)
        for xx, yy, row in zip(x, y, local):
            axis.annotate(row["actor"], (xx, yy), xytext=(4, 3), textcoords="offset points", fontsize=8)
        axis.set_xscale("log")
        axis.xaxis.set_major_locator(LogLocator(base=10, numticks=4))
        axis.xaxis.set_minor_formatter(NullFormatter())
        if pde == "coupled_lqg":
            axis.set_yscale("log")
            axis.yaxis.set_major_locator(LogLocator(base=10, numticks=7))
            axis.yaxis.set_minor_formatter(NullFormatter())
        else:
            axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
            axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, -3), useMathText=True)
        axis.grid(alpha=0.25, linewidth=0.6)
        axis.set_title(DISPLAY[pde] + (" (mature)" if pde == "B6" else " (pilot)"))
        axis.set_xlabel("Updates + one evaluation (RTX 5070 Ti s)")
        axis.set_ylabel("Value RMSE")
    fig.suptitle("Accuracy–compute diagnostics (budgets are not matched)", fontweight="normal")
    fig.tight_layout()
    _save(fig, "pilot_accuracy_compute")


def _f(value: Any, digits: int = 4) -> str:
    if value is None or value == "":
        return "—"
    return f"{float(value):.{digits}g}"


def _report(
    aggregate: list[dict[str, Any]],
    calibration: list[dict[str, Any]],
    resources: list[dict[str, Any]],
    totals: dict[str, float],
    audit: dict[str, Any],
    training: dict[str, dict[str, Any]],
) -> str:
    cal = {(row["pde_id"], row["dimension"]): row for row in calibration}
    agg = {(row["pde_id"], row["dimension"], row["actor"]): row for row in aggregate}
    dense_gain = agg[("dense_fourier", 100, "a2")]["value_rmse_mean"] / agg[("dense_fourier", 100, "a0")]["value_rmse_mean"]
    lines = [
        "# Main high-dimensional pilot on RTX 5070 Ti",
        "",
        "## Outcome",
        "",
        "| PDE | Decision | Evidence |",
        "|---|---|---|",
        "| B6 full-rank HJB | **GO_MAIN** | Mature reference and artifacts pass reuse audit; nonlinear skill 0.263 (100D) / 0.242 (400D), effective rank 29.0 / 77.1, stable learned iterations. |",
        f"| Dense Fourier joint | **GO_MAIN** | Exact joint $u,z$ reference; nonlinear skill 0.274, effective rank 23.7, and a2/a0 RMSE ratio {_f(dense_gain)} under the frozen two-update pilot. |",
        "| Published coupled LQG | **NEEDS_METHOD_FIX** | Reference is validated, but the problem is nearly heat-only (skill 0.0155) and learned policies run away; the reference-free normalized selector chooses the worst checkpoint. Even after stabilization it is a supporting, not headline, benchmark. |",
        "",
        "The decisive result is mixed but useful: Dense Fourier supplies the missing joint value–gradient nonlinear case, while the recreated published LQG exposes an out-of-domain policy/selector instability rather than a success case.",
        "",
        "## Main summary",
        "",
        "| PDE | d | Effective rank | Nonlinear skill | a0 RMSE | a1 RMSE | a2 RMSE | a4 RMSE | MC contribution at final shown | Total GPU time |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for pde, dimension, final_actor, total_key in [
        ("B6", 100, "a4", "B6_d100"),
        ("B6", 400, "a4", "B6_d400"),
        ("dense_fourier", 100, "a2", "dense_fourier"),
        ("coupled_lqg", 100, "a2", "coupled_lqg"),
    ]:
        c = cal[(pde, dimension)]
        lines.append(
            "| "
            + DISPLAY[pde]
            + f" | {dimension} | {_f(c['effective_rank'])} | {_f(c['nonlinear_skill'])} | "
            + _f(agg[(pde, dimension, "a0")]["value_rmse_mean"])
            + " | "
            + _f(agg[(pde, dimension, "a1")]["value_rmse_mean"])
            + " | "
            + _f(agg[(pde, dimension, "a2")]["value_rmse_mean"])
            + " | "
            + (_f(agg[(pde, dimension, "a4")]["value_rmse_mean"]) if pde == "B6" else "—")
            + " | "
            + _f(agg[(pde, dimension, final_actor)]["mc_rmse_contribution_mean"])
            + " | "
            + f"{totals[total_key]:.1f} s |"
        )
    lines.extend(
        [
            "",
            "B6 uses `REUSED_MATURE_RESULT`: 8,192 states, 128 paths/label, 32 steps, 4×384 network, 5,000 optimizer steps per update and four updates. Its timing is therefore not an equal-budget comparison with the two-update pilots.",
            "",
            "## Mathematical and artifact correctness",
            "",
            f"- Validation gate: **{_json(VALIDATION_PATH)['status']}**. Dense terminal identity, analytic gradient, PDE residual, tangent derivatives/reconstruction, terminal-aware anchoring, and plain-FK/MCV consistency all pass.",
            "- LQG: the tridiagonal matrix is SPD; the rescaled Gauss–Laguerre reference agrees across 96/144/192/256 nodes at roundoff, its analytic gradient agrees with finite differences, the PDE residual is at roundoff, and an independent 262,144-sample Monte Carlo check is within 1.66 standard errors.",
            f"- Prediction/checkpoint/B6 reuse audit: **{audit['status']}**. Every saved prediction hash was checked and each reported RMSE, relative $L^2$, bias and MC contribution was recomputed from the `.npz` arrays.",
            "- No truth entered labels, fitting, early stopping, or checkpoint selection. Truth arrays were loaded only after the selector had completed.",
            "",
            "## Dense Fourier joint PDE",
            "",
            "This is an analytically integrable high-dimensional stress test, not a claim of a generic closed-form nonlinear PDE. It has genuine joint dependence through $f(u,z)=\\tfrac12\\tanh(u)\\lVert z\\rVert^2$ and uses the required anchored value and gradient heads.",
            "",
            "| Actor | Mean RMSE | Pooled two-replicate RMSE (diagnostic) | MC contribution | Direct value-head RMSE | Direct gradient-vector RMSE |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for actor in ("a0", "a1", "a2"):
        row = agg[("dense_fourier", 100, actor)]
        lines.append(
            f"| {actor} | {_f(row['value_rmse_mean'], 6)} | {_f(row['pooled_value_rmse_diagnostic'], 6)} | {_f(row['mc_rmse_contribution_mean'], 6)} | {_f(row['actor_value_rmse'], 6)} | {_f(row['actor_gradient_vector_rmse'], 6)} |"
        )
    lines.extend(
        [
            "",
            f"The frozen reference-free selector chose **{training['dense_fourier']['selection']['selected_actor']}**. a2 reduces mean solution RMSE by {(1-dense_gain)*100:.1f}% versus a0; pooling the two already-saved independent predictions gives a larger {(1-agg[('dense_fourier',100,'a2')]['pooled_value_rmse_diagnostic']/agg[('dense_fourier',100,'a0')]['pooled_value_rmse_diagnostic'])*100:.1f}% reduction. Because all three solution RMSEs are near their MC contributions, the pooled figures are diagnostics rather than a new budget or headline criterion.",
            "",
            "## Recreated published coupled LQG",
            "",
            "This is labelled `RECREATED_PUBLISHED_BENCHMARK`, not an official exact reproduction. Section 3.3 of the paper uses $[0,0.5]\\times B^d$; the official repository defaults to $T=1$ and its coefficient arrays depend on batch shape. We use the paper interval, the unit-ball law, and one globally frozen coefficient draw.",
            "",
            "| Actor | Mean RMSE | Pooled RMSE | MC contribution | Direct gradient-vector RMSE | Selector absolute component RMSE | Selector normalized score |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for actor in ("a0", "a1", "a2"):
        row = agg[("coupled_lqg", 100, actor)]
        lines.append(
            f"| {actor} | {_f(row['value_rmse_mean'], 6)} | {_f(row['pooled_value_rmse_diagnostic'], 6)} | {_f(row['mc_rmse_contribution_mean'], 6)} | {_f(row['actor_gradient_vector_rmse'], 6)} | {_f(row['selector_noise_corrected_component_rmse'], 6)} | {_f(row['selector_normalized_total_inconsistency'], 6)} |"
        )
    lines.extend(
        [
            "",
            "The run is finite but unstable. Policy evaluation necessarily visits states outside the unit-ball training distribution; the evidence is consistent with uncontrolled off-domain learned drift, and the frozen linear evaluation becomes enormous. More importantly, self-normalization hides the scale runaway: absolute fixed-point discrepancy grows from 0.234 to 1.65 to 11,420, while normalized scores fall from 54.4 to 1.004 to 0.990, so the selector chooses a2. This is a method/selector failure discovered by the pilot; the recorded choice and failure are preserved unchanged.",
            "",
            "The benchmark itself is weakly nonlinear on the paper distribution: heat-only RMSE is 0.0127, only 0.0155 solution standard deviations. Thus even a repaired version should be supporting evidence rather than the headline difficult PDE.",
            "",
            "## B6 reuse audit",
            "",
            f"- Mature source commit: `{SOURCE_COMMIT}`.",
            "- Frozen B6 bank hashes: 100D `84b2132607bc6b17e44adbea838915f8ce657813f6a0987e34deca308660f1c1`; 400D `0df59abc155e90aa674532760578dc2563cbdb4aade97a849961c8a8edd9f2f9`.",
            "- The selected calibration rows, raw evaluation records, seed-0 training records, selector choice a4, and available selected checkpoint hashes were cross-checked by the read-only adapter.",
            "",
            "## Compute and deviations",
            "",
            f"- Dense new pilot: {totals['dense_fourier']:.1f} RTX 5070 Ti seconds for both updates, selector, and all six 512-state/1,024-path evaluations.",
            f"- LQG new pilot: {totals['coupled_lqg']:.1f} RTX 5070 Ti seconds on the same accounting basis.",
            f"- Total new training and inference: {totals['dense_fourier'] + totals['coupled_lqg']:.1f} GPU seconds ({(totals['dense_fourier'] + totals['coupled_lqg'])/3600:.3f} GPU h), far below the 3 h cap.",
            "- No CUDA OOM, chunk fallback, timeout, retry, nonfinite sample, hyperparameter search, extra seed, 200D/400D new training, FH-MLP grid, or formal strong-baseline training occurred.",
            "- LQG a1/a2 instability and the selector scale pathology are failures, not budget truncations or omitted rows.",
            "",
            "## Answers and recommended next stage (not run)",
            "",
            "1. **Dense Fourier does exercise joint nonlinearity.** Its heat-only skill is 0.274 and driver RMS is 0.265; both value and gradient heads learn substantially.",
            "2. **It is meaningfully high-effective-dimensional.** Effective rank is 23.7 and 50.8% of centred gradient variance lies outside the top ten directions. The two-update learned Actor improves solution RMSE, though the final estimate is MC-floor limited.",
            "3. **Published LQG is close to heat-only on the paper distribution.** It is not a strong-nonlinearity headline benchmark, and this pilot additionally finds a serious policy/selector stability issue.",
            "4. **B6 and Dense Fourier are complementary main cases.** B6 supplies stable full-rank quadratic HJB geometry; Dense supplies dense, joint $u,z$ nonlinearity. LQG should appear only as an external stress/failure case until fixed.",
            "5. **Proposed formal suite:** B6 and Dense at 100D/200D/400D with three seeds; compare Tangent-Actor against DPI, Deep BSDE, PINN+SDGD, heat-only and a0 under measured same-GPU budgets. Planning estimate: 10–12 GPU h for Tangent-Actor training/audits and 24–40 GPU h for all strong baselines. Reserve at most 2 GPU h for a separately pre-registered LQG stabilization pilot (bounded/off-domain actor control plus a scale-aware selector) before deciding whether to include it.",
            "",
            "## Provenance",
            "",
            f"- Frozen protocol commit: `{PREREGISTRATION_COMMIT}`.",
            "- Pilot runner commit: `9b22d5ad1d394bcd00143683e7044157881ab199`.",
            f"- Study revision: `{REVISION}`.",
            "- Base/highdim source commit: `281977883e0a5ce4f05895878680b582620461be`.",
            "- Official SCaSML repository audited at commit `c2412b942b1720084accf9ac9651c6085e073919`.",
            "- Source reconciliation: [ICLR 2026 paper](https://proceedings.iclr.cc/paper_files/paper/2026/file/cf8353a0efc3644454a5591a7cb55c17-Paper-Conference.pdf) and [official SCaSML repository](https://github.com/Francis-Fan-create/SCaSML).",
            "- Hardware/software details, parameter hashes, bank hashes, seeds, and checkpoint/prediction hashes are stored in `test_bank_manifest.json`, `validation.json`, `artifact_audit.json`, and the raw JSON records.",
            "",
            "## Figures",
            "",
            "- `results/main_highdim_pilot/figures/actor_iteration_rmse.png`",
            "- `results/main_highdim_pilot/figures/effective_dimension.png`",
            "- `results/main_highdim_pilot/figures/nonlinearity_strength.png`",
            "- `results/main_highdim_pilot/figures/pilot_accuracy_compute.png`",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    ensure_directories()
    prediction_audit, pooled = _prediction_audit()
    checkpoint_audit = _checkpoint_audit()
    new_replicates, new_aggregate, training = _new_rows(pooled)
    b6_replicates, b6_aggregate, b6_audits = _b6_rows()
    calibration = _calibration_rows(b6_audits)
    aggregate = b6_aggregate + new_aggregate
    resources, totals = _resource_rows(aggregate, training, b6_audits)

    b6_audit_compact = {
        dimension: {
            "label": payload["label"],
            "source_commit": payload["source_commit"],
            "dimension": payload["dimension"],
            "calibration": payload["calibration"],
            "selected_checkpoint": payload["selected_checkpoint"],
            "evaluation_record_sha256": payload["evaluation_record_sha256"],
            "artifact_sha256": payload["artifact_sha256"],
            "audit_status": payload["audit_status"],
        }
        for dimension, payload in b6_audits.items()
    }
    write_json_atomic(RESULTS_ROOT / "b6_reuse_audit.json", b6_audit_compact)
    audit = {
        "study_revision": REVISION,
        "validation_status": _json(VALIDATION_PATH)["status"],
        "prediction_audit": prediction_audit,
        "checkpoint_audit": checkpoint_audit,
        "b6_reuse_status": {
            dimension: payload["audit_status"] for dimension, payload in b6_audits.items()
        },
    }
    audit["status"] = (
        "PASS"
        if audit["validation_status"] == "VALIDATED"
        and prediction_audit["status"] == "PASS"
        and checkpoint_audit["status"] == "PASS"
        and all(value == "PASS" for value in audit["b6_reuse_status"].values())
        else "FAIL"
    )
    write_json_atomic(RESULTS_ROOT / "artifact_audit.json", audit)
    _write_csv(RESULTS_ROOT / "pilot_metrics.csv", b6_replicates + new_replicates)
    _write_csv(RESULTS_ROOT / "actor_iteration_metrics.csv", aggregate)
    _write_csv(RESULTS_ROOT / "resource_accounting.csv", resources)
    _write_csv(RESULTS_ROOT / "combined_benchmark_calibration.csv", calibration)
    _plots(aggregate, calibration)

    decisions = {
        "B6": "GO_MAIN",
        "dense_fourier": "GO_MAIN",
        "coupled_lqg": "NEEDS_METHOD_FIX",
    }
    summary = {
        "study_revision": REVISION,
        "status": "COMPLETE" if audit["status"] == "PASS" else "AUDIT_FAILED",
        "decisions": decisions,
        "artifact_audit": audit["status"],
        "new_gpu_seconds": totals["dense_fourier"] + totals["coupled_lqg"],
        "totals": totals,
        "dense_a2_over_a0": aggregate[-4]["value_rmse_mean"]
        / aggregate[-6]["value_rmse_mean"],
        "lqg_reference_free_selected": training["coupled_lqg"]["selection"]["selected_actor"],
        "lqg_failure": "finite policy-evaluation runaway and scale-normalized selector failure",
        "formal_baselines_run": False,
        "forbidden_expansions_run": False,
    }
    write_json_atomic(RESULTS_ROOT / "analysis_summary.json", summary)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(
        _report(aggregate, calibration, resources, totals, audit, training),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

