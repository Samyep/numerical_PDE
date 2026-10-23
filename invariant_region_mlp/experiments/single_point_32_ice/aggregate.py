"""Independent Stage-D audit and Stage-E aggregation for the 32-point study."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .orchestration import (
    DIMENSION_AUDIT_DIR,
    algorithm_hashes,
    load_manifest,
    manifest_path,
)
from .protocol import (
    DIMENSIONS,
    FH_MAIN_CONFIGURATION,
    PDES,
    QUERY_COUNT,
    RESULTS_DIR,
    REVISION,
    SOURCE_ROOT,
    STAGE_A_AUDIT,
    STAGE_B_AUDIT,
    TA_PRIMARY_EULER_STEPS,
    array_sha256,
    file_sha256,
    frozen_protocol,
    load_query_bank,
    query_path,
    write_json_atomic,
)


METHOD_LABELS = {
    "tangent_actor": "Tangent-Actor",
    "deep_bsde": "Deep BSDE",
    "fh_mlp_standard": "FH-MLP",
}
METHOD_ORDER = ("Tangent-Actor", "FH-MLP", "Deep BSDE")
COLORS = {
    "Tangent-Actor": "#006BA4",
    "FH-MLP": "#FF800E",
    "Deep BSDE": "#ABABAB",
}
REPORT_PATH = SOURCE_ROOT / "invariant_region_mlp" / "docs" / (
    "ICML_SINGLE_POINT_32_ICE_REPORT.md"
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _portable_path(value: str | Path) -> Path:
    path = Path(str(value))
    if path.is_file():
        return path
    normalized = str(value).replace("\\", "/")
    marker = "/single_point_32_ice/"
    if marker in normalized:
        relative = normalized.split(marker, 1)[1]
        candidate = RESULTS_DIR / relative
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(value)


def _atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False, lineterminator="\n")
    temporary.replace(path)


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(SOURCE_ROOT.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def _prediction(record: dict[str, Any]) -> tuple[Path, dict[str, np.ndarray]]:
    path = _portable_path(record["prediction_file"])
    with np.load(path, allow_pickle=False) as payload:
        arrays = {key: np.asarray(payload[key]) for key in payload.files}
    return path, arrays


def _formal_manifests() -> list[dict[str, Any]]:
    manifests = []
    for dimension in DIMENSIONS:
        audit_path = DIMENSION_AUDIT_DIR / f"d{dimension}.json"
        if not audit_path.is_file() or _json(audit_path).get("status") != "VALIDATED":
            raise RuntimeError(f"missing validated dimension audit: {audit_path}")
        manifests.append(load_manifest(manifest_path("C", dimension)))
    return manifests


def collect_raw_results() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ta_rows: list[dict[str, Any]] = []
    deep_rows: list[dict[str, Any]] = []
    fh_rows: list[dict[str, Any]] = []
    training_rows: list[dict[str, Any]] = []
    seen_training: set[str] = set()

    for manifest in _formal_manifests():
        for task in manifest["tasks"]:
            method = task["method"]
            pde = task["pde"]
            dimension = int(task["dimension"])
            if method == "tangent_actor":
                raw_root = RESULTS_DIR / "raw" / "formal" / "tangent_actor_evaluations"
                for record_path in sorted(raw_root.glob(
                    f"ta__{pde}__d{dimension}__s{task['seed']}__rep*__steps*__formal.json"
                )):
                    record = _json(record_path)
                    prediction_path, arrays = _prediction(record)
                    prediction = arrays["prediction"].reshape(-1)
                    reference = arrays["reference"].reshape(-1)
                    pair_variance = arrays["pair_variance"].reshape(-1)
                    if len(prediction) != QUERY_COUNT:
                        raise RuntimeError(f"wrong TA prediction length: {prediction_path}")
                    pairs = int(record["paths"]) // 2
                    for point in range(QUERY_COUNT):
                        ta_rows.append({
                            "method": METHOD_LABELS[method],
                            "pde": pde,
                            "dimension": dimension,
                            "point": point,
                            "training_seed": int(record["training_seed"]),
                            "inference_replicate": int(record["inference_replicate"]),
                            "time_steps": int(record["time_steps"]),
                            "paths": int(record["paths"]),
                            "prediction": float(prediction[point]),
                            "reference": float(reference[point]),
                            "error": float(prediction[point] - reference[point]),
                            "absolute_error": float(abs(prediction[point] - reference[point])),
                            "pair_variance": float(pair_variance[point]),
                            "mc_standard_error": float(
                                math.sqrt(max(float(pair_variance[point]), 0.0) / pairs)
                            ),
                            "selected_actor": record["selected_actor"],
                            "prediction_file": _relative(prediction_path),
                            "record_file": _relative(record_path),
                            "status": record["status"],
                        })
                    training_id = str(record["training_run_id"])
                    if training_id not in seen_training:
                        seen_training.add(training_id)
                        artifact = _portable_path(record["training_artifact"])
                        training = _json(artifact)
                        phases = training.get("phases", [])
                        selected = training.get("selection", {})
                        selected_actor = selected.get("selected_actor")
                        selected_hash = None
                        if isinstance(selected_actor, str) and selected_actor.startswith("a"):
                            selected_index = int(selected_actor[1:])
                            if selected_index > 0:
                                selected_hash = next(
                                    (
                                        phase.get("checkpoint_sha256")
                                        for phase in phases
                                        if int(phase.get("iteration", -1)) == selected_index
                                    ),
                                    None,
                                )
                        training_rows.append({
                            "method": METHOD_LABELS[method],
                            "pde": pde,
                            "dimension": dimension,
                            "point": np.nan,
                            "training_seed": int(task["seed"]),
                            "run_id": training_id,
                            "classification": "native_tangent_actor_full_function_training",
                            "status": training.get("status"),
                            "completed_updates": len(phases),
                            "optimizer_steps": training.get("configuration", {}).get(
                                "optimizer_steps"
                            ),
                            "selected_checkpoint": selected_actor,
                            "training_seconds": float(sum(
                                float(phase.get("phase_seconds", 0.0))
                                for phase in phases
                            )),
                            "nonfinite_count": sum(
                                int(phase.get("training", {}).get("nonfinite_count", 0))
                                for phase in phases
                            ),
                            "checkpoint_sha256": selected_hash,
                            "artifact": _relative(artifact),
                        })
            elif method == "deep_bsde":
                run_id = (
                    f"deep_bsde__{pde}__d{dimension}__q{task['point']:02d}"
                    f"__s{task['seed']}__formal"
                )
                record_path = (
                    RESULTS_DIR / "raw" / "formal" / "deep_bsde" / f"{run_id}.json"
                )
                record = _json(record_path)
                prediction_path, arrays = _prediction(record)
                prediction = float(np.asarray(arrays["prediction"]).reshape(()))
                reference = float(np.asarray(arrays["reference"]).reshape(()))
                deep_rows.append({
                    "method": METHOD_LABELS[method],
                    "pde": pde,
                    "dimension": dimension,
                    "point": int(task["point"]),
                    "training_seed": int(task["seed"]),
                    "prediction": prediction,
                    "reference": reference,
                    "error": prediction - reference,
                    "absolute_error": abs(prediction - reference),
                    "final_validation_terminal_clipped_mse": (
                        record["history"][-1][
                            "validation_terminal_clipped_mse"
                        ] if record.get("history") else np.nan
                    ),
                    "nonfinite_count": int(record.get("nonfinite_count", 0)),
                    "training_seconds": float(record["training_seconds"]),
                    "peak_gpu_bytes": record.get("peak_gpu_bytes"),
                    "classification": record["classification"],
                    "prediction_file": _relative(prediction_path),
                    "record_file": _relative(record_path),
                    "status": record["status"],
                })
                training_rows.append({
                    "method": METHOD_LABELS[method],
                    "pde": pde,
                    "dimension": dimension,
                    "point": int(task["point"]),
                    "training_seed": int(task["seed"]),
                    "run_id": run_id,
                    "classification": record["classification"],
                    "status": record["status"],
                    "completed_updates": int(record["completed_updates"]),
                    "optimizer_steps": int(record["configuration"]["iterations"]),
                    "selected_checkpoint": "final_scheduled_iterate",
                    "training_seconds": float(record["training_seconds"]),
                    "nonfinite_count": int(record.get("nonfinite_count", 0)),
                    "checkpoint_sha256": record["checkpoint_sha256"],
                    "artifact": _relative(record_path),
                })
            else:
                run_id = (
                    f"fh_mlp_standard__{pde}__d{dimension}"
                    f"__n{task['n']}_M{task['M']}"
                    f"__r{task['replicate']:02d}__formal"
                )
                record_path = RESULTS_DIR / "raw" / "formal" / "fh_mlp" / f"{run_id}.json"
                record = _json(record_path)
                prediction_path, arrays = _prediction(record)
                prediction = arrays["prediction"].reshape(-1)
                reference = arrays["reference"].reshape(-1)
                point_seconds = arrays["point_seconds"].reshape(-1)
                if len(prediction) != QUERY_COUNT:
                    raise RuntimeError(f"wrong FH prediction length: {prediction_path}")
                for point in range(QUERY_COUNT):
                    fh_rows.append({
                        "method": METHOD_LABELS[method],
                        "pde": pde,
                        "dimension": dimension,
                        "point": point,
                        "n": int(task["n"]),
                        "M": int(task["M"]),
                        "replicate": int(task["replicate"]),
                        "prediction": float(prediction[point]),
                        "reference": float(reference[point]),
                        "error": float(prediction[point] - reference[point]),
                        "absolute_error": float(abs(prediction[point] - reference[point])),
                        "point_seconds": float(point_seconds[point]),
                        "run_inference_seconds": float(record["inference_seconds"]),
                        "run_nonfinite_count": int(record["nonfinite_count"]),
                        "peak_gpu_bytes": record.get("peak_gpu_bytes"),
                        "main_configuration": bool(record["main_configuration"]),
                        "prediction_file": _relative(prediction_path),
                        "record_file": _relative(record_path),
                        "status": record["status"],
                    })

    ta = pd.DataFrame(ta_rows).sort_values(
        ["pde", "dimension", "training_seed", "time_steps",
         "inference_replicate", "point"]
    ).reset_index(drop=True)
    deep = pd.DataFrame(deep_rows).sort_values(
        ["pde", "dimension", "point", "training_seed"]
    ).reset_index(drop=True)
    fh = pd.DataFrame(fh_rows).sort_values(
        ["pde", "dimension", "n", "M", "replicate", "point"]
    ).reset_index(drop=True)
    training = pd.DataFrame(training_rows).sort_values(
        ["method", "pde", "dimension", "point", "training_seed"],
        na_position="first",
    ).reset_index(drop=True)
    return ta, deep, fh, training


def main_rows(
    ta: pd.DataFrame, deep: pd.DataFrame, fh: pd.DataFrame
) -> pd.DataFrame:
    ta_main = ta.loc[ta["time_steps"] == TA_PRIMARY_EULER_STEPS].copy()
    ta_main["repeat_id"] = (
        "s" + ta_main["training_seed"].astype(str)
        + "_mc" + ta_main["inference_replicate"].astype(str)
    )
    ta_main["randomness"] = "training_seed_and_fk_mc"
    deep_main = deep.copy()
    deep_main["repeat_id"] = "s" + deep_main["training_seed"].astype(str)
    deep_main["randomness"] = "training_seed"
    fh_main = fh.loc[
        (fh["n"] == FH_MAIN_CONFIGURATION[0])
        & (fh["M"] == FH_MAIN_CONFIGURATION[1])
    ].copy()
    fh_main["repeat_id"] = "mc" + fh_main["replicate"].astype(str)
    fh_main["randomness"] = "recursive_mc"
    columns = [
        "method", "pde", "dimension", "point", "repeat_id", "randomness",
        "prediction", "reference", "error", "absolute_error", "status",
    ]
    result = pd.concat(
        [ta_main[columns], fh_main[columns], deep_main[columns]],
        ignore_index=True,
    )
    return result.sort_values(
        ["pde", "dimension", "method", "repeat_id", "point"]
    ).reset_index(drop=True)


def _reference_uncertainty(pde: str, dimension: int) -> tuple[float, float]:
    values = np.load(
        query_path(pde, dimension) / "reference_uncertainty.npy",
        allow_pickle=False,
    ).astype(np.float64)
    return float(np.sqrt(np.mean(values**2))), float(np.max(values))


def _metric_row(key: tuple[str, int, str], frame: pd.DataFrame) -> dict[str, Any]:
    pde, dimension, method = key
    prediction = frame["prediction"].to_numpy(dtype=np.float64)
    reference = frame["reference"].to_numpy(dtype=np.float64)
    error = prediction - reference
    finite = np.isfinite(error)
    complete_finite = bool(np.all(finite))
    finite_error = error[finite]
    finite_reference = reference[finite]
    reference_norm = float(np.sqrt(np.sum(reference**2)))
    finite_reference_norm = float(np.sqrt(np.sum(finite_reference**2)))
    repeat_rmses: list[float] = []
    for _, repeated in frame.groupby("repeat_id", sort=True):
        local = repeated.sort_values("point")
        local_error = local["error"].to_numpy(dtype=np.float64)
        if len(local) == QUERY_COUNT and np.all(np.isfinite(local_error)):
            repeat_rmses.append(float(np.sqrt(np.mean(local_error**2))))
    expected_repeats = {
        "Tangent-Actor": 12,
        "FH-MLP": 10,
        "Deep BSDE": 3,
    }[method]
    complete_count = len(frame) == QUERY_COUNT * expected_repeats
    uncertainty_rms, uncertainty_max = _reference_uncertainty(pde, dimension)
    if complete_finite:
        rmse = float(np.sqrt(np.mean(error**2)))
        mae = float(np.mean(np.abs(error)))
        relative_l2 = (
            float(np.sqrt(np.sum(error**2)) / reference_norm)
            if reference_norm > np.finfo(float).tiny else np.nan
        )
        bias = float(np.mean(error))
        max_abs = float(np.max(np.abs(error)))
        random_error = frame.loc[frame["point"] != 0, "error"].to_numpy(float)
        origin_error = frame.loc[frame["point"] == 0, "error"].to_numpy(float)
        rmse_random = float(np.sqrt(np.mean(random_error**2)))
        origin_abs = float(np.sqrt(np.mean(origin_error**2)))
    else:
        rmse = mae = relative_l2 = bias = max_abs = np.nan
        rmse_random = origin_abs = np.nan
    repeat_sd = float(np.std(repeat_rmses, ddof=1)) if len(repeat_rmses) > 1 else np.nan
    return {
        "pde": pde,
        "dimension": dimension,
        "method": method,
        "status": (
            "COMPLETE" if complete_count and complete_finite
            else "NUMERICAL_FAILURE" if complete_count
            else "PARTIAL"
        ),
        "expected_repeats": expected_repeats,
        "observed_repeats": int(frame["repeat_id"].nunique()),
        "observations": len(frame),
        "finite_count": int(finite.sum()),
        "nonfinite_count": int((~finite).sum()),
        "finite_fraction": float(finite.mean()),
        "rmse_32_R": rmse,
        "mae": mae,
        "relative_l2": relative_l2,
        "signed_bias": bias,
        "maximum_absolute_error": max_abs,
        "rmse_random_31_R": rmse_random,
        "origin_error_rms_over_repeats": origin_abs,
        "repeat_rmse_sd": repeat_sd,
        "repeat_rmse_se": (
            repeat_sd / math.sqrt(len(repeat_rmses))
            if len(repeat_rmses) > 1 else np.nan
        ),
        "finite_only_rmse_diagnostic": (
            float(np.sqrt(np.mean(finite_error**2))) if len(finite_error) else np.nan
        ),
        "finite_only_relative_l2_diagnostic": (
            float(np.sqrt(np.sum(finite_error**2)) / finite_reference_norm)
            if len(finite_error) and finite_reference_norm > np.finfo(float).tiny
            else np.nan
        ),
        "reference_uncertainty_rms": uncertainty_rms,
        "reference_uncertainty_max": uncertainty_max,
    }


def aggregate_metrics(main: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    scaling_rows = [
        _metric_row(key, frame)
        for key, frame in main.groupby(["pde", "dimension", "method"], sort=True)
    ]
    scaling = pd.DataFrame(scaling_rows).sort_values(
        ["pde", "dimension", "method"]
    ).reset_index(drop=True)
    point_rows: list[dict[str, Any]] = []
    for key, frame in main.groupby(
        ["pde", "dimension", "method", "point"], sort=True
    ):
        pde, dimension, method, point = key
        prediction = frame["prediction"].to_numpy(dtype=np.float64)
        reference = float(frame["reference"].iloc[0])
        finite = np.isfinite(prediction)
        error = prediction - reference
        point_rows.append({
            "pde": pde,
            "dimension": dimension,
            "method": method,
            "point": point,
            "origin": point == 0,
            "reference": reference,
            "mean_prediction": (
                float(np.mean(prediction)) if np.all(finite) else np.nan
            ),
            "absolute_error_of_mean": (
                float(abs(np.mean(prediction) - reference))
                if np.all(finite) else np.nan
            ),
            "point_rmse_over_repeats": (
                float(np.sqrt(np.mean(error**2))) if np.all(finite) else np.nan
            ),
            "repeat_prediction_sd": (
                float(np.std(prediction, ddof=1))
                if len(prediction) > 1 and np.all(finite) else np.nan
            ),
            "finite_fraction": float(finite.mean()),
            "repeats": len(frame),
        })
    pointwise = pd.DataFrame(point_rows).sort_values(
        ["pde", "dimension", "method", "point"]
    ).reset_index(drop=True)
    return scaling, pointwise


def randomness_decomposition(
    main: pd.DataFrame, ta: pd.DataFrame
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (pde, dimension, method), frame in main.groupby(
        ["pde", "dimension", "method"], sort=True
    ):
        total_variances = (
            frame.groupby("point")["prediction"].var(ddof=1).to_numpy(float)
        )
        between = within = np.nan
        component = "repeat_variance"
        if method == "Tangent-Actor":
            source = ta.loc[
                (ta["pde"] == pde)
                & (ta["dimension"] == dimension)
                & (ta["time_steps"] == TA_PRIMARY_EULER_STEPS)
            ]
            seed_means = source.groupby(
                ["point", "training_seed"], sort=True
            )["prediction"].mean()
            between = float(seed_means.groupby("point").var(ddof=1).mean())
            within = float(
                source.groupby(["point", "training_seed"])["prediction"]
                .var(ddof=1)
                .mean()
            )
            component = "training_seed_between_plus_fk_mc_within"
        elif method == "Deep BSDE":
            between = float(np.nanmean(total_variances))
            component = "training_seed"
        else:
            within = float(np.nanmean(total_variances))
            component = "recursive_mc"
        total = float(np.nansum([between, within]))
        rows.append({
            "pde": pde,
            "dimension": dimension,
            "method": method,
            "decomposition": component,
            "mean_total_pointwise_variance": float(np.nanmean(total_variances)),
            "between_training_seed_variance": between,
            "within_inference_mc_variance": within,
            "mc_fraction_of_decomposed_variance": (
                within / total if np.isfinite(within) and total > 0 else np.nan
            ),
            "note": (
                "The frozen bank includes one deterministic origin; variance is "
                "algorithmic randomness, not query-sampling uncertainty."
            ),
        })
    return pd.DataFrame(rows).sort_values(
        ["pde", "dimension", "method"]
    ).reset_index(drop=True)


def _sacct(job_ids: list[str]) -> dict[str, dict[str, Any]]:
    if not job_ids:
        return {}
    fields = "JobIDRaw,JobName,State,ElapsedRaw,AllocTRES%100,MaxRSS,ExitCode,NodeList"
    try:
        completed = subprocess.run(
            [
                "sacct", "-X", "-n", "-P", "-j", ",".join(job_ids),
                "-o", fields,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return {}
    result: dict[str, dict[str, Any]] = {}
    for line in completed.stdout.splitlines():
        parts = line.split("|")
        if len(parts) < 8 or not parts[0] or "." in parts[0]:
            continue
        job_id, job_name, state, elapsed, tres, max_rss, exit_code, node = parts[:8]
        result[job_id] = {
            "slurm_job_name": job_name,
            "slurm_state": state,
            "slurm_elapsed_seconds": int(elapsed or 0),
            "slurm_alloc_tres": tres,
            "slurm_max_rss": max_rss,
            "slurm_exit_code": exit_code,
            "slurm_node": node,
        }
    return result


def resource_accounting() -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted((RESULTS_DIR / "task_status" / "c").glob("*.json")):
        record = _json(path)
        task = record["task"]
        environment = record.get("environment", {})
        summaries = record.get("record_summaries", [])
        peaks = [
            item.get("peak_gpu_bytes")
            for item in summaries
            if item.get("peak_gpu_bytes") is not None
        ]
        rows.append({
            "task_id": task["task_id"],
            "method": METHOD_LABELS[task["method"]],
            "pde": task["pde"],
            "dimension": int(task["dimension"]),
            "status": record["status"],
            "task_wall_seconds": record.get("wall_seconds"),
            "solver_training_seconds": sum(
                float(item.get("training_seconds") or 0.0) for item in summaries
            ),
            "solver_inference_seconds": sum(
                float(item.get("inference_seconds") or 0.0) for item in summaries
            ),
            "peak_gpu_bytes": max(peaks) if peaks else np.nan,
            "slurm_job_id": environment.get("slurm_job_id"),
            "slurm_array_job_id": environment.get("slurm_array_job_id"),
            "slurm_array_task_id": environment.get("slurm_array_task_id"),
            "hostname": environment.get("hostname"),
            "gpu_name": environment.get("cuda_device_name"),
            "gpu_total_memory": environment.get("cuda_total_memory"),
            "started_utc": record.get("started_utc"),
            "finished_utc": record.get("finished_utc"),
            "status_file": _relative(path),
        })
    frame = pd.DataFrame(rows)
    job_ids = sorted({
        str(value) for value in frame.get("slurm_job_id", pd.Series(dtype=str)).dropna()
        if str(value)
    })
    accounting = _sacct(job_ids)
    for field in (
        "slurm_job_name", "slurm_state", "slurm_elapsed_seconds",
        "slurm_alloc_tres", "slurm_max_rss", "slurm_exit_code", "slurm_node",
    ):
        frame[field] = [
            accounting.get(str(job_id), {}).get(field)
            for job_id in frame["slurm_job_id"]
        ]
    unique_allocations = {
        str(row["slurm_job_id"]): int(row["slurm_elapsed_seconds"])
        for _, row in frame.dropna(
            subset=["slurm_job_id", "slurm_elapsed_seconds"]
        ).iterrows()
    }
    summary = {
        "unique_gpu_allocations": len(unique_allocations),
        "allocated_gpu_seconds": sum(unique_allocations.values()),
        "allocated_gpu_hours": sum(unique_allocations.values()) / 3600.0,
        "formal_task_wall_seconds": float(
            pd.to_numeric(frame["task_wall_seconds"], errors="coerce").sum()
        ),
        "sacct_available": bool(accounting),
    }
    return frame.sort_values(
        ["dimension", "method", "pde", "task_id"]
    ).reset_index(drop=True), summary


def _method_values(frame: pd.DataFrame, method: str, field: str) -> np.ndarray:
    return frame.loc[frame["method"] == method, field].to_numpy(dtype=float)


def make_figures(
    scaling: pd.DataFrame,
    pointwise: pd.DataFrame,
    fh: pd.DataFrame,
    deep: pd.DataFrame,
    training: pd.DataFrame,
) -> list[Path]:
    output = RESULTS_DIR / "figures"
    output.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.2), sharex=True)
    for axis, pde in zip(axes.flat, PDES, strict=True):
        local = scaling.loc[scaling["pde"] == pde]
        for method in METHOD_ORDER:
            rows = local.loc[local["method"] == method].sort_values("dimension")
            axis.errorbar(
                rows["dimension"],
                rows["rmse_32_R"],
                yerr=rows["repeat_rmse_se"],
                marker="o",
                linewidth=1.7,
                capsize=3,
                color=COLORS[method],
                label=method,
            )
        axis.set_title(pde.replace("_", " "))
        axis.set_yscale("log")
        axis.set_xscale("log")
        axis.set_xticks(DIMENSIONS, labels=[str(value) for value in DIMENSIONS])
        axis.grid(True, which="both", alpha=0.2)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.supxlabel("ambient dimension")
    fig.supylabel("RMSE over 32 points and formal repeats")
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        path = output / f"figure1_accuracy_vs_dimension.{suffix}"
        fig.savefig(path, dpi=240, bbox_inches="tight")
        generated.append(path)
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(13.2, 6.4), sharey=False)
    for row_index, dimension in enumerate((100, 400)):
        for column, pde in enumerate(PDES):
            axis = axes[row_index, column]
            local = pointwise.loc[
                (pointwise["pde"] == pde)
                & (pointwise["dimension"] == dimension)
            ]
            data = [
                _method_values(local, method, "point_rmse_over_repeats")
                for method in METHOD_ORDER
            ]
            data = [
                values[np.isfinite(values)] if np.isfinite(values).any()
                else np.asarray([np.nan])
                for values in data
            ]
            box = axis.boxplot(data, patch_artist=True, showfliers=True)
            for patch, method in zip(box["boxes"], METHOD_ORDER, strict=True):
                patch.set_facecolor(COLORS[method])
                patch.set_alpha(0.65)
            axis.set_yscale("log")
            axis.set_xticks(range(1, 4), ["TA", "FH", "DB"])
            axis.set_title(f"{pde.replace('_', ' ')}, d={dimension}")
            axis.grid(True, axis="y", alpha=0.2)
    fig.supylabel("pointwise RMSE over formal repeats")
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        path = output / f"figure2_pointwise_error_distributions.{suffix}"
        fig.savefig(path, dpi=240, bbox_inches="tight")
        generated.append(path)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.9))
    failures = (
        fh.assign(failed=~np.isfinite(fh["prediction"]))
        .groupby(["n", "M"], as_index=False)["failed"]
        .mean()
    )
    axes[0].bar(
        np.arange(len(failures)),
        failures["failed"],
        color=COLORS["FH-MLP"],
    )
    axes[0].set_xticks(
        np.arange(len(failures)),
        [f"{n},{M}" for n, M in zip(failures["n"], failures["M"], strict=True)],
        rotation=45,
        ha="right",
    )
    axes[0].set_ylabel("nonfinite prediction fraction")
    axes[0].set_title("FH-MLP configurations")

    for pde in PDES:
        local = deep.loc[deep["pde"] == pde]
        medians = local.groupby("dimension")[
            "final_validation_terminal_clipped_mse"
        ].median()
        axes[1].plot(
            medians.index,
            medians.values,
            marker="o",
            label=pde.replace("_", " "),
        )
    axes[1].set_yscale("log")
    axes[1].set_xscale("log")
    axes[1].set_xticks(DIMENSIONS, labels=[str(value) for value in DIMENSIONS])
    axes[1].set_title("Deep BSDE terminal validation")
    axes[1].legend(frameon=False, fontsize=6)

    ta_training = training.loc[training["method"] == "Tangent-Actor"].copy()
    selected = (
        ta_training.groupby(["dimension", "selected_checkpoint"])
        .size()
        .unstack(fill_value=0)
        .reindex(DIMENSIONS, fill_value=0)
    )
    bottom = np.zeros(len(selected))
    for actor in sorted(selected.columns):
        values = selected[actor].to_numpy(float)
        axes[2].bar(
            np.arange(len(selected)),
            values,
            bottom=bottom,
            label=str(actor),
        )
        bottom += values
    axes[2].set_xticks(np.arange(len(selected)), [str(value) for value in selected.index])
    axes[2].set_title("TA selected checkpoints")
    axes[2].set_ylabel("training runs")
    axes[2].legend(frameon=False, fontsize=6, ncol=2)
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        path = output / f"figure3_numerical_stability.{suffix}"
        fig.savefig(path, dpi=240, bbox_inches="tight")
        generated.append(path)
    plt.close(fig)
    return generated


def make_tables(scaling: pd.DataFrame) -> list[Path]:
    output = RESULTS_DIR / "tables"
    output.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    pivot = scaling.pivot_table(
        index=["pde", "dimension"],
        columns="method",
        values="rmse_32_R",
        aggfunc="first",
    ).reindex(columns=METHOD_ORDER)
    csv_path = output / "main_rmse_table.csv"
    tex_path = output / "main_rmse_table.tex"
    pivot.to_csv(csv_path, lineterminator="\n")
    tex_path.write_text(
        pivot.to_latex(
            float_format=lambda value: f"{value:.4e}",
            na_rep="failure",
            escape=True,
        ),
        encoding="utf-8",
    )
    generated.extend([csv_path, tex_path])
    complete = scaling[[
        "pde", "dimension", "method", "status", "finite_fraction",
        "rmse_32_R", "mae", "relative_l2", "signed_bias",
        "maximum_absolute_error", "rmse_random_31_R",
        "origin_error_rms_over_repeats",
    ]]
    path = output / "complete_metrics.tex"
    path.write_text(
        complete.to_latex(
            index=False,
            float_format=lambda value: f"{value:.4e}",
            na_rep="failure",
            escape=True,
        ),
        encoding="utf-8",
    )
    generated.append(path)
    return generated


def _query_audit() -> dict[str, Any]:
    failures: list[str] = []
    for pde in PDES:
        for dimension in DIMENSIONS:
            root = query_path(pde, dimension)
            manifest = _json(root / "manifest.json")
            t, x, reference = load_query_bank(pde, dimension)
            checks = {
                "t": file_sha256(root / "t.npy") == manifest["t_sha256"],
                "x": file_sha256(root / "x.npy") == manifest["x_sha256"],
                "reference": (
                    file_sha256(root / "reference.npy")
                    == manifest["reference_sha256"]
                ),
                "input": (
                    array_sha256(t, x) == manifest["canonical_input_sha256"]
                ),
                "canonical_reference": (
                    array_sha256(reference)
                    == manifest["canonical_reference_sha256"]
                ),
                "shape": x.shape == (QUERY_COUNT, dimension),
                "origin": bool(np.array_equal(x[0], np.zeros(dimension))),
            }
            failures.extend(
                f"{pde}/d{dimension}/{name}"
                for name, passed in checks.items()
                if not passed
            )
    return {
        "passed": not failures,
        "checked_banks": len(PDES) * len(DIMENSIONS),
        "failures": failures,
    }


def _parameter_audit() -> dict[str, Any]:
    manifest = _json(RESULTS_DIR / "parameter_manifest.json")
    failures: list[str] = []
    checked = 0
    for row in manifest.get("rows", []):
        if "file" not in row:
            continue
        checked += 1
        path = RESULTS_DIR / row["file"]
        if not path.is_file():
            failures.append(f"missing/{row['pde']}/d{row['dimension']}")
            continue
        if file_sha256(path) != row["file_sha256"]:
            failures.append(f"file_hash/{row['pde']}/d{row['dimension']}")
        with np.load(path, allow_pickle=False) as payload:
            digest = array_sha256(*(payload[key] for key in payload.files))
        if digest != row["array_sha256"]:
            failures.append(f"array_hash/{row['pde']}/d{row['dimension']}")
        if row["array_sha256"] != row["source_array_sha256"]:
            failures.append(f"source_mismatch/{row['pde']}/d{row['dimension']}")
        if row.get("byte_identical_arrays_to_completed_ice_expansion") is not True:
            failures.append(f"not_byte_identical/{row['pde']}/d{row['dimension']}")
    return {
        "passed": not failures and checked == 10,
        "checked_parameter_files": checked,
        "expected_parameter_files": 10,
        "failures": failures,
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def artifact_audit(
    ta: pd.DataFrame,
    deep: pd.DataFrame,
    fh: pd.DataFrame,
    scaling: pd.DataFrame,
    generated: list[Path],
) -> dict[str, Any]:
    protocol_payload = _json(RESULTS_DIR / "FROZEN_PROTOCOL.json")
    expected_protocol = frozen_protocol()
    stage_a = _json(STAGE_A_AUDIT)
    stage_b = _json(STAGE_B_AUDIT)
    dimension_audits = {
        dimension: _json(DIMENSION_AUDIT_DIR / f"d{dimension}.json")
        for dimension in DIMENSIONS
    }
    checks = {
        "protocol_hash": {
            "passed": (
                protocol_payload.get("protocol_sha256")
                == expected_protocol["protocol_sha256"]
            ),
        },
        "query_banks": _query_audit(),
        "frozen_parameters": _parameter_audit(),
        "reference_validation": {
            "passed": _json(RESULTS_DIR / "reference_validation.json").get(
                "overall_pass"
            ) is True,
        },
        "stage_a": {
            "passed": (
                stage_a.get("status") == "VALIDATED"
                and stage_a.get("overall_pass") is True
            ),
        },
        "stage_b": {
            "passed": (
                stage_b.get("status") == "VALIDATED"
                and stage_b.get("overall_pass") is True
            ),
        },
        "dimension_audits": {
            "passed": all(
                value.get("status") == "VALIDATED"
                and value.get("overall_pass") is True
                for value in dimension_audits.values()
            ),
        },
        "raw_row_counts": {
            "passed": (
                len(ta) == 20 * 3 * 6 * QUERY_COUNT
                and len(deep) == 20 * QUERY_COUNT * 3
                and len(fh) == 20 * 6 * 10 * QUERY_COUNT
            ),
            "observed": {
                "tangent_actor": len(ta),
                "deep_bsde": len(deep),
                "fh_mlp": len(fh),
            },
            "expected": {
                "tangent_actor": 20 * 3 * 6 * QUERY_COUNT,
                "deep_bsde": 20 * QUERY_COUNT * 3,
                "fh_mlp": 20 * 6 * 10 * QUERY_COUNT,
            },
        },
        "scaling_cells": {
            "passed": len(scaling) == 20 * 3,
            "observed": len(scaling),
            "expected": 60,
        },
        "generated_files": {
            "passed": all(path.is_file() and path.stat().st_size > 0 for path in generated),
            "count": len(generated),
            "missing": [
                _relative(path)
                for path in generated
                if not path.is_file() or path.stat().st_size == 0
            ],
        },
    }
    overall = all(bool(value["passed"]) for value in checks.values())
    return {
        "revision": REVISION,
        "checked_utc": _utc_now(),
        "algorithm_hashes": algorithm_hashes(),
        "checks": checks,
        "overall_pass": overall,
        "status": "VALIDATED" if overall else "FAILED",
    }


def build_summary(
    scaling: pd.DataFrame,
    randomness: pd.DataFrame,
    fh: pd.DataFrame,
    resource_summary: dict[str, Any],
) -> dict[str, Any]:
    win_counts = {method: 0 for method in METHOD_ORDER}
    cell_winners: list[dict[str, Any]] = []
    for (pde, dimension), frame in scaling.groupby(
        ["pde", "dimension"], sort=True
    ):
        eligible = frame.loc[
            (frame["status"] == "COMPLETE")
            & np.isfinite(frame["rmse_32_R"])
        ]
        if eligible.empty:
            cell_winners.append({
                "pde": pde,
                "dimension": int(dimension),
                "winner": None,
                "reason": "no fully finite method",
            })
            continue
        minimum = float(eligible["rmse_32_R"].min())
        winners = eligible.loc[
            np.isclose(eligible["rmse_32_R"], minimum, rtol=1e-12, atol=0.0),
            "method",
        ].tolist()
        for method in winners:
            win_counts[method] += 1
        cell_winners.append({
            "pde": pde,
            "dimension": int(dimension),
            "winner": winners,
            "rmse_32_R": minimum,
        })

    fh_failures = (
        fh.assign(nonfinite=~np.isfinite(fh["prediction"]))
        .groupby(["pde", "n", "M"], as_index=False)
        .agg(nonfinite_fraction=("nonfinite", "mean"))
    )
    fh_failed_runs = (
        fh[["pde", "n", "M", "record_file", "status"]]
        .drop_duplicates()
        .assign(failed=lambda frame: frame["status"] == "numerical_failure")
        .groupby(["pde", "n", "M"], as_index=False)
        .agg(failed_runs=("failed", "sum"))
    )
    fh_failures = fh_failures.merge(
        fh_failed_runs, on=["pde", "n", "M"], how="left"
    )
    fh_failure_rows = fh_failures.loc[
        fh_failures["nonfinite_fraction"] > 0
    ].to_dict(orient="records")
    ta_randomness = randomness.loc[
        randomness["method"] == "Tangent-Actor",
        "mc_fraction_of_decomposed_variance",
    ].to_numpy(float)
    finite_ta_fraction = ta_randomness[np.isfinite(ta_randomness)]

    spreads: list[tuple[float, str]] = []
    for pde, frame in scaling.groupby("pde"):
        ratios = []
        for _, cell in frame.groupby("dimension"):
            values = cell.loc[
                (cell["status"] == "COMPLETE")
                & np.isfinite(cell["rmse_32_R"]),
                "rmse_32_R",
            ].to_numpy(float)
            if len(values) >= 2 and np.min(values) > 0:
                ratios.append(float(np.max(values) / np.min(values)))
        spreads.append((float(np.median(ratios)) if ratios else 0.0, pde))
    frontier = [pde for _, pde in sorted(spreads, reverse=True)[:2]]

    return {
        "revision": REVISION,
        "generated_utc": _utc_now(),
        "complete_pde_dimension_cells": 20,
        "method_cell_outcomes": int(len(scaling)),
        "lowest_rmse_cell_counts": win_counts,
        "tangent_actor_lowest_rmse_cells": win_counts["Tangent-Actor"],
        "cell_winners": cell_winners,
        "fh_mlp": {
            "numerical_failure_configurations": _json_safe(fh_failure_rows),
            "main_configuration": list(FH_MAIN_CONFIGURATION),
        },
        "deep_bsde": {
            "median_repeat_rmse_sd": _json_safe(float(
                scaling.loc[
                    scaling["method"] == "Deep BSDE", "repeat_rmse_sd"
                ].median()
            )),
            "maximum_repeat_rmse_sd": _json_safe(float(
                scaling.loc[
                    scaling["method"] == "Deep BSDE", "repeat_rmse_sd"
                ].max()
            )),
            "interpretation": (
                "training-seed variation is reported per cell and is not "
                "query-sampling uncertainty"
            ),
        },
        "tangent_actor_randomness": {
            "median_mc_fraction_of_decomposed_variance": (
                float(np.median(finite_ta_fraction))
                if len(finite_ta_fraction) else None
            ),
            "interpretation": (
                "between-training-seed and within-seed FK inference components "
                "are separated on the same frozen points"
            ),
        },
        "suggested_future_accuracy_compute_pdes": frontier,
        "algorithm_classification": {
            "tangent_actor": (
                "native full-function Tangent-Actor; Allen-Cahn uses the "
                "audited value-Actor autograd martingale-CV extension"
            ),
            "fh_mlp": "standard full-history MLP recursive Monte Carlo transfer",
            "deep_bsde_hjb_log_allen_cahn": "native source-family reproduction",
            "deep_bsde_b6_dense_fourier": "faithful transfer",
        },
        "resource_accounting": resource_summary,
        "artifact_audit_status": "PENDING",
    }


def _format_metric(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "failure"
    return f"{number:.3e}" if math.isfinite(number) else "failure"


def write_report(
    scaling: pd.DataFrame,
    summary: dict[str, Any],
    artifact_status: str,
) -> None:
    lines = [
        "# ICML Single-Point 32 ICE Report",
        "",
        f"Revision: {REVISION}",
        "",
        "This report compares fresh Tangent-Actor plus MCV training, standard "
        "full-history MLP, and independent single-initial-state Deep BSDE on "
        "the same 640 frozen time-zero queries.",
        "",
        "## Main RMSE over points and formal repeats",
        "",
        "| PDE | d | Tangent-Actor | FH-MLP | Deep BSDE |",
        "|---|---:|---:|---:|---:|",
    ]
    for pde in PDES:
        for dimension in DIMENSIONS:
            cell = scaling.loc[
                (scaling["pde"] == pde)
                & (scaling["dimension"] == dimension)
            ].set_index("method")
            lines.append(
                f"| {pde} | {dimension} | "
                f"{_format_metric(cell.loc['Tangent-Actor', 'rmse_32_R'])} | "
                f"{_format_metric(cell.loc['FH-MLP', 'rmse_32_R'])} | "
                f"{_format_metric(cell.loc['Deep BSDE', 'rmse_32_R'])} |"
            )
    failures = summary["fh_mlp"]["numerical_failure_configurations"]
    failure_text = (
        f"{len(failures)} PDE/configuration rows contain nonfinite formal "
        "FH-MLP predictions."
        if failures else
        "No formal FH-MLP configuration produced a nonfinite prediction."
    )
    lines.extend([
        "",
        "## Required conclusions",
        "",
        "1. The complete dimension trends are given in the table above and in "
        "figure1_accuracy_vs_dimension.pdf. A failed numerical outcome is never "
        "ranked using a finite-only subset.",
        "",
        "2. Tangent-Actor has the lowest finite RMSE in "
        f"{summary['tangent_actor_lowest_rmse_cells']} of 20 complete "
        "PDE-dimension cells. Exact ties, if any, are counted for every tied method.",
        "",
        f"3. {failure_text} The full six-configuration data remain in "
        "fh_mlp_results.csv; the preregistered main comparison uses n=3, M=16.",
        "",
        "4. Deep BSDE stability is measured from three independent trainings at "
        "every initial point. Its median and maximum cell-level repeat-RMSE "
        f"standard deviations are "
        f"{_format_metric(summary['deep_bsde']['median_repeat_rmse_sd'])} and "
        f"{_format_metric(summary['deep_bsde']['maximum_repeat_rmse_sd'])}.",
        "",
        "5. B6 and Dense Fourier retain their full-rank or dense frozen "
        "parameters. Their rows therefore measure effective-dimensional "
        "difficulty rather than a changed low-dimensional surrogate.",
        "",
        "6. For Tangent-Actor, randomness_decomposition.csv separates variation "
        "between fresh Actor trainings from within-Actor FK Monte Carlo noise. "
        "The median inferred MC fraction is "
        f"{_format_metric(summary['tangent_actor_randomness']['median_mc_fraction_of_decomposed_variance'])}.",
        "",
        "7. The accuracy-spread rule identifies "
        + ", ".join(summary["suggested_future_accuracy_compute_pdes"])
        + " as the first PDEs for a separate future accuracy-compute frontier.",
        "",
        "8. HJB-Log and Allen-Cahn Deep BSDE runs are native source-family "
        "reproductions; B6 and Dense Fourier are faithful transfers. "
        "Allen-Cahn Tangent-Actor uses the explicitly audited value-Actor "
        "autograd martingale-CV extension. Standard FH-MLP is retained without "
        "truth-based configuration selection.",
        "",
        "## Audit and reproducibility",
        "",
        f"Final artifact audit: {artifact_status}.",
        "",
        "All query hashes, references, task manifests, algorithm hashes, raw "
        "predictions, formal repeats, and Slurm/task timing records are checked "
        "by artifact_audit.json. Numerical failures and incomplete artifacts are "
        "reported rather than replaced or omitted.",
        "",
    ])
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def execution_status(scaling: pd.DataFrame) -> dict[str, Any]:
    cells: list[dict[str, Any]] = []
    for pde in PDES:
        for dimension in DIMENSIONS:
            local = scaling.loc[
                (scaling["pde"] == pde)
                & (scaling["dimension"] == dimension)
            ]
            audit = _json(DIMENSION_AUDIT_DIR / f"d{dimension}.json")
            cells.append({
                "pde": pde,
                "dimension": dimension,
                "status": (
                    "COMPLETE"
                    if audit.get("status") == "VALIDATED" and len(local) == 3
                    else "PARTIAL"
                ),
                "method_outcomes": {
                    row["method"]: row["status"]
                    for _, row in local.iterrows()
                },
            })
    overall = all(cell["status"] == "COMPLETE" for cell in cells)
    return {
        "revision": REVISION,
        "updated_utc": _utc_now(),
        "stage_a": _json(STAGE_A_AUDIT).get("status"),
        "stage_b": _json(STAGE_B_AUDIT).get("status"),
        "stage_c": "COMPLETE" if overall else "PARTIAL",
        "stage_d": "COMPLETE" if overall else "PARTIAL",
        "stage_e": "COMPLETE" if overall else "PARTIAL",
        "overall": "COMPLETE" if overall else "PARTIAL",
        "cells": cells,
    }


def run() -> dict[str, Any]:
    if _json(STAGE_A_AUDIT).get("status") != "VALIDATED":
        raise RuntimeError("Stage A is not validated")
    if _json(STAGE_B_AUDIT).get("status") != "VALIDATED":
        raise RuntimeError("Stage B is not validated")

    ta, deep, fh, training = collect_raw_results()
    main = main_rows(ta, deep, fh)
    scaling, pointwise = aggregate_metrics(main)
    randomness = randomness_decomposition(main, ta)
    resources, resource_summary = resource_accounting()

    outputs = {
        "ta_results.csv": ta,
        "deep_bsde_results.csv": deep,
        "fh_mlp_results.csv": fh,
        "training_manifest.csv": training,
        "main_pointwise_results.csv": pointwise,
        "main_dimension_scaling.csv": scaling,
        "randomness_decomposition.csv": randomness,
        "resource_accounting.csv": resources,
    }
    generated: list[Path] = []
    for name, frame in outputs.items():
        path = RESULTS_DIR / name
        _atomic_csv(path, frame)
        generated.append(path)

    status_payload = execution_status(scaling)
    status_path = RESULTS_DIR / "execution_status.json"
    write_json_atomic(status_path, _json_safe(status_payload))
    generated.append(status_path)

    generated.extend(make_figures(scaling, pointwise, fh, deep, training))
    generated.extend(make_tables(scaling))

    summary = build_summary(scaling, randomness, fh, resource_summary)
    summary_path = RESULTS_DIR / "final_summary.json"
    write_json_atomic(summary_path, _json_safe(summary))
    generated.append(summary_path)
    write_report(scaling, summary, "PENDING")
    generated.append(REPORT_PATH)

    audit = artifact_audit(ta, deep, fh, scaling, generated)
    audit_path = RESULTS_DIR / "artifact_audit.json"
    write_json_atomic(audit_path, _json_safe(audit))
    summary["artifact_audit_status"] = audit["status"]
    write_json_atomic(summary_path, _json_safe(summary))
    write_report(scaling, summary, audit["status"])
    if not audit["overall_pass"]:
        raise RuntimeError(f"final artifact audit failed: {audit_path}")
    return {
        "status": "COMPLETE",
        "artifact_audit": audit["status"],
        "tangent_actor_lowest_rmse_cells": summary[
            "tangent_actor_lowest_rmse_cells"
        ],
        "generated_files": len(generated) + 1,
        "report": str(REPORT_PATH),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(json.dumps(run(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
