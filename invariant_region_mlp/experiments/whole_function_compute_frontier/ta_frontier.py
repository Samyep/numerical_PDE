"""Anytime Tangent-Actor trajectories using the mature frozen implementation."""

from __future__ import annotations

from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.main_highdim_expansion import run_b6_actor as b6
from invariant_region_mlp.experiments.main_highdim_expansion import run_dense_actor as dense
from invariant_region_mlp.experiments.main_highdim_expansion.protocol import (
    GRADIENT_ACTOR_CONFIG,
    JOINT_ACTOR_CONFIG,
    REVISION as HISTORICAL_REVISION,
    stable_seed as historical_seed,
)

from .common import (
    file_sha256,
    plateau_reached,
    save_npz_atomic,
    value_metrics,
    write_json_atomic,
)
from .prepare import historical_results_root
from .protocol import (
    CHECKPOINT_DIR,
    CONFIGURATION_GRID,
    DIMENSION,
    PREDICTION_DIR,
    RAW_DIR,
    REVISION,
    STATUS_DIR,
    TA_EVALUATION_PATHS,
    TA_EVALUATION_STEPS,
    TA_INFERENCE_REPLICATES,
    TRAINING_SEED,
    ensure_directories,
)


def _canonical_run_id(pde: str) -> str:
    return f"{pde}_d{DIMENSION}_mature_s{TRAINING_SEED}"


def _trajectory_id(pde: str, configuration: str) -> str:
    return f"ta_{pde}_d{DIMENSION}_{configuration.lower()}_s{TRAINING_SEED}"


def _record_path(pde: str, configuration: str) -> Path:
    return RAW_DIR / "trajectories" / f"{_trajectory_id(pde, configuration)}.json"


def _checkpoint_path(pde: str, configuration: str, iteration: int) -> Path:
    return (
        CHECKPOINT_DIR / "tangent_actor" / pde / configuration
        / f"a{iteration}.pt"
    )


def _cache_dir(pde: str, configuration: str) -> Path:
    return RAW_DIR.parent / "cache" / "tangent_actor" / pde / configuration


def _historical_training(pde: str) -> dict[str, Any]:
    path = (
        historical_results_root() / "raw" / "actor"
        / f"training_{_canonical_run_id(pde)}.json"
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "complete" or len(payload.get("phases", [])) != 4:
        raise RuntimeError(f"historical TA-C0 is incomplete: {path}")
    return payload


def _historical_checkpoint(pde: str, iteration: int) -> Path:
    return (
        historical_results_root() / "checkpoints"
        / f"{_canonical_run_id(pde)}_a{iteration}.pt"
    )


def _configs(pde: str, configuration: str):
    values = CONFIGURATION_GRID["tangent_actor"][configuration]
    if pde == "b6_fullrank":
        return replace(
            GRADIENT_ACTOR_CONFIG,
            paths_per_label=int(values["paths_per_label"]),
            optimizer_steps=int(values["optimizer_steps"]),
        )
    if pde == "dense_fourier":
        return replace(
            JOINT_ACTOR_CONFIG,
            paths_per_label=int(values["paths_per_label"]),
            optimizer_steps=int(values["optimizer_steps"]),
        )
    raise KeyError(pde)


def _initial_and_benchmark(pde: str, config: Any, device: torch.device):
    canonical = _canonical_run_id(pde)
    if pde == "b6_fullrank":
        benchmark = b6.make_benchmark(pde, DIMENSION)
        actor = b6._initial_actor(canonical, benchmark, config, device)
    elif pde == "dense_fourier":
        dense.JOINT_ACTOR_CONFIG = config
        benchmark = dense.make_dense_fourier(DIMENSION)
        actor = dense._initial_actor(canonical, benchmark, device)
    else:
        raise KeyError(pde)
    return benchmark, actor


def _load_actor(
    pde: str, path: Path, benchmark: Any, device: torch.device
):
    return (
        b6._load_actor(path, benchmark, device)
        if pde == "b6_fullrank"
        else dense._load_actor(path, benchmark, device)
    )


def _save_actor(
    pde: str,
    path: Path,
    actor: Any,
    benchmark: Any,
    config: Any,
    iteration: int,
    fit_meta: dict[str, Any],
) -> None:
    if pde == "b6_fullrank":
        b6._save_actor(path, actor, benchmark, config, iteration)
    else:
        dense._save_actor(path, actor, benchmark, iteration, fit_meta)


def _selector_row(
    pde: str,
    actor_name: str,
    actor: Any,
    benchmark: Any,
    device: torch.device,
) -> tuple[dict[str, Any], float]:
    started = time.perf_counter()
    if pde == "b6_fullrank":
        row = b6._fixed_point(
            _canonical_run_id(pde),
            pde,
            benchmark,
            actor_name,
            actor,
            TRAINING_SEED,
            device,
        )
    else:
        rows, _ = dense._selector(
            _canonical_run_id(pde), benchmark, [actor], TRAINING_SEED, device
        )
        row = dict(rows[0])
        row["actor"] = actor_name
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    return row, time.perf_counter() - started


def _new_phase(
    pde: str,
    configuration: str,
    iteration: int,
    source_actor: Any,
    benchmark: Any,
    config: Any,
    device: torch.device,
) -> tuple[Any, dict[str, Any]]:
    cache = _cache_dir(pde, configuration)
    cache.mkdir(parents=True, exist_ok=True)
    canonical = _canonical_run_id(pde)
    phase_started = time.perf_counter()
    if pde == "b6_fullrank":
        b6.CACHE_DIR = cache
        arrays, label_meta = b6._labels(
            canonical, pde, benchmark, source_actor, iteration,
            config, TRAINING_SEED, device,
        )
        fallback_used = False
        first_failure = None
        fit_config = config
        try:
            fitted, fit_meta, trace = b6.fit_actor(
                benchmark, fit_config, source_actor,
                arrays["train_t"], arrays["train_x"], arrays["train_y"],
                arrays["valid_t"], arrays["valid_x"], arrays["valid_y"],
                seed=historical_seed(canonical, iteration, "fit"), device=device,
            )
            if not b6._fit_is_finite(fitted, fit_meta, trace):
                raise FloatingPointError("non-finite fit result")
        except (FloatingPointError, RuntimeError) as error:
            first_failure = repr(error)
            fallback_used = True
            fit_config = replace(config, learning_rate=config.learning_rate / 2.0)
            fitted, fit_meta, trace = b6.fit_actor(
                benchmark, fit_config, source_actor,
                arrays["train_t"], arrays["train_x"], arrays["train_y"],
                arrays["valid_t"], arrays["valid_x"], arrays["valid_y"],
                seed=historical_seed(canonical, iteration, "fit_fallback"), device=device,
            )
            if not b6._fit_is_finite(fitted, fit_meta, trace):
                raise FloatingPointError("half-LR fallback is also non-finite")
    else:
        dense.CACHE_DIR = cache
        dense.JOINT_ACTOR_CONFIG = config
        arrays, label_meta = dense._labels(
            canonical, benchmark, source_actor, iteration, TRAINING_SEED, device
        )
        fitted, fit_meta, trace = dense.fit_joint_actor(
            benchmark, config, source_actor,
            arrays["train_t"], arrays["train_x"], arrays["train_r"], arrays["train_q"],
            arrays["valid_t"], arrays["valid_x"], arrays["valid_r"], arrays["valid_q"],
            seed=historical_seed(canonical, iteration, "fit"), device=device,
        )
        if not all(torch.isfinite(parameter).all().item() for parameter in fitted.parameters()):
            raise FloatingPointError("nonfinite Dense Fourier actor")
        fallback_used = False
        first_failure = None
        fit_config = config
    checkpoint = _checkpoint_path(pde, configuration, iteration)
    _save_actor(pde, checkpoint, fitted, benchmark, fit_config, iteration, fit_meta)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    phase = {
        "iteration": iteration,
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": file_sha256(checkpoint),
        "phase_seconds": time.perf_counter() - phase_started,
        "labels": label_meta,
        "training": fit_meta,
        "trace": trace,
        "fallback_used": fallback_used,
        "first_failure": first_failure,
        "effective_learning_rate": fit_config.learning_rate,
        "reused_historical": False,
    }
    return fitted, phase


def train(pde: str, configuration: str, device: torch.device) -> dict[str, Any]:
    if configuration not in ("C0", "C1", "C2"):
        raise KeyError(configuration)
    ensure_directories()
    output = _record_path(pde, configuration)
    existing: dict[str, Any] | None = None
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") in {"TRAINED_NO_TRUTH", "COMPLETE"}:
            return existing
    config = _configs(pde, configuration)
    benchmark, a0 = _initial_and_benchmark(pde, config, device)
    actors: list[Any] = [a0]
    phases: list[dict[str, Any]] = []
    checkpoint_paths: list[Path | None] = [None]
    selectors: list[dict[str, Any]] = []
    cumulative_selection = 0.0
    if existing is not None and existing.get("status") == "TRAINING_IN_PROGRESS":
        if existing.get("configuration") != asdict(config):
            raise RuntimeError("TA resume configuration drift")
        phases = list(existing.get("phases", []))
        selectors = list(existing.get("selector_rows", []))
        cumulative_selection = float(
            selectors[-1].get("cumulative_selection_seconds", 0.0)
        ) if selectors else 0.0
        for phase in phases:
            checkpoint = Path(phase["checkpoint"])
            if file_sha256(checkpoint) != phase["checkpoint_sha256"]:
                raise RuntimeError(f"TA resume checkpoint drift: {checkpoint}")
            checkpoint_paths.append(checkpoint)
            actors.append(_load_actor(pde, checkpoint, benchmark, device))
    elif configuration == "C0":
        historical = _historical_training(pde)
        if historical["configuration"] != asdict(config):
            raise RuntimeError("TA-C0 historical configuration drift")
        for phase in historical["phases"]:
            iteration = int(phase["iteration"])
            checkpoint = _historical_checkpoint(pde, iteration)
            if not checkpoint.is_file() or file_sha256(checkpoint) != phase["checkpoint_sha256"]:
                raise RuntimeError(f"historical checkpoint drift: {checkpoint}")
            copied = dict(phase)
            copied["checkpoint"] = str(checkpoint)
            copied["reused_historical"] = True
            phases.append(copied)
            checkpoint_paths.append(checkpoint)
            actors.append(_load_actor(pde, checkpoint, benchmark, device))

    if not selectors:
        for iteration, actor in enumerate(actors):
            row, seconds = _selector_row(pde, f"a{iteration}", actor, benchmark, device)
            cumulative_selection += seconds
            row.update({
                "selection_seconds": seconds,
                "cumulative_selection_seconds": cumulative_selection,
            })
            selectors.append(row)

    values = CONFIGURATION_GRID["tangent_actor"][configuration]
    requested_updates = int(values["actor_updates"])
    maximum_updates = int(values["continuation_actor_updates_max"])
    if configuration != "C0":
        requested_updates = maximum_updates
    stop_reason = "configured_actor_updates"
    while len(phases) < maximum_updates:
        completed = len(phases)
        training_seconds = float(sum(float(row["phase_seconds"]) for row in phases))
        if completed >= requested_updates:
            metrics = [float(row["normalized_total_inconsistency"]) for row in selectors[1:]]
            if plateau_reached(training_seconds, metrics):
                stop_reason = "predeclared_validation_plateau"
                break
            if training_seconds >= 30 * 60:
                stop_reason = "active_training_cap"
                break
        iteration = completed + 1
        fitted, phase = _new_phase(
            pde, configuration, iteration, actors[-1], benchmark, config, device
        )
        phases.append(phase)
        actors.append(fitted)
        checkpoint_paths.append(Path(phase["checkpoint"]))
        selector, seconds = _selector_row(
            pde, f"a{iteration}", fitted, benchmark, device
        )
        cumulative_selection += seconds
        selector.update({
            "selection_seconds": seconds,
            "cumulative_selection_seconds": cumulative_selection,
        })
        selectors.append(selector)
        write_json_atomic(output, {
            "revision": REVISION,
            "trajectory_id": _trajectory_id(pde, configuration),
            "method": "tangent_actor",
            "pde": pde,
            "dimension": DIMENSION,
            "configuration_id": configuration,
            "configuration": asdict(config),
            "phases": phases,
            "selector_rows": selectors,
            "checkpoint_paths": [str(path) if path else None for path in checkpoint_paths],
            "status": "TRAINING_IN_PROGRESS",
            "truth_loaded": False,
        })
    training_seconds = float(sum(float(row["phase_seconds"]) for row in phases))
    metrics = [float(row["normalized_total_inconsistency"]) for row in selectors[1:]]
    if len(phases) >= maximum_updates and configuration == "C0" and not plateau_reached(training_seconds, metrics):
        stop_reason = "maximum_predeclared_actor_updates"
    best_so_far = "a0"
    best_metric = float(selectors[0]["normalized_total_inconsistency"])
    cumulative_training = 0.0
    checkpoints: list[dict[str, Any]] = []
    for iteration, selector in enumerate(selectors):
        if iteration:
            cumulative_training += float(phases[iteration - 1]["phase_seconds"])
        metric = float(selector["normalized_total_inconsistency"])
        if selector.get("nonfinite_count", 0) == 0 and metric < best_metric:
            best_metric, best_so_far = metric, f"a{iteration}"
        checkpoint = checkpoint_paths[iteration]
        checkpoints.append({
            "checkpoint_id": f"a{iteration}",
            "algorithm_iteration": iteration,
            "checkpoint": str(checkpoint) if checkpoint else None,
            "checkpoint_sha256": file_sha256(checkpoint) if checkpoint else None,
            "training_seconds": cumulative_training,
            "selection_seconds": float(selector["cumulative_selection_seconds"]),
            "main_axis_seconds": cumulative_training + float(selector["cumulative_selection_seconds"]),
            "validation_metric": metric,
            "causal_selected_checkpoint": best_so_far,
            "is_original_budget": iteration == 4,
            "reused_historical": bool(iteration and phases[iteration - 1]["reused_historical"]),
        })
    result = {
        "revision": REVISION,
        "trajectory_id": _trajectory_id(pde, configuration),
        "method": "tangent_actor",
        "pde": pde,
        "dimension": DIMENSION,
        "seed": TRAINING_SEED,
        "configuration_id": configuration,
        "configuration": asdict(config),
        "phases": phases,
        "selector_rows": selectors,
        "checkpoints": checkpoints,
        "training_seconds": training_seconds,
        "selection_seconds": cumulative_selection,
        "stop_reason": stop_reason,
        "status": "TRAINED_NO_TRUTH",
        "ground_truth_used_for_training_selection_or_stopping": False,
    }
    write_json_atomic(output, result)
    return result


def _bank(pde: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Path]:
    path = historical_results_root() / "inputs" / f"{pde}_d{DIMENSION}_test_bank.npz"
    with np.load(path, allow_pickle=False) as payload:
        return (
            np.asarray(payload["t"]), np.asarray(payload["x"]),
            np.asarray(payload["truth_u"]), np.asarray(payload["truth_z"]), path,
        )


def evaluate(pde: str, configuration: str, device: torch.device) -> list[dict[str, Any]]:
    gate = STATUS_DIR / "all_training_complete.json"
    if not gate.is_file():
        raise RuntimeError("global truth gate is closed: all 18 trajectories must stop first")
    record_path = _record_path(pde, configuration)
    trajectory = json.loads(record_path.read_text(encoding="utf-8"))
    if trajectory.get("status") not in {"TRAINED_NO_TRUTH", "COMPLETE"}:
        raise RuntimeError("trajectory must stop before test truth is opened")
    config = _configs(pde, configuration)
    benchmark, a0 = _initial_and_benchmark(pde, config, device)
    t, x, truth, truth_z, bank_path = _bank(pde)
    rows: list[dict[str, Any]] = []
    for checkpoint in trajectory["checkpoints"]:
        actor_name = checkpoint["checkpoint_id"]
        actor = a0 if actor_name == "a0" else _load_actor(
            pde, Path(checkpoint["checkpoint"]), benchmark, device
        )
        for replicate in range(TA_INFERENCE_REPLICATES):
            output_path = (
                PREDICTION_DIR / "tangent_actor" / pde / configuration
                / f"{actor_name}_rep{replicate:02d}.npz"
            )
            metric_path = output_path.with_suffix(".json")
            if metric_path.is_file():
                rows.append(json.loads(metric_path.read_text(encoding="utf-8")))
                continue
            if pde == "b6_fullrank":
                seed = historical_seed(
                    HISTORICAL_REVISION, pde, DIMENSION, TRAINING_SEED,
                    TA_EVALUATION_STEPS, replicate,
                )
            else:
                seed = historical_seed(
                    HISTORICAL_REVISION, "dense_inference", DIMENSION,
                    TRAINING_SEED, TA_EVALUATION_STEPS, replicate,
                )
            started = time.perf_counter()
            evaluator = b6.evaluate_policy if pde == "b6_fullrank" else dense.evaluate_policy
            output, meta = evaluator(
                benchmark, actor, t, x,
                paths=TA_EVALUATION_PATHS,
                time_steps=TA_EVALUATION_STEPS,
                seed=seed,
                device=device,
                parent_batch=8,
                path_chunk=64,
            )
            prediction = np.asarray(output["mcv"], dtype=np.float64)
            pair_variance = np.asarray(output["mcv_pair_variance"], dtype=np.float64)
            metrics = value_metrics(prediction, truth)
            save_npz_atomic(
                output_path,
                t=t, x=x, truth=truth, truth_z=truth_z,
                prediction=prediction, pair_variance=pair_variance,
            )
            row = {
                "revision": REVISION,
                "trajectory_id": trajectory["trajectory_id"],
                "method": "tangent_actor",
                "pde": pde,
                "dimension": DIMENSION,
                "configuration_id": configuration,
                "seed": TRAINING_SEED,
                "checkpoint_id": actor_name,
                "algorithm_iteration": checkpoint["algorithm_iteration"],
                "inference_replicate": replicate,
                "training_seconds": checkpoint["training_seconds"],
                "selection_seconds": checkpoint["selection_seconds"],
                "main_axis_seconds": checkpoint["main_axis_seconds"],
                "inference_seconds": float(meta.get("wall_seconds", time.perf_counter() - started)),
                "validation_metric": checkpoint["validation_metric"],
                "causal_selected_checkpoint": checkpoint["causal_selected_checkpoint"],
                "is_original_budget": checkpoint["is_original_budget"],
                "checkpoint": checkpoint["checkpoint"],
                "checkpoint_sha256": checkpoint["checkpoint_sha256"],
                "test_bank": str(bank_path),
                "test_bank_sha256": file_sha256(bank_path),
                "prediction_file": str(output_path),
                "prediction_sha256": file_sha256(output_path),
                "paths": TA_EVALUATION_PATHS,
                "euler_steps": TA_EVALUATION_STEPS,
                "mc_rmse_contribution": float(
                    np.sqrt(np.mean(pair_variance / (TA_EVALUATION_PATHS // 2)))
                ),
                "nonfinite_count": int(meta.get("nonfinite_count", 0)),
                "status": "COMPLETE",
                **metrics,
            }
            write_json_atomic(metric_path, row)
            rows.append(row)
    trajectory["status"] = "COMPLETE"
    trajectory["truth_opened_only_after_stop"] = True
    trajectory["evaluation_rows"] = len(rows)
    write_json_atomic(record_path, trajectory)
    return rows


def execute(pde: str, configuration: str, device: torch.device) -> list[dict[str, Any]]:
    train(pde, configuration, device)
    return evaluate(pde, configuration, device)


__all__ = ["evaluate", "execute", "train"]
