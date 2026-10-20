"""Resumable two-update Tangent-Actor pilots and paired inference."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.models import (
    TerminalAwareActor,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    evaluate_policy as evaluate_gradient_policy,
    fit_actor as fit_gradient_actor,
    generate_labels as generate_gradient_labels,
)
from invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.solver import (
    evaluate_policy as evaluate_joint_policy,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.protocol import (
    stable_seed,
)
from invariant_region_mlp.experiments.tangent_actor_mcv_training.labels import (
    mcv_tangent_labels,
)

from .actors import JointTerminalActor
from .common import (
    environment_record,
    load_bank,
    load_dense,
    load_lqg,
    synchronize,
)
from .equations import (
    CoupledLQGEquation,
    DenseFourierDriver,
    LQGActorBenchmark,
    array_sha256,
)
from .prepare import file_sha256
from .protocol import (
    CACHE_ROOT,
    CHECKPOINT_ROOT,
    GRADIENT_CONFIG,
    INFERENCE_PATHS,
    INFERENCE_REPLICATES,
    INFERENCE_STEPS,
    JOINT_CONFIG,
    PREDICTION_ROOT,
    RAW_ROOT,
    REFERENCE_QUADRATURE_ORDER,
    REVISION,
    TRAINING_SEED,
    UPDATES,
    VALIDATION_PATH,
    ensure_directories,
    write_json_atomic,
)
from .training import fit_joint_actor


PDE_GPU_SECONDS_CAP = 3600.0


def _atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def _atomic_torch(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save(payload, temporary)
    temporary.replace(path)


def _checkpoint(pde: str, iteration: int) -> Path:
    return CHECKPOINT_ROOT / f"{pde}_seed{TRAINING_SEED}_a{iteration}.pt"


def _training_path(pde: str) -> Path:
    return RAW_ROOT / f"training_{pde}_seed{TRAINING_SEED}.json"


def _evaluation_path(pde: str, iteration: int, replicate: int) -> Path:
    return RAW_ROOT / "evaluations" / pde / f"a{iteration}_rep{replicate:02d}.json"


def _prediction_path(pde: str, iteration: int, replicate: int) -> Path:
    return PREDICTION_ROOT / pde / f"a{iteration}_rep{replicate:02d}.npz"


def _validation_gate() -> dict[str, Any]:
    if not VALIDATION_PATH.is_file():
        raise RuntimeError("run validate.py and freeze its artifacts before training")
    payload = json.loads(VALIDATION_PATH.read_text(encoding="utf-8"))
    if payload.get("status") != "VALIDATED" or not payload.get("overall_pass"):
        raise RuntimeError("REFERENCE_NOT_VALIDATED")
    return payload


def _save_dense_actor(
    actor: JointTerminalActor, iteration: int, metadata: dict[str, Any]
) -> Path:
    path = _checkpoint("dense_fourier", iteration)
    _atomic_torch(
        path,
        {
            "state_dict": {
                name: value.detach().cpu().clone()
                for name, value in actor.state_dict().items()
            },
            "pde": "dense_fourier",
            "iteration": iteration,
            "configuration": asdict(JOINT_CONFIG),
            "actor_specification": actor.specification(),
            "metadata": metadata,
        },
    )
    return path


def _load_dense_actor(
    driver: DenseFourierDriver, iteration: int, device: torch.device
) -> JointTerminalActor:
    payload = torch.load(
        _checkpoint("dense_fourier", iteration),
        map_location=device,
        weights_only=False,
    )
    actor = JointTerminalActor(
        driver, width=JOINT_CONFIG.width, depth=JOINT_CONFIG.depth
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor


def _save_lqg_actor(
    actor: TerminalAwareActor, iteration: int, metadata: dict[str, Any]
) -> Path:
    path = _checkpoint("coupled_lqg", iteration)
    _atomic_torch(
        path,
        {
            "state_dict": {
                name: value.detach().cpu().clone()
                for name, value in actor.state_dict().items()
            },
            "pde": "coupled_lqg",
            "iteration": iteration,
            "configuration": asdict(GRADIENT_CONFIG),
            "width": GRADIENT_CONFIG.width,
            "depth": GRADIENT_CONFIG.depth,
            "activation": GRADIENT_CONFIG.activation,
            "metadata": metadata,
        },
    )
    return path


def _load_lqg_actor(
    benchmark: LQGActorBenchmark, iteration: int, device: torch.device
) -> TerminalAwareActor:
    payload = torch.load(
        _checkpoint("coupled_lqg", iteration),
        map_location=device,
        weights_only=False,
    )
    actor = TerminalAwareActor(
        benchmark,
        width=int(payload["width"]),
        depth=int(payload["depth"]),
        activation=str(payload["activation"]),
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor


def _initial_dense(driver: DenseFourierDriver, device: torch.device) -> JointTerminalActor:
    seed_everything(stable_seed(REVISION, "dense_fourier", "a0", TRAINING_SEED))
    actor = JointTerminalActor(
        driver, width=JOINT_CONFIG.width, depth=JOINT_CONFIG.depth
    ).to(device)
    actor.eval()
    return actor


def _initial_lqg(
    benchmark: LQGActorBenchmark, device: torch.device
) -> TerminalAwareActor:
    seed_everything(stable_seed(REVISION, "coupled_lqg", "a0", TRAINING_SEED))
    actor = TerminalAwareActor(
        benchmark,
        width=GRADIENT_CONFIG.width,
        depth=GRADIENT_CONFIG.depth,
        activation=GRADIENT_CONFIG.activation,
    ).to(device)
    actor.eval()
    return actor


def _dense_labels(
    driver: DenseFourierDriver,
    actor: JointTerminalActor,
    iteration: int,
    bank: dict[str, np.ndarray],
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    cache = CACHE_ROOT / f"dense_fourier_a{iteration}_labels.npz"
    meta_path = CACHE_ROOT / f"dense_fourier_a{iteration}_labels.json"
    if cache.is_file() and meta_path.is_file():
        with np.load(cache) as data:
            return {key: data[key] for key in data.files}, json.loads(
                meta_path.read_text(encoding="utf-8")
            )
    split = f"train_a{iteration}"
    train_t, train_x = bank[f"{split}_t"], bank[f"{split}_x"]
    valid_t, valid_x = bank["validation_t"], bank["validation_x"]
    kwargs = {
        "paths": JOINT_CONFIG.paths_per_label,
        "time_steps": JOINT_CONFIG.label_time_steps,
        "device": device,
        "parent_batch": JOINT_CONFIG.parent_batch,
        "path_chunk": JOINT_CONFIG.path_chunk,
    }
    train_r, train_q, train_meta = mcv_tangent_labels(
        driver,
        actor,
        train_t,
        train_x,
        seed=stable_seed(REVISION, "dense", iteration, "train"),
        **kwargs,
    )
    valid_r, valid_q, valid_meta = mcv_tangent_labels(
        driver,
        actor,
        valid_t,
        valid_x,
        seed=stable_seed(REVISION, "dense", iteration, "validation"),
        **kwargs,
    )
    replicate_r, replicate_q, replicate_meta = mcv_tangent_labels(
        driver,
        actor,
        valid_t,
        valid_x,
        seed=stable_seed(REVISION, "dense", iteration, "replicate"),
        **kwargs,
    )
    if train_q is None or valid_q is None or replicate_q is None:
        raise RuntimeError("joint labels omitted q")
    arrays = {
        "train_t": train_t,
        "train_x": train_x,
        "train_r": train_r,
        "train_q": train_q,
        "valid_t": valid_t,
        "valid_x": valid_x,
        "valid_r": valid_r,
        "valid_q": valid_q,
    }
    _atomic_npz(cache, **arrays)
    metadata = {
        "iteration": iteration,
        "source": f"a{iteration - 1}",
        "train": train_meta,
        "validation": valid_meta,
        "replication": replicate_meta,
        "value_label_replication_rmse": float(
            np.sqrt(np.mean(np.square(valid_r - replicate_r)))
        ),
        "gradient_label_replication_component_rmse": float(
            np.sqrt(np.mean(np.square(valid_q - replicate_q)))
        ),
        "gradient_label_replication_vector_rmse": float(
            np.sqrt(np.mean(np.sum(np.square(valid_q - replicate_q), axis=1)))
        ),
        "train_bank_sha256": array_sha256(train_t, train_x),
        "validation_bank_sha256": array_sha256(valid_t, valid_x),
        "value_projection_used": False,
        "ground_truth_used": False,
    }
    write_json_atomic(meta_path, metadata)
    return arrays, metadata


def _lqg_labels(
    benchmark: LQGActorBenchmark,
    actor: TerminalAwareActor,
    iteration: int,
    bank: dict[str, np.ndarray],
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    cache = CACHE_ROOT / f"coupled_lqg_a{iteration}_labels.npz"
    meta_path = CACHE_ROOT / f"coupled_lqg_a{iteration}_labels.json"
    if cache.is_file() and meta_path.is_file():
        with np.load(cache) as data:
            return {key: data[key] for key in data.files}, json.loads(
                meta_path.read_text(encoding="utf-8")
            )
    split = f"train_a{iteration}"
    train_t, train_x = bank[f"{split}_t"], bank[f"{split}_x"]
    valid_t, valid_x = bank["validation_t"], bank["validation_x"]
    kwargs = {
        "paths": GRADIENT_CONFIG.paths_per_label,
        "time_steps": GRADIENT_CONFIG.label_time_steps,
        "device": device,
        "parent_batch": GRADIENT_CONFIG.parent_batch,
        "path_chunk": GRADIENT_CONFIG.path_chunk,
        "estimator": "mcv",
    }
    train_q, train_meta = generate_gradient_labels(
        benchmark,
        actor,
        train_t,
        train_x,
        seed=stable_seed(REVISION, "lqg", iteration, "train"),
        **kwargs,
    )
    valid_q, valid_meta = generate_gradient_labels(
        benchmark,
        actor,
        valid_t,
        valid_x,
        seed=stable_seed(REVISION, "lqg", iteration, "validation"),
        **kwargs,
    )
    replicate_q, replicate_meta = generate_gradient_labels(
        benchmark,
        actor,
        valid_t,
        valid_x,
        seed=stable_seed(REVISION, "lqg", iteration, "replicate"),
        **kwargs,
    )
    arrays = {
        "train_t": train_t,
        "train_x": train_x,
        "train_q": train_q,
        "valid_t": valid_t,
        "valid_x": valid_x,
        "valid_q": valid_q,
    }
    _atomic_npz(cache, **arrays)
    metadata = {
        "iteration": iteration,
        "source": f"a{iteration - 1}",
        "train": train_meta,
        "validation": valid_meta,
        "replication": replicate_meta,
        "gradient_label_replication_component_rmse": float(
            np.sqrt(np.mean(np.square(valid_q - replicate_q)))
        ),
        "gradient_label_replication_vector_rmse": float(
            np.sqrt(np.mean(np.sum(np.square(valid_q - replicate_q), axis=1)))
        ),
        "train_bank_sha256": array_sha256(train_t, train_x),
        "validation_bank_sha256": array_sha256(valid_t, valid_x),
        "ground_truth_used": False,
    }
    write_json_atomic(meta_path, metadata)
    return arrays, metadata


def _predict_dense(
    actor: JointTerminalActor,
    driver: DenseFourierDriver,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray]:
    values, gradients = [], []
    actor.eval()
    with torch.no_grad():
        for begin in range(0, len(t), 512):
            end = min(begin + 512, len(t))
            r, q = actor.predict(
                torch.as_tensor(t[begin:end], dtype=torch.float32, device=device),
                torch.as_tensor(x[begin:end], dtype=torch.float32, device=device),
                driver.T,
            )
            values.append(r.double().cpu().numpy())
            gradients.append(q.double().cpu().numpy())
    return np.concatenate(values), np.concatenate(gradients)


def _predict_lqg(
    actor: TerminalAwareActor,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> np.ndarray:
    gradients = []
    actor.eval()
    with torch.no_grad():
        for begin in range(0, len(t), 512):
            end = min(begin + 512, len(t))
            _r, q = actor.predict(
                torch.as_tensor(t[begin:end], dtype=torch.float32, device=device),
                torch.as_tensor(x[begin:end], dtype=torch.float32, device=device),
            )
            if q is None:
                raise RuntimeError("LQG actor omitted q")
            gradients.append(q.double().cpu().numpy())
    return np.concatenate(gradients)


def _dense_selector(
    driver: DenseFourierDriver,
    actors: list[JointTerminalActor],
    bank: dict[str, np.ndarray],
    device: torch.device,
) -> tuple[list[dict[str, Any]], float]:
    t, x = bank["selector_t"], bank["selector_x"]
    rows = []
    started = time.perf_counter()
    for iteration, actor in enumerate(actors):
        labels = []
        metadata = []
        for replicate in range(2):
            r, q, meta = mcv_tangent_labels(
                driver,
                actor,
                t,
                x,
                paths=128,
                time_steps=16,
                seed=stable_seed(REVISION, "dense_selector", iteration, replicate),
                device=device,
                parent_batch=4,
                path_chunk=32,
            )
            if q is None:
                raise RuntimeError("selector labels omitted q")
            labels.append((r, q))
            metadata.append(meta)
        actor_r, actor_q = _predict_dense(actor, driver, t, x, device)
        target_r = 0.5 * (labels[0][0] + labels[1][0])
        target_q = 0.5 * (labels[0][1] + labels[1][1])
        raw_r = float(np.mean(np.square(actor_r - target_r)))
        raw_q = float(np.mean(np.square(actor_q - target_q)))
        noise_r = float(np.mean(np.square(labels[0][0] - labels[1][0])) / 4.0)
        noise_q = float(np.mean(np.square(labels[0][1] - labels[1][1])) / 4.0)
        corrected_r = max(raw_r - noise_r, 0.0)
        corrected_q = max(raw_q - noise_q, 0.0)
        r_scale = max(float(np.std(target_r)), 1e-2)
        q_scale = max(float(np.sqrt(np.mean(np.square(target_q)))), 1e-3)
        rows.append(
            {
                "actor": f"a{iteration}",
                "value_raw_rmse": math.sqrt(raw_r),
                "value_noise_corrected_rmse": math.sqrt(corrected_r),
                "gradient_raw_component_rmse": math.sqrt(raw_q),
                "gradient_noise_corrected_component_rmse": math.sqrt(corrected_q),
                "normalized_value_inconsistency": corrected_r / r_scale**2,
                "normalized_gradient_inconsistency": corrected_q / q_scale**2,
                "normalized_total_inconsistency": (
                    corrected_r / r_scale**2 + corrected_q / q_scale**2
                ),
                "minimum_weight_ess_fraction": min(
                    meta["weight_ess_fraction"] for meta in metadata
                ),
                "nonfinite_count": int(
                    sum(meta["nonfinite_count"] for meta in metadata)
                ),
                "ground_truth_used": False,
            }
        )
    return rows, time.perf_counter() - started


def _lqg_selector(
    benchmark: LQGActorBenchmark,
    actors: list[TerminalAwareActor],
    bank: dict[str, np.ndarray],
    device: torch.device,
) -> tuple[list[dict[str, Any]], float]:
    t, x = bank["selector_t"], bank["selector_x"]
    rows = []
    started = time.perf_counter()
    for iteration, actor in enumerate(actors):
        labels, metadata = [], []
        for replicate in range(2):
            q, meta = generate_gradient_labels(
                benchmark,
                actor,
                t,
                x,
                paths=128,
                time_steps=16,
                seed=stable_seed(REVISION, "lqg_selector", iteration, replicate),
                device=device,
                parent_batch=4,
                path_chunk=32,
                estimator="mcv",
            )
            labels.append(q)
            metadata.append(meta)
        actor_q = _predict_lqg(actor, t, x, device)
        target = 0.5 * (labels[0] + labels[1])
        raw = float(np.mean(np.square(actor_q - target)))
        noise = float(np.mean(np.square(labels[0] - labels[1])) / 4.0)
        corrected = max(raw - noise, 0.0)
        scale = max(float(np.sqrt(np.mean(np.square(target)))), 1e-3)
        rows.append(
            {
                "actor": f"a{iteration}",
                "gradient_raw_component_rmse": math.sqrt(raw),
                "gradient_noise_corrected_component_rmse": math.sqrt(corrected),
                "normalized_total_inconsistency": corrected / scale**2,
                "minimum_weight_ess_fraction": min(
                    meta["weight_ess_fraction"] for meta in metadata
                ),
                "nonfinite_count": int(
                    sum(meta["nonfinite_count"] for meta in metadata)
                ),
                "ground_truth_used": False,
            }
        )
    return rows, time.perf_counter() - started


def train_dense(device: torch.device) -> dict[str, Any]:
    _validation_gate()
    ensure_directories()
    output = _training_path("dense_fourier")
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") == "complete":
            return existing
    driver, bank = load_dense(), load_bank("dense_fourier")
    actor = _initial_dense(driver, device)
    _save_dense_actor(actor, 0, {"kind": "untrained_terminal_tangent"})
    phases: list[dict[str, Any]] = []
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        phases = existing.get("phases", [])
        if phases:
            actor = _load_dense_actor(driver, len(phases), device)
    started = time.perf_counter()
    for iteration in range(len(phases) + 1, UPDATES + 1):
        if time.perf_counter() - started > PDE_GPU_SECONDS_CAP:
            raise TimeoutError("dense Fourier per-PDE GPU budget exhausted")
        phase_start = time.perf_counter()
        arrays, label_meta = _dense_labels(driver, actor, iteration, bank, device)
        fitted, fit_meta, trace = fit_joint_actor(
            driver,
            JOINT_CONFIG,
            actor,
            arrays["train_t"],
            arrays["train_x"],
            arrays["train_r"],
            arrays["train_q"],
            arrays["valid_t"],
            arrays["valid_x"],
            arrays["valid_r"],
            arrays["valid_q"],
            seed=stable_seed(REVISION, "dense_fit", iteration, TRAINING_SEED),
            device=device,
        )
        checkpoint = _save_dense_actor(fitted, iteration, fit_meta)
        phases.append(
            {
                "iteration": iteration,
                "checkpoint": str(checkpoint),
                "checkpoint_sha256": file_sha256(checkpoint),
                "phase_seconds": time.perf_counter() - phase_start,
                "labels": label_meta,
                "training": fit_meta,
                "trace": trace,
                "oom_chunk_fallback_used": False,
            }
        )
        actor = fitted
        write_json_atomic(
            output,
            {
                "study_revision": REVISION,
                "pde": "dense_fourier",
                "seed": TRAINING_SEED,
                "configuration": asdict(JOINT_CONFIG),
                "phases": phases,
                "status": "training_in_progress",
                "ground_truth_used_for_training_or_selection": False,
            },
        )
        print(
            f"dense a{iteration}: labels={label_meta['train']['wall_seconds']:.1f}s "
            f"fit={fit_meta['wall_seconds']:.1f}s "
            f"val={fit_meta['best_validation_normalized_loss']:.5g}",
            flush=True,
        )
    actors = [_load_dense_actor(driver, i, device) for i in range(UPDATES + 1)]
    selector_rows, selector_seconds = _dense_selector(driver, actors, bank, device)
    selected = min(
        (row for row in selector_rows if row["nonfinite_count"] == 0),
        key=lambda row: row["normalized_total_inconsistency"],
    )
    record = {
        "study_revision": REVISION,
        "pde": "dense_fourier",
        "seed": TRAINING_SEED,
        "configuration": asdict(JOINT_CONFIG),
        "driver": driver.metadata(),
        "phases": phases,
        "selection": {
            "criterion": "minimum independent noise-corrected fixed-point inconsistency",
            "rows": selector_rows,
            "selected_actor": selected["actor"],
            "ground_truth_used": False,
            "wall_seconds": selector_seconds,
        },
        "timing": {
            "label_generation_seconds": float(
                sum(
                    phase["labels"][split]["wall_seconds"]
                    for phase in phases
                    for split in ("train", "validation", "replication")
                )
            ),
            "fitting_seconds": float(
                sum(phase["training"]["wall_seconds"] for phase in phases)
            ),
            "selector_seconds": selector_seconds,
            "training_total_invocation_seconds": time.perf_counter() - started,
        },
        "ground_truth_used_for_training_or_selection": False,
        "status": "complete",
    }
    write_json_atomic(output, record)
    return record


def train_lqg(device: torch.device) -> dict[str, Any]:
    _validation_gate()
    ensure_directories()
    output = _training_path("coupled_lqg")
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") == "complete":
            return existing
    equation, bank = load_lqg(), load_bank("coupled_lqg")
    benchmark = LQGActorBenchmark(equation)
    actor = _initial_lqg(benchmark, device)
    _save_lqg_actor(actor, 0, {"kind": "untrained_terminal_tangent"})
    phases: list[dict[str, Any]] = []
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        phases = existing.get("phases", [])
        if phases:
            actor = _load_lqg_actor(benchmark, len(phases), device)
    started = time.perf_counter()
    for iteration in range(len(phases) + 1, UPDATES + 1):
        if time.perf_counter() - started > PDE_GPU_SECONDS_CAP:
            raise TimeoutError("coupled LQG per-PDE GPU budget exhausted")
        phase_start = time.perf_counter()
        arrays, label_meta = _lqg_labels(benchmark, actor, iteration, bank, device)
        fitted, fit_meta, trace = fit_gradient_actor(
            benchmark,
            GRADIENT_CONFIG,
            actor,
            arrays["train_t"],
            arrays["train_x"],
            arrays["train_q"],
            arrays["valid_t"],
            arrays["valid_x"],
            arrays["valid_q"],
            seed=stable_seed(REVISION, "lqg_fit", iteration, TRAINING_SEED),
            device=device,
        )
        checkpoint = _save_lqg_actor(fitted, iteration, fit_meta)
        phases.append(
            {
                "iteration": iteration,
                "checkpoint": str(checkpoint),
                "checkpoint_sha256": file_sha256(checkpoint),
                "phase_seconds": time.perf_counter() - phase_start,
                "labels": label_meta,
                "training": fit_meta,
                "trace": trace,
                "oom_chunk_fallback_used": False,
            }
        )
        actor = fitted
        write_json_atomic(
            output,
            {
                "study_revision": REVISION,
                "pde": "coupled_lqg",
                "seed": TRAINING_SEED,
                "configuration": asdict(GRADIENT_CONFIG),
                "phases": phases,
                "status": "training_in_progress",
                "ground_truth_used_for_training_or_selection": False,
            },
        )
        print(
            f"lqg a{iteration}: labels={label_meta['train']['wall_seconds']:.1f}s "
            f"fit={fit_meta['wall_seconds']:.1f}s "
            f"val={fit_meta['selected_validation_normalized_mse']:.5g}",
            flush=True,
        )
    actors = [_load_lqg_actor(benchmark, i, device) for i in range(UPDATES + 1)]
    selector_rows, selector_seconds = _lqg_selector(benchmark, actors, bank, device)
    selected = min(
        (row for row in selector_rows if row["nonfinite_count"] == 0),
        key=lambda row: row["normalized_total_inconsistency"],
    )
    record = {
        "study_revision": REVISION,
        "pde": "coupled_lqg",
        "benchmark_label": "RECREATED_PUBLISHED_BENCHMARK",
        "seed": TRAINING_SEED,
        "configuration": asdict(GRADIENT_CONFIG),
        "equation": equation.metadata(),
        "phases": phases,
        "selection": {
            "criterion": "minimum independent noise-corrected fixed-point inconsistency",
            "rows": selector_rows,
            "selected_actor": selected["actor"],
            "ground_truth_used": False,
            "wall_seconds": selector_seconds,
        },
        "timing": {
            "label_generation_seconds": float(
                sum(
                    phase["labels"][split]["wall_seconds"]
                    for phase in phases
                    for split in ("train", "validation", "replication")
                )
            ),
            "fitting_seconds": float(
                sum(phase["training"]["wall_seconds"] for phase in phases)
            ),
            "selector_seconds": selector_seconds,
            "training_total_invocation_seconds": time.perf_counter() - started,
        },
        "ground_truth_used_for_training_or_selection": False,
        "status": "complete",
    }
    write_json_atomic(output, record)
    return record


def _metrics(
    prediction: np.ndarray,
    truth: np.ndarray,
    pair_variance: np.ndarray,
    *,
    reference_uncertainty: float,
) -> dict[str, float]:
    error = prediction - truth
    rmse = float(np.sqrt(np.mean(np.square(error))))
    mc = float(
        np.sqrt(np.mean(pair_variance / (INFERENCE_PATHS // 2)))
    )
    return {
        "value_rmse": rmse,
        "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
        "bias": float(np.mean(error)),
        "mae": float(np.mean(np.abs(error))),
        "mc_rmse_contribution": mc,
        "non_mc_rmse_proxy": math.sqrt(max(rmse * rmse - mc * mc, 0.0)),
        "reference_uncertainty_rmse": reference_uncertainty,
    }


def evaluate_dense(device: torch.device) -> list[dict[str, Any]]:
    training = train_dense(device)
    driver, bank = load_dense(), load_bank("dense_fourier")
    t, x = bank["test_t"], bank["test_x"]
    truth_state = driver.exact_state(t, x)
    truth, truth_q = truth_state[:, 0], truth_state[:, 1:]
    rows = []
    for iteration in range(UPDATES + 1):
        actor = _load_dense_actor(driver, iteration, device)
        actor_r, actor_q = _predict_dense(actor, driver, t, x, device)
        direct_r_error = actor_r - truth
        direct_q_error = actor_q - truth_q
        for replicate in range(INFERENCE_REPLICATES):
            output_path = _evaluation_path("dense_fourier", iteration, replicate)
            if output_path.is_file():
                rows.append(json.loads(output_path.read_text(encoding="utf-8")))
                continue
            common_seed = stable_seed(REVISION, "dense_eval", replicate)
            output, meta = evaluate_joint_policy(
                driver,
                actor,
                t,
                x,
                paths=INFERENCE_PATHS,
                time_steps=INFERENCE_STEPS,
                seed=common_seed,
                device=device,
                parent_batch=16,
                path_chunk=128,
            )
            prediction = output["mcv"]
            prediction_path = _prediction_path("dense_fourier", iteration, replicate)
            _atomic_npz(
                prediction_path,
                t=t,
                x=x,
                truth=truth,
                truth_q=truth_q,
                actor_r=actor_r,
                actor_q=actor_q,
                plain=output["plain"],
                mcv=prediction,
                plain_pair_variance=output["plain_pair_variance"],
                mcv_pair_variance=output["mcv_pair_variance"],
            )
            row = {
                "study_revision": REVISION,
                "pde": "dense_fourier",
                "benchmark_label": "ANALYTICALLY_INTEGRABLE_HIGH_DIMENSIONAL_STRESS_TEST",
                "iteration": iteration,
                "actor": f"a{iteration}",
                "selected_by_reference_free_rule": (
                    training["selection"]["selected_actor"] == f"a{iteration}"
                ),
                "replicate": replicate,
                "paths": INFERENCE_PATHS,
                "steps": INFERENCE_STEPS,
                "seed": common_seed,
                **_metrics(
                    prediction,
                    truth,
                    output["mcv_pair_variance"],
                    reference_uncertainty=2e-13,
                ),
                "plain_value_rmse": float(
                    np.sqrt(np.mean(np.square(output["plain"] - truth)))
                ),
                "actor_value_rmse": float(
                    np.sqrt(np.mean(np.square(direct_r_error)))
                ),
                "actor_gradient_component_rmse": float(
                    np.sqrt(np.mean(np.square(direct_q_error)))
                ),
                "actor_gradient_vector_rmse": float(
                    np.sqrt(np.mean(np.sum(np.square(direct_q_error), axis=1)))
                ),
                "variance_reduction": meta["variance_reduction"],
                "weight_ess_fraction": meta["weight_ess_fraction"],
                "inference_seconds": meta["wall_seconds"],
                "peak_gpu_bytes": meta["peak_gpu_bytes"],
                "nonfinite_count": meta["nonfinite_count"],
                "test_bank_sha256": array_sha256(t, x),
                "prediction_path": str(prediction_path),
                "prediction_sha256": file_sha256(prediction_path),
                "ground_truth_loaded_only_after_selection": True,
            }
            write_json_atomic(output_path, row)
            rows.append(row)
            print(
                f"dense a{iteration} rep{replicate}: rmse={row['value_rmse']:.6g} "
                f"mc={row['mc_rmse_contribution']:.6g} time={row['inference_seconds']:.1f}s",
                flush=True,
            )
    return rows


def evaluate_lqg(device: torch.device) -> list[dict[str, Any]]:
    training = train_lqg(device)
    equation, bank = load_lqg(), load_bank("coupled_lqg")
    benchmark = LQGActorBenchmark(equation)
    t, x = bank["test_t"], bank["test_x"]
    coarse = equation.exact_state(t, x, order=192)
    truth_state = equation.exact_state(t, x, order=256)
    truth, truth_q = truth_state[:, 0], truth_state[:, 1:]
    reference_uncertainty = float(
        np.sqrt(np.mean(np.square(coarse[:, 0] - truth)))
    )
    rows = []
    for iteration in range(UPDATES + 1):
        actor = _load_lqg_actor(benchmark, iteration, device)
        actor_q = _predict_lqg(actor, t, x, device)
        direct_q_error = actor_q - truth_q
        for replicate in range(INFERENCE_REPLICATES):
            output_path = _evaluation_path("coupled_lqg", iteration, replicate)
            if output_path.is_file():
                rows.append(json.loads(output_path.read_text(encoding="utf-8")))
                continue
            common_seed = stable_seed(REVISION, "lqg_eval", replicate)
            output, meta = evaluate_gradient_policy(
                benchmark,
                actor,
                t,
                x,
                paths=INFERENCE_PATHS,
                time_steps=INFERENCE_STEPS,
                seed=common_seed,
                device=device,
                parent_batch=16,
                path_chunk=128,
            )
            prediction = output["mcv"]
            prediction_path = _prediction_path("coupled_lqg", iteration, replicate)
            _atomic_npz(
                prediction_path,
                t=t,
                x=x,
                truth=truth,
                truth_q=truth_q,
                actor_q=actor_q,
                plain=output["plain"],
                mcv=prediction,
                plain_pair_variance=output["plain_pair_variance"],
                mcv_pair_variance=output["mcv_pair_variance"],
            )
            row = {
                "study_revision": REVISION,
                "pde": "coupled_lqg",
                "benchmark_label": "RECREATED_PUBLISHED_BENCHMARK",
                "iteration": iteration,
                "actor": f"a{iteration}",
                "selected_by_reference_free_rule": (
                    training["selection"]["selected_actor"] == f"a{iteration}"
                ),
                "replicate": replicate,
                "paths": INFERENCE_PATHS,
                "steps": INFERENCE_STEPS,
                "seed": common_seed,
                **_metrics(
                    prediction,
                    truth,
                    output["mcv_pair_variance"],
                    reference_uncertainty=reference_uncertainty,
                ),
                "plain_value_rmse": float(
                    np.sqrt(np.mean(np.square(output["plain"] - truth)))
                ),
                "actor_value_rmse": None,
                "actor_gradient_component_rmse": float(
                    np.sqrt(np.mean(np.square(direct_q_error)))
                ),
                "actor_gradient_vector_rmse": float(
                    np.sqrt(np.mean(np.sum(np.square(direct_q_error), axis=1)))
                ),
                "variance_reduction": meta["variance_reduction"],
                "weight_ess_fraction": meta["weight_ess_fraction"],
                "inference_seconds": meta["wall_seconds"],
                "peak_gpu_bytes": meta["peak_gpu_bytes"],
                "nonfinite_count": meta["nonfinite_count"],
                "test_bank_sha256": array_sha256(t, x),
                "prediction_path": str(prediction_path),
                "prediction_sha256": file_sha256(prediction_path),
                "ground_truth_loaded_only_after_selection": True,
            }
            write_json_atomic(output_path, row)
            rows.append(row)
            print(
                f"lqg a{iteration} rep{replicate}: rmse={row['value_rmse']:.6g} "
                f"mc={row['mc_rmse_contribution']:.6g} time={row['inference_seconds']:.1f}s",
                flush=True,
            )
    return rows


def write_run_status() -> None:
    payload = {
        "study_revision": REVISION,
        "environment": environment_record(),
        "dense_training": str(_training_path("dense_fourier")),
        "lqg_training": str(_training_path("coupled_lqg")),
        "status": "pilot_execution_complete",
    }
    write_json_atomic(RAW_ROOT / "run_status.json", payload)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pde", choices=("dense_fourier", "coupled_lqg", "all"), default="all"
    )
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--train-only", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    args = parser.parse_args()
    if args.train_only and args.evaluate_only:
        raise ValueError("choose at most one of --train-only and --evaluate-only")
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    requested = (
        ("dense_fourier", "coupled_lqg") if args.pde == "all" else (args.pde,)
    )
    total_started = time.perf_counter()
    for pde in requested:
        if not args.evaluate_only:
            (train_dense if pde == "dense_fourier" else train_lqg)(device)
        if not args.train_only:
            (evaluate_dense if pde == "dense_fourier" else evaluate_lqg)(device)
        if time.perf_counter() - total_started > 10800.0:
            raise TimeoutError("total 3 GPU-hour pilot cap exhausted")
    write_run_status()


if __name__ == "__main__":
    main()


__all__ = [
    "evaluate_dense",
    "evaluate_lqg",
    "train_dense",
    "train_lqg",
]

