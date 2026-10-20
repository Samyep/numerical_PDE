"""Mature joint value-gradient Tangent-Actor runner for Dense Fourier."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import time
from typing import Any, Iterable

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything,
)
from invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.solver import (
    evaluate_policy,
)
from invariant_region_mlp.experiments.tangent_actor_mcv_training.labels import (
    mcv_tangent_labels,
)

from .benchmarks import (
    DenseFourierDriver,
    array_sha256,
    make_dense_fourier,
    sample_cube_states,
)
from .joint_actor import JointTerminalActor
from .joint_training import fit_joint_actor
from .prepare import file_sha256, load_test_bank
from .protocol import (
    ACTOR_UPDATES,
    AUDIT_EULER_STEPS,
    CACHE_DIR,
    CHECKPOINT_DIR,
    DIMENSIONS,
    EVALUATION_PATHS,
    INFERENCE_REPLICATES,
    JOINT_ACTOR_CONFIG,
    PREDICTION_DIR,
    PRIMARY_EULER_STEPS,
    RAW_DIR,
    REVISION,
    TRAINING_SEED,
    VALIDATION_PATH,
    ensure_directories,
    stable_seed,
    write_json_atomic,
)


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


def _run_id(d: int, seed: int) -> str:
    return f"dense_fourier_d{d}_mature_s{seed}"


def _checkpoint(run_id: str, iteration: int) -> Path:
    return CHECKPOINT_DIR / f"{run_id}_a{iteration}.pt"


def _training_path(run_id: str) -> Path:
    return RAW_DIR / "actor" / f"training_{run_id}.json"


def _evaluation_path(
    run_id: str, actor_name: str, steps: int, replicate: int
) -> Path:
    return (
        RAW_DIR
        / "actor_evaluations"
        / run_id
        / actor_name
        / f"steps{steps}_rep{replicate:02d}.json"
    )


def _validation_gate() -> None:
    if not VALIDATION_PATH.is_file():
        raise RuntimeError(f"missing mathematical validation: {VALIDATION_PATH}")
    payload = json.loads(VALIDATION_PATH.read_text(encoding="utf-8"))
    if payload.get("status") != "VALIDATED" or not payload.get("overall_pass"):
        raise RuntimeError("mathematical validation gate did not pass")


def _save_actor(
    path: Path,
    actor: JointTerminalActor,
    driver: DenseFourierDriver,
    iteration: int,
    metadata: dict[str, Any],
) -> None:
    _atomic_torch(
        path,
        {
            "state_dict": {
                key: value.detach().cpu().clone()
                for key, value in actor.state_dict().items()
            },
            "benchmark": driver.name,
            "dimension": driver.d,
            "iteration": iteration,
            "configuration": asdict(JOINT_ACTOR_CONFIG),
            "actor_specification": actor.specification(),
            "metadata": metadata,
            "oracle_projection": False,
        },
    )


def _load_actor(
    path: Path, driver: DenseFourierDriver, device: torch.device
) -> JointTerminalActor:
    payload = torch.load(path, map_location=device, weights_only=False)
    if payload["configuration"] != asdict(JOINT_ACTOR_CONFIG):
        raise RuntimeError("checkpoint configuration drift")
    actor = JointTerminalActor(
        driver, width=JOINT_ACTOR_CONFIG.width, depth=JOINT_ACTOR_CONFIG.depth
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor


def _initial_actor(
    run_id: str, driver: DenseFourierDriver, device: torch.device
) -> JointTerminalActor:
    seed_everything(stable_seed(run_id, "a0"))
    actor = JointTerminalActor(
        driver, width=JOINT_ACTOR_CONFIG.width, depth=JOINT_ACTOR_CONFIG.depth
    ).to(device)
    actor.eval()
    return actor


def _states(
    d: int, count: int, *, split: str, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    return sample_cube_states(
        "dense_fourier",
        d,
        count,
        seed=stable_seed(REVISION, "dense_states", d, split, seed),
    )


def _labels(
    run_id: str,
    driver: DenseFourierDriver,
    source: JointTerminalActor,
    iteration: int,
    seed: int,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    cache = CACHE_DIR / f"{run_id}_a{iteration}_labels.npz"
    metadata_path = CACHE_DIR / f"{run_id}_a{iteration}_labels.json"
    if cache.is_file() and metadata_path.is_file():
        with np.load(cache) as data:
            arrays = {key: data[key] for key in data.files}
        return arrays, json.loads(metadata_path.read_text(encoding="utf-8"))
    cfg = JOINT_ACTOR_CONFIG
    train_t, train_x = _states(
        driver.d, cfg.training_points, split=f"train_a{iteration}", seed=seed
    )
    valid_t, valid_x = _states(
        driver.d, cfg.validation_points, split="fixed_validation", seed=seed
    )
    kwargs = {
        "paths": cfg.paths_per_label,
        "time_steps": cfg.label_time_steps,
        "device": device,
        "parent_batch": cfg.parent_batch,
        "path_chunk": cfg.path_chunk,
    }
    train_r, train_q, train_meta = mcv_tangent_labels(
        driver,
        source,
        train_t,
        train_x,
        seed=stable_seed(run_id, iteration, "train_labels"),
        **kwargs,
    )
    valid_r, valid_q, valid_meta = mcv_tangent_labels(
        driver,
        source,
        valid_t,
        valid_x,
        seed=stable_seed(run_id, iteration, "validation_labels"),
        **kwargs,
    )
    replicate_r, replicate_q, replicate_meta = mcv_tangent_labels(
        driver,
        source,
        valid_t,
        valid_x,
        seed=stable_seed(run_id, iteration, "replicate_labels"),
        **kwargs,
    )
    if train_q is None or valid_q is None or replicate_q is None:
        raise RuntimeError("joint Dense-Fourier labels omitted the gradient")
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
    write_json_atomic(metadata_path, metadata)
    return arrays, metadata


def _predict(
    actor: JointTerminalActor,
    driver: DenseFourierDriver,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray]:
    values: list[np.ndarray] = []
    gradients: list[np.ndarray] = []
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


def _selector(
    run_id: str,
    driver: DenseFourierDriver,
    actors: list[JointTerminalActor],
    seed: int,
    device: torch.device,
) -> tuple[list[dict[str, Any]], float]:
    t, x = _states(driver.d, 128, split="fixed_point_selector", seed=seed)
    rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    for iteration, actor in enumerate(actors):
        labels: list[tuple[np.ndarray, np.ndarray]] = []
        metadata: list[dict[str, Any]] = []
        for replicate in range(2):
            r, q, meta = mcv_tangent_labels(
                driver,
                actor,
                t,
                x,
                paths=128,
                time_steps=16,
                seed=stable_seed(run_id, "selector", iteration, replicate),
                device=device,
                parent_batch=4,
                path_chunk=32,
            )
            if q is None:
                raise RuntimeError("selector labels omitted q")
            labels.append((r, q))
            metadata.append(meta)
        actor_r, actor_q = _predict(actor, driver, t, x, device)
        target_r = 0.5 * (labels[0][0] + labels[1][0])
        target_q = 0.5 * (labels[0][1] + labels[1][1])
        raw_r = float(np.mean(np.square(actor_r - target_r)))
        raw_q = float(np.mean(np.square(actor_q - target_q)))
        noise_r = float(np.mean(np.square(labels[0][0] - labels[1][0])) / 4.0)
        noise_q = float(np.mean(np.square(labels[0][1] - labels[1][1])) / 4.0)
        corrected_r, corrected_q = max(raw_r - noise_r, 0.0), max(raw_q - noise_q, 0.0)
        r_scale = max(float(np.std(target_r)), 1e-2)
        q_scale = max(float(np.sqrt(np.mean(np.square(target_q)))), 1e-3)
        rows.append(
            {
                "actor": f"a{iteration}",
                "state_bank_sha256": array_sha256(t, x),
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
                    item["weight_ess_fraction"] for item in metadata
                ),
                "nonfinite_count": int(sum(item["nonfinite_count"] for item in metadata)),
                "ground_truth_used": False,
            }
        )
    return rows, time.perf_counter() - started


def train(d: int, seed: int, device: torch.device) -> dict[str, Any]:
    _validation_gate()
    ensure_directories()
    if seed != TRAINING_SEED:
        raise ValueError("formal expansion permits only seed 0")
    driver = make_dense_fourier(d)
    run_id = _run_id(d, seed)
    output = _training_path(run_id)
    actor = _initial_actor(run_id, driver, device)
    a0_path = _checkpoint(run_id, 0)
    if not a0_path.is_file():
        _save_actor(a0_path, actor, driver, 0, {"kind": "zero-training terminal actor"})
    phases: list[dict[str, Any]] = []
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("configuration") != asdict(JOINT_ACTOR_CONFIG):
            raise RuntimeError("existing training artifact has a different configuration")
        phases = existing.get("phases", [])
        if phases:
            actor = _load_actor(_checkpoint(run_id, len(phases)), driver, device)
        if existing.get("status") == "complete":
            return existing
    invocation_started = time.perf_counter()
    for iteration in range(len(phases) + 1, ACTOR_UPDATES + 1):
        phase_started = time.perf_counter()
        arrays, label_meta = _labels(
            run_id, driver, actor, iteration, seed, device
        )
        fitted, fit_meta, trace = fit_joint_actor(
            driver,
            JOINT_ACTOR_CONFIG,
            actor,
            arrays["train_t"],
            arrays["train_x"],
            arrays["train_r"],
            arrays["train_q"],
            arrays["valid_t"],
            arrays["valid_x"],
            arrays["valid_r"],
            arrays["valid_q"],
            seed=stable_seed(run_id, iteration, "fit"),
            device=device,
        )
        if not all(torch.all(torch.isfinite(p)).item() for p in fitted.parameters()):
            raise FloatingPointError("nonfinite fitted Dense-Fourier actor")
        checkpoint = _checkpoint(run_id, iteration)
        _save_actor(checkpoint, fitted, driver, iteration, fit_meta)
        phases.append(
            {
                "iteration": iteration,
                "checkpoint": str(checkpoint),
                "checkpoint_sha256": file_sha256(checkpoint),
                "phase_seconds": time.perf_counter() - phase_started,
                "labels": label_meta,
                "training": fit_meta,
                "trace": trace,
            }
        )
        actor = fitted
        write_json_atomic(
            output,
            {
                "study_revision": REVISION,
                "run_id": run_id,
                "method": "terminal_aware_joint_tangent_actor_mcv",
                "pde": "dense_fourier",
                "dimension": d,
                "seed": seed,
                "configuration": asdict(JOINT_ACTOR_CONFIG),
                "phases": phases,
                "status": "training_in_progress",
                "ground_truth_used_for_training_or_selection": False,
            },
        )
        print(
            f"{run_id} a{iteration} labels={label_meta['train']['wall_seconds']:.1f}s "
            f"fit={fit_meta['wall_seconds']:.1f}s "
            f"validation={fit_meta['best_validation_normalized_loss']:.6g}",
            flush=True,
        )
    actors = [
        _load_actor(_checkpoint(run_id, iteration), driver, device)
        for iteration in range(ACTOR_UPDATES + 1)
    ]
    selector_rows, selector_seconds = _selector(run_id, driver, actors, seed, device)
    eligible = [row for row in selector_rows if row["nonfinite_count"] == 0]
    if not eligible:
        raise FloatingPointError("all Dense-Fourier actors failed the selector")
    selected = min(eligible, key=lambda row: row["normalized_total_inconsistency"])
    result = {
        "study_revision": REVISION,
        "run_id": run_id,
        "method": "terminal_aware_joint_tangent_actor_mcv",
        "pde": "dense_fourier",
        "dimension": d,
        "seed": seed,
        "configuration": asdict(JOINT_ACTOR_CONFIG),
        "driver": driver.metadata(),
        "phases": phases,
        "selection": {
            "criterion": "minimum independent noise-corrected joint fixed-point inconsistency",
            "rows": selector_rows,
            "selected_actor": selected["actor"],
            "selected_checkpoint": str(
                _checkpoint(run_id, int(selected["actor"][1:]))
            ),
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
            "fitting_seconds": float(sum(phase["training"]["wall_seconds"] for phase in phases)),
            "phase_total_seconds": float(sum(phase["phase_seconds"] for phase in phases)),
            "reference_free_selection_seconds": selector_seconds,
            "current_invocation_seconds": time.perf_counter() - invocation_started,
        },
        "full_ambient_input": True,
        "oracle_projection": False,
        "parameter_count": sum(parameter.numel() for parameter in actors[0].parameters()),
        "ground_truth_used_for_training_or_selection": False,
        "status": "complete",
    }
    write_json_atomic(output, result)
    return result


def _actors_for_evaluation(
    run_id: str,
    driver: DenseFourierDriver,
    training: dict[str, Any],
    names: Iterable[str],
    device: torch.device,
) -> dict[str, JointTerminalActor]:
    result: dict[str, JointTerminalActor] = {}
    for name in names:
        resolved = training["selection"]["selected_actor"] if name == "selected" else name
        result[name] = _load_actor(_checkpoint(run_id, int(resolved[1:])), driver, device)
    return result


def evaluate(
    d: int,
    seed: int,
    device: torch.device,
    *,
    actor_names: tuple[str, ...] = ("selected", "a0"),
    steps: tuple[int, ...] = (PRIMARY_EULER_STEPS,),
) -> list[dict[str, Any]]:
    run_id = _run_id(d, seed)
    training = json.loads(_training_path(run_id).read_text(encoding="utf-8"))
    if training.get("status") != "complete":
        raise RuntimeError("training and reference-free selection are incomplete")
    driver = make_dense_fourier(d)
    actors = _actors_for_evaluation(run_id, driver, training, actor_names, device)
    t, x, truth, truth_q, bank_path = load_test_bank("dense_fourier", d)
    rows: list[dict[str, Any]] = []
    for actor_name, actor in actors.items():
        actor_r, actor_q = _predict(actor, driver, t, x, device)
        direct_r_error, direct_q_error = actor_r - truth, actor_q - truth_q
        for time_steps in steps:
            for replicate in range(INFERENCE_REPLICATES):
                output_path = _evaluation_path(run_id, actor_name, time_steps, replicate)
                if output_path.is_file():
                    rows.append(json.loads(output_path.read_text(encoding="utf-8")))
                    continue
                output, meta = evaluate_policy(
                    driver,
                    actor,
                    t,
                    x,
                    paths=EVALUATION_PATHS,
                    time_steps=time_steps,
                    seed=stable_seed(REVISION, "dense_inference", d, seed, time_steps, replicate),
                    device=device,
                    parent_batch=8,
                    path_chunk=64,
                )
                prediction = output["mcv"]
                pair_variance = output["mcv_pair_variance"]
                error = prediction - truth
                rmse = float(np.sqrt(np.mean(np.square(error))))
                mc = float(np.sqrt(np.mean(pair_variance / (EVALUATION_PATHS // 2))))
                prediction_path = (
                    PREDICTION_DIR
                    / "actor"
                    / run_id
                    / actor_name
                    / f"steps{time_steps}_rep{replicate:02d}.npz"
                )
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
                    mcv_pair_variance=pair_variance,
                )
                row = {
                    "study_revision": REVISION,
                    "method": "TA" if actor_name == "selected" else f"TA_{actor_name}",
                    "actor_request": actor_name,
                    "resolved_actor": training["selection"]["selected_actor"] if actor_name == "selected" else actor_name,
                    "pde": "dense_fourier",
                    "dimension": d,
                    "training_seed": seed,
                    "inference_replicate": replicate,
                    "steps": time_steps,
                    "paths": EVALUATION_PATHS,
                    "value_rmse": rmse,
                    "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
                    "bias": float(error.mean()),
                    "mae": float(np.abs(error).mean()),
                    "actor_value_rmse": float(np.sqrt(np.mean(np.square(direct_r_error)))),
                    "actor_gradient_component_rmse": float(np.sqrt(np.mean(np.square(direct_q_error)))),
                    "actor_gradient_vector_rmse": float(
                        np.sqrt(np.mean(np.sum(np.square(direct_q_error), axis=1)))
                    ),
                    "gradient_relative_l2": float(np.linalg.norm(direct_q_error) / np.linalg.norm(truth_q)),
                    "plain_value_rmse": float(np.sqrt(np.mean(np.square(output["plain"] - truth)))),
                    "mc_rmse_contribution": mc,
                    "non_mc_rmse_proxy": math.sqrt(max(rmse * rmse - mc * mc, 0.0)),
                    "variance_reduction": meta["variance_reduction"],
                    "weight_ess_fraction": meta["weight_ess_fraction"],
                    "inference_seconds": meta["wall_seconds"],
                    "peak_gpu_bytes": meta["peak_gpu_bytes"],
                    "nonfinite_count": meta["nonfinite_count"],
                    "test_states": len(t),
                    "test_bank_sha256": array_sha256(t, x),
                    "test_bank_file": str(bank_path),
                    "test_bank_file_sha256": file_sha256(bank_path),
                    "prediction_path": str(prediction_path),
                    "prediction_sha256": file_sha256(prediction_path),
                    "ground_truth_loaded_only_after_checkpoint_selection": True,
                }
                write_json_atomic(output_path, row)
                rows.append(row)
                print(
                    f"{run_id} {actor_name} steps={time_steps} rep={replicate} "
                    f"rmse={rmse:.6g} mc={mc:.3g} time={meta['wall_seconds']:.1f}s",
                    flush=True,
                )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, choices=DIMENSIONS, required=True)
    parser.add_argument("--seed", type=int, default=TRAINING_SEED)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--train-only", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--actors", default="selected,a0")
    parser.add_argument("--steps", default=str(PRIMARY_EULER_STEPS))
    args = parser.parse_args()
    if args.train_only and args.evaluate_only:
        raise ValueError("choose at most one of --train-only and --evaluate-only")
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation required")
    if not args.evaluate_only:
        result = train(args.dimension, args.seed, device)
        print(json.dumps({"training": result["status"], "selected": result["selection"]["selected_actor"]}), flush=True)
    if not args.train_only:
        evaluate(
            args.dimension,
            args.seed,
            device,
            actor_names=tuple(value for value in args.actors.split(",") if value),
            steps=tuple(int(value) for value in args.steps.split(",") if value),
        )


if __name__ == "__main__":
    main()


__all__ = ["evaluate", "train"]
