"""Mature frozen Tangent-Actor training and inference on S05/B6."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any, Iterable

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything,
    synchronize,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.models import (
    TerminalAwareActor,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    ActorConfig,
    evaluate_policy,
    fit_actor,
    generate_labels,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.protocol import (
    stable_seed,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_standard_actor import (
    MATURE,
)

from .equations import array_sha256, fixed_points, make_actor_benchmark
from .protocol import (
    CACHE_ROOT,
    CHECKPOINT_ROOT,
    PREDICTION_ROOT,
    RAW_ROOT,
    REVISION,
    TA_AUDIT_STEPS,
    TA_PATHS,
    TA_PRIMARY_STEPS,
    TA_REPLICATES,
    ensure_directories,
    write_json_atomic,
)
from .solver import load_b6_selection


def make_benchmark(pde_id: str, d: int):
    if pde_id == "B6":
        beta, gamma, _ = load_b6_selection(d)
        return make_actor_benchmark(pde_id, d, beta=beta, gamma=gamma)
    return make_actor_benchmark(pde_id, d)


def _run_id(pde_id: str, d: int, seed: int) -> str:
    return f"{pde_id}_d{d}_mature_s{seed}"


def _checkpoint(run_id: str, iteration: int) -> Path:
    return CHECKPOINT_ROOT / f"{run_id}_a{iteration}.pt"


def _training_path(run_id: str) -> Path:
    return RAW_ROOT / "actor" / f"training_{run_id}.json"


def _evaluation_path(
    run_id: str, actor_name: str, steps: int, replicate: int
) -> Path:
    return (
        RAW_ROOT
        / "actor_evaluations"
        / run_id
        / actor_name
        / f"steps{steps}_rep{replicate:02d}.json"
    )


def _save_actor(
    path: Path,
    actor: TerminalAwareActor,
    benchmark: Any,
    config: ActorConfig,
    iteration: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save(
        {
            "state_dict": {
                key: value.detach().cpu() for key, value in actor.state_dict().items()
            },
            "benchmark": benchmark.name,
            "dimension": benchmark.d,
            "iteration": iteration,
            "width": config.width,
            "depth": config.depth,
            "activation": config.activation,
            "config": asdict(config),
            "oracle_projection": False,
        },
        temporary,
    )
    temporary.replace(path)


def _load_actor(
    path: Path, benchmark: Any, device: torch.device
) -> TerminalAwareActor:
    payload = torch.load(path, map_location=device, weights_only=False)
    actor = TerminalAwareActor(
        benchmark,
        width=int(payload["width"]),
        depth=int(payload["depth"]),
        activation=str(payload["activation"]),
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor


def _initial_actor(
    run_id: str, benchmark: Any, config: ActorConfig, device: torch.device
) -> TerminalAwareActor:
    seed_everything(stable_seed(run_id, "a0"))
    actor = TerminalAwareActor(
        benchmark,
        width=config.width,
        depth=config.depth,
        activation=config.activation,
    ).to(device)
    actor.eval()
    return actor


def _states(
    pde_id: str,
    benchmark: Any,
    count: int,
    *,
    split: str,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    return benchmark.sample_states(
        count,
        seed=stable_seed("highdim_nosym_states", pde_id, benchmark.d, split, seed),
        purpose="function",
    )


def _labels(
    run_id: str,
    pde_id: str,
    benchmark: Any,
    source: TerminalAwareActor,
    iteration: int,
    config: ActorConfig,
    seed: int,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    cache = CACHE_ROOT / f"{run_id}_a{iteration}_labels.npz"
    metadata_path = CACHE_ROOT / f"{run_id}_a{iteration}_labels.json"
    if cache.is_file() and metadata_path.is_file():
        with np.load(cache) as data:
            arrays = {key: data[key] for key in data.files}
        return arrays, json.loads(metadata_path.read_text(encoding="utf-8"))
    train_t, train_x = _states(
        pde_id,
        benchmark,
        config.training_points,
        split=f"train_a{iteration}",
        seed=seed,
    )
    valid_t, valid_x = _states(
        pde_id,
        benchmark,
        config.validation_points,
        split="fixed_validation",
        seed=seed,
    )
    kwargs = {
        "paths": config.paths_per_label,
        "time_steps": config.label_time_steps,
        "device": device,
        "parent_batch": config.parent_batch,
        "path_chunk": config.path_chunk,
        "estimator": config.label_estimator,
    }
    train_y, train_meta = generate_labels(
        benchmark,
        source,
        train_t,
        train_x,
        seed=stable_seed(run_id, iteration, "train_labels"),
        **kwargs,
    )
    valid_y, valid_meta = generate_labels(
        benchmark,
        source,
        valid_t,
        valid_x,
        seed=stable_seed(run_id, iteration, "validation_labels"),
        **kwargs,
    )
    replicate_y, replicate_meta = generate_labels(
        benchmark,
        source,
        valid_t,
        valid_x,
        seed=stable_seed(run_id, iteration, "replicate_labels"),
        **kwargs,
    )
    arrays = {
        "train_t": train_t,
        "train_x": train_x,
        "train_y": train_y,
        "valid_t": valid_t,
        "valid_x": valid_x,
        "valid_y": valid_y,
    }
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, **arrays)
    difference = valid_y - replicate_y
    metadata = {
        "iteration": iteration,
        "source": f"a{iteration - 1}",
        "train": train_meta,
        "validation": valid_meta,
        "replication": replicate_meta,
        "label_replication_component_rmse": float(
            np.sqrt(np.mean(difference**2))
        ),
        "label_replication_vector_rmse": float(
            np.sqrt(np.mean(np.sum(difference**2, axis=1)))
        ),
        "train_bank_sha256": array_sha256(train_t, train_x),
        "validation_bank_sha256": array_sha256(valid_t, valid_x),
        "ground_truth_used": False,
    }
    write_json_atomic(metadata_path, metadata)
    return arrays, metadata


def _predict_q(
    actor: TerminalAwareActor,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> np.ndarray:
    parts = []
    actor.eval()
    with torch.no_grad():
        for begin in range(0, len(t), 1024):
            end = min(begin + 1024, len(t))
            _r, q = actor.predict(
                torch.as_tensor(t[begin:end], dtype=torch.float32, device=device),
                torch.as_tensor(x[begin:end], dtype=torch.float32, device=device),
            )
            if q is None:
                raise RuntimeError("gradient actor returned no q")
            parts.append(q.double().cpu().numpy())
    return np.concatenate(parts)


def _fixed_point(
    run_id: str,
    pde_id: str,
    benchmark: Any,
    actor_name: str,
    actor: TerminalAwareActor,
    seed: int,
    device: torch.device,
) -> dict[str, Any]:
    count, paths, steps = 128, 256, 32
    t, x = _states(
        pde_id, benchmark, count, split="fixed_point", seed=seed
    )
    labels = []
    metadata = []
    for replication in (0, 1):
        value, meta = generate_labels(
            benchmark,
            actor,
            t,
            x,
            paths=paths,
            time_steps=steps,
            seed=stable_seed(run_id, "fixed_point", actor_name, replication),
            device=device,
            parent_batch=8,
            path_chunk=64,
            estimator="mcv",
        )
        labels.append(value)
        metadata.append(meta)
    prediction = _predict_q(actor, t, x, device)
    target = 0.5 * (labels[0] + labels[1])
    raw = float(np.mean((prediction - target) ** 2))
    noise = float(np.mean((labels[0] - labels[1]) ** 2) / 4.0)
    corrected = max(raw - noise, 0.0)
    scale = max(float(np.sqrt(np.mean(target**2))), 1e-3)
    return {
        "actor": actor_name,
        "states": count,
        "paths_per_replication": paths,
        "time_steps": steps,
        "state_bank_sha256": array_sha256(t, x),
        "raw_component_rmse": math.sqrt(raw),
        "noise_corrected_component_rmse": math.sqrt(corrected),
        "normalized_total_inconsistency": corrected / scale**2,
        "raw_vector_rmse": math.sqrt(benchmark.d * raw),
        "noise_corrected_vector_rmse": math.sqrt(benchmark.d * corrected),
        "minimum_weight_ess_fraction": min(
            item["weight_ess_fraction"] for item in metadata
        ),
        "nonfinite_count": int(
            sum(item["nonfinite_count"] for item in metadata)
        ),
        "ground_truth_used": False,
    }


def _fit_is_finite(
    actor: TerminalAwareActor,
    metadata: dict[str, Any],
    trace: list[dict[str, float]],
) -> bool:
    values = [
        value
        for row in trace
        for value in row.values()
        if isinstance(value, (int, float))
    ]
    values.extend(
        value
        for value in metadata.values()
        if isinstance(value, (int, float))
    )
    return all(math.isfinite(float(value)) for value in values) and all(
        bool(torch.all(torch.isfinite(parameter)).item())
        for parameter in actor.parameters()
    )


def train(
    pde_id: str, d: int, seed: int, device: torch.device
) -> dict[str, Any]:
    ensure_directories()
    config = ActorConfig(**asdict(MATURE))
    if asdict(config) != asdict(MATURE):
        raise AssertionError("mature configuration changed")
    benchmark = make_benchmark(pde_id, d)
    run_id = _run_id(pde_id, d, seed)
    output = _training_path(run_id)
    actor = _initial_actor(run_id, benchmark, config, device)
    phases: list[dict[str, Any]] = []
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing["configuration"] != asdict(config):
            raise RuntimeError("existing training artifact has a different configuration")
        phases = existing.get("phases", [])
        if phases:
            actor = _load_actor(_checkpoint(run_id, len(phases)), benchmark, device)
        if existing.get("status") == "complete":
            return existing
    invocation_started = time.perf_counter()
    for iteration in range(len(phases) + 1, 5):
        phase_started = time.perf_counter()
        arrays, label_meta = _labels(
            run_id,
            pde_id,
            benchmark,
            actor,
            iteration,
            config,
            seed,
            device,
        )
        fallback_used = False
        fit_config = config
        first_failure: str | None = None
        try:
            fitted, fit_meta, trace = fit_actor(
                benchmark,
                fit_config,
                actor,
                arrays["train_t"],
                arrays["train_x"],
                arrays["train_y"],
                arrays["valid_t"],
                arrays["valid_x"],
                arrays["valid_y"],
                seed=stable_seed(run_id, iteration, "fit"),
                device=device,
            )
            if not _fit_is_finite(fitted, fit_meta, trace):
                raise FloatingPointError("non-finite fit result")
        except (FloatingPointError, RuntimeError) as error:
            first_failure = repr(error)
            fallback_used = True
            fit_config = replace(config, learning_rate=config.learning_rate / 2.0)
            fitted, fit_meta, trace = fit_actor(
                benchmark,
                fit_config,
                actor,
                arrays["train_t"],
                arrays["train_x"],
                arrays["train_y"],
                arrays["valid_t"],
                arrays["valid_x"],
                arrays["valid_y"],
                seed=stable_seed(run_id, iteration, "fit_fallback"),
                device=device,
            )
            if not _fit_is_finite(fitted, fit_meta, trace):
                raise FloatingPointError("predeclared half-LR fallback also non-finite")
        checkpoint = _checkpoint(run_id, iteration)
        _save_actor(checkpoint, fitted, benchmark, fit_config, iteration)
        phase = {
            "iteration": iteration,
            "checkpoint": str(checkpoint),
            "phase_seconds": time.perf_counter() - phase_started,
            "labels": label_meta,
            "training": fit_meta,
            "trace": trace,
            "fallback_used": fallback_used,
            "first_failure": first_failure,
            "effective_learning_rate": fit_config.learning_rate,
        }
        phases.append(phase)
        actor = fitted
        write_json_atomic(
            output,
            {
                "study_revision": REVISION,
                "run_id": run_id,
                "method": "terminal_aware_tangent_actor_mcv",
                "pde_id": pde_id,
                "dimension": d,
                "seed": seed,
                "configuration": asdict(config),
                "phases": phases,
                "status": "training_in_progress",
                "full_ambient_input": True,
                "oracle_projection": False,
                "ground_truth_used_for_training_or_selection": False,
            },
        )
        print(
            f"{run_id} a{iteration} "
            f"labels={label_meta['train']['wall_seconds']:.1f}s "
            f"fit={fit_meta['wall_seconds']:.1f}s "
            f"validation={fit_meta['selected_validation_normalized_mse']:.6g} "
            f"fallback={fallback_used}",
            flush=True,
        )
    actors = [_initial_actor(run_id, benchmark, config, device)] + [
        _load_actor(_checkpoint(run_id, iteration), benchmark, device)
        for iteration in range(1, 5)
    ]
    selector_started = time.perf_counter()
    fixed_point_rows = [
        _fixed_point(
            run_id,
            pde_id,
            benchmark,
            f"a{iteration}",
            candidate,
            seed,
            device,
        )
        for iteration, candidate in enumerate(actors)
    ]
    selector_seconds = time.perf_counter() - selector_started
    eligible = [row for row in fixed_point_rows if row["nonfinite_count"] == 0]
    if not eligible:
        raise FloatingPointError("every actor failed the reference-free selector")
    selected = min(
        eligible, key=lambda row: row["normalized_total_inconsistency"]
    )
    selected_index = int(selected["actor"][1:])
    record = {
        "study_revision": REVISION,
        "run_id": run_id,
        "method": "terminal_aware_tangent_actor_mcv",
        "pde_id": pde_id,
        "dimension": d,
        "seed": seed,
        "configuration": asdict(config),
        "phases": phases,
        "selection": {
            "criterion": "minimum independent noise-corrected fixed-point inconsistency",
            "rows": fixed_point_rows,
            "selected_actor": selected["actor"],
            "selected_checkpoint": None
            if selected_index == 0
            else str(_checkpoint(run_id, selected_index)),
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
            "phase_total_seconds": float(
                sum(phase["phase_seconds"] for phase in phases)
            ),
            "reference_free_selection_seconds": selector_seconds,
            "current_invocation_seconds": time.perf_counter() - invocation_started,
        },
        "fallback_count": int(sum(phase["fallback_used"] for phase in phases)),
        "full_ambient_input": True,
        "oracle_projection": False,
        "ground_truth_used_for_training_or_selection": False,
        "status": "complete",
    }
    write_json_atomic(output, record)
    return record


def _actors_for_evaluation(
    run_id: str,
    benchmark: Any,
    training: dict[str, Any],
    names: Iterable[str],
    device: torch.device,
) -> dict[str, TerminalAwareActor]:
    config = ActorConfig(**training["configuration"])
    result: dict[str, TerminalAwareActor] = {}
    for name in names:
        resolved = (
            training["selection"]["selected_actor"] if name == "selected" else name
        )
        if resolved == "a0":
            actor = _initial_actor(run_id, benchmark, config, device)
        else:
            actor = _load_actor(
                _checkpoint(run_id, int(resolved[1:])), benchmark, device
            )
        result[name] = actor
    return result


def evaluate(
    pde_id: str,
    d: int,
    seed: int,
    device: torch.device,
    *,
    actor_names: tuple[str, ...] = ("selected", "a0"),
    steps: tuple[int, ...] = (TA_PRIMARY_STEPS,),
) -> list[dict[str, Any]]:
    run_id = _run_id(pde_id, d, seed)
    training = json.loads(_training_path(run_id).read_text(encoding="utf-8"))
    if training.get("status") != "complete":
        raise RuntimeError("training and reference-free selection are incomplete")
    benchmark = make_benchmark(pde_id, d)
    actors = _actors_for_evaluation(
        run_id, benchmark, training, actor_names, device
    )
    t_all, x_all, is_validation = fixed_points(pde_id, d)
    t, x = t_all[~is_validation], x_all[~is_validation]
    # The analytic reference is opened only after the selector is complete.
    truth_state = benchmark.equation.exact_state(t, x)
    truth, truth_z = truth_state[:, 0], truth_state[:, 1:]
    rows: list[dict[str, Any]] = []
    for actor_name, actor in actors.items():
        actor_q = _predict_q(actor, t, x, device)
        gradient_error = actor_q - truth_z
        for time_steps in steps:
            for replicate in range(TA_REPLICATES):
                output_path = _evaluation_path(
                    run_id, actor_name, time_steps, replicate
                )
                if output_path.is_file():
                    rows.append(
                        json.loads(output_path.read_text(encoding="utf-8"))
                    )
                    continue
                output, meta = evaluate_policy(
                    benchmark,
                    actor,
                    t,
                    x,
                    paths=TA_PATHS,
                    time_steps=time_steps,
                    # Deliberately excludes actor_name: every actor receives
                    # identical Brownian paths for the diagnostic comparison.
                    seed=stable_seed(
                        "highdim_nosym_inference",
                        pde_id,
                        d,
                        seed,
                        time_steps,
                        replicate,
                    ),
                    device=device,
                    parent_batch=8,
                    path_chunk=64,
                )
                prediction = output["mcv"]
                variance = output["mcv_pair_variance"]
                error = prediction - truth
                rmse = float(np.sqrt(np.mean(error**2)))
                mc = float(
                    np.sqrt(np.mean(variance / (TA_PATHS // 2)))
                )
                prediction_path = (
                    PREDICTION_ROOT
                    / "actor"
                    / run_id
                    / actor_name
                    / f"steps{time_steps}_rep{replicate:02d}.npz"
                )
                prediction_path.parent.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(
                    prediction_path,
                    t=t,
                    x=x,
                    truth=truth,
                    truth_z=truth_z,
                    actor_q=actor_q,
                    prediction=prediction,
                    pair_variance=variance,
                )
                row = {
                    "study_revision": REVISION,
                    "method": "TA" if actor_name == "selected" else f"TA_{actor_name}",
                    "actor_request": actor_name,
                    "resolved_actor": training["selection"]["selected_actor"]
                    if actor_name == "selected"
                    else actor_name,
                    "pde_id": pde_id,
                    "dimension": d,
                    "training_seed": seed,
                    "inference_replicate": replicate,
                    "steps": time_steps,
                    "paths": TA_PATHS,
                    "value_rmse": rmse,
                    "relative_l2": float(
                        np.linalg.norm(error) / np.linalg.norm(truth)
                    ),
                    "skill": float(rmse / np.std(truth)),
                    "bias": float(np.mean(error)),
                    "mae": float(np.mean(np.abs(error))),
                    "gradient_vector_rmse": float(
                        np.sqrt(
                            np.mean(
                                np.sum(gradient_error * gradient_error, axis=1)
                            )
                        )
                    ),
                    "gradient_relative_l2": float(
                        np.linalg.norm(gradient_error) / np.linalg.norm(truth_z)
                    ),
                    "mc_rmse_contribution": mc,
                    "non_mc_rmse_proxy": math.sqrt(max(rmse * rmse - mc * mc, 0.0)),
                    "variance_reduction": meta["variance_reduction"],
                    "weight_ess_fraction": meta["weight_ess_fraction"],
                    "inference_seconds": meta["wall_seconds"],
                    "peak_gpu_bytes": meta["peak_gpu_bytes"],
                    "nonfinite_count": meta["nonfinite_count"],
                    "test_states": len(t),
                    "test_bank_sha256": array_sha256(t_all, x_all, is_validation),
                    "prediction_path": str(prediction_path),
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
    parser.add_argument("--pde", choices=("S05", "B6"), required=True)
    parser.add_argument("--dimension", type=int, choices=(100, 400), required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--train-only", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--actors", default="selected,a0")
    parser.add_argument("--steps", default=str(TA_PRIMARY_STEPS))
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation required")
    if not args.evaluate_only:
        result = train(args.pde, args.dimension, args.seed, device)
        print(
            json.dumps(
                {
                    "training": result["status"],
                    "selected": result["selection"]["selected_actor"],
                }
            ),
            flush=True,
        )
    if not args.train_only:
        evaluate(
            args.pde,
            args.dimension,
            args.seed,
            device,
            actor_names=tuple(value for value in args.actors.split(",") if value),
            steps=tuple(int(value) for value in args.steps.split(",") if value),
        )


if __name__ == "__main__":
    main()
