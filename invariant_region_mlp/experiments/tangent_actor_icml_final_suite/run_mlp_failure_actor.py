"""Terminal-aware Tangent-Actor on the nonsmooth MLP-failure PDE.

The norm is linearized with its minimum-norm subgradient at zero.  This is an
explicitly adapted extension and is never presented as a smooth-driver theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything, synchronize,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.models import (
    TerminalAwareActor,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    ActorConfig, evaluate_policy, fit_actor, generate_labels,
)

from .benchmarks import MLPFailureBenchmark
from .prepare_suite import write_json
from .protocol import RESULTS_DIR, stable_seed


CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
CACHE_DIR = RESULTS_DIR / "cache"
PREDICTION_DIR = RESULTS_DIR / "predictions"
RAW_DIR = RESULTS_DIR / "raw"
INPUT_DIR = RESULTS_DIR / "inputs"


MATURE = ActorConfig(
    training_points=8192,
    validation_points=512,
    paths_per_label=128,
    label_time_steps=32,
    optimizer_steps=5000,
    width=384,
    depth=4,
    activation="silu",
    learning_rate=1.5e-4,
    batch_size=512,
    parent_batch=8,
    path_chunk=32,
    checkpoint_interval=250,
    gradient_clip=10.0,
    training_distribution="function",
    label_estimator="mcv",
)


SMOKE = ActorConfig(
    training_points=32, validation_points=16, paths_per_label=8,
    label_time_steps=4, optimizer_steps=4, width=16, depth=2,
    activation="silu", learning_rate=2e-4, batch_size=16,
    parent_batch=4, path_chunk=4, checkpoint_interval=2,
    gradient_clip=10.0, training_distribution="function", label_estimator="mcv",
)


def _dirs() -> None:
    for path in (CHECKPOINT_DIR, CACHE_DIR, PREDICTION_DIR, RAW_DIR):
        path.mkdir(parents=True, exist_ok=True)


def _actor_path(d: int, seed: int, iteration: int, tag: str) -> Path:
    return CHECKPOINT_DIR / f"mlp_failure_d{d}_{tag}_s{seed}_a{iteration}.pt"


def _save_actor(
    path: Path, actor: TerminalAwareActor, *, d: int, seed: int,
    iteration: int, config: ActorConfig,
) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save({
        "state_dict": {key: value.detach().cpu() for key, value in actor.state_dict().items()},
        "dimension": d, "seed": seed, "iteration": iteration,
        "width": config.width, "depth": config.depth, "activation": config.activation,
        "config": asdict(config),
        "subgradient_rule": "inactive q/norm(q) when norm>1e-12; minimum-norm zero otherwise",
    }, temporary)
    temporary.replace(path)


def _load_actor(
    path: Path, benchmark: MLPFailureBenchmark, device: torch.device,
) -> TerminalAwareActor:
    payload = torch.load(path, map_location=device, weights_only=False)
    actor = TerminalAwareActor(
        benchmark, width=int(payload["width"]), depth=int(payload["depth"]),
        activation=str(payload["activation"]),
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor


def _state_bank(
    benchmark: MLPFailureBenchmark, count: int, *, seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    return benchmark.sample_states(count, seed=seed, purpose="function")


def _labels(
    benchmark: MLPFailureBenchmark, source: TerminalAwareActor,
    config: ActorConfig, *, seed: int, iteration: int, device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    train_t, train_x = _state_bank(
        benchmark, config.training_points,
        seed=stable_seed("mlp_failure_train_states", benchmark.d, seed, iteration),
    )
    valid_t, valid_x = _state_bank(
        benchmark, config.validation_points,
        seed=stable_seed("mlp_failure_label_validation_states", benchmark.d, seed),
    )
    kwargs = dict(
        paths=config.paths_per_label, time_steps=config.label_time_steps,
        device=device, parent_batch=config.parent_batch, path_chunk=config.path_chunk,
        estimator=config.label_estimator,
    )
    train, train_meta = generate_labels(
        benchmark, source, train_t, train_x,
        seed=stable_seed("mlp_failure_labels", benchmark.d, seed, iteration, "train"),
        **kwargs,
    )
    valid, valid_meta = generate_labels(
        benchmark, source, valid_t, valid_x,
        seed=stable_seed("mlp_failure_labels", benchmark.d, seed, iteration, "valid"),
        **kwargs,
    )
    replicate, replicate_meta = generate_labels(
        benchmark, source, valid_t, valid_x,
        seed=stable_seed("mlp_failure_labels", benchmark.d, seed, iteration, "replicate"),
        **kwargs,
    )
    delta = valid - replicate
    return {
        "train_t": train_t, "train_x": train_x, "train_labels": train,
        "valid_t": valid_t, "valid_x": valid_x, "valid_labels": valid,
    }, {
        "train": train_meta, "validation": valid_meta, "replication": replicate_meta,
        "validation_replication_component_rmse": float(np.sqrt(np.mean(np.square(delta)))),
        "validation_replication_vector_rmse": float(
            np.sqrt(np.mean(np.square(delta).sum(axis=1)))
        ),
        "ground_truth_used": False,
    }


def _predict_q(
    actor: TerminalAwareActor, t: np.ndarray, x: np.ndarray, device: torch.device,
) -> np.ndarray:
    parts = []
    actor.eval()
    with torch.no_grad():
        for begin in range(0, len(t), 1024):
            end = min(begin + 1024, len(t))
            _, q = actor.predict(
                torch.as_tensor(t[begin:end], dtype=torch.float32, device=device),
                torch.as_tensor(x[begin:end], dtype=torch.float32, device=device),
            )
            if q is None:
                raise RuntimeError("Gradient Actor returned no q")
            parts.append(q.double().cpu().numpy())
    return np.concatenate(parts)


def train(
    d: int, seed: int, device: torch.device, *, smoke: bool = False,
) -> dict[str, Any]:
    _dirs()
    tag = "smoke" if smoke else "mature"
    config = SMOKE if smoke else MATURE
    output = RAW_DIR / f"training_mlp_failure_d{d}_{tag}_s{seed}.json"
    benchmark = MLPFailureBenchmark(d)
    seed_everything(stable_seed("mlp_failure_initial_actor", d, seed, tag))
    actor = TerminalAwareActor(
        benchmark, width=config.width, depth=config.depth, activation=config.activation,
    ).to(device)
    actor.eval()
    phases: list[dict[str, Any]] = []
    if output.exists():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing["configuration"] != asdict(config):
            raise RuntimeError("existing training artifact has a different configuration")
        phases = existing["phases"]
        if phases:
            actor = _load_actor(Path(phases[-1]["checkpoint"]), benchmark, device)
        if existing.get("status") == "complete":
            return existing
    invocation = time.perf_counter()
    for iteration in range(len(phases) + 1, 5):
        phase_start = time.perf_counter()
        arrays, label_meta = _labels(
            benchmark, actor, config, seed=seed, iteration=iteration, device=device,
        )
        fitted, fit_meta, trace = fit_actor(
            benchmark, config, actor,
            arrays["train_t"], arrays["train_x"], arrays["train_labels"],
            arrays["valid_t"], arrays["valid_x"], arrays["valid_labels"],
            seed=stable_seed("mlp_failure_fit", d, seed, iteration), device=device,
        )
        checkpoint = _actor_path(d, seed, iteration, tag)
        _save_actor(
            checkpoint, fitted, d=d, seed=seed, iteration=iteration, config=config,
        )
        synchronize(device)
        phases.append({
            "iteration": iteration, "checkpoint": str(checkpoint),
            "phase_seconds": time.perf_counter() - phase_start,
            "labels": label_meta, "training": fit_meta, "trace": trace,
        })
        actor = fitted
        record = {
            "method": "terminal_aware_tangent_actor_mcv",
            "pde": "mlp_failure", "dimension": d, "seed": seed, "tag": tag,
            "configuration": asdict(config), "phases": phases,
            "iterations": [f"a{i}" for i in range(5)],
            "subgradient_rule": "minimum-norm b(0)=0; normalized inactive q away from zero",
            "theoretical_status": "adapted nonsmooth extension; no smooth tangent convergence claim",
            "ground_truth_used_for_training_or_checkpoint_selection": False,
            "status": "training_in_progress",
        }
        write_json(output, record)
        print(
            f"mlp_failure d={d} seed={seed} a{iteration} "
            f"label={label_meta['train']['wall_seconds']:.1f}s "
            f"fit={fit_meta['wall_seconds']:.1f}s val={fit_meta['selected_validation_rmse']:.6g}",
            flush=True,
        )

    # A fixed, independent linear-evaluation label bank ranks a0..a4 without
    # opening the analytic solution.  Two label replications are averaged.
    select_count = 32 if smoke else 512
    select_paths = 8 if smoke else 512
    select_steps = 4 if smoke else 64
    select_t, select_x = _state_bank(
        benchmark, select_count,
        seed=stable_seed("mlp_failure_selection_states", d, seed),
    )
    seed_everything(stable_seed("mlp_failure_initial_actor", d, seed, tag))
    a0 = TerminalAwareActor(
        benchmark, width=config.width, depth=config.depth, activation=config.activation,
    ).to(device)
    a0.eval()
    selection_labels = []
    selection_meta = []
    for replicate in (0, 1):
        labels, meta = generate_labels(
            benchmark, a0, select_t, select_x,
            paths=select_paths, time_steps=select_steps,
            seed=stable_seed("mlp_failure_selection_labels", d, seed, replicate),
            device=device, parent_batch=config.parent_batch,
            path_chunk=min(config.path_chunk, select_paths), estimator="mcv",
        )
        selection_labels.append(labels)
        selection_meta.append(meta)
    target = np.mean(selection_labels, axis=0)
    actors = [a0] + [
        _load_actor(_actor_path(d, seed, iteration, tag), benchmark, device)
        for iteration in range(1, 5)
    ]
    ranking = []
    for iteration, candidate in enumerate(actors):
        prediction = _predict_q(candidate, select_t, select_x, device)
        error = prediction - target
        ranking.append({
            "actor": f"a{iteration}",
            "fixed_label_component_rmse": float(np.sqrt(np.mean(np.square(error)))),
            "fixed_label_vector_rmse": float(
                np.sqrt(np.mean(np.square(error).sum(axis=1)))
            ),
            "checkpoint": None if iteration == 0 else str(_actor_path(d, seed, iteration, tag)),
        })
    selected = min(ranking, key=lambda row: row["fixed_label_vector_rmse"])
    record = {
        "method": "terminal_aware_tangent_actor_mcv",
        "pde": "mlp_failure", "dimension": d, "seed": seed, "tag": tag,
        "configuration": asdict(config), "phases": phases,
        "selection": {
            "criterion": "minimum vector RMSE to mean of two independent frozen a0 linear-evaluation gradient-label banks",
            "ground_truth_used": False, "states": select_count,
            "paths_per_replication": select_paths, "time_steps": select_steps,
            "label_replication_vector_rmse": float(
                np.sqrt(np.mean(np.square(selection_labels[0] - selection_labels[1]).sum(axis=1)))
            ),
            "label_metadata": selection_meta, "ranking": ranking,
            "selected_actor": selected["actor"], "selected_checkpoint": selected["checkpoint"],
        },
        "subgradient_rule": "minimum-norm b(0)=0; normalized inactive q away from zero",
        "theoretical_status": "adapted nonsmooth extension; no smooth tangent convergence claim",
        "ground_truth_used_for_training_or_checkpoint_selection": False,
        "invocation_seconds": time.perf_counter() - invocation,
        "status": "complete",
    }
    write_json(output, record)
    return record


def _load_test_bank(d: int) -> tuple[np.ndarray, np.ndarray, str]:
    path = INPUT_DIR / f"test_bank_mlp_failure_d{d}.npz"
    with np.load(path) as data:
        return data["t"], data["x"], str(path)


def evaluate(
    d: int, seed: int, device: torch.device, *, smoke: bool = False,
    paths: int = 2048, steps: tuple[int, ...] = (64, 128, 256, 512),
) -> dict[str, Any]:
    tag = "smoke" if smoke else "mature"
    training_path = RAW_DIR / f"training_mlp_failure_d{d}_{tag}_s{seed}.json"
    training = json.loads(training_path.read_text(encoding="utf-8"))
    benchmark = MLPFailureBenchmark(d)
    selected_actor = training["selection"]["selected_actor"]
    if selected_actor == "a0":
        config = ActorConfig(**training["configuration"])
        seed_everything(stable_seed("mlp_failure_initial_actor", d, seed, tag))
        actor = TerminalAwareActor(
            benchmark, width=config.width, depth=config.depth, activation=config.activation,
        ).to(device)
        actor.eval()
    else:
        actor = _load_actor(Path(training["selection"]["selected_checkpoint"]), benchmark, device)
    test_t, test_x, bank_path = _load_test_bank(d)
    if smoke:
        test_t, test_x = test_t[:16], test_x[:16]
        paths, steps = 8, (4,)
    truth = benchmark.exact_value_numpy(test_t, test_x)
    exact_z = benchmark.exact_z_numpy(test_t, test_x)
    predicted_q = _predict_q(actor, test_t, test_x, device)
    gradient_error = predicted_q - exact_z
    rows = []
    predictions: dict[str, np.ndarray] = {
        "t": test_t, "x": test_x, "truth": truth, "exact_z": exact_z,
        "actor_q": predicted_q,
    }
    for time_steps in steps:
        output, meta = evaluate_policy(
            benchmark, actor, test_t, test_x,
            paths=paths, time_steps=time_steps,
            seed=stable_seed("mlp_failure_postfreeze", d, seed, paths, time_steps),
            device=device, parent_batch=8 if not smoke else 4,
            path_chunk=min(64 if not smoke else 4, paths),
        )
        for estimator in ("plain", "mcv"):
            prediction = output[estimator]
            error = prediction - truth
            pair_variance = output[f"{estimator}_pair_variance"]
            rows.append({
                "steps": time_steps, "paths": paths, "estimator": estimator,
                "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
                "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
                "bias": float(np.mean(error)), "mae": float(np.mean(np.abs(error))),
                "mc_rmse_contribution": float(
                    np.sqrt(np.mean(pair_variance / (paths // 2)))
                ),
                "variance_reduction": meta["variance_reduction"],
                "weight_ess_fraction": meta["weight_ess_fraction"],
                "inference_seconds": meta["wall_seconds"],
                "peak_gpu_bytes": meta["peak_gpu_bytes"],
                "nonfinite_count": meta["nonfinite_count"],
            })
            predictions[f"{estimator}_n{time_steps}"] = prediction
            predictions[f"{estimator}_pair_variance_n{time_steps}"] = pair_variance
    prediction_path = PREDICTION_DIR / f"mlp_failure_d{d}_{tag}_s{seed}.npz"
    np.savez_compressed(prediction_path, **predictions)
    result = {
        "method": "terminal_aware_tangent_actor_mcv", "pde": "mlp_failure",
        "dimension": d, "seed": seed, "tag": tag,
        "selected_actor": selected_actor, "rows": rows,
        "gradient_actor_component_rmse": float(np.sqrt(np.mean(np.square(gradient_error)))),
        "gradient_actor_vector_rmse": float(
            np.sqrt(np.mean(np.square(gradient_error).sum(axis=1)))
        ),
        "gradient_actor_inactive_rms": float(
            np.sqrt(np.mean(np.square(predicted_q[:, 1:])))
        ),
        "test_bank": bank_path, "states": len(test_t),
        "prediction_path": str(prediction_path),
        "ground_truth_loaded_only_after_checkpoint_selection": True,
        "status": "complete",
    }
    write_json(RAW_DIR / f"evaluation_mlp_failure_d{d}_{tag}_s{seed}.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--train-only", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--evaluation-paths", type=int, default=2048)
    parser.add_argument("--evaluation-steps", default="64,128,256,512")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation required")
    if not args.evaluate_only:
        training = train(args.dimension, args.seed, device, smoke=args.smoke)
        print(json.dumps({"training": training["status"], "selected": training["selection"]["selected_actor"]}), flush=True)
    if not args.train_only:
        result = evaluate(
            args.dimension, args.seed, device, smoke=args.smoke,
            paths=args.evaluation_paths,
            steps=tuple(int(value) for value in args.evaluation_steps.split(",")),
        )
        print(json.dumps({"evaluation": result["status"], "rows": result["rows"]}), flush=True)


if __name__ == "__main__":
    main()
