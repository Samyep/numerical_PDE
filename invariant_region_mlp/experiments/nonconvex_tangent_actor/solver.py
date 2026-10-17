"""MCV tangent labels, Actor fitting, fixed-point diagnostics, and evaluation."""

from __future__ import annotations

from dataclasses import asdict
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    antithetic_normals,
    frozen,
    make_targets,
    tangent_labels,
)
from invariant_region_mlp.experiments.tangent_actor_mcv_training.labels import (
    _weighted_rollout_chunk_mcv_labels,
    mcv_tangent_labels,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.models import (
    TerminalAwareActor,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
    train_actor,
)

from .common import (
    CACHE_DIR,
    CHECKPOINT_DIR,
    array_sha256,
    atomic_json,
    atomic_npz,
    atomic_torch_save,
    seed_everything,
    stable_seed,
    state_bank,
    synchronize,
)
from .drivers import SineQuadraticDriver


def initial_actor(
    driver: SineQuadraticDriver, config: FitConfig, *, seed: int, device: torch.device
) -> TerminalAwareActor:
    """The zero-residual terminal extension used as a0."""
    seed_everything(seed)
    actor = TerminalAwareActor(
        driver, width=config.width, depth=config.depth
    ).to(device)
    actor.eval()
    return actor


def save_actor(
    path: Path,
    actor: TerminalAwareActor,
    *,
    driver: SineQuadraticDriver,
    config: FitConfig,
    metadata: dict[str, Any],
) -> None:
    atomic_torch_save(
        path,
        {
            "state_dict": {
                key: value.detach().cpu().clone()
                for key, value in actor.state_dict().items()
            },
            "driver": driver.metadata(),
            "fit_config": asdict(config),
            "metadata": metadata,
        },
    )


def load_actor(
    path: Path, driver: SineQuadraticDriver, config: FitConfig, device: torch.device
) -> TerminalAwareActor:
    payload = torch.load(path, map_location=device, weights_only=True)
    actor = TerminalAwareActor(
        driver, width=config.width, depth=config.depth
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor


def actor_prediction(
    driver: SineQuadraticDriver,
    actor: nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
    *,
    chunk: int = 2048,
) -> tuple[np.ndarray, np.ndarray]:
    r_parts: list[np.ndarray] = []
    q_parts: list[np.ndarray] = []
    with torch.no_grad():
        for begin in range(0, len(t), chunk):
            end = min(begin + chunk, len(t))
            r, q = actor.predict(
                torch.as_tensor(t[begin:end], dtype=torch.float32, device=device),
                torch.as_tensor(x[begin:end], dtype=torch.float32, device=device),
                driver.T,
            )
            if q is None:
                raise RuntimeError("joint driver Actor omitted q")
            r_parts.append(r.double().cpu().numpy())
            q_parts.append(q.double().cpu().numpy())
    return np.concatenate(r_parts), np.concatenate(q_parts)


def _label_function(estimator: str):
    if estimator == "mcv":
        return mcv_tangent_labels
    if estimator == "plain":
        return tangent_labels
    raise KeyError(estimator)


def generate_phase_labels(
    run_id: str,
    driver: SineQuadraticDriver,
    source: nn.Module,
    iteration: int,
    config: FitConfig,
    *,
    training_seed: int,
    estimator: str,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Generate independent train/validation/replicate tangent labels."""
    stem = f"{run_id}_{estimator}_a{iteration}_labels"
    cache_path = CACHE_DIR / f"{stem}.npz"
    metadata_path = CACHE_DIR / f"{stem}.json"
    if cache_path.exists() and metadata_path.exists():
        with np.load(cache_path) as archive:
            arrays = {key: archive[key] for key in archive.files}
        import json

        return arrays, json.loads(metadata_path.read_text(encoding="utf-8"))

    train_t, train_x = state_bank(
        driver.name,
        driver.d,
        f"{run_id}_train_a{iteration}",
        config.training_points,
        T=driver.T,
        seed=training_seed,
    )
    valid_t, valid_x = state_bank(
        driver.name,
        driver.d,
        f"{run_id}_fixed_label_validation",
        config.validation_points,
        T=driver.T,
        seed=training_seed,
    )
    kwargs = {
        "paths": config.paths_per_label,
        "time_steps": config.label_time_steps,
        "device": device,
        "parent_batch": config.parent_batch,
        "path_chunk": config.path_chunk,
    }
    labeler = _label_function(estimator)
    train_value, train_q, train_meta = labeler(
        driver,
        source,
        train_t,
        train_x,
        seed=stable_seed(run_id, estimator, iteration, "train_labels"),
        **kwargs,
    )
    valid_value, valid_q, valid_meta = labeler(
        driver,
        source,
        valid_t,
        valid_x,
        seed=stable_seed(run_id, estimator, iteration, "valid_labels"),
        **kwargs,
    )
    replicate_value, replicate_q, replicate_meta = labeler(
        driver,
        source,
        valid_t,
        valid_x,
        seed=stable_seed(run_id, estimator, iteration, "replicate_labels"),
        **kwargs,
    )
    if train_q is None or valid_q is None or replicate_q is None:
        raise RuntimeError("joint Actor labels require gradients")
    train_r, train_q, train_projection = make_targets(
        driver, train_value, train_q
    )
    valid_r, valid_q, valid_projection = make_targets(
        driver, valid_value, valid_q
    )
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
    difference_r = valid_value - replicate_value
    difference_q = valid_q - replicate_q
    metadata = {
        "run_id": run_id,
        "iteration": iteration,
        "source": f"a{iteration - 1}",
        "estimator": estimator,
        "train": train_meta,
        "validation": valid_meta,
        "replication": replicate_meta,
        "train_projection": train_projection,
        "validation_projection": valid_projection,
        "value_label_replication_rmse": float(
            np.sqrt(np.mean(np.square(difference_r)))
        ),
        "gradient_label_replication_component_rmse": float(
            np.sqrt(np.mean(np.square(difference_q)))
        ),
        "gradient_label_replication_vector_rmse": float(
            np.sqrt(np.mean(np.square(difference_q).sum(axis=1)))
        ),
        "train_bank_sha256": array_sha256(train_t, train_x),
        "validation_bank_sha256": array_sha256(valid_t, valid_x),
        "full_ambient_training_states": True,
        "oracle_projection_used": False,
        "ground_truth_used": False,
    }
    atomic_npz(cache_path, **arrays)
    atomic_json(metadata_path, metadata)
    return arrays, metadata


def fit_phase(
    run_id: str,
    driver: SineQuadraticDriver,
    source: TerminalAwareActor,
    iteration: int,
    config: FitConfig,
    arrays: dict[str, np.ndarray],
    label_metadata: dict[str, Any],
    *,
    training_seed: int,
    estimator: str,
    device: torch.device,
) -> tuple[TerminalAwareActor, dict[str, Any]]:
    checkpoint = CHECKPOINT_DIR / f"{run_id}_{estimator}_a{iteration}.pt"
    metadata_path = CHECKPOINT_DIR / f"{run_id}_{estimator}_a{iteration}.json"
    if checkpoint.exists() and metadata_path.exists():
        import json

        return (
            load_actor(checkpoint, driver, config, device),
            json.loads(metadata_path.read_text(encoding="utf-8")),
        )
    actor, training, trace = train_actor(
        driver,
        config,
        "terminal_aware",
        arrays["train_t"],
        arrays["train_x"],
        arrays["train_r"],
        arrays["train_q"],
        arrays["valid_t"],
        arrays["valid_x"],
        arrays["valid_r"],
        arrays["valid_q"],
        seed=stable_seed(run_id, estimator, iteration, "actor_fit", training_seed),
        device=device,
        initial=source,
    )
    metadata = {
        "run_id": run_id,
        "iteration": iteration,
        "source": f"a{iteration - 1}",
        "estimator": estimator,
        "labels": label_metadata,
        "training": training,
        "trace": trace,
        "checkpoint_selection": "minimum held-out noisy tangent-label loss",
        "ground_truth_used": False,
    }
    save_actor(
        checkpoint,
        actor,
        driver=driver,
        config=config,
        metadata=metadata,
    )
    atomic_json(metadata_path, metadata)
    return actor, metadata


def fixed_point_diagnostic(
    run_id: str,
    actor_name: str,
    driver: SineQuadraticDriver,
    actor: nn.Module,
    *,
    states: int,
    paths: int,
    time_steps: int,
    training_seed: int,
    device: torch.device,
    value_scale_floor: float,
    gradient_component_scale_floor: float,
    audit_tag: str = "primary",
) -> dict[str, Any]:
    """Independent noise-corrected A-(v_A,z_vA) consistency diagnostic."""
    t, x = state_bank(
        driver.name,
        driver.d,
        f"{run_id}_fixed_point_states",
        states,
        T=driver.T,
        seed=20261017,
    )
    kwargs = {
        "paths": paths,
        "time_steps": time_steps,
        "device": device,
        "parent_batch": 8,
        "path_chunk": min(64, paths),
    }
    first_r, first_q, first_meta = mcv_tangent_labels(
        driver,
        actor,
        t,
        x,
        seed=stable_seed(run_id, audit_tag, "common_fp_replication_1"),
        **kwargs,
    )
    second_r, second_q, second_meta = mcv_tangent_labels(
        driver,
        actor,
        t,
        x,
        seed=stable_seed(run_id, audit_tag, "common_fp_replication_2"),
        **kwargs,
    )
    if first_q is None or second_q is None:
        raise RuntimeError("fixed-point diagnostic requires q labels")
    actor_r, actor_q = actor_prediction(driver, actor, t, x, device)
    target_r = 0.5 * (first_r + second_r)
    target_q = 0.5 * (first_q + second_q)
    value_raw_mse = float(np.mean(np.square(actor_r - target_r)))
    gradient_raw_component_mse = float(np.mean(np.square(actor_q - target_q)))
    # With two independent label replications, (Y1-Y2)^2/4 estimates the
    # sampling variance of their average and can be removed from raw MSE.
    value_noise_mse = float(np.mean(np.square(first_r - second_r)) / 4.0)
    gradient_noise_component_mse = float(
        np.mean(np.square(first_q - second_q)) / 4.0
    )
    value_corrected_mse = max(value_raw_mse - value_noise_mse, 0.0)
    gradient_corrected_component_mse = max(
        gradient_raw_component_mse - gradient_noise_component_mse, 0.0
    )
    value_scale = max(float(np.std(target_r)), value_scale_floor)
    gradient_scale = max(
        float(np.sqrt(np.mean(np.square(target_q)))),
        gradient_component_scale_floor,
    )
    normalized_value = value_corrected_mse / value_scale**2
    normalized_gradient = gradient_corrected_component_mse / gradient_scale**2
    return {
        "run_id": run_id,
        "actor": actor_name,
        "audit_tag": audit_tag,
        "states": states,
        "paths_per_replication": paths,
        "time_steps": time_steps,
        "state_bank_sha256": array_sha256(t, x),
        "value_raw_rmse": math.sqrt(value_raw_mse),
        "value_label_replication_rmse": float(
            np.sqrt(np.mean(np.square(first_r - second_r)))
        ),
        "value_noise_corrected_rmse": math.sqrt(value_corrected_mse),
        "gradient_raw_component_rmse": math.sqrt(gradient_raw_component_mse),
        "gradient_raw_vector_rmse": math.sqrt(
            driver.d * gradient_raw_component_mse
        ),
        "gradient_label_replication_component_rmse": float(
            np.sqrt(np.mean(np.square(first_q - second_q)))
        ),
        "gradient_noise_corrected_component_rmse": math.sqrt(
            gradient_corrected_component_mse
        ),
        "gradient_noise_corrected_vector_rmse": math.sqrt(
            driver.d * gradient_corrected_component_mse
        ),
        "value_scale": value_scale,
        "gradient_component_scale": gradient_scale,
        "normalized_value_inconsistency": normalized_value,
        "normalized_gradient_inconsistency": normalized_gradient,
        "normalized_total_inconsistency": normalized_value + normalized_gradient,
        "first_weight_ess_fraction": first_meta["weight_ess_fraction"],
        "second_weight_ess_fraction": second_meta["weight_ess_fraction"],
        "minimum_weight_ess_fraction": min(
            first_meta["weight_ess_fraction"], second_meta["weight_ess_fraction"]
        ),
        "nonfinite_count": int(
            first_meta["nonfinite_count"] + second_meta["nonfinite_count"]
        ),
        "wall_seconds": float(
            first_meta["wall_seconds"] + second_meta["wall_seconds"]
        ),
        "selection_eligible": bool(
            first_meta["nonfinite_count"] == 0
            and second_meta["nonfinite_count"] == 0
            and np.isfinite(normalized_value + normalized_gradient)
        ),
        "ground_truth_used": False,
    }


def evaluate_policy(
    driver: SineQuadraticDriver,
    actor: nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    time_steps: int,
    seed: int,
    device: torch.device,
    parent_batch: int = 8,
    path_chunk: int = 64,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Evaluate plain and alpha=1 MCV estimates on exactly shared paths."""
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("paths must be divisible by an even path chunk")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, driver.d)
    plain_mean = np.empty(len(times), dtype=np.float64)
    mcv_mean = np.empty(len(times), dtype=np.float64)
    plain_variance = np.empty(len(times), dtype=np.float64)
    mcv_variance = np.empty(len(times), dtype=np.float64)
    generator = torch.Generator(device=device).manual_seed(seed)
    weight_sum = weight_square_sum = 0.0
    weight_min, weight_max = float("inf"), -float("inf")
    nonfinite = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(actor):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            plain_pairs: list[torch.Tensor] = []
            mcv_pairs: list[torch.Tensor] = []
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    (end - begin, path_chunk, time_steps, driver.d),
                    generator=generator,
                    device=device,
                )
                corrected, _, weight, martingale, diagnostics = (
                    _weighted_rollout_chunk_mcv_labels(
                        driver, actor, bt, bx, normals, need_z=False
                    )
                )
                plain = corrected + martingale
                half = path_chunk // 2
                plain_pairs.append(0.5 * (plain[:, :half] + plain[:, half:]))
                mcv_pairs.append(
                    0.5 * (corrected[:, :half] + corrected[:, half:])
                )
                weight_sum += float(weight.double().sum().item())
                weight_square_sum += float(weight.double().square().sum().item())
                weight_min = min(weight_min, float(weight.min().item()))
                weight_max = max(weight_max, float(weight.max().item()))
                nonfinite += int(diagnostics["nonfinite_count"])
            local_plain = torch.cat(plain_pairs, dim=1).double()
            local_mcv = torch.cat(mcv_pairs, dim=1).double()
            plain_mean[begin:end] = local_plain.mean(dim=1).cpu().numpy()
            mcv_mean[begin:end] = local_mcv.mean(dim=1).cpu().numpy()
            plain_variance[begin:end] = local_plain.var(
                dim=1, unbiased=True
            ).cpu().numpy()
            mcv_variance[begin:end] = local_mcv.var(
                dim=1, unbiased=True
            ).cpu().numpy()
    synchronize(device)
    elapsed = time.perf_counter() - started
    count = len(times) * paths
    if nonfinite or not np.all(np.isfinite(mcv_mean)):
        raise FloatingPointError("nonfinite tangent-policy evaluation")
    return {
        "plain": plain_mean,
        "mcv": mcv_mean,
        "plain_pair_variance": plain_variance,
        "mcv_pair_variance": mcv_variance,
    }, {
        "paths": paths,
        "antithetic_pairs": paths // 2,
        "time_steps": time_steps,
        "seed": seed,
        "mcv_alpha": 1.0,
        "mean_plain_pair_variance": float(np.mean(plain_variance)),
        "mean_mcv_pair_variance": float(np.mean(mcv_variance)),
        "variance_reduction": float(
            np.mean(plain_variance) / max(np.mean(mcv_variance), 1e-300)
        ),
        "mean_plain_mc_se": float(
            np.sqrt(np.sum(plain_variance / (paths // 2))) / len(times)
        ),
        "mean_mcv_mc_se": float(
            np.sqrt(np.sum(mcv_variance / (paths // 2))) / len(times)
        ),
        "weight_ess_fraction": float(
            weight_sum**2
            / max(count * weight_square_sum, np.finfo(float).tiny)
        ),
        "weight_min": weight_min,
        "weight_max": weight_max,
        "nonfinite_count": nonfinite,
        "wall_seconds": elapsed,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
    }


__all__ = [
    "actor_prediction",
    "evaluate_policy",
    "fit_phase",
    "fixed_point_diagnostic",
    "generate_phase_labels",
    "initial_actor",
    "load_actor",
    "save_actor",
]
