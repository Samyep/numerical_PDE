"""Actor fitting with dense optimization and generalization diagnostics."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
)

from .common import seed_everything
from .models import ActorKind, count_parameters, make_actor


@dataclass(frozen=True)
class FitConfig:
    training_points: int
    validation_points: int = 256
    paths_per_label: int = 128
    label_time_steps: int = 16
    optimizer_steps: int = 2500
    width: int = 192
    depth: int = 3
    learning_rate: float = 2e-4
    batch_size: int = 512
    parent_batch: int = 8
    path_chunk: int = 32
    checkpoint_interval: int = 125
    gradient_clip: float = 10.0


def _fit_metrics(
    actor: nn.Module,
    driver: ConcaveSemilinearDriver,
    t: torch.Tensor,
    x: torch.Tensor,
    r: torch.Tensor,
    q: torch.Tensor | None,
    *,
    r_scale: float,
    q_scale: float | None,
    chunk_size: int = 1024,
) -> dict[str, float | None]:
    count = len(t)
    r_square = q_square = q_vector_square = 0.0
    lower_saturation = upper_saturation = 0
    actor.eval()
    with torch.no_grad():
        for begin in range(0, count, chunk_size):
            end = min(begin + chunk_size, count)
            predicted_r, predicted_q = actor.predict(t[begin:end], x[begin:end], driver.T)
            r_error = predicted_r - r[begin:end]
            r_square += float(r_error.double().square().sum().item())
            unit = (predicted_r - driver.r_low) / (driver.r_high - driver.r_low)
            lower_saturation += int((unit < 0.01).sum().item())
            upper_saturation += int((unit > 0.99).sum().item())
            if q is not None:
                assert predicted_q is not None
                q_error = predicted_q - q[begin:end]
                q_square += float(q_error.double().square().sum().item())
                q_vector_square += float(q_error.double().square().sum(dim=1).sum().item())
    value_mse = r_square / count
    normalized = value_mse / r_scale**2
    if q is None:
        gradient_mse = gradient_component_rmse = gradient_vector_rmse = None
    else:
        assert q_scale is not None
        gradient_mse = q_square / (count * driver.d)
        gradient_component_rmse = float(np.sqrt(gradient_mse))
        gradient_vector_rmse = float(np.sqrt(q_vector_square / count))
        normalized += gradient_mse / q_scale**2
    return {
        "normalized_loss": float(normalized),
        "value_mse": float(value_mse),
        "value_rmse": float(np.sqrt(value_mse)),
        "gradient_mse": gradient_mse,
        "gradient_component_rmse": gradient_component_rmse,
        "gradient_vector_rmse": gradient_vector_rmse,
        "lower_saturation_fraction": lower_saturation / count,
        "upper_saturation_fraction": upper_saturation / count,
    }


def _tensor(value: np.ndarray, device: torch.device) -> torch.Tensor:
    return torch.as_tensor(value, dtype=torch.float32, device=device)


def train_actor(
    driver: ConcaveSemilinearDriver,
    config: FitConfig,
    actor_kind: ActorKind,
    train_t: np.ndarray,
    train_x: np.ndarray,
    train_r: np.ndarray,
    train_q: np.ndarray | None,
    valid_t: np.ndarray,
    valid_x: np.ndarray,
    valid_r: np.ndarray,
    valid_q: np.ndarray | None,
    *,
    seed: int,
    device: torch.device,
    initial: nn.Module | None,
) -> tuple[nn.Module, dict[str, Any], list[dict[str, Any]]]:
    seed_everything(seed)
    actor = make_actor(
        actor_kind, driver, width=config.width, depth=config.depth
    ).to(device)
    warm_started = initial is not None
    if initial is not None:
        if getattr(initial, "actor_kind", None) != actor_kind:
            raise ValueError("warm start actor kind mismatch")
        actor.load_state_dict(initial.state_dict())

    tt, tx, tr = _tensor(train_t, device), _tensor(train_x, device), _tensor(train_r, device)
    vt, vx, vr = _tensor(valid_t, device), _tensor(valid_x, device), _tensor(valid_r, device)
    tq = None if train_q is None else _tensor(train_q, device)
    vq = None if valid_q is None else _tensor(valid_q, device)
    r_scale = max(float(np.std(train_r)), 1e-2)
    q_scale = (
        None
        if train_q is None
        else max(float(np.sqrt(np.mean(np.square(train_q)))), 1e-3)
    )
    optimizer = torch.optim.AdamW(
        actor.parameters(), lr=config.learning_rate, weight_decay=1e-6
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=config.optimizer_steps
    )
    generator = torch.Generator(device=device).manual_seed(seed + 1)
    checkpoints = set(range(config.checkpoint_interval, config.optimizer_steps + 1,
                            config.checkpoint_interval))
    checkpoints.update((0, 1, config.optimizer_steps))

    best_loss = float("inf")
    best_state: dict[str, torch.Tensor] | None = None
    trace: list[dict[str, Any]] = []
    gradient_norm_window: list[float] = []
    clipped_window = 0

    def record(step: int) -> None:
        nonlocal best_loss, best_state, gradient_norm_window, clipped_window
        train_metrics = _fit_metrics(
            actor, driver, tt, tx, tr, tq, r_scale=r_scale, q_scale=q_scale
        )
        valid_metrics = _fit_metrics(
            actor, driver, vt, vx, vr, vq, r_scale=r_scale, q_scale=q_scale
        )
        row: dict[str, Any] = {"step": step, "learning_rate": optimizer.param_groups[0]["lr"]}
        row.update({f"train_{key}": value for key, value in train_metrics.items()})
        row.update({f"validation_{key}": value for key, value in valid_metrics.items()})
        for key in ("normalized_loss", "value_mse", "value_rmse"):
            row[f"generalization_gap_{key}"] = (
                float(valid_metrics[key]) - float(train_metrics[key])
            )
        if train_q is not None:
            for key in (
                "gradient_mse",
                "gradient_component_rmse",
                "gradient_vector_rmse",
            ):
                row[f"generalization_gap_{key}"] = (
                    float(valid_metrics[key]) - float(train_metrics[key])
                )
        else:
            row.update(
                {
                    "generalization_gap_gradient_mse": None,
                    "generalization_gap_gradient_component_rmse": None,
                    "generalization_gap_gradient_vector_rmse": None,
                }
            )
        row["gradient_norm_mean_since_checkpoint"] = (
            float(np.mean(gradient_norm_window)) if gradient_norm_window else None
        )
        row["gradient_norm_max_since_checkpoint"] = (
            float(np.max(gradient_norm_window)) if gradient_norm_window else None
        )
        row["gradient_clip_fraction_since_checkpoint"] = (
            clipped_window / len(gradient_norm_window) if gradient_norm_window else None
        )
        trace.append(row)
        loss = float(valid_metrics["normalized_loss"])
        if loss < best_loss:
            best_loss = loss
            best_state = {
                name: parameter.detach().cpu().clone()
                for name, parameter in actor.state_dict().items()
            }
        gradient_norm_window = []
        clipped_window = 0

    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    record(0)
    for step in range(1, config.optimizer_steps + 1):
        actor.train()
        index = torch.randint(
            len(tt),
            (min(config.batch_size, len(tt)),),
            generator=generator,
            device=device,
        )
        predicted_r, predicted_q = actor.predict(tt[index], tx[index], driver.T)
        value_mse = (predicted_r - tr[index]).square().mean()
        loss = value_mse / r_scale**2
        if tq is not None:
            assert predicted_q is not None and q_scale is not None
            gradient_mse = (predicted_q - tq[index]).square().mean()
            loss = loss + gradient_mse / q_scale**2
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        norm = float(
            torch.nn.utils.clip_grad_norm_(actor.parameters(), config.gradient_clip).item()
        )
        gradient_norm_window.append(norm)
        clipped_window += int(norm > config.gradient_clip)
        optimizer.step()
        scheduler.step()
        if step in checkpoints:
            record(step)

    if best_state is None:
        raise RuntimeError("no Actor checkpoint recorded")
    actor.load_state_dict(best_state)
    actor.eval()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak_gpu_bytes = int(torch.cuda.max_memory_allocated(device))
    else:
        peak_gpu_bytes = None
    selected = min(trace, key=lambda item: item["validation_normalized_loss"])
    training = {
        "config": asdict(config),
        "actor_kind": actor_kind,
        "parameters": count_parameters(actor),
        "warm_started": warm_started,
        "warm_start_step": 0,
        "value_label_scale": r_scale,
        "gradient_label_scale": q_scale,
        "train_value_label_mean": float(np.mean(train_r)),
        "train_value_label_std": float(np.std(train_r)),
        "validation_value_label_mean": float(np.mean(valid_r)),
        "validation_value_label_std": float(np.std(valid_r)),
        "train_gradient_label_component_rms": q_scale,
        "validation_gradient_label_component_rms": None
        if valid_q is None
        else float(np.sqrt(np.mean(np.square(valid_q)))),
        "selected_step": int(selected["step"]),
        "selected_by": "minimum normalized held-out noisy tangent-label loss",
        "best_validation_normalized_loss": float(best_loss),
        "selected_train_value_rmse": selected["train_value_rmse"],
        "selected_validation_value_rmse": selected["validation_value_rmse"],
        "selected_value_generalization_gap_rmse": selected[
            "generalization_gap_value_rmse"
        ],
        "selected_train_gradient_component_rmse": selected[
            "train_gradient_component_rmse"
        ],
        "selected_validation_gradient_component_rmse": selected[
            "validation_gradient_component_rmse"
        ],
        "selected_train_gradient_vector_rmse": selected[
            "train_gradient_vector_rmse"
        ],
        "selected_validation_gradient_vector_rmse": selected[
            "validation_gradient_vector_rmse"
        ],
        "selected_gradient_generalization_gap_vector_rmse": selected[
            "generalization_gap_gradient_vector_rmse"
        ],
        "selected_lower_saturation_fraction": selected[
            "validation_lower_saturation_fraction"
        ],
        "selected_upper_saturation_fraction": selected[
            "validation_upper_saturation_fraction"
        ],
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak_gpu_bytes,
        "used_exact_reference": False,
    }
    return actor, training, trace

