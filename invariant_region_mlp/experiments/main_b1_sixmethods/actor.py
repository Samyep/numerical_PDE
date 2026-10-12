"""Noisy-label ambient-input actor training."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn

from .common import load_config, seed_everything, stable_seed, training_points
from .equation import B1Equation
from .mlp_adapter import evaluate_mlp
from .neural_common import DenseNet, clone_state_dict, count_parameters, tx_input


class ActorNet(DenseNet):
    def __init__(self, d: int, width: int = 128, depth: int = 2) -> None:
        super().__init__(d + 1, d, width=width, depth=depth, activation="silu")
        self.d = d

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float = 0.25) -> torch.Tensor:
        return self(tx_input(t, x, T))


def actor_checkpoint_path(root: Path, equation: B1Equation) -> Path:
    return root / f"actor_d{equation.d}_k{equation.k}.pt"


def train_actor(
    equation: B1Equation,
    *,
    device: torch.device,
    checkpoint_path: Path,
    force: bool = False,
) -> tuple[ActorNet, dict[str, Any]]:
    config = load_config()["actor"]
    if checkpoint_path.is_file() and not force:
        payload = torch.load(checkpoint_path, map_location=device, weights_only=False)
        model = ActorNet(
            equation.d, config["hidden_width"], config["hidden_layers"]
        ).to(device)
        model.load_state_dict(payload["state_dict"])
        model.eval()
        return model, payload["metadata"]

    seed = stable_seed("B1", "actor", equation.d, equation.k)
    seed_everything(seed)
    train_t, train_x = training_points(
        equation.d, equation.k, "actor_labels_train", seed, config["training_points"]
    )
    valid_t, valid_x = training_points(
        equation.d, equation.k, "actor_labels_valid", seed, config["validation_points"]
    )
    label_started = time.perf_counter()
    train_state, train_cost = evaluate_mlp(
        equation,
        train_t,
        train_x,
        n=config["label_n"],
        M=config["label_M"],
        method="actor_label",
        repetition=0,
        device=device,
        chunk_size=512,
    )
    valid_state, valid_cost = evaluate_mlp(
        equation,
        valid_t,
        valid_x,
        n=config["label_n"],
        M=config["label_M"],
        method="actor_label",
        repetition=1,
        device=device,
        chunk_size=512,
    )
    label_seconds = time.perf_counter() - label_started

    model = ActorNet(
        equation.d, config["hidden_width"], config["hidden_layers"]
    ).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config["learning_rate"],
        weight_decay=config["weight_decay"],
    )
    tensors = {
        "train_t": torch.as_tensor(train_t, dtype=torch.float32, device=device),
        "train_x": torch.as_tensor(train_x, dtype=torch.float32, device=device),
        "train_y": torch.as_tensor(train_state[:, 1:], dtype=torch.float32, device=device),
        "valid_t": torch.as_tensor(valid_t, dtype=torch.float32, device=device),
        "valid_x": torch.as_tensor(valid_x, dtype=torch.float32, device=device),
        "valid_y": torch.as_tensor(valid_state[:, 1:], dtype=torch.float32, device=device),
    }
    generator = torch.Generator(device="cpu").manual_seed(seed)
    best_loss = float("inf")
    best_step = 0
    best_state: dict[str, torch.Tensor] | None = None
    patience = 0
    history: list[dict[str, float | int]] = []
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    train_started = time.perf_counter()
    for step in range(1, config["max_steps"] + 1):
        index = torch.randint(
            len(train_t), (config["batch_size"],), generator=generator
        ).to(device)
        prediction = model.predict(tensors["train_t"][index], tensors["train_x"][index])
        # Mean of per-example squared Euclidean errors, matching the specification.
        loss = torch.sum((prediction - tensors["train_y"][index]) ** 2, dim=-1).mean()
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        if step % config["validation_every"] == 0 or step == config["max_steps"]:
            model.eval()
            with torch.no_grad():
                valid_prediction = model.predict(tensors["valid_t"], tensors["valid_x"])
                valid_loss = torch.sum(
                    (valid_prediction - tensors["valid_y"]) ** 2, dim=-1
                ).mean().item()
            history.append({"step": step, "training_loss": loss.item(), "validation_loss": valid_loss})
            if valid_loss < best_loss * (1.0 - 1e-5):
                best_loss = valid_loss
                best_step = step
                best_state = clone_state_dict(model)
                patience = 0
            else:
                patience += 1
            model.train()
            if patience >= config["early_stopping_patience"]:
                break
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    train_seconds = time.perf_counter() - train_started
    if best_state is None:
        raise RuntimeError("actor training produced no checkpoint")
    model.load_state_dict(best_state)
    model.eval()

    # Analytic references appear only after training/checkpoint selection.
    with torch.no_grad():
        train_prediction = model.predict(tensors["train_t"], tensors["train_x"]).cpu().numpy()
        valid_prediction = model.predict(tensors["valid_t"], tensors["valid_x"]).cpu().numpy()
    exact_train = equation.exact_z(train_t, train_x)
    exact_valid = equation.exact_z(valid_t, valid_x)
    metadata: dict[str, Any] = {
        "method": "complete_actor_mlp",
        "d": equation.d,
        "k": equation.k,
        "seed": seed,
        "architecture": f"ambient (t,x)-MLP {config['hidden_layers']}x{config['hidden_width']} SiLU -> d",
        "parameter_count": count_parameters(model),
        "label_method": "pathwise terminal + level-matched baseline, original quadratic driver",
        "label_n": config["label_n"],
        "label_M": config["label_M"],
        "training_points": len(train_t),
        "validation_points": len(valid_t),
        "best_step": best_step,
        "steps_executed": step,
        "best_noisy_label_validation_loss": best_loss,
        "label_generation_seconds": label_seconds,
        "training_seconds": train_seconds,
        "total_actor_seconds": label_seconds + train_seconds,
        "peak_gpu_bytes": peak,
        "train_label_cost": train_cost,
        "valid_label_cost": valid_cost,
        "training_label_rmse_against_exact": float(
            np.sqrt(np.mean(np.sum((train_state[:, 1:] - exact_train) ** 2, axis=1)))
        ),
        "actor_validation_gradient_rmse": float(
            np.sqrt(np.mean(np.sum((valid_prediction - exact_valid) ** 2, axis=1)))
        ),
        "actor_validation_relative_error": float(
            np.sqrt(
                np.sum((valid_prediction - exact_valid) ** 2)
                / max(np.sum(exact_valid**2), np.finfo(np.float64).tiny)
            )
        ),
        "actor_validation_bias_norm": float(np.linalg.norm(np.mean(valid_prediction - exact_valid, axis=0))),
        "actor_validation_norm_mean": float(np.mean(np.linalg.norm(valid_prediction, axis=1))),
        "history": history,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
    }
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"state_dict": best_state, "metadata": metadata}, checkpoint_path)
    return model, metadata
