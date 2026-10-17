"""A frozen original-form neural PINN adaptation for the new PDEs.

This comparator uses a plain full-ambient MLP, interior PDE residuals, a soft
terminal penalty, and Hutchinson Laplacian estimates.  It does not use an
Actor, tangent labels, martingale control variates, active coordinates, or the
terminal-aware parameterization of our method.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from .common import (
    CHECKPOINT_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    RESULTS_DIR,
    atomic_json,
    atomic_torch_save,
    ensure_dirs,
    environment_payload,
    metrics,
    seed_everything,
    stable_seed,
    synchronize,
)
from .configuration import ExperimentSpec


@dataclass(frozen=True)
class PINNConfig:
    width: int = 256
    depth: int = 4
    optimizer_steps: int = 5000
    interior_batch: int = 256
    terminal_batch: int = 256
    hutchinson_probes: int = 2
    learning_rate: float = 5e-4
    terminal_weight: float = 10.0
    validation_points: int = 256
    validation_probes: int = 8
    checkpoint_interval: int = 250
    x_std: float = 1.0


class PlainPINN(nn.Module):
    def __init__(self, d: int, width: int, depth: int) -> None:
        super().__init__()
        layers: list[nn.Module] = []
        previous = d + 1
        for _ in range(depth):
            layer = nn.Linear(previous, width)
            nn.init.xavier_uniform_(layer.weight)
            nn.init.zeros_(layer.bias)
            layers.extend((layer, nn.Tanh()))
            previous = width
        final = nn.Linear(previous, 1)
        nn.init.xavier_uniform_(final.weight)
        nn.init.zeros_(final.bias)
        layers.append(final)
        self.network = nn.Sequential(*layers)

    def forward(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        inputs = torch.cat((t.reshape(-1, 1) / T, x), dim=-1)
        return self.network(inputs).reshape(-1)


def _rademacher(
    shape: tuple[int, ...], generator: torch.Generator, device: torch.device
) -> torch.Tensor:
    values = torch.randint(
        0, 2, shape, generator=generator, device=device, dtype=torch.int64
    )
    return values.to(torch.float32).mul_(2.0).sub_(1.0)


def pde_residual(
    model: PlainPINN,
    driver,
    t: torch.Tensor,
    x: torch.Tensor,
    *,
    probes: int,
    generator: torch.Generator,
) -> torch.Tensor:
    t_var = t.detach().clone().requires_grad_(True)
    x_var = x.detach().clone().requires_grad_(True)
    value = model(t_var, x_var, driver.T)
    value_t, gradient = torch.autograd.grad(
        value.sum(), (t_var, x_var), create_graph=True
    )
    laplace = torch.zeros_like(value)
    for _ in range(probes):
        epsilon = _rademacher(x_var.shape, generator, x_var.device)
        directional = (gradient * epsilon).sum()
        hessian_vector = torch.autograd.grad(
            directional, x_var, create_graph=True
        )[0]
        laplace = laplace + (hessian_vector * epsilon).sum(dim=1)
    laplace = laplace / probes
    z = driver.sigma * gradient
    return value_t + laplace + driver.driver(value, z)


def _validation(
    model: PlainPINN,
    driver,
    config: PINNConfig,
    t: torch.Tensor,
    x: torch.Tensor,
    terminal_x: torch.Tensor,
    *,
    seed: int,
    device: torch.device,
) -> dict[str, float]:
    model.eval()
    generator = torch.Generator(device=device).manual_seed(seed)
    with torch.enable_grad():
        residual = pde_residual(
            model,
            driver,
            t,
            x,
            probes=config.validation_probes,
            generator=generator,
        )
    with torch.no_grad():
        terminal_error = (
            model(
                torch.full(
                    (len(terminal_x),), driver.T, dtype=torch.float32, device=device
                ),
                terminal_x,
                driver.T,
            )
            - driver.terminal(terminal_x)
        )
    return {
        "residual_rmse": float(residual.detach().double().square().mean().sqrt().item()),
        "terminal_rmse": float(
            terminal_error.detach().double().square().mean().sqrt().item()
        ),
    }


def train(
    spec: ExperimentSpec, baseline_seed: int, device: torch.device
) -> tuple[PlainPINN, dict[str, Any]]:
    config = PINNConfig()
    driver = spec.driver()
    seed_everything(baseline_seed)
    model = PlainPINN(driver.d, config.width, config.depth).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=config.optimizer_steps
    )
    generator = torch.Generator(device=device).manual_seed(
        stable_seed("nonconvex_pinn_train", spec.run_id, baseline_seed)
    )
    valid_generator = torch.Generator(device=device).manual_seed(
        stable_seed("nonconvex_pinn_validation_bank", spec.run_id, baseline_seed)
    )
    valid_t = torch.rand(
        config.validation_points, generator=valid_generator, device=device
    ) * driver.T
    valid_x = torch.randn(
        config.validation_points,
        driver.d,
        generator=valid_generator,
        device=device,
    )
    valid_terminal_x = torch.randn(
        config.validation_points,
        driver.d,
        generator=valid_generator,
        device=device,
    )
    trace: list[dict[str, float]] = []
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for step in range(1, config.optimizer_steps + 1):
        model.train()
        t = torch.rand(
            config.interior_batch, generator=generator, device=device
        ) * driver.T
        x = config.x_std * torch.randn(
            config.interior_batch, driver.d, generator=generator, device=device
        )
        terminal_x = config.x_std * torch.randn(
            config.terminal_batch, driver.d, generator=generator, device=device
        )
        residual = pde_residual(
            model,
            driver,
            t,
            x,
            probes=config.hutchinson_probes,
            generator=generator,
        )
        terminal_prediction = model(
            torch.full(
                (config.terminal_batch,),
                driver.T,
                dtype=torch.float32,
                device=device,
            ),
            terminal_x,
            driver.T,
        )
        terminal_error = terminal_prediction - driver.terminal(terminal_x)
        residual_loss = residual.square().mean()
        terminal_loss = terminal_error.square().mean()
        loss = residual_loss + config.terminal_weight * terminal_loss
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        optimizer.step()
        scheduler.step()
        if step == 1 or step % config.checkpoint_interval == 0:
            validation = _validation(
                model,
                driver,
                config,
                valid_t,
                valid_x,
                valid_terminal_x,
                seed=stable_seed(
                    "nonconvex_pinn_validation_probes", spec.run_id, step
                ),
                device=device,
            )
            trace.append(
                {
                    "step": step,
                    "training_residual_rmse": float(
                        residual_loss.detach().sqrt().item()
                    ),
                    "training_terminal_rmse": float(
                        terminal_loss.detach().sqrt().item()
                    ),
                    "validation_residual_rmse": validation["residual_rmse"],
                    "validation_terminal_rmse": validation["terminal_rmse"],
                    "learning_rate": optimizer.param_groups[0]["lr"],
                }
            )
            print(
                f"PINN {spec.run_id} step={step} "
                f"res={validation['residual_rmse']:.4g} "
                f"terminal={validation['terminal_rmse']:.4g}",
                flush=True,
            )
    synchronize(device)
    elapsed = time.perf_counter() - started
    checkpoint = CHECKPOINT_DIR / f"pinn_{spec.run_id}_s{baseline_seed}.pt"
    atomic_torch_save(
        checkpoint,
        {
            "state_dict": {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            },
            "config": asdict(config),
            "spec": asdict(spec),
            "baseline_seed": baseline_seed,
        },
    )
    return model, {
        "config": asdict(config),
        "trace": trace,
        "training_seconds": elapsed,
        "checkpoint": str(checkpoint),
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        "checkpoint_selection": "final optimizer step; Ground Truth never consulted",
    }


@torch.no_grad()
def evaluate(
    spec: ExperimentSpec,
    model: PlainPINN,
    training: dict[str, Any],
    *,
    baseline_seed: int,
    actor_seed: int,
    device: torch.device,
) -> dict[str, Any]:
    actor_spec = ExperimentSpec(
        family=spec.family,
        d=spec.d,
        lam=spec.lam,
        eta=spec.eta,
        seed=actor_seed,
        budget="main",
        estimator="mcv",
    )
    bank_path = PREDICTION_DIR / f"post_freeze_{actor_spec.run_id}.npz"
    if not bank_path.exists():
        raise FileNotFoundError(bank_path)
    with np.load(bank_path) as archive:
        t, x, truth = archive["test_t"], archive["test_x"], archive["truth"]
    synchronize(device)
    started = time.perf_counter()
    parts = []
    for begin in range(0, len(t), 1024):
        end = min(begin + 1024, len(t))
        parts.append(
            model(
                torch.as_tensor(t[begin:end], dtype=torch.float32, device=device),
                torch.as_tensor(x[begin:end], dtype=torch.float32, device=device),
                spec.driver().T,
            )
            .double()
            .cpu()
            .numpy()
        )
    synchronize(device)
    inference_seconds = time.perf_counter() - started
    prediction = np.concatenate(parts)
    prediction_path = PREDICTION_DIR / f"pinn_{actor_spec.run_id}_s{baseline_seed}.npz"
    np.savez_compressed(prediction_path, prediction=prediction)
    result = {
        "method": "hutchinson_pinn_original_form_adaptation",
        "spec": asdict(spec),
        "actor_test_bank_run_id": actor_spec.run_id,
        "baseline_seed": baseline_seed,
        "metrics": metrics(prediction, truth),
        "training": training,
        "training_seconds": training["training_seconds"],
        "inference_seconds": inference_seconds,
        "total_seconds": training["training_seconds"] + inference_seconds,
        "plain_full_ambient_mlp": True,
        "soft_terminal_penalty": True,
        "terminal_aware_parameterization": False,
        "tangent_actor": False,
        "mcv": False,
        "oracle_projection": False,
        "ground_truth_used_for_training_or_selection": False,
        "original_algorithm_adaptation_not_published_benchmark_reproduction": True,
        "prediction_path": str(prediction_path),
    }
    atomic_json(
        RAW_DIR / f"baseline_pinn_{actor_spec.run_id}_s{baseline_seed}.json",
        result,
    )
    print(
        json.dumps(
            {
                "method": result["method"],
                "rmse": result["metrics"]["rmse"],
                "relative_l2": result["metrics"]["relative_l2"],
                "seconds": result["total_seconds"],
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=("primary", "coupled"), default="primary")
    parser.add_argument("--d", type=int, default=100)
    parser.add_argument("--lam", type=float, default=0.2)
    parser.add_argument("--eta", type=float, default=1.0)
    parser.add_argument("--actor-seed", type=int, default=0)
    parser.add_argument("--baseline-seed", type=int, default=0)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    spec = ExperimentSpec(
        family=args.family,
        d=args.d,
        lam=args.lam,
        eta=args.eta if args.family == "coupled" else 0.0,
        seed=args.actor_seed,
        budget="main",
        estimator="mcv",
    )
    ensure_dirs()
    atomic_json(
        RAW_DIR / f"environment_pinn_{spec.run_id}_s{args.baseline_seed}.json",
        environment_payload("pinn_baseline"),
    )
    model, training = train(spec, args.baseline_seed, torch.device(args.device))
    evaluate(
        spec,
        model,
        training,
        baseline_seed=args.baseline_seed,
        actor_seed=args.actor_seed,
        device=torch.device(args.device),
    )


if __name__ == "__main__":
    main()
