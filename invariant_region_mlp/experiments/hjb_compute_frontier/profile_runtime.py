"""Final-scale GH200 microbenchmarks used only for resource forecasting."""

from __future__ import annotations

import json
import time

import numpy as np
import torch

from .actor_method import evaluate_policy_cv, policy_gradient_labels
from .baselines import (
    AmortizedDeepBSDE,
    HardTerminalPINN,
    TerminalAwareDPI,
    _dpi_dataset,
    full_pinn_residual_mse,
    sdgd_objective,
)
from .common import RESULTS_DIR, atomic_json, require_cuda, seed_everything, stable_seed
from .equations import make_equation, torch_terminal
from .models import ActorNet


def _timed(device: torch.device, callback) -> tuple[float, int, object]:
    torch.cuda.synchronize(device)
    torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    result = callback()
    torch.cuda.synchronize(device)
    return (
        time.perf_counter() - started,
        int(torch.cuda.max_memory_allocated(device)),
        result,
    )


def profile() -> dict[str, object]:
    device = require_cuda()
    equation = make_equation("b1")
    seed = stable_seed("hjb_compute_frontier", "final_scale_runtime_profile")
    seed_everything(seed)
    generator = torch.Generator(device=device).manual_seed(seed)
    output: dict[str, object] = {
        "gpu_name": torch.cuda.get_device_name(device),
        "seed": seed,
        "equation": "b1_d400_k40_representative",
    }

    pinn = HardTerminalPINN(equation.d, width=512, depth=4).to(device)
    pinn_optimizer = torch.optim.Adam(pinn.parameters(), lr=1e-3)

    def pinn_steps() -> None:
        for _ in range(10):
            t = torch.rand(128, device=device, generator=generator) * equation.T
            x = torch.randn((128, equation.d), device=device, generator=generator)
            loss, _ = sdgd_objective(
                pinn,
                equation,
                t,
                x,
                coordinate_budget=32,
                generator=generator,
            )
            pinn_optimizer.zero_grad(set_to_none=True)
            loss.backward()
            pinn_optimizer.step()

    pinn_seconds, pinn_peak, _ = _timed(device, pinn_steps)
    validation_t = torch.rand(8, device=device, generator=generator) * equation.T
    validation_x = torch.randn((8, equation.d), device=device, generator=generator)
    validation_seconds, validation_peak, validation_loss = _timed(
        device,
        lambda: full_pinn_residual_mse(
            pinn, equation, validation_t, validation_x
        ),
    )
    output["pinn_sdgd"] = {
        "measured_steps": 10,
        "step_seconds": pinn_seconds / 10,
        "projected_10000_step_seconds": pinn_seconds * 1000,
        "full_laplacian_validation_seconds": validation_seconds,
        "validation_loss_finiteness": bool(np.isfinite(validation_loss)),
        "peak_gpu_bytes": max(pinn_peak, validation_peak),
    }
    del pinn, pinn_optimizer
    torch.cuda.empty_cache()

    bsde = AmortizedDeepBSDE(equation.d, width=256, depth=3).to(device)
    bsde_optimizer = torch.optim.Adam(bsde.parameters(), lr=1e-3)

    def bsde_steps() -> None:
        for _ in range(10):
            x0 = torch.randn((128, equation.d), device=device, generator=generator)
            normal = torch.randn((128, 32, equation.d), device=device, generator=generator)
            terminal_y, terminal_x = bsde.rollout(equation, x0, normal)
            loss = torch.nn.functional.huber_loss(
                terminal_y, torch_terminal(equation, terminal_x), delta=10.0
            )
            bsde_optimizer.zero_grad(set_to_none=True)
            loss.backward()
            bsde_optimizer.step()

    bsde_seconds, bsde_peak, _ = _timed(device, bsde_steps)
    output["deep_bsde"] = {
        "measured_steps": 10,
        "step_seconds": bsde_seconds / 10,
        "projected_10000_step_seconds": bsde_seconds * 1000,
        "peak_gpu_bytes": bsde_peak,
    }
    del bsde, bsde_optimizer
    torch.cuda.empty_cache()

    previous = TerminalAwareDPI(equation.d, width=512, depth=4).to(device).eval()
    dpi_label_seconds, dpi_label_peak, _ = _timed(
        device,
        lambda: _dpi_dataset(
            equation,
            previous,
            points=32,
            paths=512,
            generator=generator,
            device=device,
        ),
    )
    dpi_model = TerminalAwareDPI(equation.d, width=512, depth=4).to(device)
    dpi_optimizer = torch.optim.Adam(dpi_model.parameters(), lr=1e-3)
    train_t = torch.rand(512, device=device, generator=generator) * equation.T
    train_x = torch.randn((512, equation.d), device=device, generator=generator)
    train_u = torch.randn(512, device=device, generator=generator)
    train_z = torch.randn((512, equation.d), device=device, generator=generator)

    def dpi_updates() -> None:
        for _ in range(8):
            predicted_u, predicted_z = dpi_model.value_z(equation, train_t, train_x)
            loss = torch.mean((predicted_u - train_u).square()) + 0.1 * torch.mean(
                (predicted_z - train_z).square()
            )
            dpi_optimizer.zero_grad(set_to_none=True)
            loss.backward()
            dpi_optimizer.step()

    dpi_train_seconds, dpi_train_peak, _ = _timed(device, dpi_updates)
    projected_round = dpi_label_seconds * (4096 + 256) / 32 + dpi_train_seconds * 16
    output["dpi"] = {
        "label_points_measured": 32,
        "paths_per_label": 512,
        "label_seconds_per_point": dpi_label_seconds / 32,
        "supervised_seconds_per_epoch": dpi_train_seconds,
        "projected_seconds_per_round": projected_round,
        "projected_20_round_seconds": 20 * projected_round,
        "peak_gpu_bytes": max(dpi_label_peak, dpi_train_peak),
    }
    del previous, dpi_model, dpi_optimizer
    torch.cuda.empty_cache()

    actor = ActorNet(equation.d, width=384, depth=4).to(device)
    label_t = np.linspace(0.0, equation.T * 0.95, 16)
    label_x = np.random.default_rng(seed).standard_normal((16, equation.d))
    actor_label_seconds, actor_label_peak, label_result = _timed(
        device,
        lambda: policy_gradient_labels(
            equation,
            actor,
            label_t,
            label_x,
            paths=512,
            time_steps=32,
            seed=seed + 1,
            device=device,
            parent_batch=8,
            path_chunk=64,
        ),
    )
    actor_optimizer = torch.optim.AdamW(actor.parameters(), lr=1e-4)
    batch_t = torch.rand(1024, device=device, generator=generator) * equation.T
    batch_x = torch.randn((1024, equation.d), device=device, generator=generator)
    batch_y = torch.randn((1024, equation.d), device=device, generator=generator)

    def actor_updates() -> None:
        for _ in range(50):
            prediction = actor.predict(batch_t, batch_x, equation.T)
            loss = torch.mean(torch.sum((prediction - batch_y).square(), dim=-1)) / equation.d
            actor_optimizer.zero_grad(set_to_none=True)
            loss.backward()
            actor_optimizer.step()

    actor_train_seconds, actor_train_peak, _ = _timed(device, actor_updates)
    eval_seconds, eval_peak, _ = _timed(
        device,
        lambda: evaluate_policy_cv(
            equation,
            actor,
            np.zeros(16),
            label_x,
            paths=256,
            time_steps=64,
            seed=seed + 2,
            device=device,
            parent_batch=16,
            path_chunk=64,
        ),
    )
    projected_phase = (
        actor_label_seconds * (8192 + 512) / 16
        + actor_train_seconds * 6000 / 50
    )
    output["actor"] = {
        "policy_label_points_measured": 16,
        "policy_label_seconds_per_point": actor_label_seconds / 16,
        "label_outputs_finite": bool(
            np.all(np.isfinite(label_result[0])) and np.all(np.isfinite(label_result[1]))
        ),
        "optimizer_seconds_per_step": actor_train_seconds / 50,
        "projected_8192_label_6000_step_phase_seconds": projected_phase,
        "evaluation_seconds_per_point_256x64": eval_seconds / 16,
        "projected_evaluation_2048_seconds": eval_seconds * 128,
        "peak_gpu_bytes": max(actor_label_peak, actor_train_peak, eval_peak),
    }
    output["projection_note"] = (
        "Linear extrapolations are conservative microbenchmark forecasts; "
        "final decisions use complete pilot accounting."
    )
    return output


def main() -> None:
    result = profile()
    atomic_json(RESULTS_DIR / "runtime_profile.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
