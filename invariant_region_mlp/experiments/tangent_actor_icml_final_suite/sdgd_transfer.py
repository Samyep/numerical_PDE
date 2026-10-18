"""Faithful PINN+SDGD transfers using the pinned HJB-Log algorithm.

The source architecture, hard terminal ansatz, two-sample stochastic trace
gradient, adversarial collocation loop, Adam schedule, and 10,000 epochs are
preserved.  PDE coefficients and terminal functions are the only scientific
changes.  The loop is checkpointed at epoch boundaries for exact restart.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import random
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    JointConcaveDriver, JointConvexOriginalDriver,
)
from invariant_region_mlp.experiments.nonconvex_tangent_actor.drivers import SineQuadraticDriver
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.references import hjb_log_reference

from .benchmarks import AllenCahnBenchmark, HJBLogBenchmark, MLPFailureBenchmark, OUGMMBenchmark
from .deep_bsde_transfer import DriverBenchmark
from .prepare_suite import write_json
from .protocol import RESULTS_DIR, SOURCE_COMMITS, stable_seed
from .run_joint_actor import reference as periodic_reference
from .run_standard_actor import _allen_reference


FAMILIES = (
    "hjb_log", "ou_gmm", "joint_concave", "joint_convex",
    "nonconvex_coupled", "mlp_failure", "allen_cahn",
)
INPUT_DIR = RESULTS_DIR / "inputs"
RAW_DIR = RESULTS_DIR / "raw"
PREDICTION_DIR = RESULTS_DIR / "predictions"
CHECKPOINT_DIR = RESULTS_DIR / "official_baseline_checkpoints"


@dataclass(frozen=True)
class Config:
    epochs: int = 10_000
    width: int = 1024
    hidden_layers: int = 3
    collocation_points: int = 100
    pgd_steps: int = 20
    pgd_step_size: float = 0.05
    dimension_samples: int = 10
    learning_rate: float = 1e-3
    reporting_points: int = 20_000


class HardTerminalNet(nn.Module):
    def __init__(self, benchmark, config: Config):
        super().__init__()
        dimensions = [benchmark.d + 1] + [config.width] * config.hidden_layers + [1]
        layers: list[nn.Module] = []
        for index in range(len(dimensions) - 1):
            layers.append(nn.Linear(dimensions[index], dimensions[index + 1]))
            if index != len(dimensions) - 2:
                layers.append(nn.Tanh())
        self.nn = nn.Sequential(*layers)
        self.benchmark = benchmark

    def forward(self, tx: torch.Tensor) -> torch.Tensor:
        x, t = tx[:, :-1], tx[:, -1:]
        return self.nn(tx) * (self.benchmark.T - t) + self.benchmark.terminal(x).reshape(-1, 1)


def make_benchmark(family: str, d: int):
    if family == "hjb_log":
        return HJBLogBenchmark(d)
    if family == "ou_gmm":
        return OUGMMBenchmark.from_npz(INPUT_DIR / f"ou_gmm_d{d}.npz")
    if family == "allen_cahn":
        return AllenCahnBenchmark(d)
    if family == "mlp_failure":
        return MLPFailureBenchmark(d)
    return DriverBenchmark(family, d)


def _driver(benchmark, y: torch.Tensor, grad: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
    z = benchmark.sigma * grad
    value = benchmark.driver(y, z)
    # OUGMMBenchmark places the known linear term in its base drift for the
    # Tangent solver; here it belongs explicitly in the PINN residual.
    drift = benchmark.base_drift(torch.zeros_like(y), x)
    return value + (drift * grad).sum(-1)


def _truth(family: str, benchmark, t: np.ndarray, x: np.ndarray, *, smoke: bool = False):
    if family == "hjb_log":
        return hjb_log_reference(t, x, order=64 if smoke else 128)
    if family == "ou_gmm":
        return benchmark.exact_value(
            torch.as_tensor(t, dtype=torch.float64), torch.as_tensor(x, dtype=torch.float64)
        ).numpy()
    if family == "mlp_failure":
        return benchmark.exact_value_numpy(t, x)
    if family == "allen_cahn":
        return _allen_reference(benchmark.d, t, x, smoke=smoke)[0]
    reference, _ = periodic_reference(family)
    coordinates = x @ benchmark.impl.B
    value = reference.value(t, coordinates)
    return -value if family == "joint_convex" else value


def _bank(family: str, benchmark, count: int, seed: int, *, smoke: bool):
    rng = np.random.default_rng(seed)
    t = rng.uniform(0.0, np.nextafter(benchmark.T, 0.0), count).astype(np.float64)
    x = rng.standard_normal((count, benchmark.d)).astype(np.float64)
    truth = _truth(family, benchmark, t, x, smoke=smoke)
    return t, x, truth


def _common_bank(family: str, d: int):
    audit = json.loads((RESULTS_DIR / "PREFLIGHT_AUDIT.json").read_text())
    item = audit["test_banks"][f"{family}_d{d}"]
    source_root = Path(__file__).resolve().parents[3]
    path = source_root / item["path"]
    with np.load(path) as data:
        if "test_indices" in data:
            index = data["test_indices"]
            t, x = data["t"][index], data["x"][index]
        elif family == "nonconvex_coupled" and "test_t" in data:
            t, x = data["test_t"], data["test_x"]
        else:
            t, x = data["t"], data["x"]
    return np.asarray(t), np.asarray(x), path, item["sha256"]


def _rng_state():
    result = {"python": random.getstate(), "numpy": np.random.get_state(), "torch": torch.get_rng_state()}
    if torch.cuda.is_available():
        result["cuda"] = torch.cuda.get_rng_state_all()
    return result


def _restore_rng(state):
    random.setstate(state["python"]); np.random.set_state(state["numpy"])
    torch.set_rng_state(state["torch"].cpu())
    if torch.cuda.is_available() and "cuda" in state:
        torch.cuda.set_rng_state_all([item.cpu() for item in state["cuda"]])


class Solver:
    def __init__(self, family: str, benchmark, config: Config, device: torch.device, seed: int, smoke: bool):
        self.family, self.benchmark, self.config, self.device = family, benchmark, config, device
        self.model = HardTerminalNet(benchmark, config).to(device=device, dtype=torch.float32)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=config.learning_rate)
        self.scheduler = torch.optim.lr_scheduler.LambdaLR(
            self.optimizer, lambda epoch: 1.0 - epoch / config.epochs
        )
        self.xf = self._resample()
        count = min(config.reporting_points, 16) if smoke else config.reporting_points
        t, x, truth = _bank(
            family, benchmark, count, stable_seed("sdgd_reporting_bank", family, benchmark.d, seed),
            smoke=smoke,
        )
        self.reporting_tx = torch.as_tensor(
            np.concatenate((x, t[:, None]), axis=1), dtype=torch.float32, device=device,
        )
        self.reporting_truth = torch.as_tensor(truth[:, None], dtype=torch.float32, device=device)

    def _resample(self):
        return torch.cat((
            torch.randn((self.config.collocation_points, self.benchmark.d), device=self.device),
            self.benchmark.T * torch.rand((self.config.collocation_points, 1), device=self.device),
        ), dim=-1).requires_grad_()

    def _base(self):
        x = self.xf[:, :-1]
        t = self.xf[:, -1:].requires_grad_()
        x = x.requires_grad_()
        value = self.model(torch.cat((x, t), dim=1))
        u_t = torch.autograd.grad(value.sum(), t, create_graph=True)[0][:, 0]
        grad = torch.autograd.grad(value.sum(), x, create_graph=True)[0]
        base = u_t + _driver(self.benchmark, value[:, 0], grad, x)
        return x, grad, base

    def stochastic_loss(self, samples: int):
        x, grad, base = self._base()
        d = self.benchmark.d
        count = min(samples, d)
        coefficient = 0.5 * self.benchmark.sigma**2
        first = np.random.choice(d, count, replace=False)
        residual = base.detach().clone()
        for index in first:
            diagonal = torch.autograd.grad(grad[:, index].sum(), x, create_graph=True)[0][:, index]
            residual = residual + coefficient * d / count * diagonal.detach()
        second = np.random.choice(d, count, replace=False)
        gradient_estimator = base
        for index in second:
            diagonal = torch.autograd.grad(grad[:, index].sum(), x, create_graph=True)[0][:, index]
            gradient_estimator = gradient_estimator + coefficient * d / count * diagonal
        loss = (2.0 * residual.detach() * gradient_estimator).mean()
        return loss, residual.square().mean()

    def pgd(self):
        for _ in range(self.config.pgd_steps):
            self.xf.requires_grad_()
            loss, _ = self.stochastic_loss(self.config.dimension_samples)
            gradient = torch.autograd.grad(loss, self.xf)[0]
            self.xf = self.xf.detach() + self.config.pgd_step_size * torch.sign(gradient.detach())
            self.xf[:, -1] = torch.clamp(self.xf[:, -1], 0.0, self.benchmark.T)

    @torch.no_grad()
    def reporting_l2(self):
        error = self.reporting_truth - self.model(self.reporting_tx)
        return float(torch.linalg.vector_norm(error) / torch.linalg.vector_norm(self.reporting_truth))

    def step(self):
        self.pgd()
        loss, reported = self.stochastic_loss(self.config.dimension_samples)
        self.optimizer.zero_grad(set_to_none=True); loss.backward(); self.optimizer.step(); self.scheduler.step()
        self.xf = self._resample()
        return float(reported.detach())


def execute(family: str, d: int, seed: int, device: torch.device, *, smoke: bool = False):
    for path in (RAW_DIR, PREDICTION_DIR, CHECKPOINT_DIR):
        path.mkdir(parents=True, exist_ok=True)
    config = Config(**({
        "epochs": 2, "width": 16, "hidden_layers": 2,
        "collocation_points": 8, "pgd_steps": 2, "pgd_step_size": 0.05,
        "dimension_samples": 2, "learning_rate": 1e-3, "reporting_points": 16,
    } if smoke else {}))
    tag = "smoke" if smoke else "official"
    run_id = f"sdgd_pinn_{family}_d{d}_{tag}_s{seed}"
    output = RAW_DIR / f"{run_id}.json"
    if output.exists() and json.loads(output.read_text()).get("status") in {"complete", "smoke_nonresult"}:
        return json.loads(output.read_text())
    full_seed = stable_seed("icml_sdgd_transfer", family, d, seed)
    random.seed(full_seed); np.random.seed(full_seed); torch.manual_seed(full_seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(full_seed)
    benchmark = make_benchmark(family, d)
    solver = Solver(family, benchmark, config, device, seed, smoke)
    resume = CHECKPOINT_DIR / f"{run_id}_continuation.pt"
    start_epoch, history, cumulative = 0, [], 0.0
    if resume.exists():
        state = torch.load(resume, map_location=device, weights_only=False)
        if state["config"] != asdict(config):
            raise RuntimeError("resume configuration mismatch")
        solver.model.load_state_dict(state["model"]); solver.optimizer.load_state_dict(state["optimizer"])
        solver.scheduler.load_state_dict(state["scheduler"]); solver.xf = state["xf"].to(device).requires_grad_()
        _restore_rng(state["rng"]); start_epoch = int(state["next_epoch"])
        history, cumulative = state["history"], float(state["cumulative_seconds"])
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(); torch.cuda.synchronize()
    started = time.perf_counter()
    for epoch in range(start_epoch, config.epochs):
        loss = solver.step()
        l2 = solver.reporting_l2()
        if epoch % 100 == 0 or smoke:
            print(f"{run_id} epoch={epoch} loss={loss:.6g} l2={l2:.6g}", flush=True)
        history.append({"epoch": epoch, "reported_residual_mse": loss, "reporting_relative_l2": l2})
        if (epoch + 1) % 100 == 0 or epoch + 1 == config.epochs:
            if device.type == "cuda":
                torch.cuda.synchronize()
            temporary = resume.with_suffix(".tmp")
            torch.save({
                "model": solver.model.state_dict(), "optimizer": solver.optimizer.state_dict(),
                "scheduler": solver.scheduler.state_dict(), "xf": solver.xf.detach().cpu(),
                "rng": _rng_state(), "next_epoch": epoch + 1, "history": history,
                "cumulative_seconds": cumulative + time.perf_counter() - started,
                "config": asdict(config),
            }, temporary); temporary.replace(resume)
    if device.type == "cuda":
        torch.cuda.synchronize()
    elapsed = time.perf_counter() - started
    checkpoint = CHECKPOINT_DIR / f"{run_id}.pt"
    torch.save(solver.model.state_dict(), checkpoint)

    t, x, bank_path, bank_hash = _common_bank(family, d)
    if smoke:
        t, x = t[:16], x[:16]
    truth = _truth(family, benchmark, t, x, smoke=smoke)
    predictions, gradients = [], []
    inference_started = time.perf_counter()
    for begin in range(0, len(t), 256):
        end = min(begin + 256, len(t))
        xx = torch.as_tensor(x[begin:end], dtype=torch.float32, device=device).requires_grad_(True)
        tt = torch.as_tensor(t[begin:end, None], dtype=torch.float32, device=device)
        value = solver.model(torch.cat((xx, tt), dim=1))[:, 0]
        grad = torch.autograd.grad(value.sum(), xx)[0]
        predictions.append(value.detach().cpu().double().numpy())
        gradients.append((benchmark.sigma * grad).detach().cpu().double().numpy())
    if device.type == "cuda":
        torch.cuda.synchronize()
    inference_seconds = time.perf_counter() - inference_started
    prediction = np.concatenate(predictions); gradient = np.concatenate(gradients)
    error = prediction - truth
    prediction_path = PREDICTION_DIR / f"{run_id}.npz"
    np.savez_compressed(prediction_path, t=t, x=x, truth=truth, prediction=prediction, gradient=gradient)
    result = {
        "run_id": run_id, "method": "sdgd_pinn_official_transfer",
        "family": family, "dimension": d, "seed": seed,
        "status": "smoke_nonresult" if smoke else "complete",
        "classification": "native_official" if family == "hjb_log" and d == 100 else "faithful_transfer",
        "official_source_commit": SOURCE_COMMITS["sdgd_pinn"],
        "configuration": asdict(config), "source_configuration": "Section_5.3/HJB_Log.py",
        "dimension_sampling_rule": "min(10,d); equals the full trace at d=10",
        "terminal_aware_feature_added_by_us": False, "mcv_added_by_us": False,
        "hard_terminal_embedding_is_source_native": True,
        "history": history, "checkpoint": str(checkpoint),
        "training_seconds_this_invocation": elapsed,
        "cumulative_training_seconds": cumulative + elapsed,
        "inference_seconds": inference_seconds,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)) if device.type == "cuda" else None,
        "parameter_count": sum(item.numel() for item in solver.model.parameters()),
        "test_bank": str(bank_path), "test_bank_sha256": bank_hash, "states": len(t),
        "value_rmse": float(np.sqrt(np.mean(error**2))),
        "value_relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
        "value_bias": float(np.mean(error)), "prediction_file": str(prediction_path),
        "test_truth_used_only_for_source_native_reporting_and_post_freeze_evaluation": True,
        "test_truth_used_in_training_loss": False,
    }
    write_json(output, result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=FAMILIES, required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation required")
    result = execute(args.family, args.dimension, args.seed, device, smoke=args.smoke)
    print(json.dumps({
        "run_id": result["run_id"], "status": result["status"],
        "value_rmse": result["value_rmse"], "training_seconds": result["training_seconds_this_invocation"],
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
