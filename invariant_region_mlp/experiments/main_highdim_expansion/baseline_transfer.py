"""Faithful DPI and PINN-HTE transfers from the pinned upstream repository.

Only the equation class, horizon, and ambient tensor dimension are changed.
The upstream data generator, PISGradNet, objectives, optimizer, sample counts,
epochs, and Picard depth remain authoritative.  ``smoke`` is an explicitly
non-result mode used by CI/preflight; ``official`` is the only result mode.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shutil
import time
from types import SimpleNamespace
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    JointConcaveDriver,
    JointConvexOriginalDriver,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.reference2d import (
    PeriodicReference2D,
)
from invariant_region_mlp.experiments.nonconvex_tangent_actor.drivers import (
    SineQuadraticDriver,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.references import (
    AllenRadialConfig,
    AllenRadialReference,
    hjb_log_reference,
)

from .benchmarks import (
    array_sha256,
    exact_state,
    make_benchmark,
    terminal_gradient_torch,
)
from .prepare import file_sha256 as sha256_file, load_test_bank
from .protocol import (
    DIMENSIONS, METHODS as EXPANSION_METHODS, PDES, RESULTS_DIR,
    SOURCE_COMMITS, TRAINING_SEED, stable_seed, write_json_atomic,
)


PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
EXTERNAL_ROOT = SOURCE_ROOT.parent / "external" / "DeepPicardIteration"
HJB_DPI_CONFIG = EXTERNAL_ROOT / "scripts/hjb/base_100d_T1.0_w0.1_0.yaml"
HJB_PINN_CONFIG = EXTERNAL_ROOT / "scripts/hjb/pinn_100d_T1.0_v16.yaml"
SEMILINEAR_DPI_CONFIG = EXTERNAL_ROOT / "scripts/burgers/base_100d_T1.0_w1.0_0.yaml"
SEMILINEAR_PINN_CONFIG = EXTERNAL_ROOT / "scripts/burgers/pinn_100d_T1.0_v16_beta10.0.yaml"
RAW_DIR = RESULTS_DIR / "raw"
PREDICTION_DIR = RESULTS_DIR / "predictions"
CHECKPOINT_DIR = RESULTS_DIR / "official_baseline_checkpoints"
WORK_DIR = RESULTS_DIR / "official_baseline_work"
INPUT_DIR = RESULTS_DIR / "inputs"
REFERENCE_DIR = RESULTS_DIR / "references"

METHODS = ("dpi", "pinn_hte")
PROFILES = ("auto", "hjb_pis", "semilinear_paper")


def _resolve_profile(pde: str, profile: str) -> str:
    if pde not in PDES:
        raise KeyError(pde)
    if profile == "auto":
        return "semilinear_paper"
    if profile not in PROFILES:
        raise KeyError(profile)
    if profile == "hjb_pis":
        raise ValueError("PISGradNet is source-specific to OU/GMM and is not transferred")
    return profile


def _dirs() -> None:
    for path in (RAW_DIR, PREDICTION_DIR, CHECKPOINT_DIR, WORK_DIR):
        path.mkdir(parents=True, exist_ok=True)


def register_equations() -> dict[str, str]:
    """Register source-compatible subclasses without altering upstream code."""
    import picard.equations as eqs

    class ExpansionEquation(eqs.SimpleDiffusionEquationWithZ):
        expansion_pde = ""

        def __init__(self, nx: int, T: float, alpha: float = 2.0, **_: Any):
            super().__init__(nx=nx, T=T, alpha=alpha)
            self._benchmark = make_benchmark(self.expansion_pde, nx)

        def g(self, x: torch.Tensor) -> torch.Tensor:
            return self._benchmark.terminal(x).unsqueeze(-1)

        def g_x(self, x: torch.Tensor) -> torch.Tensor:
            return terminal_gradient_torch(self._benchmark, x)

        def fff(self, t, x, y, z):
            del t, x
            return self._benchmark.driver(y.squeeze(-1), z).unsqueeze(-1)

        def exact_solution(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
            state = exact_state(
                self._benchmark,
                t.detach().cpu().double().numpy().reshape(-1),
                x.detach().cpu().double().numpy(),
            )
            return torch.as_tensor(state[:, :1], dtype=x.dtype, device=x.device)

        def u_x(self, t: torch.Tensor, x: torch.Tensor, retain_graph: bool = False):
            del retain_graph
            state = exact_state(
                self._benchmark,
                t.detach().cpu().double().numpy().reshape(-1),
                x.detach().cpu().double().numpy(),
            )
            gradient = state[:, 1:] / self._benchmark.sigma
            return torch.as_tensor(gradient, dtype=x.dtype, device=x.device)

        def u_u_x(self, t: torch.Tensor, x: torch.Tensor):
            return self.exact_solution(t, x), self.u_x(t, x)

    names: dict[str, str] = {}
    for key in PDES:
        cls = type(f"Expansion{key.title().replace('_', '')}Equation", (ExpansionEquation,), {"expansion_pde": key})
        cls.__module__ = eqs.__name__
        setattr(eqs, cls.__name__, cls)
        names[key] = cls.__name__
    return names


def _prepare_ou_workdir(d: int) -> Path:
    target = WORK_DIR / "ou_inputs" / f"d{d}"
    target.mkdir(parents=True, exist_ok=True)
    with np.load(INPUT_DIR / f"ou_gmm_d{d}.npz") as data:
        mean = torch.as_tensor(data["mean"], dtype=torch.float64)
        var = torch.stack([2.0 * torch.eye(d, dtype=torch.float64) for _ in range(5)])
        pi = torch.as_tensor(data["pi"], dtype=torch.float64)
    torch.save(mean, target / f"mean_{d}d_ms=1.0_vs=2.0_5.pt")
    torch.save(var, target / f"var_{d}d_ms=1.0_vs=2.0_5.pt")
    torch.save(pi, target / f"pi_{d}d_ms=1.0_vs=2.0_5.pt")
    return target


def _horizon_alpha(pde: str) -> tuple[float, float]:
    if pde == "b6_fullrank":
        return 0.25, 2.0
    if pde == "dense_fourier":
        return 0.5, 2.0
    raise KeyError(pde)


def _load_test_bank(pde: str, d: int) -> tuple[np.ndarray, np.ndarray, Path, str, str]:
    t, x, _truth_u, _truth_z, path = load_test_bank(pde, d)
    return t, x, path, array_sha256(t, x), sha256_file(path)


def _truth(equation, t: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    state = exact_state(equation._benchmark, t, x)
    return state[:, 0], state[:, 1:] / equation._benchmark.sigma


def _rng_state() -> dict[str, Any]:
    result = {
        "python": random.getstate(), "numpy": np.random.get_state(),
        "torch": torch.random.get_rng_state(),
    }
    if torch.cuda.is_available():
        result["cuda"] = torch.cuda.get_rng_state_all()
    return result


def _restore_rng(state: dict[str, Any]) -> None:
    random.setstate(state["python"])
    np.random.set_state(state["numpy"])
    torch.random.set_rng_state(state["torch"])
    if torch.cuda.is_available() and "cuda" in state:
        torch.cuda.set_rng_state_all(state["cuda"])


def _patch_upstream(capture: Path) -> None:
    import lightning.pytorch as pl
    from lightning.pytorch.callbacks import RichProgressBar
    from picard import PicardRunner
    import picard.solution as solutions

    original_pinn_train = solutions.PINNSolution.train

    def pinn_train_and_capture(self):
        result = original_pinn_train(self)
        torch.save(self.dnn.state_dict(), capture)
        return result

    solutions.PINNSolution.train = pinn_train_and_capture
    # PyTorch 2.6 changed torch.load's default to weights_only=True.  The
    # pinned upstream code predates that change and reloads its own trusted
    # Lightning checkpoints without passing weights_only.  Restore the old
    # loading behavior without editing the pinned source checkout.
    os.environ.setdefault("TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD", "1")
    original_callbacks = PicardRunner.get_reporting_callbacks

    def callbacks_without_rich(self):
        callbacks, kwargs = original_callbacks(self)
        callbacks = [item for item in callbacks if not isinstance(item, RichProgressBar)]
        kwargs["enable_progress_bar"] = False
        return callbacks, kwargs

    PicardRunner.get_reporting_callbacks = callbacks_without_rich
    # Avoid non-numerical wandb finalization contacting an external service.
    os.environ.setdefault("WANDB_MODE", "offline")
    pl.seed_everything(torch.initial_seed() % (2**32), workers=True)


def _config(
    method: str, pde: str, d: int, seed: int, mode: str,
    class_name: str, profile: str,
):
    from picard import load_cfg

    resolved = _resolve_profile(pde, profile)
    if resolved == "hjb_pis":
        config_path = HJB_DPI_CONFIG if method == "dpi" else HJB_PINN_CONFIG
    else:
        config_path = SEMILINEAR_DPI_CONFIG if method == "dpi" else SEMILINEAR_PINN_CONFIG
    run_id = f"{method}_{pde}_d{d}_{resolved}_{mode}_s{seed}"
    exp_dir = WORK_DIR / run_id
    T, alpha = _horizon_alpha(pde)
    overrides = [
        "NAME", str(exp_dir), "FORCE", "True",
        "EQUATION.cls", class_name,
        "EQUATION.kwargs.nx", str(d),
        "EQUATION.kwargs.T", str(T),
        "EQUATION.kwargs.alpha", str(alpha),
    ]
    if resolved == "semilinear_paper" and method == "dpi":
        # Table 1 of arXiv:2409.08526v2 fixes K=80 for the semilinear
        # experiment. The README's gradient-supervised YAML says K=100; the
        # paper budget is authoritative for this benchmark protocol.
        overrides += ["PICARD.N", "80"]
    if mode == "smoke":
        overrides += [
            "PICARD.N", "1", "DATA.DATA_SIZE", "32",
            "DATA.NEW_SAMPLING", "False", "DATA.N_BUFFER", "1",
            "DATA.kwargs.n_estimate_terminal", "8",
            "DATA.kwargs.n_estimate_integral", "8",
            "TRAIN.N_EPOCHS", "1", "TRAIN.BATCH_SIZE", "16",
            "NETWORK.NEURONS", "[16,16]",
            "NETWORK.ACTIVATIONS", "['ELU','ELU']", "EVAL.FREQ", "1",
        ]
    elif mode == "one_iteration":
        if method != "dpi":
            raise ValueError("one_iteration applies only to DPI")
        overrides += ["PICARD.N", "1"]
    elif mode != "official":
        raise KeyError(mode)
    old = Path.cwd()
    os.chdir(config_path.parent)
    try:
        cfg = load_cfg(str(config_path), overrides)
    finally:
        os.chdir(old)
    cfg.freeze()
    return cfg, config_path, overrides, run_id, exp_dir, resolved


def _reconstruct_model(cfg, equation, checkpoint: Path, device: torch.device):
    from picard import PicardRunner
    from picard.solution_jac import PicardSolutionGradientWrapper

    class Runner(SimpleNamespace):
        def get_solution_plain(self):
            cls = PicardRunner.get_solution_cls(self.cfg)
            return cls.construct_solution(self)

    runner = Runner(equation=equation, cfg=cfg)
    model = PicardSolutionGradientWrapper.construct_solution(runner).to(device)
    payload = torch.load(checkpoint, map_location=device, weights_only=False)
    state = payload["state_dict"] if isinstance(payload, dict) and "state_dict" in payload else payload
    model.load_state_dict(state)
    return model.eval()


def execute(
    method: str, pde: str, d: int, seed: int, mode: str,
    profile: str = "auto",
) -> dict[str, Any]:
    _dirs()
    if mode == "official" and seed != 0:
        raise ValueError("single-seed benchmark protocol permits only seed 0")
    if EXTERNAL_ROOT.joinpath(".git").exists():
        import subprocess
        revision = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=EXTERNAL_ROOT, text=True
        ).strip()
        if revision != SOURCE_COMMITS[method]:
            raise RuntimeError(f"pinned source mismatch: {revision}")
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    names = register_equations()
    cfg, config_path, overrides, run_id, exp_dir, resolved_profile = _config(
        method, pde, d, seed, mode, names[pde], profile
    )
    output = RAW_DIR / f"{run_id}.json"
    if output.exists():
        existing = json.loads(output.read_text())
        expected_status = "complete" if mode == "official" else f"{mode}_nonresult"
        if existing.get("status") == expected_status:
            return existing
    capture = CHECKPOINT_DIR / f"{run_id}_final.pt"
    _patch_upstream(capture)
    from picard import PicardRunner

    cwd = WORK_DIR
    old = Path.cwd()
    os.chdir(cwd)
    started = time.perf_counter()
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats(); torch.cuda.synchronize()
    try:
        runner = PicardRunner(cfg)
        original_run_one = runner.run_one
        rng_path = CHECKPOINT_DIR / f"{run_id}_rng.pt"

        def run_one_with_rng_checkpoint():
            outcome = original_run_one()
            torch.save({"iteration": runner.i, "rng": _rng_state()}, rng_path)
            return outcome

        runner.run_one = run_one_with_rng_checkpoint
        runner.run()
    finally:
        os.chdir(old)
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    training_seconds = time.perf_counter() - started
    checkpoint = exp_dir / f"model_{cfg.PICARD.N}.ckpt" if method == "dpi" else capture
    if not checkpoint.exists():
        raise FileNotFoundError(checkpoint)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    equation = runner.equation.to(device=device)
    model = _reconstruct_model(cfg, equation, checkpoint, device)
    t, x, bank_path, bank_hash, bank_file_hash = _load_test_bank(pde, d)
    truth, truth_gradient = _truth(equation, t, x)
    reference_metadata = {
        "kind": "analytic frozen reference exposed only after training",
        "transform": "Hopf-Cole" if pde == "b6_fullrank" else "sinh",
    }
    predictions, gradients = [], []
    inference_started = time.perf_counter()
    for begin in range(0, len(t), 128):
        end = min(begin + 128, len(t))
        tt = torch.as_tensor(t[begin:end, None], dtype=torch.float64, device=device)
        xx = torch.as_tensor(x[begin:end], dtype=torch.float64, device=device).requires_grad_(True)
        value = model(torch.cat((tt, xx), dim=-1))[:, 0]
        gradient = torch.autograd.grad(value.sum(), xx)[0]
        predictions.append(value.detach().cpu().numpy())
        gradients.append(gradient.detach().cpu().numpy())
    if device.type == "cuda":
        torch.cuda.synchronize()
    inference_seconds = time.perf_counter() - inference_started
    prediction = np.concatenate(predictions)
    gradient = np.concatenate(gradients)
    error = prediction - truth
    native_t = np.asarray([0.0], dtype=np.float64)
    native_x = np.zeros((1, d), dtype=np.float64)
    native_tt = torch.zeros((1, 1), dtype=torch.float64, device=device)
    native_xx = torch.zeros((1, d), dtype=torch.float64, device=device)
    with torch.no_grad():
        native_prediction = float(model(torch.cat((native_tt, native_xx), dim=1))[0, 0])
    native_truth = float(_truth(equation, native_t, native_x)[0][0])
    prediction_path = PREDICTION_DIR / f"{run_id}.npz"
    np.savez_compressed(
        prediction_path, t=t, x=x, prediction=prediction, gradient=gradient,
        truth=truth, truth_gradient=np.empty((0, d)) if truth_gradient is None else truth_gradient,
    )
    result: dict[str, Any] = {
        "run_id": run_id, "method": method, "pde": pde, "dimension": d,
        "seed": seed, "mode": mode,
        "status": "complete" if mode == "official" else f"{mode}_nonresult",
        "classification": "faithful_transfer",
        "benchmark_profile": resolved_profile,
        "profile_role": (
            "primary" if resolved_profile == _resolve_profile(pde, "auto")
            else "diagnostic_nonprimary"
        ),
        "profile_selection_rule": (
            "paper Section 5.2 generic semilinear profile; no hard terminal embedding"
        ),
        "official_source_commit": SOURCE_COMMITS[method],
        "official_config": str(config_path), "expanded_config": str(exp_dir / "config.yaml"),
        "overrides": overrides, "checkpoint": str(checkpoint),
        "test_bank": str(bank_path), "test_bank_sha256": bank_hash,
        "test_bank_array_sha256": bank_hash,
        "test_bank_file_sha256": bank_file_hash,
        "reference_array_sha256": array_sha256(np.asarray(truth)),
        "reference_metadata": reference_metadata,
        "test_states": len(t), "value_rmse": float(np.sqrt(np.mean(error**2))),
        "value_relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
        "value_bias": float(np.mean(error)), "value_mae": float(np.mean(np.abs(error))),
        "training_seconds": training_seconds, "inference_seconds": inference_seconds,
        "total_seconds": training_seconds + inference_seconds,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated()) if torch.cuda.is_available() else None,
        "parameter_count": sum(item.numel() for item in model.parameters()),
        "prediction_file": str(prediction_path),
        "prediction_file_sha256": sha256_file(prediction_path),
        "native": {
            "point": "t=0,x=0", "prediction": native_prediction,
            "reference": native_truth,
            "absolute_error": abs(native_prediction - native_truth),
            "full_training_cost_charged": True,
        },
        "ground_truth_used_for_training_or_selection": False,
        "checkpoint_load_compatibility": (
            "TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1 restores the pinned source's "
            "pre-PyTorch-2.6 trusted Lightning checkpoint semantics"
        ),
        "terminal_aware_feature_added_by_us": False,
        "mcv_added_by_us": False,
        "oracle_projection": False,
        "transfer_note": "PDE/dimension tensors changed; baseline architecture and source budget unchanged",
    }
    if truth_gradient is not None:
        delta = gradient - truth_gradient
        result.update({
            "gradient_component_rmse": float(np.sqrt(np.mean(delta**2))),
            "gradient_vector_rmse": float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
            "gradient_relative_l2": float(np.linalg.norm(delta) / max(np.linalg.norm(truth_gradient), 1e-300)),
        })
    write_json_atomic(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--pde", choices=PDES, required=True)
    parser.add_argument("--dimension", type=int, choices=DIMENSIONS, required=True)
    parser.add_argument("--seed", type=int, default=TRAINING_SEED)
    parser.add_argument("--mode", choices=("smoke", "one_iteration", "official"), default="official")
    parser.add_argument("--profile", choices=PROFILES, default="auto")
    args = parser.parse_args()
    result = execute(
        args.method, args.pde, args.dimension, args.seed, args.mode, args.profile,
    )
    print(json.dumps({
        key: result[key] for key in (
            "run_id", "status", "value_rmse", "value_relative_l2",
            "training_seconds", "inference_seconds",
        )
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
