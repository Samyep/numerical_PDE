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
import re
import shutil
import tempfile
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

from .benchmarks import MLPFailureBenchmark, OUGMMBenchmark
from .paired_bank_repair import array_sha256
from .prepare_suite import sha256_file, write_json
from .protocol import RESULTS_DIR, SOURCE_COMMITS, stable_seed


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

PDES = (
    "hjb_log", "ou_gmm", "joint_concave", "joint_convex",
    "nonconvex_coupled", "mlp_failure", "allen_cahn",
)
METHODS = ("dpi", "pinn_hte")
PROFILES = ("auto", "hjb_pis", "semilinear_paper")


def _resolve_profile(pde: str, profile: str) -> str:
    if profile == "auto":
        # The PISGradNet is a source-native specialization for the OU/GMM
        # score-model PDE. It is not transferred to our other benchmarks.
        return "hjb_pis" if pde == "ou_gmm" else "semilinear_paper"
    if profile not in PROFILES:
        raise KeyError(profile)
    return profile


def _dirs() -> None:
    for path in (RAW_DIR, PREDICTION_DIR, CHECKPOINT_DIR, WORK_DIR):
        path.mkdir(parents=True, exist_ok=True)


def _basis_driver(pde: str, d: int):
    if pde == "joint_concave":
        return JointConcaveDriver(d=d, k=2, T=0.5)
    if pde == "joint_convex":
        return JointConvexOriginalDriver(d=d, k=2, T=0.5)
    if pde == "nonconvex_coupled":
        return SineQuadraticDriver(d=d, T=0.5, lam=0.2, eta=1.0)
    raise KeyError(pde)


def _load_periodic_reference(pde: str) -> tuple[PeriodicReference2D, float]:
    if pde == "nonconvex_coupled":
        path = (
            RESULTS_DIR.parent / "tangent_actor_accuracy_frontier/references/"
            "sine_quadratic_coupled_d100_k2_lam0p2_eta1p0_refined.npz"
        )
        sign = 1.0
    else:
        path = REFERENCE_DIR / f"{pde}_k2_n128.npz"
        sign = -1.0 if pde == "joint_convex" else 1.0
    if not path.exists():
        # Reference creation belongs to the independent Tangent-Actor runner;
        # invoking it here only creates an independent numerical artifact and
        # never passes reference values into training.
        from .run_joint_actor import reference
        reference(pde)
    with np.load(path) as data:
        reference_value = PeriodicReference2D(
            name=str(data["name"]), times=data["times"], space=data["space"],
            values=data["values"], derivative_1=data["derivative_1"],
            derivative_2=data["derivative_2"], dt=float(data["dt"]),
            store_dt=float(data["store_dt"]), n=int(data["n"]),
        )
    return reference_value, sign


class _AllenDiagnosticSurface:
    """Independent coarse surface used only by upstream logging callbacks.

    Final metrics use the separately refined radial reference.  Precomputing
    this surface avoids re-solving the radial PDE at every upstream callback.
    """

    _cache: dict[int, "_AllenDiagnosticSurface"] = {}

    def __new__(cls, d: int):
        if d in cls._cache:
            return cls._cache[d]
        value = super().__new__(cls)
        cls._cache[d] = value
        return value

    def __init__(self, d: int):
        if hasattr(self, "values"):
            return
        self.d = int(d)
        self.ds = 0.2
        self.dt = 0.0005
        self.s_max = max(100.0, 1.25 * self.d)
        solver = AllenRadialReference(AllenRadialConfig(
            s_max=self.s_max, ds=self.ds, dt=self.dt,
            dimension=self.d, horizon=0.3,
        ))
        steps = int(round(0.3 / self.dt))
        values = np.empty((steps + 1, len(solver.s)), dtype=np.float64)
        values[0] = 0.5 / (1.0 + 0.2 * solver.s)
        for index in range(1, steps + 1):
            values[index] = solver._step(values[index - 1])
        self.s = solver.s
        self.values = values
        self.derivative = np.gradient(values, self.ds, axis=1)

    def evaluate(self, t: np.ndarray, radius2: np.ndarray, derivative: bool = False) -> np.ndarray:
        tau = np.clip(0.3 - np.asarray(t), 0.0, 0.3)
        s = np.clip(np.asarray(radius2), 0.0, self.s_max)
        ft = tau / self.dt
        it = np.minimum(np.floor(ft).astype(int), len(self.values) - 2)
        wt = ft - it
        table = self.derivative if derivative else self.values
        low = np.asarray([np.interp(s[i], self.s, table[it[i]]) for i in range(len(s))])
        high = np.asarray([np.interp(s[i], self.s, table[it[i] + 1]) for i in range(len(s))])
        return low * (1.0 - wt) + high * wt


def register_equations() -> dict[str, str]:
    """Register source-compatible equation subclasses in ``picard.equations``."""
    import picard.equations as eqs

    class SuiteSimpleEquation(eqs.SimpleDiffusionEquationWithZ):
        suite_pde = ""

        def __init__(self, nx: int, T: float, alpha: float = 2.0, **_: Any):
            super().__init__(nx=nx, T=T, alpha=alpha)
            self._driver = _basis_driver(self.suite_pde, nx) if self.suite_pde in {
                "joint_concave", "joint_convex", "nonconvex_coupled"
            } else None
            if self._driver is not None:
                self.B = torch.as_tensor(self._driver.B, dtype=torch.get_default_dtype())

        def _basis(self, x: torch.Tensor) -> torch.Tensor:
            return self.B.to(dtype=x.dtype, device=x.device)

        def g(self, x: torch.Tensor) -> torch.Tensor:
            pde = self.suite_pde
            if pde == "hjb_log":
                return torch.log((1.0 + x.square().sum(-1, keepdim=True)) / 2.0)
            if pde == "mlp_failure":
                return x[..., :1].abs()
            return self._driver.terminal(x).unsqueeze(-1)

        def g_x(self, x: torch.Tensor) -> torch.Tensor:
            pde = self.suite_pde
            if pde == "hjb_log":
                return 2.0 * x / (1.0 + x.square().sum(-1, keepdim=True))
            if pde == "mlp_failure":
                result = torch.zeros_like(x)
                result[:, 0] = torch.sign(x[:, 0])
                return result
            return self._driver.terminal_gradient_torch(x)

        def fff(self, t, x, y, z):
            del t, x
            value = y.squeeze(-1)
            if self.suite_pde == "hjb_log":
                result = -0.5 * z.square().sum(-1)
            elif self.suite_pde == "mlp_failure":
                # alpha=1 here, hence z=grad u.  The zero vector uses the
                # minimum-norm subgradient through PyTorch's norm backward.
                result = torch.linalg.vector_norm(z[..., 1:], dim=-1)
            else:
                result = self._driver.driver(value, z)
            return result.unsqueeze(-1)

        def exact_solution(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
            pde = self.suite_pde
            times = t.detach().cpu().double().numpy().reshape(-1)
            states = x.detach().cpu().double().numpy()
            if pde == "hjb_log":
                result = hjb_log_reference(times, states)
            elif pde == "mlp_failure":
                result = MLPFailureBenchmark(self.nx).exact_value_numpy(times, states)
            else:
                reference, sign = _load_periodic_reference(pde)
                result = sign * reference.value(times, states @ self._driver.B)
            return torch.as_tensor(result[:, None], dtype=x.dtype, device=x.device)

        def u_x(self, t: torch.Tensor, x: torch.Tensor, retain_graph: bool = False) -> torch.Tensor:
            del retain_graph
            pde = self.suite_pde
            times = t.detach().cpu().double().numpy().reshape(-1)
            states = x.detach().cpu().double().numpy()
            if pde == "hjb_log":
                _, z = hjb_log_reference(times, states, return_gradient=True)
                gradient = z / math.sqrt(2.0)
            elif pde == "mlp_failure":
                gradient = MLPFailureBenchmark(self.nx).exact_z_numpy(times, states)
            else:
                reference, sign = _load_periodic_reference(pde)
                reduced = sign * reference.gradient(times, states @ self._driver.B)
                gradient = reduced @ self._driver.B.T
            return torch.as_tensor(gradient, dtype=x.dtype, device=x.device)

        def u_u_x(self, t: torch.Tensor, x: torch.Tensor):
            return self.exact_solution(t, x), self.u_x(t, x)

    class ICMLHJBLog(SuiteSimpleEquation):
        suite_pde = "hjb_log"

    class ICMLJointConcave(SuiteSimpleEquation):
        suite_pde = "joint_concave"

    class ICMLJointConvex(SuiteSimpleEquation):
        suite_pde = "joint_convex"

    class ICMLNonconvexCoupled(SuiteSimpleEquation):
        suite_pde = "nonconvex_coupled"

    class ICMLMLPFailure(SuiteSimpleEquation):
        suite_pde = "mlp_failure"

    class ICMLAllenCahn(eqs.SimpleDiffusionEquationWithoutZ):
        def __init__(self, nx: int, T: float = 0.3, alpha: float = 2.0, **_: Any):
            super().__init__(nx=nx, T=T, alpha=alpha)

        def sample_x0(self, n: int):
            return torch.zeros((n, self.nx), dtype=torch.get_default_dtype(), device=self.device)

        def g(self, x):
            return 0.5 / (1.0 + 0.2 * x.square().sum(-1, keepdim=True))

        def g_x(self, x):
            denominator = 1.0 + 0.2 * x.square().sum(-1, keepdim=True)
            return -0.2 * x / denominator.square()

        def f(self, t, x, y):
            del t, x
            return y - y.pow(3)

        def exact_solution(self, t, x):
            surface = _AllenDiagnosticSurface(self.nx)
            times = t.detach().cpu().double().numpy().reshape(-1)
            states = x.detach().cpu().double().numpy()
            result = surface.evaluate(times, np.square(states).sum(1))
            return torch.as_tensor(result[:, None], dtype=x.dtype, device=x.device)

        def u_x(self, t, x, retain_graph: bool = False):
            del retain_graph
            surface = _AllenDiagnosticSurface(self.nx)
            times = t.detach().cpu().double().numpy().reshape(-1)
            states = x.detach().cpu().double().numpy()
            derivative = surface.evaluate(times, np.square(states).sum(1), derivative=True)
            result = 2.0 * states * derivative[:, None]
            return torch.as_tensor(result, dtype=x.dtype, device=x.device)

        def u_u_x(self, t, x):
            return self.exact_solution(t, x), self.u_x(t, x)

    classes = {
        "hjb_log": ICMLHJBLog,
        "joint_concave": ICMLJointConcave,
        "joint_convex": ICMLJointConvex,
        "nonconvex_coupled": ICMLNonconvexCoupled,
        "mlp_failure": ICMLMLPFailure,
        "allen_cahn": ICMLAllenCahn,
    }
    names: dict[str, str] = {"ou_gmm": "OUProcessEquation"}
    for key, cls in classes.items():
        # Stable public names are required when Lightning serializes equation
        # objects into checkpoints.
        cls.__name__ = f"ICML{key.title().replace('_', '')}Equation"
        cls.__qualname__ = cls.__name__
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
    if pde in {"joint_concave", "joint_convex", "nonconvex_coupled"}:
        return 0.5, 2.0
    if pde == "allen_cahn":
        return 0.3, 2.0
    if pde == "mlp_failure":
        return 1.0, 1.0
    if pde == "hjb_log":
        return 1.0, 2.0
    return 1.0, 1.0


def _load_test_bank(pde: str, d: int) -> tuple[np.ndarray, np.ndarray, Path, str, str]:
    audit = json.loads((RESULTS_DIR / "PREFLIGHT_AUDIT.json").read_text())
    item = audit["test_banks"][f"{pde}_d{d}"]
    path = SOURCE_ROOT / item["path"]
    with np.load(path) as data:
        if "test_indices" in data:
            index = data["test_indices"]
            t, x = data["t"][index], data["x"][index]
        elif pde == "nonconvex_coupled" and "test_t" in data:
            t, x = data["test_t"], data["test_x"]
        else:
            t, x = data["t"], data["x"]
    t = np.asarray(t).reshape(-1)
    x = np.asarray(x)
    return t, x, path, array_sha256(t, x), item["sha256"]


def _truth(equation, t: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray | None]:
    # Use the independent post-freeze references exposed by the equation;
    # none of these calls are reachable from the optimizer or model selector.
    local_t = torch.as_tensor(t[:, None], dtype=torch.float64)
    local_x = torch.as_tensor(x, dtype=torch.float64)
    original_device = equation.device
    equation.to(device="cpu")
    try:
        value = equation.exact_solution(local_t, local_x).detach().numpy().reshape(-1)
        try:
            gradient = equation.u_x(local_t, local_x).detach().numpy()
        except (AttributeError, NotImplementedError):
            gradient = None
    finally:
        equation.to(device=original_device)
    return value, gradient


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


def _atomic_torch_save(payload: dict[str, Any], path: Path) -> None:
    """Commit a round-boundary state without exposing a partial file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        torch.save(payload, temporary)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _load_dpi_resume_state(
    path: Path, *, run_id: str, exp_dir: Path, cfg,
) -> dict[str, Any] | None:
    if not path.exists():
        return None
    payload = torch.load(path, map_location="cpu", weights_only=False)
    iteration = int(payload.get("iteration", 0))
    if iteration <= 0 or iteration > int(cfg.PICARD.N):
        raise RuntimeError(f"invalid DPI resume iteration {iteration} in {path}")
    if payload.get("run_id", run_id) != run_id:
        raise RuntimeError(f"DPI resume run mismatch in {path}")
    checkpoint = exp_dir / f"model_{iteration}.ckpt"
    if not checkpoint.exists():
        raise RuntimeError(
            f"DPI resume state references missing round checkpoint {checkpoint}"
        )
    if (
        iteration > int(cfg.DATA.MEMORY.REUSE)
        and "previous_dataset_size_args" not in payload
    ):
        raise RuntimeError(
            "DPI resume state predates exact dataset-size restoration; "
            "refusing an inexact continuation"
        )
    payload["checkpoint"] = checkpoint
    return payload


def _resume_metadata(
    state: dict[str, Any] | None,
) -> tuple[int, float, int]:
    """Return counters for either a fresh invocation or an exact resume."""
    if state is None:
        return 0, 0.0, 0
    return (
        int(state["iteration"]),
        float(state.get("training_seconds", 0.0)),
        int(state.get("resume_count", 0)) + 1,
    )


def _archive_uncommitted_dpi_rounds(exp_dir: Path, committed: int) -> list[str]:
    """Move artifacts newer than the atomic resume marker out of the run."""
    candidates: list[Path] = []
    for item in exp_dir.iterdir():
        match = re.fullmatch(r"model_(\d+)\.ckpt", item.name)
        if match and int(match.group(1)) > committed:
            candidates.append(item)
            continue
        match = re.fullmatch(r"pl_iter_(\d+)", item.name)
        if match and int(match.group(1)) > committed:
            candidates.append(item)
    if not candidates:
        return []
    archive = exp_dir / "resume_orphans" / f"from_round_{committed}_{time.time_ns()}"
    archive.mkdir(parents=True, exist_ok=False)
    moved = []
    for item in candidates:
        destination = archive / item.name
        item.rename(destination)
        moved.append(str(destination))
    return moved


def _construct_resumed_dpi_runner(PicardRunner, cfg, state: dict[str, Any]):
    """Construct upstream state without letting FORCE erase prior rounds."""
    exp_dir = Path(cfg.NAME)
    _archive_uncommitted_dpi_rounds(exp_dir, int(state["iteration"]))
    with tempfile.TemporaryDirectory(prefix="dpi_resume_bootstrap_", dir=WORK_DIR) as root:
        bootstrap_cfg = cfg.clone()
        bootstrap_cfg.defrost()
        bootstrap_cfg.NAME = str(Path(root) / "run")
        bootstrap_cfg.FORCE = True
        bootstrap_cfg.freeze()
        runner = PicardRunner(bootstrap_cfg)
        runner.writer.close()

    import tensorboardX

    runner.cfg = cfg
    runner.exp_dir = exp_dir
    runner.writer = tensorboardX.SummaryWriter(
        log_dir=os.path.join(cfg.LOGGING.TENSORBOARD_DIR, cfg.NAME)
    )
    runner.i = int(state["iteration"])
    runner.previous_dataset_size_args = state.get("previous_dataset_size_args")
    previous = runner.get_solution()
    checkpoint = torch.load(state["checkpoint"], weights_only=False)
    previous.load_state_dict(checkpoint["state_dict"])
    runner.u_current = previous
    runner.u_history = [previous]
    _restore_rng(state["rng"])
    return runner


def _run_picard_with_round_checkpoints(
    runner,
    *,
    state_path: Path,
    run_id: str,
    invocation_started: float,
    prior_training_seconds: float = 0.0,
    resume_count: int = 0,
    max_rounds: int | None = None,
) -> None:
    if runner.cfg.PICARD.FORMULA == "TwoLayer":
        raise RuntimeError("round-boundary resume is not implemented for TwoLayer DPI")
    original_run_one = runner.run_one
    rounds_run = 0
    while runner.i < runner.N and (max_rounds is None or rounds_run < max_rounds):
        outcome = original_run_one()
        if not outcome:
            print("Interrupted...Other processes is terminating...  ")
            break
        runner.u_history.append(runner.u_current)
        rounds_run += 1
        _atomic_torch_save({
            "schema_version": 2,
            "run_id": run_id,
            "iteration": int(runner.i),
            "rng": _rng_state(),
            "previous_dataset_size_args": runner.previous_dataset_size_args,
            "training_seconds": (
                float(prior_training_seconds)
                + float(time.perf_counter() - invocation_started)
            ),
            "resume_count": int(resume_count),
        }, state_path)


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

    cwd = _prepare_ou_workdir(d) if pde == "ou_gmm" else WORK_DIR
    old = Path.cwd()
    os.chdir(cwd)
    started = time.perf_counter()
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats(); torch.cuda.synchronize()
    rng_path = CHECKPOINT_DIR / f"{run_id}_rng.pt"
    resume_state = (
        _load_dpi_resume_state(
            rng_path, run_id=run_id, exp_dir=exp_dir, cfg=cfg,
        )
        if method == "dpi" else None
    )
    resumed_from_iteration, prior_training_seconds, resume_count = _resume_metadata(
        resume_state
    )
    try:
        runner = (
            _construct_resumed_dpi_runner(PicardRunner, cfg, resume_state)
            if resume_state else PicardRunner(cfg)
        )
        if method == "dpi":
            _run_picard_with_round_checkpoints(
                runner,
                state_path=rng_path,
                run_id=run_id,
                invocation_started=started,
                prior_training_seconds=prior_training_seconds,
                resume_count=resume_count,
            )
        else:
            runner.run()
    finally:
        os.chdir(old)
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    training_seconds = prior_training_seconds + time.perf_counter() - started
    checkpoint = exp_dir / f"model_{cfg.PICARD.N}.ckpt" if method == "dpi" else capture
    if not checkpoint.exists():
        raise FileNotFoundError(checkpoint)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    equation = runner.equation.to(device=device)
    model = _reconstruct_model(cfg, equation, checkpoint, device)
    t, x, bank_path, bank_hash, bank_file_hash = _load_test_bank(pde, d)
    if pde == "allen_cahn" and mode == "official":
        from .run_standard_actor import _allen_reference
        truth, reference_metadata = _allen_reference(d, t, x, smoke=False)
        truth_gradient = None
    else:
        truth, truth_gradient = _truth(equation, t, x)
        reference_metadata = {
            "kind": "independent suite reference exposed only after training",
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
    if pde == "allen_cahn" and mode == "official":
        native_truth = float(_allen_reference(d, native_t, native_x, smoke=False)[0][0])
    else:
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
        "classification": (
            "native_official" if (
                pde == "ou_gmm" and d == 100 and resolved_profile == "hjb_pis"
            ) else ("adapted" if pde == "mlp_failure" else "faithful_transfer")
        ),
        "benchmark_profile": resolved_profile,
        "profile_role": (
            "primary" if resolved_profile == _resolve_profile(pde, "auto")
            else "diagnostic_nonprimary"
        ),
        "profile_selection_rule": (
            "source-native specialized HJB profile only for OU/GMM"
            if resolved_profile == "hjb_pis" else
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
        "round_boundary_resume": {
            "enabled": method == "dpi",
            "resumed_from_iteration": resumed_from_iteration,
            "resume_count": resume_count,
            "state_path": str(rng_path) if method == "dpi" else None,
            "algorithmic_settings_changed": False,
        },
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
        "nonsmooth_extension": (
            "PyTorch minimum-norm norm subgradient at zero; adapted/limited theorem"
            if pde == "mlp_failure" else None
        ),
    }
    if truth_gradient is not None:
        delta = gradient - truth_gradient
        result.update({
            "gradient_component_rmse": float(np.sqrt(np.mean(delta**2))),
            "gradient_vector_rmse": float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
            "gradient_relative_l2": float(np.linalg.norm(delta) / max(np.linalg.norm(truth_gradient), 1e-300)),
        })
    write_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--pde", choices=PDES, required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
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
