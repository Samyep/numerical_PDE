"""Paper-consistent HCFL methods for the homogeneous 1D shallow-water system.

The two retained Euler ideas are translated without changing their method:

1. HLL + a learned correction assembled in Roe characteristic coordinates.
2. Central flux + nonnegative Roe dissipation + proposal-feasibility loss.

Both are tested with symmetric four- and six-cell primitive stencils (never a
five-cell stencil), a hard Tadmor projection in the forward solve, a
conservative water-depth limiter, and a fully-discrete total-entropy limiter.
Training checkpoints are selected only by independent validation rollout
error.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from torch import nn


HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
if str(EXPERIMENTS) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS))

import swe_1d_hcfl as legacy  # noqa: E402


torch.set_num_threads(2)

G = 9.81
NREF = 1024
NCOARSE = 64
FACTOR = NREF // NCOARSE
DT_SNAPSHOT = 5.0e-4
NSNAP = 18
H_FLOOR = 1.0e-5
ENTROPY_RESIDUAL_TOLERANCE = 1.0e-5
FD_ENTROPY_TOLERANCE = 1.0e-8


def stencil_shifts(stencil_cells: int) -> tuple[int, ...]:
    """Return a stencil symmetric about the interface i+1/2.

    Four cells are (i-1, i, i+1, i+2); six cells additionally include i-2
    and i+3. ``torch.roll(state, shift)`` places cell i-shift at location i.
    """
    if stencil_cells == 4:
        return (1, 0, -1, -2)
    if stencil_cells == 6:
        return (2, 1, 0, -1, -2, -3)
    raise ValueError("Only the paper ablations with 4 or 6 cells are allowed")


@dataclass(frozen=True)
class Arm:
    name: str
    model_name: str
    feasibility_weight: float
    stencil_cells: int


ARMS = (
    Arm("hll_roe_correction_s4", "hll_roe", 0.0, 4),
    Arm("hll_roe_correction_s6", "hll_roe", 0.0, 6),
    Arm("central_nonnegative_feas_s4", "central_nonnegative", 1.0e-3, 4),
    Arm("central_nonnegative_feas_s6", "central_nonnegative", 1.0e-3, 6),
)


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Strict fine-grid HLL reference data
# ---------------------------------------------------------------------------


def np_flux(state: np.ndarray) -> np.ndarray:
    h = state[..., 0]
    q = state[..., 1]
    u = q / h
    return np.stack([q, q * u + 0.5 * G * h * h], axis=-1)


def np_hll_pair(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    flux_left = np_flux(left)
    flux_right = np_flux(right)
    h_left = left[..., 0]
    h_right = right[..., 0]
    u_left = left[..., 1] / h_left
    u_right = right[..., 1] / h_right
    c_left = np.sqrt(G * h_left)
    c_right = np.sqrt(G * h_right)
    speed_left = np.minimum(u_left - c_left, u_right - c_right)
    speed_right = np.maximum(u_left + c_left, u_right + c_right)
    middle = (
        speed_right[..., None] * flux_left
        - speed_left[..., None] * flux_right
        + (speed_left * speed_right)[..., None] * (right - left)
    ) / np.maximum(speed_right - speed_left, 1.0e-14)[..., None]
    return np.where(
        (speed_left >= 0.0)[..., None],
        flux_left,
        np.where((speed_right <= 0.0)[..., None], flux_right, middle),
    )


def np_periodic_rhs(state: np.ndarray, dx: float) -> np.ndarray:
    flux = np_hll_pair(state, np.roll(state, -1, axis=-2))
    return -(flux - np.roll(flux, 1, axis=-2)) / dx


def np_periodic_ssprk2(
    state: np.ndarray,
    dt: float,
    dx: float,
) -> np.ndarray:
    stage_one = state + dt * np_periodic_rhs(state, dx)
    if not np.isfinite(stage_one).all() or (stage_one[..., 0] <= 0.0).any():
        raise RuntimeError("HLL reference lost positive depth at SSP-RK2 stage one")
    stage_two = stage_one + dt * np_periodic_rhs(stage_one, dx)
    result = 0.5 * state + 0.5 * stage_two
    if not np.isfinite(result).all() or (result[..., 0] <= 0.0).any():
        raise RuntimeError("HLL reference lost positive depth at SSP-RK2 completion")
    return result


def primitive_to_conservative(h: np.ndarray, u: np.ndarray) -> np.ndarray:
    return np.stack([h, h * u], axis=-1)


def generate_froude_coverage_ic(
    trajectories: int,
    cells: int,
    seed: int,
) -> np.ndarray:
    """Cover subcritical, transcritical, and both supercritical directions."""
    rng = np.random.default_rng(seed)
    regimes = np.resize(
        np.array(
            [
                "dam_break",
                "near_critical_right",
                "near_critical_left",
                "supercritical_right",
                "supercritical_left",
                "counterflow",
            ],
            dtype=object,
        ),
        trajectories,
    )
    rng.shuffle(regimes)
    output = np.empty((trajectories, cells, 2), dtype=np.float64)

    for index, regime in enumerate(regimes):
        h_left, h_right = rng.uniform(0.25, 2.5, 2)
        c_left = math.sqrt(G * h_left)
        c_right = math.sqrt(G * h_right)

        if regime == "dam_break":
            u_left, u_right = rng.uniform(-0.35, 0.35, 2)
        elif regime == "near_critical_right":
            u_left = rng.uniform(0.75, 1.25) * c_left
            u_right = rng.uniform(0.75, 1.25) * c_right
        elif regime == "near_critical_left":
            u_left = -rng.uniform(0.75, 1.25) * c_left
            u_right = -rng.uniform(0.75, 1.25) * c_right
        elif regime == "supercritical_right":
            u_left = rng.uniform(1.15, 1.8) * c_left
            u_right = rng.uniform(1.15, 1.8) * c_right
        elif regime == "supercritical_left":
            u_left = -rng.uniform(1.15, 1.8) * c_left
            u_right = -rng.uniform(1.15, 1.8) * c_right
        elif regime == "counterflow":
            u_left = rng.uniform(0.35, 1.2) * c_left
            u_right = -rng.uniform(0.35, 1.2) * c_right
        else:  # pragma: no cover - exhaustive by construction
            raise ValueError(regime)

        cut = int(rng.integers(cells // 3, 2 * cells // 3))
        h = np.r_[np.full(cut, h_left), np.full(cells - cut, h_right)]
        u = np.r_[np.full(cut, u_left), np.full(cells - cut, u_right)]
        shift = int(rng.integers(0, cells))
        output[index] = primitive_to_conservative(
            np.roll(h, shift), np.roll(u, shift)
        )
    return output


def strict_periodic_reference(
    initial: np.ndarray,
    snapshots: int = NSNAP,
) -> np.ndarray:
    state = np.asarray(initial, dtype=np.float64).copy()
    trajectories, cells, _ = state.shape
    if cells % NCOARSE:
        raise ValueError("Reference cell count must be divisible by 64")
    factor = cells // NCOARSE
    dx = 1.0 / cells

    def restrict(values: np.ndarray) -> np.ndarray:
        return values.reshape(
            trajectories, NCOARSE, factor, 2
        ).mean(axis=2).astype(np.float32)

    saved = [restrict(state)]
    for _ in range(1, snapshots):
        remaining = DT_SNAPSHOT
        while remaining > 1.0e-14:
            h = state[..., 0]
            u = state[..., 1] / h
            max_speed = float(np.max(np.abs(u) + np.sqrt(G * h)))
            dt = min(remaining, 0.20 * dx / max(max_speed, 1.0e-12))
            state = np_periodic_ssprk2(state, dt, dx)
            remaining -= dt
        saved.append(restrict(state))
    return np.stack(saved, axis=1)


def make_training_data(seed: int) -> torch.Tensor:
    ordinary = strict_periodic_reference(
        legacy.generate_ic(220, NREF, 21000 + seed, ood=False)
    )
    broad = strict_periodic_reference(
        legacy.generate_ic(160, NREF, 22000 + seed, ood=True)
    )
    froude = strict_periodic_reference(
        generate_froude_coverage_ic(200, NREF, 23000 + seed)
    )
    return torch.from_numpy(np.concatenate([ordinary, broad, froude], axis=0))


def make_validation_data(seed: int) -> torch.Tensor:
    ordinary = strict_periodic_reference(
        legacy.generate_ic(44, NREF, 24000 + seed, ood=False)
    )
    broad = strict_periodic_reference(
        legacy.generate_ic(32, NREF, 25000 + seed, ood=True)
    )
    froude = strict_periodic_reference(
        generate_froude_coverage_ic(60, NREF, 26000 + seed)
    )
    return torch.from_numpy(np.concatenate([ordinary, broad, froude], axis=0))


def load_or_make_data(
    output: Path,
    seed: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Cache deterministic reference trajectories locally between long runs."""
    cache = output / f"data_cache_seed{seed}.pt"
    if cache.exists():
        payload = torch.load(cache, map_location="cpu", weights_only=True)
        train_data = payload["train"]
        validation_data = payload["validation"]
        expected_train = (580, NSNAP, NCOARSE, 2)
        expected_validation = (136, NSNAP, NCOARSE, 2)
        if (
            tuple(train_data.shape) == expected_train
            and tuple(validation_data.shape) == expected_validation
        ):
            print(f"Loading deterministic data cache: {cache}", flush=True)
            return train_data, validation_data
        print(f"Ignoring stale data cache: {cache}", flush=True)

    print("Generating strict SWE training data...", flush=True)
    train_data = make_training_data(seed)
    print("Generating independent SWE validation data...", flush=True)
    validation_data = make_validation_data(seed)
    torch.save({"train": train_data, "validation": validation_data}, cache)
    return train_data, validation_data


# ---------------------------------------------------------------------------
# Torch fluxes, Roe waves, entropy, and learned proposals
# ---------------------------------------------------------------------------


def t_flux(state: torch.Tensor) -> torch.Tensor:
    h = state[..., 0].clamp_min(H_FLOOR)
    q = state[..., 1]
    u = q / h
    return torch.stack([q, q * u + 0.5 * G * h * h], dim=-1)


def primitive(state: torch.Tensor) -> torch.Tensor:
    h = state[..., 0].clamp_min(H_FLOOR)
    return torch.stack([h, state[..., 1] / h], dim=-1)


def entropy(state: torch.Tensor) -> torch.Tensor:
    h = state[..., 0].clamp_min(H_FLOOR)
    q = state[..., 1]
    return 0.5 * q * q / h + 0.5 * G * h * h


def entropy_variables(state: torch.Tensor) -> torch.Tensor:
    h = state[..., 0].clamp_min(H_FLOOR)
    u = state[..., 1] / h
    return torch.stack([G * h - 0.5 * u * u, u], dim=-1)


def entropy_potential(state: torch.Tensor) -> torch.Tensor:
    h = state[..., 0].clamp_min(H_FLOOR)
    return 0.5 * G * h * state[..., 1]


def entropy_flux(state: torch.Tensor) -> torch.Tensor:
    h = state[..., 0].clamp_min(H_FLOOR)
    u = state[..., 1] / h
    return u * (0.5 * h * u * u + G * h * h)


def t_hll(state: torch.Tensor) -> torch.Tensor:
    left = state
    right = torch.roll(state, -1, dims=-2)
    flux_left = t_flux(left)
    flux_right = t_flux(right)
    h_left = left[..., 0].clamp_min(H_FLOOR)
    h_right = right[..., 0].clamp_min(H_FLOOR)
    u_left = left[..., 1] / h_left
    u_right = right[..., 1] / h_right
    c_left = torch.sqrt(G * h_left)
    c_right = torch.sqrt(G * h_right)
    speed_left = torch.minimum(u_left - c_left, u_right - c_right)
    speed_right = torch.maximum(u_left + c_left, u_right + c_right)
    middle = (
        speed_right[..., None] * flux_left
        - speed_left[..., None] * flux_right
        + (speed_left * speed_right)[..., None] * (right - left)
    ) / (speed_right - speed_left).clamp_min(1.0e-10)[..., None]
    return torch.where(
        (speed_left >= 0.0)[..., None],
        flux_left,
        torch.where((speed_right <= 0.0)[..., None], flux_right, middle),
    )


def roe_waves(
    state: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Return R, alpha, signed eigenvalues, and entropy-fixed magnitudes."""
    left = state
    right = torch.roll(state, -1, dims=-2)
    h_left = left[..., 0].clamp_min(H_FLOOR)
    h_right = right[..., 0].clamp_min(H_FLOOR)
    u_left = left[..., 1] / h_left
    u_right = right[..., 1] / h_right
    root_left = torch.sqrt(h_left)
    root_right = torch.sqrt(h_right)
    roe_u = (
        root_left * u_left + root_right * u_right
    ) / (root_left + root_right).clamp_min(1.0e-10)
    roe_h = 0.5 * (h_left + h_right)
    roe_c = torch.sqrt(G * roe_h)
    eigenvalues = torch.stack([roe_u - roe_c, roe_u + roe_c], dim=-1)
    ones = torch.ones_like(roe_u)
    right_minus = torch.stack([ones, roe_u - roe_c], dim=-1)
    right_plus = torch.stack([ones, roe_u + roe_c], dim=-1)
    matrix = torch.stack([right_minus, right_plus], dim=-1)
    jump = right - left
    alpha = torch.linalg.solve(matrix, jump.unsqueeze(-1)).squeeze(-1)

    absolute = eigenvalues.abs()
    entropy_fix = (0.10 * roe_c).clamp_min(1.0e-6)[..., None]
    fixed = torch.where(
        absolute < entropy_fix,
        0.5 * (absolute * absolute / entropy_fix + entropy_fix),
        absolute,
    )
    return matrix, alpha, eigenvalues, fixed


def entropy_residual(flux: torch.Tensor, state: torch.Tensor) -> torch.Tensor:
    right = torch.roll(state, -1, dims=-2)
    normal = entropy_variables(right) - entropy_variables(state)
    bound = entropy_potential(right) - entropy_potential(state)
    return (normal * flux).sum(dim=-1) - bound


def strict_entropy_projection(
    proposal: torch.Tensor,
    state: torch.Tensor,
) -> torch.Tensor:
    output = proposal
    epsilon = torch.finfo(proposal.dtype).eps
    state64 = state.double()
    right = torch.roll(state64, -1, dims=-2)
    normal = entropy_variables(right) - entropy_variables(state64)
    bound = entropy_potential(right) - entropy_potential(state64)
    norm_squared = (normal * normal).sum(dim=-1)

    for _ in range(4):
        flux64 = output.double()
        residual = (normal * flux64).sum(dim=-1) - bound
        scale = bound.abs() + (normal * flux64).abs().sum(dim=-1)
        target = -8.0 * epsilon * scale
        mask = (residual > target) & (norm_squared > 1.0e-14)
        if not bool(mask.any()):
            break
        correction = torch.zeros_like(residual)
        correction[mask] = (
            residual[mask] - target[mask]
        ) / norm_squared[mask]
        output = (flux64 - correction[..., None] * normal).to(proposal.dtype)
    return output


class HLLRoeCorrectionFlux(nn.Module):
    def __init__(
        self,
        mean: np.ndarray,
        std: np.ndarray,
        width: int = 72,
        stencil_cells: int = 4,
    ):
        super().__init__()
        self.stencil_shifts = stencil_shifts(stencil_cells)
        self.net = nn.Sequential(
            nn.Linear(2 * stencil_cells, width),
            nn.Tanh(),
            nn.Linear(width, width),
            nn.Tanh(),
            nn.Linear(width, 2),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def coefficients(self, state: torch.Tensor) -> torch.Tensor:
        values = primitive(state)
        features = torch.cat(
            [
                (torch.roll(values, shift, dims=-2) - self.mean) / self.std
                for shift in self.stencil_shifts
            ],
            dim=-1,
        )
        return torch.tanh(self.net(features))

    def forward(self, state: torch.Tensor) -> torch.Tensor:
        correction_coefficients = self.coefficients(state)
        matrix, alpha, _, speeds = roe_waves(state)
        correction = torch.einsum(
            "...ij,...j->...i",
            matrix,
            correction_coefficients * speeds * alpha,
        )
        return t_hll(state) - 0.5 * correction


class CentralNonnegativeRoeFlux(nn.Module):
    def __init__(
        self,
        mean: np.ndarray,
        std: np.ndarray,
        width: int = 72,
        stencil_cells: int = 4,
    ):
        super().__init__()
        self.stencil_shifts = stencil_shifts(stencil_cells)
        self.net = nn.Sequential(
            nn.Linear(2 * stencil_cells, width),
            nn.Tanh(),
            nn.Linear(width, width),
            nn.Tanh(),
            nn.Linear(width, 2),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def coefficients(self, state: torch.Tensor) -> torch.Tensor:
        values = primitive(state)
        features = torch.cat(
            [
                (torch.roll(values, shift, dims=-2) - self.mean) / self.std
                for shift in self.stencil_shifts
            ],
            dim=-1,
        )
        return 1.0 + torch.tanh(self.net(features))

    def forward(self, state: torch.Tensor) -> torch.Tensor:
        multipliers = self.coefficients(state)
        matrix, alpha, _, speeds = roe_waves(state)
        dissipation = torch.einsum(
            "...ij,...j->...i", matrix, multipliers * speeds * alpha
        )
        right = torch.roll(state, -1, dims=-2)
        central = 0.5 * (t_flux(state) + t_flux(right))
        return central - 0.5 * dissipation


class Solver(nn.Module):
    def __init__(
        self,
        model_name: str,
        mean: np.ndarray,
        std: np.ndarray,
        width: int = 72,
        stencil_cells: int = 4,
    ) -> None:
        super().__init__()
        if model_name == "hll_roe":
            self.flux_net = HLLRoeCorrectionFlux(
                mean, std, width, stencil_cells
            )
        elif model_name == "central_nonnegative":
            self.flux_net = CentralNonnegativeRoeFlux(
                mean, std, width, stencil_cells
            )
        else:
            raise ValueError(model_name)

    def raw_flux(self, state: torch.Tensor) -> torch.Tensor:
        return self.flux_net(state)

    def flux(self, state: torch.Tensor) -> torch.Tensor:
        return strict_entropy_projection(self.raw_flux(state), state)

    def one_step(self, state: torch.Tensor) -> torch.Tensor:
        return fv_step(state, self.flux(state), DT_SNAPSHOT * NCOARSE)


# ---------------------------------------------------------------------------
# Full periodic hard-safety stack
# ---------------------------------------------------------------------------


def fv_step(
    state: torch.Tensor,
    flux: torch.Tensor,
    lam: float,
) -> torch.Tensor:
    return state - lam * (flux - torch.roll(flux, 1, dims=-2))


def admissible(state: torch.Tensor) -> torch.Tensor:
    return torch.isfinite(state).all(dim=-1) & (state[..., 0] >= H_FLOOR)


def total_entropy(state: torch.Tensor) -> torch.Tensor:
    return entropy(state.double()).sum(dim=-1)


def local_depth_limiter(
    state: torch.Tensor,
    high_flux: torch.Tensor,
    low_flux: torch.Tensor,
    lam: float,
    max_outer: int = 12,
    bisections: int = 30,
) -> tuple[torch.Tensor, torch.Tensor]:
    delta_flux = high_flux - low_flux
    low_state = fv_step(state, low_flux, lam)
    if not bool(admissible(low_state).all()):
        raise RuntimeError("Low-order periodic SWE state is not admissible")
    batch, cells, _ = state.shape
    alpha = torch.ones(batch, cells, dtype=state.dtype, device=state.device)

    for _ in range(max_outer):
        flux = low_flux + alpha[..., None] * delta_flux
        current = fv_step(state, flux, lam)
        good = admissible(current)
        if bool(good.all()):
            return flux, alpha
        bad = ~good
        direction = current - low_state
        count = int(bad.sum())
        lo = torch.zeros(count, dtype=torch.float64, device=state.device)
        hi = torch.ones(count, dtype=torch.float64, device=state.device)
        base_state = low_state[bad].double()
        delta_state = direction[bad].double()
        for _ in range(bisections):
            mid = 0.5 * (lo + hi)
            okay = admissible(base_state + mid[:, None] * delta_state)
            lo = torch.where(okay, mid, lo)
            hi = torch.where(okay, hi, mid)
        theta = torch.ones_like(alpha)
        theta[bad] = torch.clamp(lo.to(state.dtype) - 1.0e-5, min=0.0)
        interface_factor = torch.minimum(theta, torch.roll(theta, -1, dims=-1))
        alpha = alpha * interface_factor

    flux = low_flux + alpha[..., None] * delta_flux
    current = fv_step(state, flux, lam)
    if not bool(admissible(current).all()):
        failed = (~admissible(current)).any(dim=-1)
        alpha[failed] = 0.0
        flux = low_flux + alpha[..., None] * delta_flux
    if not bool(admissible(fv_step(state, flux, lam)).all()):
        raise RuntimeError("Periodic SWE depth limiter failed")
    return flux, alpha


def global_entropy_limiter(
    state: torch.Tensor,
    high_flux: torch.Tensor,
    low_flux: torch.Tensor,
    lam: float,
    bisections: int = 40,
) -> tuple[torch.Tensor, torch.Tensor]:
    low_state = fv_step(state, low_flux, lam)
    high_state = fv_step(state, high_flux, lam)
    target = total_entropy(state)
    beta = torch.ones(state.shape[0], dtype=state.dtype, device=state.device)
    need = total_entropy(high_state) > target

    if bool(need.any()):
        count = int(need.sum())
        lo = torch.zeros(count, dtype=torch.float64, device=state.device)
        hi = torch.ones(count, dtype=torch.float64, device=state.device)
        low_selected = low_state[need].double()
        direction = (high_state - low_state)[need].double()
        target_selected = target[need]
        for _ in range(bisections):
            mid = 0.5 * (lo + hi)
            okay = total_entropy(
                low_selected + mid[:, None, None] * direction
            ) <= target_selected
            lo = torch.where(okay, mid, lo)
            hi = torch.where(okay, hi, mid)
        beta[need] = torch.clamp(lo.to(state.dtype) - 1.0e-6, min=0.0)

    delta = high_flux - low_flux
    flux = low_flux + beta[:, None, None] * delta
    rounded_bad = total_entropy(fv_step(state, flux, lam)) > (
        target + FD_ENTROPY_TOLERANCE
    )
    if bool(rounded_bad.any()):
        beta[rounded_bad] = 0.0
        flux = low_flux + beta[:, None, None] * delta
    return flux, beta


@torch.no_grad()
def advance_safe_snapshot(
    model: Solver,
    state: torch.Tensor,
    cfl: float = 0.42,
) -> tuple[torch.Tensor, dict[str, float]]:
    dx = 1.0 / state.shape[-2]
    remaining = DT_SNAPSHOT
    totals: defaultdict[str, float] = defaultdict(float)
    totals["fd_beta_min"] = 1.0
    totals["maximum_entropy_residual"] = -float("inf")
    totals["maximum_total_entropy_change"] = -float("inf")

    while remaining > 1.0e-14:
        values = primitive(state)
        max_speed = float(
            (values[..., 1].abs() + torch.sqrt(G * values[..., 0])).max()
        )
        dt = min(remaining, cfl * dx / max(max_speed, 1.0e-12))

        for _ in range(24):
            lam = dt / dx
            low_flux = strict_entropy_projection(t_hll(state), state)
            low_state = fv_step(state, low_flux, lam)
            entropy_ok = bool(
                (
                    total_entropy(low_state)
                    <= total_entropy(state) + FD_ENTROPY_TOLERANCE
                ).all()
            )
            if bool(admissible(low_state).all()) and entropy_ok:
                break
            dt *= 0.5
            totals["low_order_dt_halvings"] += 1.0
        else:
            raise RuntimeError("Could not establish periodic SWE low-order premise")

        high_flux = model.flux(state)
        local_flux, alpha = local_depth_limiter(
            state, high_flux, low_flux, lam
        )
        final_flux, beta = global_entropy_limiter(
            state, local_flux, low_flux, lam
        )
        next_state = fv_step(state, final_flux, lam)
        if not bool(admissible(next_state).all()):
            raise RuntimeError("Safe periodic SWE update lost positive depth")

        residual = entropy_residual(final_flux.double(), state.double())
        entropy_change = total_entropy(next_state) - total_entropy(state)
        totals["local_active"] += float((alpha < 1.0 - 1.0e-7).sum())
        totals["local_total"] += float(alpha.numel())
        totals["local_alpha_sum"] += float(alpha.sum())
        totals["fd_active"] += float((beta < 1.0 - 1.0e-7).sum())
        totals["fd_total"] += float(beta.numel())
        totals["fd_beta_sum"] += float(beta.sum())
        totals["fd_beta_min"] = min(totals["fd_beta_min"], float(beta.min()))
        totals["entropy_violations"] += float(
            (residual > ENTROPY_RESIDUAL_TOLERANCE).sum()
        )
        totals["entropy_total"] += float(residual.numel())
        totals["maximum_entropy_residual"] = max(
            totals["maximum_entropy_residual"], float(residual.max())
        )
        totals["maximum_total_entropy_change"] = max(
            totals["maximum_total_entropy_change"], float(entropy_change.max())
        )
        totals["substeps"] += 1.0
        state = next_state
        remaining -= dt
    return state, dict(totals)


def merge_stats(total: defaultdict[str, float], step: dict[str, float]) -> None:
    for key, value in step.items():
        if key in ("fd_beta_min",):
            total[key] = min(total[key], value)
        elif key.startswith("maximum_"):
            total[key] = max(total[key], value)
        else:
            total[key] += value


@torch.no_grad()
def evaluate_dataset(
    model: Solver,
    data: torch.Tensor,
    state_std: torch.Tensor,
    split: str,
) -> dict[str, Any]:
    state = data[:, 0].clone()
    errors: list[float] = []
    minimum_depth = float(state[..., 0].min())
    totals: defaultdict[str, float] = defaultdict(float)
    totals["fd_beta_min"] = 1.0
    totals["maximum_entropy_residual"] = -float("inf")
    totals["maximum_total_entropy_change"] = -float("inf")
    started = time.perf_counter()
    for snapshot in range(1, data.shape[1]):
        state, step = advance_safe_snapshot(model, state)
        merge_stats(totals, step)
        minimum_depth = min(minimum_depth, float(state[..., 0].min()))
        errors.append(
            float((((state - data[:, snapshot]) / state_std) ** 2).mean())
        )
    return {
        "split": split,
        "trajectory_count": int(data.shape[0]),
        "rollout_nrmse": float(np.sqrt(np.mean(errors))),
        "minimum_depth": minimum_depth,
        "entropy_violation_rate": totals["entropy_violations"]
        / max(totals["entropy_total"], 1.0),
        "maximum_entropy_residual": totals["maximum_entropy_residual"],
        "local_limiter_intervention_rate": totals["local_active"]
        / max(totals["local_total"], 1.0),
        "mean_local_alpha": totals["local_alpha_sum"]
        / max(totals["local_total"], 1.0),
        "fd_entropy_intervention_rate": totals["fd_active"]
        / max(totals["fd_total"], 1.0),
        "minimum_fd_beta": totals["fd_beta_min"],
        "mean_fd_beta": totals["fd_beta_sum"]
        / max(totals["fd_total"], 1.0),
        "maximum_total_entropy_change": totals[
            "maximum_total_entropy_change"
        ],
        "mean_substeps_per_snapshot": totals["substeps"]
        / (data.shape[1] - 1),
        "evaluation_seconds": time.perf_counter() - started,
    }


# ---------------------------------------------------------------------------
# Validation-controlled training
# ---------------------------------------------------------------------------


def clone_state_dict(model: nn.Module) -> dict[str, torch.Tensor]:
    return {
        name: value.detach().cpu().clone()
        for name, value in model.state_dict().items()
    }


@torch.no_grad()
def deterministic_one_step_loss(
    model: Solver,
    data: torch.Tensor,
    state_std: torch.Tensor,
    batch_size: int = 256,
) -> float:
    inputs = data[:, :-1].reshape(-1, NCOARSE, 2)
    targets = data[:, 1:].reshape(-1, NCOARSE, 2)
    total = 0.0
    count = 0
    was_training = model.training
    model.eval()
    for start in range(0, inputs.shape[0], batch_size):
        stop = min(start + batch_size, inputs.shape[0])
        prediction = model.one_step(inputs[start:stop])
        normalized = (prediction - targets[start:stop]) / state_std
        total += float((normalized * normalized).sum())
        count += normalized.numel()
    model.train(was_training)
    return total / count


def train_to_convergence(
    arm: Arm,
    train_data: torch.Tensor,
    validation_data: torch.Tensor,
    mean: np.ndarray,
    std: np.ndarray,
    state_std: torch.Tensor,
    args: argparse.Namespace,
    output: Path,
) -> tuple[Solver, dict[str, Any], list[dict[str, Any]]]:
    torch.manual_seed(31000 + args.seed)
    model = Solver(
        arm.model_name,
        mean,
        std,
        width=args.width,
        stencil_cells=arm.stencil_cells,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    generator = torch.Generator().manual_seed(32000 + args.seed)
    best_metric = float("inf")
    best_update = -1
    best_state: dict[str, torch.Tensor] | None = None
    plateau_anchor = float("inf")
    checks_without_improvement = 0
    converged = False
    stop_reason = "max_updates"
    curve: list[dict[str, Any]] = []
    started = time.perf_counter()

    def validate(update: int, batch_loss: float | None) -> float:
        nonlocal best_metric, best_update, best_state
        train_loss = deterministic_one_step_loss(model, train_data, state_std)
        validation_loss = deterministic_one_step_loss(
            model, validation_data, state_std
        )
        rollout = evaluate_dataset(
            model, validation_data, state_std, "validation"
        )
        metric = float(rollout["rollout_nrmse"])
        improved = metric < best_metric
        if improved:
            best_metric = metric
            best_update = update
            best_state = clone_state_dict(model)
            torch.save(
                best_state,
                output / f"{arm.name}_converged_best_seed{args.seed}.pt",
            )
        row = {
            "seed": args.seed,
            "arm": arm.name,
            "model": arm.model_name,
            "update": update,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "last_minibatch_loss": batch_loss,
            "full_train_one_step_loss": train_loss,
            "validation_one_step_loss": validation_loss,
            "validation_rollout_nrmse": metric,
            "validation_minimum_depth": rollout["minimum_depth"],
            "new_absolute_best": improved,
        }
        curve.append(row)
        print(json.dumps(row), flush=True)
        return metric

    validate(0, None)
    last_loss: float | None = None
    for update in range(1, args.max_updates + 1):
        model.train()
        indices = torch.randint(
            0, train_data.shape[0], (args.batch_size,), generator=generator
        )
        times = torch.randint(
            0, train_data.shape[1] - 1, (args.batch_size,), generator=generator
        )
        inputs = train_data[indices, times]
        targets = train_data[indices, times + 1]
        raw_flux = model.raw_flux(inputs)
        projected = strict_entropy_projection(raw_flux, inputs)
        prediction = fv_step(inputs, projected, DT_SNAPSHOT * NCOARSE)
        trajectory_loss = (((prediction - targets) / state_std) ** 2).mean()
        if arm.feasibility_weight:
            raw_residual = entropy_residual(raw_flux, inputs)
            feasibility_loss = torch.relu(raw_residual).square().mean()
            loss = trajectory_loss + arm.feasibility_weight * feasibility_loss
        else:
            loss = trajectory_loss

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        last_loss = float(loss.detach())

        if update % args.validation_interval and update != args.max_updates:
            continue
        metric = validate(update, last_loss)
        if update < args.fixed_budget_updates:
            continue
        if plateau_anchor == float("inf"):
            plateau_anchor = best_metric
            checks_without_improvement = 0
            continue
        if metric <= plateau_anchor * (1.0 - args.minimum_relative_improvement):
            plateau_anchor = metric
            checks_without_improvement = 0
        else:
            checks_without_improvement += 1
        if checks_without_improvement < args.plateau_patience:
            continue
        current_lr = float(optimizer.param_groups[0]["lr"])
        if current_lr <= args.minimum_learning_rate * (1.0 + 1.0e-12):
            converged = True
            stop_reason = "validation_plateau_at_minimum_learning_rate"
            break
        new_lr = max(
            current_lr * args.learning_rate_factor,
            args.minimum_learning_rate,
        )
        for group in optimizer.param_groups:
            group["lr"] = new_lr
        checks_without_improvement = 0
        plateau_anchor = best_metric
        print(
            json.dumps(
                {
                    "arm": arm.name,
                    "event": "reduce_learning_rate",
                    "update": update,
                    "old_learning_rate": current_lr,
                    "new_learning_rate": new_lr,
                }
            ),
            flush=True,
        )

    if best_state is None:
        raise RuntimeError("No SWE validation checkpoint was created")
    model.load_state_dict(best_state)
    model.eval()
    summary = {
        "seed": args.seed,
        "arm": arm.name,
        "model": arm.model_name,
        "stencil_cells": arm.stencil_cells,
        "parameter_count": sum(p.numel() for p in model.parameters()),
        "best_update": best_update,
        "stop_update": int(curve[-1]["update"]),
        "best_validation_rollout_nrmse": best_metric,
        "final_learning_rate": float(optimizer.param_groups[0]["lr"]),
        "converged": converged,
        "stop_reason": stop_reason,
        "training_seconds": time.perf_counter() - started,
        "proposal_feasibility_weight": arm.feasibility_weight,
    }
    return model, summary, curve


def prepare_statistics(
    train_data: torch.Tensor,
) -> tuple[np.ndarray, np.ndarray, torch.Tensor]:
    values = primitive(train_data)
    mean = values.mean(dim=(0, 1, 2)).numpy()
    std = values.std(dim=(0, 1, 2)).numpy()
    state_std = train_data.std(dim=(0, 1, 2))
    return mean, std, state_std


def checkpoint_path(output: Path, arm: Arm, seed: int) -> Path:
    return output / f"{arm.name}_converged_best_seed{seed}.pt"


def load_model(
    arm: Arm,
    checkpoint: Path,
    mean: np.ndarray,
    std: np.ndarray,
    width: int,
) -> Solver:
    model = Solver(
        arm.model_name,
        mean,
        std,
        width,
        stencil_cells=arm.stencil_cells,
    )
    model.load_state_dict(
        torch.load(checkpoint, map_location="cpu", weights_only=True)
    )
    model.eval()
    return model


def self_test() -> None:
    mean = np.array([1.0, 0.0], dtype=np.float32)
    std = np.ones(2, dtype=np.float32)
    h = torch.linspace(0.6, 1.4, 32)[None]
    u = 0.3 * torch.sin(torch.linspace(0.0, 2.0 * math.pi, 32))[None]
    state = torch.stack([h, h * u], dim=-1)
    matrix, alpha, eigenvalues, speeds = roe_waves(state)
    jump = torch.roll(state, -1, dims=-2) - state
    reconstructed = torch.einsum("...ij,...j->...i", matrix, alpha)
    if float((reconstructed - jump).abs().max()) > 2.0e-6:
        raise RuntimeError("SWE Roe reconstruction failed")
    if not bool((speeds >= 0.0).all()):
        raise RuntimeError("Entropy-fixed SWE speeds must be nonnegative")
    if not bool((eigenvalues[..., 0] < eigenvalues[..., 1]).all()):
        raise RuntimeError("SWE wave ordering failed")

    if {arm.stencil_cells for arm in ARMS} != {4, 6}:
        raise RuntimeError("SWE comparison must contain exactly 4/6-cell stencils")
    if any("s5" in arm.name or arm.stencil_cells == 5 for arm in ARMS):
        raise RuntimeError("A forbidden 5-cell SWE arm was registered")

    for arm in ARMS:
        model = Solver(
            arm.model_name,
            mean,
            std,
            width=16,
            stencil_cells=arm.stencil_cells,
        )
        raw = model.raw_flux(state)
        projected = model.flux(state)
        if raw.shape != state.shape or not bool(torch.isfinite(raw).all()):
            raise RuntimeError(f"Invalid raw flux for {arm.name}")
        if float(
            entropy_residual(projected.double(), state.double())
            .detach()
            .max()
        ) > 1.0e-8:
            raise RuntimeError(f"Hard entropy projection failed for {arm.name}")
        constant = torch.tensor([[[1.2, 0.36]]]).expand(1, 32, 2).clone()
        expected = t_flux(constant)
        consistency_error = float(
            (model.raw_flux(constant) - expected).detach().abs().max()
        )
        if consistency_error > 2.0e-6:
            raise RuntimeError(
                f"Equal-state consistency failed for {arm.name}: "
                f"{consistency_error}"
            )
        translated = model.raw_flux(torch.roll(state, 5, dims=-2))
        equivariance_error = float(
            (translated - torch.roll(raw, 5, dims=-2)).detach().abs().max()
        )
        if equivariance_error > 2.0e-6:
            raise RuntimeError(f"Translation equivariance failed for {arm.name}")
        if arm.model_name == "central_nonnegative":
            coefficients = model.flux_net.coefficients(state)
            if not bool((coefficients >= 0.0).all()):
                raise RuntimeError("Nonnegative SWE Roe multiplier failed")

    print("self-test passed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("train", "evaluate", "all"), default="all")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--lr", type=float, default=3.0e-4)
    parser.add_argument("--fixed-budget-updates", type=int, default=1100)
    parser.add_argument("--max-updates", type=int, default=50000)
    parser.add_argument("--validation-interval", type=int, default=100)
    parser.add_argument("--plateau-patience", type=int, default=10)
    parser.add_argument("--minimum-relative-improvement", type=float, default=1.0e-3)
    parser.add_argument("--learning-rate-factor", type=float, default=0.3)
    parser.add_argument("--minimum-learning-rate", type=float, default=3.0e-6)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    train_data, validation_data = load_or_make_data(output, args.seed)
    mean, std, state_std = prepare_statistics(train_data)

    models: dict[str, Solver] = {}
    summaries: list[dict[str, Any]] = []
    if args.phase in ("train", "all"):
        for arm in ARMS:
            checkpoint = checkpoint_path(output, arm, args.seed)
            summary_path = output / f"report_{arm.name}_seed{args.seed}.json"
            if args.resume and checkpoint.exists() and summary_path.exists():
                report = json.loads(summary_path.read_text(encoding="utf-8"))
                models[arm.name] = load_model(
                    arm, checkpoint, mean, std, args.width
                )
                summaries.append(report["convergence"])
                continue
            model, summary, curve = train_to_convergence(
                arm,
                train_data,
                validation_data,
                mean,
                std,
                state_std,
                args,
                output,
            )
            if not summary["converged"]:
                checkpoint.unlink(missing_ok=True)
                raise RuntimeError(
                    f"{arm.name} reached the cap without validation convergence"
                )
            report = {
                "method": (
                    "same Roe-wave HCFL construction as retained Euler methods"
                ),
                "stencil_cells": arm.stencil_cells,
                "stencil_shifts": list(stencil_shifts(arm.stencil_cells)),
                "training_data_shape": list(train_data.shape),
                "validation_data_shape": list(validation_data.shape),
                "training_trajectory_counts": {
                    "legacy_in_distribution": 220,
                    "legacy_broad": 160,
                    "froude_coverage": 200,
                },
                "validation_trajectory_counts": {
                    "legacy_in_distribution": 44,
                    "legacy_broad": 32,
                    "froude_coverage": 60,
                },
                "primitive_mean": mean.tolist(),
                "primitive_std": std.tolist(),
                "state_std": state_std.tolist(),
                "convergence": summary,
            }
            summary_path.write_text(
                json.dumps(json_ready(report), indent=2), encoding="utf-8"
            )
            write_csv(
                output / f"training_curve_{arm.name}_seed{args.seed}.csv",
                curve,
            )
            models[arm.name] = model
            summaries.append(summary)
        write_csv(output / f"convergence_seed{args.seed}.csv", summaries)
    else:
        for arm in ARMS:
            models[arm.name] = load_model(
                arm,
                checkpoint_path(output, arm, args.seed),
                mean,
                std,
                args.width,
            )
            report = json.loads(
                (output / f"report_{arm.name}_seed{args.seed}.json").read_text(
                    encoding="utf-8"
                )
            )
            summaries.append(report["convergence"])

    if args.phase in ("evaluate", "all"):
        rows = [
            {
                "arm": arm.name,
                "stencil_cells": arm.stencil_cells,
                **evaluate_dataset(
                    models[arm.name], validation_data, state_std, "validation"
                ),
            }
            for arm in ARMS
        ]
        write_csv(output / f"validation_metrics_seed{args.seed}.csv", rows)
        print(json.dumps(json_ready(rows), indent=2), flush=True)


if __name__ == "__main__":
    main()
