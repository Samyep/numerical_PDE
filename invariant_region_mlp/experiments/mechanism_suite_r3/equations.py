"""Round-3 equations, fixed test distributions, and equation factory."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.interpolate import RectBivariateSpline


HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
PACKAGE_ROOT = HERE.parents[1]
for directory in (EXPERIMENTS, EXPERIMENTS / "mechanism_suite"):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from candidate_screen.candidates import L1ControlHJB, LQGameHJB  # noqa: E402
from mechanism_suite.equations import (  # noqa: E402
    ExactEquation,
    RidgeLSEHJB,
    _as_points,
    stable_seed,
)
from mechanism_suite.norm_hjb import NormDriverHJB  # noqa: E402

from .reference_c3 import load_reference, terminal_profile


BASE_SEED = 20261207
POINT_COUNT = 1200
SQRT2 = math.sqrt(2.0)
P4_REFERENCE = (
    PACKAGE_ROOT
    / "results"
    / "mechanism_suite"
    / "reference_cache"
    / "norm_hjb_beta2_lambda1_T0p5.npz"
)
RESULTS_ROOT = PACKAGE_ROOT / "results" / "mechanism_suite_r3"
C3_REFERENCE = RESULTS_ROOT / "reference_cache" / "c3_production_reference.npz"


def make_points(
    equation: Any,
    *,
    pde_id: str,
    n_points: int = POINT_COUNT,
    seed: int = BASE_SEED,
) -> dict[str, np.ndarray]:
    """Create the pre-registered fixed points and 20% validation split."""

    rng = np.random.default_rng(
        np.random.SeedSequence(
            [seed, equation.d, n_points, stable_seed(pde_id)]
        )
    )
    t = rng.uniform(
        0.0, np.nextafter(equation.T, 0.0), n_points
    ).astype(np.float64)
    if equation.family == "lqg":
        direction = rng.standard_normal((n_points, equation.d), dtype=np.float64)
        direction /= np.linalg.norm(direction, axis=1, keepdims=True)
        radius = rng.uniform(0.0, 1.0, n_points).astype(np.float64)
        x = radius[:, None] * direction
    else:
        x = rng.uniform(-1.0, 1.0, (n_points, equation.d)).astype(np.float64)
    is_validation = np.zeros(n_points, dtype=bool)
    selected = rng.choice(n_points, size=int(round(0.2 * n_points)), replace=False)
    is_validation[selected] = True
    return {"t": t, "x": x, "is_validation": is_validation}


@lru_cache(maxsize=2)
def _load_c3_spline(path_text: str) -> tuple[RectBivariateSpline, dict[str, Any]]:
    tau, s, value, metadata = load_reference(Path(path_text))
    spline = RectBivariateSpline(tau, s, value, kx=3, ky=3, s=0.0)
    return spline, metadata


@dataclass(frozen=True)
class ProductionCubicHJB(ExactEquation):
    d: int
    reference_path: str = str(C3_REFERENCE)
    kappa: float = 0.5
    A: float = 2.0
    beta: float = 2.0
    T: float = 0.25
    sigma: float = SQRT2
    direction_seed: int = BASE_SEED
    name: str = "C3_cubic_viscous_hj"
    family: str = "cubic_hj"

    def __post_init__(self) -> None:
        rng = np.random.default_rng(
            np.random.SeedSequence([self.direction_seed, self.d, 303])
        )
        w = rng.standard_normal(self.d, dtype=np.float64)
        w /= np.linalg.norm(w)
        object.__setattr__(self, "w", w)
        object.__setattr__(self, "mu", 0.0)
        _, metadata = _load_c3_spline(self.reference_path)
        parameters = metadata["parameters"]
        for key in ("kappa", "A", "beta", "T"):
            if not math.isclose(
                float(parameters[key]), float(getattr(self, key)),
                rel_tol=0.0, abs_tol=1e-14,
            ):
                raise ValueError(f"C3 reference {key} does not match equation")

    def _tau_s(self, t: np.ndarray | float, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        time, state = _as_points(t, x, self.d)
        scalar = np.einsum("...d,d->...", state, self.w)
        tau = np.broadcast_to(self.T - time, np.shape(scalar))
        return tau, scalar

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        tau, scalar = self._tau_s(t, x)
        spline, _ = _load_c3_spline(self.reference_path)
        result = spline.ev(
            tau.reshape(-1), scalar.reshape(-1)
        ).reshape(np.shape(scalar))
        at_terminal = np.abs(tau) <= 8.0 * np.finfo(np.float64).eps
        if np.any(at_terminal):
            result = np.where(
                at_terminal,
                terminal_profile(scalar, A=self.A, beta=self.beta),
                result,
            )
        return result

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        tau, scalar = self._tau_s(t, x)
        spline, _ = _load_c3_spline(self.reference_path)
        derivative = spline.ev(
            tau.reshape(-1), scalar.reshape(-1), dx=0, dy=1
        ).reshape(np.shape(scalar))
        at_terminal = np.abs(tau) <= 8.0 * np.finfo(np.float64).eps
        if np.any(at_terminal):
            derivative = np.where(
                at_terminal,
                self.A * np.tanh(self.beta * scalar),
                derivative,
            )
        return self.sigma * derivative[..., None] * self.w

    def terminal(self, x: np.ndarray) -> np.ndarray:
        scalar = np.asarray(x, dtype=np.float64) @ self.w
        return terminal_profile(scalar, A=self.A, beta=self.beta)

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        del u
        norm = np.linalg.norm(np.asarray(z, dtype=np.float64), axis=-1)
        return -self.kappa * (norm / self.sigma) ** 3

    def fzero_reference(
        self, t: np.ndarray, x: np.ndarray, quadrature_order: int = 80
    ) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        nodes, weights = np.polynomial.hermite.hermgauss(quadrature_order)
        scalar = state @ self.w
        h = np.maximum(self.T - time, 0.0)
        shifted = scalar[:, None] + 2.0 * np.sqrt(h)[:, None] * nodes[None, :]
        terminal = terminal_profile(shifted, A=self.A, beta=self.beta)
        return np.sum(terminal * weights[None, :], axis=1) / math.sqrt(math.pi)

    @property
    def segment_coefficients(self) -> tuple[float, float]:
        return (-self.A, self.A)

    @property
    def segment_endpoints(self) -> tuple[np.ndarray, np.ndarray]:
        return -self.sigma * self.A * self.w, self.sigma * self.A * self.w

    @property
    def box_low(self) -> np.ndarray:
        first, second = self.segment_endpoints
        return np.minimum(first, second)

    @property
    def box_high(self) -> np.ndarray:
        first, second = self.segment_endpoints
        return np.maximum(first, second)

    @property
    def z_ball_radius(self) -> float:
        return self.sigma * self.A

    @property
    def z_center(self) -> np.ndarray:
        return np.zeros(self.d, dtype=np.float64)

    @property
    def dose_scale(self) -> float:
        return self.sigma * self.A

    @property
    def picard_proxy(self) -> float:
        return 3.0 * self.kappa * self.A**2 * self.T

    def certificate_labels(self) -> dict[str, str]:
        return {
            "segment": "PDE-derived: Lipschitz preservation |psi_s|<=2",
            "box": "PDE-derived: coordinatewise hull of the segment",
            "ball": "PDE-derived: Euclidean hull of the segment",
            "span": "PDE-derived: translation invariance orthogonal to w",
        }

    def certificate_parameters(self) -> dict[str, Any]:
        _, metadata = _load_c3_spline(self.reference_path)
        return {
            "w_seed": self.direction_seed,
            "w_sha256": hashlib.sha256(self.w.tobytes()).hexdigest(),
            "segment_coefficients": list(self.segment_coefficients),
            "reference_path": Path(self.reference_path).as_posix(),
            "reference_metadata": metadata,
        }

    def reference_residual(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        tau, scalar = self._tau_s(t, x)
        spline, _ = _load_c3_spline(self.reference_path)
        flat_tau = tau.reshape(-1)
        flat_s = scalar.reshape(-1)
        psi_tau = spline.ev(flat_tau, flat_s, dx=1, dy=0)
        psi_s = spline.ev(flat_tau, flat_s, dx=0, dy=1)
        psi_ss = spline.ev(flat_tau, flat_s, dx=0, dy=2)
        residual = -psi_tau + psi_ss - self.kappa * np.abs(psi_s) ** 3
        return residual.reshape(np.shape(tau))


@dataclass(frozen=True)
class LQGEquation(ExactEquation):
    d: int
    gamma: float = 1.0
    T: float = 0.5
    sigma: float = SQRT2
    name: str = "LQG_exploratory"
    family: str = "lqg"

    def __post_init__(self) -> None:
        object.__setattr__(self, "mu", 0.0)

    def _p(self, t: np.ndarray | float) -> np.ndarray:
        return self.gamma / (1.0 + 4.0 * self.gamma * (self.T - np.asarray(t, dtype=np.float64)))

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        h = self.T - time
        denominator = 1.0 + 4.0 * self.gamma * h
        return 0.5 * self.d * np.log(denominator) + self.gamma * np.sum(state**2, axis=-1) / denominator

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        return 2.0 * self.sigma * self._p(time)[..., None] * state

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return self.gamma * np.sum(np.asarray(x, dtype=np.float64) ** 2, axis=-1)

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        del u
        vector = np.asarray(z, dtype=np.float64)
        return -0.5 * np.sum(vector * vector, axis=-1)

    def fzero_reference(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        return self.gamma * (
            np.sum(state**2, axis=-1) + 2.0 * self.d * (self.T - time)
        )

    @property
    def p_interval(self) -> tuple[float, float]:
        return (
            self.gamma / (1.0 + 4.0 * self.gamma * self.T),
            self.gamma,
        )

    @property
    def z_center(self) -> np.ndarray:
        # Only used by static diagnostic fallbacks; transforms use the
        # state-dependent Riccati centre.
        return np.zeros(self.d, dtype=np.float64)

    @property
    def z_ball_radius(self) -> float:
        return 2.0 * self.sigma * self.gamma

    @property
    def box_low(self) -> np.ndarray:
        return np.full(self.d, -self.z_ball_radius, dtype=np.float64)

    @property
    def box_high(self) -> np.ndarray:
        return np.full(self.d, self.z_ball_radius, dtype=np.float64)

    @property
    def dose_scale(self) -> float:
        return self.z_ball_radius

    @property
    def picard_proxy(self) -> float:
        return 4.0 * self.gamma * self.T

    def dynamic_endpoints(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        state = np.asarray(x, dtype=np.float64)
        low, high = self.p_interval
        return (
            2.0 * self.sigma * low * state,
            2.0 * self.sigma * high * state,
        )

    def analytic_residual(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        h = self.T - time
        denominator = 1.0 + 4.0 * self.gamma * h
        p = self.gamma / denominator
        u_t = -2.0 * self.d * p + 4.0 * p**2 * np.sum(state**2, axis=-1)
        laplacian = 2.0 * self.d * p
        z = 2.0 * self.sigma * p[..., None] * state
        return u_t + laplacian + self.generator(np.zeros_like(p), z)

    def certificate_labels(self) -> dict[str, str]:
        return {
            "segment": "Riccati-derived: z=2*sqrt(2)*p*x over the monotone p interval",
            "box": "Riccati-derived: coordinatewise hull of the segment",
            "ball": "Riccati-derived: Euclidean hull using |x|<=1",
            "span": "Riccati-derived: z lies in span(x)",
        }

    def certificate_parameters(self) -> dict[str, Any]:
        return {
            "gamma": self.gamma,
            "p_interval": list(self.p_interval),
            "test_support": "x=r*theta, r~Uniform[0,1]",
        }


SPECS = {
    "P1": {"kind": "p1", "primary": "box"},
    "P4": {"kind": "p4", "primary": "tight_segment"},
    "C1": {"kind": "c1", "primary": "segment"},
    "C2-convex": {
        "kind": "c2", "a": 2.0, "b": 0.0, "frac_A": 0.5,
        "primary": "box",
    },
    "C2-cancel": {
        "kind": "c2", "a": 3.0, "b": 1.0, "frac_A": 0.25,
        "primary": "box",
    },
    "C2-flip": {
        "kind": "c2", "a": 3.0, "b": 1.0, "frac_A": 0.125,
        "primary": "box",
    },
    "C3": {"kind": "c3", "primary": "segment"},
    "LQG": {"kind": "lqg", "primary": "segment"},
}


def make_equation(pde_id: str, d: int) -> ExactEquation:
    spec = SPECS[pde_id]
    kind = spec["kind"]
    if kind == "p1":
        return RidgeLSEHJB(d=d, direction_seed=BASE_SEED)
    if kind == "p4":
        return NormDriverHJB(
            d=d, reference_path=str(P4_REFERENCE), T=0.5,
            direction_seed=BASE_SEED,
        )
    if kind == "c1":
        return L1ControlHJB(
            d=d, reference_path=str(P4_REFERENCE), T=0.5,
            direction_seed=BASE_SEED,
        )
    if kind == "c2":
        return LQGameHJB(
            d=d,
            a=float(spec["a"]),
            b=float(spec["b"]),
            frac_A=float(spec["frac_A"]),
            direction_seed=BASE_SEED,
            name=pde_id,
        )
    if kind == "c3":
        return ProductionCubicHJB(d=d)
    if kind == "lqg":
        return LQGEquation(d=d)
    raise ValueError(pde_id)
