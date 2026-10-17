"""Smooth nonconvex--nonconcave drivers used by the experiment."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
    dense_orthonormal_basis,
)


class SineQuadraticDriver(ConcaveSemilinearDriver):
    r"""Driver lambda*sin(2*pi*u) - |z + eta*u*e|^2/2.

    ``eta=0`` is the separable primary benchmark.  ``eta>0`` is the genuinely
    nonseparable benchmark.  The inherited tangent evaluator is deliberately
    used without a PDE-specific solver.
    """

    has_gradient = True
    r_low, r_high = 0.0, 1.0

    def __init__(
        self,
        *,
        d: int = 100,
        T: float = 0.5,
        lam: float = 0.2,
        eta: float = 0.0,
        basis_seed: int = 61511,
    ) -> None:
        if d < 2:
            raise ValueError("the two-direction benchmark needs d >= 2")
        if lam < 0.0 or eta < 0.0:
            raise ValueError((lam, eta))
        self.d, self.k, self.T = int(d), 2, float(T)
        self.lam, self.eta = float(lam), float(eta)
        self.basis_seed = int(basis_seed)
        self.B = dense_orthonormal_basis(d, 2, seed=basis_seed)
        self.beta_1, self.beta_2 = self.B[:, 0], self.B[:, 1]
        self.e = self.beta_1
        self.family = "sine_quadratic" if eta == 0.0 else "sine_quadratic_coupled"
        lam_tag = str(lam).replace(".", "p")
        eta_tag = str(eta).replace(".", "p")
        self.name = f"{self.family}_d{d}_k2_lam{lam_tag}_eta{eta_tag}"

    def _basis(self, x: torch.Tensor) -> torch.Tensor:
        return torch.as_tensor(self.B, dtype=x.dtype, device=x.device)

    def _e(self, z: torch.Tensor) -> torch.Tensor:
        return torch.as_tensor(self.e, dtype=z.dtype, device=z.device)

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        coordinates = x @ self._basis(x)
        return 0.5 + 0.3 * torch.sin(coordinates[..., 0]) * torch.cos(
            coordinates[..., 1]
        )

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return self.terminal_from_coordinates_numpy(np.asarray(x) @ self.B)

    def terminal_from_coordinates_numpy(self, s: np.ndarray) -> np.ndarray:
        coordinates = np.asarray(s)
        return 0.5 + 0.3 * np.sin(coordinates[..., 0]) * np.cos(
            coordinates[..., 1]
        )

    def terminal_coordinate_gradient_numpy(self, s: np.ndarray) -> np.ndarray:
        coordinates = np.asarray(s)
        return np.stack(
            (
                0.3 * np.cos(coordinates[..., 0]) * np.cos(coordinates[..., 1]),
                -0.3 * np.sin(coordinates[..., 0]) * np.sin(coordinates[..., 1]),
            ),
            axis=-1,
        )

    def terminal_gradient_torch(self, x: torch.Tensor) -> torch.Tensor:
        basis = self._basis(x)
        coordinates = x @ basis
        coefficients = torch.stack(
            (
                0.3 * torch.cos(coordinates[..., 0]) * torch.cos(
                    coordinates[..., 1]
                ),
                -0.3 * torch.sin(coordinates[..., 0]) * torch.sin(
                    coordinates[..., 1]
                ),
            ),
            dim=-1,
        )
        return coefficients @ basis.T

    def _shift(self, y: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        if self.eta == 0.0:
            return z
        return z + self.eta * y[..., None] * self._e(z)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        if z is None:
            raise ValueError("SineQuadraticDriver requires z")
        shifted = self._shift(y, z)
        return self.lam * torch.sin(2.0 * math.pi * y) - 0.5 * shifted.square().sum(
            dim=-1
        )

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        if z is None:
            raise ValueError("SineQuadraticDriver requires z")
        result = 2.0 * math.pi * self.lam * torch.cos(2.0 * math.pi * y)
        if self.eta:
            result = result - self.eta * (self._shift(y, z) * self._e(z)).sum(
                dim=-1
            )
        return result

    def driver_z(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        if z is None:
            raise ValueError("SineQuadraticDriver requires z")
        return -self._shift(y, z)

    def driver_uu(self, y: torch.Tensor) -> torch.Tensor:
        return -4.0 * math.pi**2 * self.lam * torch.sin(2.0 * math.pi * y) - self.eta**2

    def driver_reduced_numpy(self, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        reduced_z = np.asarray(z, dtype=np.float64)
        shifted = reduced_z.copy()
        if self.eta:
            shifted[..., 0] += self.eta * np.asarray(y)
        return self.lam * np.sin(2.0 * np.pi * np.asarray(y)) - 0.5 * np.square(
            shifted
        ).sum(axis=-1)

    def reduced_joint_hessian_eigenvalues(self, y: np.ndarray) -> np.ndarray:
        """Eigenvalues in (u,z1,z2); ambient orthogonal z modes equal -1."""
        values = np.asarray(y, dtype=np.float64).reshape(-1)
        f_uu = -4.0 * np.pi**2 * self.lam * np.sin(2.0 * np.pi * values) - self.eta**2
        matrices = np.zeros((len(values), 3, 3), dtype=np.float64)
        matrices[:, 0, 0] = f_uu
        matrices[:, 1, 1] = -1.0
        matrices[:, 2, 2] = -1.0
        matrices[:, 0, 1] = matrices[:, 1, 0] = -self.eta
        return np.linalg.eigvalsh(matrices)

    def barrier_audit(self) -> dict[str, Any]:
        zero = float(self.driver_reduced_numpy(np.asarray(0.0), np.zeros(2)))
        one = float(self.driver_reduced_numpy(np.asarray(1.0), np.zeros(2)))
        return {
            "forward_tau_equation": "w_tau=Delta w+f(w,sqrt(2) grad w)",
            "lower_constant": 0.0,
            "lower_driver": zero,
            "lower_is_subsolution": bool(zero >= -1e-14),
            "upper_constant": 1.0,
            "upper_driver": one,
            "upper_is_supersolution": bool(one <= 1e-14),
            "terminal_min": 0.2,
            "terminal_max": 0.8,
            "invariant_interval_verified_analytically": bool(
                zero >= -1e-14 and one <= 1e-14
            ),
        }

    def metadata(self) -> dict[str, Any]:
        result = super().metadata()
        result.update(
            {
                "lambda": self.lam,
                "eta": self.eta,
                "e_equals_beta_1": True,
                "terminal": "0.5+0.3*sin(beta1.x)*cos(beta2.x)",
                "value_curvature_changes_sign": self.lam > 0.0,
                "mixed_derivative": "-eta*e",
                "invariant_interval": [0.0, 1.0],
            }
        )
        return result


def make_driver(
    *, d: int = 100, lam: float = 0.2, coupled: bool = False, eta: float = 1.0
) -> SineQuadraticDriver:
    return SineQuadraticDriver(d=d, lam=lam, eta=eta if coupled else 0.0)


__all__ = ["SineQuadraticDriver", "make_driver"]
