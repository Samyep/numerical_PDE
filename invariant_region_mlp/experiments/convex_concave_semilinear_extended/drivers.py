"""Generic concave drivers for nonseparable and effective-dimension tests."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import torch


def dense_orthonormal_basis(
    d: int, k: int, *, seed: int = 24017, maximum_k: int = 20
) -> np.ndarray:
    """Return a deterministic dense prefix so k comparisons share directions."""
    draw_k = min(maximum_k, d)
    if not 1 <= k <= draw_k:
        raise ValueError((d, k, maximum_k))
    rng = np.random.default_rng(seed + 1009 * d)
    raw = rng.standard_normal((d, draw_k))
    basis, triangular = np.linalg.qr(raw, mode="reduced")
    signs = np.where(np.diag(triangular) < 0.0, -1.0, 1.0)
    basis = basis * signs[None, :]
    return np.asarray(basis[:, :k], dtype=np.float64)


class ConcaveSemilinearDriver:
    name: str
    family: str
    d: int
    k: int
    T: float
    has_gradient: bool
    r_low: float
    r_high: float
    reconstruction_sign: float = 1.0
    sigma: float = math.sqrt(2.0)

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def terminal_from_coordinates_numpy(self, s: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def terminal_coordinate_gradient_numpy(self, s: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def terminal_gradient_numpy(self, x: np.ndarray) -> np.ndarray:
        coordinates = np.asarray(x) @ self.B
        return self.terminal_coordinate_gradient_numpy(coordinates) @ self.B.T

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        raise NotImplementedError

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        raise NotImplementedError

    def driver_z(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        if z is None:
            return torch.zeros((*y.shape, self.d), dtype=y.dtype, device=y.device)
        return torch.zeros_like(z)

    def driver_reduced_numpy(self, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def project_r(self, value: torch.Tensor) -> torch.Tensor:
        return torch.clamp(value, self.r_low, self.r_high)

    def constrain_r(self, raw: torch.Tensor) -> torch.Tensor:
        return self.r_low + (self.r_high - self.r_low) * torch.sigmoid(raw)

    def tangent_coefficients(
        self, r: torch.Tensor, q: torch.Tensor | None
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """The single generic tangent calculation used by every experiment."""
        c = self.driver_y(r, q)
        b = self.driver_z(r, q)
        f = self.driver(r, q)
        q_term = torch.zeros_like(r) if q is None else (b * q).sum(dim=-1)
        h = f - c * r - q_term
        return c, b, h

    def tangent_coefficients_autodiff(
        self, r: torch.Tensor, q: torch.Tensor | None
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        r_var = r.detach().clone().requires_grad_(True)
        if self.has_gradient:
            if q is None:
                raise ValueError("gradient-dependent driver needs q")
            q_var = q.detach().clone().requires_grad_(True)
        else:
            q_var = None
        value = self.driver(r_var, q_var)
        derivatives = torch.autograd.grad(
            value.sum(), (r_var,) if q_var is None else (r_var, q_var)
        )
        c = derivatives[0]
        b = (
            derivatives[1]
            if q_var is not None
            else torch.zeros((*r.shape, self.d), dtype=r.dtype, device=r.device)
        )
        h = value.detach() - c.detach() * r
        if q is not None:
            h = h - (b.detach() * q).sum(dim=-1)
        return c.detach(), b.detach(), h.detach()

    def working_to_original(self, value: np.ndarray) -> np.ndarray:
        return self.reconstruction_sign * np.asarray(value)

    def bregman_divergence(
        self,
        y: torch.Tensor,
        z: torch.Tensor | None,
        r: torch.Tensor,
        q: torch.Tensor | None,
    ) -> torch.Tensor:
        c, b, _ = self.tangent_coefficients(r, q)
        result = -self.driver(y, z) + self.driver(r, q) + c * (y - r)
        if self.has_gradient:
            assert z is not None and q is not None
            result = result + (b * (z - q)).sum(dim=-1)
        return result

    def metadata(self) -> dict[str, Any]:
        gram = self.B.T @ self.B
        return {
            "name": self.name,
            "family": self.family,
            "d": self.d,
            "k": self.k,
            "T": self.T,
            "sigma": self.sigma,
            "has_gradient": self.has_gradient,
            "tangent_domain": [self.r_low, self.r_high],
            "reconstruction_sign": self.reconstruction_sign,
            "basis_seed": self.basis_seed,
            "basis_orthogonality_max_abs": float(np.max(np.abs(gram - np.eye(self.k)))),
            "actor_input": "full ambient x",
            "oracle_projection_used": False,
        }


class MultiDirectionTerminalDriver(ConcaveSemilinearDriver):
    amplitude: float

    def __init__(
        self,
        *,
        d: int,
        k: int,
        T: float,
        basis_seed: int,
        family: str,
    ) -> None:
        self.d, self.k, self.T = d, k, T
        self.basis_seed, self.family = basis_seed, family
        self.name = f"{family}_d{d}_k{k}"
        self.B = dense_orthonormal_basis(d, k, seed=basis_seed)
        self.beta = self.B[:, 0]
        self.amplitude = 0.3 / k

    def terminal_original(self, x: torch.Tensor) -> torch.Tensor:
        basis = torch.as_tensor(self.B, dtype=x.dtype, device=x.device)
        return 0.5 + self.amplitude * torch.sin(x @ basis).sum(dim=-1)

    def terminal_original_numpy(self, x: np.ndarray) -> np.ndarray:
        return 0.5 + self.amplitude * np.sin(np.asarray(x) @ self.B).sum(axis=-1)

    def terminal_from_coordinates_original_numpy(self, s: np.ndarray) -> np.ndarray:
        return 0.5 + self.amplitude * np.sin(np.asarray(s)).sum(axis=-1)

    def terminal_coordinate_gradient_original_numpy(self, s: np.ndarray) -> np.ndarray:
        return self.amplitude * np.cos(np.asarray(s))


class AllenCahnDriver(MultiDirectionTerminalDriver):
    has_gradient = False
    r_low, r_high = 0.0, 1.0

    def __init__(self, *, d: int, k: int, T: float = 0.5, basis_seed: int = 24017) -> None:
        super().__init__(d=d, k=k, T=T, basis_seed=basis_seed, family="allen_cahn")

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.terminal_original(x)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return self.terminal_original_numpy(x)

    def terminal_from_coordinates_numpy(self, s: np.ndarray) -> np.ndarray:
        return self.terminal_from_coordinates_original_numpy(s)

    def terminal_coordinate_gradient_numpy(self, s: np.ndarray) -> np.ndarray:
        return self.terminal_coordinate_gradient_original_numpy(s)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del z
        return y - y.pow(3)

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del z
        return 1.0 - 3.0 * y.square()

    def driver_reduced_numpy(self, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        del z
        return y - y**3


class JointConcaveDriver(MultiDirectionTerminalDriver):
    has_gradient = True
    r_low, r_high = 0.0, 1.0

    def __init__(self, *, d: int, k: int, T: float = 0.5, basis_seed: int = 24017) -> None:
        super().__init__(d=d, k=k, T=T, basis_seed=basis_seed, family="joint_concave")

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.terminal_original(x)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return self.terminal_original_numpy(x)

    def terminal_from_coordinates_numpy(self, s: np.ndarray) -> np.ndarray:
        return self.terminal_from_coordinates_original_numpy(s)

    def terminal_coordinate_gradient_numpy(self, s: np.ndarray) -> np.ndarray:
        return self.terminal_coordinate_gradient_original_numpy(s)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        return y - y.pow(3) - 0.5 * z.square().sum(dim=-1)

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del z
        return 1.0 - 3.0 * y.square()

    def driver_z(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del y
        assert z is not None
        return -z

    def driver_reduced_numpy(self, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        return y - y**3 - 0.5 * np.square(z).sum(axis=-1)


class JointConvexOriginalDriver(MultiDirectionTerminalDriver):
    has_gradient = True
    r_low, r_high = 0.0, 1.0

    def __init__(self, *, d: int, k: int, T: float = 0.5, basis_seed: int = 24017) -> None:
        super().__init__(d=d, k=k, T=T, basis_seed=basis_seed, family="joint_convex_original")

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.terminal_original(x)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return self.terminal_original_numpy(x)

    def terminal_from_coordinates_numpy(self, s: np.ndarray) -> np.ndarray:
        return self.terminal_from_coordinates_original_numpy(s)

    def terminal_coordinate_gradient_numpy(self, s: np.ndarray) -> np.ndarray:
        return self.terminal_coordinate_gradient_original_numpy(s)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        return y.pow(3) - y + 0.5 * z.square().sum(dim=-1)

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del z
        return 3.0 * y.square() - 1.0

    def driver_z(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del y
        assert z is not None
        return z

    def driver_reduced_numpy(self, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        return y**3 - y + 0.5 * np.square(z).sum(axis=-1)


class SignTransformedDriver(ConcaveSemilinearDriver):
    """Map a convex f to the concave working driver -f(-w,-z)."""

    reconstruction_sign = -1.0

    def __init__(self, original: JointConvexOriginalDriver) -> None:
        self.original = original
        self.family = "joint_convex"
        self.name = f"{self.family}_d{original.d}_k{original.k}"
        self.d, self.k, self.T = original.d, original.k, original.T
        self.B, self.beta = original.B, original.beta
        self.basis_seed = original.basis_seed
        self.has_gradient = original.has_gradient
        self.r_low, self.r_high = -original.r_high, -original.r_low

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return -self.original.terminal(x)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return -self.original.terminal_numpy(x)

    def terminal_from_coordinates_numpy(self, s: np.ndarray) -> np.ndarray:
        return -self.original.terminal_from_coordinates_numpy(s)

    def terminal_coordinate_gradient_numpy(self, s: np.ndarray) -> np.ndarray:
        return -self.original.terminal_coordinate_gradient_numpy(s)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        return -self.original.driver(-y, -z)

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        return self.original.driver_y(-y, -z)

    def driver_z(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        return self.original.driver_z(-y, -z)

    def driver_reduced_numpy(self, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        return -self.original.driver_reduced_numpy(-y, -z)


class NonseparableDriver(ConcaveSemilinearDriver):
    """f(u,z)=u-u^3-0.5|z+eta*u*v|^2 with v=beta_1."""

    has_gradient = True
    r_low, r_high = 0.0, 1.0

    def __init__(
        self, *, d: int = 100, T: float = 0.5, eta: float = 1.0, basis_seed: int = 73129
    ) -> None:
        self.family, self.d, self.k, self.T = "nonseparable", d, 2, T
        self.name = f"nonseparable_d{d}_k2"
        self.eta, self.basis_seed = eta, basis_seed
        self.B = dense_orthonormal_basis(d, 2, seed=basis_seed)
        self.beta, self.v = self.B[:, 0], self.B[:, 0]

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        basis = torch.as_tensor(self.B, dtype=x.dtype, device=x.device)
        s = x @ basis
        return 0.5 + 0.15 * torch.sin(s[..., 0]) + 0.15 * torch.cos(s[..., 1])

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return self.terminal_from_coordinates_numpy(np.asarray(x) @ self.B)

    def terminal_from_coordinates_numpy(self, s: np.ndarray) -> np.ndarray:
        value = np.asarray(s)
        return 0.5 + 0.15 * np.sin(value[..., 0]) + 0.15 * np.cos(value[..., 1])

    def terminal_coordinate_gradient_numpy(self, s: np.ndarray) -> np.ndarray:
        value = np.asarray(s)
        return np.stack((0.15 * np.cos(value[..., 0]), -0.15 * np.sin(value[..., 1])), axis=-1)

    def _shift(self, y: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        v = torch.as_tensor(self.v, dtype=z.dtype, device=z.device)
        return z + self.eta * y[..., None] * v

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        shifted = self._shift(y, z)
        return y - y.pow(3) - 0.5 * shifted.square().sum(dim=-1)

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        shifted = self._shift(y, z)
        v = torch.as_tensor(self.v, dtype=z.dtype, device=z.device)
        return 1.0 - 3.0 * y.square() - self.eta * (shifted * v).sum(dim=-1)

    def driver_z(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        return -self._shift(y, z)

    def driver_reduced_numpy(self, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        value = np.asarray(z)
        shifted = value.copy()
        shifted[..., 0] += self.eta * y
        return y - y**3 - 0.5 * np.square(shifted).sum(axis=-1)

    def metadata(self) -> dict[str, Any]:
        result = super().metadata()
        result.update({"eta": self.eta, "v_equals_beta_1": True, "mixed_f_uz": "-eta*v"})
        return result


def make_driver(
    family: str, *, d: int = 100, k: int = 2, T: float = 0.5
) -> ConcaveSemilinearDriver:
    if family == "allen_cahn":
        return AllenCahnDriver(d=d, k=k, T=T)
    if family == "joint_concave":
        return JointConcaveDriver(d=d, k=k, T=T)
    if family == "joint_convex":
        return SignTransformedDriver(JointConvexOriginalDriver(d=d, k=k, T=T))
    if family == "nonseparable":
        if k != 2:
            raise ValueError("nonseparable pilot is defined for k=2")
        return NonseparableDriver(d=d, T=T)
    raise KeyError(family)


def terminal_signal_diagnostics(
    driver: ConcaveSemilinearDriver, *, samples: int = 8192, seed: int = 991
) -> dict[str, Any]:
    """Audit terminal signal size and ambient gradient span without training labels."""
    rng = np.random.default_rng(seed + driver.d + 31 * driver.k)
    x = rng.standard_normal((samples, driver.d))
    coordinates = x @ driver.B
    values = driver.working_to_original(driver.terminal_numpy(x))
    active_coefficients = driver.working_to_original(
        driver.terminal_coordinate_gradient_numpy(coordinates)
    )
    gradients = active_coefficients @ driver.B.T
    # Right multiplication by the orthonormal B.T preserves all nonzero
    # singular values, so the smaller coefficient matrix gives the exact
    # ambient gradient rank without a costly samples-by-d SVD.
    singular = np.linalg.svd(active_coefficients, compute_uv=False)
    tolerance = singular[0] * max(active_coefficients.shape) * np.finfo(np.float64).eps
    rank = int(np.sum(singular > tolerance))
    return {
        "samples": samples,
        "terminal_mean": float(np.mean(values)),
        "terminal_variance": float(np.var(values)),
        "terminal_std": float(np.std(values)),
        "gradient_energy": float(np.mean(np.square(gradients).sum(axis=1))),
        "gradient_norm_mean": float(np.mean(np.linalg.norm(gradients, axis=1))),
        "ambient_gradient_matrix_rank": rank,
        "expected_rank": driver.k,
        "leading_gradient_singular_values": singular[: driver.k].tolist(),
        "smallest_to_largest_active_singular_ratio": float(singular[driver.k - 1] / singular[0]),
        "basis_max_abs_entry": float(np.max(np.abs(driver.B))),
        "basis_median_abs_entry": float(np.median(np.abs(driver.B))),
    }
