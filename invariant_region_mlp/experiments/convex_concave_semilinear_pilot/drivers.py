"""Generic concave drivers and the reusable convex-to-concave transform."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import torch


class ConcaveSemilinearDriver:
    name: str
    d: int
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

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        raise NotImplementedError

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        raise NotImplementedError

    def driver_z(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        if z is None:
            return torch.zeros((*y.shape, self.d), dtype=y.dtype, device=y.device)
        return torch.zeros_like(z)

    def driver_numpy(self, y: np.ndarray, z_scalar: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def project_r(self, value: torch.Tensor) -> torch.Tensor:
        return torch.clamp(value, self.r_low, self.r_high)

    def constrain_r(self, raw: torch.Tensor) -> torch.Tensor:
        return self.r_low + (self.r_high - self.r_low) * torch.sigmoid(raw)

    def tangent_coefficients(
        self, r: torch.Tensor, q: torch.Tensor | None
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
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
        inputs = (r_var,) if q_var is None else (r_var, q_var)
        derivatives = torch.autograd.grad(value.sum(), inputs, create_graph=False)
        c = derivatives[0]
        b = (
            derivatives[1]
            if q_var is not None
            else torch.zeros((*r.shape, self.d), dtype=r.dtype, device=r.device)
        )
        h = value.detach() - c.detach() * r - (
            0.0 if q is None else (b.detach() * q).sum(dim=-1)
        )
        return c.detach(), b.detach(), h

    def working_to_original(self, value: np.ndarray) -> np.ndarray:
        return self.reconstruction_sign * np.asarray(value)

    def bregman_divergence(
        self,
        y: torch.Tensor,
        z: torch.Tensor | None,
        r: torch.Tensor,
        q: torch.Tensor | None,
    ) -> torch.Tensor:
        """D_{-f}((y,z),(r,q)); nonnegative on the tangent domain."""
        c, b, _ = self.tangent_coefficients(r, q)
        result = -self.driver(y, z) + self.driver(r, q) + c * (y - r)
        if self.has_gradient:
            assert z is not None and q is not None
            result = result + (b * (z - q)).sum(dim=-1)
        return result

    def metadata(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "d": self.d,
            "T": self.T,
            "sigma": self.sigma,
            "has_gradient": self.has_gradient,
            "tangent_domain": [self.r_low, self.r_high],
            "reconstruction_sign": self.reconstruction_sign,
            "beta_seed": getattr(self, "beta_seed", None),
            "beta_norm": float(np.linalg.norm(self.beta)),
        }


class PeriodicTerminalDriver(ConcaveSemilinearDriver):
    def __init__(self, *, d: int, T: float, beta_seed: int, name: str) -> None:
        self.d, self.T, self.beta_seed, self.name = d, T, beta_seed, name
        rng = np.random.default_rng(beta_seed)
        beta = rng.standard_normal(d)
        self.beta = beta / np.linalg.norm(beta)

    def terminal_original(self, x: torch.Tensor) -> torch.Tensor:
        beta = torch.as_tensor(self.beta, dtype=x.dtype, device=x.device)
        return 0.5 + 0.25 * torch.sin(x @ beta)

    def terminal_original_numpy(self, x: np.ndarray) -> np.ndarray:
        return 0.5 + 0.25 * np.sin(np.asarray(x) @ self.beta)


class AllenCahnDriver(PeriodicTerminalDriver):
    has_gradient = False
    r_low, r_high = 0.0, 1.0

    def __init__(self, *, d: int = 100, T: float = 0.5, beta_seed: int = 24017) -> None:
        super().__init__(d=d, T=T, beta_seed=beta_seed, name="allen_cahn")

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.terminal_original(x)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return self.terminal_original_numpy(x)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del z
        return y - y.pow(3)

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        del z
        return 1.0 - 3.0 * y.square()

    def driver_numpy(self, y: np.ndarray, z_scalar: np.ndarray) -> np.ndarray:
        del z_scalar
        return y - y**3


class JointConcaveDriver(PeriodicTerminalDriver):
    has_gradient = True
    r_low, r_high = 0.0, 1.0

    def __init__(self, *, d: int = 100, T: float = 0.5, beta_seed: int = 24017) -> None:
        super().__init__(d=d, T=T, beta_seed=beta_seed, name="joint_concave")

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.terminal_original(x)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return self.terminal_original_numpy(x)

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

    def driver_numpy(self, y: np.ndarray, z_scalar: np.ndarray) -> np.ndarray:
        return y - y**3 - 0.5 * z_scalar**2


class JointConvexOriginalDriver(PeriodicTerminalDriver):
    """Original convex driver, wrapped below before use by the solver."""

    has_gradient = True
    r_low, r_high = 0.0, 1.0

    def __init__(self, *, d: int = 100, T: float = 0.5, beta_seed: int = 24017) -> None:
        super().__init__(d=d, T=T, beta_seed=beta_seed, name="joint_convex_original")

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.terminal_original(x)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return self.terminal_original_numpy(x)

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

    def driver_numpy(self, y: np.ndarray, z_scalar: np.ndarray) -> np.ndarray:
        return y**3 - y + 0.5 * z_scalar**2


class SignTransformedDriver(ConcaveSemilinearDriver):
    """Map a convex f to the concave driver -f(-w,-z)."""

    reconstruction_sign = -1.0

    def __init__(self, original: ConcaveSemilinearDriver) -> None:
        self.original = original
        self.name = "joint_convex"
        self.d, self.T = original.d, original.T
        self.beta, self.beta_seed = original.beta, original.beta_seed
        self.has_gradient = original.has_gradient
        self.r_low, self.r_high = -original.r_high, -original.r_low

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return -self.original.terminal(x)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return -self.original.terminal_numpy(x)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        return -self.original.driver(-y, -z)

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        return self.original.driver_y(-y, -z)

    def driver_z(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        assert z is not None
        return self.original.driver_z(-y, -z)

    def driver_numpy(self, y: np.ndarray, z_scalar: np.ndarray) -> np.ndarray:
        return -self.original.driver_numpy(-y, -z_scalar)


def make_driver(name: str, *, d: int = 100) -> ConcaveSemilinearDriver:
    if name == "allen_cahn":
        return AllenCahnDriver(d=d)
    if name == "joint_concave":
        return JointConcaveDriver(d=d)
    if name == "joint_convex":
        return SignTransformedDriver(JointConvexOriginalDriver(d=d))
    raise KeyError(name)
