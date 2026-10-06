"""Transverse-aware conservative face models for the 2-D Euler ablations.

Both variants use an oriented 3-by-6 primitive-variable patch for every face.
The ``flat18`` model is deliberately unstructured: it must learn for itself
which transverse information matters.  The ``gated18`` model separates the
existing normal six-cell mapping from a transverse residual and reduces
exactly to its normal branch when the three rows agree.

``central_nonnegative18`` is the direct 2-D extension of the retained 1-D
HCFL method: central physical flux plus entropy-fixed Roe dissipation whose
four wave multipliers are constrained to ``(0, 2)``.  It retains the same
oriented 18-cell input and the same shared network for x and y faces.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import torch
from torch import nn


BASE_EXPERIMENT = Path(__file__).resolve().parents[1] / "euler_2d_clawpack_quadrants"
sys.path.insert(0, str(BASE_EXPERIMENT))
import euler_2d_common as C  # noqa: E402


TRANSVERSE_VARIANTS = ("flat18", "gated18")
MODEL_TYPES = ("normal6wide", *TRANSVERSE_VARIANTS)
CONSISTENT_VARIANT = "central_nonnegative18"
NORMAL6_CONSISTENT_VARIANT = "central_nonnegative6"
# Backward-compatible name used by the first running flat18 process.
VARIANTS = TRANSVERSE_VARIANTS
NORMAL_SHIFTS = (-3, -2, -1, 0, 1, 2)


def _mlp(input_features: int, width: int) -> nn.Sequential:
    network = nn.Sequential(
        nn.Linear(input_features, width),
        nn.Tanh(),
        nn.Linear(width, width),
        nn.Tanh(),
        nn.Linear(width, 4),
    )
    nn.init.zeros_(network[-1].weight)
    nn.init.zeros_(network[-1].bias)
    return network


def oriented_primitive_patch(state: torch.Tensor) -> torch.Tensor:
    """Return ``[..., transverse, face, 3, 6, primitive]`` patches.

    The first patch axis is ordered ``lower, centre, upper`` in the oriented
    transverse coordinate.  Constant-extrapolation padding matches the
    boundary condition used to generate the PyClaw data.
    """

    values = C.primitive(state)
    transverse_cells = values.shape[-3]
    normal_cells = values.shape[-2]
    interfaces = torch.arange(1, normal_cells, device=state.device)
    offsets = torch.tensor(NORMAL_SHIFTS, device=state.device)
    normal_indices = (interfaces[:, None] + offsets[None, :]).clamp(
        0, normal_cells - 1
    )
    normal = values.index_select(-2, normal_indices.reshape(-1)).reshape(
        *values.shape[:-3],
        transverse_cells,
        normal_cells - 1,
        len(NORMAL_SHIFTS),
        4,
    )
    centres = torch.arange(transverse_cells, device=state.device)
    rows = [
        normal.index_select(-4, (centres + shift).clamp(0, transverse_cells - 1))
        for shift in (-1, 0, 1)
    ]
    return torch.stack(rows, dim=-3)


class DirectionalPatchHLLCRoeCorrection(nn.Module):
    """HLLC plus signed Roe correction conditioned on an 18-cell patch."""

    def __init__(
        self,
        mean: np.ndarray,
        std: np.ndarray,
        variant: str,
        flat_width: int = 72,
        normal_width: int = 72,
        transverse_width: int = 32,
    ) -> None:
        super().__init__()
        if variant not in TRANSVERSE_VARIANTS:
            raise ValueError(f"Unknown transverse variant: {variant}")
        self.variant = variant
        if variant == "flat18":
            self.net = _mlp(3 * 6 * 4, flat_width)
        else:
            self.normal_net = _mlp(6 * 4, normal_width)
            self.transverse_net = _mlp(3 * 6 * 4, transverse_width)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def normalized_patch(self, state: torch.Tensor) -> torch.Tensor:
        patch = oriented_primitive_patch(state)
        return (patch - self.mean) / self.std

    def coefficient_details(
        self, state: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Return coefficients, centre-only coefficients, and gate values.

        For ``flat18``, centre-only means replacing both transverse rows by
        the centre row before applying the same network.  It measures whether
        the unstructured model learned to use or ignore transverse context.
        """

        patch = self.normalized_patch(state)
        centre = patch[..., 1, :, :]
        if self.variant == "flat18":
            features = patch.flatten(start_dim=-3)
            repeated = torch.stack((centre, centre, centre), dim=-3)
            centre_only = torch.tanh(
                self.net(repeated.flatten(start_dim=-3))
            )
            coefficients = torch.tanh(self.net(features))
            difference = torch.cat(
                (
                    (centre - patch[..., 0, :, :]).flatten(start_dim=-2),
                    (patch[..., 2, :, :] - centre).flatten(start_dim=-2),
                ),
                dim=-1,
            )
            gate = torch.tanh(difference.abs().mean(dim=-1))
            return coefficients, centre_only, gate

        normal_logits = self.normal_net(centre.flatten(start_dim=-2))
        difference = torch.cat(
            (
                (centre - patch[..., 0, :, :]).flatten(start_dim=-2),
                (patch[..., 2, :, :] - centre).flatten(start_dim=-2),
            ),
            dim=-1,
        )
        gate = torch.tanh(difference.abs().mean(dim=-1))
        transverse_features = torch.cat(
            (centre.flatten(start_dim=-2), difference), dim=-1
        )
        transverse_logits = self.transverse_net(transverse_features)
        centre_only = torch.tanh(normal_logits)
        coefficients = torch.tanh(normal_logits + gate[..., None] * transverse_logits)
        return coefficients, centre_only, gate

    def coefficients(self, state: torch.Tensor) -> torch.Tensor:
        return self.coefficient_details(state)[0]

    def forward_oriented(self, state: torch.Tensor) -> torch.Tensor:
        base = C.hllc_faces_oriented(state)
        left = state[..., :-1, :]
        right = state[..., 1:, :]
        matrix, strengths, speeds = C.roe_waves_oriented(left, right)
        correction = torch.einsum(
            "...ij,...j->...i",
            matrix,
            self.coefficients(state) * speeds * strengths,
        )
        interior = base[..., 1:-1, :] - 0.5 * correction
        return torch.cat((base[..., :1, :], interior, base[..., -1:, :]), dim=-2)


class DirectionalPatchCentralNonnegativeRoe(nn.Module):
    """Central physical flux with nonnegative learned Roe dissipation.

    A zero network output gives multiplier one for every Roe wave, hence the
    standard entropy-fixed Roe flux rather than HLLC.  The construction is
    the four-wave 2-D analogue of ``CentralRoeUpwindFlux`` used by the
    retained 1-D Euler experiments.
    """

    def __init__(
        self,
        mean: np.ndarray,
        std: np.ndarray,
        width: int = 72,
    ) -> None:
        super().__init__()
        self.net = _mlp(3 * 6 * 4, width)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def normalized_patch(self, state: torch.Tensor) -> torch.Tensor:
        patch = oriented_primitive_patch(state)
        return (patch - self.mean) / self.std

    def multiplier_details(
        self, state: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Return full-patch, centre-only, and transverse-gate diagnostics."""

        patch = self.normalized_patch(state)
        centre = patch[..., 1, :, :]
        features = patch.flatten(start_dim=-3)
        repeated = torch.stack((centre, centre, centre), dim=-3)
        multipliers = 1.0 + torch.tanh(self.net(features))
        centre_only = 1.0 + torch.tanh(
            self.net(repeated.flatten(start_dim=-3))
        )
        difference = torch.cat(
            (
                (centre - patch[..., 0, :, :]).flatten(start_dim=-2),
                (patch[..., 2, :, :] - centre).flatten(start_dim=-2),
            ),
            dim=-1,
        )
        gate = torch.tanh(difference.abs().mean(dim=-1))
        return multipliers, centre_only, gate

    def multipliers(self, state: torch.Tensor) -> torch.Tensor:
        return self.multiplier_details(state)[0]

    @staticmethod
    def _assemble_oriented(
        state: torch.Tensor, multipliers: torch.Tensor
    ) -> torch.Tensor:
        left = state[..., :-1, :]
        right = state[..., 1:, :]
        matrix, strengths, speeds = C.roe_waves_oriented(left, right)
        dissipation = torch.einsum(
            "...ij,...j->...i",
            matrix,
            multipliers * speeds * strengths,
        )
        physical = C.physical_flux_oriented(state)
        central = 0.5 * (physical[..., :-1, :] + physical[..., 1:, :])
        interior = central - 0.5 * dissipation
        return torch.cat(
            (physical[..., :1, :], interior, physical[..., -1:, :]), dim=-2
        )

    def standard_roe_faces_oriented(self, state: torch.Tensor) -> torch.Tensor:
        shape = (*state.shape[:-2], state.shape[-2] - 1, 4)
        multipliers = torch.ones(shape, dtype=state.dtype, device=state.device)
        return self._assemble_oriented(state, multipliers)

    def forward_oriented(self, state: torch.Tensor) -> torch.Tensor:
        return self._assemble_oriented(state, self.multipliers(state))


class DirectionalNormalCentralNonnegativeRoe(nn.Module):
    """The retained central/nonnegative Roe map on a normal six-cell stencil.

    This is the dimension-consistent 2-D analogue of the six-cell 1-D model:
    only cells along the face normal enter the shared neural network.  The
    transverse numerical transport, when requested, is supplied separately by
    a fixed solver and introduces no learned parameters.
    """

    def __init__(
        self,
        mean: np.ndarray,
        std: np.ndarray,
        width: int = 72,
    ) -> None:
        super().__init__()
        self.shifts = NORMAL_SHIFTS
        self.net = _mlp(6 * 4, width)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def normalized_stencil(self, state: torch.Tensor) -> torch.Tensor:
        values = C.primitive(state)
        cells = values.shape[-2]
        interfaces = torch.arange(1, cells, device=state.device)
        offsets = torch.tensor(self.shifts, device=state.device)
        indices = (interfaces[:, None] + offsets[None, :]).clamp(0, cells - 1)
        gathered = values.index_select(-2, indices.reshape(-1)).reshape(
            *values.shape[:-2], cells - 1, len(self.shifts), 4
        )
        return (gathered - self.mean) / self.std

    def multipliers(self, state: torch.Tensor) -> torch.Tensor:
        features = self.normalized_stencil(state).flatten(start_dim=-2)
        return 1.0 + torch.tanh(self.net(features))

    @staticmethod
    def _assemble_oriented(
        state: torch.Tensor, multipliers: torch.Tensor
    ) -> torch.Tensor:
        return DirectionalPatchCentralNonnegativeRoe._assemble_oriented(
            state, multipliers
        )

    def standard_roe_faces_oriented(self, state: torch.Tensor) -> torch.Tensor:
        shape = (*state.shape[:-2], state.shape[-2] - 1, 4)
        multipliers = torch.ones(shape, dtype=state.dtype, device=state.device)
        return self._assemble_oriented(state, multipliers)

    def forward_oriented(self, state: torch.Tensor) -> torch.Tensor:
        return self._assemble_oriented(state, self.multipliers(state))


class HCFL2DTransverse(C.HCFL2DEuler):
    """HCFL wrapper retaining the original hard projections and safety path."""

    def __init__(
        self,
        mean: np.ndarray,
        std: np.ndarray,
        variant: str,
    ) -> None:
        nn.Module.__init__(self)
        self.flux_net = DirectionalPatchHLLCRoeCorrection(mean, std, variant)
        self.variant = variant
        self.stencil_cells = 18


class HCFL2DCentralNonnegative(C.HCFL2DEuler):
    """Retained 1-D HCFL flux design extended to oriented 2-D faces."""

    def __init__(self, mean: np.ndarray, std: np.ndarray) -> None:
        nn.Module.__init__(self)
        self.flux_net = DirectionalPatchCentralNonnegativeRoe(mean, std)
        self.variant = CONSISTENT_VARIANT
        self.stencil_cells = 18


class HCFL2DNormalCentralNonnegative(C.HCFL2DEuler):
    """Six normal cells, with the retained 1-D HCFL flux parameterization."""

    def __init__(self, mean: np.ndarray, std: np.ndarray) -> None:
        nn.Module.__init__(self)
        self.flux_net = DirectionalNormalCentralNonnegativeRoe(mean, std)
        self.variant = NORMAL6_CONSISTENT_VARIANT
        self.stencil_cells = 6


def make_model(
    mean: np.ndarray,
    std: np.ndarray,
    variant: str,
) -> nn.Module:
    if variant == "normal6wide":
        # w^2 + 30w + 4 parameters for a 24-w-w-4 MLP; w=90 gives
        # exactly 10,804 parameters, matching flat18 exactly.
        return C.HCFL2DEuler(mean, std, width=90, stencil_cells=6)
    if variant == CONSISTENT_VARIANT:
        return HCFL2DCentralNonnegative(mean, std)
    if variant == NORMAL6_CONSISTENT_VARIANT:
        return HCFL2DNormalCentralNonnegative(mean, std)
    return HCFL2DTransverse(mean, std, variant)
