"""Transverse-aware conservative face models for the 2-D Euler ablation.

Both variants use an oriented 3-by-6 primitive-variable patch for every face.
The ``flat18`` model is deliberately unstructured: it must learn for itself
which transverse information matters.  The ``gated18`` model separates the
existing normal six-cell mapping from a transverse residual and reduces
exactly to its normal branch when the three rows agree.
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


def make_model(
    mean: np.ndarray,
    std: np.ndarray,
    variant: str,
) -> nn.Module:
    if variant == "normal6wide":
        # w^2 + 30w + 4 parameters for a 24-w-w-4 MLP; w=90 gives
        # exactly 10,804 parameters, matching flat18 exactly.
        return C.HCFL2DEuler(mean, std, width=90, stencil_cells=6)
    return HCFL2DTransverse(mean, std, variant)
