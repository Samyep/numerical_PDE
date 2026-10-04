"""Independent grid-refinement check for the local radial-dam reference."""

from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
import torch

import swe_radial_common as C


@torch.no_grad()
def solve(cells: int, device: torch.device) -> tuple[torch.Tensor, float]:
    state = C.radial_initial_state(
        np.asarray([0.5], dtype=np.float32),
        np.asarray([2.0], dtype=np.float32),
        cells,
        device,
    )
    snapshots = [state]
    started = time.perf_counter()
    for _ in range(C.ORDINARY_STATES - 1):
        remaining = C.SNAPSHOT_DT
        while remaining > 1.0e-10:
            dt = min(remaining, 0.42 / max(C.maximum_2d_rate(state), 1.0e-12))
            for _ in range(18):
                candidate = C.reference_ssprk2(state, dt)
                if bool(torch.isfinite(candidate).all()) and float(candidate[..., 0].min()) > C.H_FLOOR:
                    break
                dt *= 0.5
            else:
                raise RuntimeError("Reference grid audit lost positive depth")
            state = candidate
            remaining -= dt
        snapshots.append(state)
    return torch.stack(snapshots, dim=1), time.perf_counter() - started


def restrict(sequence: torch.Tensor) -> torch.Tensor:
    batch, times, cells, _, channels = sequence.shape
    factor = cells // C.N_COARSE
    return sequence.reshape(
        batch,
        times,
        C.N_COARSE,
        factor,
        C.N_COARSE,
        factor,
        channels,
    ).mean(dim=(3, 5))


def main() -> None:
    output = Path(__file__).resolve().parent / "results" / "swe_radial"
    output.mkdir(parents=True, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    sequences: dict[int, torch.Tensor] = {}
    runtimes: dict[int, float] = {}
    for cells in (64, 128, 256):
        sequence, runtime = solve(cells, device)
        sequences[cells] = restrict(sequence).cpu()
        runtimes[cells] = runtime
    scale = C.primitive_channel_scale(sequences[256])
    report = {
        "case": "radius=0.5, inner_height=2.0",
        "coarse_comparison_cells": C.N_COARSE,
        "saved_states": C.ORDINARY_STATES,
        "final_time": (C.ORDINARY_STATES - 1) * C.SNAPSHOT_DT,
        "runtime_seconds": runtimes,
        "primitive_nrmse": {
            "64_vs_128_rollout": C.primitive_nrmse(
                sequences[64][:, 1:], sequences[128][:, 1:], scale
            ),
            "128_vs_256_rollout": C.primitive_nrmse(
                sequences[128][:, 1:], sequences[256][:, 1:], scale
            ),
            "128_vs_256_final": C.primitive_nrmse(
                sequences[128][:, -1:], sequences[256][:, -1:], scale
            ),
        },
        "interpretation": (
            "The chosen 128-grid reference matches PDEBench's native spatial "
            "resolution, but is not continuum-exact; the 128-vs-256 gap is "
            "reported as reference-discretization uncertainty."
        ),
    }
    (output / "reference_grid_audit.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
