"""Generate the locked 2-D Euler quadrant datasets with PyClaw 5.10.

This file is executed inside the container built from ``Dockerfile.clawpack``.
It deliberately has no dependency on the PyTorch HCFL implementation.
"""

from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from clawpack import pyclaw, riemann


GAMMA = 1.4
SOURCE_URL = "https://www.clawpack.org/gallery/pyclaw/gallery/quadrants.html"

OFFICIAL_STATES = np.asarray(
    [
        [1.5, 0.0, 0.0, 1.5],
        [0.532258064516129, 1.206045378311055, 0.0, 0.3],
        [
            0.137992831541219,
            1.206045378311055,
            1.206045378311055,
            0.029032258064516,
        ],
        [0.532258064516129, 0.0, 1.206045378311055, 0.3],
    ],
    dtype=np.float64,
)


def primitive_to_conserved(values: np.ndarray) -> np.ndarray:
    rho, u, v, pressure = np.moveaxis(values, -1, 0)
    energy = pressure / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)
    return np.stack((rho, rho * u, rho * v, energy), axis=-1)


def pressure(state: np.ndarray) -> np.ndarray:
    rho, mx, my, energy = np.moveaxis(state, -1, 0)
    return (GAMMA - 1.0) * (
        energy - 0.5 * (mx * mx + my * my) / rho
    )


def quadrant_state(
    cells: int,
    primitive_states: np.ndarray,
    x_split: float,
    y_split: float,
) -> np.ndarray:
    centers = (np.arange(cells, dtype=np.float64) + 0.5) / cells
    xx, yy = np.meshgrid(centers, centers, indexing="xy")
    left = xx < x_split
    lower = yy < y_split
    primitive = np.empty((cells, cells, 4), dtype=np.float64)
    # Parameter order: upper-right, upper-left, lower-left, lower-right.
    primitive[~left & ~lower] = primitive_states[0]
    primitive[left & ~lower] = primitive_states[1]
    primitive[left & lower] = primitive_states[2]
    primitive[~left & lower] = primitive_states[3]
    return primitive_to_conserved(primitive)


def conservative_restrict(state: np.ndarray, coarse_cells: int) -> np.ndarray:
    fine_y, fine_x, channels = state.shape
    if fine_x % coarse_cells or fine_y % coarse_cells:
        raise ValueError("Fine grid must be divisible by the coarse grid")
    ratio_y = fine_y // coarse_cells
    ratio_x = fine_x // coarse_cells
    return state.reshape(
        coarse_cells, ratio_y, coarse_cells, ratio_x, channels
    ).mean(axis=(1, 3))


def build_solver() -> pyclaw.ClawSolver2D:
    solver = pyclaw.ClawSolver2D(riemann.euler_4wave_2D)
    solver.transverse_waves = 2
    solver.all_bcs = pyclaw.BC.extrap
    return solver


def run_pyclaw(
    initial_state: np.ndarray,
    coarse_cells: int,
    final_time: float,
    output_intervals: int,
) -> tuple[np.ndarray, np.ndarray]:
    fine_y, fine_x, channels = initial_state.shape
    if fine_x != fine_y or channels != 4:
        raise ValueError("Expected a square [y,x,4] Euler state")

    x_dimension = pyclaw.Dimension(0.0, 1.0, fine_x, name="x")
    y_dimension = pyclaw.Dimension(0.0, 1.0, fine_y, name="y")
    domain = pyclaw.Domain([x_dimension, y_dimension])
    state = pyclaw.State(domain, 4)
    state.problem_data["gamma"] = GAMMA
    state.problem_data["gamma1"] = GAMMA - 1.0
    state.q[...] = np.moveaxis(np.swapaxes(initial_state, 0, 1), -1, 0)

    controller = pyclaw.Controller()
    controller.solution = pyclaw.Solution(state, domain)
    controller.solver = build_solver()
    controller.tfinal = final_time
    controller.num_output_times = output_intervals
    controller.keep_copy = True
    controller.output_format = None
    controller.verbosity = 0
    controller.run()

    frames: list[np.ndarray] = []
    times: list[float] = []
    for frame in controller.frames:
        # PyClaw stores [equation,x,y]; the repository convention is [y,x,equation].
        current = np.swapaxes(np.moveaxis(np.asarray(frame.q), 0, -1), 0, 1)
        restricted = conservative_restrict(current, coarse_cells).astype(
            np.float32, copy=False
        )
        if not np.isfinite(restricted).all():
            raise RuntimeError(f"Non-finite state at t={frame.t}")
        if float(restricted[..., 0].min()) <= 0.0:
            raise RuntimeError(f"Non-positive density at t={frame.t}")
        if float(pressure(restricted).min()) <= 0.0:
            raise RuntimeError(f"Non-positive pressure at t={frame.t}")
        frames.append(restricted.copy())
        times.append(float(frame.t))
    return np.stack(frames), np.asarray(times, dtype=np.float64)


def sample_random_quadrants(count: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    states = np.empty((count, 4, 4), dtype=np.float64)
    states[..., 0] = np.exp(
        rng.uniform(math.log(0.5), math.log(2.0), size=(count, 4))
    )
    states[..., 1] = rng.uniform(-0.6, 0.6, size=(count, 4))
    states[..., 2] = rng.uniform(-0.6, 0.6, size=(count, 4))
    states[..., 3] = np.exp(
        rng.uniform(math.log(0.4), math.log(2.0), size=(count, 4))
    )
    splits = rng.uniform(0.3, 0.7, size=(count, 2))
    return states, splits


def generate(args: argparse.Namespace) -> dict[str, Any]:
    if args.kind == "official":
        primitive_states = OFFICIAL_STATES[None]
        splits = np.asarray([[0.8, 0.8]], dtype=np.float64)
    else:
        primitive_states, splits = sample_random_quadrants(args.count, args.seed)

    trajectories: list[np.ndarray] = []
    native_coarse: list[np.ndarray] = []
    expected_times: np.ndarray | None = None
    for index, (states, split) in enumerate(zip(primitive_states, splits)):
        initial_fine = quadrant_state(
            args.fine_cells, states, float(split[0]), float(split[1])
        )
        trajectory, times = run_pyclaw(
            initial_fine,
            args.coarse_cells,
            args.tfinal,
            args.output_intervals,
        )
        if expected_times is None:
            expected_times = times
        elif not np.array_equal(times, expected_times):
            raise RuntimeError("PyClaw returned inconsistent saved times")
        trajectories.append(trajectory)

        if args.include_native_coarse:
            coarse_trajectory, coarse_times = run_pyclaw(
                trajectory[0].astype(np.float64),
                args.coarse_cells,
                args.tfinal,
                args.output_intervals,
            )
            if not np.array_equal(coarse_times, times):
                raise RuntimeError("Native coarse output times differ from reference")
            native_coarse.append(coarse_trajectory)
        print(
            json.dumps(
                {
                    "completed": index + 1,
                    "count": len(primitive_states),
                    "kind": args.kind,
                    "fine_cells": args.fine_cells,
                }
            ),
            flush=True,
        )

    metadata = {
        "source": SOURCE_URL,
        "clawpack_version": version("clawpack"),
        "solver": "classic_pyclaw_euler_4wave_2D",
        "transverse_waves": 2,
        "boundary": "constant_extrapolation",
        "gamma": GAMMA,
        "kind": args.kind,
        "count": len(primitive_states),
        "seed": args.seed if args.kind == "random" else None,
        "fine_cells": args.fine_cells,
        "coarse_cells": args.coarse_cells,
        "tfinal": args.tfinal,
        "output_intervals": args.output_intervals,
        "restriction": "conservative_block_average",
        "common_coarse_initial_state_for_native_baseline": bool(
            args.include_native_coarse
        ),
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    result: dict[str, Any] = {
        "q": np.stack(trajectories),
        "times": expected_times,
        "primitive_states": primitive_states,
        "split_locations": splits,
        "metadata_json": np.asarray(json.dumps(metadata, sort_keys=True)),
    }
    if native_coarse:
        result["native_roe_coarse"] = np.stack(native_coarse)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("random", "official"), required=True)
    parser.add_argument("--count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--fine-cells", type=int, required=True)
    parser.add_argument("--coarse-cells", type=int, default=64)
    parser.add_argument("--tfinal", type=float, required=True)
    parser.add_argument("--output-intervals", type=int, required=True)
    parser.add_argument("--include-native-coarse", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.kind == "official":
        args.count = 1
    if args.fine_cells % args.coarse_cells:
        parser.error("--fine-cells must be divisible by --coarse-cells")
    if args.count <= 0 or args.output_intervals <= 0:
        parser.error("Counts must be positive")

    result = generate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **result)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "shape": list(result["q"].shape),
                "keys": sorted(result),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
