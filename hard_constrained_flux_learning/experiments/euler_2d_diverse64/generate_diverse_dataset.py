"""Generate diverse high-resolution 2-D Euler trajectories restricted to 64^2.

This script is executed inside the pinned Clawpack 5.10 container used by the
quadrant experiment.  Fine-grid finite-volume conserved variables are block
averaged to the deployment grid.  Primitive variables are never averaged.

The six initial-condition families deliberately exercise different local wave
structures.  They are balanced by construction, and every case is a
deterministic function of ``(seed, case_index)`` so data-size studies can use
nested prefixes without silently changing earlier cases.
"""

from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "euler_2d_clawpack_quadrants"
sys.path.insert(0, str(BASE))
import generate_pyclaw_dataset as Q  # noqa: E402


FAMILIES = (
    "oblique_riemann",
    "contact_shear",
    "oblique_quadrant",
    "radial_interface",
    "colliding_waves",
    "smooth_packet",
)


def _log_uniform(rng: np.random.Generator, lower: float, upper: float) -> float:
    return float(np.exp(rng.uniform(math.log(lower), math.log(upper))))


def _random_state(
    rng: np.random.Generator,
    *,
    velocity_limit: float = 0.7,
    density_bounds: tuple[float, float] = (0.4, 2.4),
    pressure_bounds: tuple[float, float] = (0.3, 2.4),
) -> list[float]:
    return [
        _log_uniform(rng, *density_bounds),
        float(rng.uniform(-velocity_limit, velocity_limit)),
        float(rng.uniform(-velocity_limit, velocity_limit)),
        _log_uniform(rng, *pressure_bounds),
    ]


def _normal_velocity_state(
    rho: float,
    normal_velocity: float,
    tangential_velocity: float,
    pressure: float,
    theta: float,
) -> list[float]:
    normal = np.asarray((math.cos(theta), math.sin(theta)))
    tangent = np.asarray((-math.sin(theta), math.cos(theta)))
    velocity = normal_velocity * normal + tangential_velocity * tangent
    return [rho, float(velocity[0]), float(velocity[1]), pressure]


def sample_spec(family: str, seed: int, case_index: int) -> dict[str, Any]:
    """Sample one deterministic, positive primitive-state specification."""

    rng = np.random.default_rng(np.random.SeedSequence((seed, case_index)))
    theta = float(rng.uniform(0.0, math.pi))
    spec: dict[str, Any] = {
        "family": family,
        "theta": theta,
        "case_index": case_index,
    }
    if family == "oblique_riemann":
        spec.update(
            offset=float(rng.uniform(-0.14, 0.14)),
            left=_random_state(rng),
            right=_random_state(rng),
        )
    elif family == "contact_shear":
        pressure = _log_uniform(rng, 0.35, 2.3)
        normal_velocity = float(rng.uniform(-0.55, 0.55))
        tangential_centre = float(rng.uniform(-0.25, 0.25))
        tangential_jump = float(rng.uniform(0.15, 0.75))
        spec.update(
            offset=float(rng.uniform(-0.14, 0.14)),
            left=_normal_velocity_state(
                _log_uniform(rng, 0.35, 2.5),
                normal_velocity,
                tangential_centre - 0.5 * tangential_jump,
                pressure,
                theta,
            ),
            right=_normal_velocity_state(
                _log_uniform(rng, 0.35, 2.5),
                normal_velocity,
                tangential_centre + 0.5 * tangential_jump,
                pressure,
                theta,
            ),
        )
    elif family == "oblique_quadrant":
        spec.update(
            centre=[float(rng.uniform(0.38, 0.62)), float(rng.uniform(0.38, 0.62))],
            offsets=[float(rng.uniform(-0.06, 0.06)), float(rng.uniform(-0.06, 0.06))],
            states=[_random_state(rng, velocity_limit=0.65) for _ in range(4)],
        )
    elif family == "radial_interface":
        spec.update(
            centre=[float(rng.uniform(0.36, 0.64)), float(rng.uniform(0.36, 0.64))],
            radius=float(rng.uniform(0.11, 0.23)),
            inside=_random_state(rng, velocity_limit=0.45),
            outside=_random_state(rng, velocity_limit=0.45),
        )
    elif family == "colliding_waves":
        pressure = _log_uniform(rng, 0.35, 1.8)
        density = _log_uniform(rng, 0.5, 1.8)
        speed = float(rng.uniform(0.3, 0.85))
        tangential = float(rng.uniform(-0.25, 0.25))
        middle_pressure = pressure * float(rng.uniform(0.65, 1.25))
        middle_density = density * float(rng.uniform(0.7, 1.35))
        spec.update(
            offset=float(rng.uniform(-0.08, 0.08)),
            half_width=float(rng.uniform(0.09, 0.18)),
            left=_normal_velocity_state(
                density * float(rng.uniform(0.8, 1.2)),
                speed,
                tangential,
                pressure * float(rng.uniform(0.85, 1.15)),
                theta,
            ),
            middle=_normal_velocity_state(
                middle_density,
                float(rng.uniform(-0.1, 0.1)),
                tangential,
                middle_pressure,
                theta,
            ),
            right=_normal_velocity_state(
                density * float(rng.uniform(0.8, 1.2)),
                -speed,
                tangential,
                pressure * float(rng.uniform(0.85, 1.15)),
                theta,
            ),
        )
    elif family == "smooth_packet":
        spec.update(
            centre=[float(rng.uniform(0.38, 0.62)), float(rng.uniform(0.38, 0.62))],
            width=float(rng.uniform(0.14, 0.24)),
            base=_random_state(
                rng,
                velocity_limit=0.35,
                density_bounds=(0.7, 1.8),
                pressure_bounds=(0.7, 1.8),
            ),
            density_amplitude=float(rng.uniform(0.08, 0.28)),
            pressure_amplitude=float(rng.uniform(0.06, 0.24)),
            velocity_amplitude=float(rng.uniform(0.06, 0.24)),
            wave_number=int(rng.integers(1, 4)),
            phase=float(rng.uniform(0.0, 2.0 * math.pi)),
        )
    else:
        raise ValueError(f"Unknown family: {family}")
    return spec


def sample_specs(count: int, seed: int) -> list[dict[str, Any]]:
    return [
        sample_spec(FAMILIES[index % len(FAMILIES)], seed, index)
        for index in range(count)
    ]


def _coordinates(cells: int) -> tuple[np.ndarray, np.ndarray]:
    centres = (np.arange(cells, dtype=np.float64) + 0.5) / cells
    return np.meshgrid(centres, centres, indexing="xy")


def initial_state(cells: int, spec: dict[str, Any]) -> np.ndarray:
    """Evaluate a specification as fine-grid conserved cell-centre data."""

    xx, yy = _coordinates(cells)
    theta = float(spec["theta"])
    normal = np.asarray((math.cos(theta), math.sin(theta)))
    tangent = np.asarray((-math.sin(theta), math.cos(theta)))
    family = str(spec["family"])

    if family in ("oblique_riemann", "contact_shear"):
        signed = (
            (xx - 0.5) * normal[0]
            + (yy - 0.5) * normal[1]
            - float(spec["offset"])
        )
        primitive = np.empty((cells, cells, 4), dtype=np.float64)
        primitive[signed < 0.0] = np.asarray(spec["left"])
        primitive[signed >= 0.0] = np.asarray(spec["right"])
    elif family == "oblique_quadrant":
        centre = np.asarray(spec["centre"])
        first = (
            (xx - centre[0]) * normal[0]
            + (yy - centre[1]) * normal[1]
            - float(spec["offsets"][0])
        )
        second = (
            (xx - centre[0]) * tangent[0]
            + (yy - centre[1]) * tangent[1]
            - float(spec["offsets"][1])
        )
        states = np.asarray(spec["states"])
        primitive = np.empty((cells, cells, 4), dtype=np.float64)
        primitive[(first >= 0.0) & (second >= 0.0)] = states[0]
        primitive[(first < 0.0) & (second >= 0.0)] = states[1]
        primitive[(first < 0.0) & (second < 0.0)] = states[2]
        primitive[(first >= 0.0) & (second < 0.0)] = states[3]
    elif family == "radial_interface":
        centre = np.asarray(spec["centre"])
        radius = np.sqrt((xx - centre[0]) ** 2 + (yy - centre[1]) ** 2)
        primitive = np.empty((cells, cells, 4), dtype=np.float64)
        primitive[radius < float(spec["radius"])] = np.asarray(spec["inside"])
        primitive[radius >= float(spec["radius"])] = np.asarray(spec["outside"])
    elif family == "colliding_waves":
        signed = (
            (xx - 0.5) * normal[0]
            + (yy - 0.5) * normal[1]
            - float(spec["offset"])
        )
        half_width = float(spec["half_width"])
        primitive = np.empty((cells, cells, 4), dtype=np.float64)
        primitive[signed < -half_width] = np.asarray(spec["left"])
        primitive[np.abs(signed) <= half_width] = np.asarray(spec["middle"])
        primitive[signed > half_width] = np.asarray(spec["right"])
    elif family == "smooth_packet":
        centre = np.asarray(spec["centre"])
        dx = xx - centre[0]
        dy = yy - centre[1]
        normal_coordinate = dx * normal[0] + dy * normal[1]
        transverse_coordinate = dx * tangent[0] + dy * tangent[1]
        width = float(spec["width"])
        envelope = np.exp(-(dx * dx + dy * dy) / (2.0 * width * width))
        phase = (
            2.0
            * math.pi
            * float(spec["wave_number"])
            * (normal_coordinate + 0.35 * transverse_coordinate)
            + float(spec["phase"])
        )
        base = np.asarray(spec["base"], dtype=np.float64)
        rho = base[0] * np.exp(
            float(spec["density_amplitude"]) * envelope * np.sin(phase)
        )
        pressure = base[3] * np.exp(
            float(spec["pressure_amplitude"]) * envelope * np.cos(phase)
        )
        sound_speed = math.sqrt(Q.GAMMA * base[3] / base[0])
        velocity_perturbation = (
            float(spec["velocity_amplitude"])
            * sound_speed
            * envelope
            * np.sin(phase + 0.5 * math.pi)
        )
        u = base[1] + velocity_perturbation * normal[0]
        v = base[2] + velocity_perturbation * normal[1]
        primitive = np.stack((rho, u, v, pressure), axis=-1)
    else:
        raise ValueError(f"Unknown family: {family}")

    if not np.isfinite(primitive).all():
        raise RuntimeError(f"Non-finite primitive initial condition for {family}")
    if float(primitive[..., 0].min()) <= 0.0 or float(primitive[..., 3].min()) <= 0.0:
        raise RuntimeError(f"Non-positive primitive initial condition for {family}")
    return Q.primitive_to_conserved(primitive)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate(args: argparse.Namespace) -> dict[str, Any]:
    specs = sample_specs(args.count, args.seed)
    trajectories: list[np.ndarray] = []
    native_coarse: list[np.ndarray] = []
    expected_times: np.ndarray | None = None

    for index, spec in enumerate(specs):
        fine_initial = initial_state(args.fine_cells, spec)
        trajectory, times = Q.run_pyclaw(
            fine_initial,
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
            coarse_trajectory, coarse_times = Q.run_pyclaw(
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
                    "count": len(specs),
                    "family": spec["family"],
                    "fine_cells": args.fine_cells,
                }
            ),
            flush=True,
        )

    metadata = {
        "clawpack_version": version("clawpack"),
        "solver": "classic_pyclaw_euler_4wave_2D",
        "transverse_waves": 2,
        "boundary": "constant_extrapolation",
        "gamma": Q.GAMMA,
        "families": list(FAMILIES),
        "balanced_family_cycle": True,
        "count": len(specs),
        "seed": args.seed,
        "fine_cells": args.fine_cells,
        "coarse_cells": args.coarse_cells,
        "tfinal": args.tfinal,
        "output_intervals": args.output_intervals,
        "restriction": "conservative_block_average_of_conserved_variables",
        "common_coarse_initial_state_for_native_baseline": bool(
            args.include_native_coarse
        ),
        "generator_sha256": _sha256(Path(__file__)),
        "base_generator_sha256": _sha256(BASE / "generate_pyclaw_dataset.py"),
    }
    result: dict[str, Any] = {
        "q": np.stack(trajectories),
        "times": expected_times,
        "family_names": np.asarray([spec["family"] for spec in specs]),
        "specs_json": np.asarray(json.dumps(specs, sort_keys=True)),
        "metadata_json": np.asarray(json.dumps(metadata, sort_keys=True)),
    }
    if native_coarse:
        result["native_roe_coarse"] = np.stack(native_coarse)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--fine-cells", type=int, required=True)
    parser.add_argument("--coarse-cells", type=int, default=64)
    parser.add_argument("--tfinal", type=float, default=0.1)
    parser.add_argument("--output-intervals", type=int, default=20)
    parser.add_argument("--include-native-coarse", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.count <= 0 or args.output_intervals <= 0:
        parser.error("Counts must be positive")
    if args.fine_cells % args.coarse_cells:
        parser.error("--fine-cells must be divisible by --coarse-cells")

    result = generate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **result)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "shape": list(result["q"].shape),
                "families": {
                    family: int(np.sum(result["family_names"] == family))
                    for family in FAMILIES
                },
                "keys": sorted(result),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
