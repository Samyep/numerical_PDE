"""Freeze parameters, fresh query banks, and the independent 2D reference."""

from __future__ import annotations

from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.reference2d import (
    _spectral_fields,
    solve_periodic_reference_2d,
)

from .benchmark import frozen_basis, make_benchmark, reduced_to_ambient_z
from .protocol import (
    BASE_SEED,
    DIMENSIONS,
    FROZEN_PROTOCOL_JSON,
    PARAMETER_DIR,
    PARAMETER_MANIFEST,
    PDE,
    QUERY_COUNT,
    REFERENCE_VALIDATION,
    ensure_directories,
    file_sha256,
    frozen_protocol,
    query_path,
    write_json_atomic,
    write_npz_atomic,
)


def _npy_bytes(array: np.ndarray) -> bytes:
    import io

    buffer = io.BytesIO()
    np.save(buffer, np.asarray(array), allow_pickle=False)
    return buffer.getvalue()


def _freeze_npy(path: Path, array: np.ndarray) -> None:
    """Create once, or verify that an already frozen artifact is identical."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _npy_bytes(array)
    if path.exists():
        if path.read_bytes() != payload:
            raise RuntimeError(f"refusing to overwrite changed frozen artifact: {path}")
        return
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(payload)
    temporary.replace(path)


def _query_inputs(dimension: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(
        np.random.SeedSequence([BASE_SEED, int(dimension), QUERY_COUNT])
    )
    t = np.zeros(QUERY_COUNT, dtype=np.float64)
    x = np.zeros((QUERY_COUNT, int(dimension)), dtype=np.float64)
    x[1:] = rng.standard_normal((QUERY_COUNT - 1, int(dimension)))
    return t, x


def _residual_diagnostics(reference, driver) -> dict[str, float]:
    """Centered PDE residual on stored slices, evaluated one slice at a time."""
    residual_squares = 0.0
    residual_count = 0
    residual_max = 0.0
    # Use every fourth interior time slice.  This is independent of the
    # splitting update and avoids materializing several 3D derivative arrays.
    for index in range(1, len(reference.times) - 1, 4):
        center = reference.values[index]
        time_derivative = (
            reference.values[index + 1] - reference.values[index - 1]
        ) / (2.0 * reference.store_dt)
        d1, d2, laplace = _spectral_fields(center)
        z = math.sqrt(2.0) * np.stack((d1, d2), axis=-1)
        residual = time_derivative + laplace + driver.driver_reduced_numpy(center, z)
        residual_squares += float(np.square(residual).sum())
        residual_count += residual.size
        residual_max = max(residual_max, float(np.max(np.abs(residual))))
    return {
        "pde_residual_rmse": math.sqrt(residual_squares / residual_count),
        "pde_residual_max_abs": residual_max,
        "pde_residual_points": residual_count,
        "pde_residual_time_stencil": "centered over stored primary slices",
    }


def _reference_and_validation() -> tuple[Any, Any, dict[str, Any]]:
    driver = make_benchmark(2)
    audit = solve_periodic_reference_2d(
        driver, n=128, dt=2.5e-4, store_dt=0.002
    )
    primary = solve_periodic_reference_2d(
        driver, n=256, dt=1.25e-4, store_dt=0.002
    )

    query_t: list[np.ndarray] = []
    query_s: list[np.ndarray] = []
    for dimension in DIMENSIONS:
        t, x = _query_inputs(dimension)
        query_t.append(t)
        query_s.append(x @ frozen_basis(dimension))
    pooled_t = np.concatenate(query_t)
    pooled_s = np.concatenate(query_s)
    primary_value = primary.value(pooled_t, pooled_s)
    audit_value = audit.value(pooled_t, pooled_s)
    primary_gradient = primary.gradient(pooled_t, pooled_s)
    audit_gradient = audit.gradient(pooled_t, pooled_s)

    terminal_rng = np.random.default_rng(
        np.random.SeedSequence([BASE_SEED, 2, 99])
    )
    terminal_s = terminal_rng.uniform(0.0, 2.0 * np.pi, (8192, 2))
    terminal_t = np.full(len(terminal_s), driver.T, dtype=np.float64)
    terminal_delta = (
        primary.value(terminal_t, terminal_s)
        - driver.terminal_from_coordinates_numpy(terminal_s)
    )
    value_delta = primary_value - audit_value
    gradient_delta = primary_gradient - audit_gradient
    residual = _residual_diagnostics(primary, driver)
    validation: dict[str, Any] = {
        "pde": PDE,
        "solver": "independent periodic Fourier Strang splitting with explicit midpoint nonlinear step",
        "primary": {
            "n": primary.n,
            "dt": primary.dt,
            "store_dt": primary.store_dt,
        },
        "audit": {"n": audit.n, "dt": audit.dt, "store_dt": audit.store_dt},
        "comparison_points": len(pooled_t),
        "value_refinement_rmse": float(np.sqrt(np.mean(np.square(value_delta)))),
        "value_refinement_max_abs": float(np.max(np.abs(value_delta))),
        "gradient_refinement_vector_rmse": float(
            np.sqrt(np.mean(np.square(gradient_delta).sum(axis=1)))
        ),
        "gradient_refinement_max_norm": float(
            np.max(np.linalg.norm(gradient_delta, axis=1))
        ),
        "terminal_interpolation_rmse": float(
            np.sqrt(np.mean(np.square(terminal_delta)))
        ),
        "terminal_interpolation_max_abs": float(np.max(np.abs(terminal_delta))),
        "working_value_min": float(primary.values.min()),
        "working_value_max": float(primary.values.max()),
        **residual,
    }
    gates = {
        "value_refinement_rmse": validation["value_refinement_rmse"] <= 5e-5,
        "value_refinement_max_abs": validation["value_refinement_max_abs"] <= 2e-4,
        "pde_residual_rmse": validation["pde_residual_rmse"] <= 2e-6,
        "terminal_interpolation_max_abs": validation[
            "terminal_interpolation_max_abs"
        ] <= 1e-4,
        "invariant_domain": (
            validation["working_value_min"] >= -1e-8
            and validation["working_value_max"] <= 1.0 + 1e-8
        ),
    }
    validation["gates"] = gates
    validation["overall_pass"] = all(gates.values())
    validation["status"] = "VALIDATED" if validation["overall_pass"] else "FAILED"
    return primary, audit, validation


def prepare() -> dict[str, Any]:
    ensure_directories()
    protocol = frozen_protocol()
    if FROZEN_PROTOCOL_JSON.exists():
        existing = json.loads(FROZEN_PROTOCOL_JSON.read_text(encoding="utf-8"))
        if existing != protocol:
            raise RuntimeError("frozen protocol JSON differs from current protocol")
    else:
        write_json_atomic(FROZEN_PROTOCOL_JSON, protocol)

    primary, audit, validation = _reference_and_validation()
    if not validation["overall_pass"]:
        write_json_atomic(REFERENCE_VALIDATION, validation)
        raise RuntimeError("BLOCKED_REFERENCE_VALIDATION")

    parameter_rows: list[dict[str, Any]] = []
    reference_rows: list[dict[str, Any]] = []
    for dimension in DIMENSIONS:
        basis = frozen_basis(dimension)
        root = PARAMETER_DIR / f"{PDE}_d{dimension}"
        _freeze_npy(root / "basis.npy", basis)
        _freeze_npy(root / "v.npy", basis[:, 0])
        parameter_rows.append(
            {
                "dimension": dimension,
                "basis_file": str(root / "basis.npy"),
                "basis_file_sha256": file_sha256(root / "basis.npy"),
                "v_file": str(root / "v.npy"),
                "v_file_sha256": file_sha256(root / "v.npy"),
                "orthogonality_max_abs": float(
                    np.max(np.abs(basis.T @ basis - np.eye(2)))
                ),
            }
        )

        t, x = _query_inputs(dimension)
        coordinates = x @ basis
        reference = primary.value(t, coordinates)
        reference_gradient = primary.gradient(t, coordinates)
        reference_z = reduced_to_ambient_z(
            dimension, math.sqrt(2.0) * reference_gradient
        )
        audit_value = audit.value(t, coordinates)
        audit_z = reduced_to_ambient_z(
            dimension, math.sqrt(2.0) * audit.gradient(t, coordinates)
        )
        uncertainty = np.abs(reference - audit_value)
        uncertainty_z = np.linalg.norm(reference_z - audit_z, axis=1)
        bank = query_path(dimension)
        for name, value in (
            ("t.npy", t),
            ("x.npy", x),
            ("reference.npy", reference),
            ("reference_z.npy", reference_z),
            ("reference_uncertainty.npy", uncertainty),
            ("reference_z_uncertainty.npy", uncertainty_z),
        ):
            _freeze_npy(bank / name, value)
        bank_npz = bank / "bank.npz"
        if bank_npz.exists():
            with np.load(bank_npz) as existing:
                if any(
                    not np.array_equal(existing[key], value)
                    for key, value in {
                        "t": t,
                        "x": x,
                        "reference": reference,
                        "reference_z": reference_z,
                    }.items()
                ):
                    raise RuntimeError(f"refusing to overwrite changed bank: {bank_npz}")
        else:
            write_npz_atomic(
                bank_npz, t=t, x=x, reference=reference, reference_z=reference_z
            )
        manifest = {
            "pde": PDE,
            "dimension": dimension,
            "points": QUERY_COUNT,
            "origin_row": 0,
            "query_seed_sequence": [BASE_SEED, dimension, QUERY_COUNT],
            "t_sha256": file_sha256(bank / "t.npy"),
            "x_sha256": file_sha256(bank / "x.npy"),
            "reference_sha256": file_sha256(bank / "reference.npy"),
            "reference_z_sha256": file_sha256(bank / "reference_z.npy"),
            "bank_sha256": file_sha256(bank_npz),
            "reference_uncertainty_rmse": float(
                np.sqrt(np.mean(np.square(uncertainty)))
            ),
            "reference_uncertainty_max_abs": float(np.max(uncertainty)),
            "reference_z_uncertainty_vector_rmse": float(
                np.sqrt(np.mean(np.square(uncertainty_z)))
            ),
        }
        write_json_atomic(bank / "manifest.json", manifest)
        reference_rows.append(manifest)

    parameter_manifest = {
        "status": "VALIDATED",
        "overall_pass": all(row["orthogonality_max_abs"] <= 1e-12 for row in parameter_rows),
        "rows": parameter_rows,
    }
    write_json_atomic(PARAMETER_MANIFEST, parameter_manifest)
    validation["banks"] = reference_rows
    validation["parameter_manifest_sha256"] = sha256(
        json.dumps(parameter_manifest, sort_keys=True).encode("utf-8")
    ).hexdigest()
    write_json_atomic(REFERENCE_VALIDATION, validation)
    return validation


def main() -> None:
    result = prepare()
    print(json.dumps({
        "status": result["status"],
        "value_refinement_rmse": result["value_refinement_rmse"],
        "pde_residual_rmse": result["pde_residual_rmse"],
        "banks": len(result["banks"]),
    }, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["prepare"]
