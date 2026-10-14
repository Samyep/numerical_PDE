"""Build and validate independent one-dimensional numerical references."""

from __future__ import annotations

from .common import RESULTS_DIR, atomic_json, ensure_dirs
from .drivers import make_driver
from .reference import validate_periodic_reference


def main() -> None:
    ensure_dirs()
    rows = {}
    for name in ("allen_cahn", "joint_concave", "joint_convex"):
        _, diagnostics = validate_periodic_reference(make_driver(name))
        rows[name] = diagnostics
        print(name, diagnostics, flush=True)
    atomic_json(
        RESULTS_DIR / "reference_validation.json",
        {
            "method": "independent periodic Fourier spectral Strang solver with explicit midpoint reaction",
            "forward_time_variable": "tau = T - t",
            "fine": {"spatial_points": 512, "internal_dt": 0.000125, "stored_dt": 0.0005},
            "coarse": {"spatial_points": 256, "internal_dt": 0.00025, "stored_dt": 0.001},
            "problems": rows,
        },
    )


if __name__ == "__main__":
    main()
