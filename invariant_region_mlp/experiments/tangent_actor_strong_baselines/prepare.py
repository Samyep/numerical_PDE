"""Build fixed references/test banks and record source/protocol metadata."""

from __future__ import annotations

import json
import platform
import subprocess
import time

import torch

from .common import (
    BENCHMARKS,
    OFFICIAL_SOURCES,
    RESULTS_DIR,
    atomic_json,
    driver_for,
    ensure_dirs,
    evaluation_data,
    load_or_build_reference,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.reference2d import (
    validate_periodic_reference_2d,
)


def main() -> None:
    ensure_dirs()
    reference_rows = []
    for benchmark in BENCHMARKS:
        driver = driver_for(benchmark)
        reference = load_or_build_reference(benchmark)
        data = evaluation_data(benchmark)
        # Re-run the independent coarse/fine audit, kept outside method timing.
        _, audit = validate_periodic_reference_2d(driver)
        reference_rows.append(
            {
                "benchmark": benchmark,
                **audit,
                "test_value_min": float(data["reference_original"][data["test_indices"]].min()),
                "test_value_max": float(data["reference_original"][data["test_indices"]].max()),
                "reference_grid_n": reference.n,
            }
        )
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = "unknown"
    atomic_json(RESULTS_DIR / "official_source_review.json", OFFICIAL_SOURCES)
    atomic_json(RESULTS_DIR / "numerical_reference_validation.json", reference_rows)
    atomic_json(
        RESULTS_DIR / "runtime_environment.json",
        {
            "host": platform.node(),
            "python": platform.python_version(),
            "torch": torch.__version__,
            "cuda_runtime": torch.version.cuda,
            "cuda_available": torch.cuda.is_available(),
            "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
            "git_commit_at_prepare": commit,
            "unix_time": time.time(),
        },
    )
    print(json.dumps({"status": "prepared", "benchmarks": list(BENCHMARKS)}))


if __name__ == "__main__":
    main()
