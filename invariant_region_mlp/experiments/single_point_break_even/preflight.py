"""H100 source-fidelity, bank-integrity, and small numerical smoke audit."""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import inspect
import json
from pathlib import Path

import numpy as np
import torch

from invariant_region_mlp.experiments.single_point_32_ice import (
    deep_bsde as source_db,
    fh_mlp as source_fh,
    tangent_actor as source_ta,
)
from invariant_region_mlp.experiments.single_point_32_ice.benchmarks import (
    make_benchmark,
)

from .deep_bsde import CONFIG, _train_one
from .fh_mlp import _solve_point
from .protocol import (
    DIMENSION,
    RESULTS_DIR,
    REVISION,
    SOURCE_ROOT,
    array_sha256,
    bank_path,
    file_sha256,
    load_bank,
    stable_seed,
    write_json_atomic,
)


def run(device: torch.device) -> dict:
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("CUDA preflight required")
    gpu = torch.cuda.get_device_name(device)
    if "H100" not in gpu.upper():
        raise RuntimeError(f"H100 required, got {gpu}")
    calibration = load_bank("calibration", require_truth=True)
    final = load_bank("final", require_truth=True)
    timing = load_bank("timing", require_truth=False)
    if "reference" in timing:
        raise RuntimeError("timing truth leak")
    reference_audit = json.loads(
        (RESULTS_DIR / "reference_audit.json").read_text(encoding="utf-8")
    )
    if reference_audit["status"] != "VALIDATED":
        raise RuntimeError("reference gate did not pass")

    benchmark = source_fh._UnprojectedDriver(make_benchmark("allen_cahn", DIMENSION))
    fh_seed = stable_seed("preflight", "fh")
    fh_first, fh_stats_first = _solve_point(
        benchmark,
        0.0,
        calibration["x"][0],
        n=2,
        M=4,
        seed=fh_seed,
        device=device,
    )
    fh_second, fh_stats_second = _solve_point(
        benchmark,
        0.0,
        calibration["x"][0],
        n=2,
        M=4,
        seed=fh_seed,
        device=device,
    )
    fh_replay = abs(fh_first - fh_second)

    deep = _train_one(
        calibration["x"][0],
        full_seed=stable_seed("preflight", "deep_bsde"),
        capture_updates=(1,),
        max_updates=1,
        device=device,
        resume_path=None,
        checkpoint_path=None,
    )
    source_files = {
        "tangent_actor": Path(inspect.getsourcefile(source_ta) or ""),
        "fh_mlp": Path(inspect.getsourcefile(source_fh) or ""),
        "deep_bsde": Path(inspect.getsourcefile(source_db) or ""),
    }
    checks = {
        "gpu_is_h100": "H100" in gpu.upper(),
        "reference_gate": reference_audit["status"] == "VALIDATED",
        "bank_shapes": (
            calibration["x"].shape == (32, 100)
            and final["x"].shape == (64, 100)
            and timing["x"].shape == (524288, 100)
        ),
        "timing_has_no_truth": "reference" not in timing,
        "fh_draw_replay_exact": fh_replay == 0.0,
        "fh_nonfinite_zero": (
            fh_stats_first.nonfinite_count == 0 and fh_stats_second.nonfinite_count == 0
        ),
        "deep_one_update_finite": (
            deep["status"] == "complete"
            and deep["nonfinite_count"] == 0
            and np.isfinite(deep["candidates"][0]["y0"])
        ),
        "deep_config_frozen": (
            CONFIG.intervals == 20
            and CONFIG.hidden_widths == (110, 110)
            and CONFIG.batch_size == 64
            and CONFIG.learning_rate == 5e-4
        ),
    }
    payload = {
        "revision": REVISION,
        "status": "VALIDATED" if all(checks.values()) else "FAILED",
        "overall_pass": all(checks.values()),
        "checks": checks,
        "gpu": gpu,
        "torch_version": torch.__version__,
        "cuda_version": torch.version.cuda,
        "banks": {
            name: {
                "file_sha256": file_sha256(bank_path(name)),
                "input_sha256": array_sha256(bank["t"], bank["x"]),
            }
            for name, bank in (
                ("calibration", calibration), ("final", final), ("timing", timing)
            )
        },
        "source_files": {
            name: {
                "path": str(path.relative_to(SOURCE_ROOT)),
                "sha256": file_sha256(path),
            }
            for name, path in source_files.items()
        },
        "fh_replay_max_abs": fh_replay,
        "fh_value": fh_first,
        "deep_smoke": deep,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(RESULTS_DIR / "preflight_audit.json", payload)
    if not payload["overall_pass"]:
        raise RuntimeError("break-even CUDA preflight failed")
    return payload


if __name__ == "__main__":
    print(json.dumps(run(torch.device("cuda")), indent=2, sort_keys=True))

