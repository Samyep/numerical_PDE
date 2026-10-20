"""Paired high-curvature FH-MLP vs *frozen* Tangent-Actor comparison.

This is a supplemental GPU experiment, not part of the ICML 140-cell grid.
It reuses the exact post-freeze Actor test bank and numerical reference and
does not retrain any Actor or select baseline hyperparameters by test error.

Examples:
  python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.strong_nonlinearity_fh_pair --dry-run
  python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.strong_nonlinearity_fh_pair --n 3 --M 32
  python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.strong_nonlinearity_fh_pair --n 4 --M 96
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import time

import numpy as np
import torch

from invariant_region_mlp.experiments.nonconvex_tangent_actor.common import (
    PREDICTION_DIR as OLD_PREDICTIONS,
    RESULTS_DIR as OLD_RESULTS,
    array_sha256,
    metrics,
)
from invariant_region_mlp.experiments.nonconvex_tangent_actor.configuration import ExperimentSpec
from invariant_region_mlp.experiments.nonconvex_tangent_actor.full_history_baseline import RawDriver
from invariant_region_mlp.experiments.nonconvex_tangent_actor.reference2d import load_reference
from invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.full_history_frontier import (
    _predict_tiled,
)

RESULT_DIR = OLD_RESULTS.parent / "strong_nonlinearity_fh_paired"
# Configurations are prespecified BEFORE any of their new RMSEs are observed.
GRID = ((3, 32), (3, 64), (3, 96), (3, 128),
        (4, 16), (4, 32), (4, 64), (4, 96))


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def actor_endpoint(spec: ExperimentSpec) -> dict:
    path = OLD_RESULTS / f"post_freeze_{spec.run_id}.json"
    if not path.exists():
        raise FileNotFoundError(f"Frozen post-freeze Actor report missing: {path}")
    payload = json.loads(path.read_text())
    stored = payload["spec"]
    for name, expected in (
        ("family", spec.family), ("d", spec.d), ("lam", spec.lam),
        ("eta", spec.eta), ("budget", spec.budget), ("seed", spec.seed)
    ):
        if isinstance(expected, float):
            if not math.isclose(float(stored[name]), expected, rel_tol=0, abs_tol=1e-12):
                raise RuntimeError(f"Actor spec mismatch: {name}")
        elif stored[name] != expected:
            raise RuntimeError(f"Actor spec mismatch: {name}")
    row = payload["selected_row"]
    if not math.isfinite(float(row["mcv_value_rmse"])):
        raise ValueError("Actor endpoint is nonfinite")
    costs = []
    table = OLD_RESULTS / "BASELINE_COMPARISON.csv"
    if table.exists():
        with table.open(newline="") as f:
            for line in csv.DictReader(f):
                if (line["method"] == "Tangent-Actor + MCV"
                    and line["family"] == spec.family
                    and int(line["d"]) == spec.d
                    and int(line["seed"]) == spec.seed
                    and line["budget"] == spec.budget
                    and math.isclose(float(line["lam"]), spec.lam, abs_tol=1e-12)
                    and math.isclose(float(line["eta"]), spec.eta, abs_tol=1e-12)):
                    costs.append(float(line["total_seconds"]))
    return {
        "selected_actor": payload["selection"]["selected_actor"],
        "value_rmse": float(row["mcv_value_rmse"]),
        "relative_l2": float(row["mcv_relative_l2"]),
        "paths": int(row["paths"]), "euler_steps": int(row["time_steps"]),
        "inference_seconds_512": float(row["inference_seconds"]),
        "total_seconds_512": costs[0] if len(costs) == 1 else None,
        "source_json": str(path),
        "source_json_sha256": sha_file(path),
        "reference_path": payload["reference_path"],
        "source_bank_array_sha256": payload["test_bank_sha256"],
    }


def get_bank(spec: ExperimentSpec, original: dict):
    path = OLD_PREDICTIONS / f"post_freeze_{spec.run_id}.npz"
    if not path.exists():
        raise FileNotFoundError(
            f"The original frozen prediction NPZ is only on the HPC workspace: {path}"
        )
    with np.load(path, allow_pickle=False) as data:
        t = np.asarray(data["test_t"], dtype=np.float64)
        x = np.asarray(data["test_x"], dtype=np.float64)
        truth = np.asarray(data["truth"], dtype=np.float64)
    if t.shape != (512,) or x.shape != (512, spec.d) or truth.shape != (512,):
        raise ValueError(f"Frozen test bank shape unexpected: {t.shape}, {x.shape}, {truth.shape}")
    if not all(np.isfinite(item).all() for item in (t, x, truth)):
        raise ValueError("Nonfinite frozen bank/reference")
    frozen_ref_path = Path(original["reference_path"])
    if not frozen_ref_path.exists():
        raise FileNotFoundError(f"Independent frozen numerical reference missing: {frozen_ref_path}")
    independent_truth = load_reference(frozen_ref_path).value(t, x @ spec.driver().B)
    if not np.allclose(truth, independent_truth, rtol=0.0, atol=1e-9):
        raise RuntimeError("Actor bank truth disagrees with independent 2D reference")
    bank_hash = array_sha256(t, x)
    if bank_hash != original["source_bank_array_sha256"]:
        raise RuntimeError("Frozen Actor test bank array hash mismatch")
    return t, x, truth, path, bank_hash


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", choices=("primary", "coupled"), default="primary")
    parser.add_argument("--lam", type=float, default=0.5)
    parser.add_argument("--eta", type=float, default=0.0)
    parser.add_argument("--budget", choices=("pilot", "main"), default="pilot")
    parser.add_argument("--d", type=int, default=100)
    parser.add_argument("--actor-seed", type=int, default=0)
    parser.add_argument("--baseline-seed", type=int, default=0)
    parser.add_argument("--n", type=int, default=3)
    parser.add_argument("--M", type=int, default=32)
    parser.add_argument("--max-seconds", type=float, default=1800.0)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if (args.n, args.M) not in GRID:
        parser.error(f"(n, M) must be one of {GRID}; no post-hoc extensions")
    if args.family == "primary" and args.eta != 0:
        parser.error("primary family requires eta=0")
    if args.max_seconds <= 0:
        parser.error("--max-seconds must be positive")

    spec = ExperimentSpec(
        family=args.family, d=args.d, lam=args.lam,
        eta=args.eta, seed=args.actor_seed, budget=args.budget, estimator="mcv"
    )
    a = actor_endpoint(spec)
    t, x, truth, bank, bank_hash = get_bank(spec, a)
    descriptor = {
        "pde": spec.run_id, "n": args.n, "M": args.M,
        "baseline_seed": args.baseline_seed, "actor_seed": args.actor_seed,
        "states": len(t), "numerical_test_bank_sha256": bank_hash,
        "truth_array_sha256": array_sha256(truth),
        "frozen_npz_sha256": sha_file(bank), "frozen_npz": str(bank),
        "actor": a, "status": "dry_run" if args.dry_run else "RUNNING",
        "estimated_960_points_note": "Only a linear-in-states extrapolation, not a measured 960-point run",
        "reference_is_independent": True, "test_truth_used_for_selection": False,
    }
    if args.dry_run:
        print(json.dumps(descriptor, indent=2, sort_keys=True))
        return

    device = torch.device(args.device)
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("Full benchmark requires CUDA; CPU is only suitable for a separate smoke test")

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    stem = f"{spec.run_id}_fh_n{args.n}_M{args.M}_s{args.baseline_seed}"
    output_json = RESULT_DIR / f"{stem}.json"
    output_npz = RESULT_DIR / f"{stem}.npz"
    if output_json.exists() or output_npz.exists():
        raise FileExistsError(f"Refusing overwrite: {output_json} or {output_npz}")

    def _timeout(signum, frame):
        del signum, frame
        raise TimeoutError(f"Exceeded {args.max_seconds}s per-cell time cap")

    signal.signal(signal.SIGALRM, _timeout)
    signal.setitimer(signal.ITIMER_REAL, args.max_seconds)
    started = time.perf_counter()
    try:
        prediction, metadata = _predict_tiled(
            RawDriver(spec.driver()), t, x, n=args.n, M=args.M,
            seed=args.baseline_seed, device=device
        )
        runtime = time.perf_counter() - started
        if not metadata["finite"] or not np.all(np.isfinite(prediction)):
            raise FloatingPointError("Raw FH-MLP produced nonfinite values")
        result = metrics(prediction, truth)
        np.savez_compressed(output_npz, prediction=prediction)
        descriptor.update({
            "status": "COMPLETE", "fh": result,
            "fh_measured_seconds_512": float(metadata["wall_seconds"]),
            "fh_total_elapsed_seconds": runtime,
            "fh_estimated_seconds_960": float(metadata["wall_seconds"]) * 960.0 / len(t),
            "fh_960_estimate_is_measured": False,
            "output_npz": str(output_npz),
            "output_npz_sha256": sha_file(output_npz),
            "memory_tiling": True,
            "memory_bytes_peak": metadata.get("peak_gpu_bytes"),
            "raw_fh_sample_stats": metadata.get("stats"),
        })
        if a["total_seconds_512"] is not None:
            descriptor["actor_estimated_seconds_960"] = (
                a["total_seconds_512"]
                - a["inference_seconds_512"]
                + a["inference_seconds_512"] * 960.0 / len(t)
            )
    except TimeoutError as exc:
        descriptor.update({"status": "OUT_OF_TIME", "reason": str(exc),
                           "elapsed_seconds": time.perf_counter() - started})
    except Exception as exc:
        descriptor.update({"status": "FAILED", "reason": repr(exc),
                           "elapsed_seconds": time.perf_counter() - started})
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    temp = output_json.with_name(output_json.name + f".tmp.{os.getpid()}")
    temp.write_text(json.dumps(descriptor, indent=2, sort_keys=True) + "\n")
    temp.replace(output_json)
    print(json.dumps({
        "status": descriptor["status"], "output": str(output_json),
        "fh_rmse": (descriptor.get("fh") or {}).get("rmse"),
        "actor_rmse": a["value_rmse"], "seconds": descriptor.get("fh_measured_seconds_512")
    }, sort_keys=True), flush=True)
    if descriptor["status"] != "COMPLETE":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
