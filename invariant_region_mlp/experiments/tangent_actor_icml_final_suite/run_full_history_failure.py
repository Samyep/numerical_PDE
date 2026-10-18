"""Reproduce the raw full-history MLP dimension-failure curve.

This wrapper calls the literal Hutzenthaler--Nguyen recursion already audited
in ``counterexample_rescue``.  It adds the frozen six-dimension grid, resumable
per-dimension artifacts, and value-RMSE bootstrap intervals needed by Figure 1.
No projection or Tangent-Actor component enters the reported raw baseline.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
from typing import Any

import numpy as np
import pandas as pd

from invariant_region_mlp.experiments.counterexample_rescue.counterexample_rescue import (
    MODES,
    run_config,
)

from .protocol import RESULTS_DIR


RAW_DIR = RESULTS_DIR / "raw"
PREDICTION_DIR = RESULTS_DIR / "predictions"
SOURCE_ROOT = Path(__file__).resolve().parents[3]
SOURCE_ALGORITHM = (
    SOURCE_ROOT / "invariant_region_mlp/experiments/counterexample_rescue/counterexample_rescue.py"
)


def _git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT,
        check=True, capture_output=True, text=True,
    ).stdout.strip()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(SOURCE_ROOT))
    except ValueError:
        return str(path)


def _value_rmse_bootstrap(usq: np.ndarray, *, seed: int, n_boot: int = 4000) -> list[float]:
    rng = np.random.default_rng(seed)
    estimates = np.empty(n_boot, dtype=np.float64)
    for start in range(0, n_boot, 250):
        stop = min(start + 250, n_boot)
        indices = rng.integers(0, len(usq), size=(stop - start, len(usq)))
        estimates[start:stop] = np.sqrt(usq[indices].mean(axis=1))
    return np.quantile(estimates, (0.025, 0.975)).tolist()


def run_dimension(d: int, *, reps: int, batch: int, seed: int,
                  force: bool = False) -> dict[str, Any]:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PREDICTION_DIR.mkdir(parents=True, exist_ok=True)
    path = RAW_DIR / f"full_history_mlp_failure_d{d}_n3_s{seed}.json"
    arrays = PREDICTION_DIR / f"full_history_mlp_failure_d{d}_n3_s{seed}.npz"
    if path.exists() and arrays.exists() and not force:
        payload = json.loads(path.read_text())
        if payload.get("status") == "complete":
            return payload

    started = time.perf_counter()
    rows, audit, states, squared, value_squared = run_config(
        d=d, n=3, m=3, reps=reps, batch=batch, seed=seed,
    )
    raw_index = MODES.index("raw")
    raw = dict(rows[raw_index])
    raw["value_rmse_boot95"] = _value_rmse_bootstrap(
        value_squared[:, raw_index], seed=seed + 1009 * d,
    )
    np.savez_compressed(
        arrays,
        raw_state=states[:, raw_index, :],
        raw_state_squared_error=squared[:, raw_index],
        raw_value_squared_error=value_squared[:, raw_index],
    )
    payload = {
        "method": "raw_full_history_mlp",
        "pde": "mlp_failure", "dimension": d, "seed": seed,
        "level": 3, "branching": 3, "replicates": reps,
        "row": raw, "algorithm_audit": audit,
        "source_algorithm": str(SOURCE_ALGORITHM.relative_to(SOURCE_ROOT)),
        "source_algorithm_sha256": _sha256(SOURCE_ALGORITHM),
        "wrapper_git_head": _git_head(),
        "prediction_artifact": _display_path(arrays),
        "invocation_seconds": time.perf_counter() - started,
        "projection_used_by_reported_raw_method": False,
        "status": "complete",
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def aggregate(records: list[dict[str, Any]]) -> Path:
    rows = []
    for payload in records:
        item = payload["row"]
        rows.append({
            "pde": "mlp_failure", "dimension": payload["dimension"],
            "method": "raw_full_history_mlp", "seed": payload["seed"],
            "level": payload["level"], "branching": payload["branching"],
            "replicates": payload["replicates"],
            "value_rmse": item["value_rmse"], "value_bias": item["value_bias"],
            "value_rmse_ci95_low": item["value_rmse_boot95"][0],
            "value_rmse_ci95_high": item["value_rmse_boot95"][1],
            "state_rmse": item["state_rmse"],
            "gradient_rmse": item["gradient_rmse"],
            "finite_fraction": item["finite_fraction"],
            "invocation_seconds": payload["invocation_seconds"],
            "source_artifact": (
                f"invariant_region_mlp/results/tangent_actor_icml_final_suite/raw/"
                f"full_history_mlp_failure_d{payload['dimension']}_n3_s{payload['seed']}.json"
            ),
            "source_algorithm_sha256": payload["source_algorithm_sha256"],
            "status": payload["status"],
        })
    output = RESULTS_DIR / "MLP_FAILURE_FULL_HISTORY.csv"
    pd.DataFrame(rows).sort_values("dimension").to_csv(output, index=False)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimensions", default="20,50,100,200,400,1000")
    parser.add_argument("--reps", type=int, default=512)
    parser.add_argument("--batch", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20261004)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    dimensions = [int(value) for value in args.dimensions.split(",")]
    records = []
    for d in dimensions:
        record = run_dimension(
            d, reps=args.reps, batch=args.batch, seed=args.seed, force=args.force,
        )
        records.append(record)
        print(json.dumps({
            "dimension": d, "value_rmse": record["row"]["value_rmse"],
            "seconds": record["invocation_seconds"], "status": record["status"],
        }), flush=True)
    output = aggregate(records)
    print(json.dumps({"output": str(output), "rows": len(records)}), flush=True)


if __name__ == "__main__":
    main()
