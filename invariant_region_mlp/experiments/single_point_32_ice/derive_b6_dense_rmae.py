"""Derive RMSE and relative MAE for the 40 frozen B6/Dense TA runs.

This is a post-training metric derivation.  It reads the already audited
pointwise predictions in ``ta_results.csv`` and never launches a solver.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .protocol import QUERY_COUNT, RESULTS_DIR, SOURCE_ROOT


TARGET_PDES = ("b6_fullrank", "dense_fourier")
TARGET_DIMENSIONS = (20, 50, 100, 200, 400)
TARGET_REPLICATES = (0, 1, 2, 3)
TARGET_TIME_STEPS = 64
EXPECTED_RUNS = (
    len(TARGET_PDES) * len(TARGET_DIMENSIONS) * len(TARGET_REPLICATES)
)

RUN_METRICS = RESULTS_DIR / "b6_dense_ta_40_run_metrics.csv"
AGGREGATE_METRICS = RESULTS_DIR / "b6_dense_ta_aggregate_metrics.csv"
AUDIT = RESULTS_DIR / "b6_dense_ta_rmae_audit.json"


def relative_mae(prediction: np.ndarray, reference: np.ndarray) -> float:
    """Return sum(abs(error)) / sum(abs(reference))."""
    prediction = np.asarray(prediction, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    denominator = float(np.sum(np.abs(reference)))
    if denominator <= np.finfo(np.float64).tiny:
        raise ValueError("relative MAE is undefined for a zero L1 reference")
    return float(np.sum(np.abs(prediction - reference)) / denominator)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative_path(path: Path) -> str:
    return path.resolve().relative_to(SOURCE_ROOT.resolve()).as_posix()


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False, lineterminator="\n")
    temporary.replace(path)


def run() -> dict:
    source = RESULTS_DIR / "ta_results.csv"
    rows = pd.read_csv(source)
    selected = rows.loc[
        rows["pde"].isin(TARGET_PDES)
        & rows["dimension"].isin(TARGET_DIMENSIONS)
        & rows["inference_replicate"].isin(TARGET_REPLICATES)
        & (rows["time_steps"] == TARGET_TIME_STEPS)
        & (rows["method"] == "Tangent-Actor")
        & (rows["status"] == "complete")
    ].copy()

    run_rows: list[dict] = []
    aggregate_rows: list[dict] = []
    grouped = selected.groupby(
        ["pde", "dimension", "inference_replicate"], sort=True
    )
    if grouped.ngroups != EXPECTED_RUNS:
        raise RuntimeError(f"expected {EXPECTED_RUNS} runs, found {grouped.ngroups}")

    for (pde, dimension, replicate), frame in grouped:
        frame = frame.sort_values("point")
        if len(frame) != QUERY_COUNT or frame["point"].nunique() != QUERY_COUNT:
            raise RuntimeError(
                f"{pde} d={dimension} replicate={replicate} does not have "
                f"exactly {QUERY_COUNT} points"
            )
        prediction = frame["prediction"].to_numpy(dtype=np.float64)
        reference = frame["reference"].to_numpy(dtype=np.float64)
        if not np.all(np.isfinite(prediction)) or not np.all(np.isfinite(reference)):
            raise RuntimeError(
                f"nonfinite value in {pde} d={dimension} replicate={replicate}"
            )
        error = prediction - reference
        run_rows.append({
            "pde": pde,
            "dimension": int(dimension),
            "method": "Tangent-Actor",
            "inference_replicate": int(replicate),
            "query_count": QUERY_COUNT,
            "time_steps": TARGET_TIME_STEPS,
            "paths": int(frame["paths"].iloc[0]),
            "selected_actor": frame["selected_actor"].iloc[0],
            "rmse": float(np.sqrt(np.mean(np.square(error)))),
            "mae": float(np.mean(np.abs(error))),
            "rmae": relative_mae(prediction, reference),
            "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(reference)),
            "bias": float(np.mean(error)),
            "nonfinite_count": 0,
        })

    run_metrics = pd.DataFrame(run_rows).sort_values(
        ["pde", "dimension", "inference_replicate"]
    ).reset_index(drop=True)

    for (pde, dimension), frame in selected.groupby(
        ["pde", "dimension"], sort=True
    ):
        frame = frame.sort_values(["inference_replicate", "point"])
        prediction = frame["prediction"].to_numpy(dtype=np.float64)
        reference = frame["reference"].to_numpy(dtype=np.float64)
        error = prediction - reference
        repeats = run_metrics.loc[
            (run_metrics["pde"] == pde)
            & (run_metrics["dimension"] == dimension),
            "rmse",
        ].to_numpy(dtype=np.float64)
        aggregate_rows.append({
            "pde": pde,
            "dimension": int(dimension),
            "method": "Tangent-Actor",
            "formal_repeats": len(TARGET_REPLICATES),
            "observations": len(frame),
            "rmse": float(np.sqrt(np.mean(np.square(error)))),
            "mae": float(np.mean(np.abs(error))),
            "rmae": relative_mae(prediction, reference),
            "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(reference)),
            "repeat_rmse_mean": float(np.mean(repeats)),
            "repeat_rmse_sd": float(np.std(repeats, ddof=1)),
        })

    aggregate = pd.DataFrame(aggregate_rows).sort_values(
        ["pde", "dimension"]
    ).reset_index(drop=True)
    _atomic_csv(run_metrics, RUN_METRICS)
    _atomic_csv(aggregate, AGGREGATE_METRICS)

    audit = {
        "status": "VALIDATED",
        "scope": list(TARGET_PDES),
        "definition": "sum(abs(prediction-reference)) / sum(abs(reference))",
        "source": _relative_path(source),
        "source_sha256": _sha256(source),
        "checks": {
            "expected_runs": EXPECTED_RUNS,
            "observed_runs": len(run_metrics),
            "query_count_per_run": QUERY_COUNT,
            "observed_point_rows": len(selected),
            "expected_point_rows": EXPECTED_RUNS * QUERY_COUNT,
            "all_metrics_finite": bool(
                np.isfinite(
                    run_metrics[["rmse", "mae", "rmae", "relative_l2", "bias"]]
                    .to_numpy(dtype=np.float64)
                ).all()
            ),
        },
        "outputs": {
            _relative_path(RUN_METRICS): _sha256(RUN_METRICS),
            _relative_path(AGGREGATE_METRICS): _sha256(AGGREGATE_METRICS),
        },
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    temporary = AUDIT.with_suffix(AUDIT.suffix + ".tmp")
    temporary.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(AUDIT)
    return audit


def main() -> None:
    print(json.dumps(run(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["relative_mae", "run"]
