"""Common accuracy metrics for stochastic and trained solvers."""

from __future__ import annotations

import math
from typing import Any

import numpy as np


def prediction_metrics(
    value: np.ndarray,
    z: np.ndarray | None,
    reference_value: np.ndarray,
    reference_z: np.ndarray,
) -> dict[str, Any]:
    observed = np.asarray(value, dtype=np.float64).reshape(-1)
    target = np.asarray(reference_value, dtype=np.float64).reshape(-1)
    finite = np.isfinite(observed)
    value_error = observed - target
    rmse = float(np.sqrt(np.mean(value_error**2))) if np.all(finite) else float("inf")
    scale = max(float(np.std(target)), np.finfo(np.float64).tiny)
    row: dict[str, Any] = {
        "value_rmse": rmse,
        "normalized_value_rmse": rmse / scale,
        "value_bias": float(np.mean(value_error)) if np.all(finite) else float("nan"),
        "value_mae": float(np.mean(np.abs(value_error))) if np.all(finite) else float("nan"),
        "value_nonfinite": int(np.size(observed) - np.count_nonzero(finite)),
    }
    if z is None:
        row.update(
            {
                "gradient_rmse": None,
                "relative_gradient_error": None,
                "gradient_bias_norm": None,
                "gradient_nonfinite": None,
            }
        )
    else:
        observed_z = np.asarray(z, dtype=np.float64)
        target_z = np.asarray(reference_z, dtype=np.float64)
        difference = observed_z - target_z
        z_finite = np.isfinite(observed_z)
        squared_error = float(np.mean(np.sum(difference**2, axis=-1)))
        reference_energy = float(np.mean(np.sum(target_z**2, axis=-1)))
        row.update(
            {
                "gradient_rmse": math.sqrt(squared_error),
                "relative_gradient_error": math.sqrt(
                    squared_error / max(reference_energy, np.finfo(np.float64).tiny)
                ),
                "gradient_bias_norm": float(np.linalg.norm(np.mean(difference, axis=0))),
                "gradient_nonfinite": int(observed_z.size - np.count_nonzero(z_finite)),
            }
        )
    return row
