"""Shared whole-function value metrics computed on one frozen test bank."""

from __future__ import annotations

import numpy as np


VALUE_METRIC_NAMES = (
    "value_rmse",
    "value_relative_l2",
    "value_bias",
    "value_mae",
    "value_rmae",
)


def value_metric_payload(
    prediction: np.ndarray, truth: np.ndarray,
) -> dict[str, float]:
    """Return absolute and relative value errors from paired pointwise arrays.

    RMAE is the relative L1 error

        sum_i |prediction_i - truth_i| / sum_i |truth_i|,

    equivalently MAE divided by the mean absolute truth.  This aggregate
    normalization remains well-defined for banks containing truth values near
    zero, unlike a mean pointwise percentage error.
    """
    prediction = np.asarray(prediction, dtype=np.float64).reshape(-1)
    truth = np.asarray(truth, dtype=np.float64).reshape(-1)
    if prediction.shape != truth.shape or prediction.size == 0:
        raise ValueError(
            f"prediction/truth must be nonempty paired vectors, got "
            f"{prediction.shape} and {truth.shape}"
        )
    if not np.isfinite(prediction).all() or not np.isfinite(truth).all():
        raise ValueError("prediction/truth contain nonfinite values")
    error = prediction - truth
    absolute_error = np.abs(error)
    mae = float(np.mean(absolute_error))
    truth_mae = float(np.mean(np.abs(truth)))
    return {
        "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
        "value_relative_l2": float(
            np.linalg.norm(error) / max(np.linalg.norm(truth), np.finfo(float).tiny)
        ),
        "value_bias": float(np.mean(error)),
        "value_mae": mae,
        "value_rmae": float(mae / max(truth_mae, np.finfo(float).tiny)),
    }
