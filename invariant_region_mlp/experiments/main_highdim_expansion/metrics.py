"""Metrics derived directly from saved whole-function predictions."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def value_metrics_from_prediction_file(path: Path | str) -> dict[str, float]:
    """Recompute RMSE, relative L2, and rMAE from a prediction artifact.

    Tangent-Actor's B6 artifact stores the primary estimate as ``prediction``;
    the Dense Fourier artifact uses ``mcv``. Baseline artifacts use
    ``prediction``.
    """

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    with np.load(path, allow_pickle=False) as payload:
        if "prediction" in payload.files:
            prediction = np.asarray(payload["prediction"], dtype=np.float64)
        elif "mcv" in payload.files:
            prediction = np.asarray(payload["mcv"], dtype=np.float64)
        else:
            raise KeyError(f"{path} contains neither 'prediction' nor 'mcv'")
        if "truth" not in payload.files:
            raise KeyError(f"{path} contains no 'truth' array")
        truth = np.asarray(payload["truth"], dtype=np.float64)
    if prediction.shape != truth.shape:
        raise ValueError(
            f"prediction/truth shape mismatch in {path}: "
            f"{prediction.shape} != {truth.shape}"
        )
    if prediction.size == 0 or not np.isfinite(prediction).all() or not np.isfinite(truth).all():
        raise ValueError(f"empty or non-finite prediction artifact: {path}")
    error = prediction - truth
    truth_l2 = float(np.linalg.norm(truth.ravel()))
    truth_l1 = float(np.abs(truth).sum())
    if truth_l2 <= 0.0 or truth_l1 <= 0.0:
        raise ValueError(f"zero truth norm in {path}")
    return {
        "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
        "value_relative_l2": float(np.linalg.norm(error.ravel()) / truth_l2),
        "value_rmae": float(np.abs(error).sum() / truth_l1),
    }
