from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .paired_bank_repair import (
    array_sha256,
    canonical_bank_arrays,
    file_sha256,
    verify_canonical_prediction,
)
from .value_metrics import value_metric_payload


def _write_prediction(path: Path, t: np.ndarray, x: np.ndarray) -> dict[str, float]:
    truth = np.array([0.25, 0.75], dtype=np.float64)
    prediction = np.array([0.20, 0.90], dtype=np.float64)
    np.savez_compressed(path, t=t, x=x, truth=truth, prediction=prediction)
    return value_metric_payload(prediction, truth)


def test_rmae_is_relative_l1_not_mean_percentage_error(tmp_path: Path) -> None:
    t = np.array([0.0, 0.5], dtype=np.float64)
    x = np.array([[1.0], [2.0]], dtype=np.float64)
    metrics = _write_prediction(tmp_path / "prediction.npz", t, x)

    assert np.isclose(metrics["value_mae"], 0.1)
    assert np.isclose(metrics["value_rmae"], 0.2)


def test_numerical_array_hash_is_independent_of_npz_container(tmp_path: Path) -> None:
    t = np.array([0.0, 0.5], dtype=np.float64)
    x = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    first = tmp_path / "first.npz"
    second = tmp_path / "second.npz"
    np.savez(first, t=t, x=x, extra=np.array([1]))
    np.savez_compressed(second, x=x, t=t, extra=np.array([2]))

    assert file_sha256(first) != file_sha256(second)
    assert array_sha256(t, x) == array_sha256(t.copy(), x.copy())


def test_new_dimension_ou_reference_uses_one_exact_float64_path() -> None:
    import torch

    from .benchmarks import OUGMMBenchmark
    from .protocol import RESULTS_DIR

    canonical = canonical_bank_arrays("ou_gmm", 20)
    benchmark = OUGMMBenchmark.from_npz(RESULTS_DIR / "inputs/ou_gmm_d20.npz")
    expected = benchmark.exact_value(
        torch.as_tensor(canonical["t"], dtype=torch.float64),
        torch.as_tensor(canonical["x"], dtype=torch.float64),
    ).numpy()

    assert canonical["truth"].dtype == np.float64
    assert np.array_equal(canonical["truth"], expected)
    assert len(canonical["truth"]) == 960


def test_strict_pairing_audit_opens_arrays_and_detects_tampering(
    tmp_path: Path, monkeypatch,
) -> None:
    from . import verify_paired_banks as verifier

    monkeypatch.setattr(verifier, "PDES", ("toy",))
    monkeypatch.setattr(verifier, "METHODS", ("actor", "baseline"))
    t = np.array([0.0, 0.5], dtype=np.float64)
    x = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    source = tmp_path / "source.npz"
    _write_prediction(source, t, x)
    truth = np.array([0.25, 0.75], dtype=np.float64)
    monkeypatch.setattr(
        verifier,
        "canonical_bank_arrays",
        lambda pde, dimension: {"t": t, "x": x, "truth": truth},
    )
    rows = []
    for method in verifier.METHODS:
        prediction_path = tmp_path / f"{method}.npz"
        metrics = _write_prediction(prediction_path, t, x)
        actual = verify_canonical_prediction(prediction_path)
        rows.append({
            "pde": "toy", "dimension": 100, "method": method,
            "seed": 0, "profile_role": "primary", "states": actual["states"],
            "test_bank_sha256": actual["test_bank_array_sha256"],
            "test_bank_array_sha256": actual["test_bank_array_sha256"],
            "reference_array_sha256": actual["reference_array_sha256"],
            "prediction_file": str(prediction_path),
            "prediction_file_sha256": actual["prediction_file_sha256"],
            "source_prediction_file": str(source),
            "source_prediction_file_sha256": file_sha256(source),
            **metrics,
        })

    csv_path = tmp_path / "MAIN_RESULTS.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    assert verifier.audit(csv_path)["status"] == "PASS"

    # Change an actual state array while leaving every CSV metadata string intact.
    _write_prediction(Path(rows[1]["prediction_file"]), t + np.array([0.0, 0.125]), x)
    failed = verifier.audit(csv_path)
    assert failed["status"] == "FAIL"
    baseline = next(
        record for record in failed["results"]["toy"]["records"]
        if record["method"] == "baseline"
    )
    assert any("test_bank" in message or "prediction_file" in message for message in baseline["errors"])
