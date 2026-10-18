"""Estimate post-gate one-seed GPU cost from completed 100D primary runs."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path

import pandas as pd

from .protocol import METHODS, PDES, RESULTS_DIR, dimensions_for


FALLBACK_SECONDS = {
    "tangent_actor": 1_800.0,
    "dpi": 900.0,
    "pinn_hte": 1_200.0,
    "sdgd_pinn": 12_000.0,
    "deep_bsde": 100.0,
}


def _scale(method: str, dimension: int) -> float:
    ratio = dimension / 100.0
    if method in {"tangent_actor", "dpi"}:
        return 0.50 + 0.50 * ratio
    if method == "pinn_hte":
        return 0.70 + 0.30 * ratio
    if method == "sdgd_pinn":
        return 0.60 + 0.40 * ratio
    if method == "deep_bsde":
        return ((dimension + 10.0) / 110.0) ** 2
    raise KeyError(method)


def _read_rows(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    frame = pd.read_csv(path)
    if "seed" in frame:
        frame = frame[pd.to_numeric(frame["seed"], errors="coerce") == 0]
    return frame


def _base_times() -> tuple[dict[tuple[str, str], float], dict[tuple[str, str], str]]:
    main = _read_rows(RESULTS_DIR / "MAIN_RESULTS.csv")
    native = _read_rows(RESULTS_DIR / "NATIVE_POINT_RESULTS.csv")
    values: dict[tuple[str, str], float] = {}
    origins: dict[tuple[str, str], str] = {}
    for frame, source in ((main, "100D primary whole-function run"), (native, "100D primary native-point run")):
        if frame.empty:
            continue
        for _, row in frame.iterrows():
            if int(row["dimension"]) != 100:
                continue
            method, pde = str(row["method"]), str(row["pde"])
            if method == "deep_bsde" and source.startswith("100D primary whole"):
                continue
            value = pd.to_numeric(row.get("total_seconds"), errors="coerce")
            if pd.notna(value) and float(value) > 0:
                values[(pde, method)] = float(value)
                origins[(pde, method)] = source
    for pde in (item.key for item in PDES):
        for method in METHODS:
            key = (pde, method)
            if key not in values:
                method_values = [seconds for (other_pde, other_method), seconds in values.items() if other_method == method]
                values[key] = float(pd.Series(method_values).median()) if method_values else FALLBACK_SECONDS[method]
                origins[key] = "method median fallback" if method_values else "predeclared method fallback"
    return values, origins


def build() -> dict[str, object]:
    base, origins = _base_times()
    rows = []
    for pde in PDES:
        for dimension in dimensions_for(pde.key):
            if dimension == 100:
                continue
            for method in METHODS:
                seconds = base[(pde.key, method)] * _scale(method, dimension)
                wide_extrapolation = dimension >= 400 or dimension <= 20
                low_factor, high_factor = (0.55, 2.20) if wide_extrapolation else (0.65, 1.75)
                rows.append({
                    "pde": pde.key, "dimension": dimension, "method": method,
                    "training_seed": 0, "base_100d_seconds": base[(pde.key, method)],
                    "base_origin": origins[(pde.key, method)],
                    "dimension_scale": _scale(method, dimension),
                    "estimated_gpu_hours": seconds / 3600.0,
                    "low_gpu_hours": seconds * low_factor / 3600.0,
                    "high_gpu_hours": seconds * high_factor / 3600.0,
                    "review_gate_status": "pending_independent_review",
                })
    frame = pd.DataFrame(rows)
    csv_path = RESULTS_DIR / "REMAINING_GPU_HOUR_ESTIMATE.csv"
    frame.to_csv(csv_path, index=False)
    by_method = {
        method: float(frame.loc[frame["method"] == method, "estimated_gpu_hours"].sum())
        for method in METHODS
    }
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "post-gate dimension sweep only; 100D and fidelity-anchor work excluded",
        "training_seeds": [0],
        "method_cells": len(frame),
        "point_estimate_gpu_hours": float(frame["estimated_gpu_hours"].sum()),
        "low_gpu_hours": float(frame["low_gpu_hours"].sum()),
        "high_gpu_hours": float(frame["high_gpu_hours"].sum()),
        "by_method_gpu_hours": by_method,
        "model": {
            "basis": "observed seed-0 100D primary wall time, one GPU per job",
            "dimension_scaling": {
                "tangent_actor_and_dpi": "0.5 + 0.5*d/100",
                "pinn_hte": "0.7 + 0.3*d/100",
                "sdgd_pinn": "0.6 + 0.4*d/100",
                "deep_bsde": "((d+10)/110)^2",
            },
            "intervals": "0.65x-1.75x normally; 0.55x-2.20x at d<=20 or d>=400",
            "queue_wait_excluded": True,
        },
        "review_gate": "pending_independent_review",
        "details_csv": str(csv_path),
    }
    (RESULTS_DIR / "REMAINING_GPU_HOUR_ESTIMATE.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(build(), indent=2, sort_keys=True))
