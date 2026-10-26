"""Turn timing-only calibration records into the immutable formal grid."""

from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from typing import Any

from invariant_region_mlp.experiments.whole_function_compute_frontier.common import (
    file_sha256,
    write_json_atomic,
)

from .engine import trajectory_path
from .protocol import (
    CALIBRATION_ANCHORS,
    CALIBRATION_MANIFEST,
    CANDIDATE_LADDERS,
    FORMAL_PROJECTED_CAP_SECONDS,
    METHODS,
    PARAMETER_NAMES,
    PDES,
    RESULTS_DIR,
    SWEEPS,
    RunSpec,
    calibration_specs,
    settings_for,
)


TARGET_SECONDS = (300.0, 600.0, 900.0, 1200.0, 1800.0)


def _completed_specs() -> list[tuple[RunSpec, dict[str, Any]]]:
    specs = calibration_specs()
    supplemental = RESULTS_DIR / "calibration_supplemental_specs.json"
    if supplemental.is_file():
        for payload in json.loads(supplemental.read_text(encoding="utf-8")):
            payload.pop("run_id", None)
            payload["checkpoint_values"] = tuple(payload["checkpoint_values"])
            specs.append(RunSpec(**payload))
    completed: list[tuple[RunSpec, dict[str, Any]]] = []
    for spec in specs:
        path = trajectory_path(spec)
        if not path.is_file():
            continue
        record = json.loads(path.read_text(encoding="utf-8"))
        if record.get("status") not in {"TRAINED_NO_TRUTH", "COMPLETE"}:
            continue
        completed.append((spec, record))
    return completed


def _fit(anchors: list[tuple[int, float]]) -> tuple[float, float]:
    anchors = sorted(anchors)
    if len(anchors) == 1:
        p, seconds = anchors[0]
        return 0.0, seconds / p
    p0, t0 = anchors[0]
    p1, t1 = anchors[-1]
    raw_slope = (t1 - t0) / max(p1 - p0, 1)
    minimum_slope = 0.05 * min(t0 / p0, t1 / p1)
    slope = max(raw_slope, minimum_slope, 1e-9)
    intercept = max(0.0, min(t0 - slope * p0, t1 - slope * p1))
    return intercept, slope


def _round_insert(method: str, sweep: str, value: float) -> int:
    if method == "pinn_hte" and sweep == "C0":
        return max(1000, int(round(value / 1000.0)) * 1000)
    if method == "dpi" and sweep == "C0":
        return max(1, int(round(value / 10.0)) * 10)
    if method == "tangent_actor" and sweep == "C0":
        return max(1, int(round(value)))
    return max(1, int(round(value)))


def _predicted_candidates(
    method: str, sweep: str, anchors: list[tuple[int, float]]
) -> tuple[list[int], dict[int, float], tuple[float, float]]:
    intercept, slope = _fit(anchors)
    candidates = list(CANDIDATE_LADDERS[(method, sweep)])
    while intercept + slope * candidates[-1] < 1500.0:
        candidates.append(candidates[-1] * 2)
    predictions = {value: intercept + slope * value for value in candidates}
    valid_end = [value for value in candidates if 1500.0 <= predictions[value] <= 2100.0]
    if not valid_end:
        inserted = _round_insert(method, sweep, (1800.0 - intercept) / slope)
        if inserted not in candidates:
            candidates.append(inserted)
            candidates.sort()
        predictions[inserted] = intercept + slope * inserted
    return candidates, predictions, (intercept, slope)


def _supplemental_spec(
    pde: str, method: str, sweep: str, value: int
) -> RunSpec:
    shortened_work = int(CALIBRATION_ANCHORS[(method, sweep)][0][1])
    scale = float(CALIBRATION_ANCHORS[(method, sweep)][0][2])
    settings = settings_for(method, sweep, value)
    if method == "tangent_actor":
        settings["actor_updates"] = shortened_work
    elif method == "dpi":
        settings["picard_iterations"] = shortened_work
    else:
        settings["training_epochs"] = shortened_work
    return RunSpec(
        phase="calibration", pde=pde, method=method, sweep=sweep,
        parameter_name=PARAMETER_NAMES[(method, sweep)],
        parameter_value=value, seed=314159, settings=settings,
        checkpoint_values=(shortened_work,), calibration_scale=scale,
    )


def build_manifest() -> dict[str, Any]:
    completed = _completed_specs()
    completed_ids = {spec.run_id for spec, _record in completed}
    missing_initial = [
        spec.run_id for spec in calibration_specs()
        if spec.run_id not in completed_ids
    ]
    if missing_initial:
        raise RuntimeError(
            f"missing {len(missing_initial)} preregistered calibration anchors: "
            f"{missing_initial}"
        )
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for spec, record in completed:
        grouped[(spec.pde, spec.method, spec.sweep)].append({
            "spec": spec,
            "training_seconds": float(record["training_seconds"]),
            "trajectory": str(trajectory_path(spec)),
            "trajectory_sha256": file_sha256(trajectory_path(spec)),
            "peak_gpu_bytes": record.get("peak_gpu_bytes"),
        })

    missing = [
        (pde, method, sweep)
        for pde in PDES for method in METHODS for sweep in SWEEPS
        if not grouped[(pde, method, sweep)]
    ]
    if missing:
        raise RuntimeError(f"missing timing calibration lines: {missing}")

    supplemental: list[RunSpec] = []
    provisional: dict[tuple[str, str, str], dict[str, Any]] = {}
    for key, rows in grouped.items():
        pde, method, sweep = key
        anchors = sorted(
            (int(row["spec"].parameter_value),
             float(row["training_seconds"]) * float(row["spec"].calibration_scale))
            for row in rows
        )
        candidates, predictions, fit = _predicted_candidates(method, sweep, anchors)
        ends = [v for v in candidates if 1500.0 <= predictions[v] <= 2100.0]
        end = min(ends, key=lambda value: abs(predictions[value] - 1800.0))
        measured_values = {int(row["spec"].parameter_value) for row in rows}
        if end not in measured_values and end not in CANDIDATE_LADDERS[(method, sweep)]:
            supplemental.append(_supplemental_spec(pde, method, sweep, end))
        selected = sorted({
            min(candidates, key=lambda value: abs(predictions[value] - target))
            for target in TARGET_SECONDS
            if predictions[min(candidates, key=lambda value: abs(predictions[value] - target))]
            <= 2100.0
        })
        selected = [value for value in selected if value <= end]
        if end not in selected:
            selected.append(end)
        while len(selected) < 3:
            available = [v for v in candidates if v <= end and v not in selected]
            if not available:
                break
            selected.append(max(available))
            selected.sort()
        provisional[key] = {
            "pde": pde, "method": method, "sweep": sweep,
            "parameter_name": PARAMETER_NAMES[(method, sweep)],
            "parameter_values": selected,
            "predicted_seconds": {str(v): predictions[v] for v in selected},
            "fit_intercept_seconds": fit[0], "fit_slope_seconds_per_unit": fit[1],
            "calibration_anchors": [
                {
                    "parameter_value": int(row["spec"].parameter_value),
                    "measured_short_run_seconds": float(row["training_seconds"]),
                    "formal_scale": float(row["spec"].calibration_scale),
                    "scaled_formal_seconds": (
                        float(row["training_seconds"])
                        * float(row["spec"].calibration_scale)
                    ),
                    "trajectory": row["trajectory"],
                    "trajectory_sha256": row["trajectory_sha256"],
                    "peak_gpu_bytes": row["peak_gpu_bytes"],
                }
                for row in rows
            ],
        }

    supplemental_path = RESULTS_DIR / "calibration_supplemental_specs.json"
    missing_supplemental = [
        spec for spec in supplemental if not trajectory_path(spec).is_file()
    ]
    if missing_supplemental:
        write_json_atomic(
            supplemental_path, [spec.payload() for spec in missing_supplemental]
        )
        return {
            "status": "SUPPLEMENTAL_TIMING_REQUIRED",
            "supplemental_specs": [spec.payload() for spec in missing_supplemental],
            "path": str(supplemental_path),
        }

    lines = [provisional[key] for key in sorted(provisional)]
    projected = 0.0
    for line in lines:
        predictions = {int(k): float(v) for k, v in line["predicted_seconds"].items()}
        if line["sweep"] == "C0":
            projected += 3.0 * predictions[max(line["parameter_values"])]
        else:
            projected += 3.0 * sum(predictions[v] for v in line["parameter_values"])

    if projected > FORMAL_PROJECTED_CAP_SECONDS:
        for line in lines:
            if line["sweep"] == "C0" or len(line["parameter_values"]) <= 3:
                continue
            values = line["parameter_values"]
            predictions = {int(k): float(v) for k, v in line["predicted_seconds"].items()}
            keep = sorted({
                min(values, key=lambda value: abs(predictions[value] - target))
                for target in (300.0, 900.0, 1800.0)
            })
            line["parameter_values"] = keep
            line["predicted_seconds"] = {str(v): predictions[v] for v in keep}
        projected = sum(
            3.0 * (
                max(float(v) for v in line["predicted_seconds"].values())
                if line["sweep"] == "C0"
                else sum(float(v) for v in line["predicted_seconds"].values())
            )
            for line in lines
        )
    if projected > FORMAL_PROJECTED_CAP_SECONDS:
        raise RuntimeError(
            f"minimum three-point design projects to {projected / 3600:.2f} GPU-hours"
        )

    payload = {
        "status": "FROZEN_BEFORE_FORMAL_TRAINING",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection_used_accuracy": False,
        "target_seconds": list(TARGET_SECONDS),
        "projected_formal_gpu_seconds": projected,
        "projected_formal_gpu_hours": projected / 3600.0,
        "selected_lines": lines,
    }
    write_json_atomic(CALIBRATION_MANIFEST, payload)
    return payload


__all__ = ["build_manifest"]
