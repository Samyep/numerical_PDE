"""Turn timing-only calibration records into the immutable formal grid."""

from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess
from typing import Any

from invariant_region_mlp.experiments.whole_function_compute_frontier.common import (
    file_sha256,
    write_json_atomic,
)

from .engine import trajectory_path
from .protocol import (
    CALIBRATION_ANCHORS,
    CALIBRATION_DIR,
    CALIBRATION_MANIFEST,
    CANDIDATE_LADDERS,
    EXPECTED_INPUTS,
    FORMAL_PROJECTED_CAP_SECONDS,
    METHODS,
    PARAMETER_NAMES,
    PDES,
    PREREGISTRATION_COMMIT,
    RESULTS_DIR,
    SOURCE_ROOT,
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


def _monotone_anchors(anchors: list[tuple[int, float]]) -> list[tuple[int, float]]:
    """Sort anchors and replace downward timing noise by a cumulative maximum."""

    ordered = sorted((int(value), float(seconds)) for value, seconds in anchors)
    result: list[tuple[int, float]] = []
    running = 0.0
    for value, seconds in ordered:
        if value <= 0 or seconds <= 0:
            raise ValueError(f"invalid timing anchor: {(value, seconds)}")
        running = max(running, seconds)
        result.append((value, running))
    return result


def _log_slope(left: tuple[int, float], right: tuple[int, float]) -> float:
    x0, y0 = left
    x1, y1 = right
    if x1 <= x0:
        raise ValueError(f"timing parameters must be distinct: {left}, {right}")
    return max(0.0, math.log(y1 / y0) / math.log(x1 / x0))


def _predict_piecewise_log_linear(
    anchors: list[tuple[int, float]], value: int
) -> float:
    """Interpolate/extrapolate in log(parameter)--log(time) coordinates."""

    points = _monotone_anchors(anchors)
    if value <= 0:
        raise ValueError(f"parameter must be positive: {value}")
    if len(points) == 1:
        # With one timing anchor, work-proportional scaling is the unique
        # preregistered-information-only extrapolation used for C0.
        p0, t0 = points[0]
        return t0 * value / p0

    if value <= points[0][0]:
        left, right = points[0], points[1]
    elif value >= points[-1][0]:
        left, right = points[-2], points[-1]
    else:
        left, right = points[0], points[1]
        for candidate_left, candidate_right in zip(points, points[1:]):
            if candidate_left[0] <= value <= candidate_right[0]:
                left, right = candidate_left, candidate_right
                break
    slope = _log_slope(left, right)
    return left[1] * math.exp(slope * math.log(value / left[0]))


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
) -> tuple[
    list[int], dict[int, float], list[tuple[int, float]], int | None, str | None
]:
    monotone = _monotone_anchors(anchors)
    candidates = sorted({
        *CANDIDATE_LADDERS[(method, sweep)],
        *(value for value, _seconds in monotone),
    })
    predictions = {
        value: _predict_piecewise_log_linear(monotone, value)
        for value in candidates
    }
    # If the entire declared ladder remains short, time exactly one doubled
    # value.  Rebuilding after that observation implements the preregistered
    # repeated-doubling rule without extrapolating directly to an impractical
    # many-orders-of-magnitude parameter.
    if predictions[candidates[-1]] < 1500.0:
        inserted = candidates[-1] * 2
        candidates.append(inserted)
        predictions[inserted] = _predict_piecewise_log_linear(monotone, inserted)
        return candidates, predictions, monotone, inserted, "ladder_below_25_double"

    valid_end = [value for value in candidates if 1500.0 <= predictions[value] <= 2100.0]
    inserted: int | None = None
    supplemental_reason: str | None = None
    if not valid_end:
        below = [value for value in candidates if predictions[value] < 1500.0]
        above = [value for value in candidates if predictions[value] > 2100.0]
        if not below or not above:
            raise RuntimeError(
                f"cannot bracket the 25--35 minute window for {(method, sweep)}"
            )
        lower = max(below)
        upper = min(above)
        inserted = _round_insert(method, sweep, math.sqrt(lower * upper))
        if inserted <= lower:
            inserted = lower + 1
        if inserted >= upper:
            inserted = upper - 1
        if inserted not in candidates:
            candidates.append(inserted)
            candidates.sort()
        predictions[inserted] = _predict_piecewise_log_linear(monotone, inserted)
        supplemental_reason = "geometric_intermediate_for_window_jump"
    return candidates, predictions, monotone, inserted, supplemental_reason


def _projected_seconds(lines: list[dict[str, Any]]) -> float:
    return sum(
        3.0 * (
            max(float(v) for v in line["predicted_seconds"].values())
            if line["sweep"] == "C0"
            else sum(float(v) for v in line["predicted_seconds"].values())
        )
        for line in lines
    )


def _nearest_distinct(
    values: list[int], predictions: dict[int, float], targets: tuple[float, ...]
) -> list[int]:
    selected: list[int] = []
    for target in targets:
        ordered = sorted(values, key=lambda value: (abs(predictions[value] - target), value))
        choice = next((value for value in ordered if value not in selected), None)
        if choice is not None:
            selected.append(choice)
    return sorted(selected)


def _supplemental_spec(
    pde: str, method: str, sweep: str, value: int
) -> RunSpec:
    settings = settings_for(method, sweep, value)
    if sweep == "C0":
        # C0's varied parameter is the work amount itself.  A supplemental
        # endpoint therefore has to execute that exact work amount; a short
        # surrogate would merely repeat the original anchor.
        shortened_work = int(value)
        scale = 1.0
    else:
        shortened_work = int(CALIBRATION_ANCHORS[(method, sweep)][0][1])
        scale = float(CALIBRATION_ANCHORS[(method, sweep)][0][2])
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
    preflight_path = RESULTS_DIR / "preflight_manifest.json"
    if not preflight_path.is_file():
        raise RuntimeError("missing timing preflight manifest")
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
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
        (
            candidates, predictions, monotone_anchors, inserted,
            supplemental_reason,
        ) = _predicted_candidates(method, sweep, anchors)
        ends = [v for v in candidates if 1500.0 <= predictions[v] <= 2100.0]
        # A geometric intermediate must be measured before its timing can be
        # trusted.  Its provisional prediction need not itself land inside the
        # window when the bracketing jump is very large.
        end = (
            min(ends, key=lambda value: abs(predictions[value] - 1800.0))
            if ends else inserted
        )
        if end is None:
            raise RuntimeError(f"no calibrated endpoint for {key}")
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
            "time_model": "monotone_piecewise_log_linear",
            "monotone_scaled_anchors": [
                {"parameter_value": value, "scaled_formal_seconds": seconds}
                for value, seconds in monotone_anchors
            ],
            "supplemental_timing_candidate": inserted,
            "supplemental_timing_reason": supplemental_reason,
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
    missing_supplemental = []
    for spec in supplemental:
        path = trajectory_path(spec)
        if not path.is_file():
            missing_supplemental.append(spec)
            continue
        try:
            status = json.loads(path.read_text(encoding="utf-8")).get("status")
        except (OSError, json.JSONDecodeError):
            status = None
        if status not in {"TRAINED_NO_TRUTH", "COMPLETE"}:
            missing_supplemental.append(spec)
    if missing_supplemental:
        registry: dict[str, dict[str, Any]] = {}
        if supplemental_path.is_file():
            for payload in json.loads(supplemental_path.read_text(encoding="utf-8")):
                registry[str(payload["run_id"])] = payload
        for spec in missing_supplemental:
            registry[spec.run_id] = spec.payload()
        write_json_atomic(
            supplemental_path, [registry[key] for key in sorted(registry)]
        )
        return {
            "status": "SUPPLEMENTAL_TIMING_REQUIRED",
            "supplemental_specs": [spec.payload() for spec in missing_supplemental],
            "supplemental_registry_size": len(registry),
            "path": str(supplemental_path),
        }

    lines = [provisional[key] for key in sorted(provisional)]
    projected = _projected_seconds(lines)

    if projected > FORMAL_PROJECTED_CAP_SECONDS:
        original_values: dict[tuple[str, str, str], list[int]] = {}
        original_predictions: dict[tuple[str, str, str], dict[int, float]] = {}
        for line in lines:
            if line["sweep"] == "C0" or len(line["parameter_values"]) <= 3:
                continue
            key = (line["pde"], line["method"], line["sweep"])
            values = list(line["parameter_values"])
            original_values[key] = values
            predictions = {int(k): float(v) for k, v in line["predicted_seconds"].items()}
            original_predictions[key] = predictions
            keep = _nearest_distinct(values, predictions, (300.0, 900.0, 1800.0))
            line["parameter_values"] = keep
            line["predicted_seconds"] = {str(v): predictions[v] for v in keep}
        projected = _projected_seconds(lines)

        # Add the candidates nearest 10 and 20 minutes in deterministic
        # lexical round-robin order while respecting the frozen cap.
        for target in (600.0, 1200.0):
            for line in sorted(lines, key=lambda row: (row["pde"], row["method"], row["sweep"])):
                key = (line["pde"], line["method"], line["sweep"])
                if line["sweep"] == "C0" or key not in original_values:
                    continue
                all_values = original_values[key]
                all_predictions = original_predictions[key]
                remaining = [v for v in all_values if v not in line["parameter_values"]]
                if not remaining:
                    continue
                value = min(remaining, key=lambda v: (abs(all_predictions[v] - target), v))
                added = 3.0 * all_predictions[value]
                if projected + added <= FORMAL_PROJECTED_CAP_SECONDS:
                    line["parameter_values"] = sorted(line["parameter_values"] + [value])
                    line["predicted_seconds"][str(value)] = all_predictions[value]
                    line["predicted_seconds"] = {
                        str(v): line["predicted_seconds"][str(v)]
                        for v in line["parameter_values"]
                    }
                    projected += added
    if projected > FORMAL_PROJECTED_CAP_SECONDS:
        raise RuntimeError(
            f"minimum three-point design projects to {projected / 3600:.2f} GPU-hours"
        )

    payload = {
        "status": "FROZEN_BEFORE_FORMAL_TRAINING",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "preregistration_commit": PREREGISTRATION_COMMIT,
        "implementation_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True
        ).strip(),
        "timing_preflight_revision": preflight.get("checked_out_revision"),
        "timing_preflight_manifest_sha256": file_sha256(preflight_path),
        "calibration_module_sha256": file_sha256(Path(__file__)),
        "frozen_input_hashes": EXPECTED_INPUTS,
        "preserved_superseded_timing_artifacts": [
            {"path": str(path), "sha256": file_sha256(path)}
            for path in sorted(
                (CALIBRATION_DIR / "archived_pre_final_chunk").glob(
                    "*/raw/trajectories/*.json"
                )
            )
        ],
        "preserved_calibration_failure_artifacts": [
            {"path": str(path), "sha256": file_sha256(path)}
            for path in sorted([
                *CALIBRATION_DIR.glob("runs/*/failure_attempts/*.json"),
                *CALIBRATION_DIR.glob(
                    "archived_pre_final_chunk/*/failure_attempts/*.json"
                ),
            ])
        ],
        "selection_used_accuracy": False,
        "target_seconds": list(TARGET_SECONDS),
        "projected_formal_gpu_seconds": projected,
        "projected_formal_gpu_hours": projected / 3600.0,
        "selected_lines": lines,
    }
    write_json_atomic(CALIBRATION_MANIFEST, payload)
    return payload


__all__ = ["build_manifest"]
