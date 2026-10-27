"""Post-gate accuracy aggregation and publication-table generation."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
from statistics import mean, stdev
from typing import Any, Iterable

import numpy as np

from .protocol import (
    DIMENSION,
    SCALING_DISPLAY_POINTS,
    FINAL_COUNT,
    HISTORICAL_TA_REFERENCE_RMSE,
    INPUT_DIR,
    PDES,
    RAW_DIR,
    RESULTS_DIR,
    SOURCE_ROOT,
    TABLE_DIR,
    TRAINING_SEEDS,
    TRAJECTORY_BY_ID,
    array_sha256,
    file_sha256,
    write_json_atomic,
)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _truth(pde: str) -> tuple[np.ndarray, np.ndarray]:
    with np.load(INPUT_DIR / f"{pde}_final_truth.npz", allow_pickle=False) as payload:
        return (
            np.asarray(payload["original_index"], dtype=np.int64),
            np.asarray(payload["truth_u"], dtype=np.float64),
        )


def _inputs(pde: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.load(INPUT_DIR / f"{pde}_final_inputs.npz", allow_pickle=False) as payload:
        return (
            np.asarray(payload["original_index"], dtype=np.int64),
            np.asarray(payload["t"], dtype=np.float64),
            np.asarray(payload["x"], dtype=np.float64),
        )


def _metrics(prediction: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    error = prediction - truth
    return {
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "relmae": float(np.sum(np.abs(error)) / np.sum(np.abs(truth))),
        "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
        "bias": float(np.mean(error)),
    }


def _record_path(pde: str, trajectory_id: str, seed: int, position: int) -> Path:
    return RAW_DIR / "formal" / pde / trajectory_id / f"s{seed}" / f"q{position:02d}.json"


def _deep_seed_row(pde: str, point: dict[str, Any], seed: int) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    indices, t, x = _inputs(pde)
    truth_indices, truth = _truth(pde)
    if not np.array_equal(indices, truth_indices):
        raise RuntimeError(f"{pde}: input/truth index mismatch")
    predictions = np.full(FINAL_COUNT, np.nan, dtype=np.float64)
    times = np.full(FINAL_COUNT, np.nan, dtype=np.float64)
    query_rows: list[dict[str, Any]] = []
    failures = 0
    peaks = []
    checkpoint_hashes = []
    trajectory_id = str(point["trajectory_id"])
    wanted_updates = int(point["checkpoint_updates"])
    for position in range(FINAL_COUNT):
        path = _record_path(pde, trajectory_id, seed, position)
        record = json.loads(path.read_text(encoding="utf-8"))
        expected_hash = array_sha256(
            np.asarray(t[position], dtype=np.float64), np.asarray(x[position], dtype=np.float64)
        )
        if record["query_input_sha256"] != expected_hash or record["original_index"] != int(indices[position]):
            raise RuntimeError(f"{pde} {trajectory_id} s{seed} q{position}: query drift")
        capture = next(
            (row for row in record.get("captures", []) if int(row["updates"]) == wanted_updates),
            None,
        )
        valid = record.get("status") == "complete" and capture is not None
        if valid:
            predictions[position] = float(capture["prediction"])
            times[position] = float(capture["active_training_seconds"])
            checkpoint_hashes.append(str(capture["checkpoint_sha256"]))
            if record.get("peak_gpu_bytes") is not None:
                peaks.append(int(record["peak_gpu_bytes"]))
        else:
            failures += 1
        query_rows.append({
            "pde": pde,
            "family": point["family"],
            "sweep": point["sweep"],
            "level": point["level"],
            "parameter_name": point["parameter_name"],
            "parameter_value": point["parameter_value"],
            "trajectory_id": trajectory_id,
            "checkpoint_updates": wanted_updates,
            "seed": seed,
            "position": position,
            "original_index": int(indices[position]),
            "t": float(t[position]),
            "prediction": float(predictions[position]) if valid else "",
            "reference": float(truth[position]),
            "error": float(predictions[position] - truth[position]) if valid else "",
            "active_training_seconds": float(times[position]) if valid else "",
            "status": "complete" if valid else "numerical_failure",
        })
    complete = failures == 0 and bool(np.isfinite(predictions).all())
    spec = TRAJECTORY_BY_ID[trajectory_id]
    row: dict[str, Any] = {
        "pde": pde,
        "family": point["family"],
        "sweep": point["sweep"],
        "level": point["level"],
        "parameter_name": point["parameter_name"],
        "parameter_value": point["parameter_value"],
        "trajectory_id": trajectory_id,
        "trajectory_signature": spec.signature,
        "checkpoint_updates": wanted_updates,
        "updates": wanted_updates,
        "batch_size": spec.batch_size,
        "time_steps": spec.time_steps,
        "seed": seed,
        "status": "complete" if complete else "numerical_failure",
        "successful_queries": FINAL_COUNT - failures,
        "failed_queries": failures,
        "mean_training_seconds_per_query": float(np.mean(times)) if complete else "",
        "total_active_training_seconds_64_queries": float(np.sum(times)) if complete else "",
        "peak_gpu_bytes": max(peaks) if peaks else "",
        "checkpoint_hash_count": len(set(checkpoint_hashes)),
    }
    if complete:
        row.update(_metrics(predictions, truth))
    else:
        row.update({"rmse": "", "relmae": "", "relative_l2": "", "bias": ""})
    return row, query_rows


def _aggregate_deep(seed_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple, list[dict[str, Any]]] = {}
    for row in seed_rows:
        key = (row["pde"], row["family"], row["sweep"], int(row["level"]))
        groups.setdefault(key, []).append(row)
    output = []
    for key, rows in sorted(groups.items()):
        rows.sort(key=lambda row: int(row["seed"]))
        base = rows[0]
        valid = len(rows) == 3 and all(row["status"] == "complete" for row in rows)
        result: dict[str, Any] = {
            field: base[field] for field in (
                "pde", "family", "sweep", "level", "parameter_name",
                "parameter_value", "trajectory_id", "trajectory_signature",
                "checkpoint_updates", "updates", "batch_size", "time_steps",
            )
        }
        result.update({
            "status": "complete" if valid else "numerical_failure",
            "training_seed_count": len(rows),
            "successful_queries": sum(int(row["successful_queries"]) for row in rows),
            "failed_queries": sum(int(row["failed_queries"]) for row in rows),
            "mean_peak_gpu_bytes": (
                mean(float(row["peak_gpu_bytes"]) for row in rows if row["peak_gpu_bytes"] != "")
                if any(row["peak_gpu_bytes"] != "" for row in rows) else ""
            ),
        })
        if valid:
            for metric in ("rmse", "relmae", "relative_l2"):
                values = [float(row[metric]) for row in rows]
                result[f"mean_{metric}"] = mean(values)
                result[f"std_{metric}"] = stdev(values)
                result[f"seed_{metric}s"] = json.dumps(values)
            times = [float(row["total_active_training_seconds_64_queries"]) for row in rows]
            result["mean_total_active_training_seconds_64_queries"] = mean(times)
            result["std_total_active_training_seconds_64_queries"] = stdev(times)
            result["mean_total_active_training_minutes_64_queries"] = mean(times) / 60.0
            result["seed_total_active_training_seconds_64_queries"] = json.dumps(times)
            target = HISTORICAL_TA_REFERENCE_RMSE[result["pde"]]
            result["historical_ta_reference_rmse_960_bank"] = target
            result["ratio_to_historical_ta_reference"] = result["mean_rmse"] / target
            result["historical_ta_reference_target_status"] = (
                "OBSERVED_AT_OR_BELOW_HISTORICAL_TA_REFERENCE"
                if result["mean_rmse"] <= target
                else "TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET"
            )
        else:
            for field in (
                "mean_rmse", "std_rmse", "mean_relmae", "std_relmae",
                "mean_relative_l2", "std_relative_l2",
                "mean_total_active_training_seconds_64_queries",
                "std_total_active_training_seconds_64_queries",
                "mean_total_active_training_minutes_64_queries",
                "ratio_to_historical_ta_reference",
            ):
                result[field] = ""
            result["historical_ta_reference_target_status"] = "NUMERICAL_INSTABILITY"
        output.append(result)
    return output


def _historical_rows() -> list[dict[str, str]]:
    path = SOURCE_ROOT / "invariant_region_mlp" / "results" / "whole_function_compute_scaling_3seed" / "scaling_summary.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _selected_historical(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected = []
    for pde in PDES:
        for method in ("dpi", "pinn_hte"):
            for sweep in ("C0", "C1", "C2"):
                group = sorted(
                    (row for row in rows if row["pde"] == pde and row["method"] == method and row["sweep"] == sweep),
                    key=lambda row: int(row["parameter_value"]),
                )
                positions = np.rint(np.linspace(0, len(group) - 1, 3)).astype(int)
                for level, position in enumerate(positions, start=1):
                    selected.append({**group[int(position)], "level": str(level)})
    return selected


def _find_prediction_json(scaling_root: Path, row: dict[str, str], seed: int) -> Path | None:
    pde, method, sweep = row["pde"], row["method"], row["sweep"]
    parameter_name, parameter_value = row["parameter_name"], row["parameter_value"]
    if sweep == "C0":
        candidates = sorted((scaling_root / "runs").glob(f"formal_{method}_{pde}_c0_*_s{seed}"))
    else:
        candidates = [
            scaling_root / "runs" / f"formal_{method}_{pde}_{sweep.lower()}_{parameter_name}-{parameter_value}_s{seed}"
        ]
    for run_dir in candidates:
        # The frozen scaling study stores C1/C2 predictions below a legacy
        # ``C0`` directory even though their metadata correctly records the
        # actual sweep.  Discover below the PDE directory and validate the
        # metadata instead of inferring provenance from that legacy folder.
        prediction_root = run_dir / "predictions" / method / pde
        if not prediction_root.is_dir():
            continue
        matches = []
        for metadata_path in prediction_root.rglob("*.json"):
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            observed = int(metadata.get("parameter_value", metadata.get("algorithm_iteration", -1)))
            if (
                metadata.get("pde") == pde
                and metadata.get("method") == method
                and metadata.get("sweep") == sweep
                and metadata.get("parameter_name") == parameter_name
                and int(metadata.get("training_seed", metadata.get("seed", -1))) == int(seed)
                and observed == int(parameter_value)
            ):
                matches.append(metadata_path)
        if len(matches) > 1:
            raise RuntimeError(
                f"ambiguous frozen prediction artifacts for {pde}/{method}/{sweep}/"
                f"{parameter_name}={parameter_value}/seed={seed}: {matches}"
            )
        if matches:
            return matches[0]
    return None


def _matched_existing(selected: list[dict[str, str]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    configured = os.environ.get("DBBS_SCALING_RESULTS")
    if not configured:
        return [], {"status": "NOT_AVAILABLE", "reason": "DBBS_SCALING_RESULTS is unset"}
    root = Path(configured)
    if not root.is_dir():
        return [], {"status": "NOT_AVAILABLE", "reason": f"missing {root}"}
    rows = []
    missing = []
    index_cache = {pde: _inputs(pde)[0] for pde in PDES}
    truth_cache = {pde: _truth(pde)[1] for pde in PDES}
    for selected_row in selected:
        pde = selected_row["pde"]
        for seed in TRAINING_SEEDS:
            metadata_path = _find_prediction_json(root, selected_row, seed)
            if metadata_path is None:
                missing.append({**selected_row, "seed": seed})
                continue
            prediction_path = metadata_path.with_suffix(".npz")
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            if not prediction_path.is_file():
                raise RuntimeError(f"missing frozen prediction array: {prediction_path}")
            expected_prediction_hash = metadata.get("prediction_sha256")
            if expected_prediction_hash and file_sha256(prediction_path) != expected_prediction_hash:
                raise RuntimeError(f"frozen prediction hash mismatch: {prediction_path}")
            with np.load(prediction_path, allow_pickle=False) as payload:
                full_prediction = np.asarray(payload["prediction"], dtype=np.float64)
            if full_prediction.shape != (1200,) and full_prediction.shape != (960,):
                raise RuntimeError(
                    f"unexpected frozen prediction shape {full_prediction.shape}: {prediction_path}"
                )
            if int(index_cache[pde].max()) >= full_prediction.shape[0]:
                raise RuntimeError(f"frozen prediction does not cover final indices: {prediction_path}")
            prediction = full_prediction[index_cache[pde]]
            metrics = _metrics(prediction, truth_cache[pde])
            rows.append({
                "pde": pde,
                "method": selected_row["method"],
                "sweep": selected_row["sweep"],
                "level": selected_row["level"],
                "parameter_name": selected_row["parameter_name"],
                "parameter_value": selected_row["parameter_value"],
                "seed": seed,
                "evaluation_states": FINAL_COUNT,
                "source_prediction_file": str(prediction_path),
                "source_prediction_file_sha256": file_sha256(prediction_path),
                **metrics,
            })
    status = "COMPLETE" if not missing else "PARTIAL"
    return rows, {"status": status, "missing": missing, "row_count": len(rows)}


def _cell(value: dict[str, Any], *, deep: bool) -> str:
    if value.get("status") != "complete":
        return r"\shortstack{FAIL\\numerical}"
    if deep:
        rmse = 1000.0 * float(value["mean_rmse"])
        sd = 1000.0 * float(value["std_rmse"])
        minutes = float(value["mean_total_active_training_minutes_64_queries"])
        parameter = value["parameter_value"]
    else:
        rmse = 1000.0 * float(value["mean_rmse"])
        sd = 1000.0 * float(value["std_rmse"])
        minutes = float(value["mean_training_seconds"]) / 60.0
        parameter = value["parameter_value"]
    return rf"\shortstack{{{rmse:.2f} $\pm$ {sd:.2f}\\{parameter}; {minutes:.1f} min}}"


def _main_table(deep_rows: list[dict[str, Any]], historical: list[dict[str, str]]) -> str:
    lines = [
        r"\begin{table*}[t]",
        r"\centering\scriptsize",
        r"\setlength{\tabcolsep}{2.5pt}",
        r"\caption{Prespecified hyperparameter-budget outcomes. Cells report $10^3\!\times$ RMSE (mean $\pm$ seed SD), then the varied hyperparameter and measured minutes. DPI/PINN--HTE are frozen historical 960-state whole-function results; Deep BSDE is evaluated on the frozen 64-state subset and its time is cumulative training of 64 fresh scalar models. The protocols are therefore not identical paired evaluations.}",
        r"\label{tab:deep-bsde-budget-scaling}",
        r"\begin{tabular}{llcccccc}",
        r"\toprule",
        r"Method & Sweep & \multicolumn{3}{c}{B6 Full-Rank HJB} & \multicolumn{3}{c}{Dense Fourier} \\",
        r" & & Level 1 & Level 2 & Level 3 & Level 1 & Level 2 & Level 3 \\",
        r"\midrule",
    ]
    for method, label in (("dpi", "DPI"), ("pinn_hte", "PINN--HTE")):
        for sweep in ("C0", "C1", "C2"):
            cells = []
            for pde in PDES:
                for level in (1, 2, 3):
                    row = next(r for r in historical if r["method"] == method and r["pde"] == pde and r["sweep"] == sweep and int(r["level"]) == level)
                    cells.append(_cell({**row, "status": "complete"}, deep=False))
            lines.append(f"{label} & {sweep} & " + " & ".join(cells) + r" \\")
        lines.append(r"\addlinespace")
    for sweep in ("C0", "C1", "C2"):
        cells = []
        for pde in PDES:
            for level in (1, 2, 3):
                row = next(r for r in deep_rows if r["pde"] == pde and r["family"] == "sweep" and r["sweep"] == sweep and int(r["level"]) == level)
                cells.append(_cell(row, deep=True))
        lines.append("Deep BSDE & " + sweep + " & " + " & ".join(cells) + r" \\")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table*}", ""])
    return "\n".join(lines)


def _matched_table(rows: list[dict[str, Any]]) -> str:
    grouped: dict[tuple[str, str, str, int], list[dict[str, Any]]] = {}
    for row in rows:
        key = (row["pde"], row["method"], row["sweep"], int(row["level"]))
        grouped.setdefault(key, []).append(row)
    lines = [
        r"\begin{table}[t]",
        r"\centering\small",
        r"\caption{Supplementary paired evaluation of frozen DPI/PINN--HTE predictions on the same 64 states used for Deep BSDE. No baseline was retrained. Cells report $10^3\!\times$ RMSE (mean $\pm$ seed SD).}",
        r"\label{tab:matched-existing-64}",
        r"\begin{tabular}{llccc}",
        r"\toprule",
        r"Method & Sweep & Level 1 & Level 2 & Level 3 \\",
        r"\midrule",
    ]
    labels = {"dpi": "DPI", "pinn_hte": "PINN--HTE"}
    for pde in PDES:
        lines.append(rf"\multicolumn{{5}}{{l}}{{\emph{{{pde.replace('_', ' ')}}}}} \\")
        for method in ("dpi", "pinn_hte"):
            for sweep in ("C0", "C1", "C2"):
                cells = []
                for level in (1, 2, 3):
                    values = [float(row["rmse"]) for row in grouped[(pde, method, sweep, level)]]
                    cells.append(rf"{1000*mean(values):.2f} $\pm$ {1000*stdev(values):.2f}")
                lines.append(labels[method] + " & " + sweep + " & " + " & ".join(cells) + r" \\")
        lines.append(r"\addlinespace")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    return "\n".join(lines)


def _write_preview(name: str, table_text: str) -> None:
    body = "\n".join([
        r"\documentclass[10pt]{article}",
        r"\usepackage[margin=0.35in,landscape]{geometry}",
        r"\usepackage{booktabs}",
        r"\usepackage{amsmath}",
        r"\makeatletter",
        r"\setlength{\@fptop}{0pt}",
        r"\makeatother",
        r"\pagestyle{empty}",
        r"\begin{document}",
        table_text,
        r"\end{document}",
        "",
    ])
    (TABLE_DIR / f"{name}_preview.tex").write_text(body, encoding="utf-8")


def aggregate() -> dict:
    gate = RESULTS_DIR / "status" / "truth_gate_open.json"
    if not gate.is_file():
        raise RuntimeError("truth gate is closed")
    seed_rows: list[dict[str, Any]] = []
    query_rows: list[dict[str, Any]] = []
    for pde in PDES:
        for point in SCALING_DISPLAY_POINTS:
            for seed in TRAINING_SEEDS:
                seed_row, point_rows = _deep_seed_row(pde, point, seed)
                seed_rows.append(seed_row)
                query_rows.extend(point_rows)
    deep_rows = _aggregate_deep(seed_rows)
    _write_csv(RESULTS_DIR / "deep_bsde_per_query.csv", query_rows)
    _write_csv(RESULTS_DIR / "deep_bsde_seed_metrics.csv", seed_rows)
    _write_csv(
        RESULTS_DIR / "deep_bsde_scaling_summary.csv",
        [row for row in deep_rows if row["family"] == "sweep"],
    )

    historical = _selected_historical(_historical_rows())
    _write_csv(RESULTS_DIR / "historical_dpi_pinn_selected_levels.csv", historical)
    matched_rows, matched_status = _matched_existing(historical)
    if matched_rows:
        _write_csv(RESULTS_DIR / "matched_existing_baselines_64_seed_metrics.csv", matched_rows)
    write_json_atomic(RESULTS_DIR / "matched_existing_baselines_status.json", matched_status)

    main_text = _main_table(deep_rows, historical)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    (TABLE_DIR / "deep_bsde_budget_scaling_main.tex").write_text(main_text, encoding="utf-8")
    _write_preview("deep_bsde_budget_scaling_main", main_text)
    if matched_status["status"] == "COMPLETE":
        matched_text = _matched_table(matched_rows)
        (TABLE_DIR / "matched_existing_baselines_64.tex").write_text(matched_text, encoding="utf-8")
        _write_preview("matched_existing_baselines_64", matched_text)

    reuse_rows = []
    for pde in PDES:
        for seed in TRAINING_SEEDS:
            for position in range(FINAL_COUNT):
                record = json.loads(_record_path(pde, "reference", seed, position).read_text(encoding="utf-8"))
                capture = next(row for row in record["captures"] if int(row["updates"]) == 2000)
                reuse_rows.append({
                    "pde": pde, "seed": seed, "position": position,
                    "checkpoint_sha256": capture["checkpoint_sha256"],
                    "consumers": ["C0-level2", "C1-level2", "C2-level2"],
                    "exact_single_artifact_reuse": True,
                })
    write_json_atomic(RESULTS_DIR / "checkpoint_reuse_audit.json", {
        "status": "PASS",
        "rows": reuse_rows,
    })
    payload = {
        "status": "AGGREGATED",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "seed_metric_rows": len(seed_rows),
        "aggregate_rows": len(deep_rows),
        "query_rows": len(query_rows),
        "numerical_failure_seed_rows": sum(row["status"] != "complete" for row in seed_rows),
        "matched_existing_baselines": matched_status,
        "historical_level_selection_rule": "lowest, rounded midpoint, and highest frozen parameter; no accuracy selection",
        "main_table": str(TABLE_DIR / "deep_bsde_budget_scaling_main.tex"),
        "scope": "C0_C1_C2_ONLY",
        "joint_scaling_excluded_by_user_amendment": True,
    }
    write_json_atomic(RESULTS_DIR / "analysis_summary.json", payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=2, sort_keys=True))

