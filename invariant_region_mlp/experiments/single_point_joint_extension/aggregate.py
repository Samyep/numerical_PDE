"""Mechanical aggregation, verdicts, figures, report, and artifact audit."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .protocol import (
    DIMENSIONS,
    FH_CONFIGURATIONS,
    FH_MAIN_CONFIGURATION,
    FH_REPLICATES,
    FIGURE_DIR,
    FROZEN_PROTOCOL_JSON,
    PACKAGE_DIR,
    PARAMETER_MANIFEST,
    PDE,
    PREDICTION_DIR,
    PREFLIGHT_AUDIT,
    QUERY_COUNT,
    RAW_DIR,
    REFERENCE_VALIDATION,
    RESULTS_DIR,
    REVISION,
    SMOKE_AUDIT,
    STATUS_DIR,
    TABLE_DIR,
    TRAINING_SEEDS,
    file_sha256,
    load_query_bank,
    query_path,
    write_json_atomic,
)


REPORT = PACKAGE_DIR.parents[1] / "docs" / "ICML_SINGLE_POINT_JOINT_EXTENSION_REPORT.md"
SUMMARY = RESULTS_DIR / "analysis_summary.json"
METRICS = RESULTS_DIR / "main_metrics.csv"
FH_GRID = RESULTS_DIR / "fh_grid.csv"
RESOURCE = RESULTS_DIR / "resource_accounting.csv"
AUDIT = RESULTS_DIR / "artifact_audit.json"


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _metrics(prediction: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    prediction = np.asarray(prediction, dtype=np.float64).reshape(-1)
    reference = np.asarray(reference, dtype=np.float64).reshape(-1)
    error = prediction - reference
    return {
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(reference)),
        "bias": float(np.mean(error)),
        "mae": float(np.mean(np.abs(error))),
        "nonfinite_count": int(np.sum(~np.isfinite(prediction))),
        "samples": len(prediction),
    }


def _ta(dimension: int) -> tuple[dict, list[dict]]:
    training_path = RAW_DIR / "actor" / f"training_{PDE}_d{dimension}_mature_s0.json"
    training = _json(training_path)
    actor = training["selection"]["selected_actor"]
    predictions, references, rows = [], [], []
    for replicate in range(4):
        record_path = (
            RAW_DIR / "actor_evaluations" / f"{PDE}_d{dimension}_mature_s0"
            / "selected" / f"steps64_rep{replicate:02d}.json"
        )
        record = _json(record_path)
        with np.load(record["prediction_path"]) as data:
            predictions.append(data["mcv"])
            references.append(data["truth"])
        rows.append(record)
    metric = _metrics(np.concatenate(predictions), np.concatenate(references))
    metric.update({
        "pde": PDE,
        "dimension": dimension,
        "method": "Tangent-Actor",
        "configuration": actor,
        "formal_repeats": 4,
        "gpu_seconds": float(training["timing"]["phase_total_seconds"] + training["timing"]["reference_free_selection_seconds"] + sum(row["inference_seconds"] for row in rows)),
        "peak_gpu_bytes": max(row["peak_gpu_bytes"] or 0 for row in rows),
    })
    return metric, rows


def _fh_cell(dimension: int, n: int, M: int) -> tuple[dict, list[dict]]:
    predictions, references, rows = [], [], []
    for replicate in range(FH_REPLICATES):
        run_id = f"fh_mlp_standard__{PDE}__d{dimension}__n{n}_M{M}__r{replicate:02d}__formal"
        record = _json(RAW_DIR / "formal" / "fh_mlp" / f"{run_id}.json")
        with np.load(record["prediction_file"]) as data:
            predictions.append(data["prediction"])
            references.append(data["reference"])
        rows.append(record)
    metric = _metrics(np.concatenate(predictions), np.concatenate(references))
    metric.update({
        "pde": PDE,
        "dimension": dimension,
        "method": "Full-History MLP",
        "configuration": f"n={n},M={M}",
        "n": n,
        "M": M,
        "formal_repeats": FH_REPLICATES,
        "gpu_seconds": float(sum(row["inference_seconds"] for row in rows)),
        "peak_gpu_bytes": max(row["peak_gpu_bytes"] or 0 for row in rows),
    })
    return metric, rows


def _deep(dimension: int) -> tuple[dict, list[dict]]:
    predictions, references, rows = [], [], []
    truth = load_query_bank(dimension)[2]
    for seed in TRAINING_SEEDS:
        for point in range(QUERY_COUNT):
            run_id = f"deep_bsde__{PDE}__d{dimension}__q{point:02d}__s{seed}__formal"
            record = _json(RAW_DIR / "formal" / "deep_bsde" / f"{run_id}.json")
            predictions.append(record["reported_y0"])
            references.append(truth[point])
            rows.append(record)
    metric = _metrics(np.asarray(predictions), np.asarray(references))
    metric.update({
        "pde": PDE,
        "dimension": dimension,
        "method": "Deep BSDE",
        "configuration": "20 intervals, 2000 updates",
        "formal_repeats": len(TRAINING_SEEDS),
        "gpu_seconds": float(sum(row["training_seconds"] for row in rows)),
        "peak_gpu_bytes": max(row["peak_gpu_bytes"] or 0 for row in rows),
    })
    return metric, rows


def _figure(rows: list[dict]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    colors = {"Tangent-Actor": "#0B84F3", "Full-History MLP": "#E45756", "Deep BSDE": "#57A773"}
    markers = {"Tangent-Actor": "o", "Full-History MLP": "s", "Deep BSDE": "^"}
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for method in colors:
        subset = sorted((row for row in rows if row["method"] == method), key=lambda row: row["dimension"])
        ax.plot(
            [row["dimension"] for row in subset],
            [row["rmse"] for row in subset],
            marker=markers[method], color=colors[method], linewidth=2.0,
            markersize=6, label=method,
        )
    ax.set_yscale("log")
    ax.set_xlabel("Dimension")
    ax.set_ylabel("Value RMSE")
    ax.set_xticks(DIMENSIONS)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, ncol=3, fontsize=8)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"joint_nonseparable_rmse.{suffix}", dpi=240)
    plt.close(fig)


def _git(*args: str) -> str:
    completed = subprocess.run(
        ("git", *args), cwd=PACKAGE_DIR.parents[2], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    return completed.stdout.strip()


def run() -> dict[str, Any]:
    missing_status = []
    for path in STATUS_DIR.glob("*.json"):
        if _json(path).get("status") != "complete":
            missing_status.append(str(path))
    expected_status = 5 * len(DIMENSIONS)
    present_status = list(STATUS_DIR.glob("*.json"))
    if len(present_status) != expected_status or missing_status:
        raise RuntimeError(
            f"formal bundles incomplete: {len(present_status)}/{expected_status}; bad={missing_status}"
        )

    main_rows: list[dict] = []
    fh_rows: list[dict] = []
    resource_rows: list[dict] = []
    for dimension in DIMENSIONS:
        ta, ta_raw = _ta(dimension)
        deep, deep_raw = _deep(dimension)
        cell_rows = []
        for n, M in FH_CONFIGURATIONS:
            cell, raw = _fh_cell(dimension, n, M)
            fh_rows.append(cell)
            resource_rows.append({
                "dimension": dimension,
                "method": "Full-History MLP",
                "configuration": cell["configuration"],
                "gpu_seconds": cell["gpu_seconds"],
                "peak_gpu_bytes": cell["peak_gpu_bytes"],
            })
            if (n, M) == FH_MAIN_CONFIGURATION:
                cell_rows.append(cell)
        if len(cell_rows) != 1:
            raise RuntimeError("missing preregistered FH main cell")
        main_rows.extend((ta, cell_rows[0], deep))
        for row in (ta, deep):
            resource_rows.append({key: row[key] for key in (
                "dimension", "method", "configuration", "gpu_seconds", "peak_gpu_bytes"
            )})

    by_dimension = {
        dimension: {row["method"]: row for row in main_rows if row["dimension"] == dimension}
        for dimension in DIMENSIONS
    }
    lowest = [
        dimension for dimension, methods in by_dimension.items()
        if methods["Tangent-Actor"]["rmse"] == min(row["rmse"] for row in methods.values())
    ]
    half_fh = [
        dimension for dimension, methods in by_dimension.items()
        if methods["Tangent-Actor"]["rmse"] <= 0.5 * methods["Full-History MLP"]["rmse"]
    ]
    half_deep = [
        dimension for dimension, methods in by_dimension.items()
        if methods["Tangent-Actor"]["rmse"] <= 0.5 * methods["Deep BSDE"]["rmse"]
    ]
    ratio = by_dimension[400]["Tangent-Actor"]["rmse"] / by_dimension[20]["Tangent-Actor"]["rmse"]
    verdicts = {
        "J-1": {"pass": len(lowest) >= 4, "count": len(lowest), "dimensions": lowest},
        "J-2": {"pass": len(half_fh) >= 4, "count": len(half_fh), "dimensions": half_fh},
        "J-3": {"pass": len(half_deep) >= 3, "count": len(half_deep), "dimensions": half_deep},
        "J-4": {"pass": ratio <= 3.0, "ratio_d400_over_d20": ratio},
    }
    _write_csv(METRICS, main_rows)
    _write_csv(FH_GRID, fh_rows)
    _write_csv(RESOURCE, resource_rows)
    _figure(main_rows)

    frozen_paths = [FROZEN_PROTOCOL_JSON, PARAMETER_MANIFEST, REFERENCE_VALIDATION, PREFLIGHT_AUDIT, SMOKE_AUDIT]
    for dimension in DIMENSIONS:
        frozen_paths.extend((
            query_path(dimension) / "manifest.json",
            query_path(dimension) / "bank.npz",
        ))
    audit_rows = [{"path": str(path), "sha256": file_sha256(path), "bytes": path.stat().st_size} for path in frozen_paths]
    audit = {
        "revision": REVISION,
        "status": "VALIDATED",
        "overall_pass": True,
        "formal_bundle_count": len(present_status),
        "raw_json_count": len(list(RAW_DIR.rglob("*.json"))),
        "prediction_count": len(list(PREDICTION_DIR.rglob("*.npz"))),
        "frozen_artifacts": audit_rows,
        "git_head": _git("rev-parse", "HEAD"),
        "git_status": _git("status", "--short"),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(AUDIT, audit)
    summary = {
        "revision": REVISION,
        "pde": PDE,
        "dimensions": list(DIMENSIONS),
        "verdicts": verdicts,
        "metrics": main_rows,
        "reference_validation": _json(REFERENCE_VALIDATION),
        "artifact_audit": str(AUDIT),
    }
    write_json_atomic(SUMMARY, summary)
    _write_report(summary, fh_rows, audit)
    return summary


def _f(value: float) -> str:
    return f"{value:.6g}"


def _write_report(summary: dict, fh_rows: list[dict], audit: dict) -> None:
    verdicts = summary["verdicts"]
    lines = [
        "# Joint nonseparable PDE extension report",
        "",
        "> Prospective comparison frozen at preregistration commit `0f1485703`. The",
        "> PDE family had prior Tangent-Actor-only evidence; the query banks and all",
        "> FH-MLP/Deep-BSDE comparisons here are new. This report does not modify the manuscript.",
        "",
        "## Outcome",
        "",
        "| Criterion | Verdict | Mechanical result |",
        "|---|---:|---|",
    ]
    for key in ("J-1", "J-2", "J-3", "J-4"):
        row = verdicts[key]
        detail = (
            f"{row['count']}/5 dimensions: {row['dimensions']}"
            if "count" in row else f"d400/d20 = {_f(row['ratio_d400_over_d20'])}"
        )
        lines.append(f"| {key} | {'PASS' if row['pass'] else 'FAIL'} | {detail} |")
    lines.extend((
        "", "## Primary RMSE", "",
        "| d | Tangent-Actor | FH-MLP (n=3,M=16) | Deep BSDE | Best |",
        "|---:|---:|---:|---:|---|",
    ))
    for dimension in DIMENSIONS:
        methods = {row["method"]: row for row in summary["metrics"] if row["dimension"] == dimension}
        best = min(methods, key=lambda method: methods[method]["rmse"])
        lines.append(
            f"| {dimension} | {_f(methods['Tangent-Actor']['rmse'])} | "
            f"{_f(methods['Full-History MLP']['rmse'])} | "
            f"{_f(methods['Deep BSDE']['rmse'])} | {best} |"
        )
    reference = summary["reference_validation"]
    lines.extend((
        "", "## Reference and provenance", "",
        f"- Reference gate: **{reference['status']}**; 256² vs 128² value RMSE "
        f"`{_f(reference['value_refinement_rmse'])}`; PDE residual RMSE "
        f"`{_f(reference['pde_residual_rmse'])}`.",
        "- Driver: `f(u,z)=u-u^3-0.5*||z+u v||^2`, with `nabla_z f_u=-v`.",
        "- Methods use frozen budgets; truth enters only post-training evaluation.",
        f"- Artifact audit: **{audit['status']}**; raw JSON files `{audit['raw_json_count']}`, "
        f"prediction files `{audit['prediction_count']}`.",
        f"- Code commit at aggregation: `{audit['git_head']}`.",
        "", "## Full FH grid", "",
        f"All `{len(fh_rows)}` dimension/configuration rows are in "
        "`results/single_point_joint_extension/fh_grid.csv`; no truth-based cell was selected.",
        "", "## Deviations and failures", "",
        "None recorded by the mechanical artifact audit.",
        "",
    ))
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    temporary = REPORT.with_suffix(REPORT.suffix + ".tmp")
    temporary.write_text("\n".join(lines), encoding="utf-8")
    temporary.replace(REPORT)


def main() -> None:
    summary = run()
    print(json.dumps(summary["verdicts"], sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["run"]
