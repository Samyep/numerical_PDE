"""Build the canonical 140-cell whole-function benchmark deliverables."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from .paired_bank_repair import array_sha256, file_sha256
from .protocol import PDES, RESULTS_DIR, SOURCE_COMMITS
from .whole_function_status import DIMENSIONS, METHODS, attempts, key
from .whole_function_job_accounting import collect as collect_job_accounting


SOURCE_ROOT = Path(__file__).resolve().parents[3]
FIGURE_DIR = RESULTS_DIR / "whole_function_figures"
TABLE_DIR = RESULTS_DIR / "whole_function_tables"
PDE_ORDER = tuple(item.key for item in PDES)
PDE_LABELS = {item.key: item.title for item in PDES}
METHOD_LABELS = {
    "tangent_actor": "Tangent-Actor + MCV",
    "dpi": "DPI",
    "pinn_hte": "PINN-HTE",
    "sdgd_pinn": "SDGD-PINN",
}
COLORS = {
    "tangent_actor": "#C73E1D", "dpi": "#8C5AAE",
    "pinn_hte": "#2A6FBB", "sdgd_pinn": "#2E8B57",
}
MARKERS = {"tangent_actor": "o", "dpi": "D", "pinn_hte": "s", "sdgd_pinn": "^"}


def _read_main() -> pd.DataFrame:
    data = pd.read_csv(RESULTS_DIR / "MAIN_RESULTS.csv")
    return data[
        data["pde"].isin(PDE_ORDER)
        & data["dimension"].isin(DIMENSIONS)
        & data["method"].isin(METHODS)
        & (data["seed"] == 0)
        & (data["profile_role"] == "primary")
    ].copy()


def _attempt_groups() -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for item in attempts():
        try:
            grouped[key(item["method"], item["pde"], int(item["dimension"]))].append(item)
        except (KeyError, TypeError, ValueError):
            continue
    return grouped


def _terminal_status(items: list[dict[str, Any]]) -> tuple[str, bool]:
    terminal = [
        item for item in items
        if item.get("status") not in {"RUNNING", "CHECKPOINTED"}
    ]
    running = any(item.get("status") == "RUNNING" for item in items)
    if not terminal:
        return "PENDING", running
    terminal.sort(key=lambda item: (str(item.get("finished_utc") or ""), str(item.get("status_file") or "")))
    value = str(terminal[-1].get("status") or "FAILED")
    return value if value in {"COMPLETE", "OUT_OF_TIME", "OUT_OF_MEMORY", "FAILED"} else "FAILED", running


def build_results_frame() -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    main = _read_main()
    duplicate = main.duplicated(["pde", "dimension", "method", "seed"], keep=False)
    if duplicate.any():
        raise RuntimeError("duplicate primary rows: " + repr(
            main.loc[duplicate, ["pde", "dimension", "method", "seed"]].to_dict("records")
        ))
    indexed = {
        (row.pde, int(row.dimension), row.method): row._asdict()
        for row in main.itertuples(index=False)
    }
    grouped_attempts = _attempt_groups()
    rows: list[dict[str, Any]] = []
    for pde in PDE_ORDER:
        for dimension in DIMENSIONS:
            for method in METHODS:
                source = indexed.get((pde, dimension, method))
                cell_attempts = grouped_attempts.get(key(method, pde, dimension), [])
                terminal, active = _terminal_status(cell_attempts)
                status = "COMPLETE" if source is not None else terminal
                attempt_seconds = sum(
                    float(item.get("actual_gpu_seconds") or 0.0)
                    for item in cell_attempts if item.get("status") != "RUNNING"
                )
                measured_seconds = (
                    attempt_seconds if attempt_seconds > 0
                    else float(source.get("allocated_gpu_seconds") or source.get("total_seconds") or 0.0)
                    if source is not None else 0.0
                )
                row: dict[str, Any] = {
                    "pde": pde, "pde_label": PDE_LABELS[pde],
                    "dimension": dimension, "method": method,
                    "method_label": METHOD_LABELS[method], "seed": 0,
                    "status": status, "active_slurm_attempt": active,
                    "training_seed_count": 1 if source is not None else 0,
                    "attempt_count": len([x for x in cell_attempts if x.get("status") != "RUNNING"]),
                    "cell_gpu_seconds": measured_seconds,
                    "cell_gpu_hours": measured_seconds / 3600.0,
                    "cell_gpu_hour_cap": 6.0,
                }
                copied = (
                    "value_rmse", "value_relative_l2", "value_bias", "value_mae",
                    "states", "test_bank", "test_bank_sha256", "test_bank_array_sha256",
                    "test_bank_file_sha256", "reference", "reference_uncertainty",
                    "reference_array_sha256", "reference_source", "prediction_file",
                    "prediction_file_sha256", "source_commit", "source_artifact",
                    "method_status", "benchmark_profile", "profile_role", "selection_rule",
                    "training_seconds", "model_selection_seconds", "inference_seconds",
                    "total_seconds", "peak_gpu_bytes", "parameter_count", "nonfinite_count",
                    "metric_provenance", "notes",
                )
                for name in copied:
                    row[name] = source.get(name) if source is not None else None
                if status != "COMPLETE":
                    for name in ("value_rmse", "value_relative_l2", "value_bias", "value_mae"):
                        row[name] = None
                rows.append(row)
    frame = pd.DataFrame(rows)
    if len(frame) != 140:
        raise RuntimeError(f"expected 140 cells, got {len(frame)}")
    return frame, [item for group in grouped_attempts.values() for item in group]


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")


def _actual_artifact_hashes(row: pd.Series) -> dict[str, Any]:
    path = Path(str(row["prediction_file"]))
    if not path.is_absolute():
        path = SOURCE_ROOT / path
    if not path.is_file():
        return {"valid": False, "error": f"missing prediction artifact {path}"}
    try:
        with np.load(path, allow_pickle=False) as data:
            if "t" in data and "x" in data:
                t, x = np.asarray(data["t"]), np.asarray(data["x"])
            elif "test_t" in data and "test_x" in data:
                t, x = np.asarray(data["test_t"]), np.asarray(data["test_x"])
            else:
                raise KeyError("missing t/x arrays")
            truth = np.asarray(data["truth"]).reshape(-1)
        return {
            "valid": bool(len(t) == len(x) == len(truth) and np.isfinite(t).all()
                          and np.isfinite(x).all() and np.isfinite(truth).all()),
            "states": len(t), "bank": array_sha256(t, x),
            "reference": array_sha256(truth), "file": file_sha256(path),
        }
    except Exception as exc:
        return {"valid": False, "error": repr(exc)}


def audit(frame: pd.DataFrame, attempt_rows: list[dict[str, Any]]) -> dict[str, Any]:
    completed = frame[frame["status"] == "COMPLETE"].copy()
    artifact_checks: list[dict[str, Any]] = []
    for _, row in completed.iterrows():
        actual = _actual_artifact_hashes(row)
        passed = bool(
            actual.get("valid")
            and int(actual.get("states", -1)) == int(row["states"])
            and actual.get("bank") == row["test_bank_array_sha256"]
            and actual.get("reference") == row["reference_array_sha256"]
            and actual.get("file") == row["prediction_file_sha256"]
        )
        artifact_checks.append({
            "pde": row["pde"], "dimension": int(row["dimension"]),
            "method": row["method"], "status": "PASS" if passed else "FAIL",
            "actual": actual,
        })

    group_checks = []
    for (pde, dimension), rows in completed.groupby(["pde", "dimension"]):
        bank = set(rows["test_bank_array_sha256"].dropna().astype(str))
        reference = set(rows["reference_array_sha256"].dropna().astype(str))
        states = set(rows["states"].dropna().astype(int))
        passed = len(bank) == len(reference) == len(states) == 1
        group_checks.append({
            "pde": pde, "dimension": int(dimension), "completed_methods": len(rows),
            "status": "PASS" if passed else "FAIL", "bank_hashes": sorted(bank),
            "reference_hashes": sorted(reference), "state_counts": sorted(states),
        })

    status_counts = Counter(frame["status"])
    active_compute_seconds = sum(
        float(item.get("actual_gpu_seconds") or 0.0)
        for item in attempt_rows if item.get("status") != "RUNNING"
    )
    job_accounting = collect_job_accounting()
    new_gpu_seconds = float(job_accounting["allocated_gpu_seconds"])
    complete_metrics = completed[["value_rmse", "value_relative_l2"]].apply(
        pd.to_numeric, errors="coerce",
    ).to_numpy(float)
    terminal = frame[frame["status"] != "COMPLETE"]
    checks = {
        "exactly_140_unique_cells": len(frame) == 140 and not frame.duplicated(
            ["pde", "dimension", "method"]
        ).any(),
        "only_frozen_dimensions": set(frame["dimension"]) == set(DIMENSIONS),
        "only_four_whole_function_methods": set(frame["method"]) == set(METHODS),
        "seed_zero_only": set(frame["seed"]) == {0},
        "complete_metrics_finite": bool(np.isfinite(complete_metrics).all()),
        "incomplete_cells_have_no_accuracy": bool(
            terminal[["value_rmse", "value_relative_l2"]].isna().all().all()
        ),
        "prediction_artifacts_match_reported_hashes": all(
            item["status"] == "PASS" for item in artifact_checks
        ),
        "same_actual_bank_and_reference_within_each_completed_group": all(
            item["status"] == "PASS" for item in group_checks
        ),
        "source_commits_present": bool(completed["source_commit"].notna().all()),
        "no_nonfinite_completed_runs": bool(
            (pd.to_numeric(completed["nonfinite_count"], errors="coerce").fillna(0) == 0).all()
        ),
        "per_cell_new_compute_within_six_hours": all(
            sum(float(item.get("actual_gpu_seconds") or 0.0) for item in items
                if item.get("status") != "RUNNING") <= 21600
            for items in _attempt_groups().values()
        ),
        "new_compute_within_120_gpu_hours": new_gpu_seconds <= 120 * 3600,
    }
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "suite_status": "complete" if int(status_counts.get("PENDING", 0)) == 0 else "in_progress",
        "checks": checks,
        "status_counts": {name: int(status_counts.get(name, 0)) for name in
                          ("COMPLETE", "OUT_OF_TIME", "OUT_OF_MEMORY", "FAILED", "PENDING")},
        "new_gpu_hours": new_gpu_seconds / 3600.0,
        "summed_cell_active_gpu_hours": active_compute_seconds / 3600.0,
        "job_accounting": job_accounting,
        "artifact_checks": artifact_checks,
        "within_group_pairing_checks": group_checks,
    }


def _configure_plot() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.2, "pdf.fonttype": 42,
    })


def _save(fig: plt.Figure, name: str) -> list[Path]:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    paths = [FIGURE_DIR / f"{name}.pdf", FIGURE_DIR / f"{name}.png"]
    fig.savefig(paths[0], bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(paths[1], dpi=300, bbox_inches="tight")
    plt.close(fig)
    return paths


def _scaling_figure(frame: pd.DataFrame, metric: str, name: str, ylabel: str) -> list[Path]:
    done = frame[frame["status"] == "COMPLETE"]
    fig, axes = plt.subplots(3, 3, figsize=(10.5, 8.0), sharex=True)
    for axis, pde in zip(axes.flat, PDE_ORDER, strict=False):
        subset = done[done["pde"] == pde]
        for method in METHODS:
            rows = subset[subset["method"] == method].sort_values("dimension")
            if rows.empty:
                continue
            axis.plot(rows["dimension"], rows[metric], color=COLORS[method],
                      marker=MARKERS[method], linewidth=1.5, label=METHOD_LABELS[method])
        if not subset.empty:
            axis.set_yscale("log")
        else:
            axis.text(0.5, 0.5, "pending", transform=axis.transAxes,
                      ha="center", va="center", color="#667085")
        axis.set_title(PDE_LABELS[pde], fontsize=9)
        axis.set_xticks(DIMENSIONS)
    for axis in axes.flat[len(PDE_ORDER):]:
        axis.axis("off")
    for axis in axes[:, 0]:
        axis.set_ylabel(ylabel)
    for axis in axes[-1, :]:
        axis.set_xlabel("Ambient dimension $d$")
    handles = [Line2D([0], [0], color=COLORS[m], marker=MARKERS[m], label=METHOD_LABELS[m]) for m in METHODS]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False)
    fig.suptitle("Completed cells only; missing/OOT cells have no numerical marker", y=0.965, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return _save(fig, name)


def _compute_figure(frame: pd.DataFrame) -> list[Path]:
    done = frame[(frame["status"] == "COMPLETE") & (frame["cell_gpu_hours"] > 0)]
    fig, axis = plt.subplots(figsize=(7.3, 4.8))
    for method in METHODS:
        rows = done[done["method"] == method]
        axis.scatter(rows["cell_gpu_hours"], rows["value_relative_l2"],
                     c=COLORS[method], marker=MARKERS[method], s=35,
                     label=METHOD_LABELS[method], alpha=0.82)
    if not done.empty:
        axis.set_xscale("log"); axis.set_yscale("log")
    axis.set_xlabel("Measured per-cell GPU hours")
    axis.set_ylabel("Value relative $L^2$")
    axis.set_title("Whole-function accuracy–compute endpoints")
    axis.legend(frameon=False, ncol=2)
    return _save(fig, "whole_function_accuracy_compute")


def _latex(frame: pd.DataFrame) -> Path:
    indexed = frame.set_index(["pde", "dimension", "method"])
    lines = [
        r"\begin{tabular}{llcccc}", r"\toprule",
        r"PDE & $d$ & Tangent-Actor & DPI & PINN-HTE & SDGD-PINN \\", r"\midrule",
    ]
    for pde in PDE_ORDER:
        for index, dimension in enumerate(DIMENSIONS):
            cells = []
            for method in METHODS:
                row = indexed.loc[(pde, dimension, method)]
                if row["status"] == "COMPLETE":
                    cells.append(rf"\shortstack{{{row['value_relative_l2']:.3g} rel.\\{row['value_rmse']:.3g} RMSE}}")
                else:
                    cells.append(str(row["status"]).replace("_", r"\_"))
            label = PDE_LABELS[pde].replace("&", r"\&") if index == 0 else ""
            lines.append(f"{label} & {dimension} & " + " & ".join(cells) + r" \\")
        lines.append(r"\addlinespace")
    lines.extend((r"\bottomrule", r"\end{tabular}"))
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    path = TABLE_DIR / "whole_function_main_table.tex"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _rankings(frame: pd.DataFrame) -> list[dict[str, Any]]:
    records = []
    for (pde, dimension), rows in frame.groupby(["pde", "dimension"]):
        complete = rows[rows["status"] == "COMPLETE"]
        if len(complete) != len(METHODS):
            continue
        ordered = complete.sort_values("value_rmse")
        records.append({
            "pde": pde, "dimension": int(dimension),
            "ranking": ordered["method"].tolist(),
            "rmse": {row.method: float(row.value_rmse) for row in ordered.itertuples()},
            "winner": ordered.iloc[0]["method"],
        })
    return records


def _report(frame: pd.DataFrame, audit_payload: dict[str, Any]) -> Path:
    rankings = _rankings(frame)
    ours = sum(item["winner"] == "tangent_actor" for item in rankings)
    counts = audit_payload["status_counts"]
    lines = [
        "# Tangent-Actor whole-function benchmark",
        "",
        f"**Integrity audit:** {audit_payload['status']}. **Execution:** "
        f"{audit_payload['suite_status']}. Generated {audit_payload['generated_utc']}.",
        "",
        "This benchmark contains the predeclared 7 PDE x 5 dimension x 4 method grid. "
        "Every training seed is 0. Deep BSDE is excluded because it is a separate single-point task.",
        "",
        "## Execution coverage",
        "",
        f"- COMPLETE: {counts['COMPLETE']}",
        f"- OUT_OF_TIME: {counts['OUT_OF_TIME']}",
        f"- OUT_OF_MEMORY: {counts['OUT_OF_MEMORY']}",
        f"- FAILED: {counts['FAILED']}",
        f"- PENDING: {counts['PENDING']}",
        f"- New measured GPU time: {audit_payload['new_gpu_hours']:.3f} GH200 hours (120-hour cap).",
        "",
        "Incomplete cells have blank accuracy fields and are absent from all numerical curves.",
        "",
        "## Strictly complete rankings",
        "",
        f"All four methods are complete in {len(rankings)} PDE/dimension settings. "
        f"Tangent-Actor has the smallest Value RMSE in {ours}/{len(rankings)} of them.",
        "",
    ]
    for item in rankings:
        order = " < ".join(f"{name} ({item['rmse'][name]:.4g})" for name in item["ranking"])
        lines.append(f"- {PDE_LABELS[item['pde']]}, d={item['dimension']}: {order}")
    lines += [
        "", "## Interpretation limits", "",
        "The baseline rows are source-faithful method transfers unless explicitly labeled native. "
        "The table therefore supports comparison on this frozen whole-function suite, rather than "
        "a claim of universal superiority over each method's original paper task. Results are "
        "single-seed training endpoints; numerical-reference and finite-MC uncertainty remain in the CSV.",
    ]
    path = RESULTS_DIR / "WHOLE_FUNCTION_FINAL_REPORT.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def build() -> dict[str, Any]:
    _configure_plot()
    frame, attempt_rows = build_results_frame()
    results_path = RESULTS_DIR / "WHOLE_FUNCTION_RESULTS.csv"
    _write_csv(frame, results_path)
    scaling = frame[frame["status"] == "COMPLETE"].copy()
    scaling_path = RESULTS_DIR / "WHOLE_FUNCTION_DIMENSION_SCALING.csv"
    _write_csv(scaling, scaling_path)
    compute_columns = [
        "pde", "dimension", "method", "status", "attempt_count",
        "training_seconds", "model_selection_seconds", "inference_seconds",
        "total_seconds", "cell_gpu_seconds", "cell_gpu_hours", "cell_gpu_hour_cap",
    ]
    compute_path = RESULTS_DIR / "WHOLE_FUNCTION_COMPUTE.csv"
    _write_csv(frame[compute_columns], compute_path)
    audit_payload = audit(frame, attempt_rows)
    audit_path = RESULTS_DIR / "WHOLE_FUNCTION_AUDIT.json"
    audit_path.write_text(json.dumps(audit_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    status_payload = {
        "generated_utc": audit_payload["generated_utc"],
        "suite_status": audit_payload["suite_status"],
        "grid_cells": 140, "status_counts": audit_payload["status_counts"],
        "new_gpu_hours": audit_payload["new_gpu_hours"],
        "cell_gpu_hour_cap": 6.0, "suite_gpu_hour_cap": 120.0,
        "active_cells": frame.loc[frame["active_slurm_attempt"],
                                  ["pde", "dimension", "method"]].to_dict("records"),
    }
    status_path = RESULTS_DIR / "WHOLE_FUNCTION_EXECUTION_STATUS.json"
    status_path.write_text(json.dumps(status_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    figures = []
    figures += _scaling_figure(frame, "value_rmse", "whole_function_accuracy_comparison", "Value RMSE")
    figures += _scaling_figure(frame, "value_relative_l2", "whole_function_dimension_scaling", "Value relative $L^2$")
    figures += _compute_figure(frame)
    table = _latex(frame)
    report = _report(frame, audit_payload)
    return {
        "status": audit_payload["status"], "counts": audit_payload["status_counts"],
        "new_gpu_hours": audit_payload["new_gpu_hours"], "results": str(results_path),
        "figures": [str(path) for path in figures], "table": str(table), "report": str(report),
    }


if __name__ == "__main__":
    print(json.dumps(build(), indent=2, sort_keys=True))
