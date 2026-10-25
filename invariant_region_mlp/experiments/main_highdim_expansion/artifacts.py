"""Build the canonical compact deliverables for the two-PDE expansion.

The builder is deliberately useful while the suite is running: absent cells
remain explicit PENDING/terminal rows and never receive fabricated metrics.
For Tangent-Actor, the primary result is exactly the frozen selected Actor,
64 Euler steps, inference replicate 0.  Replicate 1 and the 128/256-step runs
remain separate numerical-uncertainty diagnostics.
"""

from __future__ import annotations

from collections import Counter
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

from .accounting import collect as collect_accounting
from .metrics import value_metrics_from_prediction_file
from .prepare import file_sha256
from .protocol import (
    ACTOR_UPDATES,
    DIMENSIONS,
    FIGURE_DIR,
    INFERENCE_REPLICATES,
    METHODS,
    PDES,
    PRIMARY_EULER_STEPS,
    RAW_DIR,
    REPORT_PATH,
    RESULTS_DIR,
    SOURCE_COMMITS,
    SOURCE_ROOT,
    TRAINING_SEED,
    frozen_protocol,
    write_json_atomic,
)
from .status import attempts, cell_budget


METHOD_LABELS = {
    "tangent_actor": "Tangent-Actor + MCV",
    "dpi": "DPI",
    "pinn_hte": "PINN-HTE",
    "sdgd_pinn": "SDGD-PINN",
}
PDE_LABELS = {
    "b6_fullrank": "B6 Full-Rank HJB",
    "dense_fourier": "Dense Fourier Joint",
}
COLORS = {
    "tangent_actor": "#C73E1D",
    "dpi": "#7A5195",
    "pinn_hte": "#2A6FBB",
    "sdgd_pinn": "#2E8B57",
}
MARKERS = {"tangent_actor": "o", "dpi": "D", "pinn_hte": "s", "sdgd_pinn": "^"}
TABLE_DIR = RESULTS_DIR / "tables"
MAIN_PATH = RESULTS_DIR / "main_whole_function_results.csv"
ITERATION_PATH = RESULTS_DIR / "actor_iteration_results.csv"
SUMMARY_PATH = RESULTS_DIR / "final_summary.json"
AUDIT_PATH = RESULTS_DIR / "artifact_audit.json"


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None


def _finite(value: Any) -> Any:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return value
    return number if math.isfinite(number) else None


def _json_safe(value: Any) -> Any:
    """Map non-finite tabular placeholders to JSON null without changing data."""

    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _manifest_cells() -> dict[tuple[str, int], dict[str, Any]]:
    payload = _read_json(RESULTS_DIR / "parameter_and_test_bank_manifest.json") or {}
    return {
        (str(row["pde"]), int(row["dimension"])): row
        for row in payload.get("cells", [])
    }


def _attempt_seconds(method: str, pde: str, dimension: int) -> float:
    return sum(
        float(row.get("actual_gpu_seconds") or 0.0)
        for row in attempts()
        if row.get("method") == method
        and row.get("pde") == pde
        and int(row.get("dimension", -1)) == dimension
        and row.get("status") != "RUNNING"
    )


def _training_payload(pde: str, dimension: int) -> dict[str, Any] | None:
    run_id = (
        f"b6_fullrank_d{dimension}_mature_s0"
        if pde == "b6_fullrank"
        else f"dense_fourier_d{dimension}_mature_s0"
    )
    return _read_json(RAW_DIR / "actor" / f"training_{run_id}.json")


def _implementation_commit() -> str:
    calibration = _read_json(RESULTS_DIR / "b6_calibration_freeze.json") or {}
    return str(
        calibration.get("code_commit_at_calibration")
        or SOURCE_COMMITS["icml_final_suite"]
    )


def _ta_evaluation_path(
    pde: str, dimension: int, actor: str, steps: int, replicate: int
) -> Path:
    run_id = (
        f"b6_fullrank_d{dimension}_mature_s0"
        if pde == "b6_fullrank"
        else f"dense_fourier_d{dimension}_mature_s0"
    )
    return (
        RAW_DIR
        / "actor_evaluations"
        / run_id
        / actor
        / f"steps{steps}_rep{replicate:02d}.json"
    )


def _baseline_path(method: str, pde: str, dimension: int) -> Path:
    if method in {"dpi", "pinn_hte"}:
        return RAW_DIR / f"{method}_{pde}_d{dimension}_semilinear_paper_official_s0.json"
    profile = "hjb_adversarial" if pde == "b6_fullrank" else "semilinear_sde"
    return RAW_DIR / f"sdgd_pinn_{pde}_d{dimension}_{profile}_official_s0.json"


def _prediction_fields(payload: dict[str, Any]) -> tuple[str | None, str | None]:
    path = payload.get("prediction_path") or payload.get("prediction_file")
    digest = payload.get("prediction_sha256") or payload.get("prediction_file_sha256")
    return (str(path) if path else None, str(digest) if digest else None)


def _resolved_artifact_path(raw_path: Any) -> Path | None:
    if raw_path is None or pd.isna(raw_path):
        return None
    path = Path(str(raw_path))
    return path if path.is_absolute() else SOURCE_ROOT / path


def _value_metrics_from_prediction_file(raw_path: Any) -> dict[str, float]:
    path = _resolved_artifact_path(raw_path)
    if path is None:
        raise FileNotFoundError(path)
    return value_metrics_from_prediction_file(path)


def _ta_row(
    pde: str,
    dimension: int,
    manifest: dict[str, Any] | None,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    primary = _read_json(
        _ta_evaluation_path(pde, dimension, "selected", PRIMARY_EULER_STEPS, 0)
    )
    repeat = _read_json(
        _ta_evaluation_path(pde, dimension, "selected", PRIMARY_EULER_STEPS, 1)
    )
    training = _training_payload(pde, dimension)
    extra: dict[str, Any] = {}
    if repeat:
        extra = {
            "independent_replicate_value_rmse": repeat.get("value_rmse"),
            "independent_replicate_relative_l2": repeat.get("relative_l2"),
            "independent_replicate_mc_rmse_contribution": repeat.get("mc_rmse_contribution"),
            "absolute_rmse_replicate_difference": (
                abs(float(primary["value_rmse"]) - float(repeat["value_rmse"]))
                if primary and primary.get("value_rmse") is not None
                and repeat.get("value_rmse") is not None else None
            ),
        }
    if not primary or not training or training.get("status") != "complete":
        return None, extra
    prediction, prediction_hash = _prediction_fields(primary)
    derived = _value_metrics_from_prediction_file(prediction)
    timing = training.get("timing", {})
    training_seconds = float(timing.get("phase_total_seconds") or 0.0)
    selection_seconds = float(timing.get("reference_free_selection_seconds") or 0.0)
    return {
        "value_rmse": primary.get("value_rmse"),
        "value_relative_l2": primary.get("relative_l2"),
        "value_bias": primary.get("bias"),
        "value_mae": primary.get("mae"),
        "value_rmae": derived["value_rmae"],
        "gradient_vector_rmse": primary.get("gradient_vector_rmse")
        or primary.get("actor_gradient_vector_rmse"),
        "gradient_relative_l2": primary.get("gradient_relative_l2"),
        "actor_value_rmse": primary.get("actor_value_rmse"),
        "mc_rmse_contribution": primary.get("mc_rmse_contribution"),
        "variance_reduction": primary.get("variance_reduction"),
        "weight_ess_fraction": primary.get("weight_ess_fraction"),
        "steps": primary.get("steps"),
        "paths": primary.get("paths"),
        "states": primary.get("test_states"),
        "test_bank_array_sha256": primary.get("test_bank_sha256"),
        "test_bank_file": primary.get("test_bank_file"),
        "test_bank_file_sha256": primary.get("test_bank_file_sha256"),
        "reference_array_sha256": manifest.get("truth_value_sha256") if manifest else None,
        "prediction_file": prediction,
        "prediction_file_sha256": prediction_hash,
        "training_seconds": training_seconds,
        "model_selection_seconds": selection_seconds,
        "inference_seconds": primary.get("inference_seconds"),
        "total_seconds": training_seconds + selection_seconds + float(primary.get("inference_seconds") or 0.0),
        "peak_gpu_bytes": primary.get("peak_gpu_bytes"),
        "parameter_count": training.get("parameter_count"),
        "nonfinite_count": primary.get("nonfinite_count"),
        "selected_actor": training.get("selection", {}).get("selected_actor"),
        "selection_rule": training.get("selection", {}).get("criterion"),
        "classification": "proposed_method",
        "benchmark_profile": "mature_terminal_aware_tangent_actor_mcv",
        "source_commit": _implementation_commit(),
        "ground_truth_used_for_training_or_selection": False,
        **extra,
    }, extra


def _baseline_row(
    method: str, pde: str, dimension: int, manifest: dict[str, Any] | None
) -> dict[str, Any] | None:
    payload = _read_json(_baseline_path(method, pde, dimension))
    if not payload or payload.get("status") != "complete":
        return None
    prediction, prediction_hash = _prediction_fields(payload)
    derived = _value_metrics_from_prediction_file(prediction)
    training_seconds = (
        payload.get("training_seconds")
        if method != "sdgd_pinn"
        else payload.get("cumulative_training_seconds")
    )
    inference_seconds = payload.get("inference_seconds")
    return {
        "value_rmse": payload.get("value_rmse"),
        "value_relative_l2": payload.get("value_relative_l2"),
        "value_bias": payload.get("value_bias"),
        "value_mae": payload.get("value_mae"),
        "value_rmae": derived["value_rmae"],
        "gradient_vector_rmse": payload.get("gradient_vector_rmse"),
        "gradient_relative_l2": payload.get("gradient_relative_l2"),
        "actor_value_rmse": None,
        "mc_rmse_contribution": None,
        "variance_reduction": None,
        "weight_ess_fraction": None,
        "steps": None,
        "paths": None,
        "states": payload.get("test_states") or payload.get("states"),
        "test_bank_array_sha256": payload.get("test_bank_array_sha256")
        or payload.get("test_bank_sha256"),
        "test_bank_file": payload.get("test_bank"),
        "test_bank_file_sha256": payload.get("test_bank_file_sha256"),
        "reference_array_sha256": payload.get("reference_array_sha256")
        or (manifest.get("truth_value_sha256") if manifest else None),
        "prediction_file": prediction,
        "prediction_file_sha256": prediction_hash,
        "training_seconds": training_seconds,
        "model_selection_seconds": 0.0,
        "inference_seconds": inference_seconds,
        "total_seconds": payload.get("total_seconds")
        or float(training_seconds or 0.0) + float(inference_seconds or 0.0),
        "peak_gpu_bytes": payload.get("peak_gpu_bytes"),
        "parameter_count": payload.get("parameter_count"),
        "nonfinite_count": payload.get("nonfinite_count", 0),
        "selected_actor": None,
        "selection_rule": payload.get("profile_selection_rule"),
        "classification": payload.get("classification"),
        "benchmark_profile": payload.get("benchmark_profile"),
        "source_commit": payload.get("official_source_commit"),
        "ground_truth_used_for_training_or_selection": False,
        "independent_replicate_value_rmse": None,
        "independent_replicate_relative_l2": None,
        "independent_replicate_mc_rmse_contribution": None,
        "absolute_rmse_replicate_difference": None,
    }


def build_main_frame() -> pd.DataFrame:
    manifest = _manifest_cells()
    rows: list[dict[str, Any]] = []
    for pde in PDES:
        for dimension in DIMENSIONS:
            frozen = manifest.get((pde, dimension))
            for method in METHODS:
                budget = cell_budget(method, pde, dimension)
                measured = _attempt_seconds(method, pde, dimension)
                if method == "tangent_actor":
                    payload, _ = _ta_row(pde, dimension, frozen)
                else:
                    payload = _baseline_row(method, pde, dimension, frozen)
                status = "COMPLETE" if payload is not None else budget["status"]
                if status == "COMPLETE" and payload is None:
                    status = "PENDING"
                base: dict[str, Any] = {
                    "pde": pde,
                    "pde_label": PDE_LABELS[pde],
                    "dimension": dimension,
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "seed": TRAINING_SEED,
                    "status": status,
                    "attempt_count": budget["attempt_count"],
                    "cell_gpu_seconds": measured,
                    "cell_gpu_hours": measured / 3600.0,
                    "cell_gpu_hour_cap": 6.0,
                    "test_bank_array_sha256": frozen.get("test_input_sha256") if frozen else None,
                    "test_bank_file": frozen.get("test_bank_file") if frozen else None,
                    "test_bank_file_sha256": frozen.get("test_bank_file_sha256") if frozen else None,
                    "reference_array_sha256": frozen.get("truth_value_sha256") if frozen else None,
                    "states": frozen.get("test_states") if frozen else None,
                }
                if payload:
                    base.update(payload)
                    if measured <= 0:
                        internal = float(payload.get("total_seconds") or 0.0)
                        base["cell_gpu_seconds"] = internal
                        base["cell_gpu_hours"] = internal / 3600.0
                else:
                    for name in (
                        "value_rmse", "value_relative_l2", "value_bias", "value_mae", "value_rmae",
                        "gradient_vector_rmse", "gradient_relative_l2", "actor_value_rmse",
                        "mc_rmse_contribution", "variance_reduction", "weight_ess_fraction",
                        "steps", "paths", "prediction_file", "prediction_file_sha256",
                        "training_seconds", "model_selection_seconds", "inference_seconds",
                        "total_seconds", "peak_gpu_bytes", "parameter_count", "nonfinite_count", "selected_actor",
                        "selection_rule", "classification", "benchmark_profile", "source_commit",
                        "ground_truth_used_for_training_or_selection",
                        "independent_replicate_value_rmse", "independent_replicate_relative_l2",
                        "independent_replicate_mc_rmse_contribution",
                        "absolute_rmse_replicate_difference",
                    ):
                        base[name] = None
                rows.append(base)
    frame = pd.DataFrame(rows)
    if len(frame) != len(PDES) * len(DIMENSIONS) * len(METHODS):
        raise AssertionError(len(frame))
    if frame.duplicated(["pde", "dimension", "method"]).any():
        raise RuntimeError("duplicate expansion cells")
    return frame


def build_iteration_frame() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for pde in PDES:
        for dimension in DIMENSIONS:
            training = _training_payload(pde, dimension) or {}
            selected = training.get("selection", {}).get("selected_actor")
            for actor in tuple(f"a{index}" for index in range(ACTOR_UPDATES + 1)) + ("selected",):
                step_values = (64, 128, 256) if actor == "selected" else (64,)
                for steps in step_values:
                    for replicate in range(INFERENCE_REPLICATES):
                        payload = _read_json(
                            _ta_evaluation_path(pde, dimension, actor, steps, replicate)
                        )
                        row = {
                            "pde": pde,
                            "dimension": dimension,
                            "training_seed": 0,
                            "actor_request": actor,
                            "resolved_actor": (
                                selected if actor == "selected" else actor
                            ),
                            "selected_by_reference_free_rule": actor == "selected",
                            "steps": steps,
                            "paths": 1024,
                            "inference_replicate": replicate,
                            "status": "COMPLETE" if payload else "PENDING",
                        }
                        for target, source in (
                            ("value_rmse", "value_rmse"),
                            ("value_relative_l2", "relative_l2"),
                            ("value_bias", "bias"),
                            ("mc_rmse_contribution", "mc_rmse_contribution"),
                            ("variance_reduction", "variance_reduction"),
                            ("weight_ess_fraction", "weight_ess_fraction"),
                            ("actor_value_rmse", "actor_value_rmse"),
                            ("actor_gradient_vector_rmse", "actor_gradient_vector_rmse"),
                            ("gradient_vector_rmse", "gradient_vector_rmse"),
                        ):
                            row[target] = payload.get(source) if payload else None
                        rows.append(row)
    return pd.DataFrame(rows)


def _write_benchmark_definitions() -> None:
    manifest = _read_json(RESULTS_DIR / "parameter_and_test_bank_manifest.json")
    payload = {
        "protocol": frozen_protocol(),
        "definitions": {
            "b6_fullrank": {
                "equation": "u_t + Delta u - 0.5 ||z||^2 = 0; z=sqrt(2) grad u",
                "terminal": "-log[(1/(2d)) sum_m exp(b_m^T x)]",
                "parameter_seed": "SeedSequence([20261201,d,66])",
                "exact_transform": "Hopf-Cole",
            },
            "dense_fourier": {
                "equation": "u_t + Delta u + 0.5 tanh(u) ||z||^2 = 0; z=sqrt(2) grad u",
                "terminal": "asinh[0.5 + 2/sqrt(2d) sum_j a_j sin(omega_j^T x+phi_j)]",
                "parameter_seed": "SeedSequence([20261201,d,41])",
                "omega_row_norm": 1.0,
                "exact_transform": "sinh",
                "actor": "joint value-gradient terminal-aware Actor",
            },
        },
        "frozen_manifest_available": manifest is not None,
        "cells": manifest.get("cells", []) if manifest else [],
    }
    write_json_atomic(RESULTS_DIR / "benchmark_definitions.json", payload)


def _write_baseline_audit(frame: pd.DataFrame) -> None:
    completed = frame[(frame["status"] == "COMPLETE") & (frame["method"] != "tangent_actor")]
    lines = [
        "# Baseline Fidelity Audit",
        "",
        "DPI, PINN-HTE, and SDGD-PINN are **faithful method transfers** to the two new PDEs. "
        "They are not claimed as native experiments from the source papers.",
        "",
        "- No terminal-aware Tangent-Actor feature or MCV estimator is added to a baseline.",
        "- DPI/PINN-HTE use the pinned generic semilinear paper profile and its source budget.",
        "- SDGD-PINN uses the source-native adversarial HJB profile for B6 and the source-native "
        "general semilinear SDE profile for Dense Fourier.",
        "- Only the equation, horizon, ambient dimension, and evaluation bank are adapted.",
        "- Exact solutions are opened for frozen post-training evaluation only.",
        "- Smoke runs are non-result artifacts and are excluded from every table and figure.",
        "",
        f"Completed transferred cells: {len(completed)} / {3 * len(PDES) * len(DIMENSIONS)}.",
        "",
        "| Method | Source commit | Classification |",
        "|---|---|---|",
        f"| DPI | `{SOURCE_COMMITS['dpi_pinn_hte']}` | faithful_transfer |",
        f"| PINN-HTE | `{SOURCE_COMMITS['dpi_pinn_hte']}` | faithful_transfer |",
        f"| SDGD-PINN | `{SOURCE_COMMITS['sdgd_pinn']}` | faithful_transfer |",
        "",
    ]
    (RESULTS_DIR / "baseline_fidelity_audit.md").write_text("\n".join(lines), encoding="utf-8")


def _configure_plot() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.22, "pdf.fonttype": 42,
    })


def _save_figure(fig: plt.Figure, name: str) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / f"{name}.pdf", bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(FIGURE_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _figures(frame: pd.DataFrame) -> None:
    _configure_plot()
    done = frame[frame["status"] == "COMPLETE"]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.1), sharex=True)
    for axis, pde in zip(axes, PDES, strict=True):
        for method in METHODS:
            rows = done[(done["pde"] == pde) & (done["method"] == method)].sort_values("dimension")
            if not rows.empty:
                axis.plot(rows["dimension"], rows["value_rmse"], marker=MARKERS[method], color=COLORS[method], label=METHOD_LABELS[method])
        if not done[done["pde"] == pde].empty:
            axis.set_yscale("log")
        else:
            axis.text(0.5, 0.5, "formal cells pending", transform=axis.transAxes, ha="center", va="center")
        axis.set_title(PDE_LABELS[pde]); axis.set_xticks(DIMENSIONS)
        axis.set_xlabel("Ambient dimension $d$")
    axes[0].set_ylabel("Value RMSE")
    handles = [Line2D([0], [0], marker=MARKERS[m], color=COLORS[m], label=METHOD_LABELS[m]) for m in METHODS]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    _save_figure(fig, "figure_a_dimension_scaling")

    characterization_path = RESULTS_DIR / "benchmark_characterization.csv"
    char = pd.read_csv(characterization_path) if characterization_path.is_file() else pd.DataFrame()
    fig, axes = plt.subplots(1, 3, figsize=(12.3, 3.8))
    metrics = (
        ("effective_rank", "Effective rank $r_{eff}$"),
        ("variance_outside_top10", "Variance outside top 10"),
        ("S_nonlinear", "$S_{nonlinear}$"),
    )
    for axis, (metric, label) in zip(axes, metrics, strict=True):
        for pde, color, marker in zip(PDES, ("#D1495B", "#00798C"), ("o", "s"), strict=True):
            rows = char[char["pde"] == pde].sort_values("dimension") if not char.empty else char
            if not rows.empty:
                axis.plot(rows["dimension"], rows[metric], color=color, marker=marker, label=PDE_LABELS[pde])
        if char.empty:
            axis.text(0.5, 0.5, "characterization pending", transform=axis.transAxes, ha="center")
        axis.set_xticks(DIMENSIONS); axis.set_xlabel("Ambient dimension $d$"); axis.set_ylabel(label)
    if not char.empty:
        axes[0].legend(frameon=False)
    fig.tight_layout()
    _save_figure(fig, "figure_b_high_effective_dimensionality")

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.1))
    for axis, pde in zip(axes, PDES, strict=True):
        rows = done[(done["pde"] == pde) & (done["cell_gpu_hours"] > 0)]
        for method in METHODS:
            part = rows[rows["method"] == method]
            if not part.empty:
                axis.scatter(part["cell_gpu_hours"], part["value_rmse"], color=COLORS[method], marker=MARKERS[method], s=45, label=METHOD_LABELS[method])
        if not rows.empty:
            axis.set_xscale("log"); axis.set_yscale("log")
        else:
            axis.text(0.5, 0.5, "formal cells pending", transform=axis.transAxes, ha="center")
        axis.set_title(PDE_LABELS[pde]); axis.set_xlabel("Measured GPU hours / cell")
    axes[0].set_ylabel("Value RMSE")
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    _save_figure(fig, "figure_c_accuracy_compute")


def _latex(frame: pd.DataFrame) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    indexed = frame.set_index(["pde", "dimension", "method"])
    lines = [
        r"\begin{tabular}{llcccc}", r"\toprule",
        r"PDE & $d$ & Tangent-Actor & DPI & PINN-HTE & SDGD-PINN \\", r"\midrule",
    ]
    for pde in PDES:
        for number, dimension in enumerate(DIMENSIONS):
            cells = []
            for method in METHODS:
                row = indexed.loc[(pde, dimension, method)]
                cells.append(
                    f"{float(row['value_rmse']):.3g}"
                    if row["status"] == "COMPLETE" else str(row["status"]).replace("_", r"\_")
                )
            label = PDE_LABELS[pde].replace("-", r"-") if number == 0 else ""
            lines.append(f"{label} & {dimension} & " + " & ".join(cells) + r" \\")
        lines.append(r"\addlinespace")
    lines.extend((r"\bottomrule", r"\end{tabular}"))
    (TABLE_DIR / "main_highdim_expansion.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _combined(frame: pd.DataFrame) -> None:
    selected = frame[[
        "pde", "dimension", "method", "status", "value_rmse", "value_relative_l2", "value_rmae",
        "cell_gpu_hours", "classification", "test_bank_array_sha256",
    ]].copy()
    selected.insert(0, "suite", "two_highdim_expansion")
    old_path = (
        SOURCE_ROOT.parent / "source" / "invariant_region_mlp" / "results"
        / "tangent_actor_icml_final_suite" / "WHOLE_FUNCTION_RESULTS.csv"
    )
    if old_path.is_file():
        old = pd.read_csv(old_path)
        for name in selected.columns:
            if name not in old:
                old[name] = None
        old = old[list(selected.columns)]
        old["suite"] = "original_seven_pde_suite"
        selected = pd.concat((old, selected), ignore_index=True)
    selected.to_csv(RESULTS_DIR / "combined_nine_pde_whole_function_results.csv", index=False, lineterminator="\n")


def _artifact_check(row: pd.Series) -> dict[str, Any]:
    checks: dict[str, Any] = {
        "pde": row["pde"],
        "dimension": int(row["dimension"]),
        "method": row["method"],
    }
    for label, path_name, hash_name in (
        ("bank", "test_bank_file", "test_bank_file_sha256"),
        ("prediction", "prediction_file", "prediction_file_sha256"),
    ):
        raw_path, expected = row.get(path_name), row.get(hash_name)
        path = _resolved_artifact_path(raw_path)
        actual = file_sha256(path) if path is not None and path.is_file() else None
        checks[f"{label}_path"] = str(path) if path is not None else None
        checks[f"{label}_expected_sha256"] = expected if pd.notna(expected) else None
        checks[f"{label}_actual_sha256"] = actual
        checks[f"{label}_passed"] = bool(actual and pd.notna(expected) and actual == expected)
    bank_path = _resolved_artifact_path(row.get("test_bank_file"))
    prediction_path = _resolved_artifact_path(row.get("prediction_file"))
    metric_checks_passed = False
    truth_matches_bank = False
    prediction_shape: list[int] | None = None
    try:
        derived = _value_metrics_from_prediction_file(prediction_path)
        with np.load(prediction_path, allow_pickle=False) as prediction_payload:
            prediction_truth = np.asarray(prediction_payload["truth"], dtype=np.float64)
            prediction_key = "prediction" if "prediction" in prediction_payload.files else "mcv"
            prediction_shape = list(np.asarray(prediction_payload[prediction_key]).shape)
        with np.load(bank_path, allow_pickle=False) as bank_payload:
            bank_truth = np.asarray(bank_payload["truth_u"], dtype=np.float64)
        truth_matches_bank = bool(
            prediction_truth.shape == bank_truth.shape
            and np.array_equal(prediction_truth, bank_truth)
        )
        metric_checks = {
            name: bool(
                pd.notna(row.get(name))
                and np.isclose(float(row[name]), value, rtol=1e-12, atol=1e-15)
            )
            for name, value in derived.items()
        }
        checks["recomputed_metrics"] = derived
        checks["metric_checks"] = metric_checks
        metric_checks_passed = all(metric_checks.values())
    except (FileNotFoundError, KeyError, OSError, ValueError, TypeError) as error:
        checks["metric_error"] = repr(error)
    checks["prediction_shape"] = prediction_shape
    checks["truth_matches_frozen_bank"] = truth_matches_bank
    checks["metrics_recomputed_match"] = metric_checks_passed
    checks["passed"] = bool(
        checks["bank_passed"]
        and checks["prediction_passed"]
        and truth_matches_bank
        and metric_checks_passed
    )
    return checks


def _audit(frame: pd.DataFrame) -> dict[str, Any]:
    completed = frame[frame["status"] == "COMPLETE"]
    incomplete = frame[frame["status"] != "COMPLETE"]
    artifact_checks = [_artifact_check(row) for _, row in completed.iterrows()]
    group_pairing = []
    for (pde, dimension), rows in completed.groupby(["pde", "dimension"]):
        banks = set(rows["test_bank_array_sha256"].dropna().astype(str))
        refs = set(rows["reference_array_sha256"].dropna().astype(str))
        group_pairing.append({
            "pde": pde, "dimension": int(dimension), "methods": len(rows),
            "bank_hashes": sorted(banks), "reference_hashes": sorted(refs),
            "passed": len(banks) == 1 and len(refs) == 1,
        })
    classifications = completed[completed["method"] != "tangent_actor"]["classification"].dropna()
    checks = {
        "exactly_40_unique_cells": len(frame) == 40 and not frame.duplicated(["pde", "dimension", "method"]).any(),
        "seed_zero_only": set(frame["seed"]) == {0},
        "complete_metrics_finite": bool(
            np.isfinite(
                completed[["value_rmse", "value_relative_l2", "value_rmae"]]
                .apply(pd.to_numeric, errors="coerce")
                .to_numpy()
            ).all()
        ),
        "incomplete_cells_have_no_accuracy": bool(
            incomplete[["value_rmse", "value_relative_l2", "value_rmae"]].isna().all().all()
        ),
        "same_bank_and_reference_within_completed_groups": all(row["passed"] for row in group_pairing),
        "completed_artifact_hashes_match": all(row["passed"] for row in artifact_checks),
        "baselines_marked_faithful_transfer": bool((classifications == "faithful_transfer").all()),
        "no_test_truth_in_training_or_selection": bool(
            completed["ground_truth_used_for_training_or_selection"].fillna(False).eq(False).all()
        ),
        "primary_ta_is_selected_steps64_rep0": bool(
            ((completed[completed["method"] == "tangent_actor"]["steps"] == 64).all())
        ),
        "per_cell_allocation_within_six_hours": bool((frame["cell_gpu_hours"] <= 6.0 + 1e-6).all()),
    }
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "suite_status": "complete" if len(completed) == len(frame) else "in_progress_or_terminally_incomplete",
        "checks": checks,
        "status_counts": dict(Counter(frame["status"])),
        "group_pairing": group_pairing,
        "artifact_checks": artifact_checks,
    }


def _summary(frame: pd.DataFrame, audit: dict[str, Any], accounting: dict[str, Any]) -> dict[str, Any]:
    complete_groups = 0
    wins = Counter()
    group_results: list[dict[str, Any]] = []
    for (pde, dimension), rows in frame.groupby(["pde", "dimension"]):
        done = rows[rows["status"] == "COMPLETE"]
        if len(done) != len(METHODS):
            continue
        complete_groups += 1
        best = float(done["value_rmse"].min())
        winners = sorted(done[np.isclose(done["value_rmse"], best, rtol=1e-12, atol=0.0)]["method"])
        for method in winners:
            wins[method] += 1 / len(winners)
        group_results.append({"pde": pde, "dimension": int(dimension), "best_rmse": best, "winner": winners})
    dense_iterations = pd.read_csv(ITERATION_PATH) if ITERATION_PATH.is_file() else pd.DataFrame()
    dense_primary = dense_iterations[
        (dense_iterations.get("pde") == "dense_fourier")
        & (dense_iterations.get("steps") == 64)
        & (dense_iterations.get("inference_replicate") == 0)
        & dense_iterations.get("actor_request", pd.Series(dtype=str)).isin(["a0", "a1", "a4"])
        & (dense_iterations.get("status") == "COMPLETE")
    ] if not dense_iterations.empty else dense_iterations
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "suite_status": "COMPLETE" if (frame["status"] == "COMPLETE").all() else "IN_PROGRESS_OR_TERMINALLY_INCOMPLETE",
        "status_counts": {key: int(value) for key, value in Counter(frame["status"]).items()},
        "complete_four_method_settings": complete_groups,
        "tangent_actor_lowest_rmse_settings": float(wins.get("tangent_actor", 0.0)),
        "wins": dict(wins),
        "group_results": group_results,
        "dense_a0_a1_a4": (
            _json_safe(dense_primary.to_dict("records"))
            if not dense_primary.empty
            else []
        ),
        "allocated_gpu_hours": accounting.get("allocated_gpu_hours"),
        "account_credit_hours": accounting.get("account_credit_hours"),
        "budget_cap_gpu_hours": 120.0,
        "audit_status": audit["status"],
        "claim_scope": "whole-function faithful transfers on frozen shared banks; not native source-paper benchmarks",
    }


def _report(frame: pd.DataFrame, summary: dict[str, Any], audit: dict[str, Any]) -> None:
    counts = summary["status_counts"]
    complete = frame[frame["status"] == "COMPLETE"]
    lines = [
        "# ICML Main High-Dimensional Expansion: B6 and Dense Fourier",
        "",
        f"**Execution state:** {summary['suite_status']}",
        f"**Cells:** {counts.get('COMPLETE', 0)} complete, {counts.get('PENDING', 0)} pending, "
        f"{counts.get('OUT_OF_TIME', 0)} out of time, {counts.get('OUT_OF_MEMORY', 0)} OOM, "
        f"{counts.get('FAILED', 0)} failed, "
        f"{counts.get('BUDGET_TRUNCATED', 0)} budget-truncated.",
        f"**Expansion allocation:** {float(summary.get('allocated_gpu_hours') or 0):.3f} GPU hours "
        f"under the separate 120-hour cap.",
        f"**Artifact audit:** {audit['status']}.",
        "",
        "The primary Tangent-Actor row is the reference-free selected Actor evaluated once with "
        "1,024 paths and 64 Euler steps (replicate 0). Replicate 1 is reported separately and is "
        "not averaged into a fictitious higher-path estimator. Exact solutions are used only after "
        "selection for evaluation and for the pre-frozen B6 calibration/benchmark characterization.",
        "",
        "## Current formal results",
        "",
        "| PDE | d | Method | Status | Value RMSE | Relative L2 | rMAE | GPU h |",
        "|---|---:|---|---|---:|---:|---:|---:|",
    ]
    for row in frame.itertuples(index=False):
        rmse = f"{row.value_rmse:.6g}" if row.status == "COMPLETE" else "—"
        rel = f"{row.value_relative_l2:.6g}" if row.status == "COMPLETE" else "—"
        rmae = f"{row.value_rmae:.6g}" if row.status == "COMPLETE" else "—"
        lines.append(
            f"| {PDE_LABELS[row.pde]} | {row.dimension} | {METHOD_LABELS[row.method]} | "
            f"{row.status} | {rmse} | {rel} | {rmae} | {row.cell_gpu_hours:.3f} |"
        )
    lines.extend([
        "",
        "## Interpretation guardrails",
        "",
        "DPI, PINN-HTE, and SDGD-PINN rows are faithful transfers of pinned source protocols to "
        "new equations. They support a controlled comparison in this suite; they do not establish "
        "that Tangent-Actor universally outperforms every native experiment reported in those papers.",
        "",
        f"Four-method settings currently complete: {summary['complete_four_method_settings']} / 10. "
        f"Tangent-Actor currently has the lowest RMSE in {summary['tangent_actor_lowest_rmse_settings']:.3g} "
        "of those settings (ties split fractionally). This statement will update automatically as cells finish.",
        "",
        "The effective-rank, outside-top-10 variance, nonlinear-skill, iteration, MC-noise, and "
        "accuracy-compute evidence is stored in the accompanying CSVs and Figures A–C. Missing or "
        "capped cells stay explicitly non-numerical.",
        "",
        "## Reproducibility",
        "",
        "- Frozen training seed: 0; frozen whole-function bank size: 960 per PDE/dimension.",
        "- Dimensions: 20, 50, 100, 200, 400.",
        "- Per-cell cap: 6 GPU hours; expansion cap: 120 GPU hours.",
        "- Heavy banks, checkpoints, predictions, and logs remain in the workspace; committed "
        "manifests bind them by SHA-256.",
        "- Baseline classification and source commits are recorded in `baseline_fidelity_audit.md`.",
        "",
    ])
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def build() -> dict[str, Any]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    _write_benchmark_definitions()
    frame = build_main_frame()
    frame.to_csv(MAIN_PATH, index=False, lineterminator="\n")
    iterations = build_iteration_frame()
    iterations.to_csv(ITERATION_PATH, index=False, lineterminator="\n")
    _write_baseline_audit(frame)
    _combined(frame)
    audit = _audit(frame)
    write_json_atomic(AUDIT_PATH, audit)
    try:
        accounting = collect_accounting()
    except (OSError, RuntimeError, __import__("subprocess").SubprocessError) as error:
        accounting = {"error": repr(error), "allocated_gpu_hours": 0.0, "account_credit_hours": 0.0}
    summary = _summary(frame, audit, accounting)
    write_json_atomic(SUMMARY_PATH, summary)
    write_json_atomic(
        RESULTS_DIR / "execution_status.json",
        {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "status_counts": summary["status_counts"],
            "cells": frame[["method", "pde", "dimension", "status", "attempt_count", "cell_gpu_hours"]].to_dict("records"),
        },
    )
    _figures(frame)
    _latex(frame)
    _report(frame, summary, audit)
    print(json.dumps({
        "suite_status": summary["suite_status"],
        "status_counts": summary["status_counts"],
        "audit": audit["status"],
    }, sort_keys=True))
    return summary


def main() -> None:
    build()


if __name__ == "__main__":
    main()
