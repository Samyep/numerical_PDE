"""Build canonical ICML-suite tables from audited reuse and new final runs.

This script is intentionally idempotent.  It ignores every ``smoke`` artifact,
keeps native Deep-BSDE rows separate from the amortized whole-function table,
and never converts a missing run into a numerical value.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

import numpy as np
import pandas as pd

from .protocol import RESULTS_DIR, SOURCE_COMMITS


PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
ORIGINAL_EXP = SOURCE_ROOT / "invariant_region_mlp/experiments/original_pde_benchmark_comparison"
FRONTIER = SOURCE_ROOT / "invariant_region_mlp/results/tangent_actor_accuracy_frontier"
RAW_DIR = RESULTS_DIR / "raw"

PDE_ALIASES = {
    "dpi_ou_gmm": "ou_gmm",
    "coupled": "nonconvex_coupled",
}
METHOD_ALIASES = {
    "terminal_aware_tangent_actor_mcv": "tangent_actor",
    "dpi": "dpi",
    "pinn_hte": "pinn_hte",
    "sdgd_pinn_official": "sdgd_pinn",
    "sdgd_pinn_official_transfer": "sdgd_pinn",
    "deep_bsde_official_pytorch_compatibility_port": "deep_bsde",
}
METHOD_LABELS = {
    "tangent_actor": "Tangent-Actor + MCV",
    "dpi": "DPI",
    "pinn_hte": "PINN-HTE",
    "sdgd_pinn": "SDGD-PINN",
    "deep_bsde": "Deep BSDE",
}
HISTORICAL_BANK_HASHES = {
    "hjb_log": "7f3f3b201d4f137712ffca28a546383527143b537ea16418966a909c2aac5ffb",
    "ou_gmm": "14426c76454e0ac11133090a188b782343a297702e0ca4e8f106748bbff028d4",
    "allen_cahn": "731dafc5aac912af3fa584857548e2e649b0eb82a0ea42c10fefb837f4daaad8",
}
REFERENCE_NAMES = {
    "hjb_log": "independent Hopf-Cole quadrature",
    "ou_gmm": "exact evolving Gaussian-mixture density",
    "joint_concave": "independently refined periodic 2D PDE",
    "joint_convex": "independently refined periodic 2D PDE",
    "nonconvex_coupled": "independently refined Fourier 2D PDE",
    "mlp_failure": "exact one-dimensional Gaussian expectation",
    "allen_cahn": "dimension-specific radial PDE Richardson reference",
}
REFERENCE_UNCERTAINTY_DEFAULTS = {
    "hjb_log": 5.33e-15,
    "ou_gmm": 1.19e-12,
    "joint_concave": 8.453967703870988e-05,
    "joint_convex": 7.78411804651081e-05,
    "nonconvex_coupled": 8.079690914758104e-06,
    # This closed-form reference is limited only by floating-point evaluation.
    "mlp_failure": 0.0,
    "allen_cahn": 9.644675185722928e-06,
}


def _finite(value: Any) -> Any:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, np.generic):
        return value.item()
    return value


def _compact(value: Any) -> str | None:
    if value is None:
        return None
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _max_key(value: Any, key: str) -> float | None:
    found: list[float] = []
    if isinstance(value, dict):
        for name, item in value.items():
            if name == key:
                try:
                    number = float(item)
                    if np.isfinite(number):
                        found.append(number)
                except (TypeError, ValueError):
                    pass
            nested = _max_key(item, key)
            if nested is not None:
                found.append(nested)
    elif isinstance(value, (list, tuple)):
        for item in value:
            nested = _max_key(item, key)
            if nested is not None:
                found.append(nested)
    return max(found) if found else None


def _reference_uncertainty(payload: Any, pde: str) -> float | None:
    """Read the numerical-reference floor without depending on one schema."""
    for key in (
        "reference_uncertainty_value_rmse",
        "reference_uncertainty_rmse",
        "reference_uncertainty_proxy",
        "conservative_uncertainty",
    ):
        value = _max_key(payload, key)
        if value is not None:
            return value
    return REFERENCE_UNCERTAINTY_DEFAULTS.get(pde)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _pde_from_benchmark(value: str) -> str:
    base = value.rsplit("_d", 1)[0]
    return PDE_ALIASES.get(base, base)


def _base_row(**overrides: Any) -> dict[str, Any]:
    row = {
        "pde": None, "dimension": None, "method": None, "method_label": None,
        "seed": None, "method_status": None, "execution_status": "complete",
        "benchmark_profile": None, "profile_role": "primary",
        "source_kind": None, "source_commit": None, "source_artifact": None,
        "run_id": None, "checkpoint": None, "selection_rule": None,
        "steps": None, "paths": None, "estimator": None, "states": None,
        "test_bank": None, "test_bank_sha256": None, "reference": None,
        "reference_uncertainty": None, "value_rmse": None,
        "value_relative_l2": None, "value_bias": None, "value_mae": None,
        "gradient_component_rmse": None, "gradient_vector_rmse": None,
        "gradient_relative_l2": None, "actor_value_rmse": None,
        "actor_gradient_vector_rmse": None, "mc_rmse_contribution": None,
        "variance_reduction": None, "weight_ess_fraction": None,
        "training_seconds": None, "model_selection_seconds": None,
        "inference_seconds": None, "total_seconds": None,
        "allocated_gpu_seconds": None, "peak_gpu_bytes": None,
        "parameter_count": None, "nonfinite_count": None,
        "hyperparameters": None, "notes": None,
    }
    row.update(overrides)
    return {key: _finite(value) for key, value in row.items()}


def _historical_original_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    path = ORIGINAL_EXP / "function_approximation_comparison.csv"
    frame = pd.read_csv(path)
    keep_actor = (
        (frame["method"] == "terminal_aware_tangent_actor_mcv")
        & (frame["evaluation_role"] == "frozen_headline")
    )
    keep_baseline = frame["method"].isin(("dpi", "pinn_hte", "sdgd_pinn_official"))
    selected = frame[keep_actor | keep_baseline]
    rows: list[dict[str, Any]] = []
    for _, item in selected.iterrows():
        pde = _pde_from_benchmark(str(item["benchmark"]))
        method = METHOD_ALIASES[str(item["method"])]
        source_commit = (
            "d36ab9e79d0e6ffe49c43d60cdb28a064dc9ac09"
            if method == "tangent_actor" else SOURCE_COMMITS[method]
        )
        rows.append(_base_row(
            pde=pde, dimension=100, method=method,
            method_label=METHOD_LABELS[method], seed=int(item["seed"]),
            method_status=("native_official" if method != "tangent_actor" else "final_method"),
            execution_status="complete_reused", source_kind="historical_exact_reuse",
            source_commit=source_commit, source_artifact=str(path.relative_to(SOURCE_ROOT)),
            run_id=f"historical_{item['benchmark']}_{method}_s{int(item['seed'])}",
            selection_rule=(
                "reference-free frozen Actor selection" if method == "tangent_actor"
                else "published fixed budget"
            ),
            steps=item["time_steps"], paths=(
                int(item["paths_per_replicate"] * item["replicates"])
                if pd.notna(item["paths_per_replicate"]) and pd.notna(item["replicates"])
                else None
            ),
            estimator="mcv" if method == "tangent_actor" and pde != "allen_cahn" else "direct",
            states=int(item["states"]), test_bank_sha256=HISTORICAL_BANK_HASHES[pde],
            reference=REFERENCE_NAMES[pde],
            reference_uncertainty=_reference_uncertainty({}, pde),
            value_rmse=item["value_rmse"],
            value_relative_l2=item["value_relative_l2"], value_bias=item["value_bias"],
            value_mae=item["value_mae"],
            gradient_component_rmse=item["gradient_component_rmse"],
            gradient_vector_rmse=item["gradient_vector_rmse"],
            actor_value_rmse=item["direct_value_actor_rmse"],
            mc_rmse_contribution=(
                np.sqrt(item["sampling_variance_component_mse"])
                if pd.notna(item["sampling_variance_component_mse"]) else None
            ),
            variance_reduction=item["variance_reduction"],
            weight_ess_fraction=item["weight_ess_fraction"],
            training_seconds=item["training_seconds"], inference_seconds=item["inference_seconds"],
            total_seconds=item["total_seconds"], peak_gpu_bytes=item["peak_gpu_bytes"],
            parameter_count=item["parameter_count"], nonfinite_count=item["nonfinite_count"],
            notes="byte- and protocol-audited historical 100D result",
            benchmark_profile=(
                "legacy_hjb_coordinate10" if method == "sdgd_pinn"
                else ("hjb_pis" if method in {"dpi", "pinn_hte"} else "tangent_actor")
            ),
            profile_role=(
                "diagnostic_invalid_coordinate_budget" if method == "sdgd_pinn" else (
                    "diagnostic_legacy_hjb_profile"
                    if method in {"dpi", "pinn_hte"} and pde != "ou_gmm"
                    else "primary"
                )
            ),
        ))

    native_path = ORIGINAL_EXP / "native_task_comparison.csv"
    native = pd.read_csv(native_path)
    allowed = native["method"].isin(METHOD_ALIASES)
    actor_headline = (
        (native["method"] != "terminal_aware_tangent_actor_mcv")
        | (native["variant"] == "focused/a4/n64")
    )
    native_rows: list[dict[str, Any]] = []
    for _, item in native[allowed & actor_headline].iterrows():
        method = METHOD_ALIASES[str(item["method"])]
        pde = _pde_from_benchmark(str(item["benchmark"]))
        native_rows.append({
            "pde": pde, "dimension": 100, "method": method,
            "method_label": METHOD_LABELS[method], "seed": int(item["seed"]),
            "method_status": "native_official" if method == "deep_bsde" else (
                "native_official" if method in {"dpi", "pinn_hte", "sdgd_pinn"} else "final_method"
            ),
            "execution_status": "complete_reused", "source_kind": "historical_exact_reuse",
            "source_commit": SOURCE_COMMITS.get(method, "d36ab9e79d0e6ffe49c43d60cdb28a064dc9ac09"),
            "source_artifact": str(native_path.relative_to(SOURCE_ROOT)),
            "run_id": f"historical_native_{item['benchmark']}_{method}_s{int(item['seed'])}",
            "point": "t=0,x=0", "prediction": item["prediction"],
            "reference": item["reference"], "absolute_error": item["absolute_error"],
            "reference_uncertainty": _reference_uncertainty({}, pde),
            "relative_error": item["relative_error"], "mc_standard_error": item["mc_standard_error"],
            "steps": None, "paths": None, "estimator": None,
            "training_seconds": item["training_seconds"],
            "inference_seconds": item["inference_seconds"], "total_seconds": item["total_seconds"],
            "peak_gpu_bytes": item["peak_gpu_bytes"], "nonfinite_count": item["nonfinite_count"],
            "full_training_cost_charged": method == "deep_bsde",
            "notes": "native single-point task; Deep BSDE is not included in MAIN_RESULTS.csv",
            "benchmark_profile": "published_native" if method == "deep_bsde" else None,
            "profile_role": (
                "diagnostic_invalid_coordinate_budget" if method == "sdgd_pinn" else (
                    "diagnostic_legacy_hjb_profile"
                    if method in {"dpi", "pinn_hte"} and pde != "ou_gmm"
                    else "primary"
                )
            ),
        })
    return rows, native_rows


def _historical_nonconvex_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    euler_path = FRONTIER / "euler_convergence.csv"
    euler = pd.read_csv(euler_path)
    euler = euler[euler["family"] == "coupled"]
    compute = pd.read_csv(FRONTIER / "compute_accounting.csv")
    totals = compute[
        (compute["family"] == "coupled")
        & (compute["cost_category"] == "production_training_and_reference_free_selection")
    ]
    training = {}
    for _, row in totals.iterrows():
        seed = int(str(row["run_id"]).rsplit("_s", 1)[1])
        training[seed] = float(row["seconds"])
    main, frontier = [], []
    for _, item in euler.iterrows():
        seed, steps = int(item["training_seed"]), int(item["steps"])
        row = _base_row(
            pde="nonconvex_coupled", dimension=100, method="tangent_actor",
            method_label=METHOD_LABELS["tangent_actor"], seed=seed,
            method_status="final_method", execution_status="complete_reused",
            source_kind="historical_exact_reuse",
            source_commit="eb5c23d64fe7733a8e48bb3b9b68c682035c66a4",
            source_artifact=str(euler_path.relative_to(SOURCE_ROOT)),
            run_id=f"coupled_d100_mature_s{seed}", selection_rule="reference-free fixed-point consistency",
            steps=steps, paths=int(item["paths"]), estimator="mcv", states=960,
            test_bank_sha256="85a9430c6f60762f5964a9de46122e8010cd1d96b67b286a6997a824d20c2180",
            reference=REFERENCE_NAMES["nonconvex_coupled"],
            reference_uncertainty=item["reference_uncertainty_rmse"],
            value_rmse=item["value_rmse"], value_relative_l2=item["relative_l2"],
            value_bias=item["bias"], mc_rmse_contribution=item["mc_rmse_contribution"],
            variance_reduction=item["mcv_variance_reduction"], weight_ess_fraction=item["ess"],
            training_seconds=training.get(seed), inference_seconds=item["inference_seconds"],
            total_seconds=(training.get(seed, 0.0) + float(item["inference_seconds"])),
            peak_gpu_bytes=item["peak_gpu_bytes"], nonfinite_count=item["nonfinite_count"],
            notes="post-freeze predefined Euler refinement; mature 100D historical anchor",
        )
        frontier.append(row)
        if steps == 64:
            main.append(row.copy())
    return main, frontier


def _json_files(pattern: str) -> list[tuple[Path, dict[str, Any]]]:
    found = []
    for path in sorted(RAW_DIR.glob(pattern)):
        if "smoke" in path.name:
            continue
        try:
            payload = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        if payload.get("status") != "complete":
            continue
        found.append((path, payload))
    return found


def _new_tangent_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    main, frontier, native = [], [], []
    for path, result in _json_files("evaluation_*_mature_s*.json"):
        pde = result.get("pde", result.get("family"))
        d, seed = int(result["dimension"]), int(result["seed"])
        run_id = result.get("run_id", path.stem.removeprefix("evaluation_"))
        training_path = RAW_DIR / f"training_{run_id}.json"
        training = json.loads(training_path.read_text())
        config = training.get("configuration", {})
        invocation_seconds = training.get("invocation_seconds")
        phase_seconds = sum(
            float(phase.get("phase_seconds") or 0.0)
            for phase in training.get("phases", [])
        )
        training_seconds = phase_seconds or invocation_seconds
        model_selection_seconds = (
            max(float(invocation_seconds) - phase_seconds, 0.0)
            if invocation_seconds is not None and phase_seconds else None
        )
        training_peak = _max_key(training, "peak_gpu_bytes")
        selected = result.get("selected_actor")
        checkpoint = training.get("selection", {}).get("selected_checkpoint")
        bank_hash = result.get("test_bank_sha256")
        if bank_hash is None and result.get("test_bank"):
            bank_path = Path(result["test_bank"])
            if bank_path.exists():
                bank_hash = _sha256(bank_path)
        for metric in result["rows"]:
            row = _base_row(
                pde=pde, dimension=d, method="tangent_actor",
                method_label=METHOD_LABELS["tangent_actor"], seed=seed,
                method_status="adapted" if pde in {"joint_convex", "mlp_failure"} else "final_method",
                source_kind="new_final_suite", source_commit=SOURCE_COMMITS["tangent_actor"],
                source_artifact=str(path.relative_to(SOURCE_ROOT)), run_id=run_id,
                checkpoint=checkpoint, selection_rule=training.get("selection", {}).get("criterion"),
                steps=metric["steps"], paths=metric["paths"], estimator=metric["estimator"],
                states=result.get("states"), test_bank=result.get("test_bank"),
                test_bank_sha256=bank_hash, reference=_compact(result.get("reference")),
                reference_uncertainty=_reference_uncertainty(result, pde),
                value_rmse=metric["value_rmse"], value_relative_l2=metric["relative_l2"],
                value_bias=metric["bias"], value_mae=metric["mae"],
                actor_value_rmse=result.get("value_actor_rmse"),
                actor_gradient_vector_rmse=result.get("gradient_actor_vector_rmse"),
                mc_rmse_contribution=metric.get("mc_rmse_contribution"),
                variance_reduction=metric.get("variance_reduction"),
                weight_ess_fraction=metric.get("weight_ess_fraction"),
                training_seconds=training_seconds,
                model_selection_seconds=model_selection_seconds,
                inference_seconds=metric.get("inference_seconds"),
                total_seconds=(float(invocation_seconds or training_seconds or 0.0)
                               + float(metric.get("inference_seconds") or 0.0)),
                peak_gpu_bytes=max(
                    float(training_peak or 0.0), float(metric.get("peak_gpu_bytes") or 0.0),
                ) or None,
                nonfinite_count=metric.get("nonfinite_count"),
                hyperparameters=_compact(config),
                notes=(
                    "minimum-norm nonsmooth subgradient extension" if pde == "mlp_failure"
                    else "post-freeze Tangent-Actor evaluation"
                ),
            )
            frontier.append(row)
            preferred_estimator = "plain" if pde == "allen_cahn" else "mcv"
            if int(metric["steps"]) == 64 and metric["estimator"] == preferred_estimator:
                main.append(row.copy())
        for metric in result.get("native_rows", []):
            native.append({
                "pde": pde, "dimension": d, "method": "tangent_actor",
                "method_label": METHOD_LABELS["tangent_actor"], "seed": seed,
                "method_status": "adapted" if pde in {"joint_convex", "mlp_failure"} else "final_method",
                "execution_status": "complete", "source_kind": "new_final_suite",
                "source_commit": SOURCE_COMMITS["tangent_actor"],
                "source_artifact": str(path.relative_to(SOURCE_ROOT)),
                "run_id": run_id, "point": "t=0,x=0",
                "prediction": metric["prediction"], "reference": metric["reference"],
                "absolute_error": metric["absolute_error"],
                "reference_uncertainty": _reference_uncertainty(result, pde),
                "relative_error": metric["absolute_error"] / max(abs(metric["reference"]), 1e-300),
                "mc_standard_error": metric.get("mc_standard_error"),
                "steps": metric["steps"], "paths": metric["paths"],
                "estimator": metric["estimator"], "training_seconds": training_seconds,
                "model_selection_seconds": model_selection_seconds,
                "inference_seconds": metric.get("inference_seconds"),
                "total_seconds": float(invocation_seconds or training_seconds or 0.0) + float(metric.get("inference_seconds") or 0.0),
                "peak_gpu_bytes": training_peak, "nonfinite_count": metric.get("nonfinite_count"),
                "full_training_cost_charged": False,
                "notes": "post-freeze native point; independent MC stream",
            })
    return main, frontier, native


def _new_function_baselines() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    main, native = [], []
    patterns = ("dpi_*_official_s*.json", "pinn_hte_*_official_s*.json", "sdgd_pinn_*_official_s*.json")
    for pattern in patterns:
        for path, result in _json_files(pattern):
            method = METHOD_ALIASES[result["method"]]
            pde = result.get("pde", result.get("family"))
            training = result.get("training_seconds", result.get("cumulative_training_seconds"))
            row = _base_row(
                pde=pde, dimension=int(result["dimension"]), method=method,
                method_label=METHOD_LABELS[method], seed=int(result["seed"]),
                method_status=result["classification"], source_kind="new_final_suite",
                source_commit=result["official_source_commit"],
                source_artifact=str(path.relative_to(SOURCE_ROOT)), run_id=result["run_id"],
                checkpoint=result.get("checkpoint"), selection_rule="published fixed budget",
                states=result.get("test_states", result.get("states")),
                test_bank=result.get("test_bank"), test_bank_sha256=result.get("test_bank_sha256"),
                reference=_compact(result.get("reference_metadata", REFERENCE_NAMES.get(pde))),
                reference_uncertainty=_reference_uncertainty(result, pde),
                value_rmse=result["value_rmse"], value_relative_l2=result["value_relative_l2"],
                value_bias=result.get("value_bias"), value_mae=result.get("value_mae"),
                gradient_component_rmse=result.get("gradient_component_rmse"),
                gradient_vector_rmse=result.get("gradient_vector_rmse"),
                gradient_relative_l2=result.get("gradient_relative_l2"),
                training_seconds=training, inference_seconds=result.get("inference_seconds"),
                total_seconds=result.get("total_seconds", float(training or 0.0) + float(result.get("inference_seconds") or 0.0)),
                peak_gpu_bytes=result.get("peak_gpu_bytes"), parameter_count=result.get("parameter_count"),
                nonfinite_count=result.get("nonfinite_count", 0),
                hyperparameters=_compact(result.get("configuration", result.get("overrides"))),
                notes="official-method transfer" if result["classification"] != "native_official" else "native official reproduction",
                benchmark_profile=result.get("benchmark_profile", "legacy_unprofiled"),
                profile_role=result.get("profile_role", (
                    "primary" if pde == "ou_gmm" and method in {"dpi", "pinn_hte"}
                    else (
                        "diagnostic_invalid_coordinate_budget" if method == "sdgd_pinn"
                        else "diagnostic_legacy_hjb_profile"
                    )
                )),
            )
            main.append(row)
            if result.get("native"):
                item = result["native"]
                native.append({
                    "pde": pde, "dimension": int(result["dimension"]), "method": method,
                    "method_label": METHOD_LABELS[method], "seed": int(result["seed"]),
                    "method_status": result["classification"], "execution_status": "complete",
                    "source_kind": "new_final_suite", "source_commit": result["official_source_commit"],
                    "source_artifact": str(path.relative_to(SOURCE_ROOT)), "run_id": result["run_id"],
                    "point": item.get("point", "t=0,x=0"), "prediction": item["prediction"],
                    "reference": item["reference"], "absolute_error": item["absolute_error"],
                    "reference_uncertainty": _reference_uncertainty(result, pde),
                    "relative_error": item["absolute_error"] / max(abs(item["reference"]), 1e-300),
                    "mc_standard_error": None, "steps": None, "paths": None, "estimator": "network",
                    "training_seconds": training, "inference_seconds": result.get("inference_seconds"),
                    "total_seconds": result.get("total_seconds"), "peak_gpu_bytes": result.get("peak_gpu_bytes"),
                    "nonfinite_count": result.get("nonfinite_count", 0),
                    "full_training_cost_charged": True,
                    "notes": "same frozen network at native point",
                    "benchmark_profile": row["benchmark_profile"],
                    "profile_role": row["profile_role"],
                })
    return main, native


def _new_deep_bsde_rows() -> list[dict[str, Any]]:
    rows = []
    for path, result in _json_files("deep_bsde_*_official_s*.json"):
        pde = result["family"]
        rows.append({
            "pde": pde, "dimension": int(result["dimension"]), "method": "deep_bsde",
            "method_label": METHOD_LABELS["deep_bsde"], "seed": int(result["seed"]),
            "method_status": result["classification"], "execution_status": "complete",
            "source_kind": "new_final_suite", "source_commit": result["official_source_commit"],
            "source_artifact": str(path.relative_to(SOURCE_ROOT)), "run_id": result["run_id"],
            "point": "t=0,x=0", "prediction": result["reported_y0"],
            "reference": result["native_reference"], "absolute_error": result["native_absolute_error"],
            "reference_uncertainty": _reference_uncertainty(result, pde),
            "relative_error": result["native_relative_error"], "mc_standard_error": None,
            "steps": result["configuration"]["intervals"], "paths": None,
            "estimator": "native scalar Y0", "training_seconds": result["training_seconds_this_invocation"],
            "inference_seconds": 0.0, "total_seconds": result["training_seconds_this_invocation"],
            "peak_gpu_bytes": result.get("peak_gpu_bytes"), "nonfinite_count": result.get("nonfinite_count"),
            "full_training_cost_charged": True,
            "notes": "native single-initial-state solve; excluded from whole-function table",
            "benchmark_profile": "published_native",
            "profile_role": "primary",
        })
    return rows


def _git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True).strip()


def _deduplicate(rows: list[dict[str, Any]], keys: tuple[str, ...]) -> list[dict[str, Any]]:
    # New suite executions supersede a reused row only when the scientific key
    # is identical.  Sorting puts new rows last.
    priority = {"historical_exact_reuse": 0, "new_final_suite": 1}
    result: dict[tuple[Any, ...], dict[str, Any]] = {}
    for row in sorted(rows, key=lambda item: priority.get(item.get("source_kind"), 0)):
        result[tuple(row.get(key) for key in keys)] = row
    return list(result.values())


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n")
        return
    frame = pd.DataFrame(rows)
    preferred = list(rows[0])
    extras = sorted(set(frame.columns) - set(preferred))
    frame = frame[preferred + extras]
    sort = [column for column in ("pde", "dimension", "method", "seed", "steps", "paths") if column in frame]
    frame.sort_values(sort, na_position="last").to_csv(path, index=False)


def _dimension_scaling(main: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in main:
        grouped[(row["pde"], row["dimension"], row["method"], row["method_status"])].append(row)
    output = []
    for (pde, d, method, status), rows in grouped.items():
        errors = np.asarray([row["value_rmse"] for row in rows], dtype=float)
        relative = np.asarray([row["value_relative_l2"] for row in rows], dtype=float)
        costs = np.asarray([row["total_seconds"] for row in rows if row.get("total_seconds") is not None], dtype=float)
        output.append({
            "pde": pde, "dimension": d, "method": method,
            "method_label": METHOD_LABELS[method], "method_status": status,
            "seeds": len(rows), "completed_seed_ids": ";".join(str(row["seed"]) for row in rows),
            "value_rmse_mean": float(errors.mean()),
            "value_rmse_sample_std": float(errors.std(ddof=1)) if len(errors) > 1 else None,
            "value_relative_l2_mean": float(relative.mean()),
            "value_relative_l2_sample_std": float(relative.std(ddof=1)) if len(relative) > 1 else None,
            "total_seconds_mean": float(costs.mean()) if len(costs) else None,
        })
    return output


def _ablation_rows(frontier: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    # Paired estimators share the frozen Actor, Brownian budget and test bank.
    # Require both variants explicitly: a lone MCV frontier point is not an
    # ablation and must not be presented as one.
    paired: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in frontier:
        if row["method"] == "tangent_actor" and row.get("estimator") in {"plain", "mcv"}:
            paired[(
                row["pde"], row["dimension"], row["seed"], row["steps"],
                row["paths"], row["source_artifact"],
            )].append(row)
    for candidates in paired.values():
        if {row["estimator"] for row in candidates} != {"plain", "mcv"}:
            continue
        for row in candidates:
            rows.append({
                "component": "martingale_control_variate", "pde": row["pde"],
                "dimension": row["dimension"], "seed": row["seed"],
                "variant": row["estimator"], "iteration": None,
                "steps": row["steps"], "paths": row["paths"],
                "value_rmse": row["value_rmse"],
                "value_relative_l2": row["value_relative_l2"],
                "mc_rmse_contribution": row["mc_rmse_contribution"],
                "variance_reduction": row["variance_reduction"],
                "actor_value_rmse": row["actor_value_rmse"],
                "actor_gradient_vector_rmse": row["actor_gradient_vector_rmse"],
                "selection_metric": None, "training_seconds": row["training_seconds"],
                "source_artifact": row["source_artifact"],
                "evidence_role": "paired post-freeze inference ablation",
            })

    # Original 100D package stored plain and MCV errors in the same paired row.
    original_path = ORIGINAL_EXP / "function_approximation_comparison.csv"
    original = pd.read_csv(original_path)
    original = original[
        (original["method"] == "terminal_aware_tangent_actor_mcv")
        & (original["evaluation_role"] == "frozen_headline")
    ]
    for _, item in original.iterrows():
        pde = _pde_from_benchmark(str(item["benchmark"]))
        # The Allen--Cahn implementation is deliberately value-only and has no
        # learned gradient integrand.  Its identical direct columns therefore
        # do not constitute an MCV comparison.
        if pde == "allen_cahn":
            continue
        for variant, error in (("plain", item["plain_value_rmse"]), ("mcv", item["value_rmse"])):
            rows.append({
                "component": "martingale_control_variate", "pde": pde,
                "dimension": 100, "seed": int(item["seed"]), "variant": variant,
                "iteration": "a4", "steps": item["time_steps"],
                "paths": int(item["paths_per_replicate"] * item["replicates"]),
                "value_rmse": error,
                "value_relative_l2": item["value_relative_l2"] if variant == "mcv" else None,
                "mc_rmse_contribution": (
                    np.sqrt(item["sampling_variance_component_mse"])
                    if variant == "mcv" and pd.notna(item["sampling_variance_component_mse"]) else None
                ),
                "variance_reduction": item["variance_reduction"],
                "actor_value_rmse": item["direct_value_actor_rmse"],
                "actor_gradient_vector_rmse": item["gradient_vector_rmse"],
                "selection_metric": None, "training_seconds": item["training_seconds"],
                "source_artifact": str(original_path.relative_to(SOURCE_ROOT)),
                "evidence_role": "paired post-freeze inference ablation; exact historical reuse",
            })

    # Mature nonconvex run exposes every update under a common post-freeze bank.
    for seed in (0, 1, 2):
        path = FRONTIER / f"post_freeze_coupled_d100_k2_lam0p2_eta1_main_mcv_s{seed}.json"
        if not path.exists():
            continue
        payload = json.loads(path.read_text())
        selection_path = FRONTIER / f"selection_coupled_d100_k2_lam0p2_eta1_main_mcv_s{seed}.json"
        selection = {}
        if selection_path.exists():
            selection_payload = json.loads(selection_path.read_text())
            selection = {
                row["actor"]: row["normalized_total_inconsistency"]
                for row in selection_payload.get("primary_diagnostics", [])
            }
        for item in payload["iteration_rows"]:
            for variant in ("plain", "mcv"):
                rows.append({
                    "component": "martingale_control_variate", "pde": "nonconvex_coupled",
                    "dimension": 100, "seed": seed, "variant": variant,
                    "iteration": item["iteration"], "steps": item["time_steps"],
                    "paths": item["paths"],
                    "value_rmse": item[f"{variant}_value_rmse"],
                    "value_relative_l2": item[f"{variant}_relative_l2"],
                    "mc_rmse_contribution": item[f"{variant}_mc_rmse_contribution"],
                    "variance_reduction": item["mcv_variance_reduction"],
                    "actor_value_rmse": item["actor_value_rmse"],
                    "actor_gradient_vector_rmse": item["actor_gradient_vector_rmse"],
                    "selection_metric": selection.get(item["iteration"]), "training_seconds": None,
                    "source_artifact": str(path.relative_to(SOURCE_ROOT)),
                    "evidence_role": "paired post-freeze inference ablation; mature historical reuse",
                })
            rows.append({
                "component": "learned_tangent_updates", "pde": "nonconvex_coupled",
                "dimension": 100, "seed": seed, "variant": "terminal_aware_tangent_actor",
                "iteration": item["iteration"], "steps": item["time_steps"],
                "paths": item["paths"], "value_rmse": item["mcv_value_rmse"],
                "value_relative_l2": item["mcv_relative_l2"],
                "mc_rmse_contribution": item["mcv_mc_rmse_contribution"],
                "variance_reduction": item["mcv_variance_reduction"],
                "actor_value_rmse": item["actor_value_rmse"],
                "actor_gradient_vector_rmse": item["actor_gradient_vector_rmse"],
                "selection_metric": selection.get(item["iteration"]), "training_seconds": None,
                "source_artifact": str(path.relative_to(SOURCE_ROOT)),
                "evidence_role": "mature update-by-update diagnostic; test truth not used for selection",
            })

    # New final-suite joint-PDE runs use the same post-freeze protocol for all
    # a0,...,a4 checkpoints and attach their frozen reference-free consistency.
    for path, payload in _json_files("iteration_evaluation_*.json"):
        family = PDE_ALIASES.get(payload["family"], payload["family"])
        for item in payload["rows"]:
            rows.append({
                "component": "martingale_control_variate", "pde": family,
                "dimension": payload["dimension"], "seed": payload["seed"],
                "variant": item["estimator"], "iteration": item["actor"],
                "steps": item["steps"], "paths": item["paths"],
                "value_rmse": item["value_rmse"],
                "value_relative_l2": item["relative_l2"],
                "mc_rmse_contribution": item["mc_rmse_contribution"],
                "variance_reduction": item["variance_reduction"],
                "actor_value_rmse": item["actor_value_rmse"],
                "actor_gradient_vector_rmse": item["actor_gradient_vector_rmse"],
                "selection_metric": item["reference_free_total_inconsistency"],
                "training_seconds": None,
                "source_artifact": str(path.relative_to(SOURCE_ROOT)),
                "evidence_role": "paired post-freeze all-iteration diagnostic",
            })
            if item["estimator"] == "mcv":
                rows.append({
                    "component": "learned_tangent_updates", "pde": family,
                    "dimension": payload["dimension"], "seed": payload["seed"],
                    "variant": "terminal_aware_tangent_actor",
                    "iteration": item["actor"], "steps": item["steps"],
                    "paths": item["paths"], "value_rmse": item["value_rmse"],
                    "value_relative_l2": item["relative_l2"],
                    "mc_rmse_contribution": item["mc_rmse_contribution"],
                    "variance_reduction": item["variance_reduction"],
                    "actor_value_rmse": item["actor_value_rmse"],
                    "actor_gradient_vector_rmse": item["actor_gradient_vector_rmse"],
                    "selection_metric": item["reference_free_total_inconsistency"],
                    "training_seconds": None,
                    "source_artifact": str(path.relative_to(SOURCE_ROOT)),
                    "evidence_role": (
                        "post-freeze all-iteration diagnostic; truth excluded from selection"
                    ),
                })

    # Earlier controlled pilots isolate terminal-aware parameterization. They
    # are supporting ablations and never populate the main benchmark table.
    terminal_path = SOURCE_ROOT / "invariant_region_mlp/results/tangent_actor_strong_baselines/convergence_diagnostics.csv"
    if terminal_path.exists():
        diagnostic = pd.read_csv(terminal_path)
        diagnostic = diagnostic[diagnostic["method"].isin((
            "standard_tangent_actor", "tangent_actor", "terminal_tangent_Ag",
        ))]
        for _, item in diagnostic.iterrows():
            pde = "joint_concave" if str(item["benchmark"]).startswith("joint_concave") else "nonseparable_supporting"
            rows.append({
                "component": "terminal_aware_parameterization", "pde": pde,
                "dimension": int(str(item["benchmark"]).split("_d", 1)[1].split("_", 1)[0]),
                "seed": int(item["seed"]), "variant": item["method"],
                "iteration": "last", "steps": None, "paths": None,
                "value_rmse": item["last_rmse_posthoc"], "value_relative_l2": None,
                "mc_rmse_contribution": None, "variance_reduction": None,
                "actor_value_rmse": None, "actor_gradient_vector_rmse": None,
                "selection_metric": item["validation_objective_last"],
                "training_seconds": item["last_training_seconds"],
                "source_artifact": str(terminal_path.relative_to(SOURCE_ROOT)),
                "evidence_role": "supporting earlier controlled ablation; excluded from main table",
            })
    return rows


def _update_manifest(main: list[dict[str, Any]], native: list[dict[str, Any]]) -> dict[str, Any]:
    path = RESULTS_DIR / "EXECUTION_MANIFEST.json"
    manifest = json.loads(path.read_text())
    completed: dict[tuple[str, int, str], set[int]] = defaultdict(set)
    for row in main + native:
        completed[(row["pde"], int(row["dimension"]), row["method"])].add(int(row["seed"]))
    for cell in manifest["method_cells"]:
        key = (cell["pde"], int(cell["dimension"]), cell["method"])
        primary_seeds = set(completed.get(key, set()))
        previously_recorded = {int(value) for value in cell.get("completed_seeds", [])}
        seeds = sorted(primary_seeds | previously_recorded)
        cell["completed_seeds"] = seeds
        cell["historical_completed_seeds"] = sorted(previously_recorded - {0})
        if seeds:
            target = {int(value) for value in str(cell["target_training_seeds"]).split(";")}
            matching = [
                row for row in main + native
                if (row["pde"], int(row["dimension"]), row["method"]) == key
            ]
            reused_only = all(row.get("source_kind") == "historical_exact_reuse" for row in matching)
            prefix = "completed_reused" if reused_only else "completed"
            cell["execution_status"] = (
                prefix if target.issubset(seeds) else f"{prefix}_partial_seed_target"
            )
            artifacts: dict[str, dict[str, Any]] = {}
            for artifact in cell.get("artifacts", []):
                if isinstance(artifact, dict) and artifact.get("path"):
                    artifacts[artifact["path"]] = artifact
                elif isinstance(artifact, str):
                    candidate = SOURCE_ROOT / artifact
                    artifacts[artifact] = {
                        "path": artifact,
                        "bytes": candidate.stat().st_size if candidate.exists() else None,
                        "sha256": _sha256(candidate) if candidate.exists() else None,
                    }
            for row in matching:
                artifact = row["source_artifact"]
                candidate = SOURCE_ROOT / artifact
                artifacts[artifact] = {
                    "path": artifact,
                    "bytes": candidate.stat().st_size if candidate.exists() else None,
                    "sha256": _sha256(candidate) if candidate.exists() else None,
                }
            cell["artifacts"] = [artifacts[name] for name in sorted(artifacts)]
        elif cell["execution_status"].startswith("completed"):
            # Old reuse labels are retained only if the canonical importer
            # actually found their audited numerical rows.
            cell["execution_status"] = "pending"
        elif (
            cell["execution_status"] == "pending_validity_preflight"
            and cell["pde"] == "mlp_failure"
        ):
            # Stage-A smoke coverage exercised each source-compatible
            # minimum-norm subgradient path.  This changes only the execution
            # gate; the methods remain explicitly classified as adapted.
            cell["execution_status"] = "pending"
    counts = defaultdict(int)
    for cell in manifest["method_cells"]:
        counts[cell["execution_status"]] += 1
    manifest["status"] = "stage_b_fidelity_gate_pending"
    manifest["updated_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["execution_counts"] = dict(sorted(counts.items()))
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def aggregate() -> dict[str, Any]:
    historical_main, historical_native = _historical_original_rows()
    nonconvex_main, nonconvex_frontier = _historical_nonconvex_rows()
    new_tangent_main, new_tangent_frontier, new_tangent_native = _new_tangent_rows()
    new_baseline_main, new_baseline_native = _new_function_baselines()
    new_deep_bsde_native = _new_deep_bsde_rows()

    all_main = historical_main + nonconvex_main + new_tangent_main + new_baseline_main
    all_native = historical_native + new_tangent_native + new_baseline_native + new_deep_bsde_native

    # The amended protocol has one fixed primary seed. Historical seed 1/2
    # and superseded baseline profiles remain immutable evidence but cannot
    # enter primary rankings, means, figures, or manifest completion logic.
    primary_main_candidates = [
        row for row in all_main
        if int(row["seed"]) == 0 and row.get("profile_role", "primary") == "primary"
    ]
    primary_native_candidates = [
        row for row in all_native
        if int(row["seed"]) == 0 and row.get("profile_role", "primary") == "primary"
    ]
    main = _deduplicate(
        primary_main_candidates,
        ("pde", "dimension", "method", "seed"),
    )
    native = _deduplicate(
        primary_native_candidates,
        ("pde", "dimension", "method", "seed", "steps", "estimator"),
    )
    excluded = []
    # Compare against the final de-duplicated rows, not every candidate.  This
    # keeps superseded-but-valid primary runs visible in the historical export
    # after a corrected profile replaces the same scientific cell.
    primary_main_ids = {id(row) for row in main}
    primary_native_ids = {id(row) for row in native}
    for table_kind, rows, kept in (
        ("whole_function", all_main, primary_main_ids),
        ("native_point", all_native, primary_native_ids),
    ):
        for row in rows:
            if id(row) in kept:
                continue
            excluded.append({
                **row, "table_kind": table_kind,
                "primary_exclusion_reason": (
                    "historical_nonzero_seed" if int(row["seed"]) != 0
                    else (
                        row.get("profile_role", "nonprimary_profile")
                        if row.get("profile_role", "primary") != "primary"
                        else "superseded_primary_duplicate"
                    )
                ),
            })
    frontier_candidates = [
        row for row in historical_main + nonconvex_frontier + new_tangent_frontier + new_baseline_main
        if int(row["seed"]) == 0 and row.get("profile_role", "primary") == "primary"
    ]
    frontier = _deduplicate(
        frontier_candidates,
        ("pde", "dimension", "method", "seed", "steps", "paths", "estimator"),
    )
    scaling = _dimension_scaling(main)
    ablations = _ablation_rows(frontier)
    _write_csv(RESULTS_DIR / "MAIN_RESULTS.csv", main)
    _write_csv(RESULTS_DIR / "NATIVE_POINT_RESULTS.csv", native)
    _write_csv(RESULTS_DIR / "ACCURACY_COMPUTE_FRONTIER.csv", frontier)
    _write_csv(RESULTS_DIR / "DIMENSION_SCALING.csv", scaling)
    _write_csv(RESULTS_DIR / "ABLATION_RESULTS.csv", ablations)
    _write_csv(RESULTS_DIR / "HISTORICAL_AND_NONPRIMARY_RESULTS.csv", excluded)
    review_rows = [
        {**row, "comparison_task": "whole_function"}
        for row in main if int(row["dimension"]) == 100
    ] + [
        {**row, "comparison_task": "single_initial_state"}
        for row in native if int(row["dimension"]) == 100
    ]
    _write_csv(RESULTS_DIR / "100D_RESULTS.csv", review_rows)
    manifest = _update_manifest(main, native)
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": _git_head(), "main_rows": len(main), "native_rows": len(native),
        "frontier_rows": len(frontier), "dimension_scaling_rows": len(scaling),
        "ablation_rows": len(ablations),
        "historical_and_nonprimary_rows": len(excluded),
        "review_100d_rows": len(review_rows),
        "completed_method_cells": sum(
            item["execution_status"].startswith("completed")
            for item in manifest["method_cells"]
        ),
        "total_method_cells": len(manifest["method_cells"]),
        "smoke_artifacts_included": False,
        "deep_bsde_in_main_results": False,
        "primary_training_seed": 0,
        "best_seed_selection": False,
    }
    (RESULTS_DIR / "AGGREGATION_STATUS.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    print(json.dumps(aggregate(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
