"""Aggregate the completed frontier study into tables, figures, and reports."""

from __future__ import annotations

import csv
import glob
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import (
    FIGURE_DIR,
    OLD_PREDICTION_DIR,
    OLD_REFERENCE_DIR,
    OLD_RESULTS_DIR,
    PACKAGE_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    REFERENCE_DIR,
    RESULTS_DIR,
    atomic_json,
    ensure_dirs,
    metrics,
    write_csv,
)
from .drivers import SineQuadraticDriver
from .reference2d import load_reference


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _mean_sd(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=np.float64)
    return float(array.mean()), float(array.std(ddof=1)) if len(array) > 1 else 0.0


def _mean_ci95(values: list[float]) -> tuple[float, float]:
    """Small-sample two-sided t interval for independently seeded runs."""
    mean, sd = _mean_sd(values)
    if len(values) <= 1:
        return mean, mean
    critical = {2: 12.706, 3: 4.303}.get(len(values), 1.96)
    half = critical * sd / math.sqrt(len(values))
    return max(0.0, mean - half), mean + half


def _primary_bank_truth() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    driver = SineQuadraticDriver(d=100, lam=0.2, eta=0.0)
    bank_path = OLD_PREDICTION_DIR / (
        "post_freeze_primary_d100_k2_lam0p2_eta0_main_mcv_s0.npz"
    )
    with np.load(bank_path) as bank:
        t, x = bank["test_t"].copy(), bank["test_x"].copy()
    reference = load_reference(REFERENCE_DIR / f"{driver.name}_refined.npz")
    return t, x, reference.value(t, x @ driver.B)


def _full_history_rows() -> list[dict[str, Any]]:
    _, _, primary_truth = _primary_bank_truth()
    rows: list[dict[str, Any]] = []
    for path in sorted((OLD_RESULTS_DIR / "raw").glob("baseline_full_history_*.json")):
        payload = _load(path)
        if payload["spec"]["family"] != "primary" or payload["spec"]["d"] != 100:
            continue
        prediction_path = Path(payload["prediction_path"])
        with np.load(prediction_path) as archive:
            accuracy = metrics(archive["prediction"], primary_truth)
        rows.append({
            "family": "primary", "d": 100, "source": "parent_frozen_grid",
            "role": "parent_grid", "n": payload["n"], "M": payload["M"],
            "seed": payload["baseline_seed"], "finite": payload["runtime"]["finite"],
            "value_rmse": accuracy["rmse"], "relative_l2": accuracy["relative_l2"],
            "bias": accuracy["bias"], "wall_seconds": payload["runtime"]["wall_seconds"],
            "peak_gpu_bytes": payload["runtime"]["peak_gpu_bytes"],
            "reference": "refined_256_recomputed", "memory_tiled": False,
            "nonfinite_count": payload["runtime"]["stats"]["nonfinite_count"],
        })
    for path in sorted(RAW_DIR.glob("full_history_*.json")):
        payload = _load(path)
        accuracy = payload["metrics"]
        rows.append({
            "family": payload["spec"]["family"], "d": payload["spec"]["d"],
            "source": "accuracy_frontier", "role": payload["configuration_role"],
            "n": payload["n"], "M": payload["M"], "seed": payload["baseline_seed"],
            "finite": payload["runtime"]["finite"], "value_rmse": accuracy["rmse"],
            "relative_l2": accuracy["relative_l2"], "bias": accuracy["bias"],
            "wall_seconds": payload["runtime"]["wall_seconds"],
            "peak_gpu_bytes": payload["runtime"]["peak_gpu_bytes"],
            "reference": payload["reference_path"],
            "memory_tiled": payload["runtime"].get("terminal_sample_memory_tiling", False),
            "nonfinite_count": payload["runtime"]["stats"]["nonfinite_count"],
        })
    return rows


def _precision_rows() -> list[dict[str, Any]]:
    rows = []
    payloads = []
    for path in sorted(RAW_DIR.glob("precision_primary_a4_seed*.json")):
        if path.name.startswith("environment_"):
            continue
        payloads.append(_load(path))
    # Retain every raw run locally, while the headline table uses the largest
    # completed path budget for each frozen training seed.
    best_payload = {}
    for payload in payloads:
        seed = int(payload["training_seed"])
        paths = int(payload["configuration"]["paths"])
        if seed not in best_payload or paths > int(best_payload[seed]["configuration"]["paths"]):
            best_payload[seed] = payload
    for payload in [best_payload[key] for key in sorted(best_payload)]:
        refinement = {
            item["coarse_steps"]: item for item in payload["refinement_rows"]
        }
        training = _load(
            OLD_RESULTS_DIR / f"training_primary_d100_k2_lam0p2_eta0_main_mcv_s{payload['training_seed']}.json"
        )
        for row in payload["rows"]:
            rows.append({
                "family": "primary", "d": 100, "method": "Tangent-Actor + MCV",
                "training_seed": payload["training_seed"], "actor": "a4",
                "steps": row["steps"], "paths": row["paths"],
                "value_rmse": row["value_rmse"], "relative_l2": row["relative_l2"],
                "bias": row["mean_bias"], "rmse_ci95_low": row["rmse_ci95_low"],
                "rmse_ci95_high": row["rmse_ci95_high"],
                "mc_rmse_contribution": row["mc_rmse_contribution"],
                "reference_uncertainty_rmse": row["reference_uncertainty_rmse"],
                "residual_error_proxy": row["residual_non_mc_non_reference_rmse_proxy"],
                "mcv_variance_reduction": row["mcv_variance_reduction"],
                "ess": row["weight_ess_fraction"], "nonfinite_count": row["nonfinite_count"],
                "peak_gpu_bytes": row.get("peak_gpu_bytes"),
                "training_seconds": training["wall_seconds_this_invocation"],
                "inference_seconds": row["wall_seconds"],
                "total_production_seconds": training["wall_seconds_this_invocation"] + row["wall_seconds"],
                "post_freeze_followup": True,
                "evaluation_paths_selected_as_largest_completed_budget": True,
                "coupled_prediction_difference_rms": refinement.get(row["steps"], {}).get("coupled_prediction_difference_rms"),
                "noise_corrected_weak_difference_rms": refinement.get(row["steps"], {}).get("noise_corrected_weak_difference_rms"),
                "observed_weak_order": refinement.get(row["steps"], {}).get("observed_order_using_next_difference"),
            })
    for path in sorted(RAW_DIR.glob("selected_precision_*.json")):
        payload = _load(path)
        refinement = {
            item["coarse_steps"]: item
            for item in payload.get("refinement_rows", [])
        }
        for row in payload["rows"]:
            rows.append({
                "family": payload["spec"]["family"], "d": payload["spec"]["d"],
                "method": "Tangent-Actor + MCV", "training_seed": payload["spec"]["seed"],
                "actor": payload["selected_actor"], "steps": row["steps"],
                "paths": row["paths"], "value_rmse": row["value_rmse"],
                "relative_l2": row["relative_l2"], "bias": row["mean_bias"],
                "rmse_ci95_low": row["rmse_ci95_low"],
                "rmse_ci95_high": row["rmse_ci95_high"],
                "mc_rmse_contribution": row["mc_rmse_contribution"],
                "reference_uncertainty_rmse": row["reference_uncertainty_rmse"],
                "residual_error_proxy": row["residual_non_mc_non_reference_rmse_proxy"],
                "mcv_variance_reduction": row["mcv_variance_reduction"],
                "ess": row["weight_ess_fraction"], "nonfinite_count": row["nonfinite_count"],
                "peak_gpu_bytes": row.get("peak_gpu_bytes"),
                "training_seconds": None, "inference_seconds": row["wall_seconds"],
                "total_production_seconds": None, "post_freeze_followup": True,
                "evaluation_paths_selected_as_largest_completed_budget": True,
                "coupled_prediction_difference_rms": refinement.get(row["steps"], {}).get("coupled_prediction_difference_rms"),
                "noise_corrected_weak_difference_rms": refinement.get(row["steps"], {}).get("noise_corrected_weak_difference_rms"),
                "observed_weak_order": None,
            })
    return rows


def _nonseparable_rows() -> list[dict[str, Any]]:
    rows = []
    for path in sorted(RESULTS_DIR.glob("post_freeze_coupled_d100*_main_mcv_s*.json")):
        payload = _load(path)
        selection = payload["selection"]
        training_path = RESULTS_DIR / f"training_{payload['spec']['family']}_d100_k2_lam0p2_eta1_main_mcv_s{payload['spec']['seed']}.json"
        training = _load(training_path)
        score = {
            item["actor"]: item["normalized_total_inconsistency"]
            for item in selection["ranking_audit_diagnostics"]
        }
        phases = {f"a{phase['iteration']}": phase for phase in training["phases"]}
        for item in payload["iteration_rows"]:
            phase = phases.get(item["iteration"])
            rows.append({
                "seed": payload["spec"]["seed"], "iteration": item["iteration"],
                "selected": item["iteration"] == selection["selected_actor"],
                "value_rmse_64": item["mcv_value_rmse"],
                "relative_l2_64": item["mcv_relative_l2"],
                "gradient_vector_rmse": item["actor_gradient_vector_rmse"],
                "value_actor_rmse": item["actor_value_rmse"],
                "active_gradient_rmse": item["actor_gradient_active_vector_rmse"],
                "inactive_gradient_rms": item["actor_gradient_inactive_vector_rms"],
                "fixed_point_score": score[item["iteration"]],
                "mcv_variance_reduction": item["mcv_variance_reduction"],
                "ess": item["weight_ess_fraction"],
                "label_value_replication_rmse": None if phase is None else phase["labels"]["value_label_replication_rmse"],
                "label_gradient_replication_rmse": None if phase is None else phase["labels"]["gradient_label_replication_vector_rmse"],
                "selected_optimizer_step": None if phase is None else phase["training"]["selected_step"],
                "train_gradient_fit_rmse": None if phase is None else phase["training"]["selected_train_gradient_vector_rmse"],
                "validation_gradient_fit_rmse": None if phase is None else phase["training"]["selected_validation_gradient_vector_rmse"],
                "phase_training_seconds": None if phase is None else phase["training"]["wall_seconds"],
            })
    return rows


def _ambient_rows() -> list[dict[str, Any]]:
    candidates: list[tuple[str, Path]] = []
    candidates.extend(("parent", path) for path in sorted(OLD_RESULTS_DIR.glob("post_freeze_primary_d*_lam0p2_eta0_*mcv_s*.json")))
    candidates.extend(("frontier", path) for path in sorted(RESULTS_DIR.glob("post_freeze_primary_d*_lam0p2_eta0_main_mcv_s*.json")))
    candidates.extend(("low_rank", path) for path in sorted(RESULTS_DIR.glob("post_freeze_primary_d*_lowrank_r*.json")))
    rows = []
    for source, path in candidates:
        payload = _load(path)
        d = int(payload["spec"]["d"])
        if d not in (100, 400):
            continue
        selection = payload["selection"]
        selected = payload["selected_row"]
        a0 = payload["iteration_rows"][0]
        run_id = selected["run_id"]
        training_root = OLD_RESULTS_DIR if source == "parent" else RESULTS_DIR
        training_path = training_root / f"training_{run_id}.json"
        training = _load(training_path)
        phases = training["phases"]
        threshold = _load(PACKAGE_DIR / "protocol.json")["ambient_400"][
            "persistent_generalization_gap_definition"
        ]
        gap_flags = []
        for phase in phases:
            fit = phase["training"]
            train_gradient = float(fit["selected_train_gradient_vector_rmse"])
            validation_gradient = float(fit["selected_validation_gradient_vector_rmse"])
            replication = float(phase["labels"]["gradient_label_replication_vector_rmse"])
            gap_flags.append(
                validation_gradient >= float(threshold["validation_over_train_gradient_rmse_ratio"]) * max(train_gradient, 1e-300)
                and validation_gradient - train_gradient >= float(threshold["minimum_excess_as_fraction_of_label_replication_rmse"]) * replication
            )
        selected_iteration = int(selection["selected_iteration"])
        selected_phase = None if selected_iteration == 0 else phases[selected_iteration - 1]
        selected_fit = None if selected_phase is None else selected_phase["training"]
        selected_replication = None if selected_phase is None else selected_phase["labels"]["gradient_label_replication_vector_rmse"]
        fit_config = training.get("fit_config", phases[0]["training"].get("config", {}))
        rows.append({
            "source": source, "run_id": run_id, "d": d,
            "seed": payload["spec"]["seed"],
            "architecture": payload.get("architecture", "full_ambient"),
            "rank": payload.get("rank"), "budget": payload["spec"]["budget"],
            "selected_actor": selection["selected_actor"],
            "a0_value_rmse": a0["mcv_value_rmse"],
            "selected_value_rmse": selected["mcv_value_rmse"],
            "value_improvement_fraction": (a0["mcv_value_rmse"] - selected["mcv_value_rmse"]) / a0["mcv_value_rmse"],
            "gradient_vector_rmse": selected["actor_gradient_vector_rmse"],
            "active_gradient_rmse": selected["actor_gradient_active_vector_rmse"],
            "inactive_gradient_rms": selected["actor_gradient_inactive_vector_rms"],
            "gradient_effective_rank": selected.get("actor_residual_gradient_effective_rank"),
            "gradient_rank_90pct": selected.get("actor_residual_gradient_rank_90pct"),
            "mcv_variance_reduction": selected["mcv_variance_reduction"],
            "ess": selected["weight_ess_fraction"],
            "fixed_point_score": next(
                item["normalized_total_inconsistency"]
                for item in selection["ranking_audit_diagnostics"]
                if item["actor"] == selection["selected_actor"]
            ),
            "training_points": fit_config.get("training_points"),
            "optimizer_steps": fit_config.get("optimizer_steps"),
            "selected_optimizer_step": None if selected_fit is None else selected_fit["selected_step"],
            "selected_train_gradient_fit_rmse": None if selected_fit is None else selected_fit["selected_train_gradient_vector_rmse"],
            "selected_validation_gradient_fit_rmse": None if selected_fit is None else selected_fit["selected_validation_gradient_vector_rmse"],
            "selected_train_value_fit_rmse": None if selected_fit is None else selected_fit["selected_train_value_rmse"],
            "selected_validation_value_fit_rmse": None if selected_fit is None else selected_fit["selected_validation_value_rmse"],
            "selected_validation_normalized_loss": None if selected_fit is None else selected_fit["best_validation_normalized_loss"],
            "selected_gradient_label_replication_rmse": selected_replication,
            "generalization_gap_phases": sum(gap_flags),
            "generalization_gap_total_phases": len(gap_flags),
            "persistent_generalization_gap": sum(gap_flags) >= int(threshold["minimum_update_phases"]),
        })
    return rows


def _compute_rows(full_history, precision, nonseparable, ambient) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(OLD_RESULTS_DIR.glob("training_primary_d100_k2_lam0p2_eta0_main_mcv_s*.json")):
        payload = _load(path)
        rows.append({"run_id": payload["spec"]["seed"], "family": "primary", "d": 100, "stage": "historical_training_and_selection", "cost_category": "production_training_and_reference_free_selection", "cost_scope": "top_level", "seconds": payload["wall_seconds_this_invocation"], "peak_gpu_bytes": None})
    for path in sorted(RESULTS_DIR.glob("training_*.json")):
        payload = _load(path)
        run_id = payload.get("run_id", payload["spec"].get("seed"))
        rows.append({"run_id": run_id, "family": payload["spec"]["family"], "d": payload["spec"]["d"], "stage": "training_and_selection", "cost_category": "production_training_and_reference_free_selection", "cost_scope": "top_level", "seconds": payload["wall_seconds_this_invocation"], "peak_gpu_bytes": None})
        for phase in payload["phases"]:
            rows.append({"run_id": run_id, "family": payload["spec"]["family"], "d": payload["spec"]["d"], "stage": f"a{phase['iteration']}_actor_fit", "cost_category": "training_component", "cost_scope": "component_nonadditive", "seconds": phase["training"]["wall_seconds"], "peak_gpu_bytes": phase["training"].get("peak_gpu_bytes")})
            rows.append({"run_id": run_id, "family": payload["spec"]["family"], "d": payload["spec"]["d"], "stage": f"a{phase['iteration']}_label_generation", "cost_category": "supervision_component", "cost_scope": "component_nonadditive", "seconds": sum(float(phase["labels"][name]["wall_seconds"]) for name in ("train", "validation", "replication")), "peak_gpu_bytes": max(phase["labels"][name].get("peak_gpu_bytes") or 0 for name in ("train", "validation", "replication"))})
    for row in precision:
        rows.append({"run_id": row["training_seed"], "family": row["family"], "d": row["d"], "stage": f"precision_inference_{row['steps']}", "cost_category": "post_freeze_precision_diagnostic", "cost_scope": "component", "seconds": row["inference_seconds"], "peak_gpu_bytes": row.get("peak_gpu_bytes")})
    for row in full_history:
        rows.append({"run_id": f"fh_{row['family']}_n{row['n']}_M{row['M']}_s{row['seed']}", "family": row["family"], "d": row["d"], "stage": "full_history_inference", "cost_category": "baseline_frontier_evaluation", "cost_scope": "top_level", "seconds": row["wall_seconds"], "peak_gpu_bytes": row["peak_gpu_bytes"]})
    for path in sorted(RESULTS_DIR.glob("post_freeze_*.json")):
        payload = _load(path)
        rows.append({
            "run_id": payload.get("run_id", payload["selected_row"]["run_id"]),
            "family": payload["spec"]["family"], "d": payload["spec"]["d"],
            "stage": "post_freeze_actor_and_policy_audit",
            "cost_category": "ground_truth_diagnostic",
            "cost_scope": "top_level", "seconds": payload["wall_seconds"],
            "peak_gpu_bytes": None,
        })
    for family in ("primary", "coupled"):
        path = REFERENCE_DIR / f"{family}_refinement_audit.json"
        if path.exists():
            payload = _load(path)
            rows.append({"run_id": f"reference_{family}_n{payload['refined_n']}", "family": family, "d": 2, "stage": "independent_reference_refinement", "cost_category": "diagnostic_cpu_reference", "cost_scope": "top_level_latest_refinement", "seconds": payload["wall_seconds"], "peak_gpu_bytes": None})
    return rows


def _plot_primary(precision: list[dict], full_history: list[dict]) -> None:
    primary = [row for row in precision if row["family"] == "primary"]
    if not primary:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2))
    for seed in sorted({row["training_seed"] for row in primary}):
        local = sorted((row for row in primary if row["training_seed"] == seed), key=lambda row: row["steps"])
        axes[0].plot([r["steps"] for r in local], [r["value_rmse"] for r in local], "o-", label=f"Actor seed {seed}")
    endpoint = [row for row in full_history if row["family"] == "primary" and row["source"] == "accuracy_frontier" and row["finite"]]
    for (n, M) in sorted({(row["n"], row["M"]) for row in endpoint}):
        values = [row["value_rmse"] for row in endpoint if row["n"] == n and row["M"] == M]
        axes[0].axhline(np.mean(values), linestyle="--", alpha=.55, label=f"FH n={n}, M={M}")
    axes[0].set(xscale="log", yscale="log", xlabel="Euler steps", ylabel="Value RMSE", title="Frozen Actor: time-step convergence")
    axes[0].set_xticks([64, 128, 256, 512], ["64", "128", "256", "512"])
    axes[0].grid(True, which="both", alpha=.25)
    axes[0].legend(fontsize=7)

    fh = [row for row in full_history if row["family"] == "primary" and row["finite"]]
    axes[1].scatter([row["wall_seconds"] for row in fh], [row["value_rmse"] for row in fh], marker="s", alpha=.6, label="Full-History")
    axes[1].scatter([row["total_production_seconds"] for row in primary], [row["value_rmse"] for row in primary], marker="o", alpha=.8, label="Tangent-Actor + MCV")
    axes[1].set(xscale="log", yscale="log", xlabel="Measured production GPU seconds", ylabel="Value RMSE", title="Accuracy–compute frontier")
    axes[1].grid(True, which="both", alpha=.25)
    axes[1].legend()
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"primary_accuracy_frontier.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)


def _plot_nonseparable(rows: list[dict]) -> None:
    if not rows:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.1))
    for seed in sorted({row["seed"] for row in rows}):
        local = sorted((row for row in rows if row["seed"] == seed), key=lambda row: int(row["iteration"][1:]))
        x = [int(row["iteration"][1:]) for row in local]
        axes[0].plot(x, [row["value_rmse_64"] for row in local], "o-", label=f"seed {seed}")
        axes[1].plot(x, [row["gradient_vector_rmse"] for row in local], "o-", label=f"seed {seed}")
    for axis, ylabel, title in zip(axes, ("Value RMSE", "Gradient-vector RMSE"), ("Nonseparable policy values", "Nonseparable Actor gradients"), strict=True):
        axis.set(xlabel="Actor iteration", ylabel=ylabel, title=title, yscale="log")
        axis.grid(True, which="both", alpha=.25)
        axis.legend()
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"nonseparable_high_budget.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)


def _plot_ambient(rows: list[dict]) -> None:
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    labels, values, colors = [], [], []
    for row in rows:
        if row["source"] == "parent" and not (row["d"] == 400 and row["budget"] == "pilot"):
            continue
        labels.append(f"d={row['d']}\n{row['architecture']}" + ("\n" + str(row["rank"]) if row["rank"] else ""))
        values.append(row["selected_value_rmse"])
        colors.append("#4472C4" if row["architecture"] == "full_ambient" else "#ED7D31")
    ax.bar(range(len(values)), values, color=colors)
    ax.set_yscale("log")
    ax.set_xticks(range(len(values)), labels)
    ax.set_ylabel("Selected Value RMSE")
    ax.set_title("Ambient-dimension generalization audit")
    ax.grid(True, axis="y", which="both", alpha=.25)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"ambient_400d_generalization.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)


def run() -> dict[str, Any]:
    ensure_dirs()
    frozen_protocol = _load(PACKAGE_DIR / "protocol.json")
    full_history = _full_history_rows()
    precision = _precision_rows()
    nonseparable = _nonseparable_rows()
    ambient = _ambient_rows()
    if not precision:
        raise RuntimeError("precision evaluation results are not complete")
    write_csv(RESULTS_DIR / "full_history_extended_grid.csv", full_history)
    write_csv(RESULTS_DIR / "euler_convergence.csv", precision)
    write_csv(RESULTS_DIR / "nonseparable_three_seed_results.csv", nonseparable)
    write_csv(RESULTS_DIR / "ambient_scaling_comparison.csv", ambient)
    accuracy = []
    for row in precision:
        accuracy.append(dict(row))
    for row in full_history:
        accuracy.append({
            "family": row["family"], "d": row["d"], "method": "Full-History MLP",
            "training_seed": row["seed"], "actor": None, "steps": None,
            "paths": None, "value_rmse": row["value_rmse"], "relative_l2": row["relative_l2"],
            "bias": row["bias"], "rmse_ci95_low": None, "rmse_ci95_high": None,
            "mc_rmse_contribution": None, "reference_uncertainty_rmse": None,
            "residual_error_proxy": None, "mcv_variance_reduction": None,
            "ess": None, "nonfinite_count": row["nonfinite_count"],
            "peak_gpu_bytes": row["peak_gpu_bytes"],
            "training_seconds": 0.0, "inference_seconds": row["wall_seconds"],
            "total_production_seconds": row["wall_seconds"], "post_freeze_followup": row["source"] == "accuracy_frontier",
            "n": row["n"], "M": row["M"], "configuration_role": row["role"],
        })
    write_csv(RESULTS_DIR / "accuracy_frontier.csv", accuracy)
    compute = _compute_rows(full_history, precision, nonseparable, ambient)
    write_csv(RESULTS_DIR / "compute_accounting.csv", compute)
    _plot_primary(precision, full_history)
    _plot_nonseparable(nonseparable)
    _plot_ambient(ambient)

    primary512 = [row for row in precision if row["family"] == "primary" and row["steps"] == 512]
    primary_actor_mean, primary_actor_sd = _mean_sd([row["value_rmse"] for row in primary512]) if primary512 else (math.nan, math.nan)
    endpoint_groups = {}
    for row in full_history:
        if row["family"] == "primary" and row["source"] == "accuracy_frontier" and row["role"] == "endpoint" and row["finite"]:
            endpoint_groups.setdefault((row["n"], row["M"]), []).append(row["value_rmse"])
    endpoint_summary = [
        {"n": key[0], "M": key[1], "seeds": len(values), "rmse_mean": _mean_sd(values)[0], "rmse_sd": _mean_sd(values)[1], "rmse_ci95_low": _mean_ci95(values)[0], "rmse_ci95_high": _mean_ci95(values)[1]}
        for key, values in sorted(endpoint_groups.items())
    ]
    strongest = min(endpoint_summary, key=lambda row: row["rmse_mean"]) if endpoint_summary else None
    coupled_endpoint_groups = {}
    for row in full_history:
        if row["family"] == "coupled" and row["source"] == "accuracy_frontier" and row["role"] == "endpoint" and row["finite"]:
            coupled_endpoint_groups.setdefault((row["n"], row["M"]), []).append(row["value_rmse"])
    coupled_endpoint_summary = [
        {"n": key[0], "M": key[1], "seeds": len(values), "rmse_mean": _mean_sd(values)[0], "rmse_sd": _mean_sd(values)[1], "rmse_ci95_low": _mean_ci95(values)[0], "rmse_ci95_high": _mean_ci95(values)[1]}
        for key, values in sorted(coupled_endpoint_groups.items())
    ]
    strongest_coupled = min(coupled_endpoint_summary, key=lambda row: row["rmse_mean"]) if coupled_endpoint_summary else None
    coupled_selected = [row for row in nonseparable if row["selected"]]
    coupled_precision = [row for row in precision if row["family"] == "coupled"]
    coupled_finest_steps = max((row["steps"] for row in coupled_precision), default=None)
    coupled_finest = [row for row in coupled_precision if row["steps"] == coupled_finest_steps]
    coupled_actor_finest_mean = None if not coupled_finest else _mean_sd([row["value_rmse"] for row in coupled_finest])[0]
    coupled_actor_finest_sd = None if not coupled_finest else _mean_sd([row["value_rmse"] for row in coupled_finest])[1]
    d400 = [row for row in ambient if row["source"] == "frontier" and row["d"] == 400]
    low_rank = [row for row in ambient if row["source"] == "low_rank"]
    primary_advantage = None if strongest is None or not primary512 else strongest["rmse_mean"] - primary_actor_mean
    primary_by_steps = []
    for steps in (64, 128, 256, 512):
        local = [row for row in precision if row["family"] == "primary" and row["steps"] == steps]
        if local:
            primary_by_steps.append({
                "steps": steps,
                "seeds": len(local),
                "rmse_mean": _mean_sd([row["value_rmse"] for row in local])[0],
                "rmse_sd": _mean_sd([row["value_rmse"] for row in local])[1],
                "mc_rmse_mean": float(np.mean([row["mc_rmse_contribution"] for row in local])),
                "residual_proxy_mean": float(np.mean([row["residual_error_proxy"] for row in local])),
                "inference_seconds_mean": float(np.mean([row["inference_seconds"] for row in local])),
                "total_production_seconds_mean": float(np.mean([row["total_production_seconds"] for row in local])),
            })
    threshold = float(frozen_protocol["reporting_accuracy_threshold_rmse"])
    actor_threshold = [row for row in primary_by_steps if row["rmse_mean"] <= threshold]
    fh_threshold = []
    for key, local in {
        (row["source"], row["n"], row["M"]): [
            item for item in full_history
            if item["family"] == "primary" and item["finite"]
            and (item["source"], item["n"], item["M"]) == (row["source"], row["n"], row["M"])
        ]
        for row in full_history if row["family"] == "primary" and row["finite"]
    }.items():
        rmse_mean = float(np.mean([item["value_rmse"] for item in local]))
        if rmse_mean <= threshold:
            fh_threshold.append({
                "source": key[0], "n": key[1], "M": key[2],
                "rmse_mean": rmse_mean,
                "seconds_mean": float(np.mean([item["wall_seconds"] for item in local])),
            })
    robust_primary = bool(
        primary_advantage is not None
        and primary_advantage > 0
        and primary_advantage > 2.0 * max(row["reference_uncertainty_rmse"] for row in primary512)
        and max(row["rmse_ci95_high"] for row in primary512)
        < strongest["rmse_ci95_low"] - max(row["reference_uncertainty_rmse"] for row in primary512)
    )
    nonseparable_support = bool(
        len(coupled_selected) >= 3
        and all(row["value_rmse_64"] < 0.00198 for row in coupled_selected)
    )
    d400_improvements = [row["value_improvement_fraction"] for row in d400]
    d400_seed0_meaningful = bool(
        any(row["seed"] == 0 and row["value_improvement_fraction"] >= 0.10 for row in d400)
    )
    d400_fixed = bool(
        len(d400) >= 3
        and np.mean(d400_improvements) >= 0.10
        and min(d400_improvements) > 0.0
    )
    if robust_primary and nonseparable_support and d400_fixed:
        decision = "GO"
    elif robust_primary and nonseparable_support:
        decision = "CONDITIONAL GO"
    elif strongest is not None and strongest["rmse_mean"] <= primary_actor_mean:
        decision = "NO-GO for the broad accuracy claim"
    else:
        decision = "CONDITIONAL GO"
    summary = {
        "decision": decision,
        "primary_actor_512_rmse_mean": primary_actor_mean,
        "primary_actor_512_rmse_sd": primary_actor_sd,
        "strongest_full_history_endpoint": strongest,
        "primary_actor_advantage_positive": None if primary_advantage is None else primary_advantage > 0,
        "primary_actor_advantage_rmse": primary_advantage,
        "robust_primary_advantage": robust_primary,
        "primary_euler_summary": primary_by_steps,
        "reporting_accuracy_threshold_rmse": threshold,
        "actor_minimum_cost_to_threshold": None if not actor_threshold else min(actor_threshold, key=lambda row: row["total_production_seconds_mean"]),
        "full_history_minimum_cost_to_threshold": None if not fh_threshold else min(fh_threshold, key=lambda row: row["seconds_mean"]),
        "nonseparable_selected_runs": len(coupled_selected),
        "nonseparable_support": nonseparable_support,
        "nonseparable_all_selected_improve_over_a0": bool(
            len(coupled_selected) >= 3
            and all(
                row["value_rmse_64"] < next(
                    initial["value_rmse_64"] for initial in nonseparable
                    if initial["seed"] == row["seed"] and initial["iteration"] == "a0"
                )
                for row in coupled_selected
            )
        ),
        "coupled_actor_finest_steps": coupled_finest_steps,
        "coupled_actor_finest_rmse_mean": coupled_actor_finest_mean,
        "coupled_actor_finest_rmse_sd": coupled_actor_finest_sd,
        "strongest_coupled_full_history_endpoint": strongest_coupled,
        "coupled_actor_beats_full_history": None if strongest_coupled is None or coupled_actor_finest_mean is None else coupled_actor_finest_mean < strongest_coupled["rmse_mean"],
        "d400_fixed": d400_fixed,
        "d400_seed0_meaningful": d400_seed0_meaningful,
        "d400_completed_seeds": len(d400),
        "d400_value_improvement_mean": None if not d400 else float(np.mean(d400_improvements)),
        "d400_value_improvement_sd": None if not d400 else _mean_sd(d400_improvements)[1],
        "low_rank_runs": len(low_rank),
        "counts": {
            "precision_rows": len(precision), "full_history_rows": len(full_history),
            "nonseparable_rows": len(nonseparable), "ambient_rows": len(ambient),
        },
    }
    atomic_json(RESULTS_DIR / "FINAL_SUMMARY.json", summary)
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
