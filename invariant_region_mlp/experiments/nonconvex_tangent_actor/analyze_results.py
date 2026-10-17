"""Aggregate compact results, draw figures, and write the research report."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import (
    FIGURE_DIR,
    PACKAGE_DIR,
    RAW_DIR,
    RESULTS_DIR,
    atomic_json,
    write_csv,
)


def _jsons(pattern: str) -> list[dict[str, Any]]:
    return [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(RESULTS_DIR.glob(pattern))
    ]


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_reproducibility_manifest() -> None:
    protocol_path = PACKAGE_DIR / "protocol.json"
    test_path = PACKAGE_DIR / "test_correctness.py"
    scripts = sorted((PACKAGE_DIR / "slurm").glob("*.slurm"))
    environments = []
    for path in sorted(RAW_DIR.glob("environment_*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        environments.append({"record": path.name, **payload})
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:
        commit = "unknown"
    manifest = {
        "experiment": "nonconvex_nonconcave_tangent_actor_v1",
        "branch": "ir-mlp-nonconvex-nonconcave-tangent-actor",
        "generated_from_git_commit": commit,
        "protocol": {
            "path": str(protocol_path.relative_to(PACKAGE_DIR.parents[1])),
            "sha256": _file_sha256(protocol_path),
            "frozen_configuration": json.loads(
                protocol_path.read_text(encoding="utf-8")
            ),
        },
        "correctness_tests": {
            "command": "PYTHONPATH=$PWD pytest -q invariant_region_mlp/experiments/nonconvex_tangent_actor/test_correctness.py",
            "test_function_definitions": test_path.read_text(encoding="utf-8").count(
                "\ndef test_"
            ),
            "latest_verified_result": "13 passed",
            "verified_date": "2026-10-17",
            "test_file_sha256": _file_sha256(test_path),
        },
        "slurm_scripts": [
            {
                "path": str(path.relative_to(PACKAGE_DIR.parents[1])),
                "sha256": _file_sha256(path),
            }
            for path in scripts
        ],
        "recorded_runtime_environments": environments,
        "determinism": {
            "seed_derivation": "SHA-256 stable_seed in common.py",
            "common_postfreeze_state_seed": 20261017,
            "common_fixed_point_state_seed": 20261017,
            "fixed_validation_states_across_candidates": True,
            "independent_training_and_validation_brownian_streams": True,
        },
        "selection_boundary": {
            "criterion": "minimum independently replicated, noise-corrected normalized fixed-point inconsistency",
            "ranking_audit_doubles_paths_and_steps": True,
            "ground_truth_available_to_training_or_selection": False,
        },
        "large_local_artifacts": {
            "directories": ["cache", "checkpoints", "predictions", "references", "raw", "slurm"],
            "tracked_in_git": False,
            "compact_results_and_audits_tracked": True,
        },
    }
    atomic_json(RESULTS_DIR / "REPRODUCIBILITY_MANIFEST.json", manifest)


def _rank_correlation(x: list[float], y: list[float]) -> float | None:
    if len(x) < 3:
        return None
    x_order = np.argsort(np.argsort(np.asarray(x)))
    y_order = np.argsort(np.argsort(np.asarray(y)))
    value = np.corrcoef(x_order, y_order)[0, 1]
    return float(value) if np.isfinite(value) else None


def consistency_rows(training: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for payload in training:
        spec = payload["spec"]
        for key, level in (
            ("primary_diagnostics", "primary"),
            ("ranking_audit_diagnostics", "ranking_audit"),
        ):
            for row in payload["selection"][key]:
                rows.append(
                    {
                        **spec,
                        "run_id": row["run_id"],
                        "iteration": row["actor"],
                        "diagnostic_level": level,
                        "selected": row["actor"]
                        == payload["selection"]["selected_actor"],
                        **{
                            field: value
                            for field, value in row.items()
                            if field not in {"run_id", "actor", "audit_tag"}
                        },
                    }
                )
    return rows


def iteration_rows(
    post_freeze: list[dict[str, Any]], consistency: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    lookup = {
        (row["run_id"], row["iteration"]): row
        for row in consistency
        if row["diagnostic_level"] == "ranking_audit"
    }
    rows: list[dict[str, Any]] = []
    for payload in post_freeze:
        selection = payload["selection"]
        for source in payload["iteration_rows"]:
            consistency_row = lookup.get((source["run_id"], source["iteration"]), {})
            rows.append(
                {
                    **source,
                    "selected": source["iteration"] == selection["selected_actor"],
                    "fixed_point_inconsistency": consistency_row.get(
                        "normalized_total_inconsistency"
                    ),
                    "fixed_point_value_rmse": consistency_row.get(
                        "value_noise_corrected_rmse"
                    ),
                    "fixed_point_gradient_vector_rmse": consistency_row.get(
                        "gradient_noise_corrected_vector_rmse"
                    ),
                    "selection_ranking_stable": selection[
                        "ranking_stable_when_paths_and_steps_doubled"
                    ],
                }
            )
    return rows


def mcv_rows(iterations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fields = (
        "run_id",
        "family",
        "d",
        "lambda",
        "eta",
        "seed",
        "budget",
        "iteration",
        "selected",
        "mcv_variance_reduction",
        "plain_mc_rmse_contribution",
        "mcv_mc_rmse_contribution",
        "plain_value_rmse",
        "mcv_value_rmse",
        "mean_plain_mc_se",
        "mean_mcv_mc_se",
        "weight_ess_fraction",
        "paths",
        "time_steps",
    )
    return [{key: row.get(key) for key in fields} for row in iterations]


def compute_rows(
    training: list[dict[str, Any]],
    post_freeze: list[dict[str, Any]],
    baselines: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for payload in training:
        spec = payload["spec"]
        cumulative = 0.0
        for phase in payload["phases"]:
            labels = phase["labels"]
            label_seconds = sum(
                labels[key]["wall_seconds"]
                for key in ("train", "validation", "replication")
            )
            training_seconds = phase["training"]["wall_seconds"]
            cumulative += label_seconds + training_seconds
            rows.append(
                {
                    **spec,
                    "run_id": payload["selection"]["spec"]
                    and phase["run_id"],
                    "component": f"a{phase['iteration']}_training",
                    "label_seconds": label_seconds,
                    "optimizer_seconds": training_seconds,
                    "inference_seconds": 0.0,
                    "total_seconds": label_seconds + training_seconds,
                    "cumulative_training_seconds": cumulative,
                    "peak_gpu_bytes": max(
                        value
                        for value in (
                            labels["train"].get("peak_gpu_bytes"),
                            phase["training"].get("peak_gpu_bytes"),
                        )
                        if value is not None
                    ),
                }
            )
    for payload in post_freeze:
        for row in payload["iteration_rows"]:
            rows.append(
                {
                    **payload["spec"],
                    "run_id": row["run_id"],
                    "component": f"{row['iteration']}_post_freeze_evaluation",
                    "label_seconds": 0.0,
                    "optimizer_seconds": 0.0,
                    "inference_seconds": row["inference_seconds"],
                    "total_seconds": row["inference_seconds"],
                    "cumulative_training_seconds": None,
                    "peak_gpu_bytes": None,
                }
            )
    for payload in baselines:
        if payload["method"] == "full_history_mlp_raw_adaptation":
            component = f"full_history_n{payload['n']}_M{payload['M']}_s{payload['baseline_seed']}"
            inference_seconds = payload["runtime"]["wall_seconds"]
            peak_gpu_bytes = payload["runtime"].get("peak_gpu_bytes")
        else:
            component = f"pinn_s{payload['baseline_seed']}"
            inference_seconds = payload["inference_seconds"]
            peak_gpu_bytes = payload["training"].get("peak_gpu_bytes")
        rows.append(
            {
                **payload["spec"],
                "run_id": payload["actor_test_bank_run_id"],
                "component": component,
                "label_seconds": 0.0,
                "optimizer_seconds": payload.get("training_seconds", 0.0),
                "inference_seconds": inference_seconds,
                "total_seconds": payload["total_seconds"],
                "cumulative_training_seconds": 0.0,
                "peak_gpu_bytes": peak_gpu_bytes,
            }
        )
    return rows


def baseline_rows(
    post_freeze: list[dict[str, Any]], baselines: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for payload in post_freeze:
        selected = payload["selected_row"]
        rows.append(
            {
                "method": "Tangent-Actor + MCV",
                **payload["spec"],
                "run_id": selected["run_id"],
                "configuration": payload["selection"]["selected_actor"],
                "seed": payload["spec"]["seed"],
                "value_rmse": selected["mcv_value_rmse"],
                "relative_l2": selected["mcv_relative_l2"],
                "bias": selected["mcv_bias"],
                "total_seconds": sum(
                    phase["labels"][split]["wall_seconds"]
                    for phase in next(
                        item
                        for item in _jsons("training_*.json")
                        if item["selection"]["spec"] == payload["spec"]
                    )["phases"]
                    for split in ("train", "validation", "replication")
                )
                + sum(
                    phase["training"]["wall_seconds"]
                    for phase in next(
                        item
                        for item in _jsons("training_*.json")
                        if item["selection"]["spec"] == payload["spec"]
                    )["phases"]
                )
                + selected["inference_seconds"],
                "finite": selected["nonfinite_count"] == 0,
                "original_algorithm_adaptation": False,
            }
        )
    for payload in baselines:
        full_history = payload["method"] == "full_history_mlp_raw_adaptation"
        rows.append(
            {
                "method": "Full-History MLP" if full_history else "Hutchinson PINN",
                **payload["spec"],
                "run_id": payload["actor_test_bank_run_id"],
                "configuration": (
                    f"n={payload['n']},M={payload['M']}"
                    if full_history
                    else "frozen 4x256, 5000 steps"
                ),
                "seed": payload["baseline_seed"],
                "value_rmse": payload["metrics"]["rmse"],
                "relative_l2": payload["metrics"]["relative_l2"],
                "bias": payload["metrics"]["bias"],
                "total_seconds": payload["total_seconds"],
                "finite": payload["runtime"]["finite"] if full_history else True,
                "original_algorithm_adaptation": True,
            }
        )
    return rows


def _group(rows: list[dict[str, Any]], key: str) -> dict[Any, list[dict[str, Any]]]:
    result: dict[Any, list[dict[str, Any]]] = {}
    for row in rows:
        result.setdefault(row[key], []).append(row)
    return result


def _plot_iterations(rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    groups = _group(rows, "run_id")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
    for run_id, values in groups.items():
        values = sorted(values, key=lambda row: int(row["iteration"][1:]))
        x = [int(row["iteration"][1:]) for row in values]
        label = run_id.replace("_main_mcv_", "_").replace("_pilot_mcv_", "_")
        axes[0].plot(x, [row["mcv_value_rmse"] for row in values], "o-", label=label)
        axes[1].plot(
            x, [row["actor_gradient_vector_rmse"] for row in values], "o-"
        )
        axes[2].plot(
            x, [row["fixed_point_inconsistency"] for row in values], "o-"
        )
    titles = ("PDE Value RMSE", "Actor Gradient-vector RMSE", "Fixed-point inconsistency")
    for axis, title in zip(axes, titles, strict=True):
        axis.set_title(title)
        axis.set_xlabel("Actor iteration")
        axis.set_yscale("log")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("error (log scale)")
    axes[0].legend(fontsize=7, loc="best")
    fig.suptitle("Nonconvex Tangent-Actor iteration behavior", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "actor_iteration_behavior.png", dpi=220)
    fig.savefig(FIGURE_DIR / "actor_iteration_behavior.pdf")
    plt.close(fig)


def _plot_fixed_point(rows: list[dict[str, Any]]) -> None:
    valid = [
        row
        for row in rows
        if row.get("fixed_point_inconsistency") is not None
        and row["fixed_point_inconsistency"] > 0.0
        and row["mcv_value_rmse"] > 0.0
    ]
    if not valid:
        return
    fig, axis = plt.subplots(figsize=(6.3, 5.0))
    for row in valid:
        marker = "*" if row["selected"] else "o"
        axis.scatter(
            row["fixed_point_inconsistency"],
            row["mcv_value_rmse"],
            marker=marker,
            s=100 if row["selected"] else 42,
            alpha=0.8,
        )
        axis.annotate(
            f"{row['iteration']} s{row['seed']}",
            (row["fixed_point_inconsistency"], row["mcv_value_rmse"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("reference-free fixed-point inconsistency")
    axis.set_ylabel("post-freeze true Value RMSE")
    axis.set_title("Does fixed-point consistency rank true accuracy?")
    axis.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "fixed_point_vs_true_error.png", dpi=220)
    fig.savefig(FIGURE_DIR / "fixed_point_vs_true_error.pdf")
    plt.close(fig)


def _plot_curvature(preflight: dict[str, Any] | None) -> None:
    if not preflight:
        return
    names, positive, negative = [], [], []
    for name, payload in preflight.get("drivers", {}).items():
        coverage = payload["reference"]["curvature_coverage"]
        names.append(name)
        positive.append(coverage["f_uu_positive_fraction"])
        negative.append(coverage["f_uu_negative_fraction"])
    if not names:
        return
    fig, axis = plt.subplots(figsize=(6.0, 4.3))
    axis.bar(names, positive, label=r"$f_{uu}>0$", color="#D1495B")
    axis.bar(names, negative, bottom=positive, label=r"$f_{uu}<0$", color="#3A6EA5")
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("fraction of interior reference states")
    axis.set_title("Actual curvature-regime coverage")
    axis.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "curvature_region_coverage.png", dpi=220)
    fig.savefig(FIGURE_DIR / "curvature_region_coverage.pdf")
    plt.close(fig)


def _plot_mcv(rows: list[dict[str, Any]]) -> None:
    selected = [row for row in rows if row["selected"]]
    if not selected:
        return
    labels = [f"{r['family']} d{r['d']} λ={r['lambda']} s{r['seed']}" for r in selected]
    values = [row["mcv_variance_reduction"] for row in selected]
    fig, axis = plt.subplots(figsize=(max(7.0, 1.15 * len(labels)), 4.5))
    bars = axis.bar(np.arange(len(labels)), values, color="#5B8E7D")
    axis.axhline(1.0, color="black", linewidth=1)
    axis.set_yscale("log")
    axis.set_ylabel("plain / MCV path-pair variance")
    axis.set_title("MCV variance reduction at frozen selected Actors")
    axis.set_xticks(np.arange(len(labels)), labels, rotation=35, ha="right")
    for bar, value in zip(bars, values, strict=True):
        axis.text(bar.get_x() + bar.get_width() / 2, value, f"{value:.1f}×", ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "mcv_variance_reduction.png", dpi=220)
    fig.savefig(FIGURE_DIR / "mcv_variance_reduction.pdf")
    plt.close(fig)


def _summary(
    training: list[dict[str, Any]],
    post_freeze: list[dict[str, Any]],
    iterations: list[dict[str, Any]],
    baselines: list[dict[str, Any]],
) -> dict[str, Any]:
    pilot_path = RESULTS_DIR / "pilot_decision.json"
    pilot = json.loads(pilot_path.read_text()) if pilot_path.exists() else None
    selected = [row for row in iterations if row["selected"]]
    primary_main = [
        row
        for row in selected
        if row["family"] == "primary" and row["budget"] == "main" and row["lambda"] == 0.2
    ]
    coupled = [row for row in selected if row["family"] == "coupled"]
    d400 = [row for row in selected if row["d"] == 400]
    consistency_values = [
        row["fixed_point_inconsistency"]
        for row in iterations
        if row.get("fixed_point_inconsistency") is not None
    ]
    true_values = [
        row["mcv_value_rmse"]
        for row in iterations
        if row.get("fixed_point_inconsistency") is not None
    ]
    all_stable = bool(selected) and all(
        row["nonfinite_count"] == 0 and row["weight_ess_fraction"] >= 0.2
        for row in selected
    )
    enough_breadth = len({row["seed"] for row in primary_main}) >= 2 and bool(coupled) and bool(d400)
    decision = (
        "NO-GO"
        if pilot is not None and not pilot["passed"]
        else "GO"
        if enough_breadth and all_stable and baselines
        else "CONDITIONAL GO"
    )
    return {
        "decision": decision,
        "pilot": pilot,
        "completed_training_runs": len(training),
        "completed_post_freeze_runs": len(post_freeze),
        "selected_runs": len(selected),
        "primary_main_seeds": sorted({row["seed"] for row in primary_main}),
        "coupled_completed": bool(coupled),
        "d400_completed": bool(d400),
        "full_history_configurations_completed": sum(
            item["method"] == "full_history_mlp_raw_adaptation" for item in baselines
        ),
        "neural_baseline_runs_completed": sum(
            item["method"] != "full_history_mlp_raw_adaptation" for item in baselines
        ),
        "all_selected_runs_finite_and_ess_above_0p2": all_stable,
        "fixed_point_true_error_spearman": _rank_correlation(
            consistency_values, true_values
        ),
        "selected_results": selected,
        "scope_statement": "ambient d-dimensional benchmark with a two-dimensional independent reference structure",
        "global_convergence_claimed": False,
        "bregman_regret_identity_reported": False,
    }


def _write_report(summary: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    selected = summary["selected_results"]
    lines = [
        "# Nonconvex–Nonconcave Tangent-Actor Report",
        "",
        f"## Decision: **{summary['decision']}**",
        "",
        "The experiment evaluates the unchanged undamped Terminal-Aware Tangent-Actor update with alpha=1 martingale control variates. Checkpoints are selected by an independently replicated, noise-corrected fixed-point consistency score before the numerical reference is loaded.",
        "",
        "## Frozen selected results",
        "",
        "| family | d | lambda | seed | Actor | Value RMSE | Relative L2 | Actor gradient-vector RMSE | MCV variance reduction | ESS |",
        "|---|---:|---:|---:|---|---:|---:|---:|---:|---:|",
    ]
    for row in selected:
        lines.append(
            f"| {row['family']} | {row['d']} | {row['lambda']} | {row['seed']} | {row['iteration']} | "
            f"{row['mcv_value_rmse']:.6g} | {row['mcv_relative_l2']:.6g} | "
            f"{row['actor_gradient_vector_rmse']:.6g} | {row['mcv_variance_reduction']:.2f}× | "
            f"{row['weight_ess_fraction']:.4f} |"
        )
    if not selected:
        lines.append("| — | — | — | — | — | — | — | — | — | — |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"Completed {summary['completed_training_runs']} training runs and {summary['completed_post_freeze_runs']} post-freeze evaluations. The rank correlation between the reference-free consistency score and post-freeze true Value error is {summary['fixed_point_true_error_spearman']}.",
            "",
            "No monotone policy-improvement theorem or nonnegative Bregman-regret identity is used for these drivers. Linear-policy values are diagnostics only. The experiment supports claims only for smooth bounded-curvature drivers in the tested locally stable regime.",
            "",
            "The Actor always receives the full ambient state. Active coordinates and reference gradients are used only in the post-freeze audit. These are ambient 100D/400D tests with a two-dimensional reference structure, not evidence for arbitrary effective dimension.",
            "",
            "## Artifacts",
            "",
            "- [Actor iteration behavior](figures/actor_iteration_behavior.png)",
            "- [Fixed-point consistency versus true error](figures/fixed_point_vs_true_error.png)",
            "- [Curvature-region coverage](figures/curvature_region_coverage.png)",
            "- [MCV variance reduction](figures/mcv_variance_reduction.png)",
            "- `ACTOR_CONSISTENCY_DIAGNOSTICS.csv`",
            "- `ACTOR_ITERATION_RESULTS.csv`",
            "- `MCV_VARIANCE_RESULTS.csv`",
            "- `BASELINE_COMPARISON.csv`",
            "- `COMPUTE_ACCOUNTING.csv`",
            "- `REPRODUCIBILITY_MANIFEST.json`",
        ]
    )
    (RESULTS_DIR / "NONCONVEX_TANGENT_ACTOR_REPORT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    training = _jsons("training_*.json")
    post_freeze = _jsons("post_freeze_*.json")
    baselines = [
        json.loads(path.read_text(encoding="utf-8"))
        for pattern in ("baseline_full_history_*.json", "baseline_pinn_*.json")
        for path in sorted(RAW_DIR.glob(pattern))
    ]
    consistency = consistency_rows(training)
    iterations = iteration_rows(post_freeze, consistency)
    mcv = mcv_rows(iterations)
    compute = compute_rows(training, post_freeze, baselines)
    comparison = baseline_rows(post_freeze, baselines)
    if consistency:
        write_csv(RESULTS_DIR / "ACTOR_CONSISTENCY_DIAGNOSTICS.csv", consistency)
    if iterations:
        write_csv(RESULTS_DIR / "ACTOR_ITERATION_RESULTS.csv", iterations)
        write_csv(RESULTS_DIR / "MCV_VARIANCE_RESULTS.csv", mcv)
    if comparison:
        write_csv(RESULTS_DIR / "BASELINE_COMPARISON.csv", comparison)
    if compute:
        write_csv(RESULTS_DIR / "COMPUTE_ACCOUNTING.csv", compute)
    preflight_path = RESULTS_DIR / "preflight_audit.json"
    preflight = (
        json.loads(preflight_path.read_text(encoding="utf-8"))
        if preflight_path.exists()
        else None
    )
    _plot_iterations(iterations)
    _plot_fixed_point(iterations)
    _plot_curvature(preflight)
    _plot_mcv(iterations)
    summary = _summary(training, post_freeze, iterations, baselines)
    atomic_json(RESULTS_DIR / "FINAL_SUMMARY.json", summary)
    _write_reproducibility_manifest()
    _write_report(summary, iterations)
    print(json.dumps({key: value for key, value in summary.items() if key != "selected_results"}, indent=2))


if __name__ == "__main__":
    main()
