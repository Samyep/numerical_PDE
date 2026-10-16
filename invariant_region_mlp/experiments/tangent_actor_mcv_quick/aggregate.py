"""Aggregate MCV replications, audit metrics, plot, and render the report."""

from __future__ import annotations

from collections import defaultdict
import csv
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
    BENCHMARKS,
    FIGURE_DIR,
    PACKAGE_DIR,
    RAW_DIR,
    RESULTS_DIR,
    atomic_json,
    ensure_dirs,
    metrics,
    write_csv,
)


BENCHMARK_LABELS = {
    "nonseparable_d100_k2": "Nonseparable d=100, k=2",
    "joint_concave_d400_k2": "Joint concave d=400, k=2",
}
VARIANT_LABELS = {
    "no_mcv": "No MCV",
    "direct": r"Direct $\alpha=1$",
    "calibrated": r"Calibrated $\alpha$",
}
VARIANT_COLORS = {
    "no_mcv": "#7A7A7A",
    "direct": "#006BA4",
    "calibrated": "#FF800E",
}


def _json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def _mean(rows: list[dict[str, Any]], key: str) -> float:
    values = [float(row[key]) for row in rows if row.get(key) is not None]
    return float(np.mean(values)) if values else math.nan


def _mean_std(rows: list[dict[str, Any]], key: str, digits: int = 4) -> str:
    values = np.asarray(
        [float(row[key]) for row in rows if row.get(key) is not None], dtype=np.float64
    )
    if not len(values):
        return "—"
    if len(values) == 1:
        return f"{values[0]:.{digits}g}"
    return f"{values.mean():.{digits}g} ± {values.std(ddof=1):.{digits}g}"


def load_records() -> tuple[
    list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]
]:
    fragments = sorted(RAW_DIR.glob("mcv_*.json"))
    if not fragments:
        raise FileNotFoundError(f"no raw MCV fragments in {RAW_DIR}")
    variance: dict[tuple[Any, ...], dict[str, Any]] = {}
    runtime: dict[tuple[Any, ...], dict[str, Any]] = {}
    calibration: dict[tuple[Any, ...], dict[str, Any]] = {}
    runs: list[dict[str, Any]] = []
    for path in fragments:
        run = _json(path)
        runs.append({**run, "raw_path": str(path)})
        for row in run["variance_rows"]:
            key = (
                row["benchmark"], int(row["seed"]), row["checkpoint_key"],
                int(row["paths"]), row["variant"],
            )
            if key in variance:
                raise ValueError(f"duplicate variance row {key}")
            variance[key] = row
        for row in run["runtime_rows"]:
            key = (
                row["benchmark"], int(row["seed"]), row["checkpoint_key"],
                int(row["paths"]),
            )
            if key in runtime:
                raise ValueError(f"duplicate runtime row {key}")
            runtime[key] = row
        for row in run["calibration"]:
            key = (row["benchmark"], int(row["seed"]), row["checkpoint_key"])
            if key in calibration:
                raise ValueError(f"duplicate calibration row {key}")
            calibration[key] = row
    return list(variance.values()), list(runtime.values()), list(calibration.values()), runs


def _coverage_audit(variance: list[dict[str, Any]]) -> dict[str, Any]:
    expected = {
        (benchmark, seed, checkpoint, paths, variant)
        for benchmark in BENCHMARKS
        for seed in (0, 1)
        for checkpoint in ("a1", "a3", "late")
        for paths in (64, 128, 512)
        for variant in ("no_mcv", "direct", "calibrated")
    }
    observed = {
        (
            row["benchmark"], int(row["seed"]), row["checkpoint_key"],
            int(row["paths"]), row["variant"],
        )
        for row in variance
    }
    return {
        "expected_rows": len(expected),
        "observed_rows": len(observed),
        "missing": [list(item) for item in sorted(expected - observed)],
        "unexpected": [list(item) for item in sorted(observed - expected)],
        "passed": observed == expected,
    }


def _metric_audit(
    variance: list[dict[str, Any]], runs: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    run_by_row: dict[tuple[str, int, str], dict[str, Any]] = {}
    for run in runs:
        for checkpoint in run["checkpoint_keys"]:
            key = (run["benchmark"], int(run["seed"]), checkpoint)
            run_by_row[key] = run
    audit_rows = []
    for row in variance:
        run = run_by_row[(row["benchmark"], int(row["seed"]), row["checkpoint_key"])]
        with np.load(run["prediction_path"]) as archive:
            predictions = archive[row["prediction_key"]]
            reference = archive["reference"]
        recomputed_mean = metrics(predictions.mean(axis=0), reference)
        replicate_rmse = float(np.mean([metrics(item, reference)["rmse"] for item in predictions]))
        differences = {
            "rmse_of_replicate_mean": abs(
                recomputed_mean["rmse"] - float(row["rmse_of_replicate_mean"])
            ),
            "relative_l2_of_replicate_mean": abs(
                recomputed_mean["relative_l2"]
                - float(row["relative_l2_of_replicate_mean"])
            ),
            "bias_of_replicate_mean": abs(
                recomputed_mean["bias"] - float(row["bias_of_replicate_mean"])
            ),
            "mean_replicate_rmse": abs(
                replicate_rmse - float(row["mean_replicate_rmse"])
            ),
        }
        audit_rows.append(
            {
                "benchmark": row["benchmark"],
                "seed": row["seed"],
                "checkpoint_key": row["checkpoint_key"],
                "paths": row["paths"],
                "variant": row["variant"],
                **{f"abs_difference_{key}": value for key, value in differences.items()},
                "passed": max(differences.values()) < 1e-12,
            }
        )
    maximum = max(
        value
        for row in audit_rows
        for key, value in row.items()
        if key.startswith("abs_difference_")
    )
    summary = {
        "rows_checked": len(audit_rows),
        "rows_passed": sum(bool(row["passed"]) for row in audit_rows),
        "maximum_absolute_difference": maximum,
    }
    summary["passed"] = summary["rows_checked"] == summary["rows_passed"]
    return audit_rows, summary


def _correctness_audit(
    variance: list[dict[str, Any]], runtime: list[dict[str, Any]], runs: list[dict[str, Any]]
) -> dict[str, Any]:
    fingerprints: dict[str, set[str]] = defaultdict(set)
    for run in runs:
        fingerprints[run["benchmark"]].add(run["state_bank_fingerprint"])
    martingale_rows = [row for row in variance if row["variant"] == "direct"]
    mean_consistency_rows = [row for row in variance if row["variant"] != "no_mcv"]
    return {
        "all_runs_inference_only": all(not run["training_performed"] for run in runs),
        "original_solver_unmodified": all(not run["ordinary_solver_modified"] for run in runs),
        "reference_opened_after_evaluation": all(
            run["reference_opened_after_all_stochastic_evaluation"] for run in runs
        ),
        "ground_truth_used_for_alpha_or_selection": any(
            run["ground_truth_used_for_alpha_or_selection"] for run in runs
        ),
        "fixed_state_bank_consistent_across_seeds_and_fragments": all(
            len(values) == 1 for values in fingerprints.values()
        ) and set(fingerprints) == set(BENCHMARKS),
        "ordinary_fk_already_uses_mcv": any(
            row["ordinary_fk_already_uses_mcv"] for row in runtime
        ),
        "ordinary_fk_match_max_abs": max(
            float(row["ordinary_fk_match_max_abs"]) for row in runtime
        ),
        "martingale_zero_in_ci_fraction": float(np.mean([
            bool(row["martingale_ci95_contains_zero"]) for row in martingale_rows
        ])),
        "mcv_plain_mean_agreement_ci_fraction": float(np.mean([
            bool(row["prediction_difference_ci95_contains_zero"])
            for row in mean_consistency_rows
        ])),
        "nonfinite_rows": sum(int(row.get("nonfinite_count", 0)) for row in variance),
        "passed": (
            all(not run["training_performed"] for run in runs)
            and all(not run["ordinary_solver_modified"] for run in runs)
            and all(run["reference_opened_after_all_stochastic_evaluation"] for run in runs)
            and not any(run["ground_truth_used_for_alpha_or_selection"] for run in runs)
            and all(len(values) == 1 for values in fingerprints.values())
            and set(fingerprints) == set(BENCHMARKS)
            and not any(row["ordinary_fk_already_uses_mcv"] for row in runtime)
            and max(float(row["ordinary_fk_match_max_abs"]) for row in runtime) < 2e-6
        ),
    }


def _plot_variance(variance: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), sharey=True)
    positions = {"a1": 0, "a3": 1, "late": 2}
    for axis, benchmark in zip(axes, BENCHMARKS):
        for variant in ("direct", "calibrated"):
            means, lows, highs = [], [], []
            for checkpoint in ("a1", "a3", "late"):
                local = [
                    row for row in variance
                    if row["benchmark"] == benchmark
                    and row["checkpoint_key"] == checkpoint
                    and int(row["paths"]) == 128
                    and row["variant"] == variant
                ]
                values = np.asarray(
                    [float(row["variance_ratio_replicates"]) for row in local]
                )
                means.append(float(values.mean()))
                lows.append(float(values.min()))
                highs.append(float(values.max()))
            x = np.arange(3)
            axis.plot(
                x, means, marker="o", linewidth=2, color=VARIANT_COLORS[variant],
                label=VARIANT_LABELS[variant],
            )
            axis.fill_between(x, lows, highs, color=VARIANT_COLORS[variant], alpha=0.16)
        axis.axhline(1.0, color="black", linestyle="--", linewidth=1)
        axis.set_xticks(list(positions.values()), ["a1", "a3", "late"])
        axis.set_yscale("log")
        axis.set_xlabel("Frozen Actor checkpoint")
        axis.set_title(BENCHMARK_LABELS[benchmark])
        axis.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("Variance ratio to ordinary FK (M=128)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.98))
    fig.suptitle("Martingale control-variate effectiveness across Actor iterations", y=1.07)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"mcv_variance_ratio_by_iteration.{suffix}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_accuracy_compute(variance: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), sharey=False)
    for axis, benchmark in zip(axes, BENCHMARKS):
        for variant in ("no_mcv", "direct", "calibrated"):
            points = []
            for paths in (64, 128, 512):
                local = [
                    row for row in variance
                    if row["benchmark"] == benchmark
                    and row["checkpoint_key"] == "late"
                    and int(row["paths"]) == paths
                    and row["variant"] == variant
                ]
                points.append((
                    _mean(local, "amortized_seconds_per_replicate"),
                    _mean(local, "mean_replicate_rmse"),
                ))
            axis.plot(
                [item[0] for item in points], [item[1] for item in points],
                marker="o", linewidth=2, color=VARIANT_COLORS[variant],
                label=VARIANT_LABELS[variant],
            )
            for paths, point in zip((64, 128, 512), points):
                axis.annotate(str(paths), point, xytext=(4, 3), textcoords="offset points", fontsize=8)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Inference seconds for 128 states (calibration amortized)")
        axis.set_ylabel("Mean single-replicate Value RMSE")
        axis.set_title(BENCHMARK_LABELS[benchmark])
        axis.grid(True, which="both", alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.98))
    fig.suptitle("Accuracy versus inference compute for the latest Actor", y=1.07)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"mcv_accuracy_vs_inference_compute.{suffix}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def _variance_table(variance: list[dict[str, Any]]) -> str:
    lines = [
        "| Benchmark | Actor | M | Direct R | Calibrated R | Direct reduction | Calibrated reduction |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for benchmark in BENCHMARKS:
        for checkpoint in ("a1", "a3", "late"):
            for paths in (64, 128, 512):
                base = [
                    row for row in variance
                    if row["benchmark"] == benchmark
                    and row["checkpoint_key"] == checkpoint
                    and int(row["paths"]) == paths
                ]
                direct = [row for row in base if row["variant"] == "direct"]
                calibrated = [row for row in base if row["variant"] == "calibrated"]
                lines.append(
                    f"| {BENCHMARK_LABELS[benchmark]} | {checkpoint} | {paths} | "
                    f"{_mean_std(direct, 'variance_ratio_replicates', 4)} | "
                    f"{_mean_std(calibrated, 'variance_ratio_replicates', 4)} | "
                    f"{_mean_std(direct, 'effective_variance_reduction', 4)}× | "
                    f"{_mean_std(calibrated, 'effective_variance_reduction', 4)}× |"
                )
    return "\n".join(lines)


def _runtime_table(runtime: list[dict[str, Any]]) -> str:
    lines = [
        "| Benchmark | Actor | M | Plain s | MCV s | Overhead | Extra peak MiB |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for benchmark in BENCHMARKS:
        for checkpoint in ("a1", "a3", "late"):
            for paths in (128, 512):
                local = [
                    row for row in runtime
                    if row["benchmark"] == benchmark
                    and row["checkpoint_key"] == checkpoint
                    and int(row["paths"]) == paths
                ]
                lines.append(
                    f"| {BENCHMARK_LABELS[benchmark]} | {checkpoint} | {paths} | "
                    f"{_mean_std(local, 'plain_mean_seconds', 4)} | "
                    f"{_mean_std(local, 'mcv_mean_seconds', 4)} | "
                    f"{100*_mean(local, 'mcv_runtime_overhead_fraction'):.2f}% | "
                    f"{_mean(local, 'additional_peak_gpu_bytes') / 2**20:.3g} |"
                )
    return "\n".join(lines)


def _render_report(
    variance: list[dict[str, Any]], runtime: list[dict[str, Any]],
    calibration: list[dict[str, Any]], coverage: dict[str, Any],
    metric_audit: dict[str, Any], correctness: dict[str, Any],
) -> None:
    manifest_path = RESULTS_DIR / "execution_manifest.json"
    manifest = _json(manifest_path) if manifest_path.exists() else {"jobs": []}
    completed_jobs = [
        job for job in manifest.get("jobs", []) if job.get("state") == "COMPLETED"
    ]
    total_gpu_hours = sum(
        float(job.get("allocated_gpu_hours", 0.0)) for job in completed_jobs
    )
    late = [row for row in variance if row["checkpoint_key"] == "late"]
    reproducible_by_benchmark: dict[str, bool] = {}
    best_variant_by_benchmark: dict[str, str] = {}
    for benchmark in BENCHMARKS:
        candidates = {}
        for variant in ("direct", "calibrated"):
            local = [
                row for row in late
                if row["benchmark"] == benchmark and row["variant"] == variant
            ]
            candidates[variant] = _mean(local, "variance_ratio_replicates")
        best = min(candidates, key=candidates.get)
        best_variant_by_benchmark[benchmark] = best
        relevant = [
            row for row in late
            if row["benchmark"] == benchmark and row["variant"] == best
        ]
        reproducible_by_benchmark[benchmark] = (
            _mean(relevant, "variance_ratio_replicates") < 0.8
            and sum(float(row["variance_ratio_ci95_high"]) < 1.0 for row in relevant)
            >= math.ceil(0.75 * len(relevant))
        )
    consistency = correctness["mcv_plain_mean_agreement_ci_fraction"] >= 0.9
    best_late = [
        row for row in late
        if row["variant"] == best_variant_by_benchmark[row["benchmark"]]
    ]
    efficiency = _mean(best_late, "effective_variance_reduction_per_gpu_second")
    decision = (
        "GO" if all(reproducible_by_benchmark.values()) and consistency and efficiency > 1.1
        else "NO-GO"
    )
    direct_late = [row for row in late if row["variant"] == "direct"]
    calibrated_late = [row for row in late if row["variant"] == "calibrated"]
    direct_ratio = _mean(direct_late, "variance_ratio_replicates")
    calibrated_ratio = _mean(calibrated_late, "variance_ratio_replicates")
    direct_suffices = direct_ratio <= 1.1 * calibrated_ratio and direct_ratio < 0.8
    early = [
        row for row in variance
        if row["checkpoint_key"] == "a1" and row["variant"] == "direct"
    ]
    helps_more_late = direct_ratio < _mean(early, "variance_ratio_replicates")
    chosen = direct_late if direct_suffices else calibrated_late
    path_fraction = _mean(chosen, "path_saving_fraction")
    path_factor = _mean(chosen, "effective_variance_reduction")
    overhead = _mean(
        [row for row in runtime if row["checkpoint_key"] == "late"],
        "mcv_runtime_overhead_fraction",
    )
    late_diagnostic_rows = [
        row for row in late if row["variant"] == "direct"
    ]
    payoff_antithetic_ratio = _mean(
        late_diagnostic_rows, "payoff_antithetic_variance_ratio"
    )
    martingale_antithetic_ratio = _mean(
        late_diagnostic_rows, "martingale_antithetic_variance_ratio"
    )
    ess_min = min(float(row["weight_ess_fraction"]) for row in variance)
    ess_max = max(float(row["weight_ess_fraction"]) for row in variance)
    alpha_table = [
        "| Benchmark | Actor | Calibrated alpha | Calibration s | Pair corr(P,H) |",
        "|---|---|---:|---:|---:|",
    ]
    for benchmark in BENCHMARKS:
        for checkpoint in ("a1", "a3", "late"):
            local = [
                row for row in calibration
                if row["benchmark"] == benchmark and row["checkpoint_key"] == checkpoint
            ]
            alpha_table.append(
                f"| {BENCHMARK_LABELS[benchmark]} | {checkpoint} | "
                f"{_mean_std(local, 'alpha', 4)} | {_mean_std(local, 'calibration_seconds', 4)} | "
                f"{_mean_std(local, 'payoff_martingale_pair_correlation', 4)} |"
            )
    lines = [
        "# Rapid Martingale Control-Variate Validation for Tangent-Actor",
        "",
        f"**Decision: {decision}.** The decision uses independently calibrated coefficients, eight fresh MC replications, paired Brownian streams and antithetic pair averages as the variance units. Ground Truth is used only for the posthoc accuracy columns.",
        "",
        "## Frozen protocol",
        "",
        r"No Actor was retrained and the original FK solver was not modified. For every Euler step the auxiliary estimator accumulates $H\leftarrow H+C_mq(t_m,X_m)^\top\Delta W_m$ at the left endpoint before the state and exponential weight are advanced. The ordinary and MCV predictions use the same Brownian paths. A separate 512-path stream estimates one scalar $\alpha$ per Actor after within-state centering.",
        "",
        "The ordinary implementation reports `martingale_control_variate=False`; identical-seed ordinary values and the $\alpha=0$ augmented values agree to the tolerance recorded in the correctness audit.",
        "",
        "## Variance ratio across iterations",
        "",
        "The 95% intervals below are paired state bootstraps of the replicate-variance ratio; all variants share Brownian streams within each replicate.",
        "",
        _variance_table(variance),
        "",
        "![Variance ratio versus Actor iteration](figures/mcv_variance_ratio_by_iteration.png)",
        "",
        "## Calibration",
        "",
        "\n".join(alpha_table),
        "",
        "## Runtime and memory",
        "",
        _runtime_table(runtime),
        "",
        "![Accuracy versus inference compute](figures/mcv_accuracy_vs_inference_compute.png)",
        "",
        "Calibration is a one-time per-Actor cost. Calibrated efficiency includes one eighth of that cost per reported replicate; direct MCV has no calibration cost.",
        "",
        "## Sampling and resource diagnostics",
        "",
        f"For late Actors, antithetic pairing alone has payoff and martingale variance ratios **{payoff_antithetic_ratio:.3g}** and **{martingale_antithetic_ratio:.3g}** relative to ordinary independent-path averaging. Thus antithetic sampling already removes **{100*(1-payoff_antithetic_ratio):.1f}%** of payoff variance by this diagnostic, and the reported MCV gain is additional. Exponential-weight ESS fractions range from **{ess_min:.3f}** to **{ess_max:.3f}**.",
        "",
        f"The completed Slurm allocations used **{total_gpu_hours:.3f} GH200 GPU-hours** across **{len(completed_jobs)}** jobs. Exact job IDs, nodes, elapsed times and launcher commits are in `execution_manifest.json`.",
        "",
        "## Correctness and interpretation",
        "",
        f"The martingale-zero 95% interval contains zero in **{100*correctness['martingale_zero_in_ci_fraction']:.1f}%** of direct-MCV rows. The paired MCV-minus-ordinary mean interval contains zero in **{100*correctness['mcv_plain_mean_agreement_ci_fraction']:.1f}%** of corrected rows. These are multiple 95% checks, so an occasional miss is expected. Metric replay checked **{metric_audit['rows_checked']}** rows with maximum discrepancy **{metric_audit['maximum_absolute_difference']:.3g}**. Coverage passed: **{coverage['passed']}**; correctness audit passed: **{correctness['passed']}**.",
        "",
        "MCV changes sampling variance only. Differences in posthoc Value RMSE across variants at finite $M$ are MC realization effects; the frozen Actor and Euler bias are unchanged.",
        "",
        "## Answers",
        "",
        f"1. **Does MCV help?** {'Yes' if decision == 'GO' else 'Not reliably enough'}. The best late-Actor mean effective variance reduction is **{_mean(best_late, 'effective_variance_reduction'):.3g}×**.",
        f"2. **Does it help more for later Actors?** {'Yes' if helps_more_late else 'No'}. Direct-MCV mean variance ratio changes from **{_mean(early, 'variance_ratio_replicates'):.3g}** at a1 to **{direct_ratio:.3g}** for the late Actors.",
        f"3. **Paths saved for the same standard error:** approximately **{100*path_fraction:.1f}%**, or a **{path_factor:.3g}×** path reduction, for the preferred late-Actor variant under the observed $1/M$ scaling.",
        f"4. **Does $\alpha=1$ suffice?** {'Yes; calibration gives less than a 10% relative variance advantage.' if direct_suffices else f'No; calibrated mean ratio {calibrated_ratio:.3g} is materially below direct {direct_ratio:.3g}.'}",
        f"5. **Accuracy per GPU second:** the preferred estimator gives **{efficiency:.3g}×** effective variance reduction per GPU second, with mean kernel overhead **{100*overhead:.2f}%**. Its total RMSE benefit can be smaller when Actor/Euler bias dominates.",
        f"6. **Optional inference enhancement:** {'Yes. Keep it inference-only and default to the direct coefficient when calibration is unnecessary.' if decision == 'GO' else 'No under the current evidence; retain the code as a diagnostic only.'}",
    ]
    (RESULTS_DIR / "MCV_QUICK_VALIDATION.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    ensure_dirs()
    variance, runtime, calibration, runs = load_records()
    coverage = _coverage_audit(variance)
    if not coverage["passed"]:
        raise RuntimeError(f"incomplete experiment coverage: {coverage['missing'][:5]}")
    audit_rows, audit_summary = _metric_audit(variance, runs)
    correctness = _correctness_audit(variance, runtime, runs)
    write_csv(RESULTS_DIR / "mcv_variance_by_iteration.csv", variance)
    write_csv(RESULTS_DIR / "mcv_runtime.csv", runtime)
    write_csv(RESULTS_DIR / "mcv_calibration.csv", calibration)
    write_csv(RESULTS_DIR / "metric_recomputation_audit.csv", audit_rows)
    atomic_json(RESULTS_DIR / "metric_recomputation_summary.json", audit_summary)
    atomic_json(RESULTS_DIR / "coverage_audit.json", coverage)
    atomic_json(RESULTS_DIR / "correctness_audit.json", correctness)
    protocol = _json(PACKAGE_DIR / "protocol.json")
    atomic_json(
        RESULTS_DIR / "reproducibility_metadata.json",
        {
            "branch": "ir-mlp-tangent-actor-mcv-quick",
            "aggregation_parent_commit": subprocess.check_output(
                ["git", "rev-parse", "HEAD"],
                cwd=PACKAGE_DIR.parents[2], text=True,
            ).strip(),
            "protocol": protocol,
            "raw_fragments": [run["raw_path"] for run in runs],
            "hardware": "NVIDIA GH200 120GB",
            "software": "Python 3.12.9; PyTorch 2.10.0+cu129",
        },
    )
    _plot_variance(variance)
    _plot_accuracy_compute(variance)
    _render_report(
        variance, runtime, calibration, coverage, audit_summary, correctness
    )
    print(json.dumps({
        "variance_rows": len(variance),
        "runtime_rows": len(runtime),
        "coverage_passed": coverage["passed"],
        "metric_audit_passed": audit_summary["passed"],
        "correctness_passed": correctness["passed"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
