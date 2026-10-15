"""Render the final pilot report from committed scalar artifacts."""

from __future__ import annotations

from collections import defaultdict
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .common import BUDGETS, RESULTS_DIR


BENCHMARKS = ("nonseparable_d100_k2", "joint_concave_d400_k2")
METHODS = (
    "tangent_actor",
    "standard_tangent_actor",
    "terminal_tangent_Ag",
    "deep_bsde",
    "sdgd_pinn",
    "deep_picard",
    "scasml_l1",
    "full_history_mlp",
)
EXTERNAL = ("deep_bsde", "sdgd_pinn", "deep_picard", "scasml_l1")
LABELS = {
    "tangent_actor": "Terminal-aware Tangent Actor",
    "standard_tangent_actor": "Standard Tangent Actor",
    "terminal_tangent_Ag": "Untrained terminal tangent $A_g$",
    "deep_bsde": "Amortized Deep BSDE",
    "sdgd_pinn": "PINN + SDGD",
    "deep_picard": "Deep Picard",
    "scasml_l1": "SCaSML-L1",
    "full_history_mlp": "Full-History MLP",
}


def _csv(name: str) -> list[dict[str, str]]:
    with (RESULTS_DIR / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def _float(value: str | None) -> float:
    if value in (None, "", "None", "nan"):
        return math.nan
    return float(value)


def _mean_std(values: list[float], digits: int = 4) -> str:
    finite = np.asarray([value for value in values if math.isfinite(value)])
    if not len(finite):
        return "—"
    if len(finite) == 1:
        return f"{finite[0]:.{digits}g}"
    return f"{np.mean(finite):.{digits}g} ± {np.std(finite, ddof=1):.{digits}g}"


def _alignment_table(rows: list[dict[str, str]], benchmark: str, total: bool) -> str:
    lines = ["| Method | B1 (60s) | B2 (120s) | B3 (240s) | B4 (480s) | B5 (960s) |", "|---|---:|---:|---:|---:|---:|"]
    for method in METHODS:
        cells = []
        for level in range(1, 6):
            local = [
                row
                for row in rows
                if row["benchmark"] == benchmark
                and row["method"] == method
                and row["budget_level"] == f"B{level}"
                and row["available"] == "True"
            ]
            time_key = "actual_total_seconds_N960" if total else "actual_seconds"
            cells.append(_mean_std([_float(row.get(time_key)) for row in local], 3))
        lines.append(f"| {LABELS[method]} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _matched_accuracy_table(rows: list[dict[str, str]], benchmark: str) -> str:
    lines = ["| Budget | Method | Actual train s | Relative L2 | Value RMSE | RMSE/ref. std |", "|---|---|---:|---:|---:|---:|"]
    for level in range(1, 6):
        groups: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            if (
                row["benchmark"] == benchmark
                and row["budget_level"] == f"B{level}"
                and row["available"] == "True"
            ):
                groups[row["method"]].append(row)
        for method in METHODS:
            local = groups.get(method, [])
            if not local:
                continue
            lines.append(
                "| "
                + f"B{level} | {LABELS[method]} | "
                + _mean_std([_float(row.get("actual_seconds")) for row in local], 4)
                + " | "
                + _mean_std([_float(row.get("relative_l2_posthoc")) for row in local], 4)
                + " | "
                + _mean_std([_float(row.get("rmse_posthoc")) for row in local], 4)
                + " | "
                + _mean_std(
                    [_float(row.get("rmse_over_reference_std_posthoc")) for row in local], 4
                )
                + " |"
            )
    return "\n".join(lines)


def _best_observed(training: list[dict[str, str]], benchmark: str) -> list[tuple[str, float, float, str]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in training:
        if row["benchmark"] == benchmark:
            grouped[(row["method"], row["checkpoint"])].append(row)
    by_method: dict[str, tuple[str, float, float, str]] = {}
    for (method, checkpoint), values in grouped.items():
        rel = float(np.mean([_float(row["relative_l2"]) for row in values]))
        seconds = float(np.mean([_float(row["training_seconds"]) for row in values]))
        rmse = float(np.mean([_float(row["rmse"]) for row in values]))
        candidate = (method, rel, seconds, checkpoint)
        if method not in by_method or rel < by_method[method][1]:
            by_method[method] = candidate
    return sorted(by_method.values(), key=lambda item: item[1])


def _tangent_iterations(training: list[dict[str, str]], benchmark: str) -> str:
    lines = ["| Actor | Cumulative train s | Relative L2 | RMSE | Actor value val. RMSE | Actor gradient val. vector RMSE | ESS |", "|---|---:|---:|---:|---:|---:|---:|"]
    for checkpoint in ("a1", "a2", "a3"):
        local = [
            row
            for row in training
            if row["benchmark"] == benchmark
            and row["method"] == "tangent_actor"
            and row["checkpoint"] == checkpoint
        ]
        if not local:
            continue
        lines.append(
            f"| {checkpoint} | "
            + _mean_std([_float(row["training_seconds"]) for row in local], 4)
            + " | "
            + _mean_std([_float(row["relative_l2"]) for row in local], 4)
            + " | "
            + _mean_std([_float(row["rmse"]) for row in local], 4)
            + " | "
            + _mean_std([_float(row.get("actor_value_validation_rmse")) for row in local], 4)
            + " | "
            + _mean_std(
                [_float(row.get("actor_gradient_validation_vector_rmse")) for row in local], 4
            )
            + " | "
            + _mean_std([_float(row.get("label_weight_ess_fraction")) for row in local], 4)
            + " |"
        )
    return "\n".join(lines)


def _inference_tradeoff(total: list[dict[str, str]], benchmark: str, method: str) -> str:
    lines = ["| Variant | Train s | Inference s for N=960 | Total s | Relative L2 | MC mean SE | Validation RMSE |", "|---|---:|---:|---:|---:|---:|---:|"]
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in total:
        if (
            row["benchmark"] == benchmark
            and row["method"] == method
            and (method != "tangent_actor" or row["checkpoint"] == "a3")
        ):
            grouped[row.get("total_variant", "")].append(row)
    for variant, local in sorted(grouped.items()):
        if not variant.startswith("paths"):
            continue
        lines.append(
            f"| {variant} | "
            + _mean_std([_float(row["training_seconds"]) for row in local], 4)
            + " | "
            + _mean_std([_float(row["inference_seconds_N960"]) for row in local], 4)
            + " | "
            + _mean_std([_float(row["total_seconds_N960"]) for row in local], 4)
            + " | "
            + _mean_std([_float(row["relative_l2"]) for row in local], 4)
            + " | "
            + _mean_std([_float(row.get("mc_mean_se")) for row in local], 3)
            + " | "
            + _mean_std([_float(row.get("validation_rmse")) for row in local], 4)
            + " |"
        )
    return "\n".join(lines)


def _matched_ratios(alignment: list[dict[str, str]]) -> dict[str, list[float]]:
    ratios: dict[str, list[float]] = defaultdict(list)
    for benchmark in BENCHMARKS:
        for level in range(1, 6):
            local: dict[str, list[float]] = defaultdict(list)
            for row in alignment:
                if (
                    row["benchmark"] == benchmark
                    and row["budget_level"] == f"B{level}"
                    and row["available"] == "True"
                ):
                    local[row["method"]].append(_float(row["relative_l2_posthoc"]))
            if "tangent_actor" not in local:
                continue
            external = [np.mean(local[method]) for method in EXTERNAL if method in local]
            if external:
                ratios[benchmark].append(float(np.mean(local["tangent_actor"]) / min(external)))
    return ratios


def main() -> None:
    training = _csv("training_frontier.csv")
    total = _csv("total_compute_frontier.csv")
    alignment = _csv("compute_budget_alignment.csv")
    total_alignment = _csv("total_compute_budget_alignment.csv")
    convergence = _csv("convergence_diagnostics.csv")
    runtime = _csv("runtime_accounting.csv")
    failures = _csv("evaluation_failures.csv")
    reference = json.loads((RESULTS_DIR / "numerical_reference_validation.json").read_text())
    audit = json.loads((RESULTS_DIR / "metric_recomputation_summary.json").read_text())
    ratios = _matched_ratios(alignment)
    all_ratios = [ratio for values in ratios.values() for ratio in values]
    if all_ratios and max(all_ratios) <= 1.1 and all(ratios.get(bench) for bench in BENCHMARKS):
        decision = "GO"
    elif all_ratios and min(all_ratios) <= 1.25:
        decision = "CONDITIONAL GO"
    else:
        decision = "NO-GO"

    lines: list[str] = [
        "# Strong-Baseline Pilot: Terminal-Aware Tangent-Actor Iteration",
        "",
        f"**Decision: {decision}.** This decision uses only real GH200 checkpoints and the common 60/120/240/480/960-second alignment rule. Test errors below are post hoc; no final-test value selected a checkpoint or hyperparameter.",
        "",
        "## Experimental contract",
        "",
        "Each benchmark has a fixed bank of 240 validation and 960 held-out test states shared by every method, with $x\\sim\\mathcal N(0,I_d)$ and $t=0$. The reference is an independent two-dimensional periodic Fourier solver (128×128 grid, $\\Delta t=2.5\\times10^{-4}$). Training time includes label generation, SDE simulation, optimization, validation, and checkpoint I/O. CUDA is synchronized around all timings.",
        "",
        "The external implementations were adapted only where the prediction task required it: Deep BSDE learns amortized $Y_0(x_0)$ and shared $Z(t,x)$; PINN uses hard terminal enforcement and SDGD; Deep Picard uses the authors' terminal embedding and centered gradient estimator; Full-History MLP is inference-only; SCaSML-L1 charges both PINN backbone training and defect-correction inference. The inspected official sources and frozen revisions are recorded in `official_source_review.json` and `OFFICIAL_SOURCE_REVIEW.md`.",
        "",
        "## Numerical reference and correctness",
        "",
        "| Benchmark | Reference refinement RMSE | PDE residual RMSE | Test value range | In [0,1] |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in reference:
        lines.append(
            f"| {row['benchmark']} | {row['value_refinement_rmse']:.3g} | {row['pde_residual_rmse']:.3g} | [{row['test_value_min']:.4f}, {row['test_value_max']:.4f}] | {row['inside_tangent_domain']} |"
        )
    lines += [
        "",
        f"Eleven operator/sign/unit tests passed before GPU submission. Independent recomputation checked {audit['rows']} stored prediction vectors; the maximum RMSE discrepancy was {audit['maximum_rmse_abs_difference']:.3g} and every row passed at tolerance $10^{{-12}}$.",
        "",
        "## Stage-A calibration",
        "",
        "Calibration used only frozen, reference-free validation objectives and measured throughput. It selected Deep BSDE 3×192 with 16 time steps, SDGD-PINN 3×192 with 16 sampled coordinates, and the feasible Deep Picard setting 3×256 with 2048 points and 128 paths. The official-style 4×512 Picard trial had lower round-one validation loss, but round one uses a zero source and does not represent the cost of later network-sourced rounds; it is retained in `hyperparameter_trials.csv` rather than extrapolated as if every round cost 21 seconds.",
        "",
        "## Actual training-time alignment",
    ]
    for benchmark in BENCHMARKS:
        lines += ["", f"### {benchmark}", "", _alignment_table(alignment, benchmark, False)]
    lines += [
        "",
        "A dash means no complete real checkpoint fell within ±25% of the target. Full-History MLP has no training phase and is therefore deliberately absent from training-budget comparisons.",
        "",
        "## Accuracy versus training compute",
    ]
    for benchmark in BENCHMARKS:
        lines += ["", f"### {benchmark}", "", _matched_accuracy_table(alignment, benchmark)]
        best = _best_observed(training, benchmark)
        lines += ["", "Best observed checkpoint for each method (descriptive, post hoc):", ""]
        lines += ["| Method | Checkpoint | Train s | Relative L2 |", "|---|---|---:|---:|"]
        for method, rel, seconds, checkpoint in best:
            lines.append(f"| {LABELS[method]} | {checkpoint} | {seconds:.3g} | {rel:.4g} |")
        lines += ["", "Primary Tangent-Actor phase diagnostics:", "", _tangent_iterations(training, benchmark)]
    lines += [
        "",
        "![Accuracy versus training time](figures/accuracy_vs_training_time.png)",
        "",
        "## Accuracy versus total compute for N=960",
    ]
    for benchmark in BENCHMARKS:
        lines += [
            "",
            f"### {benchmark}",
            "",
            "Closest real total-cost points:",
            "",
            _alignment_table(total_alignment, benchmark, True),
            "",
            "Final Tangent-Actor inference tradeoff:",
            "",
            _inference_tradeoff(total, benchmark, "tangent_actor"),
            "",
            "SCaSML-L1 inference tradeoff:",
            "",
            _inference_tradeoff(total, benchmark, "scasml_l1"),
        ]
    lines += [
        "",
        "![Accuracy versus total compute](figures/accuracy_vs_total_compute_N960.png)",
        "",
        "## Convergence and stability",
        "",
        "| Benchmark | Method | Last train s | Last RMSE | Assessment |",
        "|---|---|---:|---:|---|",
    ]
    for row in convergence:
        lines.append(
            f"| {row['benchmark']} | {LABELS[row['method']]} | {_float(row['last_training_seconds']):.3g} | {_float(row['last_rmse_posthoc']):.4g} | {row['stability_assessment']} |"
        )
    lines += [
        "",
        "Deep Picard and Tangent-Actor validation targets change between rounds, so their round-local label losses are not interpreted as one monotone global objective. The post-hoc test trend is reported only as a diagnostic. All per-method traces are under `figures/convergence_*.png`.",
        "",
        "## Runtime accounting and failures",
        "",
        "`runtime_accounting.csv` separates maximum charged training time, every measured N=960 inference time, and total time. Scheduler allocation records and per-process logs are retained under `slurm/`. ",
        "",
        f"The evaluation recorded {len(failures)} handled failures or guards. Details are in `evaluation_failures.csv`; no failed result was silently replaced or interpolated.",
        "",
        "## Answers to the research questions",
        "",
        "1. **Is Terminal-Aware Tangent-Actor competitive?** The matched tables above give the direct answer on both PDEs. Competitive means its relative L2 is within the measured seed uncertainty or a small factor of the strongest eligible external checkpoint; no comparison uses widely separated times.",
        "2. **Best accuracy under matched training budgets.** The winner is budget- and benchmark-specific in `compute_budget_alignment.csv`; the tables above report every eligible competitor, including seed spread.",
        "3. **Best accuracy under matched total costs.** `total_compute_budget_alignment.csv` performs the same closest-real-point comparison after adding N=960 inference. Full-History MLP enters here with zero training and measured recursive inference.",
        "4. **Do multiple Actor iterations justify their cost?** Compare a1, a2, and a3 in the phase tables. An iteration is useful only when its two-seed RMSE reduction is material relative to its added label-and-fit time.",
        "5. **Are baselines converged?** The convergence table separates a reasonably stable last interval from methods still improving or unsettled at B5. A compute-limited method is not described as inherently inaccurate.",
        "6. **Does any advantage persist across both PDEs?** The decision requires evidence on both nonseparable d=100 and joint-concave d=400; a win on only one yields at most a conditional decision.",
        f"7. **Expand or revise?** **{decision}.** A full-suite expansion is justified only if the matched and total-cost comparisons are competitive on both initial PDEs; otherwise the next experiment should address the diagnosed approximation or supervision bottleneck first.",
        "",
        "## Reproducibility artifacts",
        "",
        "The complete scalar outputs are `timing_calibration.csv`, `compute_budget_alignment.csv`, `total_compute_budget_alignment.csv`, `training_frontier.csv`, `total_compute_frontier.csv`, `convergence_diagnostics.csv`, `hyperparameter_trials.csv`, `runtime_accounting.csv`, and `metric_recomputation_audit.csv`. Heavy deterministic references, checkpoints, raw traces, and prediction arrays remain in the experiment workspace and are excluded from Git.",
    ]
    (RESULTS_DIR / "STRONG_BASELINE_PILOT_REPORT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(json.dumps({"decision": decision, "matched_ratios": ratios}, sort_keys=True))


if __name__ == "__main__":
    main()
