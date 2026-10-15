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
    lines = [
        "| Actor | Cumulative train s | Relative L2 | RMSE | Actor value val. RMSE | Actor gradient val. vector RMSE | Value-label repl. RMSE | Gradient-label repl. vector RMSE | Weight ESS fraction |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
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
            + _mean_std(
                [_float(row.get("value_label_replication_rmse")) for row in local], 4
            )
            + " | "
            + _mean_std(
                [
                    _float(row.get("gradient_label_replication_vector_rmse"))
                    for row in local
                ],
                4,
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
        final_checkpoint = "a3" if method == "tangent_actor" else "PINN_B5"
        if (
            row["benchmark"] == benchmark
            and row["method"] == method
            and row["checkpoint"] == final_checkpoint
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


def _available_rows(
    rows: list[dict[str, str]], benchmark: str, level: int, method: str
) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if row["benchmark"] == benchmark
        and row["budget_level"] == f"B{level}"
        and row["method"] == method
        and row["available"] == "True"
    ]


def _mean_field(rows: list[dict[str, str]], field: str) -> float:
    values = [_float(row.get(field)) for row in rows]
    return float(np.mean([value for value in values if math.isfinite(value)]))


def _shared_training_summary(alignment: list[dict[str, str]]) -> str:
    lines = [
        "| Benchmark | Budget | Ours: rel. L2 (train s) | Strongest external: rel. L2 (train s) | Ours reduction |",
        "|---|---|---:|---:|---:|",
    ]
    for benchmark in BENCHMARKS:
        for level in range(1, 6):
            ours = _available_rows(alignment, benchmark, level, "tangent_actor")
            external: list[tuple[float, str, list[dict[str, str]]]] = []
            for method in EXTERNAL:
                local = _available_rows(alignment, benchmark, level, method)
                if local:
                    external.append(
                        (_mean_field(local, "relative_l2_posthoc"), method, local)
                    )
            if not ours or not external:
                continue
            ext_rel, ext_method, ext_rows = min(external)
            our_rel = _mean_field(ours, "relative_l2_posthoc")
            gain = 100.0 * (1.0 - our_rel / ext_rel)
            lines.append(
                f"| {benchmark} | B{level} | {our_rel:.6f} "
                f"({_mean_field(ours, 'actual_seconds'):.1f}) | "
                f"{LABELS[ext_method]} {ext_rel:.6f} "
                f"({_mean_field(ext_rows, 'actual_seconds'):.1f}) | {gain:.1f}% |"
            )
    return "\n".join(lines)


def _total_compute_winners(total_alignment: list[dict[str, str]]) -> str:
    lines = [
        "| Benchmark | Budget | Best eligible method | Rel. L2 | Actual total s | Seeds |",
        "|---|---|---|---:|---:|---:|",
    ]
    for benchmark in BENCHMARKS:
        for level in range(1, 6):
            candidates: list[tuple[float, str, list[dict[str, str]]]] = []
            for method in METHODS:
                local = _available_rows(total_alignment, benchmark, level, method)
                if local:
                    candidates.append(
                        (_mean_field(local, "relative_l2_posthoc"), method, local)
                    )
            if not candidates:
                continue
            rel, method, local = min(candidates)
            lines.append(
                f"| {benchmark} | B{level} | {LABELS[method]} | {rel:.6f} | "
                f"{_mean_field(local, 'actual_total_seconds_N960'):.1f} | {len(local)} |"
            )
    return "\n".join(lines)


def _full_history_summary(total: list[dict[str, str]]) -> str:
    lines = [
        "| Benchmark | Best tested configuration | Rel. L2 | Inference s for N=960 |",
        "|---|---|---:|---:|",
    ]
    for benchmark in BENCHMARKS:
        groups: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in total:
            if row["benchmark"] == benchmark and row["method"] == "full_history_mlp":
                groups[row["total_variant"]].append(row)
        candidates = [
            (_mean_field(local, "relative_l2"), variant, local)
            for variant, local in groups.items()
        ]
        rel, variant, local = min(candidates)
        lines.append(
            f"| {benchmark} | {variant} | {rel:.6f} | "
            f"{_mean_field(local, 'inference_seconds_N960'):.3f} |"
        )
    return "\n".join(lines)


def _scheduler_table(rows: list[dict[str, str]]) -> str:
    lines = [
        "| Job | Stage | State | Node | Wall s | GPUs | Allocated GPU h | Productive |",
        "|---:|---|---|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['job_id']} | {row['stage']} | {row['state']} | {row['node']} | "
            f"{_float(row['elapsed_seconds']):.0f} | {_float(row['allocated_gpus']):.0f} | "
            f"{_float(row['allocated_gpu_hours']):.3f} | {row['productive']} |"
        )
    return "\n".join(lines)


def main() -> None:
    training = _csv("training_frontier.csv")
    total = _csv("total_compute_frontier.csv")
    alignment = _csv("compute_budget_alignment.csv")
    total_alignment = _csv("total_compute_budget_alignment.csv")
    convergence = _csv("convergence_diagnostics.csv")
    runtime = _csv("runtime_accounting.csv")
    scheduler = _csv("scheduler_accounting.csv")
    failures = _csv("evaluation_failures.csv")
    reference = json.loads((RESULTS_DIR / "numerical_reference_validation.json").read_text())
    audit = json.loads((RESULTS_DIR / "metric_recomputation_summary.json").read_text())
    ratios = _matched_ratios(alignment)
    # The matched signal is positive on both PDEs, but only B2 contains a
    # complete primary-Actor checkpoint on both benchmarks.  Two k=2 PDEs and
    # two seeds are not enough evidence for an unconditional full-suite claim,
    # especially while Deep BSDE is still improving at B5.
    decision = "CONDITIONAL GO"

    lines: list[str] = [
        "# Strong-Baseline Pilot: Terminal-Aware Tangent-Actor Iteration",
        "",
        f"**Decision: {decision}.** Terminal-Aware Tangent-Actor is clearly competitive at moderate compute on both pilot PDEs, but the evidence does not yet support an unconditional full-suite claim. The primary method has a shared matched training comparison only at B2, the pilot has two seeds and two effective-dimension-two problems, and amortized Deep BSDE is still improving when it obtains the lowest B5 error. Test errors below are post hoc; no final-test value selected a checkpoint, hyperparameter, or inference budget.",
        "",
        "## Experimental contract",
        "",
        "Each benchmark has a fixed bank of 240 validation and 960 held-out test states shared by every method, with $x\\sim\\mathcal N(0,I_d)$ and $t=0$. The reference is an independent two-dimensional periodic Fourier solver (128×128 grid, $\\Delta t=2.5\\times10^{-4}$). Training time includes label generation, SDE simulation, optimization, validation, and checkpoint I/O. CUDA is synchronized around all timings.",
        "",
        "The external implementations were adapted only where the prediction task required it: Deep BSDE learns amortized $Y_0(x_0)$ and shared $Z(t,x)$; PINN uses hard terminal enforcement and SDGD; Deep Picard uses the authors' terminal embedding and centered gradient estimator; Full-History MLP is inference-only; SCaSML-L1 charges both PINN backbone training and defect-correction inference. The inspected official sources and frozen revisions are recorded in `official_source_review.json` and `OFFICIAL_SOURCE_REVIEW.md`.",
        "",
        "## Main measured findings",
        "",
        "The only training budget populated by the primary Tangent-Actor and external methods on both benchmarks is B2. Within the predeclared ±25% window, ours has the lowest mean error in both cases:",
        "",
        _shared_training_summary(alignment),
        "",
        "This is a medium-budget advantage rather than uniform dominance. At B5, Deep BSDE reaches relative $L^2$ 0.001861 on the nonseparable problem and 0.003207 on the joint-concave problem, below the best tested Tangent-Actor MC values 0.004547 and 0.004894. Those Tangent-Actor points cost only 389 and 280 total seconds, respectively, while Deep BSDE costs about 961 seconds.",
        "",
        "The terminal-aware construction is material. Standard Tangent-Actor remains at relative $L^2$ 0.007949 on the nonseparable problem and 0.028118 at $d=400$, whereas the primary method reaches 0.006173 and 0.006904 before increasing inference MC. The untrained terminal tangent is also weaker on both benchmarks.",
        "",
        "![Compute-aligned strong-baseline summary](figures/strong_baseline_accuracy_compute_summary.png)",
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
        f"Twelve operator/sign/unit tests passed. Independent recomputation checked {audit['rows']} stored prediction vectors; the maximum RMSE discrepancy was {audit['maximum_rmse_abs_difference']:.3g} and every row passed at tolerance $10^{{-12}}$.",
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
        "",
        "Stochastic inference budgets are chosen with the independent validation bank within each fixed trained checkpoint; final-test error never ranks inference variants. Across trained checkpoints, the selected row is the closest real total time within the ±25% budget window. The best eligible method at each total-cost budget is:",
        "",
        _total_compute_winners(total_alignment),
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
        "Full-History MLP has no training cost and all tested configurations finish far below B1, so it appears on the continuous total-compute plot but not in a 60-second matched cell. Its best tested results are descriptive post-hoc values:",
        "",
        _full_history_summary(total),
        "",
        "The level-3 recursion is unstable for these cubic value-gradient drivers: increasing its sample base lowers error, but it remains worse than level 2 and becomes expensive. Full-History MLP is therefore an extremely fast, lower-accuracy option in this pilot rather than the accuracy leader.",
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
        "Deep BSDE and SDGD-PINN are still improving at B5 on both benchmarks and must be described as compute-limited. SCaSML is mixed across seeds. Deep Picard reaches only about 238–249 seconds in 40 complete rounds, so B4 and B5 are unavailable and its result is also compute-limited. Deep Picard and Tangent-Actor validation targets change between rounds, so their round-local label losses are not interpreted as one monotone global objective. The post-hoc test trend is reported only as a diagnostic. All per-method traces are under `figures/convergence_*.png`.",
        "",
        "## Runtime accounting and failures",
        "",
        "`runtime_accounting.csv` separates maximum charged training time, every measured N=960 inference time, and total time. Scheduler-level allocation is:",
        "",
        _scheduler_table(scheduler),
        "",
        f"Productive jobs used {sum(_float(row['allocated_gpu_hours']) for row in scheduler if row['productive'] == 'True'):.3f} allocated GPU-hours. Including the canceled infrastructure launch gives {sum(_float(row['allocated_gpu_hours']) for row in scheduler):.3f} GPU-hours. Method frontiers use synchronized single-GPU elapsed times rather than four-GPU allocation wall time. Per-process logs remain under `slurm/` in the workspace.",
        "",
        f"The evaluation recorded {len(failures)} handled failures or guards. Details are in `evaluation_failures.csv`; no failed result was silently replaced or interpolated.",
        "",
        "## Answers to the research questions",
        "",
        "1. **Is Terminal-Aware Tangent-Actor competitive?** Yes at moderate compute. At the common B2 training budget it lowers relative $L^2$ by 33.2% versus Deep BSDE on the nonseparable problem and by 50.6% versus SCaSML-L1 at $d=400$. It does not have a B5 training checkpoint and therefore does not displace Deep BSDE's best long-budget accuracy.",
        "2. **Best accuracy under matched training budgets.** Ours wins the shared B2 comparison on both PDEs. Where ours has no eligible checkpoint, Deep BSDE is the strongest external method at B3–B5 on both problems; SCaSML-L1 is strongest at joint-concave B1, while Standard Tangent-Actor is the single-seed nonseparable B1 leader.",
        "3. **Best accuracy under matched total costs.** Ours has the lowest eligible mean error at nonseparable B1–B4 and joint-concave B1–B3. The nonseparable B1 and B2 Tangent-Actor cells contain only one seed because the other seed lies just outside the time window. Deep BSDE wins B5 on both PDEs; joint-concave B4 has no eligible Tangent-Actor point.",
        "4. **Do multiple Actor iterations justify their cost?** The move from a1 to a2 reduces mean relative error by 17.5% on the nonseparable PDE and 14.4% at $d=400$, so a2 is justified. The move from a2 to a3 improves the nonseparable error by only 5.5% and worsens the $d=400$ mean by 0.23%. A third update is not robustly justified; use a2 by default or a reference-free validation stopping rule.",
        "5. **Are baselines converged?** No for Deep BSDE or SDGD-PINN: both improve through B5. Deep Picard is limited by the 40-round run ending near B3. SCaSML is mixed, and Full-History MLP's tested level-3 recursion remains unstable. These outcomes are compute-limited diagnostics rather than claims of inherent solver inaccuracy.",
        "6. **Does the advantage persist across both PDEs?** Yes at the one shared matched training budget and across all moderate total-cost cells where ours is eligible. The scope is still narrow: two seeds, two $k=2$ PDEs, and one independent reference construction per PDE.",
        f"7. **Expand or revise?** **{decision}.** Keep the method and proceed to a focused confirmation: add a third seed, include higher effective dimension or another independently referenced PDE, retain Deep BSDE through B5 (and a common B6 only if every continuously trained method is extended), and predeclare an a2/validation stopping rule. A broad ICML suite should follow only if that confirmation preserves the advantage.",
        "",
        "## Reproducibility artifacts",
        "",
        "The complete scalar outputs are `timing_calibration.csv`, `compute_budget_alignment.csv`, `total_compute_budget_alignment.csv`, `training_frontier.csv`, `total_compute_frontier.csv`, `convergence_diagnostics.csv`, `hyperparameter_trials.csv`, `runtime_accounting.csv`, `scheduler_accounting.csv`, and `metric_recomputation_audit.csv`. `execution_manifest.json`, `runtime_environment.json`, and `protocol.json` record the jobs, hardware, software, seeds, and frozen configurations. Heavy deterministic references, checkpoints, raw traces, and prediction arrays remain in the experiment workspace and are excluded from Git.",
    ]
    (RESULTS_DIR / "STRONG_BASELINE_PILOT_REPORT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(json.dumps({"decision": decision, "matched_ratios": ratios}, sort_keys=True))


if __name__ == "__main__":
    main()
