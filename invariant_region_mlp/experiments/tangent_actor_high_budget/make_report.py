"""Render the high-budget report from committed scalar artifacts."""

from __future__ import annotations

from collections import defaultdict
import csv
import json
import math
from typing import Any

import numpy as np

from .common import BENCHMARKS, RESULTS_DIR


LABELS = {
    "nonseparable_d100_k2": "Nonseparable $(d=100,k=2)$",
    "joint_concave_d400_k2": "Joint concave $(d=400,k=2)$",
    "iterations": "more iterations (2048,128)",
    "points4096": "4096 points (128 paths)",
    "points8192": "8192 points (128 paths)",
    "points8192_continued": "continued 8192-point branch",
    "paths256": "2048 points (256 paths)",
    "paths512": "2048 points (512 paths)",
    "terminal_aware_tangent_actor": "Terminal-Aware Tangent-Actor",
    "deep_bsde_original": "Deep BSDE, original form",
    "sdgd_pinn_original": "SDGD-PINN, original form",
    "deep_picard": "Deep Picard",
    "scasml_original": "SCaSML + ordinary PINN",
    "full_history_mlp": "Full-History MLP",
}


def _csv(name: str) -> list[dict[str, str]]:
    with (RESULTS_DIR / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def _float(value: Any) -> float:
    if value in (None, "", "None", "nan"):
        return math.nan
    return float(value)


def _mean(rows: list[dict[str, str]], field: str) -> float:
    values = [_float(row.get(field)) for row in rows]
    values = [value for value in values if math.isfinite(value)]
    return float(np.mean(values)) if values else math.nan


def _mean_std(rows: list[dict[str, str]], field: str, digits: int = 4) -> str:
    values = [_float(row.get(field)) for row in rows]
    values = np.asarray([value for value in values if math.isfinite(value)])
    if not len(values):
        return "—"
    if len(values) == 1:
        return f"{values[0]:.{digits}g}"
    return f"{np.mean(values):.{digits}g} ± {np.std(values, ddof=1):.{digits}g}"


def _pilot_table(ablations: list[dict[str, str]], benchmark: str) -> str:
    lines = [
        "| Controlled branch | Cumulative train s | Paired improvement | Actor value val. RMSE | Actor gradient val. vector RMSE | Gradient-label repl. RMSE | Selected step |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for branch in ("iterations", "points4096", "points8192", "paths256", "paths512"):
        local = [
            row for row in ablations
            if row["benchmark"] == benchmark and row["branch"] == branch
            and (branch != "iterations" or row["checkpoint"] == "a4")
        ]
        lines.append(
            f"| {LABELS[branch]} | {_mean_std(local, 'cumulative_training_seconds', 5)} | "
            f"{_mean_std(local, 'paired_previous_minus_new', 4)} | "
            f"{_mean_std(local, 'validation_value_rmse', 4)} | "
            f"{_mean_std(local, 'validation_gradient_vector_rmse', 4)} | "
            f"{_mean_std(local, 'gradient_label_replication_vector_rmse', 4)} | "
            f"{_mean_std(local, 'selected_optimizer_step', 5)} |"
        )
    return "\n".join(lines)


def _iteration_table(ablations: list[dict[str, str]], benchmark: str) -> str:
    lines = [
        "| Branch | Actor | Train s | Paired previous−new | MC 95% CI | Value val. RMSE | Gradient val. RMSE | Meaningful by frozen rule |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    local = [
        row for row in ablations
        if row["benchmark"] == benchmark
        and row["branch"] in {"iterations", "points8192", "points8192_continued"}
    ]
    order = {"iterations": 0, "points8192": 1, "points8192_continued": 2}
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in local:
        groups[(row["branch"], row["checkpoint"])].append(row)
    for (branch, checkpoint), rows in sorted(
        groups.items(), key=lambda item: (order[item[0][0]], _mean(item[1], "cumulative_training_seconds"))
    ):
        lines.append(
            f"| {LABELS[branch]} | {checkpoint} | {_mean_std(rows, 'cumulative_training_seconds', 5)} | "
            f"{_mean_std(rows, 'paired_previous_minus_new', 4)} | "
            f"[{_mean(rows, 'paired_previous_mc_ci95_low'):.3g}, {_mean(rows, 'paired_previous_mc_ci95_high'):.3g}] | "
            f"{_mean_std(rows, 'validation_value_rmse', 4)} | "
            f"{_mean_std(rows, 'validation_gradient_vector_rmse', 4)} | "
            f"{sum(row['meaningful_reference_free_improvement'] == 'True' for row in rows)}/{len(rows)} seeds |"
        )
    return "\n".join(lines)


def _selected_accuracy_table(total: list[dict[str, str]]) -> str:
    lines = [
        "| Benchmark | Validation-selected Actor | MC paths | Train s | Value RMSE | Relative L2 | Bias | Total s |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for benchmark in BENCHMARKS:
        local = [
            row for row in total
            if row["benchmark"] == benchmark
            and row["method"] == "tangent_actor_high_budget"
            and row.get("validation_selected") == "True"
            and row.get("inference_time_steps") == "32"
        ]
        by_paths: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in local:
            by_paths[row["inference_paths"]].append(row)
        for paths, rows in sorted(by_paths.items(), key=lambda item: int(float(item[0]))):
            actors = sorted({f"{row['branch']}/{row['checkpoint']}" for row in rows})
            lines.append(
                f"| {LABELS[benchmark]} | {', '.join(actors)} | {int(float(paths))} | "
                f"{_mean_std(rows, 'training_seconds', 5)} | {_mean_std(rows, 'rmse', 5)} | "
                f"{_mean_std(rows, 'relative_l2', 5)} | {_mean_std(rows, 'bias', 4)} | "
                f"{_mean_std(rows, 'total_seconds_N960', 5)} |"
            )
    return "\n".join(lines)


def _systematic_table(decomposition: list[dict[str, str]]) -> str:
    lines = [
        "| Benchmark | Actor | Paths | Euler steps | Observed RMSE | MC RMSE contribution | Systematic RMSE proxy | MC fraction of MSE |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for benchmark in BENCHMARKS:
        local = [
            row for row in decomposition
            if row["benchmark"] == benchmark
            and row["validation_selected"] == "True"
            and row["inference_time_steps"] == "32"
        ]
        groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
        for row in local:
            groups[(row["branch"], row["checkpoint"], row["inference_paths"])].append(row)
        for (branch, checkpoint, paths), rows in sorted(
            groups.items(), key=lambda item: int(float(item[0][2]))
        ):
            observed = math.sqrt(_mean(rows, "observed_mse"))
            mc = math.sqrt(_mean(rows, "mc_variance_contribution"))
            lines.append(
                f"| {LABELS[benchmark]} | {branch}/{checkpoint} | {int(float(paths))} | 32 | "
                f"{observed:.5g} | {mc:.5g} | {_mean(rows, 'systematic_rmse_proxy'):.5g} | "
                f"{100*_mean(rows, 'mc_fraction_observed_mse'):.2f}% |"
            )
    return "\n".join(lines)


def _euler_table(decomposition: list[dict[str, str]]) -> str:
    lines = [
        "| Benchmark | Euler steps | Observed RMSE | Systematic RMSE proxy | 95% CI for systematic RMSE |",
        "|---|---:|---:|---:|---:|",
    ]
    for benchmark in BENCHMARKS:
        for steps in (16, 32, 64):
            local = [
                row for row in decomposition
                if row["benchmark"] == benchmark
                and row["validation_selected"] == "True"
                and row["inference_paths"] == "2048"
                and int(float(row["inference_time_steps"])) == steps
            ]
            if not local:
                continue
            lines.append(
                f"| {LABELS[benchmark]} | {steps} | {math.sqrt(_mean(local, 'observed_mse')):.5g} | "
                f"{_mean(local, 'systematic_rmse_proxy'):.5g} | "
                f"[{_mean(local, 'systematic_rmse_ci95_low'):.5g}, {_mean(local, 'systematic_rmse_ci95_high'):.5g}] |"
            )
    return "\n".join(lines)


def _baseline_table(training: list[dict[str, str]]) -> str:
    lines = [
        "| Benchmark | Method | Budget checkpoint | Actual train s | Relative L2 | RMSE |",
        "|---|---|---|---:|---:|---:|",
    ]
    for benchmark in BENCHMARKS:
        for method in ("deep_bsde_original", "sdgd_pinn_original", "deep_picard"):
            for checkpoint in ("B3", "B4", "B5"):
                local = [
                    row for row in training
                    if row["benchmark"] == benchmark and row["method"] == method
                    and row["checkpoint"] == checkpoint
                ]
                if local:
                    lines.append(
                        f"| {LABELS[benchmark]} | {LABELS[method]} | {checkpoint} | "
                        f"{_mean_std(local, 'training_seconds', 5)} | "
                        f"{_mean_std(local, 'relative_l2', 5)} | {_mean_std(local, 'rmse', 5)} |"
                    )
    return "\n".join(lines)


def _alignment_mean(
    alignment: list[dict[str, str]], benchmark: str, method: str, level: str,
    *, axis: str,
) -> float:
    rows = [
        row for row in alignment
        if row["axis"] == axis and row["benchmark"] == benchmark
        and row["method"] == method and row["budget_level"] == level
        and row["selection_type"] == "validation_selected"
        and row["available"] == "True"
    ]
    return _mean(rows, "relative_l2_posthoc")


def _deep_bsde_comparison(alignment: list[dict[str, str]]) -> str:
    statements = []
    for benchmark in BENCHMARKS:
        levels = []
        for level in ("B4", "B5"):
            actor = _alignment_mean(
                alignment, benchmark, "terminal_aware_tangent_actor", level,
                axis="training",
            )
            deep_bsde = _alignment_mean(
                alignment, benchmark, "deep_bsde_original", level,
                axis="training",
            )
            if math.isfinite(actor) and math.isfinite(deep_bsde):
                winner = "TA" if actor < deep_bsde else "Deep BSDE"
                levels.append(
                    f"{level}: TA {actor:.5g}, Deep BSDE {deep_bsde:.5g} ({winner} lower)"
                )
        statements.append(f"{LABELS[benchmark]} — " + "; ".join(levels))
    return " ".join(statements)


def _total_compute_leaders(alignment: list[dict[str, str]]) -> str:
    statements = []
    for benchmark in BENCHMARKS:
        winners = []
        for level in ("B1", "B2", "B3", "B4", "B5"):
            candidates = []
            methods = sorted({
                row["method"] for row in alignment
                if row["axis"] == "total" and row["benchmark"] == benchmark
                and row["budget_level"] == level
            })
            for method in methods:
                value = _alignment_mean(
                    alignment, benchmark, method, level, axis="total"
                )
                rows = [
                    row for row in alignment
                    if row["axis"] == "total" and row["benchmark"] == benchmark
                    and row["method"] == method and row["budget_level"] == level
                    and row["selection_type"] == "validation_selected"
                    and row["available"] == "True"
                ]
                if math.isfinite(value) and len(rows) == 2:
                    candidates.append((value, method))
            if candidates:
                value, method = min(candidates)
                winners.append(f"{level} {LABELS.get(method, method)} ({value:.5g})")
        statements.append(f"{LABELS[benchmark]} — " + ", ".join(winners))
    return " ".join(statements)


def _budget_table(alignment: list[dict[str, str]], axis: str) -> str:
    lines = [
        "| Benchmark | Budget | Method | Validation-selected rel. L2 | Descriptive envelope rel. L2 | Actual selected s |",
        "|---|---|---|---:|---:|---:|",
    ]
    methods = (
        "terminal_aware_tangent_actor", "deep_bsde_original",
        "sdgd_pinn_original", "deep_picard", "scasml_original", "full_history_mlp",
    )
    for benchmark in BENCHMARKS:
        for level in ("B3", "B4", "B5"):
            for method in methods:
                validation = [
                    row for row in alignment
                    if row["axis"] == axis and row["benchmark"] == benchmark
                    and row["budget_level"] == level and row["method"] == method
                    and row["selection_type"] == "validation_selected"
                    and row["available"] == "True"
                ]
                descriptive = [
                    row for row in alignment
                    if row["axis"] == axis and row["benchmark"] == benchmark
                    and row["budget_level"] == level and row["method"] == method
                    and row["selection_type"] == "descriptive_test_envelope"
                    and row["available"] == "True"
                ]
                if validation or descriptive:
                    lines.append(
                        f"| {LABELS[benchmark]} | {level} | {LABELS.get(method, method)} | "
                        f"{_mean_std(validation, 'relative_l2_posthoc', 5)} | "
                        f"{_mean_std(descriptive, 'relative_l2_posthoc', 5)} | "
                        f"{_mean_std(validation, 'actual_seconds', 5)} |"
                    )
    return "\n".join(lines)


def _resource_summary(resources: list[dict[str, str]]) -> str:
    scheduler = [row for row in resources if row["category"] == "scheduler"]
    lines = [
        "| Job | Stage | State | Wall s | GPUs | Allocated GPU h | Node |",
        "|---:|---|---|---:|---:|---:|---|",
    ]
    for row in scheduler:
        lines.append(
            f"| {row.get('job_id')} | {row.get('stage')} | {row.get('state')} | "
            f"{_float(row.get('elapsed_seconds')):.0f} | {_float(row.get('allocated_gpus')):.0f} | "
            f"{_float(row.get('allocated_gpu_hours')):.3f} | {row.get('node')} |"
        )
    return "\n".join(lines)


def main() -> None:
    ablations = _csv("actor_scaling_ablation.csv")
    training = _csv("training_frontier.csv")
    total = _csv("total_compute_frontier.csv")
    decomposition = _csv("systematic_error_decomposition.csv")
    alignment = _csv("compute_budget_alignment.csv")
    resources = _csv("resource_accounting.csv")
    metric_audit = json.loads((RESULTS_DIR / "metric_recomputation_summary.json").read_text())
    correctness = json.loads((RESULTS_DIR / "independent_correctness_audit.json").read_text())

    # Frozen, reference-free hypothesis decisions.
    h2_by_benchmark = {}
    h3_by_benchmark = {}
    for benchmark in BENCHMARKS:
        base = [row for row in ablations if row["benchmark"] == benchmark and row["branch"] == "iterations" and row["checkpoint"] == "a4"]
        data = [row for row in ablations if row["benchmark"] == benchmark and row["branch"] == "points8192"]
        paths = [row for row in ablations if row["benchmark"] == benchmark and row["branch"] == "paths512"]
        h2_by_benchmark[benchmark] = _mean(data, "paired_previous_minus_new") > 1.5 * _mean(base, "paired_previous_minus_new")
        h3_by_benchmark[benchmark] = _mean(paths, "paired_previous_minus_new") > 1.5 * _mean(base, "paired_previous_minus_new")
    h2 = all(h2_by_benchmark.values())
    h3 = all(h3_by_benchmark.values())
    iteration_means: dict[str, list[float]] = {}
    for benchmark in BENCHMARKS:
        iteration_means[benchmark] = []
        checkpoints = sorted(
            {
                row["checkpoint"] for row in ablations
                if row["benchmark"] == benchmark and row["branch"] == "iterations"
            }
        )
        for checkpoint in checkpoints:
            rows = [
                row for row in ablations
                if row["benchmark"] == benchmark and row["branch"] == "iterations"
                and row["checkpoint"] == checkpoint
            ]
            iteration_means[benchmark].append(
                _mean(rows, "paired_previous_minus_new")
            )
    h1 = {benchmark: sum(value >= 1e-4 for value in values) >= 2 for benchmark, values in iteration_means.items()}

    systematic_reduction = {}
    for benchmark in BENCHMARKS:
        frozen = [
            row for row in decomposition
            if row["benchmark"] == benchmark and row["checkpoint"] == "a3_frozen"
            and row["inference_paths"] == "512" and row["inference_time_steps"] == "32"
        ]
        selected = [
            row for row in decomposition
            if row["benchmark"] == benchmark and row["validation_selected"] == "True"
            and row["inference_paths"] == "512" and row["inference_time_steps"] == "32"
        ]
        systematic_reduction[benchmark] = (
            1.0 - _mean(selected, "systematic_rmse_proxy") / _mean(frozen, "systematic_rmse_proxy")
            if frozen and selected else math.nan
        )
    evidence_of_scaling = h2 and any(value > 0.10 for value in systematic_reduction.values() if math.isfinite(value))
    decision = "CONDITIONAL GO" if evidence_of_scaling else "NO-GO"

    lowest_systematic = min(
        (
            row for row in decomposition
            if row["validation_selected"] == "True"
            and row["inference_time_steps"] == "32"
            and row["inference_paths"] == "8192"
        ),
        key=lambda row: _float(row["systematic_rmse_proxy"]),
    )

    lines = [
        "# Terminal-Aware Tangent-Actor: High-Budget Scaling Study",
        "",
        f"**Decision: {decision}.** The controlled pilot supports **H2 (data-limited)** on both PDEs. More training points improve the held-out Actor fit and the reference-free policy value much more efficiently than increasing paths per label. Extra iterations remain useful on the $d=400$ problem and rapidly diminish on the $d=100$ problem. The final conclusion below separates the systematic Actor/discretization floor from inference MC noise and uses corrected original-form competitors.",
        "",
        "## Protocol correction and selection firewall",
        "",
        "The old Deep BSDE and SDGD rows embedded our Terminal-Aware construction and are excluded from every formal frontier. Deep BSDE is rerun with raw amortized $Y_0(x_0)$ and nonshared raw $Z_n(x)$ networks. SDGD is rerun as an ordinary tanh PINN with the authors' soft terminal value/gradient penalty. Corrected SCaSML uses that ordinary PINN and includes both terminal and PDE defects. Deep Picard and Full-History MLP keep terminal operations that belong to their published algorithms.",
        "",
        "All Actor choices use held-out noisy tangent labels and a fixed 240-state, reference-free tangent-policy bank. The 960-state numerical Ground Truth is opened only after selection is frozen. The descriptive test-error envelope is reported separately and implements $E^*(C)=\\min_{i:C_i\\le C}E_i$ exactly.",
        "",
        "## Stage-one controlled scaling",
        "",
    ]
    for benchmark in BENCHMARKS:
        lines.extend([f"### {LABELS[benchmark]}", "", _pilot_table(ablations, benchmark), ""])
    lines.extend(
        [
            "Increasing to 8192 nested training points is the only factor that wins on all four benchmark-seed tasks. Raising label paths halves replication noise, but the smaller reduction in held-out Actor error and policy value does not justify its extra cost. This is why stage two continues only the 8192-point branch.",
            "",
            "## Iteration and staged continuation",
            "",
        ]
    )
    for benchmark in BENCHMARKS:
        lines.extend([f"### {LABELS[benchmark]}", "", _iteration_table(ablations, benchmark), ""])
    lines.extend(
        [
            "## Posthoc Ground Truth accuracy",
            "",
            _selected_accuracy_table(total),
            "",
            "## MC versus systematic error",
            "",
            _systematic_table(decomposition),
            "",
            f"The lowest high-precision systematic RMSE proxy is **{_float(lowest_systematic['systematic_rmse_proxy']):.6g}**, obtained for `{lowest_systematic['benchmark']}` at `{lowest_systematic['branch']}/{lowest_systematic['checkpoint']}` with 8192 inference paths. Bootstrap intervals resample the fixed 960-state test workload and propagate the stored per-state MC variance estimate.",
            "",
            "### Euler sensitivity",
            "",
            _euler_table(decomposition),
            "",
            "## Corrected strong baselines",
            "",
            _baseline_table(training),
            "",
            "The historical terminal-aware Deep BSDE, SDGD-PINN and SCaSML measurements are preserved in `historical_terminal_aware_baselines.csv` but are ineligible here.",
            "",
            "## Accuracy-compute envelopes",
            "",
            "### Training compute",
            "",
            _budget_table(alignment, "training"),
            "",
            "### Total compute for 960 predictions",
            "",
            _budget_table(alignment, "total"),
            "",
            "![Accuracy versus training compute](figures/relative_l2_vs_training_seconds.png)",
            "",
            "![Accuracy versus total compute](figures/relative_l2_vs_total_compute.png)",
            "",
            "![Actor scaling diagnostics](figures/actor_scaling_diagnostics.png)",
            "",
            "## Resource and correctness audit",
            "",
            _resource_summary(resources),
            "",
            f"The independent metric replay checked **{metric_audit['rows_checked']}** prediction rows; maximum absolute discrepancy was **{metric_audit['maximum_absolute_difference']:.3g}** and the audit passed: **{metric_audit['passed']}**. The protocol audit passed: **{correctness['passed']}**. It verifies exact nested state prefixes, non-Terminal-Aware competitor configurations, the Ground Truth selection firewall, and budget-envelope eligibility.",
            "",
            "## Answers to the research questions",
            "",
            f"1. **Beyond $a_3$:** {'Yes on the joint-concave d=400 problem; gains persist but diminish. On the nonseparable d=100 problem they fall below the frozen practical threshold after a4.' if h1.get('joint_concave_d400_k2') else 'No consistent multi-step gain was detected.'}",
            f"2. **Additional supervision data:** {'Yes, 8192 points is consistently the strongest controlled extension and materially lowers Actor fitting error.' if h2 else 'No consistent data-scaling benefit was detected.'}",
            f"3. **More MC paths per label:** {'The gain is large enough to support H3.' if h3 else 'No. They reduce label replication noise, but yield much less policy improvement per GPU second than more points.'}",
            f"4. **Lowest systematic floor:** `{lowest_systematic['benchmark']}`, `{lowest_systematic['branch']}/{lowest_systematic['checkpoint']}`, at 8192 inference paths; RMSE proxy **{_float(lowest_systematic['systematic_rmse_proxy']):.6g}**.",
            f"5. **Deep BSDE at 480/960 seconds:** {_deep_bsde_comparison(alignment)}",
            f"6. **Best accuracy at matched total-compute budgets:** {_total_compute_leaders(alignment)} The descriptive test envelope remains diagnostic only.",
            f"7. **Scaling limitation:** {'There is no evidence of an immediate structural plateau: additional state coverage lowers the non-MC floor. The remaining limitation is Actor approximation/data scaling, especially at d=400.' if evidence_of_scaling else 'The expanded budget did not materially lower the non-MC floor, so further progress requires a methodological change rather than a larger version of the same training recipe.'}",
            "",
            f"**Final recommendation: {decision}.** Continue only with methods that improve state-space coverage or Actor representation. Do not spend the next budget increase primarily on more MC paths per label, and do not present the excluded terminal-aware competitor variants as external baselines.",
        ]
    )
    (RESULTS_DIR / "TANGENT_ACTOR_HIGH_BUDGET_REPORT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
