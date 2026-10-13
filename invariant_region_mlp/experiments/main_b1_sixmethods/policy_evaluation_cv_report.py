"""Validate and report the martingale control-variate policy-evaluation pilot."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import PROJECT_DIR, RESULTS_DIR, atomic_json


CELLS = ((100, 10), (400, 40), (400, 1))
OUTPUT_DIR = RESULTS_DIR / "policy_evaluation_cv_pilot"
ACTOR_ITERATION_DIR = RESULTS_DIR / "actor_iteration_pilot"
FIGURE_DIR = OUTPUT_DIR / "figures"
REPORT_PATH = PROJECT_DIR / "docs" / "MAIN_B1_POLICY_EVALUATION_CV.md"


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _format(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}g}"


def _rmse(value: float, ci95: float) -> str:
    return f"{_format(value)} ± {_format(ci95, 2)}"


def _finite(value: Any) -> bool:
    if value is None or isinstance(value, (str, bool)):
        return True
    if isinstance(value, (int, float)):
        return math.isfinite(float(value))
    if isinstance(value, list):
        return all(_finite(item) for item in value)
    if isinstance(value, dict):
        return all(_finite(item) for item in value.values())
    return True


def _row_lookup(rows: list[dict[str, Any]]) -> dict[tuple[int, int, int, int], dict[str, Any]]:
    return {
        (
            int(row["d"]),
            int(row["k"]),
            int(row["actor_iteration"]),
            int(row["time_steps"]),
        ): row
        for row in rows
    }


def _baselines() -> dict[tuple[int, int], dict[str, Any]]:
    output: dict[tuple[int, int], dict[str, Any]] = {}
    for d, k in CELLS:
        result = _load(ACTOR_ITERATION_DIR / f"result_d{d}_k{k}.json")
        output[d, k] = result["baselines"]
    return output


def _validate(result: dict[str, Any]) -> dict[str, Any]:
    failures: list[str] = []
    rows = result.get("rows", [])
    if result.get("status") != "complete":
        failures.append("result status is not complete")
    if len(rows) != 18:
        failures.append(f"expected 18 rows, found {len(rows)}")
    if not _finite(result):
        failures.append("result contains a nonfinite number")
    if result.get("protocol", {}).get("actor_retraining") is not False:
        failures.append("actor retraining flag is invalid")
    if result.get("protocol", {}).get("same_brownian_for_raw_and_corrected") is not True:
        failures.append("raw/CV Brownian pairing flag is invalid")
    if result.get("protocol", {}).get("additional_actor_calls_for_control_variate") != 0:
        failures.append("control variate used additional actor calls")

    actor_results = {
        (d, k): _load(ACTOR_ITERATION_DIR / f"result_d{d}_k{k}.json")
        for d, k in CELLS
    }
    test_hashes: dict[str, str] = {}
    for cell in result.get("cells", []):
        key = (int(cell["d"]), int(cell["k"]))
        name = f"d{key[0]}_k{key[1]}"
        observed = cell["test_point_sha256"]
        expected = actor_results[key]["test_point_sha256"]
        test_hashes[name] = observed
        if observed != expected:
            failures.append(f"{name}: test-point hash mismatch")

    legacy_differences = [
        float(row["legacy_n16_max_abs_difference"])
        for row in rows
        if row["time_steps"] == 16
    ]
    if len(legacy_differences) != 6 or max(legacy_differences) != 0.0:
        failures.append("16-step raw predictions do not exactly match legacy")
    for row in rows:
        if len(row.get("seed_rows", [])) != 3:
            failures.append("a row does not contain three seed records")
        costs = row.get("costs", [])
        if [cost.get("repetition") for cost in costs] != [0, 1, 2]:
            failures.append("evaluation repetitions are not 0,1,2")
        if not all(cost.get("raw_and_corrected_share_paths") for cost in costs):
            failures.append("a raw/CV result did not share paths")
        if any(cost.get("additional_actor_calls_for_control_variate") != 0 for cost in costs):
            failures.append("a CV result used extra actor calls")
        expected_nested = int(row["time_steps"]) > 16
        if any(cost.get("nested_brownian_refinement") != expected_nested for cost in costs):
            failures.append("Brownian nesting metadata mismatch")

    prediction_path = OUTPUT_DIR / "predictions" / "paired_predictions.npz"
    prediction_arrays = 0
    prediction_values = 0
    nonfinite_predictions = 0
    if not prediction_path.is_file():
        failures.append("paired prediction archive is missing")
    else:
        with np.load(prediction_path, allow_pickle=False) as payload:
            prediction_arrays = len(payload.files)
            if prediction_arrays != 36:
                failures.append(
                    f"expected 36 prediction arrays, found {prediction_arrays}"
                )
            for name in payload.files:
                array = np.asarray(payload[name])
                prediction_values += array.size
                if array.shape != (3, 960):
                    failures.append(f"{name}: unexpected shape {array.shape}")
                if not np.all(np.isfinite(array)):
                    nonfinite_predictions += int(np.size(array) - np.isfinite(array).sum())
    if nonfinite_predictions:
        failures.append("prediction archive contains nonfinite values")

    checkpoints: dict[str, dict[str, str]] = {}
    for d, k in CELLS:
        name = f"d{d}_k{k}"
        checkpoints[name] = {}
        for iteration in (0, 1):
            path = (
                ACTOR_ITERATION_DIR
                / "checkpoints"
                / f"actor_iteration_{name}_q{iteration}.pt"
            )
            if not path.is_file():
                failures.append(f"{name}: q{iteration} checkpoint missing")
            else:
                checkpoints[name][f"q{iteration}"] = _sha256(path)
        if len(set(checkpoints[name].values())) != 2:
            failures.append(f"{name}: q0 and q1 checkpoints are not distinct")

    return {
        "schema_version": 1,
        "status": "PASS" if not failures else "FAIL",
        "rows_checked": len(rows),
        "seed_records_checked": sum(len(row.get("seed_rows", [])) for row in rows),
        "prediction_arrays_checked": prediction_arrays,
        "prediction_values_checked": prediction_values,
        "nonfinite_predictions": nonfinite_predictions,
        "legacy_n16_max_abs_difference": max(legacy_differences, default=None),
        "test_point_sha256": test_hashes,
        "checkpoint_sha256": checkpoints,
        "failures": failures,
    }


def _summarize(
    result: dict[str, Any],
    baselines: dict[tuple[int, int], dict[str, Any]],
    validation: dict[str, Any],
) -> dict[str, Any]:
    lookup = _row_lookup(result["rows"])
    cells: list[dict[str, Any]] = []
    for d, k in CELLS:
        row16 = lookup[d, k, 1, 16]
        row64 = lookup[d, k, 1, 64]
        baseline = baselines[d, k]
        cells.append(
            {
                "d": d,
                "k": k,
                "q1_n16_raw_rmse": row16["raw_value_rmse"],
                "q1_n16_corrected_rmse": row16["corrected_value_rmse"],
                "q1_n64_raw_rmse": row64["raw_value_rmse"],
                "q1_n64_corrected_rmse": row64["corrected_value_rmse"],
                "q1_n64_corrected_bias": row64["corrected_value_bias"],
                "q1_n64_bias_squared_fraction_of_mse": (
                    float(row64["corrected_value_bias"]) ** 2
                    / float(row64["corrected_value_rmse"]) ** 2
                ),
                "q1_n64_variance_reduction_factor": row64[
                    "antithetic_variance_reduction_factor"
                ],
                "q1_raw_rmse_change_n64_minus_n16": (
                    row64["raw_value_rmse"] - row16["raw_value_rmse"]
                ),
                "q1_corrected_rmse_change_n64_minus_n16": (
                    row64["corrected_value_rmse"]
                    - row16["corrected_value_rmse"]
                ),
                "q1_corrected_bias_change_n64_minus_n16": (
                    row64["corrected_value_bias"]
                    - row16["corrected_value_bias"]
                ),
                "converged_actor_mlp_rmse": baseline[
                    "converged_actor_mlp_value_rmse"
                ],
                "pinn_sdgd_rmse": baseline["pinn_sdgd_value_rmse"],
                "q1_n64_corrected_beats_converged_actor_mlp": (
                    row64["corrected_value_rmse"]
                    < baseline["converged_actor_mlp_value_rmse"]
                ),
                "q1_n64_corrected_beats_pinn_sdgd": (
                    row64["corrected_value_rmse"]
                    < baseline["pinn_sdgd_value_rmse"]
                ),
            }
        )
    scheduler = _load(OUTPUT_DIR / "scheduler_accounting.json")
    return {
        "schema_version": 1,
        "status": validation["status"],
        "cells": cells,
        "q1_n64_corrected_beats_raw_cells": sum(
            lookup[d, k, 1, 64]["corrected_value_rmse"]
            < lookup[d, k, 1, 64]["raw_value_rmse"]
            for d, k in CELLS
        ),
        "q1_n64_corrected_beats_converged_actor_mlp_cells": sum(
            cell["q1_n64_corrected_beats_converged_actor_mlp"] for cell in cells
        ),
        "q1_n64_corrected_beats_pinn_sdgd_cells": sum(
            cell["q1_n64_corrected_beats_pinn_sdgd"] for cell in cells
        ),
        "cv_rmse_improved_rows": sum(
            row["corrected_value_rmse"] < row["raw_value_rmse"]
            for row in result["rows"]
        ),
        "variance_reduction_factor_range": [
            min(row["antithetic_variance_reduction_factor"] for row in result["rows"]),
            max(row["antithetic_variance_reduction_factor"] for row in result["rows"]),
        ],
        "program_elapsed_seconds": result["elapsed_seconds"],
        "scheduler_accounting": scheduler,
        "actor_retraining": False,
        "actor_structure_prior": result["protocol"]["actor_structure_prior"],
        "nonlinear_mlp_teacher_used": False,
    }


def _plot(
    result: dict[str, Any], baselines: dict[tuple[int, int], dict[str, Any]]
) -> None:
    lookup = _row_lookup(result["rows"])
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(2, 3, figsize=(13.4, 7.2))
    steps = np.asarray((16, 32, 64))
    styles = {0: ("--", 0.65), 1: ("-", 1.0)}
    for column, (d, k) in enumerate(CELLS):
        top = axes[0, column]
        bottom = axes[1, column]
        for iteration in (0, 1):
            linestyle, alpha = styles[iteration]
            raw = np.asarray(
                [lookup[d, k, iteration, n]["raw_value_rmse"] for n in steps]
            )
            raw_ci = np.asarray(
                [lookup[d, k, iteration, n]["raw_value_ci95"] for n in steps]
            )
            corrected = np.asarray(
                [
                    lookup[d, k, iteration, n]["corrected_value_rmse"]
                    for n in steps
                ]
            )
            corrected_ci = np.asarray(
                [
                    lookup[d, k, iteration, n]["corrected_value_ci95"]
                    for n in steps
                ]
            )
            top.errorbar(
                steps,
                raw,
                yerr=raw_ci,
                color="tab:blue",
                marker="o",
                linestyle=linestyle,
                alpha=alpha,
                capsize=3,
                label=f"q{iteration} raw",
            )
            top.errorbar(
                steps,
                corrected,
                yerr=corrected_ci,
                color="tab:orange",
                marker="s",
                linestyle=linestyle,
                alpha=alpha,
                capsize=3,
                label=f"q{iteration} martingale CV",
            )
            reduction = np.asarray(
                [
                    lookup[d, k, iteration, n][
                        "antithetic_variance_reduction_factor"
                    ]
                    for n in steps
                ]
            )
            bottom.plot(
                steps,
                reduction,
                marker="o" if iteration == 0 else "s",
                linewidth=2,
                label=f"q{iteration}",
            )
        baseline = baselines[d, k]
        top.axhline(
            baseline["converged_actor_mlp_value_rmse"],
            color="black",
            linestyle=":",
            label="Converged Actor--MLP",
        )
        top.axhline(
            baseline["pinn_sdgd_value_rmse"],
            color="tab:green",
            linestyle="-.",
            label="PINN+SDGD",
        )
        top.set_title(f"d={d}, k={k}")
        top.set_yscale("log")
        top.grid(True, which="both", alpha=0.25)
        bottom.axhline(1.0, color="black", linewidth=1, alpha=0.6)
        bottom.set_xlabel("Euler time steps")
        bottom.set_xticks(steps)
        bottom.set_yscale("log")
        bottom.grid(True, which="both", alpha=0.25)
    axes[0, 0].set_ylabel("Value RMSE (log scale)")
    axes[1, 0].set_ylabel("Estimator variance reduction (raw / CV)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper center", ncol=6, frameon=False)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    figure.savefig(FIGURE_DIR / "policy_evaluation_cv.png", dpi=180)
    figure.savefig(FIGURE_DIR / "policy_evaluation_cv.pdf")
    plt.close(figure)


def _report(summary: dict[str, Any], rows: list[dict[str, Any]]) -> str:
    lookup = _row_lookup(rows)
    cell_lookup = {(cell["d"], cell["k"]): cell for cell in summary["cells"]}
    main_table = [
        "| (d,k) | q1 raw, n=16 | q1 CV, n=16 | q1 raw, n=64 | q1 CV, n=64 | "
        "Variance reduction, n=64 | CV bias, n=64 |",
        "|:---:|---:|---:|---:|---:|---:|---:|",
    ]
    step_table = [
        "| (d,k) | Estimator | n=16 RMSE | n=32 RMSE | n=64 RMSE | n64-n16 |",
        "|:---:|:---|---:|---:|---:|---:|",
    ]
    for d, k in CELLS:
        row16 = lookup[d, k, 1, 16]
        row64 = lookup[d, k, 1, 64]
        cell = cell_lookup[d, k]
        main_table.append(
            f"| ({d},{k}) | {_rmse(row16['raw_value_rmse'], row16['raw_value_ci95'])} | "
            f"{_rmse(row16['corrected_value_rmse'], row16['corrected_value_ci95'])} | "
            f"{_rmse(row64['raw_value_rmse'], row64['raw_value_ci95'])} | "
            f"{_rmse(row64['corrected_value_rmse'], row64['corrected_value_ci95'])} | "
            f"{_format(row64['antithetic_variance_reduction_factor'])}x | "
            f"{_format(row64['corrected_value_bias'])} |"
        )
        for estimator, label in (("raw", "Raw"), ("corrected", "Martingale CV")):
            values = [lookup[d, k, 1, n][f"{estimator}_value_rmse"] for n in (16, 32, 64)]
            step_table.append(
                f"| ({d},{k}) | {label} | {_format(values[0])} | "
                f"{_format(values[1])} | {_format(values[2])} | "
                f"{_format(values[2] - values[0])} |"
            )

    d100 = cell_lookup[100, 10]
    d400k40 = cell_lookup[400, 40]
    d400k1 = cell_lookup[400, 1]
    scheduler = summary["scheduler_accounting"]
    return f"""# Martingale control-variate Policy Evaluation pilot

## Numerical answer

The proposed control variate works and costs no extra Actor forward calls. It
reduces the antithetic-pair estimator variance for every tested actor and time
resolution. For `q1` at 64 Euler steps the reduction is
**{_format(d100['q1_n64_variance_reduction_factor'])}x** at `(100,10)`,
**{_format(d400k40['q1_n64_variance_reduction_factor'])}x** at `(400,40)`, and
**{_format(d400k1['q1_n64_variance_reduction_factor'])}x** at `(400,1)`.

The practical conclusion is conditional:

- `(100,10)` is a strong success. `q1` RMSE falls from 0.05701 for the
  published 16-step raw estimator to **0.02811** for 64-step martingale CV.
- `(400,40)` remains a clear failure. The corresponding RMSE is **0.14614**
  with bias **+0.12915**, versus 0.08700 for PINN+SDGD and 0.08087 for the
  converged Actor--MLP. Refining 16 to 64 steps changes CV RMSE by only
  {_format(d400k40['q1_corrected_rmse_change_n64_minus_n16'])}.
- `(400,1)` improves modestly from 0.05061 raw at 16 steps to **0.04818** CV at
  16 steps. Increasing to 64 steps gives 0.04900, so extra time resolution does
  not help this cell under the matched Brownian paths.

Thus the `(400,40)` gap is not mainly finite-path variance or 16-step Euler
error. The leading remaining explanation is policy error from the frozen small
`q1` label/training budget. No additional Actor iteration is justified by this
experiment.

## Matched q1 results

{chr(10).join(main_table)}

RMSE intervals are 95% normal approximations across three independent seed
banks. Raw and corrected values use exactly the same 256 antithetic paths and
the same Actor evaluations. Variance reduction is computed from independent
antithetic-pair means, so it includes covariance inside each pair. Empirical
cross-seed variance gives closely matching ratios.

![Policy evaluation control-variate results](../results/main_b1_sixmethods/policy_evaluation_cv_pilot/figures/policy_evaluation_cv.png)

## Time-discretization check

{chr(10).join(step_table)}

The 32- and 64-step paths are Brownian-bridge refinements of the original
16-step paths. For `(100,10)`, refining the grid lowers corrected bias from
0.02752 to 0.02059, so time discretization is material alongside Monte Carlo
variance. For `(400,40)`, corrected bias changes only from 0.13072 to 0.12915;
its large residual error cannot be removed by increasing paths or by this grid
refinement. At `(400,1)`, the small RMSE increase with refinement is consistent
with loss of favorable error cancellation, rather than a variance problem.

These errors are measured against the optimal analytic value `u*`. Therefore
the reported bias combines policy suboptimality and residual Euler bias; it is
not bias of the Monte Carlo estimator for the discretized policy value. The
martingale term has zero expectation and changes only estimator variance.

## One improvement remains enough

At 64 steps with the control variate, q0-to-q1 RMSE changes are:

- `(100,10)`: {_format(lookup[100, 10, 0, 64]['corrected_value_rmse'])} to
  {_format(lookup[100, 10, 1, 64]['corrected_value_rmse'])}.
- `(400,40)`: {_format(lookup[400, 40, 0, 64]['corrected_value_rmse'])} to
  {_format(lookup[400, 40, 1, 64]['corrected_value_rmse'])}.
- `(400,1)`: {_format(lookup[400, 1, 0, 64]['corrected_value_rmse'])} to
  {_format(lookup[400, 1, 1, 64]['corrected_value_rmse'])}.

One policy improvement helps all three cells. This supports the simple pipeline
`pathwise initialization -> one policy improvement -> accurate linear policy
evaluation` and gives no evidence that more Actor Iteration rounds are useful.

## Scope, structure prior, and correctness

No actor was retrained or selected in this experiment. The actors remain
teacher-free with respect to nonlinear MLP labels, but they retain the
terminal-gradient-derived rank-`k` output subspace and receive `k` explicitly.
This is a structural prior and must be stated when comparing with general HJB
methods.

The full suite has 30 passing tests. A linear-terminal exact-policy test reaches
numerical zero variance and fixes the stochastic-integral sign and scaling.
Brownian refinement preserves every coarse increment and antithetic pair. All
six 16-step raw prediction banks match the published pilot **bit for bit**
(maximum absolute difference 0). Artifact validation checked 18 result rows,
54 seed records, 36 prediction arrays, and 103,680 finite predictions.

Job `{scheduler['job_id']}` completed `0:0` on one GH200 in
{scheduler['elapsed_seconds'] / 60:.2f} GPU-minutes. There were no failed runs,
NaNs, or memory failures.

## Reproduction

From the repository root inside a single-GH200 allocation:

```bash
python -m pytest -q invariant_region_mlp/experiments/main_b1_sixmethods
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_policy_evaluation_cv
python -m invariant_region_mlp.experiments.main_b1_sixmethods.policy_evaluation_cv_report
```

## Artifacts

- `result.json`: full protocol, seed rows, costs, and nested-grid diagnostics.
- `comparison.csv`: compact 18-row q0/q1 and 16/32/64-step comparison.
- `summary.json`: machine-readable conclusions and baseline comparisons.
- `validation.json`: parity, hashes, shapes, and finite-value checks.
- `scheduler_accounting.json`: exact Slurm allocation record.
- `manifest.json`: SHA-256 hashes for the report and compact artifacts.
"""


def build_report() -> dict[str, Any]:
    result = _load(OUTPUT_DIR / "result.json")
    baselines = _baselines()
    validation = _validate(result)
    atomic_json(OUTPUT_DIR / "validation.json", validation)
    if validation["status"] != "PASS":
        raise RuntimeError("policy-evaluation CV validation failed")
    summary = _summarize(result, baselines, validation)
    atomic_json(OUTPUT_DIR / "summary.json", summary)
    _plot(result, baselines)
    REPORT_PATH.write_text(_report(summary, result["rows"]))

    artifacts = [
        OUTPUT_DIR / "result.json",
        OUTPUT_DIR / "comparison.csv",
        OUTPUT_DIR / "summary.json",
        OUTPUT_DIR / "validation.json",
        OUTPUT_DIR / "scheduler_accounting.json",
        FIGURE_DIR / "policy_evaluation_cv.png",
        FIGURE_DIR / "policy_evaluation_cv.pdf",
        REPORT_PATH,
    ]
    manifest = {
        "schema_version": 1,
        "artifacts": [
            {
                "path": str(path.relative_to(PROJECT_DIR)),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in artifacts
        ],
    }
    atomic_json(OUTPUT_DIR / "manifest.json", manifest)
    for item in manifest["artifacts"]:
        path = PROJECT_DIR / item["path"]
        if path.stat().st_size != item["bytes"] or _sha256(path) != item["sha256"]:
            raise RuntimeError(f"manifest verification failed for {path}")
    return summary


def main() -> None:
    print(json.dumps(build_report(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
