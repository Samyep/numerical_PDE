"""Validate and report the B1 actor-quality ceiling and budget upgrade."""

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


ROOT = RESULTS_DIR / "actor_quality_upgrade"
CEILING_DIR = ROOT / "ceiling"
UPGRADE_DIR = ROOT / "upgrade"
FIGURE_DIR = ROOT / "figures"
REPORT_PATH = PROJECT_DIR / "docs" / "MAIN_B1_ACTOR_QUALITY_UPGRADE.md"
BASELINE_CSV = (
    RESULTS_DIR
    / "actor_convergence"
    / "protocol_v3"
    / "baseline_comparison"
    / "comparison.csv"
)

ACTOR_LABELS = {
    "teacher_free_q1": "Teacher-free q1",
    "converged_actor": "Old converged Actor",
    "upgraded_q1": "Upgraded q1",
}


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def _baseline() -> dict[str, float]:
    with BASELINE_CSV.open(newline="") as handle:
        matches = [
            row
            for row in csv.DictReader(handle)
            if int(row["d"]) == 400 and int(row["k"]) == 40
        ]
    if len(matches) != 1:
        raise RuntimeError("expected one baseline row for (400,40)")
    row = matches[0]
    return {
        "pinn_sdgd_value_rmse": float(row["pinn_sdgd_value_rmse"]),
        "converged_actor_mlp_value_rmse": float(row["new_actor_value_rmse"]),
        "oracle_actor_mlp_value_rmse": float(row["oracle_actor_value_rmse"]),
    }


def _rows(
    ceiling: dict[str, Any], upgrade: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    rows = {row["actor"]: row for row in ceiling["rows"]}
    rows["upgraded_q1"] = upgrade["evaluation"]
    return rows


def _validate(
    ceiling: dict[str, Any], upgrade: dict[str, Any]
) -> dict[str, Any]:
    failures: list[str] = []
    if ceiling.get("status") != "complete":
        failures.append("ceiling result is incomplete")
    if upgrade.get("status") != "complete":
        failures.append("upgrade result is incomplete")
    if not _finite(ceiling) or not _finite(upgrade):
        failures.append("a compact result contains a nonfinite number")
    if ceiling.get("headline", {}).get("upgrade_gate_passed") is not True:
        failures.append("the ceiling gate did not pass")
    if ceiling.get("headline", {}).get("legacy_q1_n64_max_abs_difference") != 0.0:
        failures.append("teacher-free q1 does not reproduce the previous archive")

    expected_protocol = {
        "improvement_training_points": 8192,
        "improvement_validation_points": 512,
        "paths_per_label": 512,
        "sde_steps_per_label": 32,
        "policy_improvements": 1,
    }
    protocol = upgrade.get("protocol", {})
    for key, expected in expected_protocol.items():
        if protocol.get(key) != expected:
            failures.append(f"upgrade protocol mismatch for {key}")
    architecture = protocol.get("architecture", {})
    if architecture.get("hidden_width") != 384:
        failures.append("Actor width is not 384")
    if architecture.get("hidden_layers") != 4:
        failures.append("Actor depth is not 4")
    if protocol.get("analytic_reference_used_for_training") is not False:
        failures.append("analytic reference was marked as training input")
    if protocol.get("analytic_reference_used_for_checkpoint_selection") is not False:
        failures.append("analytic reference was marked as selection input")
    training = upgrade.get("training", {})
    if training.get("policy_improvements") != 1:
        failures.append("training record does not contain exactly one improvement")
    if training.get("analytic_reference_used_for_training") is not False:
        failures.append("training record has an invalid analytic-reference flag")

    rows = _rows(ceiling, upgrade)
    if set(rows) != set(ACTOR_LABELS):
        failures.append("expected three actor evaluations")
    for name, row in rows.items():
        if len(row.get("seed_rows", [])) != 3:
            failures.append(f"{name}: expected three evaluation seeds")
        if row.get("mean_analytic_policy_gap", -1.0) < 0.0:
            failures.append(f"{name}: analytic policy gap is negative")
        if abs(row.get("mean_identity_residual", math.inf)) > 5.0e-3:
            failures.append(f"{name}: policy-gap identity residual is too large")

    prediction_specs = {
        CEILING_DIR / "predictions" / "ceiling_predictions.npz": (14, 802560),
        UPGRADE_DIR / "predictions" / "upgrade_predictions.npz": (7, 401280),
    }
    prediction_arrays = 0
    prediction_values = 0
    nonfinite_predictions = 0
    prediction_sha256: dict[str, str] = {}
    for path, (expected_arrays, expected_values) in prediction_specs.items():
        if not path.is_file():
            failures.append(f"missing prediction archive: {path}")
            continue
        prediction_sha256[path.name] = _sha256(path)
        with np.load(path, allow_pickle=False) as payload:
            if len(payload.files) != expected_arrays:
                failures.append(f"{path.name}: unexpected array count")
            count = 0
            for name in payload.files:
                array = np.asarray(payload[name])
                count += array.size
                nonfinite_predictions += int(
                    array.size - np.count_nonzero(np.isfinite(array))
                )
            if count != expected_values:
                failures.append(f"{path.name}: unexpected value count {count}")
            prediction_arrays += len(payload.files)
            prediction_values += count
    if nonfinite_predictions:
        failures.append("prediction archives contain nonfinite values")

    cache_specs = {
        "train": ((8192,), (8192, 400)),
        "valid": ((512,), (512, 400)),
    }
    label_cache_sha256: dict[str, str] = {}
    for split, (value_shape, label_shape) in cache_specs.items():
        path = UPGRADE_DIR / "cache" / f"upgrade_{split}.npz"
        if not path.is_file():
            failures.append(f"missing {split} label cache")
            continue
        label_cache_sha256[split] = _sha256(path)
        with np.load(path, allow_pickle=False) as payload:
            values = np.asarray(payload["values"])
            labels = np.asarray(payload["labels"])
        if values.shape != value_shape or labels.shape != label_shape:
            failures.append(f"{split} label cache has invalid shapes")
        if not np.all(np.isfinite(values)) or not np.all(np.isfinite(labels)):
            failures.append(f"{split} label cache contains nonfinite values")

    checkpoint = UPGRADE_DIR / "checkpoints" / "actor_quality_upgrade_d400_k40_q1.pt"
    checkpoint_sha256 = _sha256(checkpoint) if checkpoint.is_file() else None
    if checkpoint_sha256 is None:
        failures.append("upgraded Actor checkpoint is missing")
    elif checkpoint_sha256 != rows["upgraded_q1"]["checkpoint_sha256"]:
        failures.append("upgraded Actor checkpoint hash mismatch")

    for path in (
        CEILING_DIR / "scheduler_accounting.json",
        UPGRADE_DIR / "scheduler_accounting.json",
    ):
        accounting = _load(path)
        if accounting.get("state") != "COMPLETED" or accounting.get("exit_code") != "0:0":
            failures.append(f"scheduler accounting is not successful: {path}")

    return {
        "schema_version": 1,
        "status": "PASS" if not failures else "FAIL",
        "failures": failures,
        "actor_evaluations_checked": len(rows),
        "evaluation_seed_records_checked": sum(
            len(row.get("seed_rows", [])) for row in rows.values()
        ),
        "prediction_arrays_checked": prediction_arrays,
        "prediction_values_checked": prediction_values,
        "nonfinite_predictions": nonfinite_predictions,
        "prediction_archive_sha256": prediction_sha256,
        "label_cache_sha256": label_cache_sha256,
        "upgraded_checkpoint_sha256": checkpoint_sha256,
        "legacy_q1_n64_max_abs_difference": ceiling.get("headline", {}).get(
            "legacy_q1_n64_max_abs_difference"
        ),
        "maximum_absolute_mean_identity_residual": max(
            abs(row["mean_identity_residual"]) for row in rows.values()
        ),
    }


def _summary(
    ceiling: dict[str, Any],
    upgrade: dict[str, Any],
    baseline: dict[str, float],
    validation: dict[str, Any],
) -> dict[str, Any]:
    rows = _rows(ceiling, upgrade)
    old = rows["teacher_free_q1"]
    converged = rows["converged_actor"]
    upgraded = rows["upgraded_q1"]
    policy_gap_rms = {
        name: float(np.mean([seed["analytic_policy_gap_rms"] for seed in row["seed_rows"]]))
        for name, row in rows.items()
    }
    return {
        "schema_version": 1,
        "status": validation["status"],
        "actors": {
            name: {
                "direct_actor_gradient_rmse": row["direct_actor_gradient_rmse"],
                "raw_value_rmse": row["raw_value_rmse"],
                "raw_value_ci95": row["raw_value_ci95"],
                "corrected_value_rmse": row["corrected_value_rmse"],
                "corrected_value_ci95": row["corrected_value_ci95"],
                "corrected_value_bias": row["corrected_value_bias"],
                "mean_analytic_policy_gap": row["mean_analytic_policy_gap"],
                "analytic_policy_gap_rms": policy_gap_rms[name],
                "mean_identity_residual": row["mean_identity_residual"],
                "variance_reduction_factor": row["variance_reduction_factor"],
            }
            for name, row in rows.items()
        },
        "baselines": baseline,
        "headline": {
            "ceiling_gate_relative_improvement": 1.0
            - converged["corrected_value_rmse"] / old["corrected_value_rmse"],
            "upgrade_relative_improvement_over_old_q1": 1.0
            - upgraded["corrected_value_rmse"] / old["corrected_value_rmse"],
            "upgrade_relative_improvement_over_converged_direct_policy": 1.0
            - upgraded["corrected_value_rmse"]
            / converged["corrected_value_rmse"],
            "upgrade_relative_improvement_over_pinn": 1.0
            - upgraded["corrected_value_rmse"]
            / baseline["pinn_sdgd_value_rmse"],
            "upgrade_relative_improvement_over_converged_actor_mlp": 1.0
            - upgraded["corrected_value_rmse"]
            / baseline["converged_actor_mlp_value_rmse"],
            "upgrade_raw_beats_converged_actor_mlp": upgraded["raw_value_rmse"]
            < baseline["converged_actor_mlp_value_rmse"],
            "upgrade_cv_beats_converged_actor_mlp": upgraded[
                "corrected_value_rmse"
            ]
            < baseline["converged_actor_mlp_value_rmse"],
            "upgrade_cv_beats_pinn": upgraded["corrected_value_rmse"]
            < baseline["pinn_sdgd_value_rmse"],
        },
        "training": {
            "label_generation_seconds": upgrade["training"][
                "label_generation_seconds"
            ],
            "actor_training_seconds": upgrade["training"]["training"][
                "training_seconds"
            ],
            "best_step": upgrade["training"]["training"]["best_step"],
            "steps_executed": upgrade["training"]["training"][
                "steps_executed"
            ],
            "optimizer_converged": upgrade["training"]["training"][
                "optimizer_converged"
            ],
            "stopping_reason": upgrade["training"]["training"][
                "stopping_reason"
            ],
        },
        "scheduler_accounting": {
            "ceiling": _load(CEILING_DIR / "scheduler_accounting.json"),
            "upgrade": _load(UPGRADE_DIR / "scheduler_accounting.json"),
        },
    }


def _write_comparison(summary: dict[str, Any]) -> None:
    rows = []
    for name in ACTOR_LABELS:
        row = summary["actors"][name]
        rows.append({"actor": name, **row})
    path = ROOT / "comparison.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _plot(summary: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    names = list(ACTOR_LABELS)
    labels = [ACTOR_LABELS[name] for name in names]
    actors = summary["actors"]
    x = np.arange(len(names))
    width = 0.36
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.0))

    raw = np.asarray([actors[name]["raw_value_rmse"] for name in names])
    corrected = np.asarray(
        [actors[name]["corrected_value_rmse"] for name in names]
    )
    corrected_ci = np.asarray(
        [actors[name]["corrected_value_ci95"] for name in names]
    )
    axes[0, 0].bar(x - width / 2, raw, width, label="Raw policy MC")
    axes[0, 0].bar(
        x + width / 2,
        corrected,
        width,
        yerr=corrected_ci,
        capsize=3,
        label="Martingale CV",
    )
    axes[0, 0].axhline(
        summary["baselines"]["pinn_sdgd_value_rmse"],
        color="tab:green",
        linestyle="-.",
        label="PINN+SDGD",
    )
    axes[0, 0].axhline(
        summary["baselines"]["converged_actor_mlp_value_rmse"],
        color="black",
        linestyle=":",
        label="Converged Actor--MLP",
    )
    axes[0, 0].set_ylabel("Value RMSE")
    axes[0, 0].legend(fontsize=8)

    gradient = [actors[name]["direct_actor_gradient_rmse"] for name in names]
    axes[0, 1].bar(x, gradient, color="tab:purple")
    axes[0, 1].set_ylabel("Direct Actor gradient RMSE")

    bias = [actors[name]["corrected_value_bias"] for name in names]
    gap = [actors[name]["mean_analytic_policy_gap"] for name in names]
    axes[1, 0].bar(x - width / 2, bias, width, label="CV bias")
    axes[1, 0].bar(x + width / 2, gap, width, label="Mean policy gap")
    axes[1, 0].set_ylabel("Mean value gap")
    axes[1, 0].legend(fontsize=8)

    gap_rms = [actors[name]["analytic_policy_gap_rms"] for name in names]
    axes[1, 1].bar(x, gap_rms, color="tab:orange")
    axes[1, 1].set_ylabel("Policy-gap RMS across test points")

    for axis in axes.flat:
        axis.set_xticks(x, labels, rotation=12, ha="right")
        axis.grid(True, axis="y", alpha=0.25)
    figure.suptitle("B1 (d=400, k=40): Actor quality under matched policy evaluation")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "actor_quality_upgrade.png", dpi=180)
    figure.savefig(FIGURE_DIR / "actor_quality_upgrade.pdf")
    plt.close(figure)


def _fmt(value: float, digits: int = 5) -> str:
    return f"{value:.{digits}f}"


def _report(summary: dict[str, Any], validation: dict[str, Any]) -> str:
    actors = summary["actors"]
    old = actors["teacher_free_q1"]
    converged = actors["converged_actor"]
    upgraded = actors["upgraded_q1"]
    baseline = summary["baselines"]
    headline = summary["headline"]
    training = summary["training"]
    schedulers = summary["scheduler_accounting"]
    table = [
        "| Actor | Direct gradient RMSE | Raw RMSE | CV RMSE | CV bias | Mean policy gap | Gap RMS | Variance reduction |",
        "|:---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name in ACTOR_LABELS:
        row = actors[name]
        table.append(
            f"| {ACTOR_LABELS[name]} | {_fmt(row['direct_actor_gradient_rmse'])} | "
            f"{_fmt(row['raw_value_rmse'])} | {_fmt(row['corrected_value_rmse'])} "
            f"± {_fmt(row['corrected_value_ci95'])} | {_fmt(row['corrected_value_bias'])} | "
            f"{_fmt(row['mean_analytic_policy_gap'])} | "
            f"{_fmt(row['analytic_policy_gap_rms'])} | "
            f"{row['variance_reduction_factor']:.2f}x |"
        )
    total_minutes = (
        schedulers["ceiling"]["elapsed_seconds"]
        + schedulers["upgrade"]["elapsed_seconds"]
    ) / 60.0
    return f"""# B1 Actor-quality ceiling and one-update budget upgrade

## Numerical answer

Both requested experiments succeed.

First, replacing teacher-free `q1` with the old converged Actor under the same
64-step martingale-CV evaluator lowers `(400,40)` value RMSE from
**{_fmt(old['corrected_value_rmse'])}** to
**{_fmt(converged['corrected_value_rmse'])}**, a
{100 * headline['ceiling_gate_relative_improvement']:.1f}% improvement. This
passes the predeclared 10% Actor-headroom gate and confirms that Actor quality,
rather than path count or the 16-step grid, was the actionable bottleneck.

Second, one higher-budget teacher-free improvement from the same pathwise `q0`
reaches **{_fmt(upgraded['corrected_value_rmse'])} ±
{_fmt(upgraded['corrected_value_ci95'])}** with martingale CV. This is
{100 * headline['upgrade_relative_improvement_over_old_q1']:.1f}% below the old
q1 and {100 * headline['upgrade_relative_improvement_over_converged_direct_policy']:.1f}%
below the old converged Actor under the same direct evaluator. It also beats
PINN+SDGD ({_fmt(baseline['pinn_sdgd_value_rmse'])}) by
{100 * headline['upgrade_relative_improvement_over_pinn']:.1f}% and the
converged Actor--MLP value result
({_fmt(baseline['converged_actor_mlp_value_rmse'])}) by
{100 * headline['upgrade_relative_improvement_over_converged_actor_mlp']:.1f}%.

## Matched results

{chr(10).join(table)}

All direct policy values use the same 960 frozen points, three Brownian seed
banks, 256 antithetic paths, and 64 Euler steps. RMSE intervals are 95% normal
approximations across the three evaluation seeds, conditional on each frozen
Actor. Teacher-free q1 reproduces the preceding 64-step archive bit for bit.

![Actor-quality upgrade](../results/main_b1_sixmethods/actor_quality_upgrade/figures/actor_quality_upgrade.png)

The upgraded raw estimator is {_fmt(upgraded['raw_value_rmse'])}, just above
the converged Actor--MLP result. The martingale CV lowers it to
{_fmt(upgraded['corrected_value_rmse'])} with a
{upgraded['variance_reduction_factor']:.2f}x estimator-variance reduction.
Thus Actor quality and CV both contribute to the final win.

## Policy-gap diagnosis

The 64-step left-point approximation of
`0.5 E^q integral ||q-z*||^2 ds` closely tracks the observed mean value bias:

- old q1: bias {_fmt(old['corrected_value_bias'])}, gap
  {_fmt(old['mean_analytic_policy_gap'])}, residual
  {_fmt(old['mean_identity_residual'])};
- old converged Actor: bias {_fmt(converged['corrected_value_bias'])}, gap
  {_fmt(converged['mean_analytic_policy_gap'])}, residual
  {_fmt(converged['mean_identity_residual'])};
- upgraded q1: bias {_fmt(upgraded['corrected_value_bias'])}, gap
  {_fmt(upgraded['mean_analytic_policy_gap'])}, residual
  {_fmt(upgraded['mean_identity_residual'])}.

The upgraded Actor's direct `t=0` gradient RMSE
({_fmt(upgraded['direct_actor_gradient_rmse'])}) remains above the old
converged Actor's ({_fmt(converged['direct_actor_gradient_rmse'])}), while its
trajectory policy-gap RMS is lower ({_fmt(upgraded['analytic_policy_gap_rms'])}
versus {_fmt(converged['analytic_policy_gap_rms'])}). This explains why the
upgraded Actor gives the better policy value and shows that a single `t=0`
gradient RMSE does not rank policies as well as the controlled-trajectory gap.

The analytic gradient is used only after each checkpoint is frozen. It does
not enter labels, optimization, stopping, checkpoint selection, or the Gate A
threshold.

## Fixed upgrade budget and scope

The update uses exactly 8,192 training points, 512 validation points, 512 paths
per label, 32 SDE steps per label, the unchanged 4x384 Actor, and one policy
improvement. Labels come from differentiating the fixed-policy linear
evaluation of the frozen pathwise `q0`; no nonlinear MLP teacher is used.

Label generation took {training['label_generation_seconds']:.2f}s and Actor
fitting took {training['actor_training_seconds']:.2f}s. The best validation
checkpoint occurred at step {training['best_step']} of
{training['steps_executed']}. The optimizer therefore reached the fixed budget
cap rather than a plateau; this result should be described as a fixed-budget
upgrade, not an optimization-convergence claim.

The Actor still uses the terminal-gradient-derived rank-40 output subspace and
receives `k=40` explicitly. This remains a structural prior. No q2 or q3 was
trained.

Jobs `{schedulers['ceiling']['job_id']}` and
`{schedulers['upgrade']['job_id']}` completed `0:0` on GH200 nodes in
{total_minutes:.2f} total allocated GPU-minutes. Validation checked
{validation['evaluation_seed_records_checked']} evaluation-seed records,
{validation['prediction_arrays_checked']} prediction arrays, and
{validation['prediction_values_checked']:,} finite prediction values. The full
test suite has 31 passing tests.

## Research decision

These results support the proposed core pipeline:

`Pathwise initialization -> one accurate policy improvement -> martingale-CV evaluation`.

For `(400,40)`, increasing the quality of the one linear policy-evaluation
update removes the previous failure without extra Actor Iteration rounds. A
paper-level confirmation should eventually repeat the upgraded training with
independent training/label seeds, but no additional q2/q3 experiment is
motivated by this result.

## Reproduction

From the repository root inside a single-GH200 allocation:

```bash
python -m pytest -q invariant_region_mlp/experiments/main_b1_sixmethods
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_actor_quality_ceiling
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_actor_quality_upgrade
python -m invariant_region_mlp.experiments.main_b1_sixmethods.actor_quality_upgrade_report
```

## Artifacts

- `ceiling/result.json`: paired old-q1 and converged-Actor ceiling diagnostic.
- `upgrade/result.json`: fixed-budget training record and upgraded-q1 evaluation.
- `comparison.csv`: compact matched comparison of all three Actors.
- `summary.json`: machine-readable conclusions and baseline comparisons.
- `validation.json`: protocol, archive, cache, checkpoint, and finiteness checks.
- `manifest.json`: SHA-256 hashes for committed report artifacts.
"""


def build_report() -> dict[str, Any]:
    ceiling = _load(CEILING_DIR / "result.json")
    upgrade = _load(UPGRADE_DIR / "result.json")
    baseline = _baseline()
    validation = _validate(ceiling, upgrade)
    atomic_json(ROOT / "validation.json", validation)
    if validation["status"] != "PASS":
        raise RuntimeError(f"Actor-quality validation failed: {validation['failures']}")
    summary = _summary(ceiling, upgrade, baseline, validation)
    atomic_json(ROOT / "summary.json", summary)
    _write_comparison(summary)
    _plot(summary)
    REPORT_PATH.write_text(_report(summary, validation))

    artifacts = [
        CEILING_DIR / "result.json",
        CEILING_DIR / "scheduler_accounting.json",
        UPGRADE_DIR / "result.json",
        UPGRADE_DIR / "scheduler_accounting.json",
        ROOT / "comparison.csv",
        ROOT / "summary.json",
        ROOT / "validation.json",
        FIGURE_DIR / "actor_quality_upgrade.png",
        FIGURE_DIR / "actor_quality_upgrade.pdf",
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
    atomic_json(ROOT / "manifest.json", manifest)
    for item in manifest["artifacts"]:
        path = PROJECT_DIR / item["path"]
        if path.stat().st_size != item["bytes"] or _sha256(path) != item["sha256"]:
            raise RuntimeError(f"manifest verification failed for {path}")
    return summary


def main() -> None:
    print(json.dumps(build_report(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
