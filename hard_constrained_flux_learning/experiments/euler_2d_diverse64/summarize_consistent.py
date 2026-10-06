"""Summarize and audit the method-consistent fixed-64 comparison."""

from __future__ import annotations

from collections import defaultdict
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FAMILIES = (
    "oblique_riemann",
    "contact_shear",
    "oblique_quadrant",
    "radial_interface",
    "colliding_waves",
    "smooth_packet",
)
LABELS = {
    "all": "All",
    "oblique_riemann": "Oblique\nRiemann",
    "contact_shear": "Contact /\nshear",
    "oblique_quadrant": "Oblique\nquadrant",
    "radial_interface": "Radial\ninterface",
    "colliding_waves": "Colliding\nwaves",
    "smooth_packet": "Smooth\npacket",
}
ROE = "PyClaw Roe-64"
HLLC = "HLLC-64"
SIGNED = "HLLC + signed Roe-18"
CONSISTENT = "Central + nonnegative Roe-18"
METHODS = (ROE, HLLC, SIGNED, CONSISTENT)
METHOD_LABELS = {
    ROE: "PyClaw Roe-64",
    HLLC: "HLLC-64",
    SIGNED: "HLLC + signed Roe",
    CONSISTENT: "Central + nonnegative Roe",
}
COLORS = {
    ROE: "#3B6FB6",
    HLLC: "#E28E2C",
    SIGNED: "#8E63B6",
    CONSISTENT: "#2A9D72",
}
NEURAL = (SIGNED, CONSISTENT)
COMPONENT_MAE_METRICS = (
    "density_mae",
    "x_velocity_mae",
    "y_velocity_mae",
    "pressure_mae",
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")


def load_rows() -> list[dict[str, str]]:
    with (RESULTS / "consistent_test_case_metrics.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        return list(csv.DictReader(handle))


def selected(
    rows: list[dict[str, str]], method: str, family: str
) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if row["method"] == method
        and (family == "all" or row["family"] == family)
    ]


def method_family_mean(
    rows: list[dict[str, str]], method: str, family: str, metric: str
) -> tuple[float, float]:
    values = selected(rows, method, family)
    if method not in NEURAL:
        data = np.asarray([float(row[metric]) for row in values])
        return float(data.mean()), 0.0
    seed_means = []
    for seed in (0, 1, 2):
        data = np.asarray(
            [float(row[metric]) for row in values if int(row["seed"]) == seed]
        )
        seed_means.append(float(data.mean()))
    return float(np.mean(seed_means)), float(np.std(seed_means))


def plot_validation() -> None:
    figure, axis = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    seed_colors = ("#176B4C", "#2A9D72", "#70C59D")
    for seed, color in zip((0, 1, 2), seed_colors):
        for prefix, linestyle, label_prefix, alpha in (
            ("flat18", "--", "HLLC + signed", 0.7),
            ("central_nonnegative18", "-", "Central + nonnegative", 1.0),
        ):
            with (RESULTS / f"{prefix}_curve_seed{seed}.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                rows = list(csv.DictReader(handle))
            updates = [int(row["update"]) for row in rows]
            values = [float(row["rollout_nmae_completed"]) for row in rows]
            axis.plot(
                updates,
                values,
                color=color,
                linestyle=linestyle,
                linewidth=1.7,
                alpha=alpha,
                label=f"{label_prefix}, seed {seed}",
            )
            best = int(np.argmin(values))
            axis.scatter(
                updates[best], values[best], color=color, s=30, zorder=4,
                marker="o" if prefix == "central_nonnegative18" else "x",
            )
    axis.set_xlabel("Optimizer update")
    axis.set_ylabel("Validation rollout NMAE")
    axis.set_title("Controlled validation convergence: only flux parameterization changes")
    axis.grid(alpha=0.22)
    axis.legend(ncol=2, fontsize=7)
    figure.savefig(RESULTS / "consistent_validation_convergence.png", dpi=220)
    plt.close(figure)


def plot_metric(
    rows: list[dict[str, str]],
    metric: str,
    output: str,
    ylabel: str,
    title: str,
    factor: float = 1.0,
) -> None:
    families = ("all", *FAMILIES)
    positions = np.arange(len(families), dtype=float)
    width = 0.19
    figure, axis = plt.subplots(figsize=(11.2, 4.8), constrained_layout=True)
    offsets = np.arange(len(METHODS), dtype=float) - (len(METHODS) - 1) / 2
    for method_index, method in enumerate(METHODS):
        means = []
        errors = []
        for family in families:
            mean, error = method_family_mean(rows, method, family, metric)
            means.append(factor * mean)
            errors.append(factor * error)
        axis.bar(
            positions + offsets[method_index] * width,
            means,
            width,
            yerr=errors if method in NEURAL else None,
            capsize=2.5,
            color=COLORS[method],
            label=METHOD_LABELS[method],
        )
    if metric == "density_tv_relative_error_signed":
        axis.axhline(0.0, color="0.2", linewidth=0.9)
    axis.set_xticks(positions, [LABELS[family] for family in families])
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    axis.grid(axis="y", alpha=0.22)
    axis.legend(ncol=4, fontsize=8)
    figure.savefig(RESULTS / output, dpi=220)
    plt.close(figure)


def case_means(
    rows: list[dict[str, str]], method: str, metric: str
) -> dict[int, float]:
    values: defaultdict[int, list[float]] = defaultdict(list)
    for row in rows:
        if row["method"] == method:
            values[int(row["case"])].append(float(row[metric]))
    return {case: float(np.mean(items)) for case, items in values.items()}


def paired_comparison(
    rows: list[dict[str, str]], candidate: str, baseline: str
) -> dict[str, float | int | list[float]]:
    candidate_case = case_means(rows, candidate, "nmae")
    baseline_case = case_means(rows, baseline, "nmae")
    cases = sorted(candidate_case)
    difference = np.asarray(
        [candidate_case[case] - baseline_case[case] for case in cases]
    )
    rng = np.random.default_rng(20261005)
    bootstrap = np.empty(20000, dtype=np.float64)
    for index in range(bootstrap.size):
        sample = rng.integers(0, difference.size, size=difference.size)
        bootstrap[index] = float(difference[sample].mean())
    return {
        "mean_paired_nmae_difference": float(difference.mean()),
        "bootstrap_95_percent_interval": [
            float(np.quantile(bootstrap, 0.025)),
            float(np.quantile(bootstrap, 0.975)),
        ],
        "candidate_case_wins_out_of_36": int((difference < 0.0).sum()),
    }


def audit(rows: list[dict[str, str]]) -> dict:
    summary = load_json(RESULTS / "consistent_test_summary.json")
    data_audit = load_json(RESULTS / "data_audit.json")
    reports = [
        load_json(RESULTS / f"central_nonnegative18_report_seed{seed}.json")
        for seed in (0, 1, 2)
    ]
    seed_means: dict[str, list[float]] = {}
    for method in NEURAL:
        seed_means[method] = [
            float(
                np.mean(
                    [
                        float(row["nmae"])
                        for row in rows
                        if row["method"] == method and int(row["seed"]) == seed
                    ]
                )
            )
            for seed in (0, 1, 2)
        ]
    means = {
        method: method_family_mean(rows, method, "all", "nmae")[0]
        for method in METHODS
    }
    component_mae = {
        method: {
            metric: method_family_mean(rows, method, "all", metric)[0]
            for metric in COMPONENT_MAE_METRICS
        }
        for method in METHODS
    }
    safety = summary["safety"]
    consistent_safety = [
        safety[f"{CONSISTENT} seed {seed}"] for seed in (0, 1, 2)
    ]
    batch_substeps = sum(
        float(item["batch_substeps"]) for item in consistent_safety
    )
    entropy_activations = sum(
        float(item["entropy_limiter_active"]) for item in consistent_safety
    )
    positivity_activations = sum(
        float(item["positivity_limiter_active"]) for item in consistent_safety
    )
    learned = summary["learned_flux_audit"]
    checks = {
        "data_audit_passed": bool(data_audit["passed"]),
        "all_training_runs_validation_stopped": all(
            report["stop_reason"]
            == "validation_plateau_at_minimum_learning_rate"
            for report in reports
        ),
        "all_consistent_test_trajectories_completed": summary["accuracy"][
            f"{CONSISTENT}/all"
        ]["completion_rate"]
        == 1.0,
        "positive_density_and_pressure": min(
            summary["accuracy"][f"{CONSISTENT}/all"]["minimum_density"],
            summary["accuracy"][f"{CONSISTENT}/all"]["minimum_pressure"],
        )
        > 0.0,
        "fully_discrete_entropy_within_tolerance": max(
            float(item["maximum_entropy_balance"])
            for item in consistent_safety
        )
        <= 5.0e-7,
        "interface_projection_residual_below_5e-6": max(
            float(item["maximum_interface_residual"])
            for item in consistent_safety
        )
        <= 5.0e-6,
        "relative_conservation_closure_below_1e-7": max(
            float(item["maximum_conservation_closure"])
            for item in consistent_safety
        )
        <= 1.0e-7,
        "roe_multipliers_nonnegative_and_bounded": min(
            float(item["minimum_roe_multiplier"])
            for item in learned.values()
        )
        >= 0.0
        and max(
            float(item["maximum_roe_multiplier"])
            for item in learned.values()
        )
        <= 2.0,
        "learned_multiplier_change_is_nonzero": min(
            float(item["raw_proposal_change_from_standard_roe_l1_ratio"])
            for item in learned.values()
        )
        > 1.0e-4,
    }
    result = {
        "passed": all(checks.values()),
        "checks": checks,
        "validation_best_update": [
            int(report["best_update"]) for report in reports
        ],
        "validation_best_nmae": [
            float(report["best_validation_nmae"]) for report in reports
        ],
        "test_nmae": {
            **means,
            f"{SIGNED} seed means": seed_means[SIGNED],
            f"{SIGNED} between-seed std": float(np.std(seed_means[SIGNED])),
            f"{CONSISTENT} seed means": seed_means[CONSISTENT],
            f"{CONSISTENT} between-seed std": float(
                np.std(seed_means[CONSISTENT])
            ),
        },
        "test_component_mae": component_mae,
        "relative_nmae_change_consistent_vs_signed": (
            means[CONSISTENT] - means[SIGNED]
        )
        / means[SIGNED],
        "relative_nmae_change_consistent_vs_pyclaw_roe": (
            means[CONSISTENT] - means[ROE]
        )
        / means[ROE],
        "paired_comparisons": {
            SIGNED: paired_comparison(rows, CONSISTENT, SIGNED),
            ROE: paired_comparison(rows, CONSISTENT, ROE),
            HLLC: paired_comparison(rows, CONSISTENT, HLLC),
        },
        "maximum_interface_projection_residual": max(
            float(item["maximum_interface_residual"])
            for item in consistent_safety
        ),
        "maximum_relative_conservation_closure": max(
            float(item["maximum_conservation_closure"])
            for item in consistent_safety
        ),
        "maximum_fully_discrete_entropy_balance": max(
            float(item["maximum_entropy_balance"])
            for item in consistent_safety
        ),
        "positivity_fallback_activations": positivity_activations,
        "entropy_fallback_activations": entropy_activations,
        "entropy_fallback_rate_per_batch_substep": entropy_activations
        / max(batch_substeps, 1.0),
        "minimum_entropy_fallback_beta": min(
            float(item["minimum_entropy_beta"])
            for item in consistent_safety
        ),
        "learned_flux_audit": learned,
    }
    save_json(RESULTS / "consistent_AUDIT.json", result)
    return result


def main() -> None:
    rows = load_rows()
    plot_validation()
    plot_metric(
        rows,
        "nmae",
        "consistent_test_nmae_by_family.png",
        "Rollout NMAE (lower is better)",
        "Held-out fixed-64 accuracy by initial-condition family",
    )
    plot_metric(
        rows,
        "density_tv_relative_error_signed",
        "consistent_test_density_tv_by_family.png",
        "Signed final density-TV error (%)",
        "Held-out roughness: positive values indicate excess TV",
        factor=100.0,
    )
    result = audit(rows)
    print(json.dumps(result, indent=2), flush=True)
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
