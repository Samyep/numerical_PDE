"""Compare original and enriched training data with all method choices frozen."""

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
BASELINE_RESULTS = HERE.parent / "euler_2d_normal6_fixed_transverse" / "results"
PREFIX = "central_nonnegative6_fixed_transverse"
MODEL_LABEL = "HCFL-6 + fixed transverse Roe"
ORIGINAL = "Original 192 trajectories"
ENRICHED = "Enriched 384 trajectories"
ROE = "PyClaw Roe-64"
FAMILIES = (
    "oblique_riemann",
    "contact_shear",
    "oblique_quadrant",
    "radial_interface",
    "colliding_waves",
    "smooth_packet",
)
DISPLAY = {
    "all": "All",
    "oblique_riemann": "Oblique\nRiemann",
    "contact_shear": "Contact /\nshear",
    "oblique_quadrant": "Oblique\nquadrant",
    "radial_interface": "Radial\ninterface",
    "colliding_waves": "Colliding\nwaves",
    "smooth_packet": "Smooth\npacket",
}
MAE_METRICS = {
    "density_mae": r"Density $\rho$",
    "x_velocity_mae": r"Velocity $u$",
    "y_velocity_mae": r"Velocity $v$",
    "pressure_mae": r"Pressure $p$",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def model_rows(path: Path) -> list[dict[str, str]]:
    return [row for row in read_rows(path) if row["method"] == MODEL_LABEL]


def values(
    rows: list[dict[str, str]], family: str, metric: str
) -> np.ndarray:
    selected = [
        float(row[metric])
        for row in rows
        if family == "all" or row["family"] == family
    ]
    return np.asarray(selected, dtype=np.float64)


def seed_means(
    rows: list[dict[str, str]], family: str, metric: str
) -> np.ndarray:
    result = []
    for seed in (0, 1, 2):
        selected = [
            float(row[metric])
            for row in rows
            if int(row["seed"]) == seed
            and (family == "all" or row["family"] == family)
        ]
        result.append(float(np.mean(selected)))
    return np.asarray(result)


def case_means(rows: list[dict[str, str]], metric: str) -> dict[int, float]:
    grouped: defaultdict[int, list[float]] = defaultdict(list)
    for row in rows:
        grouped[int(row["case"])].append(float(row[metric]))
    return {case: float(np.mean(items)) for case, items in grouped.items()}


def paired_difference(
    enriched: list[dict[str, str]], original: list[dict[str, str]], metric: str
) -> dict[str, float | int | list[float]]:
    candidate = case_means(enriched, metric)
    baseline = case_means(original, metric)
    cases = sorted(candidate)
    difference = np.asarray(
        [candidate[case] - baseline[case] for case in cases]
    )
    rng = np.random.default_rng(20261007)
    bootstrap = np.empty(20000, dtype=np.float64)
    for index in range(bootstrap.size):
        sample = rng.integers(0, difference.size, size=difference.size)
        bootstrap[index] = float(difference[sample].mean())
    return {
        "mean_difference_enriched_minus_original": float(difference.mean()),
        "bootstrap_95_percent_interval": [
            float(np.quantile(bootstrap, 0.025)),
            float(np.quantile(bootstrap, 0.975)),
        ],
        "enriched_case_wins_out_of_36": int((difference < 0.0).sum()),
    }


def plot_validation() -> None:
    figure, axis = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    colors = ("#176B4C", "#2A9D72", "#70C59D")
    for seed, color in zip((0, 1, 2), colors):
        for directory, linestyle, label in (
            (BASELINE_RESULTS, "--", ORIGINAL),
            (RESULTS, "-", ENRICHED),
        ):
            path = directory / f"{PREFIX}_curve_seed{seed}.csv"
            with path.open(newline="", encoding="utf-8") as handle:
                curve = list(csv.DictReader(handle))
            update = [int(row["update"]) for row in curve]
            nmae = [float(row["rollout_nmae_completed"]) for row in curve]
            axis.plot(
                update,
                nmae,
                color=color,
                linestyle=linestyle,
                linewidth=1.7,
                label=f"{label}, seed {seed}",
            )
            best = int(np.argmin(nmae))
            axis.scatter(
                update[best],
                nmae[best],
                color=color,
                marker="o" if linestyle == "-" else "x",
                s=30,
                zorder=4,
            )
    axis.set_xlabel("Optimizer update")
    axis.set_ylabel("Validation rollout NMAE")
    axis.set_title("Data-only ablation: validation convergence")
    axis.grid(alpha=0.22)
    axis.legend(ncol=2, fontsize=7)
    figure.savefig(RESULTS / "validation_data_ablation.png", dpi=220)
    plt.close(figure)


def plot_family_metric(
    original: list[dict[str, str]],
    enriched: list[dict[str, str]],
    metric: str,
    filename: str,
    ylabel: str,
    factor: float = 1.0,
) -> None:
    families = ("all", *FAMILIES)
    position = np.arange(len(families), dtype=float)
    width = 0.34
    figure, axis = plt.subplots(figsize=(10.5, 4.8), constrained_layout=True)
    for offset, rows, label, color in (
        (-0.5, original, ORIGINAL, "#8E63B6"),
        (0.5, enriched, ENRICHED, "#2A9D72"),
    ):
        means = []
        errors = []
        for family in families:
            per_seed = factor * seed_means(rows, family, metric)
            means.append(float(per_seed.mean()))
            errors.append(float(per_seed.std()))
        axis.bar(
            position + offset * width,
            means,
            width,
            yerr=errors,
            capsize=2.5,
            color=color,
            label=label,
        )
    if metric == "density_tv_relative_error_signed":
        axis.axhline(0.0, color="0.2", linewidth=0.9)
    axis.set_xticks(position, [DISPLAY[item] for item in families])
    axis.set_ylabel(ylabel)
    axis.set_title("Same HCFL method and test set; training data only changes")
    axis.grid(axis="y", alpha=0.22)
    axis.legend()
    figure.savefig(RESULTS / filename, dpi=220)
    plt.close(figure)


def plot_component_mae(
    original: list[dict[str, str]], enriched: list[dict[str, str]]
) -> None:
    """Plot dimensional primitive-variable MAEs without mixing their units."""
    position = np.arange(len(MAE_METRICS), dtype=float)
    width = 0.34
    figure, axis = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    for offset, rows, label, color in (
        (-0.5, original, ORIGINAL, "#8E63B6"),
        (0.5, enriched, ENRICHED, "#2A9D72"),
    ):
        means = []
        errors = []
        for metric in MAE_METRICS:
            per_seed = seed_means(rows, "all", metric)
            means.append(float(per_seed.mean()))
            errors.append(float(per_seed.std()))
        axis.bar(
            position + offset * width,
            means,
            width,
            yerr=errors,
            capsize=2.5,
            color=color,
            label=label,
        )
    axis.set_xticks(position, list(MAE_METRICS.values()))
    axis.set_ylabel("Held-out rollout MAE (lower is better)")
    axis.set_title("Primitive-variable MAE; no cross-variable aggregation")
    axis.grid(axis="y", alpha=0.22)
    axis.legend()
    figure.savefig(RESULTS / "test_mae_data_ablation.png", dpi=220)
    plt.close(figure)


def main() -> None:
    original = model_rows(BASELINE_RESULTS / "test_case_metrics.csv")
    enriched = model_rows(RESULTS / "test_case_metrics.csv")
    data_audit = load_json(RESULTS / "data_audit.json")
    test_summary = load_json(RESULTS / "test_summary.json")
    reports = [
        load_json(RESULTS / f"{PREFIX}_report_seed{seed}.json")
        for seed in (0, 1, 2)
    ]
    baseline_reports = [
        load_json(BASELINE_RESULTS / f"{PREFIX}_report_seed{seed}.json")
        for seed in (0, 1, 2)
    ]

    plot_validation()
    plot_family_metric(
        original,
        enriched,
        "nmae",
        "test_nmae_data_ablation.png",
        "Rollout NMAE (lower is better)",
    )
    plot_family_metric(
        original,
        enriched,
        "density_tv_relative_error_signed",
        "test_tv_data_ablation.png",
        "Signed final density-TV error (%)",
        factor=100.0,
    )
    plot_component_mae(original, enriched)

    original_nmae = float(values(original, "all", "nmae").mean())
    enriched_nmae = float(values(enriched, "all", "nmae").mean())
    original_tv = float(
        values(original, "all", "density_tv_relative_error_signed").mean()
    )
    enriched_tv = float(
        values(enriched, "all", "density_tv_relative_error_signed").mean()
    )
    original_mae = {
        metric: float(values(original, "all", metric).mean())
        for metric in MAE_METRICS
    }
    enriched_mae = {
        metric: float(values(enriched, "all", metric).mean())
        for metric in MAE_METRICS
    }
    expected_statistics_hash = baseline_reports[0]["train_data_sha256"]
    safety = [
        test_summary["safety"][f"{MODEL_LABEL} seed {seed}"]
        for seed in (0, 1, 2)
    ]
    checks = {
        "data_audit_passed": bool(data_audit["passed"]),
        "all_runs_validation_plateau_stopped": all(
            report["stop_reason"]
            == "validation_plateau_at_minimum_learning_rate"
            for report in reports
        ),
        "normalization_and_metric_scales_frozen": all(
            report["normalization_statistics_sha256"]
            == expected_statistics_hash
            for report in reports
        )
        and test_summary["metric_scale_data_sha256"]
        == expected_statistics_hash,
        "validation_hash_unchanged": all(
            report["validation_data_sha256"]
            == baseline_reports[index]["validation_data_sha256"]
            for index, report in enumerate(reports)
        ),
        "all_test_trajectories_completed": test_summary["accuracy"][
            f"{MODEL_LABEL}/all"
        ]["completion_rate"]
        == 1.0,
        "positivity_fallback_unused": all(
            float(item["positivity_limiter_active"]) == 0.0 for item in safety
        ),
        "all_test_states_positive": test_summary["accuracy"][
            f"{MODEL_LABEL}/all"
        ]["minimum_density"]
        > 0.0
        and test_summary["accuracy"][f"{MODEL_LABEL}/all"][
            "minimum_pressure"
        ]
        > 0.0,
        "fully_discrete_entropy_balance_within_5e-7_tolerance": all(
            float(item["maximum_entropy_balance"]) <= 5.0e-7
            for item in safety
        ),
        "hard_projection_residual_below_1e-6": all(
            float(item["maximum_interface_residual"]) <= 1.0e-6
            for item in safety
        ),
        "conservation_closure_below_1e-7": all(
            float(item["maximum_conservation_closure"]) <= 1.0e-7
            for item in safety
        ),
    }
    result = {
        "passed": all(checks.values()),
        "checks": checks,
        "original_training_trajectories": 192,
        "enriched_training_trajectories": 384,
        "validation_and_test_unchanged": True,
        "original_test_nmae": original_nmae,
        "enriched_test_nmae": enriched_nmae,
        "relative_test_nmae_change": enriched_nmae / original_nmae - 1.0,
        "original_signed_density_tv_error": original_tv,
        "enriched_signed_density_tv_error": enriched_tv,
        "original_primitive_mae": original_mae,
        "enriched_primitive_mae": enriched_mae,
        "relative_primitive_mae_change": {
            metric: enriched_mae[metric] / original_mae[metric] - 1.0
            for metric in MAE_METRICS
        },
        "paired_nmae": paired_difference(enriched, original, "nmae"),
        "paired_primitive_mae": {
            metric: paired_difference(enriched, original, metric)
            for metric in MAE_METRICS
        },
        "best_validation_nmae": [
            float(report["best_validation_nmae"]) for report in reports
        ],
        "best_update": [int(report["best_update"]) for report in reports],
        "data_audit": data_audit,
    }
    (RESULTS / "AUDIT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(result, indent=2), flush=True)
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
