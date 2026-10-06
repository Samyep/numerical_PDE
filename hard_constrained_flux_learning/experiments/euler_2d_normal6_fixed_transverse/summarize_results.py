"""Create figures and a machine-readable audit for the transverse experiment."""

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
OLD_RESULTS = HERE.parent / "euler_2d_diverse64" / "results"
PREFIX = "central_nonnegative6_fixed_transverse"
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
LEARNED = "HCFL-18 learned transverse"
FIXED = "HCFL-6 + fixed transverse Roe"
DISABLED = "HCFL-6 transverse disabled after training"
METHODS = (ROE, HLLC, LEARNED, FIXED)
NEURAL = (LEARNED, FIXED, DISABLED)
COLORS = {
    ROE: "#3B6FB6",
    HLLC: "#E28E2C",
    LEARNED: "#8E63B6",
    FIXED: "#2A9D72",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")


def load_rows() -> list[dict[str, str]]:
    with (RESULTS / "test_case_metrics.csv").open(
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
    figure, axis = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    seed_colors = ("#176B4C", "#2A9D72", "#70C59D")
    for seed, color in zip((0, 1, 2), seed_colors):
        for path, linestyle, label_prefix, alpha in (
            (
                OLD_RESULTS / f"central_nonnegative18_curve_seed{seed}.csv",
                "--",
                "18-cell learned transverse",
                0.7,
            ),
            (
                RESULTS / f"{PREFIX}_curve_seed{seed}.csv",
                "-",
                "6-cell fixed transverse",
                1.0,
            ),
        ):
            with path.open(newline="", encoding="utf-8") as handle:
                curve = list(csv.DictReader(handle))
            updates = [int(row["update"]) for row in curve]
            values = [float(row["rollout_nmae_completed"]) for row in curve]
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
                updates[best],
                values[best],
                color=color,
                s=28,
                zorder=4,
                marker="o" if linestyle == "-" else "x",
            )
    axis.set_xlabel("Optimizer update")
    axis.set_ylabel("Validation rollout NMAE")
    axis.set_title("Validation convergence: learned versus fixed transverse context")
    axis.grid(alpha=0.22)
    axis.legend(ncol=2, fontsize=7)
    figure.savefig(RESULTS / "validation_convergence.png", dpi=220)
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
    offsets = np.arange(len(METHODS), dtype=float) - (len(METHODS) - 1) / 2
    figure, axis = plt.subplots(figsize=(11.3, 4.9), constrained_layout=True)
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
            label=method,
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
    differences = np.asarray(
        [candidate_case[case] - baseline_case[case] for case in cases]
    )
    rng = np.random.default_rng(20261006)
    bootstrap = np.empty(20000, dtype=np.float64)
    for index in range(bootstrap.size):
        sample = rng.integers(0, differences.size, size=differences.size)
        bootstrap[index] = float(differences[sample].mean())
    return {
        "mean_paired_nmae_difference": float(differences.mean()),
        "bootstrap_95_percent_interval": [
            float(np.quantile(bootstrap, 0.025)),
            float(np.quantile(bootstrap, 0.975)),
        ],
        "candidate_case_wins_out_of_36": int((differences < 0.0).sum()),
    }


def audit(rows: list[dict[str, str]]) -> dict:
    summary = load_json(RESULTS / "test_summary.json")
    reports = [
        load_json(RESULTS / f"{PREFIX}_report_seed{seed}.json")
        for seed in (0, 1, 2)
    ]
    fixed_safety = [
        summary["safety"][f"{FIXED} seed {seed}"] for seed in (0, 1, 2)
    ]
    learned = summary["learned_flux_audit"]
    batch_substeps = sum(float(item["batch_substeps"]) for item in fixed_safety)
    entropy_activations = sum(
        float(item["entropy_limiter_active"]) for item in fixed_safety
    )
    positivity_activations = sum(
        float(item["positivity_limiter_active"]) for item in fixed_safety
    )
    checks = {
        "all_training_runs_validation_stopped": all(
            report["stop_reason"]
            == "validation_plateau_at_minimum_learning_rate"
            for report in reports
        ),
        "all_test_trajectories_completed": summary["accuracy"][f"{FIXED}/all"][
            "completion_rate"
        ]
        == 1.0,
        "fully_discrete_entropy_within_tolerance": max(
            float(item["maximum_entropy_balance"]) for item in fixed_safety
        )
        <= 5.0e-7,
        "interface_projection_residual_below_5e-6": max(
            float(item["maximum_interface_residual"]) for item in fixed_safety
        )
        <= 5.0e-6,
        "relative_conservation_closure_below_1e-7": max(
            float(item["maximum_conservation_closure"]) for item in fixed_safety
        )
        <= 1.0e-7,
        "transverse_solver_has_zero_trainable_parameters": all(
            report["transverse_trainable_parameters"] == 0 for report in reports
        ),
        "roe_multipliers_nonnegative_and_bounded": min(
            float(item["minimum_roe_multiplier"]) for item in learned.values()
        )
        >= 0.0
        and max(float(item["maximum_roe_multiplier"]) for item in learned.values())
        <= 2.0,
        "learned_normal_flux_is_not_standard_roe": min(
            float(item["normal_proposal_change_from_standard_roe_l1_ratio"])
            for item in learned.values()
        )
        > 1.0e-4,
    }
    means = {
        method: method_family_mean(rows, method, "all", "nmae")[0]
        for method in (*METHODS, DISABLED)
    }
    result = {
        "passed": all(checks.values()),
        "checks": checks,
        "validation_best_update": [int(report["best_update"]) for report in reports],
        "validation_best_nmae": [
            float(report["best_validation_nmae"]) for report in reports
        ],
        "test_nmae": means,
        "relative_nmae_change_fixed_vs_learned_transverse": (
            means[FIXED] - means[LEARNED]
        )
        / means[LEARNED],
        "relative_nmae_change_fixed_vs_pyclaw_roe": (means[FIXED] - means[ROE])
        / means[ROE],
        "relative_nmae_change_fixed_vs_same_weights_transverse_disabled": (
            means[FIXED] - means[DISABLED]
        )
        / means[DISABLED],
        "paired_comparisons": {
            baseline: paired_comparison(rows, FIXED, baseline)
            for baseline in (LEARNED, ROE, HLLC, DISABLED)
        },
        "maximum_interface_projection_residual": max(
            float(item["maximum_interface_residual"]) for item in fixed_safety
        ),
        "maximum_relative_conservation_closure": max(
            float(item["maximum_conservation_closure"]) for item in fixed_safety
        ),
        "maximum_fully_discrete_entropy_balance": max(
            float(item["maximum_entropy_balance"]) for item in fixed_safety
        ),
        "positivity_fallback_activations": positivity_activations,
        "entropy_fallback_activations": entropy_activations,
        "entropy_fallback_rate_per_batch_substep": entropy_activations
        / max(batch_substeps, 1.0),
        "minimum_entropy_fallback_beta": min(
            float(item["minimum_entropy_beta"]) for item in fixed_safety
        ),
        "learned_flux_audit": learned,
    }
    save_json(RESULTS / "AUDIT.json", result)
    return result


def main() -> None:
    rows = load_rows()
    plot_validation()
    plot_metric(
        rows,
        "nmae",
        "test_nmae_by_family.png",
        "Rollout NMAE (lower is better)",
        "Held-out fixed-64 accuracy by initial-condition family",
    )
    plot_metric(
        rows,
        "density_tv_relative_error_signed",
        "test_density_tv_by_family.png",
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
