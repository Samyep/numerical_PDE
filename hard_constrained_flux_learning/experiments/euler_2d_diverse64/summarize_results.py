"""Create the fixed-64 result figures and an executable result audit."""

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
METHODS = ("PyClaw Roe-64", "HLLC-64", "flat18")
METHOD_LABELS = {
    "PyClaw Roe-64": "Roe-64",
    "HLLC-64": "HLLC-64",
    "flat18": "HCFL-64",
}
COLORS = {
    "PyClaw Roe-64": "#3B6FB6",
    "HLLC-64": "#E28E2C",
    "flat18": "#2A9D72",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")


def load_rows() -> list[dict[str, str]]:
    with (RESULTS / "test_case_metrics.csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def selected(rows: list[dict[str, str]], method: str, family: str) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if row["method"] == method and (family == "all" or row["family"] == family)
    ]


def method_family_mean(
    rows: list[dict[str, str]], method: str, family: str, metric: str
) -> tuple[float, float]:
    values = selected(rows, method, family)
    if method != "flat18":
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
    figure, axis = plt.subplots(figsize=(7.4, 4.2), constrained_layout=True)
    for seed, color in zip((0, 1, 2), ("#2A9D72", "#4DAF88", "#79C6A5")):
        with (RESULTS / f"flat18_curve_seed{seed}.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        updates = [int(row["update"]) for row in rows]
        values = [float(row["rollout_nmae_completed"]) for row in rows]
        axis.plot(updates, values, color=color, linewidth=1.7, label=f"HCFL seed {seed}")
        best = int(np.argmin(values))
        axis.scatter(updates[best], values[best], color=color, s=34, zorder=4)
    axis.axhline(0.028396874330503186, color="#E28E2C", linestyle="--", label="Untrained HLLC")
    axis.set_xlabel("Optimizer update")
    axis.set_ylabel("Validation rollout NMAE")
    axis.set_title("Validation-selected convergence")
    axis.grid(alpha=0.22)
    axis.legend(ncol=2, fontsize=8)
    figure.savefig(RESULTS / "validation_convergence.png", dpi=220)
    plt.close(figure)


def plot_accuracy(rows: list[dict[str, str]]) -> None:
    families = ("all", *FAMILIES)
    positions = np.arange(len(families), dtype=float)
    width = 0.24
    figure, axis = plt.subplots(figsize=(10.0, 4.6), constrained_layout=True)
    for method_index, method in enumerate(METHODS):
        means = []
        errors = []
        for family in families:
            mean, error = method_family_mean(rows, method, family, "nmae")
            means.append(mean)
            errors.append(error)
        axis.bar(
            positions + (method_index - 1) * width,
            means,
            width,
            yerr=errors if method == "flat18" else None,
            capsize=3,
            color=COLORS[method],
            label=METHOD_LABELS[method],
        )
    axis.set_xticks(positions, [LABELS[family] for family in families])
    axis.set_ylabel("Rollout NMAE (lower is better)")
    axis.set_title("Held-out fixed-64 accuracy by initial-condition family")
    axis.grid(axis="y", alpha=0.22)
    axis.legend(ncol=3)
    figure.savefig(RESULTS / "test_nmae_by_family.png", dpi=220)
    plt.close(figure)


def plot_tv(rows: list[dict[str, str]]) -> None:
    families = ("all", *FAMILIES)
    positions = np.arange(len(families), dtype=float)
    width = 0.24
    figure, axis = plt.subplots(figsize=(10.0, 4.6), constrained_layout=True)
    for method_index, method in enumerate(METHODS):
        means = []
        errors = []
        for family in families:
            mean, error = method_family_mean(
                rows, method, family, "density_tv_relative_error_signed"
            )
            means.append(100.0 * mean)
            errors.append(100.0 * error)
        axis.bar(
            positions + (method_index - 1) * width,
            means,
            width,
            yerr=errors if method == "flat18" else None,
            capsize=3,
            color=COLORS[method],
            label=METHOD_LABELS[method],
        )
    axis.axhline(0.0, color="0.2", linewidth=0.9)
    axis.set_xticks(positions, [LABELS[family] for family in families])
    axis.set_ylabel("Signed final density-TV error (%)")
    axis.set_title("Held-out roughness: positive values indicate excess TV")
    axis.grid(axis="y", alpha=0.22)
    axis.legend(ncol=3)
    figure.savefig(RESULTS / "test_density_tv_by_family.png", dpi=220)
    plt.close(figure)


def audit(rows: list[dict[str, str]]) -> dict:
    summary = load_json(RESULTS / "test_summary.json")
    data_audit = load_json(RESULTS / "data_audit.json")
    reports = [load_json(RESULTS / f"flat18_report_seed{seed}.json") for seed in (0, 1, 2)]
    seed_nmae = []
    for seed in (0, 1, 2):
        values = [
            float(row["nmae"])
            for row in rows
            if row["method"] == "flat18" and int(row["seed"]) == seed
        ]
        seed_nmae.append(float(np.mean(values)))
    roe_nmae, _ = method_family_mean(rows, "PyClaw Roe-64", "all", "nmae")
    hllc_nmae, _ = method_family_mean(rows, "HLLC-64", "all", "nmae")
    hcfl_nmae = float(np.mean(seed_nmae))

    flat_by_case: defaultdict[int, list[float]] = defaultdict(list)
    for row in rows:
        if row["method"] == "flat18":
            flat_by_case[int(row["case"])].append(float(row["nmae"]))
    flat_case = {case: float(np.mean(values)) for case, values in flat_by_case.items()}
    baseline_case = {
        method: {
            int(row["case"]): float(row["nmae"])
            for row in rows
            if row["method"] == method
        }
        for method in ("PyClaw Roe-64", "HLLC-64")
    }
    wins = {
        method: sum(flat_case[case] < baseline_case[method][case] for case in flat_case)
        for method in baseline_case
    }
    safety = summary["safety"]
    hcfl_safety = [safety[f"flat18 seed {seed}"] for seed in (0, 1, 2)]
    batch_substeps = sum(float(item["batch_substeps"]) for item in hcfl_safety)
    entropy_activations = sum(float(item["entropy_limiter_active"]) for item in hcfl_safety)
    learned = summary["learned_flux_nonzero_audit"]
    checks = {
        "data_audit_passed": bool(data_audit["passed"]),
        "all_training_runs_validation_stopped": all(
            report["stop_reason"] == "validation_plateau_at_minimum_learning_rate"
            for report in reports
        ),
        "all_hcfl_test_trajectories_completed": summary["accuracy"]["flat18/all"]["completion_rate"] == 1.0,
        "positive_density_and_pressure": min(
            summary["accuracy"]["flat18/all"]["minimum_density"],
            summary["accuracy"]["flat18/all"]["minimum_pressure"],
        ) > 0.0,
        "fully_discrete_entropy_balance_nonpositive": max(
            float(item["maximum_entropy_balance"]) for item in hcfl_safety
        ) <= 0.0,
        "interface_projection_residual_below_5e-6": max(
            float(item["maximum_interface_residual"]) for item in hcfl_safety
        ) <= 5.0e-6,
        "relative_conservation_closure_below_1e-7": max(
            float(item["maximum_conservation_closure"]) for item in hcfl_safety
        ) <= 1.0e-7,
        "no_positivity_fallback": sum(
            float(item["positivity_limiter_active"]) for item in hcfl_safety
        ) == 0.0,
        "learned_correction_is_nonzero": min(
            float(item["raw_correction_to_hllc_flux_l1_ratio"])
            for item in learned.values()
        ) > 1.0e-4,
        "hcfl_improves_hllc_mean_nmae": hcfl_nmae < hllc_nmae,
    }
    result = {
        "passed": all(checks.values()),
        "checks": checks,
        "validation_best_update": [int(report["best_update"]) for report in reports],
        "validation_best_nmae": [float(report["best_validation_nmae"]) for report in reports],
        "test_nmae": {
            "PyClaw Roe-64": roe_nmae,
            "HLLC-64": hllc_nmae,
            "HCFL-64_mean": hcfl_nmae,
            "HCFL-64_seed_means": seed_nmae,
            "HCFL-64_between_seed_std": float(np.std(seed_nmae)),
        },
        "hcfl_relative_nmae_improvement_over_hllc": (hllc_nmae - hcfl_nmae) / hllc_nmae,
        "hcfl_relative_nmae_gap_above_roe": (hcfl_nmae - roe_nmae) / roe_nmae,
        "hcfl_case_wins_out_of_36": wins,
        "maximum_interface_projection_residual": max(
            float(item["maximum_interface_residual"]) for item in hcfl_safety
        ),
        "maximum_relative_conservation_closure": max(
            float(item["maximum_conservation_closure"]) for item in hcfl_safety
        ),
        "maximum_fully_discrete_entropy_balance": max(
            float(item["maximum_entropy_balance"]) for item in hcfl_safety
        ),
        "positivity_fallback_activations": sum(
            float(item["positivity_limiter_active"]) for item in hcfl_safety
        ),
        "entropy_fallback_activations": entropy_activations,
        "entropy_fallback_rate_per_batch_substep": entropy_activations
        / max(batch_substeps, 1.0),
        "minimum_entropy_fallback_beta": min(
            float(item["minimum_entropy_beta"]) for item in hcfl_safety
        ),
        "learned_flux_nonzero_audit": learned,
    }
    save_json(RESULTS / "AUDIT.json", result)
    return result


def main() -> None:
    rows = load_rows()
    plot_validation()
    plot_accuracy(rows)
    plot_tv(rows)
    result = audit(rows)
    print(json.dumps(result, indent=2), flush=True)
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
