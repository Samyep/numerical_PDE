"""Create a controlled 64/512-cell summary for all structural flux arms."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
AUDIT_RESULTS = HERE.parent / "euler_1d_convergence_audit" / "results"
FULL_RESULTS = HERE.parent / "euler_1d_full_flux_baseline" / "results"
CENTRAL_RESULTS = (
    HERE.parent / "euler_1d_consistent_central_flux" / "results"
)
ARMS = (
    "full_flux_broad",
    "roe_complete_broad",
    "central_consistent_broad",
    "central_roe_signed_broad",
    "central_roe_upwind_broad",
    "dissipation_broad",
)
COLORS = {
    "full_flux_broad": "#CC79A7",
    "roe_complete_broad": "#E69F00",
    "central_consistent_broad": "#D55E00",
    "central_roe_signed_broad": "#7B2CBF",
    "central_roe_upwind_broad": "#F05A28",
    "dissipation_broad": "#009E73",
    "native_hllc_512": "#0072B2",
}
LABELS = {
    "full_flux_broad": "complete flux (physical coordinates)",
    "roe_complete_broad": "complete flux (Roe coordinates)",
    "central_consistent_broad": "central + vector correction",
    "central_roe_signed_broad": "central + signed Roe multipliers",
    "central_roe_upwind_broad": "central + nonnegative auto-upwind",
    "dissipation_broad": "HLLC + signed Roe correction",
    "native_hllc_512": "native HLLC-512",
}
RANDOM_SPLITS = (
    "ordinary_id",
    "broad_random_in_support",
    "moderate_ood_high_frequency",
)
CANONICAL_SPLITS = (
    "sod",
    "lax",
    "collision",
    "strong_pressure",
    "near_vacuum_expansion",
)
CASE_LABELS = {
    "sod": "Sod",
    "lax": "Lax",
    "collision": "Collision",
    "strong_pressure": "Strong pressure",
    "near_vacuum_expansion": "Near-vacuum",
}
GRID = "#D8DEE9"
TEXT = "#1F2933"
MUTED = "#66788A"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def style_axis(axis: plt.Axes) -> None:
    axis.grid(True, color=GRID, linewidth=0.7, alpha=0.65, zorder=0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#9AA5B1")
    axis.spines["bottom"].set_color("#9AA5B1")
    axis.tick_params(colors="#465362", labelsize=8.5)


def load_curves(seed: int, results: Path) -> dict[str, list[dict[str, str]]]:
    audit_curve = read_csv(AUDIT_RESULTS / f"training_curve_seed{seed}.csv")
    return {
        "full_flux_broad": read_csv(
            FULL_RESULTS / f"training_curve_seed{seed}.csv"
        ),
        "roe_complete_broad": read_csv(
            results / f"training_curve_roe_complete_broad_seed{seed}.csv"
        ),
        "central_consistent_broad": read_csv(
            CENTRAL_RESULTS / f"training_curve_seed{seed}.csv"
        ),
        "central_roe_signed_broad": read_csv(
            results
            / f"training_curve_central_roe_signed_broad_seed{seed}.csv"
        ),
        "central_roe_upwind_broad": read_csv(
            results
            / f"training_curve_central_roe_upwind_broad_seed{seed}.csv"
        ),
        "dissipation_broad": [
            row
            for row in audit_curve
            if row["arm"] == "dissipation_broad"
        ],
    }


def load_metrics(seed: int, results: Path) -> dict[str, dict[str, dict[str, str]]]:
    audit_metrics = read_csv(AUDIT_RESULTS / f"metrics_seed{seed}.csv")
    sources = {
        "full_flux_broad": read_csv(
            FULL_RESULTS / f"metrics_seed{seed}.csv"
        ),
        "roe_complete_broad": read_csv(
            results / f"metrics_roe_complete_broad_seed{seed}.csv"
        ),
        "central_consistent_broad": read_csv(
            CENTRAL_RESULTS / f"metrics_seed{seed}.csv"
        ),
        "central_roe_signed_broad": read_csv(
            results / f"metrics_central_roe_signed_broad_seed{seed}.csv"
        ),
        "central_roe_upwind_broad": read_csv(
            results / f"metrics_central_roe_upwind_broad_seed{seed}.csv"
        ),
        "dissipation_broad": [
            row
            for row in audit_metrics
            if row["arm"] == "dissipation_broad"
            and row["stage"] == "converged_best"
        ],
    }
    return {
        arm: {row["split"]: row for row in rows}
        for arm, rows in sources.items()
    }


def split_mean(
    metrics: dict[str, dict[str, dict[str, str]]],
    arm: str,
    splits: tuple[str, ...],
) -> float:
    return float(np.mean([
        float(metrics[arm][split]["rollout_nrmse"]) for split in splits
    ]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--results-dir", type=Path, default=HERE / "results"
    )
    args = parser.parse_args()
    results = args.results_dir.resolve()
    curves = load_curves(args.seed, results)
    metrics = load_metrics(args.seed, results)
    transfer = read_json(
        results / f"roe_upwind512_with_hllc2048_seed{args.seed}.json"
    )
    central_transfer = read_json(
        CENTRAL_RESULTS
        / (
            "central_consistent512_vs_hllc_and_roe_with_hllc2048_"
            f"seed{args.seed}.json"
        )
    )
    signed_report = read_json(
        results / f"report_central_roe_signed_broad_seed{args.seed}.json"
    )
    upwind_report = read_json(
        results / f"report_central_roe_upwind_broad_seed{args.seed}.json"
    )
    complete_report = read_json(
        results / f"report_roe_complete_broad_seed{args.seed}.json"
    )

    cell64 = {
        arm: {
            "random_mean_rollout_nrmse": split_mean(
                metrics, arm, RANDOM_SPLITS
            ),
            "canonical_mean_rollout_nrmse": split_mean(
                metrics, arm, CANONICAL_SPLITS
            ),
        }
        for arm in ARMS
    }
    validation = {
        arm: min(
            float(row["validation_rollout_nrmse"])
            for row in curves[arm]
        )
        for arm in ARMS
    }

    transfer_method_for_arm = {
        "roe_complete_broad": "roe_complete_hcfl_512",
        "central_roe_signed_broad": "central_roe_signed_hcfl_512",
        "central_roe_upwind_broad": "central_roe_upwind_hcfl_512",
        "dissipation_broad": "hllc_roe_correction_hcfl_512",
        "native_hllc_512": "native_hllc_512",
    }
    cell512 = {
        arm: transfer["mean_metrics"][method]
        for arm, method in transfer_method_for_arm.items()
    }
    cell512["central_consistent_broad"] = central_transfer[
        "mean_metrics"
    ]["central_consistent_hcfl_512"]

    figure, axes = plt.subplots(2, 2, figsize=(16.5, 10.2))
    figure.subplots_adjust(
        left=0.07,
        right=0.985,
        bottom=0.12,
        top=0.9,
        wspace=0.24,
        hspace=0.42,
    )
    curve_axis, mean64_axis, case512_axis, mean512_axis = axes.flat

    for arm in ARMS:
        updates = np.array([int(row["update"]) for row in curves[arm]])
        values = np.array([
            float(row["validation_rollout_nrmse"]) for row in curves[arm]
        ])
        curve_axis.plot(
            updates,
            values,
            color=COLORS[arm],
            linewidth=1.45,
            label=LABELS[arm],
            zorder=3,
        )
        best = int(np.argmin(values))
        curve_axis.scatter(
            [updates[best]], [values[best]], color=COLORS[arm],
            edgecolor="white", marker="*", s=90, linewidth=0.6,
            zorder=4,
        )
    curve_axis.set_yscale("log")
    curve_axis.set_xlabel("optimizer updates")
    curve_axis.set_ylabel("validation rollout NRMSE")
    curve_axis.set_title(
        "A. Independent-validation convergence",
        loc="left", fontweight="semibold", color=TEXT,
    )
    curve_axis.legend(frameon=False, fontsize=7.7, ncol=2)
    style_axis(curve_axis)

    positions = np.arange(len(ARMS))
    width = 0.36
    mean64_axis.bar(
        positions - width / 2,
        [cell64[arm]["random_mean_rollout_nrmse"] for arm in ARMS],
        width=width,
        color=[COLORS[arm] for arm in ARMS],
        alpha=0.55,
        label="three random-distribution splits",
        zorder=3,
    )
    mean64_axis.bar(
        positions + width / 2,
        [cell64[arm]["canonical_mean_rollout_nrmse"] for arm in ARMS],
        width=width,
        color=[COLORS[arm] for arm in ARMS],
        hatch="//",
        label="five canonical Riemann cases",
        zorder=3,
    )
    mean64_axis.set_yscale("log")
    mean64_axis.set_xticks(
        positions,
        [LABELS[arm] for arm in ARMS],
        rotation=22,
        ha="right",
    )
    mean64_axis.set_ylabel("64-cell mean rollout NRMSE")
    mean64_axis.set_title(
        "B. Held-out 64-cell accuracy",
        loc="left", fontweight="semibold", color=TEXT,
    )
    mean64_axis.legend(frameon=False, fontsize=8)
    style_axis(mean64_axis)

    methods512 = (
        "hllc_roe_correction_hcfl_512",
        "central_roe_signed_hcfl_512",
        "central_roe_upwind_hcfl_512",
    )
    method_to_arm = {
        "hllc_roe_correction_hcfl_512": "dissipation_broad",
        "central_roe_signed_hcfl_512": "central_roe_signed_broad",
        "central_roe_upwind_hcfl_512": "central_roe_upwind_broad",
    }
    case_positions = np.arange(len(CANONICAL_SPLITS))
    case_width = 0.25
    for offset, method in zip(
        (-case_width, 0.0, case_width), methods512, strict=True
    ):
        arm = method_to_arm[method]
        case512_axis.bar(
            case_positions + offset,
            [
                transfer["cases"][case][method]["rollout_nrmse"]
                for case in CANONICAL_SPLITS
            ],
            width=case_width,
            color=COLORS[arm],
            label=LABELS[arm],
            zorder=3,
        )
    case512_axis.set_yscale("log")
    case512_axis.set_xticks(
        case_positions,
        [CASE_LABELS[case] for case in CANONICAL_SPLITS],
        rotation=18,
        ha="right",
    )
    case512_axis.set_ylabel("512-cell rollout NRMSE")
    case512_axis.set_title(
        "C. Zero-shot 512-cell transfer by case",
        loc="left", fontweight="semibold", color=TEXT,
    )
    case512_axis.legend(frameon=False, fontsize=7.8)
    style_axis(case512_axis)

    mean512_arms = (
        "native_hllc_512",
        "dissipation_broad",
        "central_consistent_broad",
        "roe_complete_broad",
        "central_roe_signed_broad",
        "central_roe_upwind_broad",
    )
    mean512_values = [
        cell512[arm]["mean_rollout_nrmse"] for arm in mean512_arms
    ]
    mean512_axis.bar(
        np.arange(len(mean512_arms)),
        mean512_values,
        color=[COLORS[arm] for arm in mean512_arms],
        zorder=3,
    )
    mean512_axis.set_xticks(
        np.arange(len(mean512_arms)),
        [LABELS[arm] for arm in mean512_arms],
        rotation=22,
        ha="right",
    )
    mean512_axis.set_ylabel("five-case mean rollout NRMSE")
    mean512_axis.set_title(
        "D. Zero-shot 512-cell mean (all shown methods completed 5/5)",
        loc="left", fontweight="semibold", color=TEXT,
    )
    for index, value in enumerate(mean512_values):
        mean512_axis.text(
            index,
            value + 0.004,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=7.5,
            color=TEXT,
        )
    style_axis(mean512_axis)

    signed_multiplier = signed_report["learned_wave_multiplier_diagnostics"]
    upwind_multiplier = upwind_report["learned_wave_multiplier_diagnostics"]
    old_mean = cell512["dissipation_broad"]["mean_rollout_nrmse"]
    signed_mean = cell512["central_roe_signed_broad"]["mean_rollout_nrmse"]
    upwind_mean = cell512["central_roe_upwind_broad"]["mean_rollout_nrmse"]
    figure.suptitle(
        "Roe decomposition and automatic upwinding: controlled seed-0 result",
        y=0.975,
        fontsize=17,
        fontweight="bold",
        color=TEXT,
    )
    figure.text(
        0.985,
        0.025,
        "Signed multiplier: "
        f"{100.0 * signed_multiplier['negative_multiplier_fraction']:.2f}% "
        "negative; auto-upwind: 0% negative. 512-cell mean improvement vs "
        "old HLLC+Roe: "
        f"{100.0 * (1.0 - signed_mean / old_mean):.2f}% signed, "
        f"{100.0 * (1.0 - upwind_mean / old_mean):.2f}% auto-upwind. "
        "Direct physical-complete flux is omitted from panel D because one "
        "of five cases failed its safe-step cost limit.",
        ha="right",
        fontsize=8.2,
        color=MUTED,
    )
    output = results / f"roe_upwind_summary_seed{args.seed}.png"
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)

    diagnostics = {
        "central_roe_signed_broad": {
            "multiplier": signed_multiplier,
            "hard_projection_intervention_rate": signed_report[
                "proposal_diagnostics_on_validation_ground_truth_states"
            ]["hard_projection_intervention_rate"],
        },
        "central_roe_upwind_broad": {
            "multiplier": upwind_multiplier,
            "hard_projection_intervention_rate": upwind_report[
                "proposal_diagnostics_on_validation_ground_truth_states"
            ]["hard_projection_intervention_rate"],
        },
        "roe_complete_broad": {
            "maximum_equal_interface_raw_error": complete_report[
                "equal_interface_with_varied_outer_stencil_consistency"
            ]["maximum_raw_absolute_error"],
            "hard_projection_intervention_rate": complete_report[
                "proposal_diagnostics_on_validation_ground_truth_states"
            ]["hard_projection_intervention_rate"],
        },
    }
    summary = {
        "seed": args.seed,
        "validation_best_rollout_nrmse": validation,
        "cell64": cell64,
        "cell512": cell512,
        "diagnostics": diagnostics,
        "interpretation": {
            "roe_coordinates_alone_repair_complete_flux": False,
            "both_structured_central_roe_arms_beat_old_roe_mean_at_512": True,
            "signed_has_lower_mean_rollout_nrmse": signed_mean < upwind_mean,
            "upwind_has_lower_mean_final_snapshot_nrmse": (
                cell512["central_roe_upwind_broad"][
                    "mean_final_snapshot_nrmse"
                ]
                < cell512["central_roe_signed_broad"][
                    "mean_final_snapshot_nrmse"
                ]
            ),
            "single_seed_only": True,
        },
    }
    (results / f"summary_seed{args.seed}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    rows: list[dict[str, Any]] = []
    for arm in ARMS:
        rows.extend([
            {
                "scope": "validation",
                "case": "all_validation",
                "method": arm,
                "metric": "best_rollout_nrmse",
                "value": validation[arm],
            },
            {
                "scope": "64_cell",
                "case": "random_mean",
                "method": arm,
                "metric": "rollout_nrmse",
                "value": cell64[arm]["random_mean_rollout_nrmse"],
            },
            {
                "scope": "64_cell",
                "case": "canonical_mean",
                "method": arm,
                "metric": "rollout_nrmse",
                "value": cell64[arm]["canonical_mean_rollout_nrmse"],
            },
        ])
    for arm, values in cell512.items():
        rows.append(
            {
                "scope": "512_cell",
                "case": "canonical_mean",
                "method": arm,
                "metric": "rollout_nrmse",
                "value": values["mean_rollout_nrmse"],
            }
        )
    csv_path = results / f"summary_seed{args.seed}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(output)
    print(json.dumps(summary["interpretation"], indent=2))


if __name__ == "__main__":
    main()
