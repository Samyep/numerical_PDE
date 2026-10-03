"""Plot selected nonperiodic profiles for 4/5/6-cell Euler flux models.

This is an inference-only visualization.  The four- and six-cell checkpoints
come from the symmetric-stencil experiment, while the five-cell checkpoints
are the previously converged legacy asymmetric models.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch


HERE = Path(__file__).resolve().parent
HCFL_ROOT = HERE.parents[1]
LOCAL_RUN = HCFL_ROOT / "local_run"
NONPERIODIC_EXPERIMENT = (
    HCFL_ROOT / "experiments" / "euler_1d_nonperiodic_audit"
)
for module_path in (LOCAL_RUN, NONPERIODIC_EXPERIMENT):
    if str(module_path) not in sys.path:
        sys.path.insert(0, str(module_path))

import euler_ablation_runner as base  # noqa: E402
import run_nonperiodic_audit as nonperiodic  # noqa: E402


CASES = ("lax", "collision", "contact_right_exit")
ZOOM_WINDOWS = {
    "lax": (0.40, 0.65),
    "collision": (0.42, 0.58),
    "contact_right_exit": (0.90, 1.00),
}
VARIABLES = (r"Density $\rho$", r"Velocity $u$", r"Pressure $p$")
STENCIL_LABELS = {
    4: "4-cell symmetric input",
    5: "5-cell legacy asymmetric input",
    6: "6-cell symmetric input",
}
MODEL_SPECS = {
    4: {
        "hllc_roe": (
            "dissipation_4",
            HERE / "results" / "hllc_roe_s4_converged_best_seed0.pt",
        ),
        "central_feas": (
            "central_roe_upwind_4",
            HERE
            / "results"
            / "central_nonnegative_feas_s4_converged_best_seed0.pt",
        ),
    },
    5: {
        "hllc_roe": (
            "dissipation",
            HCFL_ROOT
            / "experiments"
            / "euler_1d_convergence_audit"
            / "results"
            / "dissipation_broad_converged_best_seed0.pt",
        ),
        "central_feas": (
            "central_roe_upwind",
            HCFL_ROOT
            / "experiments"
            / "euler_1d_roe_upwind_ablation"
            / "results"
            / "central_roe_upwind_feas_broad_converged_best_seed0.pt",
        ),
    },
    6: {
        "hllc_roe": (
            "dissipation_6",
            HERE / "results" / "hllc_roe_s6_converged_best_seed0.pt",
        ),
        "central_feas": (
            "central_roe_upwind_6",
            HERE
            / "results"
            / "central_nonnegative_feas_s6_converged_best_seed0.pt",
        ),
    },
}
METHOD_STYLES = {
    "hllc_roe": {
        "label": "HLLC + Roe correction",
        "color": "#0072B2",
        "linestyle": "-",
    },
    "central_feas": {
        "label": "Central + nonnegative Roe + feasibility",
        "color": "#D55E00",
        "linestyle": "--",
    },
}


def load_model(model_name: str, checkpoint: Path) -> base.Solver:
    if not checkpoint.exists():
        raise FileNotFoundError(checkpoint)
    model = base.Solver(
        model_name,
        np.zeros(3, dtype=np.float32),
        np.ones(3, dtype=np.float32),
        width=72,
    )
    state = torch.load(checkpoint, map_location="cpu", weights_only=True)
    model.load_state_dict(state)
    model.eval()
    return model


def primitive_final(trajectory: torch.Tensor) -> np.ndarray:
    return base.primitive(trajectory)[0, -1].detach().cpu().numpy()


def local_profile(
    values: np.ndarray,
    cells: int,
    window: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    x = (np.arange(cells) + 0.5) / cells
    mask = (x >= window[0]) & (x <= window[1])
    return x[mask], values[mask]


def displayed_component(
    values: np.ndarray,
    case: str,
    variable: int,
) -> np.ndarray:
    """Expose tiny post-exit velocity/pressure errors without axis offsets."""
    if case == "contact_right_exit" and variable in (1, 2):
        return values - 1.0
    return values


def total_variation(values: np.ndarray) -> float:
    return float(np.abs(np.diff(values)).sum())


def extrema_count(values: np.ndarray, scale: float) -> int:
    if values.size < 3:
        return 0
    left = values[1:-1] - values[:-2]
    right = values[2:] - values[1:-1]
    threshold = 1.0e-3 * max(scale, 1.0e-8)
    return int(
        (
            (left * right < 0.0)
            & (np.minimum(np.abs(left), np.abs(right)) > threshold)
        ).sum()
    )


def profile_diagnostics(
    profile: np.ndarray,
    cells: int,
    window: tuple[float, float],
) -> list[dict[str, float | int]]:
    _, local = local_profile(profile, cells, window)
    rows: list[dict[str, float | int]] = []
    for variable in range(3):
        values = local[:, variable]
        scale = max(float(np.ptp(values)), float(np.max(np.abs(values))), 1.0)
        rows.append(
            {
                "local_total_variation": total_variation(values),
                "local_significant_extrema": extrema_count(values, scale),
                "local_minimum": float(values.min()),
                "local_maximum": float(values.max()),
            }
        )
    return rows


@torch.no_grad()
def generate(args: argparse.Namespace) -> None:
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    models = {
        stencil: {
            method: load_model(model_name, checkpoint)
            for method, (model_name, checkpoint) in methods.items()
        }
        for stencil, methods in MODEL_SPECS.items()
    }

    references: dict[str, torch.Tensor] = {}
    native: dict[str, torch.Tensor] = {}
    predictions: dict[int, dict[str, dict[str, torch.Tensor]]] = {
        stencil: {method: {} for method in METHOD_STYLES}
        for stencil in MODEL_SPECS
    }
    diagnostics: list[dict[str, Any]] = []

    for case in CASES:
        print(json.dumps({"stage": "reference", "case": case}), flush=True)
        references[case], _ = nonperiodic.strict_hllc_rollout(case, 2048)
        native[case], _ = nonperiodic.strict_hllc_rollout(case, 512)
        reference_profile = primitive_final(references[case])
        for variable, values in enumerate(
            profile_diagnostics(reference_profile, 2048, ZOOM_WINDOWS[case])
        ):
            diagnostics.append(
                {
                    "case": case,
                    "stencil_size": "reference_2048",
                    "method": "hllc_2048",
                    "variable": variable,
                    **values,
                }
            )
        for stencil, stencil_models in models.items():
            for method, model in stencil_models.items():
                print(
                    json.dumps(
                        {
                            "stage": "learned_rollout",
                            "case": case,
                            "stencil": stencil,
                            "method": method,
                        }
                    ),
                    flush=True,
                )
                trajectory, _ = nonperiodic.learned_rollout(model, case, 512)
                predictions[stencil][method][case] = trajectory
                profile = primitive_final(trajectory)
                for variable, values in enumerate(
                    profile_diagnostics(profile, 512, ZOOM_WINDOWS[case])
                ):
                    diagnostics.append(
                        {
                            "case": case,
                            "stencil_size": stencil,
                            "method": method,
                            "variable": variable,
                            **values,
                        }
                    )

    limits: dict[tuple[str, int], tuple[float, float]] = {}
    for case in CASES:
        window = ZOOM_WINDOWS[case]
        profiles: list[tuple[np.ndarray, int]] = [
            (primitive_final(references[case]), 2048),
            (primitive_final(native[case]), 512),
        ]
        profiles.extend(
            (primitive_final(predictions[stencil][method][case]), 512)
            for stencil in MODEL_SPECS
            for method in METHOD_STYLES
        )
        for variable in range(3):
            values = np.concatenate(
                [
                    displayed_component(
                        local_profile(profile, cells, window)[1][:, variable],
                        case,
                        variable,
                    )
                    for profile, cells in profiles
                ]
            )
            low = float(values.min())
            high = float(values.max())
            padding = 0.06 * max(high - low, 1.0e-8)
            limits[(case, variable)] = (low - padding, high + padding)

    for stencil in MODEL_SPECS:
        figure, axes = plt.subplots(3, 3, figsize=(15.4, 9.4), squeeze=False)
        for column, case in enumerate(CASES):
            window = ZOOM_WINDOWS[case]
            reference_profile = primitive_final(references[case])
            native_profile = primitive_final(native[case])
            x_reference, reference_local = local_profile(
                reference_profile, 2048, window
            )
            x_native, native_local = local_profile(native_profile, 512, window)
            for row in range(3):
                axis = axes[row, column]
                axis.plot(
                    x_reference,
                    displayed_component(reference_local[:, row], case, row),
                    color="#CBD5E1",
                    linewidth=3.0,
                    label="native HLLC-2048" if row == 0 else None,
                    zorder=1,
                )
                axis.plot(
                    x_native,
                    displayed_component(native_local[:, row], case, row),
                    color="#374151",
                    linewidth=1.1,
                    linestyle=":",
                    label="native HLLC-512" if row == 0 else None,
                    zorder=2,
                )
                for method, style in METHOD_STYLES.items():
                    profile = primitive_final(
                        predictions[stencil][method][case]
                    )
                    x_learned, local = local_profile(profile, 512, window)
                    axis.plot(
                        x_learned,
                        displayed_component(local[:, row], case, row),
                        color=style["color"],
                        linewidth=1.55,
                        linestyle=style["linestyle"],
                        label=style["label"] if row == 0 else None,
                        zorder=4,
                    )
                axis.set_xlim(*window)
                axis.set_ylim(*limits[(case, row)])
                axis.grid(True, alpha=0.2, linewidth=0.6)
                axis.spines["top"].set_visible(False)
                axis.spines["right"].set_visible(False)
                if case == "contact_right_exit":
                    axis.axvline(1.0, color="#111827", linewidth=0.7, alpha=0.5)
                if column == 0:
                    axis.set_ylabel(VARIABLES[row])
                elif case == "contact_right_exit" and row == 1:
                    axis.set_ylabel(r"Velocity error $u-1$")
                elif case == "contact_right_exit" and row == 2:
                    axis.set_ylabel(r"Pressure error $p-1$")
                if row == 0:
                    axis.set_title(
                        f"{nonperiodic.CASES[case]['display']}\n"
                        f"zoom x in [{window[0]:.2f}, {window[1]:.2f}]",
                        fontsize=11,
                        fontweight="bold",
                    )
                if row == 2:
                    axis.set_xlabel("x (cell centers joined by straight lines)")

        handles, labels = axes[0, 0].get_legend_handles_labels()
        figure.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.965),
            ncol=4,
            frameon=False,
            fontsize=9.5,
        )
        figure.suptitle(
            f"Nonperiodic transmissive smoothness check: "
            f"{STENCIL_LABELS[stencil]} (seed 0, t = 0.0252)",
            y=0.995,
            fontsize=15,
            fontweight="bold",
        )
        figure.tight_layout(rect=(0, 0, 1, 0.92))
        figure.savefig(
            output / f"selected_nonperiodic_smoothness_s{stencil}_seed0.png",
            dpi=220,
            bbox_inches="tight",
            facecolor="white",
        )
        plt.close(figure)

    payload = {
        "seed": 0,
        "training": "no training; inference-only plot",
        "boundary_condition": "transmissive nonperiodic",
        "cases": list(CASES),
        "zoom_windows": ZOOM_WINDOWS,
        "display_note": (
            "For contact_right_exit, velocity and pressure are plotted as "
            "u-1 and p-1; diagnostics remain in the original variables."
        ),
        "model_sources": {
            str(stencil): {
                method: str(checkpoint.resolve())
                for method, (_, checkpoint) in methods.items()
            }
            for stencil, methods in MODEL_SPECS.items()
        },
        "local_profile_diagnostics": diagnostics,
    }
    (output / "selected_nonperiodic_smoothness_seed0.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


def left_exit_displayed_component(
    values: np.ndarray,
    variable: int,
) -> np.ndarray:
    """Plot post-exit velocity/pressure as errors from (-1, 1)."""
    if variable == 1:
        return values + 1.0
    if variable == 2:
        return values - 1.0
    return values


@torch.no_grad()
def generate_contact_left(args: argparse.Namespace) -> None:
    """Compare 4/5/6 cells on the left-going boundary-interaction case."""
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    case = "contact_left_exit"
    window = (0.0, 0.10)
    models = {
        stencil: {
            method: load_model(model_name, checkpoint)
            for method, (model_name, checkpoint) in methods.items()
        }
        for stencil, methods in MODEL_SPECS.items()
    }

    print(json.dumps({"stage": "reference", "case": case}), flush=True)
    reference, _ = nonperiodic.strict_hllc_rollout(case, 2048)
    native, _ = nonperiodic.strict_hllc_rollout(case, 512)
    predictions: dict[int, dict[str, torch.Tensor]] = {
        stencil: {} for stencil in MODEL_SPECS
    }
    diagnostics: list[dict[str, Any]] = []
    reference_profile = primitive_final(reference)
    for variable, values in enumerate(
        profile_diagnostics(reference_profile, 2048, window)
    ):
        diagnostics.append(
            {
                "stencil_size": "reference_2048",
                "method": "hllc_2048",
                "variable": variable,
                **values,
            }
        )

    for stencil, stencil_models in models.items():
        for method, model in stencil_models.items():
            print(
                json.dumps(
                    {
                        "stage": "learned_rollout",
                        "case": case,
                        "stencil": stencil,
                        "method": method,
                    }
                ),
                flush=True,
            )
            trajectory, _ = nonperiodic.learned_rollout(model, case, 512)
            predictions[stencil][method] = trajectory
            profile = primitive_final(trajectory)
            for variable, values in enumerate(
                profile_diagnostics(profile, 512, window)
            ):
                diagnostics.append(
                    {
                        "stencil_size": stencil,
                        "method": method,
                        "variable": variable,
                        **values,
                    }
                )

    native_profile = primitive_final(native)
    profiles: list[tuple[np.ndarray, int]] = [
        (reference_profile, 2048),
        (native_profile, 512),
    ]
    profiles.extend(
        (primitive_final(predictions[stencil][method]), 512)
        for stencil in MODEL_SPECS
        for method in METHOD_STYLES
    )
    limits: dict[int, tuple[float, float]] = {}
    for variable in range(3):
        values = np.concatenate(
            [
                left_exit_displayed_component(
                    local_profile(profile, cells, window)[1][:, variable],
                    variable,
                )
                for profile, cells in profiles
            ]
        )
        low = float(values.min())
        high = float(values.max())
        padding = 0.06 * max(high - low, 1.0e-8)
        limits[variable] = (low - padding, high + padding)

    figure, axes = plt.subplots(3, 3, figsize=(15.4, 9.4), squeeze=False)
    x_reference, reference_local = local_profile(
        reference_profile, 2048, window
    )
    x_native, native_local = local_profile(native_profile, 512, window)
    row_labels = (
        r"Density $\rho$",
        r"Velocity error $u+1$",
        r"Pressure error $p-1$",
    )
    for column, stencil in enumerate((4, 5, 6)):
        for row in range(3):
            axis = axes[row, column]
            axis.plot(
                x_reference,
                left_exit_displayed_component(
                    reference_local[:, row], row
                ),
                color="#CBD5E1",
                linewidth=3.0,
                label="native HLLC-2048" if row == 0 else None,
                zorder=1,
            )
            axis.plot(
                x_native,
                left_exit_displayed_component(native_local[:, row], row),
                color="#374151",
                linewidth=1.1,
                linestyle=":",
                label="native HLLC-512" if row == 0 else None,
                zorder=2,
            )
            for method, style in METHOD_STYLES.items():
                profile = primitive_final(predictions[stencil][method])
                x_learned, local = local_profile(profile, 512, window)
                axis.plot(
                    x_learned,
                    left_exit_displayed_component(local[:, row], row),
                    color=style["color"],
                    linewidth=1.55,
                    linestyle=style["linestyle"],
                    label=style["label"] if row == 0 else None,
                    zorder=4,
                )
            axis.set_xlim(*window)
            axis.set_ylim(*limits[row])
            axis.axvline(0.0, color="#111827", linewidth=0.7, alpha=0.5)
            axis.grid(True, alpha=0.2, linewidth=0.6)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
            if column == 0:
                axis.set_ylabel(row_labels[row])
            if row == 0:
                axis.set_title(
                    STENCIL_LABELS[stencil],
                    fontsize=11,
                    fontweight="bold",
                )
            if row == 2:
                axis.set_xlabel("x (cell centers joined by straight lines)")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=4,
        frameon=False,
        fontsize=9.5,
    )
    figure.suptitle(
        "Contact exits left: nonperiodic transmissive 4/5/6-cell "
        "smoothness check (seed 0, t = 0.0252, zoom x in [0.00, 0.10])",
        y=0.995,
        fontsize=14.5,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    figure.savefig(
        output / "contact_left_exit_smoothness_456_seed0.png",
        dpi=220,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(figure)

    payload = {
        "seed": 0,
        "training": "no training; inference-only plot",
        "boundary_condition": "transmissive nonperiodic",
        "case": case,
        "zoom_window": window,
        "display_note": (
            "Velocity and pressure are plotted as u+1 and p-1; local "
            "diagnostics remain in the original primitive variables."
        ),
        "model_sources": {
            str(stencil): {
                method: str(checkpoint.resolve())
                for method, (_, checkpoint) in methods.items()
            }
            for stencil, methods in MODEL_SPECS.items()
        },
        "local_profile_diagnostics": diagnostics,
    }
    (output / "contact_left_exit_smoothness_456_seed0.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=HERE / "results")
    parser.add_argument("--contact-left-only", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    if arguments.contact_left_only:
        generate_contact_left(arguments)
    else:
        generate(arguments)
