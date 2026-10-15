"""Stage-A correctness/throughput calibration using validation objectives only."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import torch

from .common import (
    BUDGETS,
    CHECKPOINT_DIR,
    RESULTS_DIR,
    driver_for,
    ensure_dirs,
    evaluation_data,
    write_csv,
)
from .deep_picard import train_deep_picard
from .full_history_mlp import predict_full_history_mlp
from .neural_baselines import (
    load_sdgd_pinn,
    train_deep_bsde_timed,
    train_sdgd_pinn_timed,
)
from .scasml import predict_scasml_level1
from .tangent_baseline import train_tangent_method


BENCHMARK = "nonseparable_d100_k2"
SEED = 97


def _continuous_rows(method: str, result: dict[str, Any], trial: str) -> list[dict[str, Any]]:
    checkpoint = result["checkpoints"][-1]
    steps = int(checkpoint["step"])
    elapsed = float(checkpoint["training_seconds"])
    rate = steps / elapsed
    return [
        {
            "method": method,
            "trial": trial,
            "calibration_seconds": elapsed,
            "calibration_units": steps,
            "units_name": "optimizer_steps",
            "units_per_second": rate,
            **{f"estimated_units_B{index + 1}": round(rate * budget) for index, budget in enumerate(BUDGETS)},
            "validation_objective": checkpoint["validation_objective"],
        }
    ]


def main() -> None:
    ensure_dirs()
    if not torch.cuda.is_available():
        raise RuntimeError("GH200 CUDA allocation required for timing calibration")
    device = torch.device("cuda")
    timing: list[dict[str, Any]] = []
    trials: list[dict[str, Any]] = []

    deep_a = train_deep_bsde_timed(
        BENCHMARK,
        SEED,
        device,
        calibration_only_seconds=30.0,
        width=192,
        depth=3,
        time_steps=16,
        batch_size=128,
        learning_rate=5e-4,
        trial_tag="_A",
    )
    deep_b = train_deep_bsde_timed(
        BENCHMARK,
        SEED,
        device,
        calibration_only_seconds=30.0,
        width=256,
        depth=3,
        time_steps=20,
        batch_size=64,
        learning_rate=1e-3,
        trial_tag="_B",
    )
    deep_results = {"A": deep_a, "B": deep_b}
    deep_selected = min(
        deep_results,
        key=lambda key: deep_results[key]["checkpoints"][-1]["validation_objective"],
    )
    for label, result in deep_results.items():
        timing.extend(_continuous_rows("deep_bsde", result, label))
        checkpoint = result["checkpoints"][-1]
        trials.append(
            {
                "method": "deep_bsde",
                "trial": label,
                "configuration": json.dumps(result["configuration"], sort_keys=True),
                "validation_metric": "frozen terminal Huber",
                "validation_objective": checkpoint["validation_objective"],
                "seconds": checkpoint["training_seconds"],
                "selected": label == deep_selected,
                "test_metric_seen": False,
            }
        )

    pinn_a = train_sdgd_pinn_timed(
        BENCHMARK,
        SEED,
        device,
        calibration_only_seconds=30.0,
        width=192,
        depth=3,
        batch_size=128,
        coordinate_budget=16,
        learning_rate=5e-4,
        trial_tag="_A",
    )
    pinn_b = train_sdgd_pinn_timed(
        BENCHMARK,
        SEED,
        device,
        calibration_only_seconds=30.0,
        width=256,
        depth=4,
        batch_size=64,
        coordinate_budget=8,
        learning_rate=1e-3,
        trial_tag="_B",
    )
    pinn_results = {"A": pinn_a, "B": pinn_b}
    pinn_selected = min(
        pinn_results,
        key=lambda key: pinn_results[key]["checkpoints"][-1]["validation_sampled_residual_rms"],
    )
    for label, result in pinn_results.items():
        timing.extend(_continuous_rows("sdgd_pinn", result, label))
        checkpoint = result["checkpoints"][-1]
        trials.append(
            {
                "method": "sdgd_pinn",
                "trial": label,
                "configuration": json.dumps(result["configuration"], sort_keys=True),
                "validation_metric": "frozen sampled residual RMS",
                "validation_objective": checkpoint["validation_sampled_residual_rms"],
                "seconds": checkpoint["training_seconds"],
                "selected": label == pinn_selected,
                "test_metric_seen": False,
            }
        )

    dpi_a = train_deep_picard(
        BENCHMARK,
        SEED,
        device,
        width=512,
        depth=4,
        training_points=4096,
        paths=512,
        epochs=16,
        batch_size=512,
        learning_rate=1e-3,
        calibration_rounds=1,
        trial_tag="_A",
    )
    dpi_b = train_deep_picard(
        BENCHMARK,
        SEED,
        device,
        width=256,
        depth=3,
        training_points=2048,
        paths=128,
        epochs=20,
        batch_size=512,
        learning_rate=5e-4,
        calibration_rounds=1,
        trial_tag="_B",
    )
    dpi_results = {"A": dpi_a, "B": dpi_b}
    # Round one uses the zero Picard source, so its throughput does not measure
    # the network-sourced rounds that dominate the main run.  Trial A has the
    # lower one-round label loss but cannot provide a useful multi-round ladder
    # inside B1--B5; retain it as an attempted official-style configuration and
    # select the calibrated feasible trial B without consulting test error.
    dpi_selected = "B"
    for label, result in dpi_results.items():
        row = result["history"][-1]
        seconds = row["cumulative_training_seconds"]
        timing.append(
            {
                "method": "deep_picard",
                "trial": label,
                "calibration_seconds": seconds,
                "calibration_units": 1,
                "units_name": "complete_rounds",
                "units_per_second": 1.0 / seconds,
                **{f"estimated_units_B{index + 1}": None for index, _ in enumerate(BUDGETS)},
                "validation_objective": row["validation_normalized_loss"],
                "calibration_scope": "zero-source round one only; later complete rounds are timed directly",
            }
        )
        trials.append(
            {
                "method": "deep_picard",
                "trial": label,
                "configuration": json.dumps(result["configuration"], sort_keys=True),
                "validation_metric": "held-out normalized Picard-label loss",
                "validation_objective": row["validation_normalized_loss"],
                "seconds": seconds,
                "selected": label == dpi_selected,
                "test_metric_seen": False,
            }
        )

    tangent = train_tangent_method("tangent_actor", BENCHMARK, SEED, device, iterations=1)
    phase = tangent["phases"][0]
    timing.append(
        {
            "method": "tangent_actor",
            "trial": "predeclared_TA_A",
            "calibration_seconds": phase["cumulative_training_seconds"],
            "calibration_units": 1,
            "units_name": "complete_actor_phases",
            "units_per_second": 1.0 / phase["cumulative_training_seconds"],
            **{f"estimated_units_B{index + 1}": None for index, _ in enumerate(BUDGETS)},
            "validation_objective": phase["training"]["best_validation_normalized_loss"],
            "calibration_scope": "constant-source phase one only; a1--a3 are timed as complete phases",
        }
    )
    trials.append(
        {
            "method": "tangent_actor",
            "trial": "predeclared_TA_A",
            "configuration": json.dumps(tangent["configuration"], sort_keys=True),
            "validation_metric": "held-out normalized tangent-label loss",
            "validation_objective": phase["training"]["best_validation_normalized_loss"],
            "seconds": phase["cumulative_training_seconds"],
            "selected": True,
            "test_metric_seen": False,
        }
    )

    validation = evaluation_data(BENCHMARK)
    vx = validation["x"][validation["validation_indices"][:64]]
    for n, M in ((1, 64), (2, 8), (3, 4)):
        _, metadata = predict_full_history_mlp(
            BENCHMARK, vx, device, n=n, M=M, seed=SEED
        )
        timing.append(
            {
                "method": "full_history_mlp",
                "trial": f"n{n}_M{M}",
                "calibration_seconds": metadata["inference_seconds"],
                "calibration_units": len(vx),
                "units_name": "validation_state_predictions",
                "units_per_second": len(vx) / metadata["inference_seconds"],
                "estimated_960_state_seconds": metadata["inference_seconds"] * 960 / len(vx),
                "validation_objective": None,
            }
        )

    selected_pinn = pinn_results[pinn_selected]["checkpoints"][-1]["checkpoint"]
    pinn_model, _ = load_sdgd_pinn(Path(selected_pinn), device)
    _, _, scasml_meta = predict_scasml_level1(
        pinn_model,
        BENCHMARK,
        vx,
        device,
        paths=8,
        probes=4,
        seed=SEED,
    )
    timing.append(
        {
            "method": "scasml_l1",
            "trial": f"pinn_{pinn_selected}_paths8_probes4",
            "calibration_seconds": scasml_meta["inference_seconds"],
            "calibration_units": len(vx),
            "units_name": "validation_state_predictions",
            "units_per_second": len(vx) / scasml_meta["inference_seconds"],
            "estimated_960_state_seconds": scasml_meta["inference_seconds"] * 960 / len(vx),
            "validation_objective": scasml_meta["mean_mc_se_across_workload"],
        }
    )

    write_csv(RESULTS_DIR / "timing_calibration.csv", timing)
    write_csv(RESULTS_DIR / "hyperparameter_trials.csv", trials)
    selections = {
        "deep_bsde": deep_selected,
        "sdgd_pinn": pinn_selected,
        "deep_picard": dpi_selected,
        "selection_basis": "independent reference-free validation only",
        "test_metrics_seen": False,
    }
    (RESULTS_DIR / "calibration_selection.json").write_text(
        json.dumps(selections, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(selections, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
