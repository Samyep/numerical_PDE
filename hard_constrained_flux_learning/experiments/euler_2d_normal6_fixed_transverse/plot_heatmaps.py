"""Plot deterministic held-out final-density fields for the transverse audit."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

import torch


HERE = Path(__file__).resolve().parent
DIVERSE = HERE.parent / "euler_2d_diverse64"
sys.path.insert(0, str(DIVERSE))
import plot_final_heatmaps as P  # noqa: E402
import run_consistent64 as S  # noqa: E402
import run_diverse64 as R  # noqa: E402


RUNNER_SPEC = importlib.util.spec_from_file_location(
    "normal6_fixed_transverse_runner", HERE / "run_experiment.py"
)
if RUNNER_SPEC is None or RUNNER_SPEC.loader is None:
    raise ImportError("Could not load the fixed-transverse experiment runner")
X = importlib.util.module_from_spec(RUNNER_SPEC)
RUNNER_SPEC.loader.exec_module(X)


@torch.no_grad()
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    data = args.data_dir.resolve()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    archive = R.load_npz(data / "test.npz")
    families = archive["family_names"]
    indices = P.representative_indices(families)
    times = archive["times"]
    reference = torch.from_numpy(archive["q"][indices]).float()
    roe = torch.from_numpy(archive["native_roe_coarse"][indices]).float()
    initial = reference[:, 0].to(device)
    learned_transverse, _ = R.B.rollout_hcfl(
        S.load_model(args.seed, device), initial, times
    )
    fixed_transverse, _ = R.B.rollout_hcfl(
        X.load_model(args.seed, device), initial, times
    )

    payload: dict[str, object] = {
        "indices": indices,
        "families": [str(families[index]) for index in indices],
        "time": float(times[-1]),
        "Reference 512→64": reference[:, -1, ..., 0].numpy(),
        "PyClaw Roe": roe[:, -1, ..., 0].numpy(),
        "HCFL learned transverse": learned_transverse[:, -1, ..., 0].numpy(),
        "HCFL fixed transverse": fixed_transverse[:, -1, ..., 0].numpy(),
    }
    X.RESULTS.mkdir(parents=True, exist_ok=True)
    P.plot_solution_fields(payload, X.RESULTS / "final_time_density_heatmaps.png")
    metrics = P.plot_errors(
        payload, X.RESULTS / "final_time_density_error_heatmaps.png"
    )
    metrics.update(
        {
            "checkpoint_seed": args.seed,
            "learned_stencil_cells": 6,
            "fixed_transverse_trainable_parameters": 0,
        }
    )
    (X.RESULTS / "final_time_heatmap_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()
