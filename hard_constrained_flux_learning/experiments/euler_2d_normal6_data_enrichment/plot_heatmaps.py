"""Plot held-out final-density fields for the data-only ablation."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

import torch


HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "euler_2d_normal6_fixed_transverse"
DIVERSE = HERE.parent / "euler_2d_diverse64"
sys.path.insert(0, str(DIVERSE))
import plot_final_heatmaps as P  # noqa: E402
import run_diverse64 as R  # noqa: E402


RUNNER_SPEC = importlib.util.spec_from_file_location(
    "normal6_data_enrichment_heatmap_runner", BASE / "run_experiment.py"
)
if RUNNER_SPEC is None or RUNNER_SPEC.loader is None:
    raise ImportError("Could not load the fixed-transverse experiment runner")
E = importlib.util.module_from_spec(RUNNER_SPEC)
sys.modules[RUNNER_SPEC.name] = E
RUNNER_SPEC.loader.exec_module(E)


def load_from(results: Path, seed: int, device: torch.device) -> torch.nn.Module:
    """Load the same method from an explicitly selected result directory."""
    E.RESULTS = results
    return E.load_model(seed, device)


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

    original_model = load_from(BASE / "results", args.seed, device)
    original, _ = R.B.rollout_hcfl(original_model, initial, times)
    enriched_model = load_from(HERE / "results", args.seed, device)
    enriched, _ = R.B.rollout_hcfl(enriched_model, initial, times)

    payload: dict[str, object] = {
        "indices": indices,
        "families": [str(families[index]) for index in indices],
        "time": float(times[-1]),
        "Reference 512→64": reference[:, -1, ..., 0].numpy(),
        "PyClaw Roe": roe[:, -1, ..., 0].numpy(),
        "HCFL original data": original[:, -1, ..., 0].numpy(),
        "HCFL enriched data": enriched[:, -1, ..., 0].numpy(),
    }
    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    P.plot_solution_fields(payload, results / "final_time_density_heatmaps.png")
    metrics = P.plot_errors(
        payload, results / "final_time_density_error_heatmaps.png"
    )
    metrics.update(
        {
            "checkpoint_seed": args.seed,
            "comparison": "training-data-only enrichment",
            "original_training_trajectories": 192,
            "enriched_training_trajectories": 384,
        }
    )
    (results / "final_time_heatmap_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()
