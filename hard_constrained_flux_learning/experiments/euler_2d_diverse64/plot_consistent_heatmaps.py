"""Plot held-out final fields for the method-consistent 2-D comparison."""

from __future__ import annotations

import json

import numpy as np
import torch

import plot_final_heatmaps as P
import run_consistent64 as X
import run_diverse64 as R


@torch.no_grad()
def main() -> None:
    seed = 0
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    archive = R.load_npz(R.DATA / "test.npz")
    families = archive["family_names"]
    indices = P.representative_indices(families)
    times = archive["times"]
    reference = torch.from_numpy(archive["q"][indices]).float()
    roe = torch.from_numpy(archive["native_roe_coarse"][indices]).float()
    initial = reference[:, 0].to(device)
    hllc, _ = R.B.rollout_hllc(initial, times)
    signed, _ = R.B.rollout_hcfl(R.load_model(seed, device), initial, times)
    consistent, _ = R.B.rollout_hcfl(X.load_model(seed, device), initial, times)

    payload: dict[str, object] = {
        "indices": indices,
        "families": [str(families[index]) for index in indices],
        "time": float(times[-1]),
        "Reference 512→64": reference[:, -1, ..., 0].numpy(),
        "PyClaw Roe-64": roe[:, -1, ..., 0].numpy(),
        "HLLC-64": hllc[:, -1, ..., 0].numpy(),
        "HLLC + signed Roe-18": signed[:, -1, ..., 0].numpy(),
        "Central + nonnegative Roe-18": consistent[:, -1, ..., 0].numpy(),
    }
    R.RESULTS.mkdir(parents=True, exist_ok=True)
    P.plot_solution_fields(
        payload, R.RESULTS / "consistent_final_time_density_heatmaps.png"
    )
    metrics = P.plot_errors(
        payload, R.RESULTS / "consistent_final_time_density_error_heatmaps.png"
    )
    metrics["comparison"] = (
        "method-consistent central+nonnegative Roe versus old HLLC+signed Roe"
    )
    (R.RESULTS / "consistent_final_time_heatmap_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()
