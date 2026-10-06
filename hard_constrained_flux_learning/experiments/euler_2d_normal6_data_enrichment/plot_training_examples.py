"""Plot deterministic examples from the added training-only families."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
SHARDS = HERE / "data" / "shards"
RESULTS = HERE / "results"
FAMILIES = (
    "crossed_riemann",
    "corrugated_contact",
    "shock_contact_triple",
    "shock_vortex",
    "double_blast",
    "four_way_collision",
)
DISPLAY = {
    "crossed_riemann": "Crossed Riemann",
    "corrugated_contact": "Corrugated contact",
    "shock_contact_triple": "Shock-contact triple",
    "shock_vortex": "Shock-vortex",
    "double_blast": "Double blast",
    "four_way_collision": "Four-way collision",
}


def pressure(state: np.ndarray) -> np.ndarray:
    rho = state[..., 0]
    return 0.4 * (
        state[..., 3]
        - 0.5 * (state[..., 1] ** 2 + state[..., 2] ** 2) / rho
    )


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(
        len(FAMILIES), 4, figsize=(11.0, 15.0), constrained_layout=True
    )
    metrics: list[dict[str, float | str]] = []
    for row, family in enumerate(FAMILIES):
        with np.load(SHARDS / f"{family}.npz", allow_pickle=False) as archive:
            initial = archive["q"][0, 0]
            final = archive["q"][0, -1]
            final_time = float(archive["times"][-1])
        fields = (
            (initial[..., 0], r"initial $\rho$"),
            (final[..., 0], rf"$t={final_time:.2f}$ $\rho$"),
            (pressure(initial), r"initial $p$"),
            (pressure(final), rf"$t={final_time:.2f}$ $p$"),
        )
        for column, (field, title) in enumerate(fields):
            axis = axes[row, column]
            image = axis.imshow(
                field,
                origin="lower",
                extent=(0.0, 1.0, 0.0, 1.0),
                interpolation="nearest",
                cmap="viridis",
            )
            axis.set_xticks([])
            axis.set_yticks([])
            if row == 0:
                axis.set_title(title)
            if column == 0:
                axis.set_ylabel(DISPLAY[family])
            figure.colorbar(image, ax=axis, fraction=0.046, pad=0.02)
        metrics.append(
            {
                "family": family,
                "case_index": 0,
                "final_time": final_time,
                "minimum_density": float(final[..., 0].min()),
                "minimum_pressure": float(pressure(final).min()),
                "maximum_density": float(final[..., 0].max()),
                "maximum_pressure": float(pressure(final).max()),
            }
        )
    figure.suptitle(
        "Added 512→64 training trajectories (first deterministic case per family)"
    )
    figure.savefig(RESULTS / "enriched_training_examples.png", dpi=220)
    plt.close(figure)
    (RESULTS / "enriched_training_examples.json").write_text(
        json.dumps(metrics, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()
