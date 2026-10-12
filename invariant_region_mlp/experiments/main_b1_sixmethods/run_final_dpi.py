"""Archive the preliminary DPI pass and resumably run the frozen final profile."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
from typing import Any

from .aggregate import aggregate
from .common import CHECKPOINT_DIR, PACKAGE_DIR, RESULTS_DIR, load_config
from .dpi_adapter import dpi_checkpoint_path, dpi_training_signature
from .equation import B1Equation
from .runner import run_method


PROFILE_PATH = PACKAGE_DIR / "dpi_final_config.json"
ARCHIVE_ROOT = RESULTS_DIR / "superseded" / "dpi_initial_stage_b"


def _load_profile() -> dict[str, Any]:
    with PROFILE_PATH.open() as handle:
        return json.load(handle)


def _matches_final(path: Path, signature: dict[str, Any], seeds: int) -> bool:
    if not path.is_file():
        return False
    with path.open() as handle:
        payload = json.load(handle)
    shared = payload.get("shared", {})
    return all(
        shared.get(f"seed_{seed}", {}).get("training_config") == signature
        for seed in range(seeds)
    )


def _archive_once(source: Path, relative: Path) -> None:
    destination = ARCHIVE_ROOT / relative
    if source.is_file() and not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def main() -> None:
    profile = _load_profile()["resource_adapted_final_configuration"]
    config = load_config()
    actual = config["dpi"]
    for key, expected in profile.items():
        if actual.get(key) != expected:
            raise RuntimeError(
                f"config.json DPI setting {key!r} is {actual.get(key)!r}, expected {expected!r}"
            )

    signature = dpi_training_signature(actual, "stage_b")
    raw_root = RESULTS_DIR / "raw" / "stage_b"
    seeds = int(config["stage_b"]["neural_seeds"])
    scale = float(config["pde"]["frozen_coefficient_scale"])
    for d in config["pde"]["dimensions"]:
        for k in config["pde"]["effective_dimensions"]:
            output = raw_root / f"dpi_d{d}_k{k}.json"
            prediction = raw_root / f"dpi_d{d}_k{k}_predictions.npz"
            if _matches_final(output, signature, seeds):
                print(f"final DPI already complete: d={d}, k={k}")
                continue

            _archive_once(output, Path("raw") / output.name)
            _archive_once(prediction, Path("raw") / prediction.name)
            equation = B1Equation(d=d, k=k, coefficient_scale=scale)
            for seed in range(seeds):
                checkpoint = dpi_checkpoint_path(
                    CHECKPOINT_DIR, equation, seed, "stage_b"
                )
                _archive_once(checkpoint, Path("checkpoints") / checkpoint.name)

            output.unlink(missing_ok=True)
            prediction.unlink(missing_ok=True)
            run_method("stage_b", "dpi", d, k, force=False)

    aggregate("stage_b")


if __name__ == "__main__":
    main()
