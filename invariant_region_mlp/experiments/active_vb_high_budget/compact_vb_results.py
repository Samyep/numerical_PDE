"""Create Git-sized, per-repetition raw artifacts from detailed local NPZ files.

The detailed files retain every predicted z coordinate and remain available
locally.  The compact files preserve float64 value predictions, nonlinear
value corrections, and pointwise gradient-error norms, while the companion
CSV preserves every scalar diagnostic and work counter.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
RESULT_ROOT = PROJECT_ROOT / "results" / "active_vb_high_budget"


def compact_stage(stage: str, *, overwrite: bool = False) -> list[Path]:
    source = RESULT_ROOT / "repetitions" / stage
    paths = sorted(source.glob("*.npz"))
    if not paths:
        raise FileNotFoundError(f"no files under {source}")
    by_dimension: dict[int, list[Path]] = {}
    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"]))
        by_dimension.setdefault(int(metadata["dimension"]), []).append(path)

    outputs = []
    compact_root = RESULT_ROOT / "compact_repetitions"
    compact_root.mkdir(parents=True, exist_ok=True)
    for dimension, dimension_paths in sorted(by_dimension.items()):
        output = compact_root / f"{stage}_d{dimension:03d}.npz"
        if output.exists() and not overwrite:
            raise FileExistsError(f"refusing to overwrite {output}")
        records = []
        value_predictions = []
        nonlinear_corrections = []
        gradient_error_norms = []
        for path in dimension_paths:
            with np.load(path, allow_pickle=False) as data:
                prediction = np.asarray(data["prediction"], dtype=np.float64)
                truth = np.asarray(data["truth"], dtype=np.float64)
                nonlinear = np.asarray(data["nonlinear_u_correction"], dtype=np.float64)
                metadata = json.loads(str(data["metadata_json"]))
            value_predictions.append(prediction[:, 0])
            nonlinear_corrections.append(nonlinear)
            gradient_error_norms.append(np.linalg.norm(prediction[:, 1:] - truth[:, 1:], axis=1))
            records.append(
                {
                    "source": str(path.relative_to(RESULT_ROOT)).replace("\\", "/"),
                    "n": metadata["n"],
                    "M": metadata["M"],
                    "method": metadata["method"],
                    "repetition": metadata["repetition"],
                    "paired_seed": metadata["paired_seed"],
                    "wall_clock_seconds": metadata["wall_clock_seconds"],
                    "work": metadata["work"],
                    "metrics": metadata["metrics"],
                }
            )
        temporary = output.with_name(output.stem + f".tmp-{os.getpid()}.npz")
        np.savez_compressed(
            temporary,
            value_prediction=np.stack(value_predictions),
            nonlinear_u_correction=np.stack(nonlinear_corrections),
            gradient_error_l2_per_point=np.stack(gradient_error_norms),
            records_json=np.asarray(json.dumps(records, sort_keys=True, allow_nan=True)),
        )
        os.replace(temporary, output)
        outputs.append(output)
        print(f"wrote {output} ({output.stat().st_size / 2**20:.2f} MiB)", flush=True)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", default="main")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    compact_stage(args.stage, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
