"""Make the requested real-GH200 compute-time vs RMSE curve from FH result JSONs.

Never draws fictitious FH points, never interpolates/extrapolates an
unobserved target-crossing. Requires original frozen 100D bank and Actor row.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .runner import D, DEPTH, M_GRID, TEST_POINTS, _array_hash, load_comparison, sha_file, write_json_atomic


def read_points(comparison: dict, folder: Path):
    points, incomplete = [], []
    for M in M_GRID:
        stem = f"ou_gmm_d100_fh_n{DEPTH}_M{M}_seed0"
        path = folder / (stem + ".json")
        if not path.exists():
            continue
        row = json.loads(path.read_text())
        if row.get("status") != "COMPLETE":
            incomplete.append({
                "M": M, "status": row.get("status", "UNKNOWN"),
                "completed_states": row.get("completed_states", 0),
                "gpu_seconds": row.get("gpu_seconds", 0.0),
            })
            continue
        if row.get("depth") != DEPTH or row.get("states") != TEST_POINTS:
            raise RuntimeError(f"Invalid configuration/shape: {path}")
        a = row.get("actor", {})
        if (a.get("bank_array_sha256") != comparison["bank_array_sha256"]
                or abs(a.get("actor_rmse", float("inf")) - comparison["actor_rmse"]) > 1e-9):
            raise RuntimeError(f"Frozen comparison mismatch: {path}")
        if a.get("truth_array_sha256") != comparison["truth_array_sha256"]:
            raise RuntimeError(f"Reference changed across M: {path}")

        npz_path = Path(row["prediction_file"])
        if not npz_path.exists():
            npz_path = folder / (stem + ".npz")
        if not npz_path.exists() or sha_file(npz_path) != row["prediction_sha256"]:
            raise RuntimeError(f"Saved prediction SHA mismatch: {npz_path}")
        with np.load(npz_path, allow_pickle=False) as data:
            t, x = np.array(data["t"]), np.array(data["x"])
            y, v = np.array(data["truth"]), np.array(data["prediction"])
        if t.shape != (TEST_POINTS,) or x.shape != (TEST_POINTS, D):
            raise RuntimeError(f"Bad saved prediction shape: {npz_path}")
        if (_array_hash(t, x) != comparison["bank_array_sha256"]
            or not np.array_equal(y, comparison["truth"])):
            raise RuntimeError(f"Actual paired test arrays differ: {npz_path}")
        rmse = float(np.sqrt(np.mean((v-y)**2)))
        if abs(rmse-float(row["value_rmse"]))>1e-10:
            raise RuntimeError(f"RMSE changed from measured saved prediction: {npz_path}")
        point = {
            "M": M, "n": DEPTH, "states": TEST_POINTS,
            "gpu_seconds": float(row["gpu_seconds"]),
            "rmse": rmse, "relative_l2": float(row["relative_l2"]),
            "status": "COMPLETE"
        }
        if not math.isfinite(point["gpu_seconds"]) or point["gpu_seconds"] <= 0:
            raise RuntimeError(f"Invalid measured runtime: {path}")
        points.append(point)
    return sorted(points, key=lambda p:p["gpu_seconds"]), incomplete


def make_plot(comparison: dict, points: list[dict], folder: Path):
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.set_xscale("log")
    ax.set_yscale("log")
    target_rmse = comparison["actor_rmse"]
    target_seconds = comparison["actor_total_seconds"]
    ax.axhline(target_rmse, color="#AD5C26", ls="--", lw=1.3,
               label=f"Tangent-Actor target RMSE {target_rmse:.4g}")
    ax.plot(target_seconds, target_rmse, color="#AD5C26", marker="D", ms=7,
            label=f"Tangent-Actor measured total {target_seconds:.0f} s")
    if points:
        ax.plot([p["gpu_seconds"] for p in points],
                [p["rmse"] for p in points], marker="o", lw=1.9,
                color="#3063A2", label=f"FH-MLP (n={DEPTH}, one seed)")
        for p in points:
            ax.annotate(f"M={p['M']}", (p["gpu_seconds"], p["rmse"]),
                        xytext=(5, 6), textcoords="offset points", fontsize=8)
    hits = [p for p in points if p["rmse"] <= target_rmse]
    if hits:
        candidate = min(hits, key=lambda p:p["gpu_seconds"])
        ax.axvline(candidate["gpu_seconds"], color="#3B8B70", ls=":", lw=1.2,
                   label=f"First observed FH target hit: {candidate['gpu_seconds']:.1f} s")
    ax.set_xlabel("Measured GH200 inference time for 960 points (seconds)")
    ax.set_ylabel("Whole-function Value RMSE (lower is better)")
    ax.set_title("OU/GMM HJB, 100D: Full-History MLP accuracy–compute")
    ax.grid(True, which="both", alpha=0.23)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(folder / "ou_gmm_fh_time_accuracy.png", dpi=200)
    fig.savefig(folder / "ou_gmm_fh_time_accuracy.pdf")
    plt.close(fig)
    return min(hits, key=lambda p:p["gpu_seconds"]) if hits else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("invariant_region_mlp/results/ou_gmm_fh_frontier"))
    args = parser.parse_args()
    comp = load_comparison(args.source_root.resolve())
    folder = args.out_dir
    folder.mkdir(parents=True, exist_ok=True)
    points, incomplete = read_points(comp, folder)
    if not points:
        print("NO_MEASURED_FH_DATA: no curve fabricated. Run the GH200 jobs first.")
        return

    hit = make_plot(comp, points, folder)
    report = {
        "benchmark": "ou_gmm_d100", "test_points": TEST_POINTS, "n_fixed": DEPTH,
        "single_sweep_parameter": "M", "prespecified_grid": list(M_GRID),
        "actor_rmse": comp["actor_rmse"],
        "actor_total_seconds": comp["actor_total_seconds"],
        "MLP_points": points, "incomplete_cells": incomplete,
        "target_status": "REACHED_IN_MEASURED_GRID" if hit else "NOT_REACHED_IN_MEASURED_GRID",
        "FH_first_observed_target_seconds": hit["gpu_seconds"] if hit else None,
        "FH_first_observed_target_M": hit["M"] if hit else None,
        "actor_vs_FH_time_ratio": (comp["actor_total_seconds"]/hit["gpu_seconds"]
                                   if hit else None),
        "interpretation": "First measured threshold-reaching completed fixed-grid M, not a proof of optimal M/time; no extrapolated crossing."
    }
    write_json_atomic(folder / "ou_gmm_fh_time_accuracy_summary.json", report)
    csv_path = folder / "ou_gmm_fh_time_accuracy.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "n", "M", "states", "gpu_seconds", "rmse", "relative_l2", "status"])
        writer.writeheader()
        writer.writerows(points)
    print(json.dumps({
        "status":report["target_status"], "completed_points":len(points),
        "FH_first_observed_target_seconds":report["FH_first_observed_target_seconds"],
        "actor_seconds":report["actor_total_seconds"],
        "png":str(folder/"ou_gmm_fh_time_accuracy.png"),
        "pdf":str(folder/"ou_gmm_fh_time_accuracy.pdf")
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
