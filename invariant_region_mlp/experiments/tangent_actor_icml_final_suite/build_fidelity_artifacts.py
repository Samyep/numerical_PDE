"""Build convergence evidence for the independent baseline-fidelity review."""

from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd

from .protocol import RESULTS_DIR


SOURCE_ROOT = Path(__file__).resolve().parents[3]
ORIGINAL = SOURCE_ROOT / "invariant_region_mlp/results/original_pde_benchmark_comparison"
FIGURE_DIR = RESULTS_DIR / "fidelity_figures"


def _event_rows(path: Path, tags: Iterable[str], method: str, anchor: str) -> list[dict[str, object]]:
    from tensorboard.backend.event_processing.event_accumulator import EventAccumulator

    accumulator = EventAccumulator(str(path)); accumulator.Reload()
    rows = []
    for tag in tags:
        grouped: dict[int, list[object]] = defaultdict(list)
        if tag not in accumulator.Tags().get("scalars", []):
            continue
        for item in accumulator.Scalars(tag):
            grouped[int(item.step)].append(item)
        for step, candidates in sorted(grouped.items()):
            # DPI logs multiple validation callbacks inside each Picard round.
            # The last callback is the fully trained network for that round.
            item = candidates[-1]
            rows.append({
                "method": method, "anchor": anchor, "step": step,
                "metric": tag, "value": float(item.value),
                "wall_time": float(item.wall_time), "source": str(path),
            })
    return rows


def _history_rows(path: Path, method: str, anchor: str) -> list[dict[str, object]]:
    if not path.exists():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = []
    for item in payload.get("history", []):
        step = int(item.get("epoch", item.get("step", 0)))
        for key in (
            "reporting_relative_error", "reporting_relative_l2",
            "reported_residual_mse", "validation_terminal_clipped_mse",
        ):
            value = item.get(key)
            if value is not None:
                rows.append({
                    "method": method, "anchor": anchor, "step": step,
                    "metric": key, "value": float(value),
                    "wall_time": None, "source": str(path),
                })
    return rows


def _latest(pattern: str) -> Path | None:
    paths = sorted(Path().glob(pattern))
    return paths[-1] if paths else None


def build() -> list[Path]:
    rows: list[dict[str, object]] = []
    dpi_event = ORIGINAL / "official_work/dpi_runs/dpi_s0_official/events.out.tfevents.1792202615.gh040.hsn.cm.delta.internal.ncsa.edu"
    pinn_events = sorted((ORIGINAL / "official_work/dpi_inputs/pinn").glob(
        "OUProcessEquation_dim=100_T=1.0_v=16_tw=0.0_lr=0.001_*/events.out.tfevents.*"
    ))
    if dpi_event.exists():
        rows += _event_rows(dpi_event, ("rMAE", "rMAEg"), "DPI", "OU/GMM 100D")
    if pinn_events:
        rows += _event_rows(pinn_events[-1], ("rMAE", "rMAE_gradient", "loss_total"), "PINN-HTE", "OU/GMM 100D")
    for family, label in (("hjb_log", "HJB-Log 100D"), ("allen_cahn", "Allen-Cahn 100D")):
        rows += _history_rows(
            ORIGINAL / f"raw/deep_bsde_{family}_d100_s0.json",
            "Deep BSDE", label,
        )
    sdgd = RESULTS_DIR / "raw/sdgd_pinn_hjb_log_d250_hjb_adversarial_official_s0.json"
    rows += _history_rows(sdgd, "SDGD-PINN", "HJB-Log 250D")

    frame = pd.DataFrame(rows)
    output_csv = RESULTS_DIR / "BASELINE_CONVERGENCE.csv"
    frame.to_csv(output_csv, index=False)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.2))
    panels = (
        (axes[0, 0], "DPI", "rMAE", "Picard round", "relative MAE"),
        (axes[0, 1], "PINN-HTE", "rMAE", "optimizer step", "relative MAE"),
        (axes[1, 0], "Deep BSDE", "validation_terminal_clipped_mse", "optimizer step", "validation terminal MSE"),
        (axes[1, 1], "SDGD-PINN", "reporting_relative_error", "optimizer epoch", "relative L2"),
    )
    for axis, method, metric, xlabel, ylabel in panels:
        subset = frame[(frame["method"] == method) & (frame["metric"] == metric)] if not frame.empty else frame
        if subset.empty:
            axis.text(0.5, 0.5, "pending", ha="center", va="center", transform=axis.transAxes)
        else:
            for anchor, group in subset.groupby("anchor"):
                group = group.sort_values("step")
                axis.plot(group["step"], group["value"], label=anchor, linewidth=1.5)
            axis.set_yscale("log")
            axis.legend(frameon=False, fontsize=8)
        axis.set_title(method)
        axis.set_xlabel(xlabel); axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25)
    fig.suptitle("Published-method native-anchor convergence (seed 0)")
    fig.tight_layout()
    pdf = FIGURE_DIR / "baseline_native_convergence.pdf"
    png = FIGURE_DIR / "baseline_native_convergence.png"
    fig.savefig(pdf, bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return [output_csv, pdf, png]


if __name__ == "__main__":
    print("\n".join(str(path) for path in build()))
