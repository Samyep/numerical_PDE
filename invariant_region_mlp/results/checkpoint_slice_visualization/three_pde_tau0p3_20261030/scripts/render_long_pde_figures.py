#!/usr/bin/env python3
"""Render one publication-style long comparison figure per completed PDE.

The figure layout mirrors the requested visual reference:

* a large reference surface and reference value/gradient slice on the left;
* four method columns containing surface, value/gradient, value-error, and
  gradient-error panels;
* one reference-derived solution scale and two method-shared logarithmic
  discrepancy scales.

Only frozen arrays already written by ``evaluate_slice.py`` are read.  The
renderer performs no smoothing, clipping, interpolation, or model selection.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import LogFormatterSciNotation, MaxNLocator
import numpy as np


# The script lives in ``<bundle>/scripts`` both in the source workspace and in
# the portable export, so the parent directory is the data root on any machine.
ROOT = Path(__file__).resolve().parents[1]
PDES = ("joint_concave", "nonconvex_coupled", "ou_gmm")
PDE_TITLES = {
    "joint_concave": "400D Jointly Concave",
    "nonconvex_coupled": "100D Nonconvex–Nonconcave",
    "ou_gmm": "100D OU/GMM HJB",
}
METHODS = ("tangent_actor", "dpi", "sdgd_pinn", "pinn_hte")
METHOD_LABELS = {
    "tangent_actor": "Tangent-Actor",
    "dpi": "DPI",
    "sdgd_pinn": "SDGD-PINN",
    "pinn_hte": "PINN-HTE",
}
TIME_FRACTION = 0.3
TIME_TAG = "tau0p3"


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 7.2,
            "axes.titlesize": 7.3,
            "axes.labelsize": 6.5,
            "xtick.labelsize": 5.7,
            "ytick.labelsize": 5.7,
            "axes.linewidth": 0.58,
            "xtick.major.width": 0.48,
            "ytick.major.width": 0.48,
            "xtick.major.size": 2.2,
            "ytick.major.size": 2.2,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.035,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def load_arrays(pde: str) -> tuple[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]:
    reference_path = ROOT / "references" / f"{pde}_{TIME_TAG}.npz"
    if not reference_path.exists():
        raise FileNotFoundError(reference_path)
    with np.load(reference_path) as archive:
        reference = {key: np.asarray(archive[key]) for key in archive.files}

    predictions: dict[str, dict[str, np.ndarray]] = {}
    for method in METHODS:
        path = ROOT / "predictions" / f"{pde}_{TIME_TAG}_{method}.npz"
        if not path.exists():
            raise FileNotFoundError(path)
        with np.load(path) as archive:
            predictions[method] = {
                key: np.asarray(archive[key]) for key in archive.files
            }
    return reference, predictions


def validate_arrays(
    pde: str,
    reference: dict[str, np.ndarray],
    predictions: dict[str, dict[str, np.ndarray]],
) -> dict[str, object]:
    expected_reference = {
        "axis": (81,),
        "s1": (81, 81),
        "s2": (81, 81),
        "truth": ((6561,), (81, 81)),
        "truth_gradient": ((6561, 2), (81, 81, 2)),
    }
    for key, expected in expected_reference.items():
        if key not in reference:
            raise KeyError(f"{pde}: missing reference key {key}")
        allowed = expected if isinstance(expected[0], tuple) else (expected,)
        if reference[key].shape not in allowed:
            raise ValueError(
                f"{pde}: {key} has shape {reference[key].shape}, expected {allowed}"
            )
        if not np.isfinite(reference[key]).all():
            raise ValueError(f"{pde}: non-finite reference values in {key}")

    method_shapes = {
        "prediction": (81, 81),
        "gradient": (81, 81, 2),
        "value_error": (81, 81),
        "gradient_error": (81, 81, 2),
        "value_discrepancy": (81, 81),
        "gradient_discrepancy": (81, 81),
    }
    audit: dict[str, object] = {"pde": pde, "time_fraction": TIME_FRACTION, "methods": {}}
    for method, arrays in predictions.items():
        for key, shape in method_shapes.items():
            if key not in arrays:
                raise KeyError(f"{pde}/{method}: missing {key}")
            if arrays[key].shape != shape:
                raise ValueError(
                    f"{pde}/{method}: {key} has shape {arrays[key].shape}, expected {shape}"
                )
            if not np.isfinite(arrays[key]).all():
                raise ValueError(f"{pde}/{method}: non-finite values in {key}")

        value_rmse = float(np.sqrt(np.mean(np.square(arrays["value_error"]))))
        gradient_rmse = float(
            np.sqrt(np.mean(np.sum(np.square(arrays["gradient_error"]), axis=-1)))
        )
        audit["methods"][method] = {
            "value_rmse": value_rmse,
            "reduced_gradient_vector_rmse": gradient_rmse,
        }
    return audit


def reshape_reference(reference: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    return (
        reference["truth"].reshape(81, 81),
        reference["truth_gradient"].reshape(81, 81, 2),
    )


def positive_log_norm(values: list[np.ndarray]) -> LogNorm:
    flattened = np.concatenate([np.asarray(value).ravel() for value in values])
    positive = flattened[np.isfinite(flattened) & (flattened > 0)]
    if positive.size == 0:
        return LogNorm(vmin=1e-12, vmax=1e-11)
    vmax = float(np.max(positive))
    vmin = max(float(np.quantile(positive, 0.01)), vmax * 1e-5, 1e-12)
    return LogNorm(vmin=vmin, vmax=max(vmax, vmin * 1.001))


def quiver_scale(
    truth_gradient: np.ndarray,
    predictions: dict[str, dict[str, np.ndarray]],
    axis: np.ndarray,
) -> float:
    gradients = [truth_gradient] + [predictions[m]["gradient"] for m in METHODS]
    q95 = max(
        float(np.quantile(np.linalg.norm(gradient, axis=-1), 0.95))
        for gradient in gradients
    )
    target_length = 0.105 * float(axis[-1] - axis[0])
    return max(q95 / max(target_length, 1e-12), 1e-12)


def finish_2d_axis(ax: mpl.axes.Axes, axis: np.ndarray) -> None:
    ticks = [float(axis[0]), 0.0, float(axis[-1])]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlabel(r"$s_1$", labelpad=0.3)
    ax.set_ylabel(r"$s_2$", labelpad=0.3)
    ax.tick_params(pad=1.0)
    for spine in ax.spines.values():
        spine.set_color("#27313c")
        spine.set_linewidth(0.62)


def draw_quiver(
    ax: mpl.axes.Axes,
    axis: np.ndarray,
    gradient: np.ndarray,
    scale: float,
) -> None:
    indices = np.arange(0, 81, 10)
    xx, yy = np.meshgrid(axis[indices], axis[indices], indexing="xy")
    local = gradient[np.ix_(indices, indices)]
    ax.quiver(
        xx,
        yy,
        local[..., 0],
        local[..., 1],
        color="#17212b",
        alpha=0.74,
        angles="xy",
        scale_units="xy",
        scale=scale,
        width=0.00410,
        headwidth=3.0,
        headlength=3.8,
        headaxislength=3.35,
        pivot="mid",
        zorder=5,
    )


def draw_heatmap(
    ax: mpl.axes.Axes,
    axis: np.ndarray,
    values: np.ndarray,
    norm: mpl.colors.Normalize,
    cmap: str,
) -> mpl.image.AxesImage:
    return ax.imshow(
        values,
        origin="lower",
        extent=(axis[0], axis[-1], axis[0], axis[-1]),
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        aspect="equal",
        rasterized=True,
    )


def draw_surface(
    ax: mpl.axes.Axes,
    s1: np.ndarray,
    s2: np.ndarray,
    values: np.ndarray,
    norm: mpl.colors.Normalize,
    large: bool,
) -> None:
    stride = 2
    ax.plot_surface(
        s1[::stride, ::stride],
        s2[::stride, ::stride],
        values[::stride, ::stride],
        cmap="RdBu_r",
        norm=norm,
        linewidth=0,
        antialiased=True,
        alpha=0.94,
        rasterized=True,
    )
    ax.view_init(elev=27, azim=-126)
    ax.set_xlabel(r"$s_1$", labelpad=-6.5 if large else -8.8)
    ax.set_ylabel(r"$s_2$", labelpad=-6.5 if large else -8.8)
    ax.set_zlabel(r"$u$", labelpad=-6.0 if large else -8.0)
    s1_min, s1_max = float(np.min(s1)), float(np.max(s1))
    s2_min, s2_max = float(np.min(s2)), float(np.max(s2))
    ax.set_xticks([s1_min, 0.0, s1_max])
    ax.set_yticks([s2_min, 0.0, s2_max])
    ax.zaxis.set_major_locator(MaxNLocator(nbins=3))
    ax.tick_params(pad=-1.0 if large else -2.4, labelsize=5.0 if large else 4.4)
    ax.grid(True, linewidth=0.35, color="#9da7b2", alpha=0.72)
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.set_facecolor((0.97, 0.975, 0.98, 1.0))
        pane.set_edgecolor((0.48, 0.52, 0.57, 0.76))
        pane.set_alpha(0.20)
    try:
        ax.set_box_aspect((1.25, 1.0, 0.62 if large else 0.54))
    except AttributeError:
        pass


def discrepancy_label(kind: str, rmse: float) -> str:
    return f"{kind}\nRMSE {rmse:.2e}"


def draw_discrepancy(
    ax: mpl.axes.Axes,
    axis: np.ndarray,
    values: np.ndarray,
    norm: LogNorm,
    label: str,
) -> mpl.image.AxesImage:
    image = draw_heatmap(ax, axis, np.maximum(values, norm.vmin), norm, "magma")
    ax.text(
        0.025,
        0.965,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.25,
        linespacing=0.95,
        color="white",
        bbox={
            "boxstyle": "round,pad=0.18,rounding_size=0.05",
            "facecolor": (0.08, 0.09, 0.11, 0.70),
            "edgecolor": (1.0, 1.0, 1.0, 0.18),
            "linewidth": 0.35,
        },
        zorder=7,
    )
    finish_2d_axis(ax, axis)
    return image


def add_colorbar(
    fig: mpl.figure.Figure,
    cax: mpl.axes.Axes,
    mappable: mpl.cm.ScalarMappable,
    label: str,
    logarithmic: bool,
) -> None:
    formatter = LogFormatterSciNotation(base=10) if logarithmic else None
    colorbar = fig.colorbar(mappable, cax=cax, orientation="vertical", format=formatter)
    colorbar.ax.tick_params(labelsize=5.0, width=0.42, length=1.9, pad=1.2)
    colorbar.outline.set_linewidth(0.52)
    colorbar.set_label(label, fontsize=5.8, labelpad=2.8)


def render_pde(pde: str) -> dict[str, object]:
    reference, predictions = load_arrays(pde)
    audit = validate_arrays(pde, reference, predictions)
    axis = reference["axis"]
    truth, truth_gradient = reshape_reference(reference)

    value_min = float(np.min(truth))
    value_max = float(np.max(truth))
    if np.isclose(value_min, value_max):
        delta = max(abs(value_min), 1.0) * 1e-3
        value_min, value_max = value_min - delta, value_max + delta
    solution_norm = Normalize(vmin=value_min, vmax=value_max, clip=False)
    value_error_norm = positive_log_norm(
        [predictions[method]["value_discrepancy"] for method in METHODS]
    )
    gradient_error_norm = positive_log_norm(
        [predictions[method]["gradient_discrepancy"] for method in METHODS]
    )
    arrow_scale = quiver_scale(truth_gradient, predictions, axis)

    fig = plt.figure(figsize=(13.25, 7.35), facecolor="#fbfbfa")
    grid = fig.add_gridspec(
        4,
        6,
        left=0.035,
        right=0.972,
        bottom=0.066,
        top=0.915,
        width_ratios=(1.62, 1.0, 1.0, 1.0, 1.0, 0.045),
        height_ratios=(1.0, 1.0, 1.0, 1.0),
        wspace=0.25,
        hspace=0.24,
    )
    fig.suptitle(
        rf"{PDE_TITLES[pde]}  $\;\mid\;$  $t=0.3T$",
        x=0.50,
        y=0.970,
        fontsize=10.0,
        fontweight="semibold",
    )

    reference_surface = fig.add_subplot(grid[0:2, 0], projection="3d")
    draw_surface(
        reference_surface,
        reference["s1"],
        reference["s2"],
        truth,
        solution_norm,
        large=True,
    )
    reference_surface.set_title("Reference surface", pad=2.8, fontweight="semibold")

    reference_heatmap = fig.add_subplot(grid[2:4, 0])
    draw_heatmap(reference_heatmap, axis, truth, solution_norm, "RdBu_r")
    draw_quiver(reference_heatmap, axis, truth_gradient, arrow_scale)
    reference_heatmap.set_title(
        "Reference 2D + gradient", pad=3.2, fontweight="semibold"
    )
    finish_2d_axis(reference_heatmap, axis)

    for column, method in enumerate(METHODS, start=1):
        arrays = predictions[method]
        surface = fig.add_subplot(grid[0, column], projection="3d")
        draw_surface(
            surface,
            reference["s1"],
            reference["s2"],
            arrays["prediction"],
            solution_norm,
            large=False,
        )
        surface.set_title(METHOD_LABELS[method], pad=2.0, fontweight="semibold")

        prediction = fig.add_subplot(grid[1, column])
        draw_heatmap(
            prediction, axis, arrays["prediction"], solution_norm, "RdBu_r"
        )
        draw_quiver(prediction, axis, arrays["gradient"], arrow_scale)
        prediction.set_title("2D + gradient", pad=2.8, fontsize=6.3, fontweight="semibold")
        finish_2d_axis(prediction, axis)

        method_audit = audit["methods"][method]
        value_error = fig.add_subplot(grid[2, column])
        draw_discrepancy(
            value_error,
            axis,
            arrays["value_discrepancy"],
            value_error_norm,
            discrepancy_label("Value discrepancy", method_audit["value_rmse"]),
        )

        gradient_error = fig.add_subplot(grid[3, column])
        draw_discrepancy(
            gradient_error,
            axis,
            arrays["gradient_discrepancy"],
            gradient_error_norm,
            discrepancy_label(
                "Gradient discrepancy",
                method_audit["reduced_gradient_vector_rmse"],
            ),
        )

    solution_cax = fig.add_subplot(grid[0:2, 5])
    value_error_cax = fig.add_subplot(grid[2, 5])
    gradient_error_cax = fig.add_subplot(grid[3, 5])
    add_colorbar(
        fig,
        solution_cax,
        mpl.cm.ScalarMappable(norm=solution_norm, cmap="RdBu_r"),
        r"$u(t,x)$",
        logarithmic=False,
    )
    add_colorbar(
        fig,
        value_error_cax,
        mpl.cm.ScalarMappable(norm=value_error_norm, cmap="magma"),
        r"$|\widehat{u}-u|$",
        logarithmic=True,
    )
    add_colorbar(
        fig,
        gradient_error_cax,
        mpl.cm.ScalarMappable(norm=gradient_error_norm, cmap="magma"),
        r"$\|\widehat{\nabla_s u}-\nabla_s u\|_2$",
        logarithmic=True,
    )

    stem = ROOT / "figures" / f"{pde}_{TIME_TAG}_long_comparison"
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), dpi=300, facecolor=fig.get_facecolor())
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor=fig.get_facecolor())
    fig.savefig(stem.with_suffix(".svg"), facecolor=fig.get_facecolor())
    plt.close(fig)
    audit.update(
        {
            "output_pdf": str(stem.with_suffix(".pdf")),
            "output_png": str(stem.with_suffix(".png")),
            "output_svg": str(stem.with_suffix(".svg")),
            "solution_limits": [value_min, value_max],
            "value_discrepancy_limits": [
                value_error_norm.vmin,
                value_error_norm.vmax,
            ],
            "gradient_discrepancy_limits": [
                gradient_error_norm.vmin,
                gradient_error_norm.vmax,
            ],
        }
    )
    return audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pde",
        action="append",
        choices=PDES,
        help="PDE to render; repeat as needed. Defaults to all three.",
    )
    args = parser.parse_args()
    configure_style()
    selected = tuple(args.pde) if args.pde else PDES
    audits = [render_pde(pde) for pde in selected]
    audit_path = ROOT / "figures" / f"long_comparison_{TIME_TAG}_audit.json"
    audit_path.write_text(json.dumps(audits, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"rendered": list(selected), "audit": str(audit_path)}, sort_keys=True))


if __name__ == "__main__":
    main()
