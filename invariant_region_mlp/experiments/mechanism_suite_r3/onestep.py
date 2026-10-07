"""Registered C2 terminal-block check of the one-step Jensen prediction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .equations import BASE_SEED, RESULTS_ROOT, make_equation
from .run import write_json_atomic


OUTPUT_PATH = RESULTS_ROOT / "s2_onestep_check.json"
DIMENSIONS = (10, 50, 200, 1000)
SAMPLE_SIZES = (4, 16)
CONFIGS = ("C2-convex", "C2-cancel", "C2-flip")
SQRT2 = math.sqrt(2.0)
NODES, WEIGHTS = np.polynomial.hermite_e.hermegauss(120)
WEIGHTS = WEIGHTS / np.sum(WEIGHTS)


def terminal_phi(s: np.ndarray) -> np.ndarray:
    lambdas = np.asarray([0.5, 1.5, 3.0], dtype=np.float64)
    logits = np.asarray(s, dtype=np.float64)[..., None] * lambdas
    maximum = np.max(logits, axis=-1, keepdims=True)
    return -(
        maximum[..., 0]
        + np.log(np.sum(np.exp(logits - maximum), axis=-1))
        - math.log(3.0)
    )


def run_case(
    pde_id: str,
    d: int,
    M: int,
    *,
    s0: float = 0.3,
    h: float = 0.2,
) -> dict[str, Any]:
    equation = make_equation(pde_id, d)
    repetitions = 20_000 if d <= 200 else 6_000
    chunk = 1000 if d <= 200 else 250
    rng = np.random.default_rng(
        np.random.SeedSequence(
            [BASE_SEED, d, M, {"C2-convex": 1, "C2-cancel": 2, "C2-flip": 3}[pde_id]]
        )
    )
    w = equation.w
    x = s0 * w
    delta_q = terminal_phi(s0 + SQRT2 * math.sqrt(h) * NODES) - terminal_phi(s0)
    c_h = float(np.sum(WEIGHTS * delta_q**2) / h)
    predicted_factor = float(equation.predicted_orthogonal_bias_factor())
    predicted_gap = c_h * predicted_factor / M

    sum_gap = 0.0
    sum_gap2 = 0.0
    sum_raw = 0.0
    count = 0
    for begin in range(0, repetitions, chunk):
        size = min(chunk, repetitions - begin)
        normal = rng.standard_normal((size, M, d), dtype=np.float64)
        shifted = (x + SQRT2 * math.sqrt(h) * normal) @ w
        delta = terminal_phi(shifted) - terminal_phi(s0)
        zhat = np.einsum("rm,rmd->rd", delta, normal) / (M * math.sqrt(h))
        parallel_scalar = zhat @ w
        parallel = parallel_scalar[:, None] * w
        raw_f = equation.generator(np.zeros(size), zhat)
        parallel_f = equation.generator(np.zeros(size), parallel)
        gap = raw_f - parallel_f
        sum_gap += float(np.sum(gap))
        sum_gap2 += float(np.sum(gap**2))
        sum_raw += float(np.sum(raw_f))
        count += size
    measured = sum_gap / count
    variance = max(sum_gap2 / count - measured**2, 0.0)
    standard_error = math.sqrt(variance / count)
    ratio = measured / predicted_gap if predicted_gap != 0.0 else None
    return {
        "pde": pde_id,
        "d": d,
        "M": M,
        "repetitions": repetitions,
        "chunk_size": chunk,
        "s": s0,
        "h": h,
        "c_h": c_h,
        "K": predicted_factor,
        "predicted_c_over_M_times_K": predicted_gap,
        "measured_orthogonal_jensen_gap": measured,
        "measured_standard_error": standard_error,
        "measured_over_predicted": ratio,
        "mean_raw_generator": sum_raw / count,
        "criterion_applicable": abs(predicted_factor) >= 5.0,
        "criterion_passed": (
            bool(0.9 <= ratio <= 1.1)
            if abs(predicted_factor) >= 5.0 and ratio is not None
            else None
        ),
    }


def run_all(path: Path = OUTPUT_PATH) -> dict[str, Any]:
    rows = []
    for pde_id in CONFIGS:
        for d in DIMENSIONS:
            for M in SAMPLE_SIZES:
                row = run_case(pde_id, d, M)
                rows.append(row)
                print(
                    f"one-step {pde_id} d={d} M={M} "
                    f"ratio={row['measured_over_predicted']:.4f}",
                    flush=True,
                )
    applicable = [row for row in rows if row["criterion_applicable"]]
    payload = {
        "schema_version": 1,
        "base_seed": BASE_SEED,
        "estimator": "brute-force R^d centred terminal-block gradient estimator",
        "measured_gap": "E[f(z_hat)-f((w.z_hat)w)], isolating the orthogonal Jensen term",
        "prediction": "(c_h/M) K",
        "rows": rows,
        "T3": {
            "threshold": "ratio in [0.9,1.1] for every row with |K|>=5",
            "applicable_rows": len(applicable),
            "passed": bool(applicable and all(row["criterion_passed"] for row in applicable)),
        },
    }
    write_json_atomic(path, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run_all(args.output)


if __name__ == "__main__":
    main()

