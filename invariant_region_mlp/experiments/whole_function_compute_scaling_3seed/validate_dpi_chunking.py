"""Direct numerical validation of the estimator-preserving DPI chunker."""

from __future__ import annotations

import json
from types import MethodType, SimpleNamespace

import torch

from invariant_region_mlp.experiments.whole_function_compute_frontier import (
    baseline_frontier as base,
)
from invariant_region_mlp.experiments.whole_function_compute_frontier.common import (
    write_json_atomic,
)

from .engine import _chunked_dpi_integral_with_gradients
from .protocol import RESULTS_DIR


def main() -> None:
    base._external_root()
    from picard.data import OnlineDataGenerator
    from picard.equations import Cha

    torch.set_default_dtype(torch.float64)
    torch.manual_seed(918273)
    equation = Cha(nx=2, alpha=1.0, T=1.0)
    solution = torch.nn.Sequential(
        torch.nn.Linear(3, 8), torch.nn.Tanh(), torch.nn.Linear(8, 1)
    ).double()
    generator = OnlineDataGenerator(
        equation, solution, N=2, i=1, device="cpu",
        n_estimate_terminal=4, n_estimate_integral=16,
        hessian_approximation=SimpleNamespace(method=None, kwargs={}),
    )
    tx = torch.tensor([[0.15, -0.2, 0.3], [0.55, 0.1, -0.4]])
    original = OnlineDataGenerator.estimate_integral_with_gradients

    # When the nominal count is at most 8,192, the wrapper deliberately uses
    # one chunk.  It must therefore reproduce the upstream implementation down
    # to floating-point roundoff under an identical RNG state.
    torch.manual_seed(112233)
    reference_small = original(generator, tx.clone())
    torch.manual_seed(112233)
    chunked_small = _chunked_dpi_integral_with_gradients(generator, tx.clone())
    small_max_abs = float((reference_small - chunked_small).abs().max())
    if small_max_abs > 1e-12:
        raise AssertionError(f"single-chunk mismatch: {small_max_abs}")

    # For an actual multi-chunk call, replace only the Brownian draw routine by
    # a deterministic indexed tape.  The upstream one-shot sum and chunked sum
    # then consume exactly the same 2,500 contributions, so their outputs must
    # agree up to reduction-order roundoff.
    generator.n_estimate_integral = 2500
    calls: list[int] = []

    def deterministic_tape(self, points, count, return_dW=False):
        start = int(getattr(self, "_validation_offset", 0))
        self._validation_offset = start + count
        calls.append(count)
        batch = points.size(0)
        t0 = points[:, :1]
        x0 = points[:, 1:]
        fraction = (
            torch.arange(start, start + count, dtype=points.dtype)
            .add(0.5).div(self.n_estimate_integral)
        )
        fraction = fraction.repeat(batch).view(-1, 1)
        t = torch.repeat_interleave(t0, count, dim=0)
        x = torch.repeat_interleave(x0, count, dim=0)
        s = t + (self.T - t) * fraction
        index = torch.arange(start, start + count, dtype=points.dtype)
        dw_one = torch.stack((torch.sin(index + 0.1), torch.cos(index + 0.2)), dim=1)
        dw = dw_one.repeat(batch, 1)
        xs = x + torch.sqrt(s - t) * self.equation.alpha_sqrt * dw
        sxs = torch.cat([s, xs], dim=-1)
        if return_dW:
            return t, s, xs, sxs, dw
        return t, s, xs, sxs

    generator.generate_sx_for_integral = MethodType(deterministic_tape, generator)
    generator._validation_offset = 0
    calls.clear()
    reference_large = original(generator, tx.clone())
    reference_calls = list(calls)
    generator._validation_offset = 0
    calls.clear()
    chunked_large = _chunked_dpi_integral_with_gradients(generator, tx.clone())
    chunked_calls = list(calls)
    large_max_abs = float((reference_large - chunked_large).abs().max())
    if sum(chunked_calls) != generator.n_estimate_integral:
        raise AssertionError(f"sample-count drift: {chunked_calls}")
    if large_max_abs > 2e-11:
        raise AssertionError(f"multi-chunk mismatch: {large_max_abs}")

    payload = {
        "status": "PASS",
        "single_chunk_max_abs": small_max_abs,
        "multi_chunk_max_abs": large_max_abs,
        "upstream_calls": reference_calls,
        "chunked_calls": chunked_calls,
        "nominal_samples": int(generator.n_estimate_integral),
        "all_nominal_samples_used": sum(chunked_calls) == generator.n_estimate_integral,
        "finite": bool(torch.isfinite(chunked_large).all()),
    }
    write_json_atomic(RESULTS_DIR / "dpi_chunking_validation.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
