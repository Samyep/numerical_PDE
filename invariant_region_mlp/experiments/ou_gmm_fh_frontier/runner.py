"""Frozen 100D OU/GMM HJB: raw Full-History MLP time-accuracy sweep.

Mathematically necessary OU change relative to Brownian-only FH:
    dX = X ds + dW,  X_{t+h} = exp(h)*x + sqrt((exp(2h)-1)/2)*Z
    z = grad_x u; terminal derivative = exp(h) grad g(X_{t+h})
    BEL score = exp(h) / sqrt(var(h)) * Z.

No terminal truth, GMM responsibilities of the *solution*, or learned
Actor enters any MLP update. Only the known terminal g and grad(g) are used.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.benchmarks import OUGMMBenchmark
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.preflight import _array_hash
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.protocol import stable_seed

# Fix Picard depth; scan ONLY the Monte Carlo base M.
DEPTH = 4
M_GRID = (4, 8, 12, 16, 24, 32, 48, 64)
SEED = 0
D = 100
TEST_POINTS = 960
TARGET = "tangent_actor"
ALPHA = 0.5  # same beta time sampler as the existing FH-MLP control
MAX_TERMINAL_FLOATS = 1_600_000


@dataclass
class Counters:
    terminal_samples: int = 0
    transition_samples: int = 0
    driver_evaluations: int = 0
    recursive_states: int = 0
    recursive_calls: int = 0


def sha_file(path: Path) -> str:
    obj = sha256()
    with path.open("rb") as h:
        for block in iter(lambda: h.read(1 << 20), b""):
            obj.update(block)
    return obj.hexdigest()


def write_json_atomic(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(tmp, path)


def save_npz_atomic(path: Path, **arrays) -> None:
    tmp = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp, path)


def ou_variance(h: torch.Tensor) -> torch.Tensor:
    return 0.5 * torch.expm1(2.0 * h)


def ou_state(
    x: torch.Tensor, h: torch.Tensor, z: torch.Tensor
) -> torch.Tensor:
    """Broadcast x: (*,d), h: (*), z: (*,d) for OU dX=Xds+dW."""
    return torch.exp(h)[..., None] * x + torch.sqrt(ou_variance(h).clamp_min(0))[..., None] * z


def weighted_ou_bel(
    remaining: torch.Tensor, fraction: torch.Tensor, delta: torch.Tensor,
    alpha: float = ALPHA
) -> torch.Tensor:
    """Importance weight times OU BEL score, WITHOUT its standard normal Z.

    For beta(alpha,1) time sampling, importance=(T-t)*fraction^(1-alpha)/alpha.
    OU gradient score=e^delta Z/sqrt(var(delta)).
    The algebraically equivalent stable product is
    sqrt(T-t)/alpha * fraction^(1/2-alpha) * e^delta
        * sqrt(delta/var(delta)) * Z.
    In particular alpha=1/2 cancels the small-fraction singularity.
    """
    # delta / var(delta) -> 1 as delta -> 0; don't form 0/0.
    ratio = torch.where(
        delta.abs() < 1e-6, torch.ones_like(delta),
        torch.sqrt((delta / ou_variance(delta).clamp_min(1e-30)).clamp_min(0)),
    )
    return (
        torch.sqrt(remaining.clamp_min(0))[:, None] / alpha
        * fraction.pow(0.5 - alpha)
        * torch.exp(delta) * ratio
    )


class RawOUFullHistory:
    """Same telescoping FH recursion as GenericFullHistoryMLP, OU kernel.

    U_0=0, U_1=E[g(X_T)] and for n>=2, differences
    f(U_l)-1[l>1] f(U_{l-1}) at levels l=1,...,n-1.
    We integrate the known constant -d source exactly once for n>=2,
    without sampling its zero gradient. The quadratic part remains sampled.
    The base f(u,z) = -0.5||z||^2-d is evaluated without projection.
    """

    def __init__(
        self, benchmark: OUGMMBenchmark, M: int, device: torch.device,
        *, seed: int, dtype: torch.dtype = torch.float32,
    ):
        self.benchmark, self.M, self.device, self.dtype = benchmark, M, device, dtype
        self.d = benchmark.d
        self.T = benchmark.T
        self.generator = torch.Generator(device=device).manual_seed(int(seed))
        self.stats = Counters()
        self.mean = benchmark.mean.to(device=device, dtype=dtype)
        self.log_pi = benchmark.pi.to(device=device, dtype=dtype).log()
        self.mean2 = self.mean.square().sum(-1)
        self.log_normalizer = -0.5 * self.d * math.log(4.0 * math.pi)

    def normal(self, shape) -> torch.Tensor:
        return torch.randn(shape, generator=self.generator, dtype=self.dtype, device=self.device)

    def terminal_value_gradient(self, x: torch.Tensor):
        """Known terminal GMM g=-log p and exact analytic grad(g).

        Variance at terminal is 2*I, five components; stable logsumexp.
        This is terminal data, NOT use of the exact nonterminal HJB solution.
        """
        x2 = x.square().sum(dim=-1, keepdim=True)
        distance2 = (x2 + self.mean2 - 2.0 * (x @ self.mean.T)).clamp_min(0)
        logits = self.log_pi + self.log_normalizer - distance2 / 4.0
        log_density = torch.logsumexp(logits, dim=-1)
        weights = torch.softmax(logits, dim=-1)
        gradient = 0.5 * (x - weights @ self.mean)
        return -log_density, gradient

    def terminal_estimate(self, n: int, t: torch.Tensor, x: torch.Tensor):
        b = len(t)
        count = self.M ** n
        remaining = (self.T - t).clamp_min(0)
        exp_h = torch.exp(remaining)
        std_h = torch.sqrt(ou_variance(remaining))
        parent, _ = self.terminal_value_gradient(x)
        # Accumulate sums in float64 even when paths run in float32.
        value_sum = torch.zeros(b, dtype=torch.float64, device=self.device)
        grad_sum = torch.zeros((b, self.d), dtype=torch.float64, device=self.device)
        block = max(1, min(count, MAX_TERMINAL_FLOATS // max(b * self.d, 1)))
        for begin in range(0, count, block):
            size = min(block, count - begin)
            noise = self.normal((b, size, self.d))
            end = x[:, None, :] * exp_h[:, None, None] + std_h[:, None, None] * noise
            values, gradients = self.terminal_value_gradient(end.reshape(-1, self.d))
            values = values.reshape(b, size)
            gradients = gradients.reshape(b, size, self.d)
            value_sum += (values - parent[:, None]).double().sum(dim=1)
            grad_sum += gradients.double().sum(dim=1)
        result = torch.empty((b, self.d + 1), dtype=self.dtype, device=self.device)
        result[:, 0] = parent + (value_sum / count).to(self.dtype)
        result[:, 1:] = exp_h[:, None] * (grad_sum / count).to(self.dtype)
        self.stats.terminal_samples += b * count
        return result

    def driver(self, state: torch.Tensor):
        z = state[..., 1:]
        # Treat the known spatially constant -d source deterministically.
        # This removes pure MC noise in its gradient; no oracle value is used.
        return -0.5 * z.square().sum(dim=-1)

    def solve(self, n: int, t: torch.Tensor, x: torch.Tensor):
        t = t.to(dtype=self.dtype, device=self.device).reshape(-1)
        x = x.to(dtype=self.dtype, device=self.device).reshape(-1, self.d)
        b = len(t)
        self.stats.recursive_calls += 1
        self.stats.recursive_states += b
        if n == 0:
            return torch.zeros((b, self.d+1), dtype=self.dtype, device=self.device)

        output = self.terminal_estimate(n, t, x)
        remaining = (self.T - t).clamp_min(0)
        if n >= 2:
            # Source -d contributes exactly -d*(T-t) to value, zero to z.
            # n=1 is the heat/OU semigroup initialization (no f yet).
            output[:, 0] += -self.d * remaining
        tiny = torch.finfo(self.dtype).tiny

        for level in range(1, n):
            siblings = self.M ** (n-level)
            uniform = torch.rand(
                (b, siblings), generator=self.generator,
                dtype=self.dtype, device=self.device,
            ).clamp_min(tiny)
            fraction = uniform.pow(1.0 / ALPHA)
            delta = remaining[:, None] * fraction
            noise = self.normal((b, siblings, self.d))
            child_x = (
                torch.exp(delta)[..., None] * x[:, None, :]
                + torch.sqrt(ou_variance(delta))[..., None] * noise
            )
            child_t = (t[:, None] + delta).reshape(-1)
            child_x = child_x.reshape(-1, self.d)

            high = self.solve(level, child_t, child_x).reshape(b, siblings, self.d+1)
            difference = self.driver(high)
            self.stats.driver_evaluations += b * siblings
            if level > 1:
                low = self.solve(level-1, child_t, child_x).reshape(b, siblings, self.d+1)
                difference = difference - self.driver(low)
                self.stats.driver_evaluations += b * siblings

            importance = remaining[:, None] * fraction.pow(1.0 - ALPHA) / ALPHA
            output[:, 0] += (importance * difference).mean(dim=1)
            bel_times_importance = weighted_ou_bel(remaining, fraction, delta)
            output[:, 1:] += (
                (bel_times_importance * difference)[..., None] * noise
            ).mean(dim=1)
            self.stats.transition_samples += b * siblings
        return output


def load_comparison(source_root: Path):
    root = source_root / "invariant_region_mlp/results/tangent_actor_icml_final_suite"
    main_path = root / "MAIN_RESULTS.csv"
    parameters = root / "inputs/ou_gmm_d100.npz"
    predictions = root / "paired_predictions/main_ou_gmm_d100_tangent_actor_s0_n64p0_p1024_mcv.npz"
    for path in (main_path, parameters, predictions):
        if not path.exists():
            raise FileNotFoundError(f"Frozen original suite artifact is missing: {path}")
    with main_path.open(newline="") as fp:
        rows = [r for r in csv.DictReader(fp) if
                r["pde"] == "ou_gmm" and r["method"] == TARGET
                and r["dimension"] == "100" and r["seed"] == "0"
                and r["estimator"] == "mcv" and r["scope"] == "primary"]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one frozen OU Actor row, found {len(rows)}")
    row = rows[0]
    with np.load(predictions, allow_pickle=False) as data:
        t = np.array(data["t"], dtype=np.float64)
        x = np.array(data["x"], dtype=np.float64)
        truth = np.array(data["truth"], dtype=np.float64)
        actor_prediction = np.array(data["prediction"], dtype=np.float64)
    if t.shape != (TEST_POINTS,) or x.shape != (TEST_POINTS, D):
        raise RuntimeError("Unexpected OU/GMM whole-function test shape")
    if truth.shape != t.shape or actor_prediction.shape != t.shape:
        raise RuntimeError("Paired prediction/reference shape mismatch")
    if not all(np.isfinite(v).all() for v in (t, x, truth, actor_prediction)):
        raise RuntimeError("Nonfinite frozen benchmark artifacts")
    bank_hash = _array_hash(t, x)
    if bank_hash != row["test_bank_array_sha256"]:
        raise RuntimeError("Frozen OU/GMM test bank differs from MAIN_RESULTS.csv")
    actor_rmse = float(np.sqrt(np.mean((actor_prediction-truth)**2)))
    recorded_rmse = float(row["value_rmse"])
    if abs(actor_rmse - recorded_rmse) > 1e-8:
        raise RuntimeError(f"Frozen Actor RMSE mismatch: {actor_rmse} vs {recorded_rmse}")
    benchmark = OUGMMBenchmark.from_npz(parameters)
    with torch.no_grad():
        fresh_truth = benchmark.exact_value(
            torch.tensor(t, dtype=torch.float64),
            torch.tensor(x, dtype=torch.float64)
        ).numpy()
    if float(np.max(np.abs(truth-fresh_truth))) > 1e-7:
        raise RuntimeError("Stored truth is not the exact OU/GMM reference")
    return {
        "benchmark": benchmark, "t": t, "x": x, "truth": truth,
        "actor_rmse": actor_rmse,
        "actor_total_seconds": float(row["total_seconds"]),
        "actor_inference_seconds": float(row["inference_seconds"]),
        "bank_array_sha256": bank_hash,
        "truth_array_sha256": _array_hash(truth),
        "parameter_npz_sha256": sha_file(parameters),
        "canonical_actor_npz_sha256": sha_file(predictions),
        "actor_main_results_sha256": sha_file(main_path),
        "source_root": str(source_root),
    }


def static_metadata(comp: dict) -> dict:
    return {name: comp[name] for name in (
        "actor_rmse", "actor_total_seconds", "actor_inference_seconds",
        "bank_array_sha256", "truth_array_sha256",
        "parameter_npz_sha256", "canonical_actor_npz_sha256",
        "actor_main_results_sha256", "source_root",
    )}


def eval_chunk_size(M: int) -> int:
    return max(1, min(16, 2_000_000 // (D * M**DEPTH)))


def run_one(
    comp: dict, M: int, output_dir: Path, device: torch.device, *,
    max_seconds: float, seed: int = SEED,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"ou_gmm_d100_fh_n{DEPTH}_M{M}_seed{seed}"
    final_path = output_dir / (stem + ".json")
    progress_path = output_dir / (stem + ".progress.npz")
    pred_path = output_dir / (stem + ".npz")

    if final_path.exists():
        completed = json.loads(final_path.read_text())
        if completed.get("status") == "COMPLETE":
            print(f"SKIP completed M={M}: {final_path}", flush=True)
            return completed

    chunk = eval_chunk_size(M)
    prediction = np.full(TEST_POINTS, np.nan, dtype=np.float64)
    finished = 0
    accumulated = 0.0
    if progress_path.exists():
        with np.load(progress_path, allow_pickle=False) as old:
            if str(old["bank_hash"]) != comp["bank_array_sha256"]:
                raise RuntimeError("Resume bank hash mismatch")
            if int(old["chunk_size"]) != chunk:
                raise RuntimeError("Resume chunk size mismatch")
            prediction = np.array(old["prediction"], dtype=np.float64)
            finished = int(old["finished"])
            accumulated = float(old["elapsed"])
        if finished % chunk != 0 and finished != TEST_POINTS:
            raise RuntimeError("Progress not at an evaluation chunk boundary")

    state = {
        "status": "RUNNING", "method": "raw_full_history_mlp_ou_kernel",
        "depth": DEPTH, "M": M, "seed": seed, "dtype": "float32",
        "states": TEST_POINTS, "time_beta_alpha": ALPHA,
        "terminal_gradient": "analytic_known_g_pathwise",
        "known_constant_source": "-d*(T-t) integrated exactly for n>=2",
        "OU_kernel": "X_s=exp(h)x+sqrt((exp(2h)-1)/2)Z",
        "OU_BEL": "exp(h)/sqrt((exp(2h)-1)/2)*Z",
        "training_seconds": 0.0, "value_rmse": None,
        "completed_states": finished, "gpu_seconds": accumulated,
        "actor": static_metadata(comp),
        "single_parameter_sweep": True, "all_grid_M": list(M_GRID),
        "truth_used_for_training": False, "oracle_actor_used": False,
        "terminal_reference_used_inside_MLP": False,
    }
    write_json_atomic(final_path, state)

    try:
        for chunk_index, begin in enumerate(range(0, TEST_POINTS, chunk)):
            end = min(begin + chunk, TEST_POINTS)
            if end <= finished:
                continue
            if accumulated >= max_seconds:
                state["status"] = "OUT_OF_TIME"
                break
            model = RawOUFullHistory(
                comp["benchmark"], M, device,
                seed=stable_seed("ou_gmm_full_history_fixed_depth", seed, DEPTH, M, chunk_index),
            )
            tt = torch.from_numpy(comp["t"][begin:end]).to(device=device, dtype=torch.float32)
            xx = torch.from_numpy(comp["x"][begin:end]).to(device=device, dtype=torch.float32)
            if device.type == "cuda":
                torch.cuda.synchronize(device)
            started = time.perf_counter()
            with torch.no_grad():
                result = model.solve(DEPTH, tt, xx)
            if device.type == "cuda":
                torch.cuda.synchronize(device)
            seconds = time.perf_counter() - started
            value = result[:, 0].double().cpu().numpy()
            accumulated += seconds
            if not np.all(np.isfinite(value)):
                state["status"] = "NONFINITE"
                break
            prediction[begin:end] = value
            finished = end
            save_npz_atomic(
                progress_path, finished=np.asarray(finished), chunk_size=np.asarray(chunk),
                elapsed=np.asarray(accumulated), prediction=prediction,
                bank_hash=np.asarray(comp["bank_array_sha256"]),
            )
            state["completed_states"] = finished
            state["gpu_seconds"] = accumulated
            write_json_atomic(final_path, state)
            if chunk_index == 0 or end == TEST_POINTS or chunk_index % max(1, (32//chunk)) == 0:
                print(json.dumps({
                    "M": M, "done": end, "all": TEST_POINTS,
                    "gpu_s": round(accumulated, 3), "status": state["status"]
                }), flush=True)
        if finished == TEST_POINTS and state["status"] == "RUNNING":
            e = prediction - comp["truth"]
            state["value_rmse"] = float(np.sqrt(np.mean(e**2)))
            state["relative_l2"] = float(
                np.linalg.norm(e)/np.linalg.norm(comp["truth"]))
            state["bias"] = float(e.mean())
            state["mae"] = float(np.abs(e).mean())
            state["status"] = "COMPLETE"
            save_npz_atomic(
                pred_path, t=comp["t"], x=comp["x"],
                truth=comp["truth"], prediction=prediction
            )
            state["prediction_file"] = str(pred_path)
            state["prediction_sha256"] = sha_file(pred_path)
            state["hit_actor_threshold"] = state["value_rmse"] <= comp["actor_rmse"]
    except (KeyboardInterrupt, Exception) as err:
        state["status"] = "INTERRUPTED" if isinstance(err, KeyboardInterrupt) else "FAILED"
        state["error"] = repr(err)
        write_json_atomic(final_path, state)
        raise
    write_json_atomic(final_path, state)
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True,
                        help="Original HPC source tree with immutable full-suite input NPZs")
    parser.add_argument("--out-dir", type=Path, default=Path("invariant_region_mlp/results/ou_gmm_fh_frontier"))
    parser.add_argument("--M", type=int, default=None, help="One of the prespecified M values")
    parser.add_argument("--all", action="store_true", help="Evaluate M grid in ascending order")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--max-cell-seconds", type=float, default=3600)
    parser.add_argument("--max-sweep-seconds", type=float, default=7200)
    args = parser.parse_args()
    if (args.M is None) == (not args.all):
        parser.error("specify exactly one of --M or --all")
    if args.M is not None and args.M not in M_GRID:
        parser.error(f"Only fixed M values permitted: {M_GRID}")
    comp = load_comparison(args.source_root.resolve())
    summary = {"target_rmse":comp["actor_rmse"], "target_seconds":comp["actor_total_seconds"],
               "depth":DEPTH, "M_grid":list(M_GRID), "states":TEST_POINTS,
               "frozen_bank_hash":comp["bank_array_sha256"], "status":"INPUTS_AUDITED"}
    print(json.dumps(summary, indent=2), flush=True)
    if args.dry_run:
        return
    device = torch.device(args.device)
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("Run the timed 960-point benchmark on GH200 CUDA, not CPU")
    scan_start = time.perf_counter()
    for M in M_GRID if args.all else (args.M,):
        used = time.perf_counter() - scan_start
        remaining = min(args.max_cell_seconds, args.max_sweep_seconds - used)
        if remaining <= 0:
            print("SWEEP_OUT_OF_TIME", flush=True)
            break
        item = run_one(comp, M, args.out_dir, device, max_seconds=remaining)
        if item["status"] != "COMPLETE":
            print(f"STOP at M={M}, status={item['status']}; no fabricated RMSE", flush=True)
            break
        # The first threshold hit is all that the user asks to locate.
        if args.all and item["hit_actor_threshold"]:
            print(f"FIRST_OBSERVED_HIT M={M}: {item['gpu_seconds']:.2f} sec", flush=True)
            break


if __name__ == "__main__":
    main()
