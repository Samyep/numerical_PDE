"""Exact, checkpointable execution of the official monolithic SDGD script.

The upstream HJB-Log file exposes no resume hook and a 10,000-epoch run is
longer than this cluster's two-hour interactive limit.  This launcher executes
the official source verbatim through its class definitions and replaces only
the final uncheckpointed ``model.train_adam()`` call with the same loop split
at epoch boundaries.  Model, optimizer, scheduler, collocation state, and all
Torch/NumPy RNG states are persisted, so a split run follows the same update
sequence as an uninterrupted run.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np
import torch

from .common import OFFICIAL_WORK_DIR, RAW_DIR, RESULTS_DIR, atomic_json, ensure_dirs


SOURCE = Path(
    "/work/nvme/bgvf/ywang151/original_pde_source_audit/SDGD_PINN/Section_5.3/HJB_Log.py"
)
SOURCE_COMMIT = "c597cc4fe9caac0b3ca118ea84453405f25ff833"


def _load_official_definitions(seed: int, total_epochs: int) -> dict[str, Any]:
    text = SOURCE.read_text()
    marker = "\nmodel = PINN()\n"
    if text.count(marker) != 1:
        raise RuntimeError("official SDGD source layout changed")
    definitions = text.split(marker, 1)[0]
    namespace: dict[str, Any] = {
        "__file__": str(SOURCE),
        "__name__": "__main__",
        "__package__": None,
    }
    previous_argv = sys.argv
    sys.argv = [
        str(SOURCE),
        "--SEED", str(seed),
        "--epochs", str(total_epochs),
        "--device", "cuda",
    ]
    try:
        exec(compile(definitions, str(SOURCE), "exec"), namespace)
    finally:
        sys.argv = previous_argv
    return namespace


def _rng_state() -> dict[str, Any]:
    return {
        "numpy": np.random.get_state(),
        "torch_cpu": torch.get_rng_state(),
        "torch_cuda": torch.cuda.get_rng_state_all(),
    }


def _restore_rng(state: dict[str, Any]) -> None:
    np.random.set_state(state["numpy"])
    # ``map_location='cuda'`` is used for the continuation payload so model
    # and optimizer tensors resume on the allocated GPU.  RNG APIs, however,
    # require their serialized ByteTensor states on CPU, including CUDA
    # generator states.  Moving only these opaque state tensors back to CPU
    # changes no random bits.
    torch.set_rng_state(state["torch_cpu"].detach().cpu())
    torch.cuda.set_rng_state_all(
        [value.detach().cpu() for value in state["torch_cuda"]]
    )


def _move_optimizer_state(optimizer: torch.optim.Optimizer, device: torch.device) -> None:
    for state in optimizer.state.values():
        for key, value in state.items():
            if torch.is_tensor(value):
                state[key] = value.to(device)


def run_chunk(
    *,
    seed: int,
    start_epoch: int,
    end_epoch: int,
    total_epochs: int,
    tag: str,
) -> dict[str, Any]:
    ensure_dirs()
    if not 0 <= start_epoch < end_epoch <= total_epochs:
        raise ValueError("require 0 <= start < end <= total_epochs")
    if not torch.cuda.is_available():
        raise RuntimeError("a CUDA allocation is required")
    work = OFFICIAL_WORK_DIR / "sdgd"
    work.mkdir(parents=True, exist_ok=True)
    (work / "saved_loss_l2").mkdir(exist_ok=True)
    continuation = work / f"sdgd_hjb_log_s{seed}_{tag}_continuation.pt"
    invocation_started = time.perf_counter()
    previous = Path.cwd()
    os.chdir(work)
    try:
        namespace = _load_official_definitions(seed, total_epochs)
        args = namespace["args"]
        model = namespace["PINN"]()
        optimizer = torch.optim.Adam(model.net_params_pinn, lr=model.adam_lr)
        lr_lambda = lambda epoch: 1 - epoch / total_epochs
        scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)
        previous_loop_seconds = 0.0
        previous_process_seconds = 0.0
        previous_peak_gpu_bytes = 0
        current_l2: float | None = None
        if start_epoch:
            if not continuation.exists():
                raise FileNotFoundError(continuation)
            payload = torch.load(continuation, map_location="cuda", weights_only=False)
            if int(payload["next_epoch"]) != start_epoch:
                raise RuntimeError(
                    f"checkpoint expects epoch {payload['next_epoch']}, requested {start_epoch}"
                )
            if int(payload["total_epochs"]) != total_epochs or int(payload["seed"]) != seed:
                raise RuntimeError("continuation metadata differs")
            model.u_net.load_state_dict(payload["model_state_dict"])
            optimizer.load_state_dict(payload["optimizer_state_dict"])
            _move_optimizer_state(optimizer, torch.device("cuda"))
            scheduler.load_state_dict(payload["scheduler_state_dict"])
            model.xf = payload["xf"].to("cuda").requires_grad_()
            model.saved_loss = list(payload["saved_loss"])
            model.saved_l2 = list(payload["saved_l2"])
            current_l2 = payload["current_l2"]
            previous_loop_seconds = float(payload["cumulative_loop_seconds"])
            previous_process_seconds = float(payload["cumulative_process_seconds"])
            previous_peak_gpu_bytes = int(payload["peak_gpu_bytes"])
            _restore_rng(payload["rng_state"])
        elif continuation.exists():
            continuation.unlink()

        torch.cuda.reset_peak_memory_stats()
        torch.cuda.synchronize()
        loop_started = time.perf_counter()
        for epoch in range(start_epoch, end_epoch):
            # This is a line-for-line numerical transcription of the official
            # PINN.train_adam body.  Progress-bar display is intentionally not
            # part of the update sequence.
            model.pgd()
            if args.method == 0:
                assert args.batch_size == args.dim
                loss = model.get_loss_pinn()
                saved_loss = loss
            elif args.method == 1:
                assert args.batch_size < args.dim
                loss, saved_loss = getattr(
                    model, "get_loss_pinn_stochastic_" + args.dataset
                )(args.batch_size)
            elif args.method == 2:
                assert args.batch_size < args.dim
                loss, saved_loss = getattr(
                    model, "get_loss_pinn_stochastic2_" + args.dataset
                )(args.batch_size)
            else:
                raise ValueError(args.method)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            if args.use_sch:
                scheduler.step()
            current_loss = saved_loss.item()
            model.Resample_HJB_Log()
            if epoch % 100 == 0:
                current_l2 = model.L2_pinn()
                print(
                    "epoch %d, loss: %e, l2: %e"
                    % (epoch, current_loss, current_l2),
                    flush=True,
                )
            if args.save_loss:
                model.saved_loss.append(current_loss)
                if epoch % 100 != 0:
                    current_l2 = model.L2_pinn()
                model.saved_l2.append(current_l2)
        torch.cuda.synchronize()
        chunk_loop_seconds = time.perf_counter() - loop_started
        cumulative_loop_seconds = previous_loop_seconds + chunk_loop_seconds
        peak_gpu_bytes = max(previous_peak_gpu_bytes, torch.cuda.max_memory_allocated())
        process_seconds = time.perf_counter() - invocation_started
        cumulative_process_seconds = previous_process_seconds + process_seconds
        payload = {
            "seed": seed,
            "tag": tag,
            "next_epoch": end_epoch,
            "total_epochs": total_epochs,
            "model_state_dict": {
                key: value.detach().cpu().clone()
                for key, value in model.u_net.state_dict().items()
            },
            "optimizer_state_dict": optimizer.state_dict(),
            "scheduler_state_dict": scheduler.state_dict(),
            "xf": model.xf.detach().cpu(),
            "saved_loss": list(model.saved_loss),
            "saved_l2": list(model.saved_l2),
            "current_l2": current_l2,
            "rng_state": _rng_state(),
            "cumulative_loop_seconds": cumulative_loop_seconds,
            "cumulative_process_seconds": cumulative_process_seconds,
            "peak_gpu_bytes": peak_gpu_bytes,
        }
        torch.save(payload, continuation)

        result: dict[str, Any] = {
            "method": "sdgd_pinn_official_resumable_compatibility_port",
            "benchmark": "hjb_log_d100",
            "seed": seed,
            "tag": tag,
            "start_epoch": start_epoch,
            "end_epoch": end_epoch,
            "total_epochs": total_epochs,
            "chunk_loop_seconds": chunk_loop_seconds,
            "chunk_process_seconds": process_seconds,
            "cumulative_loop_seconds": cumulative_loop_seconds,
            "cumulative_process_seconds": cumulative_process_seconds,
            "peak_gpu_bytes": peak_gpu_bytes,
            "continuation_checkpoint": str(continuation),
            "official_source": str(SOURCE),
            "official_source_commit": SOURCE_COMMIT,
            "numerical_loop_change": False,
            "split_boundary_state": (
                "model, Adam, LambdaLR, collocation points, NumPy RNG, Torch CPU RNG, and all CUDA RNG states"
            ),
            "terminal_aware_feature_added_by_us": False,
            "mcv_added_by_us": False,
            "status": "complete",
        }
        atomic_json(
            RAW_DIR / f"official_sdgd_resumable_s{seed}_{tag}_e{end_epoch}.json",
            result,
        )
        if end_epoch == total_epochs:
            checkpoint = work / f"sdgd_hjb_log_s{seed}_official.pt"
            torch.save(model.u_net.state_dict(), checkpoint)
            if args.save_loss:
                frame = namespace["pd"].DataFrame(
                    {
                        "loss": np.asarray(model.saved_loss),
                        "L2": np.asarray(model.saved_l2),
                    }
                )
                frame.to_excel(
                    work
                    / "saved_loss_l2"
                    / (
                        f"{args.dataset}_Linf_dim={args.dim}"
                        f"_batch_pgd={args.batch_size_pgd}"
                        f"_method_pgd={args.method_adv}"
                        f"_batch={args.batch_size}"
                        f"_method={args.method}_SEED={args.SEED}.xlsx"
                    ),
                    index=False,
                )
            result.update(
                {
                    "method": "sdgd_pinn_official",
                    "calibration": False,
                    "epochs": total_epochs,
                    "training_seconds_including_official_per_epoch_test_reporting": cumulative_process_seconds,
                    "training_loop_seconds": cumulative_loop_seconds,
                    "checkpoint": str(checkpoint),
                    "parameter_count": sum(
                        parameter.numel() for parameter in model.u_net.parameters()
                    ),
                    "published_defaults": {
                        "width": 1024,
                        "hidden_layers": 3,
                        "activation": "tanh",
                        "collocation_points": 100,
                        "pgd_steps_per_epoch": 20,
                        "pgd_step_size": 0.05,
                        "dimension_samples": 10,
                        "optimizer": "Adam",
                        "learning_rate": 0.001,
                        "linear_lr_schedule": True,
                        "hard_terminal_ansatz": True,
                        "method": 2,
                        "adversarial_method": 2,
                    },
                }
            )
            atomic_json(RAW_DIR / f"official_sdgd_s{seed}_official.json", result)
    finally:
        os.chdir(previous)
    return result


def compare_validation(seed: int, tag: str) -> dict[str, Any]:
    work = OFFICIAL_WORK_DIR / "sdgd"
    continuation = torch.load(
        work / f"sdgd_hjb_log_s{seed}_{tag}_continuation.pt",
        map_location="cpu",
        weights_only=False,
    )
    official = torch.load(
        work / f"sdgd_hjb_log_s{seed}_calibration.pt",
        map_location="cpu",
        weights_only=False,
    )
    differences = {
        key: float((continuation["model_state_dict"][key] - official[key]).abs().max())
        for key in official
    }
    result = {
        "seed": seed,
        "tag": tag,
        "compared_epochs": int(continuation["next_epoch"]),
        "maximum_parameter_absolute_difference": max(differences.values()),
        "all_parameters_bitwise_equal": all(value == 0.0 for value in differences.values()),
        "per_tensor_maximum_absolute_difference": differences,
    }
    atomic_json(RESULTS_DIR / "sdgd_resumable_equivalence_audit.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--total-epochs", type=int, default=10_000)
    parser.add_argument("--tag", default="official10000")
    parser.add_argument("--compare-validation", action="store_true")
    args = parser.parse_args()
    if args.compare_validation:
        result = compare_validation(args.seed, args.tag)
    else:
        if args.start is None or args.end is None:
            parser.error("--start and --end are required for a training chunk")
        result = run_chunk(
            seed=args.seed,
            start_epoch=args.start,
            end_epoch=args.end,
            total_epochs=args.total_epochs,
            tag=args.tag,
        )
    print(json.dumps(result, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
