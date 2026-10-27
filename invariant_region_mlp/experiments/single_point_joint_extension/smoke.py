"""Small end-to-end infrastructure gate before formal GPU dispatch."""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import torch

from invariant_region_mlp.experiments.main_highdim_expansion.joint_actor import (
    JointTerminalActor,
)
from invariant_region_mlp.experiments.main_highdim_expansion.joint_training import (
    fit_joint_actor,
)
from invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.solver import (
    evaluate_policy,
)
from invariant_region_mlp.experiments.tangent_actor_mcv_training.labels import (
    mcv_tangent_labels,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
)

from .benchmark import make_benchmark, sample_states
from .deep_bsde import execute as deep_execute
from .fh_mlp import execute as fh_execute
from .protocol import (
    PREFLIGHT_AUDIT,
    REVISION,
    SMOKE_AUDIT,
    require_preflight,
    stable_seed,
    write_json_atomic,
)


def _ta_smoke(device: torch.device) -> dict:
    driver = make_benchmark(100)
    torch.manual_seed(stable_seed("smoke", "ta_actor"))
    source = JointTerminalActor(driver, width=16, depth=2).to(device)
    source.eval()
    train_t, train_x = sample_states(100, 16, seed=stable_seed("smoke", "train"))
    valid_t, valid_x = sample_states(100, 8, seed=stable_seed("smoke", "valid"))
    kwargs = {
        "paths": 4,
        "time_steps": 2,
        "device": device,
        "parent_batch": 4,
        "path_chunk": 4,
    }
    train_r, train_q, train_meta = mcv_tangent_labels(
        driver, source, train_t, train_x, seed=stable_seed("smoke", "train_labels"), **kwargs
    )
    replay_r, replay_q, _ = mcv_tangent_labels(
        driver, source, train_t, train_x, seed=stable_seed("smoke", "train_labels"), **kwargs
    )
    valid_r, valid_q, _ = mcv_tangent_labels(
        driver, source, valid_t, valid_x, seed=stable_seed("smoke", "valid_labels"), **kwargs
    )
    config = FitConfig(
        training_points=16,
        validation_points=8,
        paths_per_label=4,
        label_time_steps=2,
        optimizer_steps=2,
        width=16,
        depth=2,
        learning_rate=1.5e-4,
        batch_size=8,
        parent_batch=4,
        path_chunk=4,
        checkpoint_interval=1,
        gradient_clip=10.0,
    )
    actor, fit_meta, _ = fit_joint_actor(
        driver,
        config,
        source,
        train_t,
        train_x,
        train_r,
        train_q,
        valid_t,
        valid_x,
        valid_r,
        valid_q,
        seed=stable_seed("smoke", "fit"),
        device=device,
    )
    output, eval_meta = evaluate_policy(
        driver,
        actor,
        valid_t[:4],
        valid_x[:4],
        paths=4,
        time_steps=2,
        seed=stable_seed("smoke", "evaluate"),
        device=device,
        parent_batch=4,
        path_chunk=4,
    )
    replay_exact = np.array_equal(train_r, replay_r) and np.array_equal(train_q, replay_q)
    finite = all(np.all(np.isfinite(value)) for value in (train_r, train_q, output["mcv"]))
    return {
        "label_replay_bitwise_equal": replay_exact,
        "finite": finite,
        "fit_selected_step": fit_meta["selected_step"],
        "label_nonfinite_count": train_meta["nonfinite_count"],
        "evaluation_nonfinite_count": eval_meta["nonfinite_count"],
        "passed": replay_exact and finite and eval_meta["nonfinite_count"] == 0,
    }


def run(device: torch.device) -> dict:
    require_preflight()
    started = time.perf_counter()
    checks = {
        "tangent_actor": _ta_smoke(device),
        "fh_mlp": fh_execute(100, 2, 4, 0, device, pilot=True),
        "deep_bsde": deep_execute(100, 0, 0, device, pilot=True),
    }
    passed = bool(
        checks["tangent_actor"]["passed"]
        and checks["fh_mlp"]["status"] == "pilot_complete"
        and checks["deep_bsde"]["status"] == "pilot_complete"
        and checks["fh_mlp"]["nonfinite_count"] == 0
        and checks["deep_bsde"]["nonfinite_count"] == 0
    )
    payload = {
        "revision": REVISION,
        "status": "VALIDATED" if passed else "FAILED",
        "overall_pass": passed,
        "device": str(device),
        "preflight_sha256": __import__("hashlib").sha256(
            PREFLIGHT_AUDIT.read_bytes()
        ).hexdigest(),
        "checks": checks,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(SMOKE_AUDIT, payload)
    if not passed:
        raise RuntimeError("joint-extension smoke gate failed")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    result = run(device)
    print(json.dumps({
        "status": result["status"],
        "elapsed_seconds": result["elapsed_seconds"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["run"]
