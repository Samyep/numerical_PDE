"""Execute one frozen expansion method/PDE/dimension cell."""

from __future__ import annotations

import argparse
import json
import torch

from .protocol import DIMENSIONS, METHODS, PDES, TRAINING_SEED


def execute(method: str, pde: str, dimension: int) -> dict:
    if method not in METHODS or pde not in PDES or dimension not in DIMENSIONS:
        raise ValueError((method, pde, dimension))
    device = torch.device("cuda")
    if not torch.cuda.is_available():
        raise RuntimeError("formal cells require a CUDA allocation")
    if method == "tangent_actor":
        actors = tuple(f"a{index}" for index in range(5)) + ("selected",)
        if pde == "b6_fullrank":
            from .run_b6_actor import evaluate, train

            training = train(pde, dimension, TRAINING_SEED, device)
            primary = evaluate(
                pde,
                dimension,
                TRAINING_SEED,
                device,
                actor_names=actors,
                steps=(64,),
            )
            audit = evaluate(
                pde,
                dimension,
                TRAINING_SEED,
                device,
                actor_names=("selected",),
                steps=(128, 256),
            )
        else:
            from .run_dense_actor import evaluate, train

            training = train(dimension, TRAINING_SEED, device)
            primary = evaluate(
                dimension,
                TRAINING_SEED,
                device,
                actor_names=actors,
                steps=(64,),
            )
            audit = evaluate(
                dimension,
                TRAINING_SEED,
                device,
                actor_names=("selected",),
                steps=(128, 256),
            )
        return {
            "status": "complete",
            "method": method,
            "pde": pde,
            "dimension": dimension,
            "selected_actor": training["selection"]["selected_actor"],
            "primary_rows": len(primary),
            "audit_rows": len(audit),
        }
    if method in {"dpi", "pinn_hte"}:
        from .baseline_transfer import execute as execute_baseline

        result = execute_baseline(
            method, pde, dimension, TRAINING_SEED, "official", "auto"
        )
    else:
        from .sdgd_transfer import execute as execute_sdgd

        result = execute_sdgd(
            pde, dimension, TRAINING_SEED, device, smoke=False, profile="auto"
        )
    return {
        "status": result["status"],
        "method": method,
        "pde": pde,
        "dimension": dimension,
        "value_rmse": result["value_rmse"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--pde", choices=PDES, required=True)
    parser.add_argument("--dimension", type=int, choices=DIMENSIONS, required=True)
    args = parser.parse_args()
    print(json.dumps(execute(args.method, args.pde, args.dimension), sort_keys=True))


if __name__ == "__main__":
    main()
