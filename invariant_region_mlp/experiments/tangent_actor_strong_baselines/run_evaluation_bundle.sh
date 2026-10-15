#!/bin/bash
set -euo pipefail

BENCHMARK=$1
SEED=$2
PYTHON_BIN=/sw/user/python/miniforge3-pytorch-2.10.0/bin/python

export OMP_NUM_THREADS=12
"$PYTHON_BIN" -c "import torch; print('evaluation bundle', '$BENCHMARK', '$SEED', 'device', torch.cuda.get_device_name(0), flush=True)"
"$PYTHON_BIN" -u -m invariant_region_mlp.experiments.tangent_actor_strong_baselines.evaluate_run \
  --benchmark "$BENCHMARK" --seed "$SEED"
