#!/bin/bash
set -euo pipefail

BENCHMARK=$1
SEED=$2
PYTHON_BIN=/sw/user/python/miniforge3-pytorch-2.10.0/bin/python

export OMP_NUM_THREADS=12
"$PYTHON_BIN" -c "import torch; print('evaluation bundle', '$BENCHMARK', '$SEED', 'device', torch.cuda.get_device_name(0), flush=True)"

# The original fourth Slurm step encountered a node-local VNI launch failure
# and was restarted 75 seconds late in the same allocation.  It completed all
# expensive continuous baselines but the allocation ended during its final
# Deep-Picard run.  Repair exactly that missing raw record before evaluation;
# the deterministic seed/configuration and method-internal timing are unchanged.
DPI_RAW=invariant_region_mlp/results/tangent_actor_strong_baselines/raw/deep_picard_${BENCHMARK}_s${SEED}_main.json
if [[ "$BENCHMARK" == "nonseparable_d100_k2" && "$SEED" == "1" && ! -f "$DPI_RAW" ]]; then
  echo "repairing missing Deep Picard training record for $BENCHMARK seed=$SEED"
  "$PYTHON_BIN" -u -m invariant_region_mlp.experiments.tangent_actor_strong_baselines.run_training \
    --method deep_picard --benchmark "$BENCHMARK" --seed "$SEED"
fi

"$PYTHON_BIN" -u -m invariant_region_mlp.experiments.tangent_actor_strong_baselines.evaluate_run \
  --benchmark "$BENCHMARK" --seed "$SEED"
