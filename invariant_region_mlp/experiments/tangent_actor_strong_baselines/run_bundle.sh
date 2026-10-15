#!/bin/bash
set -euo pipefail

BENCHMARK=$1
SEED=$2
PYTHON_BIN=/sw/user/python/miniforge3-pytorch-2.10.0/bin/python
METHODS=(tangent_actor standard_tangent_actor deep_bsde sdgd_pinn deep_picard)

export OMP_NUM_THREADS=12
"$PYTHON_BIN" -c "import torch; print('bundle', '$BENCHMARK', '$SEED', 'device', torch.cuda.get_device_name(0), flush=True)"
for METHOD in "${METHODS[@]}"; do
  echo "starting method=$METHOD benchmark=$BENCHMARK seed=$SEED"
  "$PYTHON_BIN" -u -m invariant_region_mlp.experiments.tangent_actor_strong_baselines.run_training \
    --method "$METHOD" --benchmark "$BENCHMARK" --seed "$SEED"
done
