#!/bin/bash
set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
TASK_FILE=${1:?usage: run_queue_worker.sh TASK_FILE SLOT}
SLOT=${2:?usage: run_queue_worker.sh TASK_FILE SLOT}

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-12}

tail -n +2 "$TASK_FILE" | while IFS=',' read -r slot sequence pde variant seed estimate extra; do
  [[ "$slot" == "$SLOT" ]] || continue
  if [[ -n "${extra:-}" ]]; then
    echo "invalid queue row slot=$slot sequence=$sequence" >&2
    exit 2
  fi
  rep0="$SOURCE_ROOT/invariant_region_mlp/results/final_component_ablation/raw/evaluation_${pde}_d100_${variant}_s${seed}_rep0.json"
  rep1="$SOURCE_ROOT/invariant_region_mlp/results/final_component_ablation/raw/evaluation_${pde}_d100_${variant}_s${seed}_rep1.json"
  if [[ -s "$rep0" && -s "$rep1" ]]; then
    echo "CELL_SKIP slot=$SLOT seq=$sequence pde=$pde variant=$variant seed=$seed reason=complete"
    continue
  fi
  echo "CELL_BEGIN slot=$SLOT seq=$sequence pde=$pde variant=$variant seed=$seed estimate=$estimate utc=$(date -u +%FT%TZ)"
  common=(--pde "$pde" --variant "$variant" --seed "$seed" --device cuda)
  if [[ "$variant" == "a3_no_updates" ]]; then
    python -m invariant_region_mlp.experiments.final_component_ablation.run_ablation \
      "${common[@]}" --replicate 0 --evaluate-only
  else
    python -m invariant_region_mlp.experiments.final_component_ablation.run_ablation \
      "${common[@]}" --replicate 0
  fi
  python -m invariant_region_mlp.experiments.final_component_ablation.run_ablation \
    "${common[@]}" --replicate 1 --evaluate-only
  echo "CELL_END slot=$SLOT seq=$sequence pde=$pde variant=$variant seed=$seed utc=$(date -u +%FT%TZ)"
done
