#!/bin/bash
#SBATCH --job-name=sp32_stagea
#SBATCH --account=isye
#SBATCH --partition=ice-gpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from the experiment checkout root}}
SCRIPT_DIR=$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_32_ice/slurm
source "$SCRIPT_DIR/common_env.sh"

nvidia-smi --query-gpu=index,uuid,name,memory.total,driver_version --format=csv,noheader
python -m pytest -q invariant_region_mlp/experiments/single_point_32_ice/test_single_point_32_ice.py
python -m invariant_region_mlp.experiments.single_point_32_ice.preflight --device cuda

# Small end-to-end CUDA smoke runs. They are namespaced as pilot and are never
# eligible for a formal table.
python -m invariant_region_mlp.experiments.single_point_32_ice.tangent_actor \
  train --pde hjb_log --dimension 20 --seed 0 --device cuda --pilot
python -m invariant_region_mlp.experiments.single_point_32_ice.tangent_actor \
  evaluate --pde hjb_log --dimension 20 --seed 0 --replicate 0 --device cuda --pilot
python -m invariant_region_mlp.experiments.single_point_32_ice.tangent_actor \
  train --pde allen_cahn --dimension 20 --seed 0 --device cuda --pilot
python -m invariant_region_mlp.experiments.single_point_32_ice.tangent_actor \
  evaluate --pde allen_cahn --dimension 20 --seed 0 --replicate 0 --device cuda --pilot
python -m invariant_region_mlp.experiments.single_point_32_ice.tangent_actor \
  train --pde b6_fullrank --dimension 20 --seed 0 --device cuda --pilot
python -m invariant_region_mlp.experiments.single_point_32_ice.tangent_actor \
  evaluate --pde b6_fullrank --dimension 20 --seed 0 --replicate 0 --device cuda --pilot
python -m invariant_region_mlp.experiments.single_point_32_ice.tangent_actor \
  train --pde dense_fourier --dimension 20 --seed 0 --device cuda --pilot
python -m invariant_region_mlp.experiments.single_point_32_ice.tangent_actor \
  evaluate --pde dense_fourier --dimension 20 --seed 0 --replicate 0 --device cuda --pilot
python -m invariant_region_mlp.experiments.single_point_32_ice.deep_bsde \
  --pde hjb_log --dimension 20 --point 0 --seed 0 --device cuda --pilot
python -m invariant_region_mlp.experiments.single_point_32_ice.fh_mlp \
  --pde dense_fourier --dimension 20 --n 2 --M 16 --replicate 0 --device cuda --pilot

python - <<'PY'
from datetime import datetime, timezone
import os
from invariant_region_mlp.experiments.single_point_32_ice.protocol import RESULTS_DIR, write_json_atomic
write_json_atomic(RESULTS_DIR / "stage_a_complete.json", {
    "status": "COMPLETE",
    "finished_utc": datetime.now(timezone.utc).isoformat(),
    "job_id": os.environ["SLURM_JOB_ID"],
})
PY
