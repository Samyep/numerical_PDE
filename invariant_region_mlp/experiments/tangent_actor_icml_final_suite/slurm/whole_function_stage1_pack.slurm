#!/bin/bash
#SBATCH --job-name=icml_whole_pilot
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=2
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=220G
#SBATCH --time=02:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
DIMENSION=${DIMENSION:?submit with DIMENSION=20 or 50}
case "$DIMENSION" in 20|50) ;; *) echo "pilot dimension must be 20 or 50" >&2; exit 2;; esac

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.review_gate || exit $?
python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_job_accounting \
  --reserve-job-id "$SLURM_JOB_ID" || exit $?
nvidia-smi --query-gpu=index,name,memory.total,driver_version --format=csv,noheader
python -c 'import torch; print("torch.cuda.is_available=", torch.cuda.is_available(), "count=", torch.cuda.device_count())'
pilot_cell_cap=$(python - "$SLURM_JOB_ID" <<'PY'
import re
import subprocess
import sys
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_job_accounting import _duration_seconds
info = subprocess.check_output(["scontrol", "show", "job", "-o", sys.argv[1]], text=True)
match = re.search(r"\bTimeLimit=(\S+)", info)
if not match:
    raise SystemExit("missing Slurm TimeLimit")
print(max(300, _duration_seconds(match.group(1)) - 300))
PY
) || exit $?
echo "pilot cell cap from allocation: ${pilot_cell_cap}s"

methods=(tangent_actor dpi pinn_hte sdgd_pinn)
pids=()
for slot in 0 1 2 3; do
  method=${methods[$slot]}
  task_out="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_ID}_${slot}_${method}_joint_concave_d${DIMENSION}_step.out"
  task_err="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_ID}_${slot}_${method}_joint_concave_d${DIMENSION}_step.err"
  srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    --mem=110G --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
    --output="$task_out" --error="$task_err" \
    env WHOLE_CELL_TASK_ID="$slot" CELL_CAP_SECONDS="$pilot_cell_cap" \
    "$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/run_whole_function_cell.sh" \
    "$method" joint_concave "$DIMENSION" &
  pids+=("$!")
done

status=0
for slot in 0 1 2 3; do
  if wait "${pids[$slot]}"; then code=0; else code=$?; fi
  if [[ $code -ne 0 && $code -ne 124 ]]; then
    echo "pilot slot failed: slot=$slot method=${methods[$slot]} dimension=$DIMENSION exit=$code" >&2
    status=1
  fi
done
exit "$status"
