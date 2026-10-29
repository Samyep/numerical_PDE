#!/bin/bash
#SBATCH --job-name=ta_iter_diag_2g
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=2
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/final_component_ablation/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/final_component_ablation/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
RESULT_ROOT="$SOURCE_ROOT/invariant_region_mlp/results/final_component_ablation"
REGULAR_JOB_ID=3484104

# Keep whichever allocation starts first.  Cancel only a still-pending regular
# twin; never interrupt it after it has begun useful work.
regular_state=$(squeue -h -j "$REGULAR_JOB_ID" -o '%T' || true)
if [[ -n "$regular_state" && "$regular_state" != "PENDING" ]]; then
  echo "regular twin $REGULAR_JOB_ID is already $regular_state; exiting"
  exit 0
fi
if [[ "$regular_state" == "PENDING" ]]; then
  scancel --state=PENDING "$REGULAR_JOB_ID" || true
fi
regular_state=$(squeue -h -j "$REGULAR_JOB_ID" -o '%T' || true)
if [[ -n "$regular_state" ]]; then
  echo "regular twin $REGULAR_JOB_ID remains $regular_state; exiting"
  exit 0
fi

cd "$SOURCE_ROOT"
module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

echo "git_head=$(git rev-parse HEAD)"
echo "job_id=$SLURM_JOB_ID nodes=$SLURM_JOB_NUM_NODES tasks=$SLURM_NTASKS"
python -V
nvidia-smi --query-gpu=index,name,memory.total,driver_version --format=csv,noheader
python -c 'import torch; print("torch", torch.__version__, "cuda", torch.cuda.is_available(), "count", torch.cuda.device_count())'
python -m invariant_region_mlp.experiments.final_component_ablation.preflight
python -m invariant_region_mlp.experiments.final_component_ablation.run_iteration_diagnostics --audit-only

chains=("hjb_log nonconvex_coupled" "joint_concave joint_convex")
pids=()
for slot in 0 1; do
  slot_out="$RESULT_ROOT/slurm/${SLURM_JOB_ID}_iteration_slot${slot}.out"
  slot_err="$RESULT_ROOT/slurm/${SLURM_JOB_ID}_iteration_slot${slot}.err"
  srun --exclusive --exact --nodes=1 --ntasks=1 \
    --cpus-per-task="$SLURM_CPUS_PER_TASK" --mem=60G \
    --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
    --output="$slot_out" --error="$slot_err" \
    bash -lc '
      set -euo pipefail
      for pde in $1; do
        python -m invariant_region_mlp.experiments.final_component_ablation.run_iteration_diagnostics \
          --pde "$pde" --device cuda
      done
    ' _ "${chains[$slot]}" &
  pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
  wait "$pid" || status=1
done
if (( status != 0 )); then
  exit "$status"
fi

python -c '
import json
from pathlib import Path
root = Path("invariant_region_mlp/results/final_component_ablation/raw")
paths = sorted(root.glob("iteration_diagnostic_*_rep0.json"))
complete = sum(json.loads(path.read_text()).get("status") == "complete" for path in paths)
print({"iteration_diagnostics": len(paths), "complete": complete})
if len(paths) != 60 or complete != 60:
    raise SystemExit(1)
'
