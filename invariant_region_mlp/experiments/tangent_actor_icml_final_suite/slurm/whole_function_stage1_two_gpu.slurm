#!/bin/bash
#SBATCH --job-name=icml_whole_pilot
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=2
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=220G
#SBATCH --time=01:30:00
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
python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.review_gate
python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_job_accounting \
  --reserve-job-id "$SLURM_JOB_ID"
nvidia-smi --query-gpu=index,name,memory.total,driver_version --format=csv,noheader
python -c 'import torch; print("torch.cuda.is_available=", torch.cuda.is_available(), "count=", torch.cuda.device_count())'

slow_out="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_ID}_slow_sdgd_d${DIMENSION}.out"
slow_err="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_ID}_slow_sdgd_d${DIMENSION}.err"
fast_out="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_ID}_fast_lane_d${DIMENSION}.out"
fast_err="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_ID}_fast_lane_d${DIMENSION}.err"
helper="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/run_whole_function_cell.sh"
fast_worker="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/run_stage1_fast_lane.sh"

srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
  --mem=110G --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
  --output="$slow_out" --error="$slow_err" \
  env WHOLE_CELL_TASK_ID=slow CELL_CAP_SECONDS=5100 \
  "$helper" sdgd_pinn joint_concave "$DIMENSION" &
slow_pid=$!
srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
  --mem=110G --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
  --output="$fast_out" --error="$fast_err" "$fast_worker" "$DIMENSION" &
fast_pid=$!

status=0
wait "$slow_pid" || status=1
wait "$fast_pid" || status=1
exit "$status"
