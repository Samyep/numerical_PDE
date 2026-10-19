#!/bin/bash
#SBATCH --job-name=icml_whole_resume2
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=2
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=220G
#SBATCH --time=01:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
DIMENSION=${DIMENSION:?submit with DIMENSION}
METHOD_A=${METHOD_A:-tangent_actor}
METHOD_B=${METHOD_B:-sdgd_pinn}
CAP_A_SECONDS=${CAP_A_SECONDS:-3300}
CAP_B_SECONDS=${CAP_B_SECONDS:-3300}

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

helper="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/run_whole_function_cell.sh"
out_a="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_ID}_resume_a_${METHOD_A}__joint_concave__d${DIMENSION}.out"
err_a=${out_a%.out}.err
out_b="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_ID}_resume_b_${METHOD_B}__joint_concave__d${DIMENSION}.out"
err_b=${out_b%.out}.err

srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
  --mem=110G --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
  --output="$out_a" --error="$err_a" \
  env WHOLE_CELL_TASK_ID=resume_a CELL_CAP_SECONDS="$CAP_A_SECONDS" \
  "$helper" "$METHOD_A" joint_concave "$DIMENSION" &
pid_a=$!
srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
  --mem=110G --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
  --output="$out_b" --error="$err_b" \
  env WHOLE_CELL_TASK_ID=resume_b CELL_CAP_SECONDS="$CAP_B_SECONDS" \
  "$helper" "$METHOD_B" joint_concave "$DIMENSION" &
pid_b=$!

status=0
wait "$pid_a" || status=1
wait "$pid_b" || status=1
exit "$status"
