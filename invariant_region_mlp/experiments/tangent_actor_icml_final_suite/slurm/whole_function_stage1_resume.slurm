#!/bin/bash
#SBATCH --job-name=icml_whole_resume
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=120G
#SBATCH --time=00:30:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
METHOD=${METHOD:?submit with METHOD}
PDE=${PDE:?submit with PDE}
DIMENSION=${DIMENSION:?submit with DIMENSION}
RESUME_CAP_SECONDS=${RESUME_CAP_SECONDS:-1500}

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
exec srun --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
  --mem=120G --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
  env WHOLE_CELL_TASK_ID=resume CELL_CAP_SECONDS="$RESUME_CAP_SECONDS" \
  "$helper" "$METHOD" "$PDE" "$DIMENSION"
