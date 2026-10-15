#!/bin/bash
#SBATCH --job-name=ta_baseline_main
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=01:00:00
#SBATCH --array=0-19%4
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_strong_baselines/slurm/%x_%A_%a.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_strong_baselines/slurm/%x_%A_%a.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
PYTHON_BIN=/sw/user/python/miniforge3-pytorch-2.10.0/bin/python
METHODS=(tangent_actor standard_tangent_actor deep_bsde sdgd_pinn deep_picard)
BENCHMARKS=(nonseparable_d100_k2 joint_concave_d400_k2)

METHOD_INDEX=$((SLURM_ARRAY_TASK_ID / 4))
REMAINDER=$((SLURM_ARRAY_TASK_ID % 4))
BENCHMARK_INDEX=$((REMAINDER / 2))
SEED=$((REMAINDER % 2))
METHOD=${METHODS[$METHOD_INDEX]}
BENCHMARK=${BENCHMARKS[$BENCHMARK_INDEX]}

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

"$PYTHON_BIN" -V
git rev-parse HEAD
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
echo "task=$SLURM_ARRAY_TASK_ID method=$METHOD benchmark=$BENCHMARK seed=$SEED"
"$PYTHON_BIN" -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'device', torch.cuda.get_device_name(0))"

srun "$PYTHON_BIN" -u -m invariant_region_mlp.experiments.tangent_actor_strong_baselines.run_training \
  --method "$METHOD" --benchmark "$BENCHMARK" --seed "$SEED"
