#!/bin/bash
#SBATCH --job-name=orig_actor_eval
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=02:00:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
: "${BENCHMARK:?set BENCHMARK to hjb_log_d100, dpi_ou_gmm_d100, or allen_cahn_d100}"
: "${SEEDS:=0 1 2}"
: "${DIAGNOSTIC_STEPS=32 128}"
: "${RUN_SCREEN:=1}"
DECISION=$SOURCE_ROOT/invariant_region_mlp/results/original_pde_benchmark_comparison/final_actor_selection_${BENCHMARK}.json

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

read -r TRIAL ITERATION < <(python -c 'import json,sys; d=json.load(open(sys.argv[1]))["selected"]; print(d["trial"], d["iteration"])' "$DECISION")

python -V
git rev-parse HEAD
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"

read -r -a SEED_LIST <<< "$SEEDS"
for SEED in "${SEED_LIST[@]}"; do
  srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.evaluate_actor \
    --benchmark "$BENCHMARK" --trial "$TRIAL" --seed "$SEED" --iteration "$ITERATION" \
    --paths 256 --steps 64 --replicates 4 --native-paths 8192
done

read -r -a DIAGNOSTIC_STEP_LIST <<< "$DIAGNOSTIC_STEPS"
for STEPS in "${DIAGNOSTIC_STEP_LIST[@]}"; do
  srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.evaluate_actor \
    --benchmark "$BENCHMARK" --trial "$TRIAL" --seed 0 --iteration "$ITERATION" \
    --paths 256 --steps "$STEPS" --replicates 4 --native-paths 8192
done

if [[ "$RUN_SCREEN" == "1" ]]; then
  srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.evaluate_actor_screen \
    --benchmark "$BENCHMARK" --paths 128 --steps 64
fi
