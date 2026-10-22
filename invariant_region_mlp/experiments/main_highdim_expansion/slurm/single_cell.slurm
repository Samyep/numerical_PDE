#!/bin/bash
#SBATCH --job-name=icml_exp_cell
#SBATCH --account=isye
#SBATCH --partition=ice-gpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=96G
#SBATCH --time=06:10:00
#SBATCH --output=invariant_region_mlp/results/main_highdim_expansion/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/main_highdim_expansion/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from the repository root}}
SCRIPT_DIR=$SOURCE_ROOT/invariant_region_mlp/experiments/main_highdim_expansion/slurm
source "$SCRIPT_DIR/common_env.sh"

method=${METHOD:?submit with METHOD}
pde=${PDE:?submit with PDE}
dimension=${DIMENSION:?submit with DIMENSION}

python -m invariant_region_mlp.experiments.main_highdim_expansion.review_gate
python -m invariant_region_mlp.experiments.main_highdim_expansion.accounting \
  --reserve-job-id "$SLURM_JOB_ID"
nvidia-smi --query-gpu=index,uuid,name,memory.total,driver_version --format=csv,noheader
python -c 'import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available(), torch.cuda.get_device_name(0))'

exec "$SCRIPT_DIR/run_cell.sh" "$method" "$pde" "$dimension"
