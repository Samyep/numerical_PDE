#!/bin/bash
#SBATCH --partition=ice-gpu
#SBATCH --account=isye
#SBATCH --qos=coe-ice
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --job-name=wcf_eval
set -euo pipefail
source "$(dirname "$0")/slurm_env.sh"

pdes=(b6_fullrank b6_fullrank b6_fullrank b6_fullrank b6_fullrank b6_fullrank b6_fullrank b6_fullrank b6_fullrank dense_fourier dense_fourier dense_fourier dense_fourier dense_fourier dense_fourier dense_fourier dense_fourier dense_fourier)
methods=(tangent_actor tangent_actor tangent_actor dpi dpi dpi pinn_hte pinn_hte pinn_hte tangent_actor tangent_actor tangent_actor dpi dpi dpi pinn_hte pinn_hte pinn_hte)
configs=(C0 C1 C2 C0 C1 C2 C0 C1 C2 C0 C1 C2 C0 C1 C2 C0 C1 C2)
index="${SLURM_ARRAY_TASK_ID:?array task id is required}"
python -m invariant_region_mlp.experiments.whole_function_compute_frontier.run_frontier evaluate \
  --method "${methods[$index]}" --pde "${pdes[$index]}" \
  --configuration "${configs[$index]}" --device cuda
