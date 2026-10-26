#!/bin/bash
#SBATCH --partition=ice-gpu
#SBATCH --account=isye
#SBATCH --qos=coe-ice
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=01:15:00
set -euo pipefail
source "${WCF_SCRIPT_DIR:-${SLURM_SUBMIT_DIR}/invariant_region_mlp/experiments/whole_function_compute_frontier}/slurm_env.sh"

stage="${1:?stage must be c0 or sensitivity}"
if [[ "${stage}" == "c0" ]]; then
  pdes=(b6_fullrank b6_fullrank b6_fullrank dense_fourier dense_fourier dense_fourier)
  methods=(tangent_actor dpi pinn_hte tangent_actor dpi pinn_hte)
  configs=(C0 C0 C0 C0 C0 C0)
elif [[ "${stage}" == "sensitivity" ]]; then
  pdes=(b6_fullrank b6_fullrank b6_fullrank b6_fullrank b6_fullrank b6_fullrank dense_fourier dense_fourier dense_fourier dense_fourier dense_fourier dense_fourier)
  methods=(tangent_actor tangent_actor dpi dpi pinn_hte pinn_hte tangent_actor tangent_actor dpi dpi pinn_hte pinn_hte)
  configs=(C1 C2 C1 C2 C1 C2 C1 C2 C1 C2 C1 C2)
else
  echo "unknown stage: ${stage}" >&2
  exit 2
fi
index="${SLURM_ARRAY_TASK_ID:?array task id is required}"
python -m invariant_region_mlp.experiments.whole_function_compute_frontier.run_frontier train \
  --method "${methods[$index]}" --pde "${pdes[$index]}" \
  --configuration "${configs[$index]}" --device cuda
