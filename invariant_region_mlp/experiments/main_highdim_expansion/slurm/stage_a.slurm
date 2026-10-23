#!/bin/bash
#SBATCH --job-name=icml_exp_stagea
#SBATCH --account=isye
#SBATCH --partition=ice-gpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=invariant_region_mlp/results/main_highdim_expansion/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/main_highdim_expansion/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from the repository root}}
SCRIPT_DIR=$SOURCE_ROOT/invariant_region_mlp/experiments/main_highdim_expansion/slurm
source "$SCRIPT_DIR/common_env.sh"

python -m invariant_region_mlp.experiments.main_highdim_expansion.accounting --reserve-job-id "$SLURM_JOB_ID"
nvidia-smi --query-gpu=index,uuid,name,memory.total,driver_version --format=csv,noheader || true
python -m invariant_region_mlp.experiments.main_highdim_expansion.infrastructure \
  probe --context stage_a
python -m invariant_region_mlp.experiments.main_highdim_expansion.machine_info

python -m invariant_region_mlp.experiments.main_highdim_expansion.calibrate --device cuda
python -m invariant_region_mlp.experiments.main_highdim_expansion.prepare
python -m invariant_region_mlp.experiments.main_highdim_expansion.reuse_audit
python -m pytest -q invariant_region_mlp/experiments/main_highdim_expansion
python -m invariant_region_mlp.experiments.main_highdim_expansion.validate --device cuda

for method in dpi pinn_hte; do
  for pde in b6_fullrank dense_fourier; do
    python -m invariant_region_mlp.experiments.main_highdim_expansion.baseline_transfer \
      --method "$method" --pde "$pde" --dimension 100 --seed 0 --mode smoke --profile auto
  done
done
for pde in b6_fullrank dense_fourier; do
  python -m invariant_region_mlp.experiments.main_highdim_expansion.sdgd_transfer \
    --family "$pde" --dimension 100 --seed 0 --device cuda --smoke --profile auto
done

python -m invariant_region_mlp.experiments.main_highdim_expansion.characterize \
  --device cuda --dense-heat-pairs 4096
python -m invariant_region_mlp.experiments.main_highdim_expansion.runtime_profile
python -m invariant_region_mlp.experiments.main_highdim_expansion.review_gate

python - <<'PY'
from datetime import datetime, timezone
from invariant_region_mlp.experiments.main_highdim_expansion.protocol import RESULTS_DIR, write_json_atomic
write_json_atomic(RESULTS_DIR / "stage_a_complete.json", {
    "status": "COMPLETE",
    "finished_utc": datetime.now(timezone.utc).isoformat(),
    "job_id": __import__("os").environ["SLURM_JOB_ID"],
})
PY
