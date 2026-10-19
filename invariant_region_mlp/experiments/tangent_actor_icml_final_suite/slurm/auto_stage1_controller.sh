#!/bin/bash
# Lightweight login-node controller for the two approved Stage-1 pilot packs.
set -uo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
PILOT_SCRIPT="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/whole_function_stage1_pack.slurm"
RESULTS_DIR="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite"
first_job=${1:?usage: auto_stage1_controller.sh FIRST_JOB_ID}

source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"

wait_for_job() {
  local job=$1
  while squeue -h -j "$job" | grep -q .; do
    echo "$(date -u +%FT%TZ) waiting job=$job $(squeue -h -j "$job" -o '%T|%M|%S|%R')"
    sleep 30
  done
  for _ in $(seq 1 10); do
    state=$(sacct -X -n -j "$job" -o State | awk 'NF {print $1; exit}')
    [[ -n "${state:-}" ]] && break
    sleep 10
  done
  echo "$(date -u +%FT%TZ) finished job=$job state=${state:-UNKNOWN}"
}

refresh_artifacts() {
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.aggregate_results
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_job_accounting
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_artifacts
  python - <<'PY'
import json
from pathlib import Path
p = Path("invariant_region_mlp/results/tangent_actor_icml_final_suite/WHOLE_FUNCTION_AUDIT.json")
d = json.loads(p.read_text())
if d["status"] != "PASS" or not d["checks"]["new_compute_within_120_gpu_hours"]:
    raise SystemExit("whole-function audit did not pass")
PY
}

require_four_complete() {
  local job=$1 dimension=$2
  python - "$job" "$dimension" <<'PY'
import sys
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_status import attempts
job, dimension = sys.argv[1], int(sys.argv[2])
rows = [r for r in attempts() if str(r.get("job_id")) == job and int(r.get("dimension", -1)) == dimension]
statuses = {r.get("method"): r.get("status") for r in rows}
expected = {"tangent_actor", "dpi", "pinn_hte", "sdgd_pinn"}
print({"job": job, "dimension": dimension, "statuses": statuses})
if set(statuses) != expected or any(statuses[name] != "COMPLETE" for name in expected):
    raise SystemExit("pilot did not produce four COMPLETE cells; stopping staged dispatch")
PY
}

wait_for_job "$first_job"
refresh_artifacts || exit $?
require_four_complete "$first_job" 20 || exit $?

if squeue -h -u "$USER" -o '%j' | grep -q '^icml_whole'; then
  echo "duplicate protection: another icml_whole job is active" >&2
  exit 6
fi
sbatch --test-only --time=01:30:00 --export=ALL,DIMENSION=50 "$PILOT_SCRIPT" || exit $?
second_job=$(sbatch --parsable --time=01:30:00 --export=ALL,DIMENSION=50 "$PILOT_SCRIPT") || exit $?
echo "$second_job" > "$RESULTS_DIR/STAGE1_D50_JOB_ID.txt"
echo "$(date -u +%FT%TZ) submitted d50 job=$second_job"
wait_for_job "$second_job"
refresh_artifacts || exit $?
require_four_complete "$second_job" 50 || exit $?
pytest -q invariant_region_mlp/experiments/tangent_actor_icml_final_suite
echo "$(date -u +%FT%TZ) Stage 1 complete jobs=$first_job,$second_job"
