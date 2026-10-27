#!/bin/bash
set -euo pipefail
script_dir=$(cd "$(dirname "$0")" && pwd)
source "$script_dir/common_env.sh"
python - <<'PY'
import json
from invariant_region_mlp.experiments.deep_bsde_budget_scaling_3seed.protocol import RESULTS_DIR
path = RESULTS_DIR / "preflight_allocation_gate.json"
payload = json.loads(path.read_text())
if not payload.get("formal_submission_authorized", False):
    raise SystemExit("formal submission blocked by allocation gate")
print(payload["projected_total_gpu_hours_with_contingency"])
PY
mkdir -p invariant_region_mlp/results/deep_bsde_budget_scaling_3seed/slurm
job_id=$(sbatch --parsable "$script_dir/formal.slurm")
python - "$job_id" <<'PY'
import json, sys
from datetime import datetime, timezone
from invariant_region_mlp.experiments.deep_bsde_budget_scaling_3seed.protocol import RESULTS_DIR, write_json_atomic
write_json_atomic(RESULTS_DIR / "submitted_jobs.json", {
    "formal_array_job_id": sys.argv[1],
    "submitted_utc": datetime.now(timezone.utc).isoformat(),
})
PY
echo "$job_id"

