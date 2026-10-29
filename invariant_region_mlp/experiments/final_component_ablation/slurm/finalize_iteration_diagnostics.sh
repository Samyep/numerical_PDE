#!/bin/bash
set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
RESULT_ROOT="$SOURCE_ROOT/invariant_region_mlp/results/final_component_ablation"
PYTHON=/sw/user/python/miniforge3-pytorch-2.10.0/bin/python
ACTIVE_JOB_IDS=(3484104 3484142)
ACCOUNTING_JOB_IDS=(
  3474284 3474354 3474773 3474794 3474893 3474277 3476461
  3480678 3480679 3484104 3484120 3484142
)

cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export MPLBACKEND=Agg

while true; do
  complete=$(
    "$PYTHON" - <<'PY'
import json
from pathlib import Path
root = Path("invariant_region_mlp/results/final_component_ablation/raw")
count = 0
for path in root.glob("iteration_diagnostic_*_d100_a0_full_s*_a*_rep0.json"):
    try:
        count += json.loads(path.read_text()).get("status") == "complete"
    except (OSError, json.JSONDecodeError):
        pass
print(count)
PY
  )
  active=$(
    (squeue -h -j "$(IFS=,; echo "${ACTIVE_JOB_IDS[*]}")" \
      -o '%i:%T' 2>/dev/null || true) | paste -sd, -
  )
  printf '%s complete=%s/60 active=%s\n' \
    "$(date --iso-8601=seconds)" "$complete" "${active:-none}"
  if [[ "$complete" == "60" && -z "$active" ]]; then
    break
  fi
  if [[ -z "$active" && "$complete" != "60" ]]; then
    echo "all diagnostic jobs are terminal but only $complete/60 results exist" >&2
    sacct -j "$(IFS=,; echo "${ACTIVE_JOB_IDS[*]}")" \
      --format=JobID,JobName%24,State,Elapsed,ExitCode,AllocTRES -P
    exit 1
  fi
  sleep 60
done

job_ids=$(IFS=,; echo "${ACCOUNTING_JOB_IDS[*]}")
"$PYTHON" -m invariant_region_mlp.experiments.final_component_ablation.build_artifacts \
  --job-ids "$job_ids"
"$PYTHON" -m invariant_region_mlp.experiments.final_component_ablation.test_correctness
"$PYTHON" -m py_compile \
  invariant_region_mlp/experiments/final_component_ablation/build_artifacts.py \
  invariant_region_mlp/experiments/final_component_ablation/run_iteration_diagnostics.py
git diff --check
echo "FINALIZE_COMPLETE $(date --iso-8601=seconds)"
