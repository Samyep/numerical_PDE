#!/bin/bash
set -euo pipefail
script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(git -C "${script_dir}" rev-parse --show-toplevel)"
result_dir="${repo_root}/invariant_region_mlp/results/whole_function_compute_frontier"
mkdir -p "${result_dir}/slurm"
source "${script_dir}/slurm_env.sh"
python -m invariant_region_mlp.experiments.whole_function_compute_frontier.run_frontier freeze
start="$(date '+%Y-%m-%dT%H:%M:%S')"
if [[ "${WCF_REUSE_SMOKE:-0}" == "1" ]] && python - "${result_dir}/raw/smoke_nonresults.json" <<'PY'
import json, pathlib, sys
path = pathlib.Path(sys.argv[1])
raise SystemExit(0 if path.is_file() and json.loads(path.read_text()).get("status") == "PASS" else 1)
PY
then
  smoke="reused_verified_pass"
  c0="$(sbatch --parsable --array=0-5%1 --job-name=wcf_c0 --output="${result_dir}/slurm/c0_%A_%a.out" "${script_dir}/slurm_train.sh" c0)"
else
  smoke="$(sbatch --parsable --output="${result_dir}/slurm/smoke_%j.out" "${script_dir}/slurm_smoke.sh")"
  c0="$(sbatch --parsable --dependency=afterok:${smoke} --array=0-5%1 --job-name=wcf_c0 --output="${result_dir}/slurm/c0_%A_%a.out" "${script_dir}/slurm_train.sh" c0)"
fi
sens="$(sbatch --parsable --dependency=afterok:${c0} --array=0-11%1 --job-name=wcf_sens --output="${result_dir}/slurm/sens_%A_%a.out" "${script_dir}/slurm_train.sh" sensitivity)"
gate="$(sbatch --parsable --dependency=afterok:${sens} --output="${result_dir}/slurm/gate_%j.out" "${script_dir}/slurm_gate.sh")"
eval_job="$(sbatch --parsable --dependency=afterok:${gate} --array=0-17%1 --output="${result_dir}/slurm/eval_%A_%a.out" "${script_dir}/slurm_evaluate.sh")"
final="$(sbatch --parsable --dependency=afterok:${eval_job} --export=ALL,WCF_SLURM_START="${start}" --output="${result_dir}/slurm/final_%j.out" "${script_dir}/slurm_final.sh")"
python - "${result_dir}/submitted_jobs.json" "${start}" "${smoke}" "${c0}" "${sens}" "${gate}" "${eval_job}" "${final}" <<'PY'
import json, pathlib, sys
path, start, smoke, c0, sens, gate, evaluate, final = sys.argv[1:]
pathlib.Path(path).write_text(json.dumps({
    "submitted_at": start,
    "jobs": {"smoke": smoke, "c0": c0, "sensitivity": sens,
             "gate": gate, "evaluate": evaluate, "final": final},
    "no_sdgd_jobs": True,
}, indent=2, sort_keys=True) + "\n")
PY
echo "smoke=${smoke} c0=${c0} sensitivity=${sens} gate=${gate} evaluate=${eval_job} final=${final}"
