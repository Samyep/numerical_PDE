"""Summarise defect_probe.jsonl: mean skill per (k, d, cell, method)."""
import json
import sys
from collections import defaultdict

import numpy as np

rows = [json.loads(l) for l in open(sys.argv[1])]
cell = defaultdict(lambda: defaultdict(list))
sur = defaultdict(dict)
for r in rows:
    cell[(r["k"], r["d"], r["n"], r["M"])][r["method"]].append(r["skill"])
    if "surrogate_skill" in r:
        sur[(r["k"], r["d"])][r["method"].split("_e")[-1]] = r["surrogate_skill"]
cols = ["raw", "centre", "oracle_state", "path", "path_cv_0", "path_cv_gradg"]
eps = ["0.03", "0.1", "0.3"]
dcols = [f"defect_{t}_e{e}" for e in eps for t in ("bismut", "path")]
print("| k | d | (n,M) | " + " | ".join(cols + [c.replace("defect_", "") for c in dcols]) + " |")
print("|---|---|---|" + "---|" * (len(cols) + len(dcols)))
f = lambda v: f"{np.mean(v):.3g}" if v else "-"
for key in sorted(cell):
    c = cell[key]
    print(f"| {key[0]} | {key[1]} | ({key[2]},{key[3]}) | " + " | ".join(f(c[m]) for m in cols + dcols) + " |")
print("\nsurrogate-alone skill (eps -> skill):")
for kd, v in sorted(sur.items()):
    print(kd, {e: round(v[e], 3) for e in eps if e in v})
