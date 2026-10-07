"""Summarise screen.jsonl: per-cell mean skill, raw generator bias, gap closed, equal-cost comparison."""
import json
import sys
from collections import defaultdict

import numpy as np

rows = [json.loads(l) for l in open(sys.argv[1] if len(sys.argv) > 1 else "../../results/candidate_screen/screen.jsonl")]
cell = defaultdict(lambda: defaultdict(list))
meta = {}
for r in rows:
    k = (r["label"], r["d"], r["n"], r["M"])
    cell[k][r["method"]].append(r["skill"])
    if r["method"] == "raw":
        cell[k]["_rawbias"].append(r["generator_bias"] if r["generator_bias"] is not None else np.nan)
    cell[k]["_f"].append(r["f_calls"] or 0)
    if "predicted_bias_factor" in r:
        meta[(r["label"], r["d"])] = r["predicted_bias_factor"]

m = lambda v: float(np.mean(v)) if len(v) else float("nan")
lines = ["| candidate | d | (n,M) | raw | segment | box | oracle | centre | f=0 | raw gen. bias | Gc(best cert.) |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
eq_cost = []
for label in dict.fromkeys(k[0] for k in cell):
    for d in sorted({k[1] for k in cell if k[0] == label}):
        best = {}
        for (n, M) in sorted({(k[2], k[3]) for k in cell if k[0] == label and k[1] == d}, key=lambda x: (x[0], x[1])):
            c = cell[(label, d, n, M)]
            sk = {meth: m(c[meth]) for meth in ("raw", "segment", "box", "oracle_state", "centre", "f_zero")}
            cert = min(sk["segment"], sk["box"])
            gc = (sk["raw"] - cert) / (sk["raw"] - sk["oracle_state"]) if sk["raw"] > sk["oracle_state"] else float("nan")
            lines.append(f"| {label} | {d} | ({n},{M}) | {sk['raw']:.3g} | {sk['segment']:.3f} | {sk['box']:.3f} | "
                         f"{sk['oracle_state']:.3f} | {sk['centre']:.3f} | {sk['f_zero']:.3f} | {m(c['_rawbias']):.3g} | {gc:.2f} |")
            f = m(c["_f"])
            for meth in ("raw", "segment", "box", "centre"):
                best.setdefault(meth, []).append((sk[meth], f, (n, M)))
        rb = min(best["raw"]); cb = min(best["segment"] + best["box"]); ce = min(best["centre"])
        cheaper = [v for v in best["segment"] + best["box"] if v[1] <= rb[1] * 1.0001]
        cc = min(cheaper) if cheaper else (float("nan"), 0, None)
        eq_cost.append(f"| {label} | {d} | {rb[0]:.3f} @{rb[2]} | {cc[0]:.3f} @{cc[2]} | {cb[0]:.3f} @{cb[2]} | {ce[0]:.3f} @{ce[2]} | "
                       f"{meta.get((label, d), float('nan')):.1f} |")
print("\n".join(lines))
print("\n| candidate | d | best raw | certified at <= raw cost | best certified | best centre | predicted C2 bias factor |")
print("|---|---|---|---|---|---|---|")
print("\n".join(eq_cost))
