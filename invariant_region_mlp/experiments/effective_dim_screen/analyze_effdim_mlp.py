"""Summarise effdim_mlp.jsonl: per-cell mean skill and raw generator bias; equal-cost comparison per (k, d)."""
import json
import sys
from collections import defaultdict

import numpy as np

rows = [json.loads(l) for l in open(sys.argv[1])]
cell = defaultdict(lambda: defaultdict(list))
for r in rows:
    cell[(r["k"], r["d"], r["n"], r["M"])][r["method"]].append(r)
meths = ["raw", "sub_box", "subspace", "box", "ball", "oracle_state", "centre", "f_zero"]
mean = lambda v, key="skill": float(np.mean([x[key] for x in v])) if v else float("nan")
print("| k | d | (n,M) | " + " | ".join(meths) + " | raw gen. bias | sub_box gen. bias |")
print("|---|---|---|" + "---|" * (len(meths) + 2))
for key in sorted(cell):
    c = cell[key]
    vals = " | ".join(f"{mean(c[m]):.3g}" for m in meths)
    print(f"| {key[0]} | {key[1]} | ({key[2]},{key[3]}) | {vals} | {mean(c['raw'], 'generator_bias'):.3g} | "
          f"{mean(c['sub_box'], 'generator_bias'):.3g} |")
print("\nEqual cost (generator calls): best raw over cells vs best certified at <= that cost")
for k, d in sorted({(kk[0], kk[1]) for kk in cell}):
    pts = {m: [(mean(cell[kk][m]), mean(cell[kk][m], "f_calls"), kk[2:]) for kk in cell
               if kk[:2] == (k, d) and cell[kk][m]] for m in meths}
    rb = min(pts["raw"])
    best = {m: min([p for p in pts[m] if p[1] <= rb[1] * 1.0001], default=(float("nan"), 0, None))
            for m in ("sub_box", "subspace", "box", "ball", "centre")}
    print(f"k={k} d={d}: raw {rb[0]:.3f}@{rb[2]} | " + " | ".join(f"{m} {v[0]:.3f}@{v[2]}" for m, v in best.items()))
