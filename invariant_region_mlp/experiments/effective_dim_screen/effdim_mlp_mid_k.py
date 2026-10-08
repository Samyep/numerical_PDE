"""Fill in intermediate effective dimensions k in {20, 30} for the MLP screen (same code as effdim_mlp.py).
d = 100: cells (2,32), (3,6), (4,6); d = 400: cell (2,32) only. Settings are the screened passing
configurations (k=20: s=4, T=0.1; k=30: s=5, T=0.1). 240 points, 2 paired repetitions.
Usage: python effdim_mlp_mid_k.py OUT.jsonl [workers]
"""
import json, sys, time
from concurrent.futures import ProcessPoolExecutor
import effdim_mlp as E

SETTINGS = [dict(k=20, scale=4.0, T=0.1), dict(k=30, scale=5.0, T=0.1)]
jobs = [(st, 100, n, M, m, r) for st in SETTINGS for (n, M) in E.CELLS for m in E.METHODS for r in range(2)]
jobs += [(st, 400, 2, 32, m, r) for st in SETTINGS for m in E.METHODS for r in range(2)]
if __name__ == "__main__":
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=int(sys.argv[2]) if len(sys.argv) > 2 else 2) as pool, open(sys.argv[1], "a") as fh:
        for row in pool.map(E.task, jobs, chunksize=1):
            fh.write(json.dumps(row) + "\n"); fh.flush()
    print(f"done {len(jobs)} jobs in {time.time() - t0:.0f}s", flush=True)
