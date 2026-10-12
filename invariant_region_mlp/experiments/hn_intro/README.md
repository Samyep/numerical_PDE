# HN (2025) introduction figure

Independent CPU NumPy implementation of Hutzenthaler--Nguyen arXiv:2506.23969 Eq. (6) with the prescribed Beta(1/2,1) random-time rule and centered BEL terminal estimator. The pathwise comparison changes only the terminal gradient. No actor/network is used.

From repo root (with Python dependencies installed):

```bash
python invariant_region_mlp/experiments/hn_intro/run_hn_intro.py \
  --dims 50 100 200 400 --n 3 --m 3 --reps 64 \
  --outdir invariant_region_mlp/results/hn_intro
```

Output `results.csv`, `summary.csv`, figure PDF/PNG, and report in the selected output directory. The archived report is `invariant_region_mlp/docs/HN_INTRO_REPORT.md`.

Exact pathwise theoretical RMSE for root terminal N=27 samples: value $\sqrt{(1-2/\pi)/27}$, full state $\sqrt{(2-2/\pi)/27}$. Confirm the code's built-in assertions pass.

**Caveat:** small finite-depth numerical test does not prove HN's asymptotic complexity theorem and does not establish neural actor necessity.
