# Invariant-Region MLP (IR-MLP)

Certified projection inside multilevel Picard recursions for high-dimensional semilinear PDEs.

This directory contains the current IR-MLP research prototype. The complete source snapshot, including all experiment scripts and JSON results, is stored in `archive/invariant-region-mlp-source.tar.gz`.

## Current results

- SCaSML-style HJB, n=2, M=10, 1,200 points x 10 repetitions: heuristic rel-L2 1.53-1.57 versus 0.79-0.85 with the exact-matrix joint projection.
- A non-oracle uniform HJB bound derived only from the public coefficient ranges gives rel-L2 0.83/0.85/0.86/0.88 in 100/120/140/160D.
- 100D nonlinear funding: at n=4,M=3, baseline MAE 3.025, joint hard projection 0.877, coordinate box 2.708.
- 100D counterparty-credit-risk negative control: the certified value interval is almost never violated and hard projection is effectively identical to baseline.

See `docs/current_round_report.md` for the compact tables and `integration/scasml/` for the upstream patch design.
