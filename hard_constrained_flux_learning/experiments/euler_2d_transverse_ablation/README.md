# 2-D Euler transverse-stencil ablation

This experiment tests the diagnosed limitation of the original 2-D HCFL
model while keeping its conservative face update, signed Roe correction,
hard interface-entropy projection, proposal-feasibility loss, and deployment
safety wrapper unchanged.

Each face receives an oriented `3 x 6 = 18` cell patch.  Two nearly
parameter-matched models are compared:

- `flat18`: one MLP sees all cells and must learn whether transverse context
  should be ignored or used;
- `gated18`: a normal six-cell branch is augmented by a transverse residual
  that is exactly zero when the three rows coincide.

An additional `normal6wide` capacity control uses only the original six
normal cells but has exactly the same 10,804 parameters as `flat18`.  It was
added after `flat18` validation training began, but before any new test or
official-quadrant evaluation, specifically to challenge the transverse-input
attribution.

The original official quadrant result has already been observed.  It is
therefore only a post-hoc diagnostic here, never a model-selection set.  This
controlled ablation reuses the original train/validation split to isolate the
architecture change; a later confirmatory claim needs new held-out data.

Run:

```text
python test_transverse.py
python run_ablation.py train --variant normal6wide --seed 0
python run_ablation.py train --variant flat18 --seed 0
python run_ablation.py train --variant gated18 --seed 0
python run_ablation.py train --variant flat18 --seed 1
python run_ablation.py train --variant flat18 --seed 2
python run_ablation.py evaluate --variants normal6wide flat18 gated18 --seeds 0 1 2
python audit_results.py
```

After architecture selection, evaluate all three `flat18` seeds separately;
the seed-zero ablation controls remain single-seed diagnostic comparisons.
