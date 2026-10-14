# Frozen external-baseline protocols

Frozen before GPU results: 2026-10-14 (America/New_York).

This study uses exactly one preset configuration for each external baseline on both PDEs. This follows the user's final instruction that external baselines should not be tuned and should be run once where possible. Checkpoints are retained to form compute curves, but analytic Value or gradient references never choose an external checkpoint. The machine-readable source of truth is `experiments/hjb_compute_frontier/protocol_frozen.json`.

## PINN + SDGD

The official [SDGD HJB-Log implementation](https://github.com/zheyuanhu01/SDGD_PINN/blob/main/Section_5.3/HJB_Log.py) uses a hard terminal ansatz, a 4-layer width-1024 Tanh network, 10,000 Adam epochs, learning-rate decay, coordinate-sampled Laplacians, and adversarial collocation.

Our frozen scaled setting uses a 4-layer width-512 Tanh correction network with the exact terminal condition hard-enforced, batch size 128, two independent sets of 32 sampled coordinates for the unbiased squared-residual gradient, 10,000 Adam steps, initial learning rate `1e-3`, and decay by `0.2` at steps 5,000 and 8,000. Checkpoints are saved at 2,000, 5,000, and 10,000 steps. Eight frozen validation points use all 400 diagonal Hessian entries, so checkpoint selection minimizes an exact full-Laplacian residual MSE rather than a noisy stochastic residual product. The full preset trajectory runs even if an intermediate residual stalls.

Width 512 halves the official hidden width while retaining its depth and full optimizer horizon. We omit a separate adversarial-collocation candidate because the user's final instruction freezes one baseline configuration rather than tuning variants; this limitation will be reported with the results.

## Amortized Deep BSDE

The official [DeepBSDE HJB configuration](https://github.com/frankhan91/DeepBSDE/blob/master/configs/hjb_lq_d100.json) is a single-start solver with dimension 100, 20 time intervals, two width-110 hidden layers, batch size 64, learning rate `1e-2`, 2,000 iterations, and float64 arithmetic.

Our task requires predictions for arbitrary (x_0\sim N(0,I_{400})), so the frozen adaptation learns an amortized scalar (u_0(x_0)) network and a shared time-dependent (z(t,x)) network. Both have three width-256 SiLU layers. It uses 32 Euler steps, batch size 128, 10,000 Adam steps, initial learning rate `1e-3`, decay by `0.2` at 5,000 and 8,000, terminal Huber training loss, and gradient clipping. Checkpoints are saved at 2,000, 5,000, and 10,000 steps. A fixed independent bank of 256 initial states and Brownian paths selects the checkpoint with the lowest terminal mismatch MSE. This is explicitly an amortized adaptation, not a claim to reproduce the official single-point architecture exactly.

## Deep Picard Iteration

The representative [official DPI HJB configuration](https://github.com/DeepOptimalControl/DeepPicardIteration/blob/main/scripts/hjb/base_100d_T1.0_w0.1_0.yaml) uses 40 Picard rounds, 4,096 training points, 4,096 terminal and integral MC samples, 16 epochs, batch size 512, a 4-layer width-512 ELU model, gradient supervision with weight 0.1, warm starts, and a specialized `PISGRADNET`.

Our frozen scaled setting retains 4,096 training points per round, 16 epochs, batch size 512, a 4-layer width-512 ELU network, value plus (z) supervision with weight 0.1, Adam at `1e-3`, and warm starts. It uses 20 preset Picard rounds and 512 MC paths per value/gradient label. A joint value/(z) network blends both outputs with the analytic terminal value and terminal derivative, providing the terminal-aware behavior of the published specialized model without using any nonterminal reference. Checkpoints are saved after rounds 1, 2, 5, 10, and 20. Because different rounds regress different Picard targets, their losses are not treated as a comparable stopping sequence; the preset 20-round ceiling selects the final model, while within-round held-out target loss audits fitting quality.

## Fairness and selection

- B1 uses (d=400,k=40,T=0.25); HJB-Log uses (d=400,T=1).
- Every method receives the same terminal condition and diffusion convention (dX=\sqrt2,dW).
- External baselines cannot access `exact_u` or `exact_z` during training or checkpoint selection. Static leakage checks enforce this at preflight.
- All methods use the same new 2,048-point final test bank per PDE. The historical 960 B1 points are excluded.
- External methods expose their whole preset checkpoint trajectory. No point on their curve is selected by test error.
- The Actor method follows the user's later override and may choose among declared Actor candidates by test RMSE. Results label this `test_rmse_oracle_allowed_by_user`; search and selection cost are disclosed separately.
