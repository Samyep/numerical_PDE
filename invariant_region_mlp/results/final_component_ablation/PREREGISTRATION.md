# Final Tangent-Actor Component Ablation: Preregistration

This document freezes the experiment before any new formal GPU result is
generated on branch `ir-mlp-final-component-ablation`.

## Primary matrix

The study contains exactly four 100-dimensional PDEs (`hjb_log`,
`joint_concave`, `joint_convex`, `nonconvex_coupled`), exactly four methods
(`a0_full`, `a1_no_mcv`, `a2_no_terminal_anchor`, `a3_no_updates`), and
training seeds 0, 1, and 2. There is no hyperparameter search.

All learned methods use 8,192 training states, 512 validation states, 128
paths per label, 32 Euler label steps, 5,000 optimizer steps per update,
width 384, depth 4, and four tangent updates. The HJB benchmark preserves its
published `mixed` training-state distribution; the joint benchmarks preserve
their audited source distributions. Learning rate, optimizer, activation,
loss normalization, clipping, batching, warm starts, and noisy-label
checkpoint rule are unchanged.

The primary whole-function evaluation uses 960 frozen states, 64 Euler steps,
and 1,024 paths. Replicate 0 is primary. Replicate 1 is an independent
inference repeat and remains separate. RMSE is primary; RMAE is the original
ratio of sums.

## Ablations

- `a0_full`: MCV labels, MCV fixed-point selection, MCV inference, hard
  terminal anchoring, and learned updates.
- `a1_no_mcv`: plain labels, plain fixed-point selection, and plain inference.
- `a2_no_terminal_anchor`: retains the exact terminal-data offset and zero
  residual initialization but removes the `(T-t)` multiplier; all value range
  constraints remain.
- `a3_no_updates`: evaluates the deterministic terminal extension `a0` using
  MCV and does no fitting.

Full checkpoints are reused only after exact source/PDE/configuration/seed/
selection verification. Reused actors are re-evaluated whenever the frozen
evaluation protocol does not match exactly.

## Frozen 960-state nonconvex clarification

The historical Oscillatory Value–Gradient artifact has only 512 states. The
task mandates 960 states for every PDE. Therefore its old RMSE is not a
primary result here. Before formal evaluation, this study freezes a 960-state
bank using the already audited ICML-final-suite rule
`stable_seed("icml_final_test_bank", "nonconvex_coupled", 100)`, with
`t ~ Uniform[0,0.5)` and `x ~ N(0,I_100)`. The existing independently refined
2D periodic reference supplies truth on those states. All four variants use
this identical bank.

## Scientific gates

Formal runs require tests for PDE identity, convex sign conversion, matching
initial actors, terminal identity and its deliberate removal, fixed-actor
plain/MCV expectation agreement, pathwise gradients, retained value bounds,
bank/reference hashes, and absence of test-reference use before selection.
Failures are reported without variant-specific stabilization. The allocated
GPU-hour ceiling is 36 hours on a single consistent GH200 model.

Ground truth is forbidden in training, validation, selection, and tuning.
Test RMSE is computed only after each run's checkpoint and reference-free
selection record have been frozen.
