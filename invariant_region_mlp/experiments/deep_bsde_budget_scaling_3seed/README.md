# Deep BSDE budget scaling, three seeds

This isolated study trains only native scalar-initial-value Deep BSDE models.
Each final whole-function query is a distinct model initialized at its frozen
`(t_i,x_i)` and integrated over `T-t_i`.  It does not train or modify
Tangent-Actor, DPI, PINN-HTE, SDGD-PINN, or FH-MLP.

Workflow:

1. `python -m ...prepare` freezes input-only time-stratified subsets.
2. `python -m ...validate_math` and `pytest .../test_correctness.py` audit the adaptation.
3. Submit `slurm/preflight.slurm`; formal submission is blocked unless its 120-H100-hour gate passes.
4. Submit `slurm/submit_formal.sh`; the 42 array tasks are resumable bundles of 64 queries.
5. Run `finalize_training`, then aggregate only after the truth gate opens.

