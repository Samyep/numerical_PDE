# Literature positioning and fair-comparison boundary

This note separates **same-task numerical baselines** from **method-family
positioning**.  Published errors are not copied into the HCFL result table
when the governing equation, initial-condition distribution, grid, time
horizon, reference solver, or training objective differs.

## Same-task numerical comparisons used here

All methods in `run_64_benchmark.py` receive the same 64 finite-volume cell
averages.  A fine reference is created only by piecewise-constant
prolongation of those averages, fine-grid evolution, and conservative
restriction back to 64 cells.  The numerical table includes:

1. native first-order HLLC-64 (Euler) or HLL-64 (SWE);
2. TVD MUSCL-HLLC-64 or MUSCL-HLL-64;
3. the common four-cell HCFL model over independent training seeds; and
4. a compact residual FNO trained on exactly the same 64-cell trajectories;
   and
5. a four-cell unconstrained learned numerical flux trained on the same data.

The FNO is a controlled architecture baseline, not a line-by-line
reproduction of a published implementation.  It directly predicts the next
state and has no finite-volume flux form, positivity layer, or entropy
projection.

The learned-flux negative control retains the finite-volume flux-difference
update and exact equal-state consistency, but has no hard Tadmor projection,
positivity limiter, or fully-discrete entropy line search.  It makes the
comparison to the learned-discretization/neural-FV family numerical as well as
conceptual, while still being labelled as a matched family baseline rather
than a reproduction of one paper.

## Representative method families

| Family | Representative primary work | What is structurally enforced | Why published numbers are not inserted into our numerical ranking |
|---|---|---|---|
| PINN | [Raissi, Perdikaris & Karniadakis, Physics Informed Deep Learning](https://arxiv.org/abs/1711.10561) | The PDE residual enters the optimization objective. | A conventional PINN is normally optimized for a particular space-time problem and does not share our amortized local numerical-flux task. Shock solutions also make a strong-form residual comparison protocol-dependent. |
| Conservative PINN | [Jagtap, Kharazmi & Karniadakis, cPINN](https://doi.org/10.1016/j.cma.2020.113028) and its [authors' implementation](https://github.com/AmeyaJagtap/Conservative_PINNs) | Strong flux continuity across neural subdomain interfaces. | This is domain-decomposed solution fitting, not a reusable 64-cell time-step operator. Flux continuity is not the same statement as positivity or a Tadmor/fully-discrete entropy inequality. |
| Operator learning | [Li et al., Fourier Neural Operator](https://arxiv.org/abs/2010.08895) | A mesh-oriented mapping between function spaces; the original architecture does not by itself impose FV conservation, positivity, or entropy admissibility. | We therefore train a same-data compact FNO locally and report its own accuracy, completion, conservation drift, and entropy diagnostics instead of importing Burgers/Darcy/Navier–Stokes numbers. |
| Physics-encoded divergence-free operator | [Khorrami et al., PeFNO](https://arxiv.org/abs/2408.15408) | A stress-potential representation makes the predicted stress field divergence-free by construction. | This is the closest conceptual comparison to a hard output constraint, but `div(sigma)=0` for quasi-static solids is not the Euler/SWE entropy condition. A 1D compressible velocity is not supposed to be divergence-free, so a numerical head-to-head baseline would be physically mismatched. |
| Learned discretization | [Bar-Sinai et al., data-driven discretizations](https://doi.org/10.1073/pnas.1814058116) | Learned low-resolution spatial discretization inside a known PDE evolution. | It is highly relevant motivation for coarse-grid accuracy, but the paper studies different 1D equations and does not supply our system-level hard positivity and entropy guarantees. |
| Conservative learned flux | [Chen, Gelb & Lee, conservative form network](https://arxiv.org/abs/2211.14375) | Flux-difference form gives discrete conservation and improves shock-speed behavior while learning an unknown physical flux. | Their principal problem is discovery of an unknown governing flux. HCFL assumes the Euler/SWE physical flux is known and learns a stable coarse numerical flux/dissipation correction. |
| Neural finite volume | [Lichtlé et al., (U)NFV](https://arxiv.org/abs/2505.23702) | Conservation is built into an extended space-time FV architecture; supervised and weak-residual variants are proposed. | This is the closest current learned-FV family, but the reported study focuses on first-order scalar conservation laws. Its published error ratios cannot be transferred to Euler/SWE systems or interpreted as hard entropy/positivity guarantees. |

## Defensible HCFL claims from this experiment

Subject to the reported protocol and tolerances, the intended claims are:

- higher 64-cell rollout accuracy than specified native first-order baselines;
- a direct comparison with a stronger TVD MUSCL baseline, whether HCFL wins
  or loses;
- finite-volume conservation up to floating-point roundoff for periodic
  rollouts;
- positive density/pressure (Euler) or depth (SWE) for every accepted HCFL
  update;
- a hard interface Tadmor inequality after projection; and
- a nonincreasing fully-discrete total mathematical entropy for accepted
  periodic updates, because the final flux is line-searched toward an
  already-checked low-order endpoint.

The last two guarantees are conditional on the solver being able to establish
the admissible low-order premise by its adaptive CFL reduction.  They are
forward-solver properties, not claims that the raw neural proposal is always
feasible.

## Claims this experiment does **not** support

- universal accuracy superiority over PINNs, neural operators, or every
  learned finite-volume method;
- a divergence-free compressible Euler/SWE velocity field;
- TVD, monotonicity, or shock-oscillation-free behavior merely from entropy
  admissibility;
- generalization to bathymetric/well-balanced SWE; or
- a theorem covering all admissible initial states and boundary conditions.
