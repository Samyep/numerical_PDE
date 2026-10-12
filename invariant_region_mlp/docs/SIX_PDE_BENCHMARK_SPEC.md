# Six PDE Benchmark Specification — Bias-Stabilized Full-History MLP

Status: **research design / proposed benchmark suite**; not a claim that all six have already been implemented or solved. This file is a specification, not experimental results.

Prepared 2026-10-11. Source baseline: branch `ir-mlp-baseline-gpu`; read-only existing implementations under `experiments/mechanism_suite/`, `experiments/candidate_screen/`, `experiments/double_estimator/`, and `experiments/gpu_mlp/`. The separately running non-ridge and actor jobs must not be modified.

## Research question and scope

**Target:** Repair *full-history MLP's own* nonlinear noise-to-bias amplification for gradient-dependent semilinear PDEs. The primary main-table comparison is **unmodified full-history MLP vs one complete repaired MLP**. Pathwise terminal estimator, level-matched source baseline, double estimator and learned actor appear primarily in mechanism/ablation experiments. Actor is a candidate principal replacement for double, **not yet empirically validated**.

For the semilinear cases, use the existing convention

\[
u_t+\Delta_x u+f(u,z)=0,\qquad z=\sqrt{2}\nabla_x u,\quad
u(T,x)=g(x),\quad x\in\mathbb R^d.
\]

Experiment axes: ambient dimension \(d\in\{100,400\}\), effective dimension \(k\), Picard depth \(n\), sampling parameter \(M\), and *the same* total-work accounting as current experiments (generator calls, Gaussian-coordinate work, and wall time, including actor training when used). Do **not** infer full-rank difficulty solely from a syntactically full-rank terminal. Before comparing solvers report the heat-only \(f=0\) discrepancy, spatial variability of \(f(u^*,z^*)\), gradient covariance effective rank, and Hessian effective rank.

A \(k\)-ridge solution means \(u(t,x)=v(t,B^\top x)\) for fixed \(B\in\mathbb R^{d\times k}\), \(B^\top B=I_k\). Multi-ridge with \(k=10\) or 20 is **still ridge**, but not one-dimensional; compare it alongside genuinely full-rank tests. Full rank need not imply a hard benchmark: high-dimensional concentration can make the nonlinear correction nearly spatially constant.

## Benchmark overview

| ID | Equation / test | Geometry | Ground truth | Present status |
| --- | --- | --- | --- | --- |
| B1 | Coupled multi-ridge quadratic HJ, \(f(z)=-\|z\|^2/2\) | \(k=1,5,10,20,40\), \(d=100,400\) | **Analytic** Hopf–Cole exponential sum, value and gradient | Extends P1/MR; new coupled \(k>1\) inputs needed |
| B2 | Negative Euclidean norm, \(f(z)=-\lambda\|z\|/\sqrt2\) | Radial \(k=10,20\), ambient \(d=100,400\) | Independently refined **1D radial** PDE, not yet produced | Existing P4 has only \(k=1\) reference |
| B3 | Smooth cubic, \(f(z)=-\kappa(\|z\|/\sqrt2)^3\) | Radial \(k=10,20\), \(d=100,400\) | Independently refined **1D radial** PDE, not yet produced | Existing C3 has only \(k=1\) reference |
| B4 | **Fully nonlinear** controlled-volatility HJB | Hessian-dependent, \(k=2\)–5 active coordinates | Independently refined **1D component** PDEs, not yet produced | **Outside** current \((u,z)\)-only MLP; conditional extension |
| B5 | Two-actor zero-sum Isaacs LQ game | Coupled \(k=10,20\) subspace with two control blocks | **Analytic** finite-horizon matrix Riccati | Existing C2 uses \(k=1\); new coupled game needed |
| B6 | Full-rank, nonseparable quadratic HJ | \(k=d=100,400\) | **Analytic** Hopf–Cole exponential sum | New terminal/gradient and cost calibration needed |

**Scope warning:** B1 and B6 have the *same differential operator* and different terminal geometries; these are six distinct **benchmark instances**, not six mathematically distinct PDE operators. B4 does **not** support a claim that existing gradient-only full-history MLP already solves general fully nonlinear HJBs.

## B1 — Coupled multi-ridge quadratic HJ (principal controlled-\(k\) benchmark)

Let \(B\in\mathbb R^{d\times k}\), \(B^\top B=I\), \(y=B^\top x\). Choose \(K\ge k+1\) affine-spanning vectors \(b_m\in\mathbb R^k\) (for \(k\ge2\), initial proposal \(K=2k+1\)). Define

\[
\begin{split}
u_t+\Delta u-\tfrac12\|z\|^2&=0,\quad z=\sqrt2\nabla u,\\
g_d(x)&=-\log\left[\frac1K\sum_{m=1}^{K}
           e^{b_m^\top B^\top x}\right],\qquad u(T,x)=g_d(x).
\end{split}
\]

**Exact solution, for any \(d,k,K\):** with \(\tau=T-t\),

\[
u(t,x)=-\log\left[\frac1K\sum_m
   \exp\{b_m^\top B^\top x+\tau\|b_m\|^2\}\right].
\]

Define \(\pi_m\propto\exp\{b_m^\top B^\top x+\tau\|b_m\|^2\}\). Then

\[
z^*(t,x)=-\sqrt2 B\sum_m\pi_m b_m,\quad
a^*(t,x)=-z^*(t,x).
\]

The Hessian is
\[
\nabla_x^2u=-B\,\mathrm{Cov}_{\pi}(b_m)\,B^\top.
\]
If \(b_m\) affinely span \(\mathbb R^k\), this has rank \(k\) at every finite \((t,x)\); nevertheless effective rank on a test distribution can be much lower. Unless the vectors admit a special decomposition, the effective coordinates are **coupled**, not additive separable.

**Link to existing experiments:** \(k=1\), \(K=3\), \(b=(0.5,1.5,3)\) reproduces the existing P1 exponential profile (given a matching direction \(B\)). The existing MR is also \(k=10\) but uses a different randomly generated terminal and time horizon; keep it as a legacy mechanism benchmark rather than silently replacing it. The new B1 permits controlled \(k=1,5,10,20,40\).

**Identity check:** \(w=e^{-u}\) satisfies \(w_t+\Delta w=0\), terminal \(w(T,x)=e^{-g_d(x)}\); Gaussian moment generating functions give the analytic sum above. Validate PDE residual, terminal equality, analytic \(z\), and exact formula using an independent heat-expectation MC test.

**Calibration before looking at method results:** seed and scale of \(b_m\), training/test distribution, and \(T\) must be frozen *after checking reference-only properties* (nontrivial heat-only correction; nontrivial spatial variation of \(f(z^*)\); gradient-covariance effective rank not collapsing to 1; no runaway exponential weights). Parameter changes motivated by MLP outcomes would violate the intended design.

## B2 — Nonsmooth concave norm driver with truly \(k\)-dimensional gradients

For \(B^\top B=I_k\), \(y=B^\top x\), \(r=\|y\|_2\), choose \(k=10,20\) and \(d=100,400\):

\[
u_t+\Delta u-\frac{\lambda}{\sqrt2}\|z\|_2=0,\quad z=\sqrt2\nabla u,
\qquad
g(x)=\frac{A}{\beta}\log\cosh(\beta r).
\]

Starting parameters may inherit P4 (\(\lambda=1\), \(A=1\), \(\beta=2\), \(T=0.25\)), **subject to independent reference calibration**. The exact reduced PDE, writing \(u(t,x)=v(t,r)\), is

\[
v_t+v_{rr}+\frac{k-1}{r}v_r-\lambda|v_r|=0,\qquad
v(T,r)=\frac A\beta\log\cosh(\beta r),\quad v_r(t,0)=0.
\]

At \(r=0\), the regular radial Laplacian is interpreted by its smooth limit \(k\,v_{rr}(t,0)\). Produce a converged monotone radial finite-volume or finite-difference solution with spatial/time/boundary refinement; **do not reuse the \(k=1\) P4 reference**. Although reduced to a 1D radial equation, at general points \(\nabla u\) points in varying directions inside \(\mathrm{span}(B)\) and has \(k\)-rank gradient covariance for an appropriate isotropic nondegenerate sampling distribution.

The Hamiltonian is nonsmooth at \(z=0\); a fixed actor may choose \(a_\theta\in\{a:\|a\|\le\lambda/\sqrt2\}\), with
\[
-\lambda\|z\|/\sqrt2=\min_{\|a\|\le\lambda/\sqrt2} a^\top z.
\]
Do **not** claim a global quadratic actor-regret bound at \(z=0\).

## B3 — Smooth, nonquadratic cubic driver with radial \(k\)-dimensional gradients

Reuse B2's \(B,y,r\), but choose

\[
u_t+\Delta u-\kappa\left(\frac{\|z\|_2}{\sqrt2}\right)^3=0,
\qquad
g(x)=\frac A\beta\log\cosh(\beta r).
\]

Initial parameters from C3: \(\kappa=0.25,A=2,\beta=2,T=0.25\); **not yet validated for \(k=10,20\)**. The reduced radial equation is

\[
v_t+v_{rr}+\frac{k-1}{r}v_r-\kappa|v_r|^3=0,\qquad v_r(t,0)=0.
\]

Independently refine a stable radial reference. This is concave and differentiable in \(z\), but its Hessian is not globally bounded on all of \(\mathbb R^d\): the second-order tangent-error bound requires a certified or experimentally controlled region for \(\|z\|\) and actor predictions. The actor slope is

\[
\nabla_z f(q)=-\frac{3\kappa}{2\sqrt2}\|q\|q.
\]

The existing C3 implementation/reference treats \(k=1\); keep it as a validation special case, not the full \(k>1\) ground truth.

## B4 — Genuinely fully nonlinear controlled-volatility HJB (conditional stretch benchmark)

This is **not** a routine input to the current \((u,z)\) semilinear MLP. Take \(k\in\{2,3,5\}\) active physical coordinates, \(\eta>0\), and

\[
u_t+\Delta u+\eta\sum_{j=1}^k\min(0,u_{x_jx_j})-\tfrac12\|z\|^2=0,
\quad z=\sqrt2\nabla u,
\quad
g(x)=\frac1{\sqrt k}\sum_{j=1}^k\gamma_j\sin(x_j+\varphi_j).
\]

The Hessian term equals \(\min_{a_j\in[0,1]}\eta\sum_j a_j u_{x_jx_j}\), corresponding to controlled per-coordinate diffusion variance; the quadratic gradient term has a separate drift-control representation. This operator is fully nonlinear in the highest-order derivative.

With the displayed *additive* terminal and coordinate-aligned controlled volatility, the solution decomposes into independently solvable **one-dimensional fully nonlinear** component equations; that makes a high-accuracy, independently refined reference feasible. This is intentionally a reference-friendly structural benchmark, not evidence for general nonseparable fully nonlinear HJB.

**Gate:** do not claim B4 solved by the existing algorithm until a mathematically valid Hessian/controlled-diffusion extension, independent reference verification, and a fair stronger baseline exist. If this gate is not met, B4 is a prospective extension, not a main-table success.

## B5 — Two-actor, zero-sum coupled Isaacs LQ game

Split ambient coordinates into control blocks \(A,B\) (rename the *basis* below \(\mathsf B\) to avoid confusion). Choose orthonormal \(\mathsf B=[\mathsf B_A,\mathsf B_B]\in\mathbb R^{d\times k}\), with the basis columns supported within their respective ambient coordinate blocks. Choose a **dense, positive-definite** \(Q\in\mathbb R^{k\times k}\), including cross-block off-diagonal terms. Let \(y=\mathsf B^\top x\), \(k=k_A+k_B=10,20\).

\[
u_t+\Delta u-\frac a2\|z_A\|^2+\frac b2\|z_B\|^2=0,\quad
z=\sqrt2\nabla u,\quad
g(x)=\tfrac12 y^\top Qy,
\]
where \(a,b>0\) and \(z_A,z_B\) denote the *full ambient gradient coordinate blocks*, so the generator still responds to off-subspace BEL noise. Under the chosen block-supported \(\mathsf B\), the true gradient is contained in \(\mathsf B\), with a nontrivially coupled \(k\)-dimensional reduced game.

Set \(D_k=\operatorname{diag}(-a I_{k_A},b I_{k_B})\). As long as the inverses below exist throughout the horizon, the exact solution is
\[
\begin{aligned}
P(t)&=[Q^{-1}-2(T-t)D_k]^{-1},\\
u(t,x)&=\tfrac12 y^\top P(t)y+c(t),\\
c(t)&=\int_t^T\operatorname{tr}P(s)\,ds,\\
z^*(t,x)&=\sqrt2\,\mathsf B P(t)y.
\end{aligned}
\]
This follows from \(P'=-2P D_kP\), \(c'=-\operatorname{tr}P\). Evaluate the scalar integral accurately (or via a validated matrix log-determinant expression). Check invertibility and PDE residual across all times/test points.

The min/max control representation is
\[
\min_{p}\max_{q}\left[
p^\top z_A+\frac{\|p\|^2}{2a}
+q^\top z_B-\frac{\|q\|^2}{2b}
\right],
\quad
p^*=-a z_A,\ \ q^*=b z_B.
\]
Compare both learned actors and note the driver is neither convex nor concave in the *combined* gradient. The legacy C2-cancel/C2-flip cases are \(k=1\) ridge variants useful as ablations only, **not** the full-coupling B5.

## B6 — Full-rank nonseparable quadratic HJ (structural stress test)

Same equation and actor as B1, but set \(k=d\), \(\mathsf B=I_d\), and \(K\ge d+1\) dense vectors \(b_m\in\mathbb R^d\) which affinely span \(\mathbb R^d\):
\[
\begin{aligned}
u_t+\Delta u-\tfrac12\|z\|^2&=0,\\
g(x)&=-\log\left[\frac1K\sum_m e^{b_m^\top x}\right],\\
u^*(t,x)&=-\log\left[\frac1K\sum_m
e^{b_m^\top x+(T-t)\|b_m\|^2}\right],\\
z^*(t,x)&=-\sqrt2\sum_m \pi_m(t,x)b_m.
\end{aligned}
\]

If \(\{b_m\}\) affinely span \(\mathbb R^d\), the Hessian \(-\mathrm{Cov}_\pi(b_m)\) has rank \(d\) at every finite input. A generic dense set has no fixed small ridge subspace and no independent coordinate decomposition. Yet an exponent sum is still a deliberately manufactured *analytically soluble* family: do not equate it with an arbitrary fully coupled nonlinear HJB.

**Compute caveat:** \(K=d+1\) already makes each naive dense evaluation \(O(d^2)\); larger \(K\) can be prohibitive inside recursive MLP. Include the terminal-evaluation cost in wall time, and calibrate the \(b_m\) scale (softmax concentration, gradient effective rank, driver spatial variability, heat-only deviation) *before* method evaluation. Do not assume full-rank automatically implies an appreciable nonlinear correction.

## Cross-benchmark checks and evaluation rules

1. **Before algorithms:** verify terminal match, PDE residual, analytic or independently refined reference, and gradient via independent finite differences. Log all seeds, \(d,k\), basis generation, terminal coefficients and time horizon. For B2/B3 radial references check convergence near \(r=0\) and at the outer boundary; report reference uncertainty separately from solver error.
2. **Avoid the easy-benchmark trap:** report \(\|u^*-u_{\mathrm{heat}}\|\), spatial standard deviation of \(f(z^*)\), \(\|\nabla u^*\|\), gradient-covariance effective rank, Hessian effective rank, and any exact symmetries. Examine whether the driver is nearly constant in space.
3. **Primary comparison:** vanilla full-history MLP versus **one complete corrected MLP** under matched \(n,M,d,k\), repeated seeds and matched wall time/operation counts. Do not fill the main benchmark table with separate path/baseline/double/actor variants. Explain actor training overhead; do not compare a trained method on inference-only cost against an untrained baseline.
4. **Ablations:** remove pathwise, remove level-matched baseline, remove actor, and compare actor versus double *when the Hamiltonian admits the double estimator*. Double is not a universal default for B2/B3. Report directly measured value bias and variance over independent repetitions, not just a slope-based inference.
5. **Targets are hypotheses:** \(n=3\) beating the finite-depth \(n=2\) accuracy floor; strong performance as \(k\) rises; dimension robustness when \(d-k\) is large; and actor's additional benefit beyond already strong path+baseline. No result is claimed ahead of experiments.
6. **Qualification:** B1/B2/B3/B5/B6 are semilinear in the highest-order term and are the primary scope. B4 is explicitly conditional pending a new fully nonlinear algorithm and should never be used as evidence for an unimplemented feature.

## Immediate implementation plan

- **First:** B1 coupled \(k\)-ridge with analytic formula; use the P1 code shape, extending vectors to a \(K\times k\) coefficient matrix and \(d\times k\) orthonormal basis. Verify \(k=1\) exactly reproduces P1 on matched inputs. Do *not* overwrite current P1/MR results.
- **Second:** B5 analytic coupled game (Riccati) and B6 full-rank nonseparable LSE, with nonlinearity-strength calibration.
- **Third:** B2/B3 stable radial reference builders for \(k=10,20\), then add new GPU PDE families without modifying frozen baselines.
- **Last / conditional:** B4 controlled-volatility / Hessian extension; do not make the current semilinear paper's correctness depend on this stretch benchmark.

Related in-repository code: `experiments/mechanism_suite/equations.py` (P1), `experiments/mechanism_suite/norm_hjb.py` (P4), `experiments/candidate_screen/candidates.py` (C2/C3), `experiments/double_estimator/equations.py` (MR/game registry), `experiments/gpu_mlp/solver.py` (GPU).
