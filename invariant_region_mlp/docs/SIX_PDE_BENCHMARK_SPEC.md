# Six PDE Benchmark Specification — Bias-Stabilized Full-History MLP

Status: **research design / proposed benchmark suite**; not a claim that all six have already been implemented or solved. This file is a specification, not experimental results.

Prepared 2026-10-11. Source baseline: branch `ir-mlp-baseline-gpu`; read-only existing implementations under `experiments/mechanism_suite/`, `experiments/candidate_screen/`, `experiments/double_estimator/`, and `experiments/gpu_mlp/`. The separately running non-ridge and actor jobs must not be modified.

## Research question and scope

**Target:** Repair *full-history MLP's own* nonlinear noise-to-bias amplification for gradient-dependent semilinear PDEs. The primary main-table comparison is **unmodified full-history MLP vs one complete repaired MLP**. Pathwise terminal estimator, level-matched source baseline, double estimator and learned actor appear primarily in mechanism/ablation experiments. Actor is a candidate principal replacement for double, **not yet empirically validated**.

For B1, B3, B4, B5, and B6, use the existing convention

\[
u_t+\Delta_x u+f(u,z)=0,\qquad z=\sqrt{2}\nabla_x u,\quad
u(T,x)=g(x),\quad x\in\mathbb R^d.
\]

**B2 deliberately reproduces the original Hutzenthaler--Nguyen (2025) counterexample unchanged**: it uses \(\sigma=1\), \(\tfrac12\Delta\), \(T=1\), and \(z=\nabla u\), rather than the \(\sigma=\sqrt2\) convention above. A code adapter may rescale it, but published-benchmark reproduction must validate mathematical equivalence and explicitly report any rescaling.

Experiment axes: ambient dimension \(d\in\{100,400\}\), effective dimension \(k\), Picard depth \(n\), sampling parameter \(M\), and *the same* total-work accounting as current experiments (generator calls, Gaussian-coordinate work, and wall time, including actor training when used). Do **not** infer full-rank difficulty solely from a syntactically full-rank terminal. Before comparing solvers report the heat-only \(f=0\) discrepancy, spatial variability of \(f(u^*,z^*)\), gradient covariance effective rank, and Hessian effective rank.

A \(k\)-ridge solution means \(u(t,x)=v(t,B^\top x)\) for fixed \(B\in\mathbb R^{d\times k}\), \(B^\top B=I_k\). Multi-ridge with \(k=10\) or 20 is **still ridge**, but not one-dimensional; compare it alongside genuinely full-rank tests. Full rank need not imply a hard benchmark: high-dimensional concentration can make the nonlinear correction nearly spatially constant.

## Benchmark overview

| ID | Equation / test | Geometry | Ground truth | Present status |
| --- | --- | --- | --- | --- |
| B1 | Coupled multi-ridge quadratic HJ, \(f(z)=-\|z\|^2/2\) | \(k=1,5,10,20,40\), \(d=100,400\) | **Analytic** Hopf–Cole exponential sum, value and gradient | Extends P1/MR; new coupled \(k>1\) inputs needed |
| B2 | **HN (2025) original counterexample + analytic \(k\)-dimensional extension**, \(H_{d,k}(z)=\|z_{k+1:d}\|_2\) | \(k=1,10,30,50\); \(d=50,100,200,400\); \(k=d=50\) negative control | **Exact** sum of Gaussian absolute-value expectations and gradients | \(k=1\) published original; \(k>1\) new extension, original HN lower bound not automatically inherited |
| B3 | Smooth cubic, \(f(z)=-\kappa(\|z\|/\sqrt2)^3\) | Radial \(k=10,20\), \(d=100,400\) | Independently refined **1D radial** PDE, not yet produced | Existing C3 has only \(k=1\) reference |
| B4 | **Fully nonlinear** controlled-volatility HJB | Hessian-dependent, \(k=2\)–5 active coordinates | Independently refined **1D component** PDEs, not yet produced | **Outside** current \((u,z)\)-only MLP; conditional extension |
| B5 | Two-actor zero-sum Isaacs LQ game | Coupled \(k=10,20\) subspace with two control blocks | **Analytic** finite-horizon matrix Riccati | Existing C2 uses \(k=1\); new coupled game needed |
| B6 | Full-rank, nonseparable quadratic HJ | \(k=d=100,400\) | **Analytic** Hopf–Cole exponential sum | New terminal/gradient and cost calibration needed |

**Scope warning:** B2 at \(k=1\) is the original published HN counterexample; \(k>1\) is our analytic generalization, not a published HN theorem. It uses a different diffusion convention from the other benchmarks. B1 and B6 have the *same differential operator* and different terminal geometries; these are six distinct **benchmark instances**, not six mathematically distinct PDE operators. B4 does **not** support a claim that existing gradient-only full-history MLP already solves general fully nonlinear HJBs.

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

## B2 — Hutzenthaler--Nguyen (2025) and an analytic multi-active-dimensional extension

**Original published source:** M. Hutzenthaler and T. A. Nguyen, *Full history recursive multilevel Picard approximations suffer from the curse of dimensionality for the Hamilton--Jacobi--Bellman equation of a stochastic control problem* (2025), arXiv:2506.23969, https://arxiv.org/abs/2506.23969.

### Frozen experiment grid

Use **effective dimensions** \(k\in\{1,10,30,50\}\) and **ambient dimensions** \(d\in\{50,100,200,400\}\). Keep *all* 16 combinations in reporting; classify the corner \(k=d=50\) separately as a **degenerate \(f\equiv0\) negative control** rather than a successful nonlinear-failure repair. There are 15 nondegenerate \(k<d\) cells.

| effective \(k\) / ambient \(d\) | 50 | 100 | 200 | 400 |
| --- | --- | --- | --- | --- |
| 1 | original HN | original HN | original HN | original HN |
| 10 | extension | extension | extension | extension |
| 30 | extension | extension | extension | extension |
| 50 | **\(k=d\), null-driver negative control** | extension | extension | extension |

### PDE and terminal condition

**Preserve HN's original normalization:** \(T=1\), diffusion generator \(\frac12\Delta\), \(\sigma=I_d\), and state \(z=\nabla_x u\). For any \(1\le k\le d\), define

\[
\boxed{
\begin{aligned}
\partial_t u_{d,k}(t,x)+\tfrac12\Delta_x u_{d,k}(t,x)
+\left(\sum_{j=k+1}^{d}|\partial_{x_j}u_{d,k}(t,x)|^2\right)^{1/2}&=0,\\
u_{d,k}(1,x)&=g_k(x)=\frac{1}{\sqrt{k}}\sum_{j=1}^k|x_j|.
\end{aligned}}
\]

The Euclidean norm of an empty coordinate block is defined to be zero (so \(k=d\) gives \(H_{d,d}\equiv0\)). When **\(k=1\)** this *exactly* reduces to the original HN PDE and terminal condition \(g(x)=|x_1|\). **Only these \(k=1\) runs are replications of the published HN counterexample/theorem.** For \(k>1\), this is our newly constructed extension; the published HN lower bound is not claimed to hold without a new proof.

The Hamiltonian is

\[
H_{d,k}(z)=\|z_{k+1:d}\|_2
=\sup_{\substack{a_{1:k}=0\\\|a\|_2\le1}}a^\top z.
\]

It is convex and globally 1-Lipschitz for \(k<d\), and nondifferentiable where \(z_{k+1:d}=0\). Any admissible actor is pointwise optimal on the true solution; in particular \(a\equiv0\) is exact. Consequently B2 tests **MLP's intrinsic noise-to-bias amplification**, not whether a learned actor is needed to solve this intrinsically reduced problem.

### Exact analytic value, gradient, and simple checks

Let \(\tau=1-t>0\), \(G\sim\mathcal N(0,1)\), and write

\[
\psi_\tau(s)=\mathbb E|s+\sqrt{\tau}G|
=\sqrt{\frac{2\tau}{\pi}}\exp\left(-\frac{s^2}{2\tau}\right)
+s\left(2\Phi\left(\frac{s}{\sqrt{\tau}}\right)-1\right).
\]

Since \(u^*_{d,k}\) depends only on \(x_1,\ldots,x_k\), the nonlinear Hamiltonian vanishes **on the true solution**, and the exact solution is

\[
u^*_{d,k}(t,x)=\frac1{\sqrt{k}}\sum_{j=1}^k\psi_{1-t}(x_j).
\]

For \(t<1\), its gradient is

\[
\partial_{x_j}u^*_{d,k}(t,x)=
\begin{cases}
\displaystyle\frac{2\Phi(x_j/\sqrt{1-t})-1}{\sqrt{k}},&1\le j\le k,\\
0,&k<j\le d.
\end{cases}
\]

The \(k\times k\) active Hessian is diagonal with strictly positive entries for \(t<1\), hence has rank \(k\). This extension is **additively separable** (not a general coupled \(k\)-dimensional solution), while the original \(k=1\) is an actual 1D ridge.

At the origin,
\[
u^*_{d,k}(0,0)=\sqrt{\frac{2k}{\pi}},
\qquad\nabla u^*_{d,k}(0,0)=0.
\]
The \(k=1\) value is \(\sqrt{2/\pi}\). At terminal time the absolute-value payoff is not differentiable at zero coordinates, but Gaussian terminal points avoid those kinks almost surely; for terminal pathwise estimates use the a.e. gradient
\(\nabla g_k(x)=k^{-1/2}(\operatorname{sgn}(x_1),\ldots,\operatorname{sgn}(x_k),0,\ldots,0)\), with an explicit convention on the null set.

### Quantities to report and numerical interpretation

1. **Reproduce the published result for \(k=1\):** compare against the source's specified full-history recursion and scaling, including its joint-state \((u,\nabla u)\) error. When checking the original theorem use its published sampling schedule (e.g. \(M=n\)), not just the convenient \(n,M\) values of our GPU runs. Do not claim the original theorem for the \(k>1\) family.
2. **Main experimental comparison:** unmodified full-history MLP against **one complete repaired MLP**, across the frozen \(4\times4\) grid at comparable depth/sampling and matched computational work, with paired RNG and repeated roots. The \(k=d=50\) run is reported as a diagnostic negative control, not a nonlinearity-repair success.
3. **Error metrics:** value absolute RMSE and value RMSE normalized by the standard deviation of the exact value over fixed test points; report the **full gradient RMSE** and a joint-state RMSE to match HN's theorem, without promising that a value-only correction repairs every output gradient component.
4. **Potential confounding when \(k\) changes:** the normalization \(k^{-1/2}\) keeps \(\|\nabla g_k\|_2=1\) almost everywhere, **but does not keep the terminal value scale fixed**: \(u^*_{d,k}(0,0)=\sqrt{2k/\pi}\). Standard BEL gradient-label variances can consequently vary with \(k\) for reasons beyond the dimension \(d-k\). Report the empirical first-layer score variance and terminal mean/variance per \(k\); do not attribute every observed \(k\) trend to recursion alone. Use the same normalized terminal family for all methods, and never retune after inspecting method results.
5. **Comparators:** the known exact active-subspace projection and \(a\equiv0\) are legitimate, extremely strong structure-aware negative controls. A learned actor achieving \(a\approx0\) here does not demonstrate neural prediction skill. Double's polynomial-specific guarantee does not apply directly to this norm Hamiltonian.
6. **Scientific scope:** distinguish *full-history MLP failure on a formally high-dimensional PDE* from *high intrinsic solution complexity*. Since \(H(\nabla u^*)=0\), the problem is an exact noise-amplification stress test rather than an example with an active nonlinear contribution along the true solution. The multi-\(k\) extension remains analytically solvable.

The B2 protocol is separate from the unchanged B1/B3/B4/B5/B6 specifications and from ongoing actor/nonridge experiments.


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

1. **Before algorithms:** verify terminal match, PDE residual, analytic or independently refined reference, and gradient via independent finite differences. Log all seeds, \(d,k\), basis generation, terminal coefficients and time horizon. For B3 radial references check convergence near \(r=0\) and at the outer boundary; report reference uncertainty separately from solver error. B2 instead uses an analytic Gaussian-absolute-value reference and requires a faithful reproduction of the original HN full-history recursion.
2. **Avoid the easy-benchmark trap:** report \(\|u^*-u_{\mathrm{heat}}\|\), spatial standard deviation of \(f(z^*)\), \(\|\nabla u^*\|\), gradient-covariance effective rank, Hessian effective rank, and any exact symmetries. Examine whether the driver is nearly constant in space.
3. **Primary comparison:** vanilla full-history MLP versus **one complete corrected MLP** under matched \(n,M,d,k\), repeated seeds and matched wall time/operation counts. Do not fill the main benchmark table with separate path/baseline/double/actor variants. Explain actor training overhead; do not compare a trained method on inference-only cost against an untrained baseline.
4. **Ablations:** remove pathwise, remove level-matched baseline, remove actor, and compare actor versus double *when the Hamiltonian admits the double estimator*. Double is not a universal default for B2/B3; B2 additionally needs a published-recurrence and full-state error comparison. Report directly measured value bias and variance over independent repetitions, not just a slope-based inference.
5. **Targets are hypotheses:** \(n=3\) beating the finite-depth \(n=2\) accuracy floor; strong performance as \(k\) rises; dimension robustness when \(d-k\) is large; and actor's additional benefit beyond already strong path+baseline. No result is claimed ahead of experiments.
6. **Qualification:** B1/B2/B3/B5/B6 are semilinear in the highest-order term and are the primary scope. B4 is explicitly conditional pending a new fully nonlinear algorithm and should never be used as evidence for an unimplemented feature.

## Immediate implementation plan

- **First:** B1 coupled \(k\)-ridge with analytic formula; use the P1 code shape, extending vectors to a \(K\times k\) coefficient matrix and \(d\times k\) orthonormal basis. Verify \(k=1\) exactly reproduces P1 on matched inputs. Do *not* overwrite current P1/MR results.
- **Second:** B5 analytic coupled game (Riccati) and B6 full-rank nonseparable LSE, with nonlinearity-strength calibration.
- **Third:** reproduce HN 2025 at \(k=1\), then evaluate its frozen \(k\in\{1,10,30,50\}\), \(d\in\{50,100,200,400\}\) analytic extension (\(k=d=50\) is a null-driver control), preserving \(T=1\), \(\tfrac12\Delta\), and joint-state error. Independently build B3's radial reference for \(k=10,20\). Add new GPU PDE families without modifying frozen baselines.
- **Last / conditional:** B4 controlled-volatility / Hessian extension; do not make the current semilinear paper's correctness depend on this stretch benchmark.

Related in-repository code: `experiments/mechanism_suite/equations.py` (P1), `experiments/mechanism_suite/norm_hjb.py` (legacy P4, **not** B2), `experiments/candidate_screen/candidates.py` (C2/C3), `experiments/double_estimator/equations.py` (MR/game registry), `experiments/gpu_mlp/solver.py` (GPU).
