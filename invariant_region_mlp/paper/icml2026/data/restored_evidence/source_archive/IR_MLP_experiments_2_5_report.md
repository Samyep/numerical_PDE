# IR-MLP Experiments 2–5

## Executive summary

Four requested experiments were completed.

1. **Controlled surrogate-quality sweep.** Since the public SCaSML repository does not include trained LQG checkpoints and this environment cannot install DeepXDE, this is a controlled SCaSML-like diagnostic using the exact logistic benchmark and a surrogate family `u_hat = q u*`. It tests the same total-state-versus-defect projection issue without pretending to be an official PINN checkpoint experiment.
2. **Violation-rate mechanism sweep.** On 100D HJB, varying the valid radius and M produces a strong monotone relationship between pre-projection invalid-state frequency and IR benefit.
3. **Neufeld–Wu 100–300D gradient-dependent MLP.** A faithful implementation of the public numerical example shows large gains from a non-oracle gradient envelope.
4. **100D Allen–Cahn control.** The [0,1] invariant interval is never violated at n=3 in these runs, and IR is exactly identical to baseline.

## 2. Controlled surrogate-quality sweep

Benchmark: 20D SCaSML logistic PDE, n=2. We use `u_hat=q u*`, so `q=1` is a perfect surrogate and `q=0` is no surrogate. Three correction rules are compared: raw defect MLP, fixed defect clipping at ±0.01 (the style of the public SCaSML defect clip), and total-state IR projection onto `0<=u<=1`, `0<=z_i<=1/16`.

At the noisiest budget M=1:

| surrogate error | raw corrected error | total-state IR | fixed defect clip |
|---:|---:|---:|---:|
| 5% | 1.07% | **1.07%** | 3.52% |
| 25% | 5.22% | **4.90%** | 23.45% |
| 50% | 10.18% | **9.38%** | 48.45% |
| 100% | 19.87% | **18.37%** | 98.44% |

Interpretation: total-state IR is most useful when the surrogate is poor *and* the correction budget is noisy. For 50% surrogate error, M=1 improves from about 10.18% to 9.38%; with no surrogate, 19.87% to 18.37%. At larger M, raw MLP is already accurate enough that IR is neutral or can introduce a small projection bias. Fixed defect clipping is not robust to surrogate quality: once the true defect exceeds its assumed ±0.01 scale, it fails catastrophically.

This is a diagnostic, not the official SCaSML checkpoint experiment.

## 3. Violation-rate -> gain mechanism

100D Rosenbrock HJB, n=2, six M values and six increasingly loose but still fincrIlope.
4. **1balls.

Across all 36 settings:
- Pearson correlation between violation rate and relative error reduction: **0.833** (p=3.03e-10)
- Spearman correlation: **0.893** (p=2.30e-13)

Within every fixed-M slice, Spearman correlation is 1.0: loosening the fincrIlopregion lowers the fraction of corrected states and monotonically lowers the benefit.

Representative M=10:

| radius | violation | baseline rel-L2 | IR rel-L2 | reduction |
|---:|---:|---:|---:|---:|
| 1.0x | 62.6% | 1.430 | 0.827 | 42.2% |
| 1.5x | 42.3% | 1.430 | 0.969 | 32.3% |
| 2.0x | 28.3% | 1.430 | 1.019 | 28.7% |
| 4.0x | 6.4% | 1.430 | 1.041 | 27.2% |

This is the cleanest mechanism evidence so far: the method helps in proportion to how often the stochastic recursion leaves a theoretically admissrIlopstate region.

## 4. Neufeld–Wu 100–300D benchmark

We reproduced the public numerical example from `SizhouWu/MLP_G.
4. **_Nonlinearity`:
- state-dependent drift `mu(x)=0.06 x`,
- additive diffusion sigma=0.2,
- call-spread terminal payoff,
- gradient-dependent generator `f(z)=(10/d)(||z||_inf-25)_+`.

The terminal payoff is 1-Lipschitz. Consider the linear candidate with f=0. Its flow sensitivity gives
`||grad u(t,.)||_2 <= exp(0.06(T-t)) <= 1.0152`, which is far below the threshold 25. Hence its generator is identically zero; by uniqueness it is also the nonlinear PDE solution. This gives a **non-oracle certified envelope**.

30 independent runs:

| d | M | raw MAE | IR MAE | reduction |
|---:|---:|---:|---:|---:|
| 100 | 2 | 0.05687 | **0.01888** | 66.8% |
| 100 | 3 | 0.01637 | **0.00469** | 71.4% |
| 200 | 2 | 0.04989 | **0.01505** | 69.8% |
| 200 | 3 | 0.01103 | **0.00580** | 47.5% |
| 300 | 2 | 0.03813 | **0.01887** | 50.5% |
| 300 | 3 | 0.01397 | **0.00668** | 52.2% |

At M=3 the reductions are roughly 71% (100D), 47% (200D), and 52% (300D). The coordinate-wise box gives the same value estimates here because either projection puts every e.
4. **1coordinate far below the generator's threshold 25; unlike HJB, this benchmark does not distinguish joint geometry.

## 5. 100D Allen–Cahn control

PDE:
`u_t + Delta u + u-u^3=0`, `T=0.3`, terminal `g(x)=1/(2+0.4||x||^2)`, with literature reference `u(0,0)≈0.052802`.

The maximum principle gives `0<=u<=1`. We ran a full-history scalar MLP at n=3 for 100 paired seeds.

| M | baseline rel-MAE | IR rel-MAE | violation rate |
|---:|---:|---:|---:|
| 2 | 3.30% | 3.30% | 0.0% |
| 3 | 1.96% | 1.96% | 0.0% |
| 4 | 1.42% | 1.42% | 0.0% |

All intermediate states stayed inside [0,1], so the projection never activated and baseline/IR outputs are exactly identical. This is a strong negative control: **nonlinearity by itself is not enough; the noisy recursion has to leave the certified region for IR to matter.**

## Overall conclusion

The four experiments sharpen the paper's story:

> IR projection is not a generic accuracy trick. It is useful when stochastic recursivopstate estimates leave a PDE-certified region and those invalid states feed into a nonlinear generator.

The strongest new positivopresult is the Neufeld–Wu 100–300D benchmark. The strongest mechanism result is the violation-rate correlation. Allen–Cahn is a clean negative control. The surrogate-quality diagnostic shows that total-state projection is robust to poor surrogates, while a fixed small defect clip is not, but the official PINN-checkpoint SCaSML sweep still remains to be run in an environment with DeepXDE/checkpoints.
