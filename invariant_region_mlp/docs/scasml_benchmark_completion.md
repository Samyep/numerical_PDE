# SCaSML benchmark completion for IR-MLP

## Goal

Complete the benchmark families used in the public SCaSML study and compare:

1. the authors' reported full-history MLP result;
2. our matched pure-MLP baseline;
3. IR-MLP with certified projection inside the recursion;
4. the authors' reported full SCaSML result.

The new standalone runs use Picard level n=2, Monte Carlo base M=10, 1,000 interior plus 200 boundary test points, and three paired repetitions. The main scientific track uses the standard Elworthy--Bismut--Li terminal-gradient normalization. The HJB row reuses the existing 10-repetition 100--160D current-round experiment and the non-oracle uniform radius sqrt(15).

## Main comparison

### Linear convection--diffusion (LCD)

| d | our MLP | IR-MLP | hard reduction | paper MLP | paper SCaSML |
|---:|---:|---:|---:|---:|---:|
| 10 | 0.2240 +/- 0.0099 | 0.2240 +/- 0.0099 | 0.0% | 0.227 | 0.0274 |
| 20 | 0.2398 +/- 0.0085 | 0.2398 +/- 0.0085 | 0.0% | 0.235 | 0.0472 |
| 30 | 0.2496 +/- 0.0218 | 0.2496 +/- 0.0218 | 0.0% | 0.238 | 0.0972 |
| 60 | 0.2581 +/- 0.0089 | 0.2581 +/- 0.0089 | 0.0% | 0.239 | 0.132 |

This is an expected mechanism control. The generator is f=0, so changing the recursive z state cannot feed back into u. The matched raw MLP is close to the paper's reported MLP values.

### Gradient-dependent nonlinear / viscous-Burgers benchmark (VB)

The public executable class currently uses sigma=0.25. We use that value in the main comparison because it is the setting that reproduces the reported MLP table; see the discrepancy notes below.

| d | corrected-EBL MLP | IR-MLP | hard reduction | paper MLP | paper SCaSML |
|---:|---:|---:|---:|---:|---:|
| 20 | 0.0531 +/- 0.0016 | **0.0270 +/- 0.0005** | 49.1% | 0.0836 | 0.00403 |
| 40 | 0.0713 +/- 0.0015 | **0.0414 +/- 0.0023** | 41.9% | 0.104 | 0.0292 |
| 60 | 0.0818 +/- 0.0036 | **0.0522 +/- 0.0011** | 36.2% | 0.117 | 0.0288 |
| 80 | 0.0928 +/- 0.0034 | **0.0629 +/- 0.0041** | 32.2% | 0.119 | 0.0564 |

For diagnosis, reproducing the public MLP_full_history.py terminal-gradient scaling gives unprojected errors 0.0732, 0.0946, 0.1030, and 0.1096 in 20/40/60/80D, respectively, which are closer to the paper MLP column. Thus the standard EBL correction already improves the baseline; the invariant-region projection then gives an additional 32--49% reduction relative to that corrected baseline.

IR-MLP does not beat the full SCaSML method here in general. The closest case is 80D: 0.0629 versus the paper SCaSML value 0.0564.

### Diffusion--reaction / oscillating-solution benchmark (DR)

| d | our MLP | IR-MLP | hard reduction | paper MLP | paper SCaSML |
|---:|---:|---:|---:|---:|---:|
| 100 | 0.08836 +/- 0.00485 | **0.06821 +/- 0.00316** | 22.8% | 0.0899 | 0.0111 |
| 120 | 0.09289 +/- 0.00174 | **0.06597 +/- 0.00096** | 29.0% | 0.0913 | 0.0103 |
| 140 | 0.09283 +/- 0.00458 | **0.06498 +/- 0.00347** | 30.0% | 0.0897 | 0.0300 |
| 160 | 0.09122 +/- 0.00492 | **0.06085 +/- 0.00305** | 33.3% | 0.0900 | 0.0322 |

This is the cleanest reproduction check: our unprojected MLP nearly matches the reported paper MLP at all four dimensions. Certified value/gradient projection reduces error by about 23--33%, but full SCaSML remains better.

### HJB--Rosenbrock (existing current-round result)

| d | our corrected MLP | uniform hard sqrt(15) | hard reduction | paper MLP | paper SCaSML |
|---:|---:|---:|---:|---:|---:|
| 100 | 1.527 | **0.830** | 45.6% | 5.63 | 0.0553 |
| 120 | 1.564 | **0.851** | 45.6% | 5.50 | 0.0666 |
| 140 | 1.568 | **0.864** | 44.9% | 5.37 | 0.0684 |
| 160 | 1.549 | **0.879** | 43.3% | 5.27 | 0.0994 |

The hard result uses only the public coefficient ranges, not the realized matrix. It strongly stabilizes pure MLP but is not a substitute for the trained-surrogate plus defect-correction SCaSML pipeline.

## Upstream reproducibility discrepancies

Three discrepancies matter for any direct numerical comparison.

1. **LCD test domain.** The paper text describes [0,0.5]^d. The public full-history LCD class does not override Equation.test_geometry(), so the executable test path samples [-0.5,0.5]^d. The latter is used here because it reproduces the paper MLP table.
2. **VB diffusion coefficient.** The current paper text prints sigma=sqrt(2), while the public Grad_Dependent_Nonlinear class returns sigma=0.25. The 0.25 setting is the one that tracks the reported MLP errors and is therefore the main comparison setting here.
3. **Terminal gradient normalization.** Public MLP_full_history.py computes the terminal z term by dividing by T-t. The standard EBL weight divides by sqrt(T-t). Main IR-MLP claims use the corrected EBL normalization; the public-code scaling is retained only as a diagnostic.

These differences mean the paper's SCaSML numbers should be treated as an external reference column, not as an assertion that every row is bit-for-bit identical to our corrected pure-MLP implementation.

## Verdict

The completed suite supports a narrower and stronger mechanism claim:

- LCD: no nonlinear feedback, so hard projection has essentially no effect.
- VB: gradient-dependent nonlinear feedback is present; hard projection reduces corrected-MLP error by 32--49%.
- DR: nonlinear value feedback is present; hard projection reduces error by 23--33%.
- HJB: quadratic gradient feedback is especially sensitive; the non-oracle hard ball reduces error by about 43--46%.

The method does **not** generally beat full SCaSML. The natural next experiment is therefore to insert the invariant-region projector into the official SCaSML defect recursion. For a surrogate state z_hat and defect z_breve, project the total state z_total=z_hat+z_breve and convert back to the projected defect. That experiment directly tests whether the two methods are complementary.
