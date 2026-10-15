# Official baseline source review

Reviewed before implementing the adaptations in this experiment:

| Method | Primary source inspected | Frozen source revision | Adaptation used here |
|---|---|---|---|
| Deep BSDE | Han, Jentzen and E; authors' `frankhan91/DeepBSDE` | `e76fed80` | The original scalar \(Y_0\) becomes an amortized \(Y_0(x_0)\); a shared terminal-aware \(Z(t,x)\) predicts all 960 states. The BSDE update is \(dY=-f(Y,Z)dt+Z\,dW\). |
| PINN + SDGD | Hu et al.; authors' `zheyuanhu01/SDGD_PINN` | `c597cc4f` | Hard terminal embedding and two independent coordinate subsets give an unbiased gradient of the squared full PDE residual. |
| Deep Picard Iteration | Han et al.; authors' `DeepOptimalControl/DeepPicardIteration` | `e8616c2f` | PISGradNet-style terminal embedding, value and gradient regression, centered Bismut control variate, warm reload, up to 40 rounds. Paths are reduced from 4096 to 512 to fit the common GH200 budget. |
| Full-History MLP | Neufeld and Wu gradient-dependent MLP paper; the public full-history implementation bundled with SCaSML | paper arXiv:2310.12545; code `c2412b94` | Batched GPU recursion with pathwise terminal gradients and domain truncation for the cubic driver. It has no learned parameters and therefore appears only on total-compute curves. |
| SCaSML | Fan et al.; authors' `Francis-Fan-create/SCaSML` | `c2412b94` | The compatible first defect/Picard correction is applied to each hard-terminal SDGD-PINN checkpoint. Its backbone training and residual-simulation inference are both charged. A Hutchinson trace estimate replaces an exact ambient Laplacian. |

All five formulations accept nonlinear dependence on both value and gradient. No implementation receives the independent spectral solution during training, validation, checkpoint selection, or inference-budget selection.
