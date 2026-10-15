# Strong-baseline protocol correction

The preceding pilot made Deep BSDE and SDGD-PINN easier by placing the new
Terminal-Aware construction inside both competitors. Those measurements remain
in the historical result directory, but they are excluded from the formal
comparison in this study.

## Deep BSDE

The pinned upstream DeepBSDE implementation at commit
`e76fed80b1995daf4a0dd6d3e2cf64a0931dabda` learns a free `y_init`, a free
`z_init`, and nonshared batch-normalized ReLU subnetworks for later time slices.
It does not add `g(x)` to `Y0` and does not add the terminal gradient to `Z`.

The common 960-state prediction task requires amortization over the initial
state. The corrected adaptation therefore replaces only the scalar `y_init`
with a raw network `Y0(x0)`. It retains nonshared BN/ReLU `Z_n(x)` networks and
the upstream division of their outputs by dimension. No network receives
`g(x)`, `grad g(x)`, or a terminal residual ansatz.

## SDGD-PINN

The pinned SDGD-PINN semilinear implementation at commit
`c597cc4fe9caac0b3ca118ea84453405f25ff833` uses an ordinary tanh MLP and adds
the soft terminal penalty

\[
20\,\operatorname{MSE}(u(T,X),g(X))
+\operatorname{MSE}(z(T,X),z_g(X)).
\]

The corrected adaptation follows that form. The network maps `(t,x)` directly
to `u`; it no longer uses `g(x)+(T-t)N(t,x)`. Terminal values and gradients
appear only as supervised boundary data in the authors' loss.

## Other baselines

Deep Picard's own terminal embedding and Full-History MLP's pathwise terminal
evaluation remain intact because they are parts of those published methods.
SCaSML is recomputed using the corrected ordinary SDGD-PINN as its backbone.

The final tables label the invalid older Deep BSDE, SDGD-PINN, and SCaSML rows
as historical and do not let them enter a formal accuracy-compute envelope.
