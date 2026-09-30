# Burgers trajectory-only flux identifiability

For a conservative two-point update
[
u_i^{n+1}=u_i^n-lambdaig(F(u_i,u_{i+1})-F(u_{i-1},u_i)ig),
]
suppose a teacher flux `G` generates the exact same one-step map for every local triplet `(a,b,c)`. Then
[
F(b,c)-F(a,b)=G(b,c)-G(a,b).
]
With `H=F-G`,
[
H(b,c)=H(a,b)
]
for all `a,b,c`. This implies `H` is constant on a connected state domain. If the learned flux is consistent,
[
F(u,u)=f(u)=G(u,u),
]
the constant is zero, so `F=G`.

The accompanying trajectory-only experiment uses Godunov trajectories as the teacher, never shows flux labels to the network, and empirically approaches the Godunov flux surface to roughly 1e-2 RMSE over `[-1.5,1.5]^2`.

This result is useful as a controlled identification check. The intended research setting is harder: fine-grid trusted trajectories are downsampled to a coarse grid, in which case the optimal learned coarse flux need not equal any classical coarse-grid flux.
