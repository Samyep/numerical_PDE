# Why the counterexample has the invariant subspace span(e1)

Consider

\[
\partial_t u_d+\tfrac12\Delta u_d+
\left(\sum_{j=2}^d |\partial_{x_j}u_d|^2\right)^{1/2}=0,
\qquad u_d(1,x)=|x_1|.
\]

The key point is that the gradient constraint is certified by the HJB/control structure, not inferred from numerical ground truth.

## 1. Control representation of the Hamiltonian

For

\[
H(p)=\left(\sum_{j=2}^d p_j^2\right)^{1/2},
\]

Cauchy--Schwarz gives

\[
H(p)=\sup_{a_1=0,\ \|a\|_2\le1} a^\top p.
\]

Hence the PDE is the HJB equation for the controlled diffusion

\[
dX_s=a_s\,ds+dW_s,
\qquad a_{s,1}=0,
\qquad \|a_s\|_2\le1,
\]

with terminal reward \(g(X_1)=|X_{1,1}|\).

## 2. Why the value can only depend on x1

Every admissible control has zero first component, so the first state coordinate is completely uncontrolled:

\[
X_{s,1}=x_1+W_1(s)-W_1(t).
\]

The terminal reward also reads only this first coordinate. Therefore every admissible control produces exactly the same payoff distribution:

\[
\mathbb E|X_{1,1}|
=
\mathbb E|x_1+W_1(1)-W_1(t)|.
\]

Taking the supremum over controls changes nothing. Thus the value function is

\[
u_d(t,x)=\mathbb E|x_1+W_1(1)-W_1(t)|,
\]

which is independent of \(x_2,\ldots,x_d\). Consequently

\[
\partial_{x_2}u_d=\cdots=\partial_{x_d}u_d=0,
\qquad
\nabla u_d(t,x)\in \operatorname{span}(e_1).
\]

This is the exact structural certificate used by IR-MLP.

## 3. Why this is the HJB solution: verification intuition

Let \(w\) be a sufficiently smooth solution and let \(X^a\) follow any admissible control. Ito's formula yields

\[
dw(s,X_s^a)=
\left(w_t+\tfrac12\Delta w+a_s^\top\nabla w\right)ds
+\nabla w^\top dW_s.
\]

Because

\[
a_s^\top\nabla w
\le
\sup_{a_1=0,\|a\|\le1}a^\top\nabla w,
\]

the PDE implies that the drift above is nonpositive. Therefore

\[
w(t,x)\ge \mathbb E g(X_1^a)
\]

for every admissible control, so \(w\) lies above the control value. If one chooses a maximizing feedback direction for the Hamiltonian, the drift becomes zero and equality is attained. Hence any classical solution agrees with the value function. For the nonsmooth terminal datum \(|x_1|\), the rigorous version is the standard viscosity-solution dynamic-programming/comparison argument for this globally Lipschitz Hamiltonian.

## 4. Consequence for IR-MLP

The certified set is

\[
\mathcal C_d=\operatorname{span}(e_1),
\qquad
Q_d z=(z_1,0,\ldots,0).
\]

For every numerical gradient estimate \(z\),

\[
f_d(Q_dz)=0.
\]

Thus stochastic gradient noise in the \(d-1\) forbidden directions is removed before it can be converted into a positive nonlinear source. This is the exact mechanism behind the counterexample rescue theorem.
