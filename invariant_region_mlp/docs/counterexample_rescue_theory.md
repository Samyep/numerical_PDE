# Structural rescue theorem and limitations

Let g_d(x)=|x_1| and f_d(v)=||v_(2:d)||_2, with u_t+(1/2)Delta u+f_d(grad u)=0 and u(1,x)=g_d(x). Work in the unique continuous at-most-linear-growth solution class. The solution is smooth for t<1.

## Lemma: certificate from symmetry
For any a with a_1=0, v(t,x)=u(t,x+a) has identical PDE and terminal data. Uniqueness gives v=u. Thus grad u belongs to S=span(e_1). The restricted PDE is the one-dimensional heat equation, with u(t,x)=E|x_1+sqrt(1-t)G|. This argument certifies the subspace before any numerical reference is evaluated.

## Proposition: collapse of the projected recurrence
Take exactly the centered terminal block and inverse-density weights in Hutzenthaler--Nguyen (2025), Eq. (6). Replace every recursive return Y=(U,V) by (U,QV), Q=diag(1,0,...,0). Initial states are zero. Since f_d(Qv)=0 for every v, induction over depth makes every nonlinear multilevel difference zero pathwise. The return is the projected terminal Monte Carlo block.

## Theorem: exact full-state MSE at the origin
For N=m^n independent G_i~N(0,1), the surviving coordinates are U_N=N^-1 sum |G_i| and V_1,N=N^-1 sum |G_i|G_i. The true state is (sqrt(2/pi),0,...,0). Therefore Var(U_N)=(1-2/pi)/N and Var(V_1,N)=3/N. All other gradients vanish. Hence full Euclidean state MSE is (4-2/pi)/N for every d,n,m>=1.

This requires projection of the returned gradient. Generator-only projection leaves d-1 terminal-gradient variances of 1/N each, hence MSE (d+3-2/pi)/N. Final-only projection cannot change U_N of the original nonlinear recursion.

## Cost corollary
After proving the corrections identically zero, prune them. N=ceil((4-2/pi)/epsilon^2) suffices, at O(dN) cost with ambient Brownian simulation or O(N+d) with active-coordinate simulation and dense output. An arbitrary N can be realized by n=1,m=N. For diagonal n=m, the smallest n with n^n>=A=(4-2/pi)/epsilon^2 has n^n<e*n*A, when the preceding diagonal budget is below A. Thus a conservative diagonal bound has an additional logarithmic factor. The unpruned matched-tree implementation is NOT charged only N samples.

## Theorem: general subspace commutation
Let B in R^(d x r) satisfy B^TB=I_r and Q=BB^T. Assume g_d(x)=g_r(B^Tx) and f_d(y,Bv)=f_r(y,v). Use identity diffusion and the same full-history recursion, but project all returned gradients with Q. Couple reduced Brownian motion W^r=B^TW^d and use the same random-time tree. Then the projected d-dimensional state equals the embedding (U^r,BV^r) of the r-dimensional MLP result at (t,B^Tx), for every finite n,m.

Proof: the centered terminal difference depends only on B^Tx and B^T(W_1-W_t). Its projected gradient weight is B[B^T(W_1-W_t)]/(1-t). Inductively the child states are embeddings of reduced states, so the restriction assumption makes each scalar driver difference equal. The projected outer Brownian weight is the embedded reduced weight. Add terms. Euclidean norms of embedded states are unchanged because B is an isometry on R^r.

For f_r=0 and L-Lipschitz g_r, h=1-t<=1, the collapsed estimator has MSE <= L^2*r*(h+r+2)/N. Indeed D=g_r(a+sqrt(h)G)-g_r(a) obeys D^2<=L^2*h*||G||^2. Bound the value variance by L^2*h*r and the trace of gradient covariance by L^2 E||G||^4=L^2*r*(r+2), then divide by N. This controls retained rank, not ambient dimension.

## What is not proved
The source lower bound remains valid for its original unprojected diagonal family: for every p>=0 and n>2p, liminf d^-p RMSE is infinite. This is not one fixed n outperforming every polynomial. Our method supplies certified structure and changes the algorithm. Reduced Monte Carlo is exactly as good on this family; the target is not intrinsically high-dimensional.

An arbitrary Euclidean Lipschitz extension is not automatically dimension-uniform in a maximum-norm MLP theorem. Nor is a nested full-history estimator automatically an unbiased-noise Picard oracle. Generic convex projection lemmas do not imply final-error dominance between two different recursive trajectories. The subspace commutation identity is the separate exact argument used here.
