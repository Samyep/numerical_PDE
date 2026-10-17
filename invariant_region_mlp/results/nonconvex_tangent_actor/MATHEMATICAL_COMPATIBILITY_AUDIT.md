# Mathematical Compatibility Audit

The implementation uses the existing generic tangent-coefficient and weighted Feynman–Kac code. No PDE-specific linear solver is introduced.

| Driver | max tangent error | max Hessian error | f_uu>0 fraction | f_uu<0 fraction | invariant [0,1] | gate |
|---|---:|---:|---:|---:|---|---|
| primary | 4.441e-16 | 0.000e+00 | 0.457 | 0.543 | True | PASS |

The interval check uses the forward-time form `w_tau = Delta w + f(w,sqrt(2) grad w)`: 0 is a subsolution and 1 is a supersolution for both drivers. The numerical reference independently remains inside this interval.

The local Taylor remainder follows from the bounded joint Hessian on the invariant interval and bounded Actor region. This is a local consistency statement; it is not used as a global convergence theorem.
