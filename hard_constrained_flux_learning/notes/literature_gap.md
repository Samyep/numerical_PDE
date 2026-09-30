# Literature / Positioning Checklist

## Current positioning

The project should **not** claim any of the following as new by themselves:

- learning numerical fluxes from data;
- conservative flux-form neural networks;
- entropy-aware learned fluxes;
- neural PDE solvers for shocks;
- hard constraints in neural numerical solvers.

The intended gap is narrower:

> an arbitrary learned flux proposal is passed through a cheap, differentiable, closed-form projection that guarantees the Tadmor interface entropy inequality exactly, while retaining the standard conservative update.

## Nearest-neighbor families to compare carefully

1. Conservative Flux-form Networks (CFN): learned local conservative schemes from trajectory data.
2. Entropy-stable CFN / structure-preserving learned flux methods: entropy stability through a specific scheme template.
3. NESCFN and other methods that jointly learn flux/entropy or encourage entropy structure numerically.
4. Neural Riemann solvers with hard physical constraints such as positivity, symmetry, consistency, invariance, or scaling.
5. Classical entropy-stable fluxes, entropy-conservative flux + dissipation constructions, algebraic flux correction, invariant-domain/positivity limiters.
6. Neural Conservation Laws / divergence-free architectures: exact conservation but a different architectural mechanism and problem emphasis.

## Questions the paper must answer

- Does direct hard projection add something that a fixed entropy-stable template does not?
- Does it preserve learned accuracy in smooth regions better than a strongly dissipative safe solver?
- How often and how strongly does the hard layer intervene?
- Can positivity/admissibility be enforced jointly and cheaply?
- Can the guarantee be made fully discrete rather than only semi-discrete?
- Does the learned proposal discover useful higher-order/local structure rather than collapsing to a classical flux?

## Literature-review standard

Before making a novelty claim, search specifically for combinations of:

- learned numerical flux + Tadmor projection;
- neural flux + half-space projection;
- certified entropy-stable neural flux;
- hard entropy constraints + finite volume + neural network;
- learned Riemann solver + entropy inequality projection.
