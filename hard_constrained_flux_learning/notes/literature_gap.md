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


## 2026 re-audit after the 2D Euler/trust experiments

A targeted search was repeated after the framework had evolved beyond the
initial hard-projection idea.

### Closest current neighbors

#### HCNRS: hard-constrained neural Riemann solver (Zhang et al., 2026)

`Learning the Exact Flux: Neural Riemann Solvers with Hard Constraints`
(arXiv:2603.30007) is now the closest hard-constraint neighbor.

It learns a surrogate for the exact Riemann solver and enforces five constraints
by construction:

- positivity;
- consistency;
- mirror symmetry;
- Galilean invariance;
- scaling invariance.

Its target is the local exact Riemann map, and it demonstrates shallow-water
and Euler rollouts including a 2D Euler implosion problem.

Important distinction for HCFL:
- HCNRS does **not** use direct projection onto the Tadmor entropy half-space as
  its defining constraint;
- its primary supervised target is the exact Riemann solver, whereas HCFL is
  trained from coarse solution trajectories without flux labels;
- HCFL learns a local correction to an approximate classical flux and retains
  a separate hard entropy/admissibility safety stack.

This paper must be a primary baseline / related-work comparison.

#### ESCFN (Liu et al., J. Sci. Comput. 2026)

ESCFN embeds neural flux learning into a second-order non-oscillatory
Kurganov--Tadmor scheme. It provides an entropy-stable learned construction but
is tied to the chosen classical numerical template.

Distinction:
HCFL starts from an arbitrary learned proposal (or learned correction around a
chosen base flux) and enforces the known face entropy condition directly by
projection. The hard layer is therefore separable from the proposal
architecture / base solver.

#### NESCFN (Liu, Zhang, Gelb, JCP 2026)

NESCFN jointly learns a conservation law and convex entropy from trajectories.
The authors explicitly state that their current construction does not provide a
formal proof of entropy stability for general learned fluxes; their
central-flux-based entropy-conservative surrogate need not satisfy Tadmor's
shuffle condition exactly and structural conditions are enforced only to
numerical tolerance.

Distinction:
HCFL assumes the physical entropy pair is known and uses that knowledge to
enforce the face condition exactly. This is a less ambitious system-identification
setting but yields a stronger certified numerical-solver statement.

### OOD trust/fallback search

The targeted search did not identify a close neural hyperbolic-flux method that
uses the following specific construction:

1. calibrate a state-space support score using training trajectories only;
2. project both learned and classical fluxes into the same hard feasible set;
3. convexly interpolate toward the classical flux as the current state leaves
   training support.

This should **not** be claimed as globally novel without a broader literature
review, but it currently appears to be a useful differentiator from the closest
learned-Riemann / learned-flux papers.

### Current safest novelty statement

Avoid:
> first hard-constrained neural solver for hyperbolic PDEs.

That is false because HCNRS and other hard-constrained methods exist.

Avoid:
> first entropy-stable neural flux.

That is false because ESCFN and related work exist.

Current defensible target:
> a trajectory-trained conservative flux learner in which an arbitrary learned
> proposal is minimally projected into the known Tadmor entropy-feasible set,
> combined with local admissibility limiting and a training-calibrated convex
> fallback toward a classical entropy-feasible flux under severe distribution
> shift.

The novelty of this exact combination still requires final publication-level
checking before submission.
